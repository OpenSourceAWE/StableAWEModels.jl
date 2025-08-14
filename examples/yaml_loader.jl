using YAML
using Logging
using SymbolicAWEModels
using LinearAlgebra
using VortexStepMethod, ControlPlots, KiteUtils
using Plots, PlotlyJS
plotlyjs()

# Default diameter for wing section segments (in mm)
DEFAULT_WING_DIAMETER_MM = 1.0

# ------------------ helpers ------------------

function parse_table(tbl)::Vector{NamedTuple}
    haskey(tbl, "headers") || throw(ArgumentError("table is missing `headers`"))
    haskey(tbl, "data")    || throw(ArgumentError("table is missing `data`"))
    headers = String.(tbl["headers"])
    rows = tbl["data"]
    out = NamedTuple[]
    for (k, row) in enumerate(rows)
        # skip empty or comment rows
        if isempty(row) || (isa(row[1], String) && startswith(row[1], "#"))
            continue
        end
        # allow missing trailing columns (fill with nothing)
        if length(row) < length(headers)
            row = vcat(row, fill(nothing, length(headers) - length(row)))
        end
        if length(row) > length(headers)
            @warn "Skipping row $k in table: has $(length(row)) values, expected $(length(headers)). Row: $row"
            continue
        end
        nt = NamedTuple{Tuple(Symbol.(headers))}(Tuple(row))
        push!(out, nt)
    end
    return out
end
"""
    neutral_transform(idx::Int; points::Vector{Point}, base_point_idx::Int=1,
                           rot_point_idx::Int=2, base_pos::Vector{Float64}=[0.0,0.0,0.0])
    This enables no transformation, but still a call to the transformations

"""
function neutral_transform(idx::Int; points, base_point_idx::Int=1,
                           rot_point_idx::Int=2, base_pos::AbstractVector{<:Real}=nothing)
    if base_pos === nothing
        base_pos = [0.0, 0.0, 0.0]
    end
    dir0 = points[rot_point_idx].pos_cad .- points[base_point_idx].pos_cad
    if norm(dir0) < 1e-12
        error("neutral_transform: base_point_idx and rot_point_idx coincide; pick different indices")
    end
    curr_elev   = KiteUtils.calc_elevation(dir0)
    curr_azim   = -KiteUtils.azimuth_east(dir0)
    return SymbolicAWEModels.Transform(idx, curr_elev, curr_azim, 0.0;
        base_pos=base_pos, base_point_idx=base_point_idx, rot_point_idx=rot_point_idx)
end


"""
        load_sys_struct_from_yaml(yaml_path::AbstractString; system_name="from_yaml", set=nothing)

Build a `SymbolicAWEModels.SystemStructure` from a YAML geometry file.

# YAML expectations (example)
bridle_nodes:
    headers: ["id","x","y","z","type"]   # `type` ∈ {"static","dynamic","knot","pulley"}
    # type mapping: static/pulley → STATIC, dynamic/knot → DYNAMIC
    data:    [[1,0,0,0,"static"], [2,1,0,0,"dynamic"]]
bridle_lines:
    headers: ["name","rest_length","diameter","material","density"]
    data:    [["A", 1.2, 0.002, "dyneema", 970.0]]
bridle_connections:
    data:    [["A", 1, 2, 0]]
    # If present, ck is validated and logged, but currently ignored in segment creation.
"""
function load_sys_struct_from_yaml(yaml_path::AbstractString; system_name="from_yaml", set=nothing())
    # Load YAML data
    data = YAML.load_file(yaml_path)

    # settings
    set === nothing && (set = SymbolicAWEModels.load_settings("base"))

    # --- Points ---
    points = SymbolicAWEModels.Point[]
    node_rows = parse_table(data["bridle_nodes"])
    for r in node_rows
        # Accept both 'static'/'dynamic' and 'knot'/'pulley' as types
        kind = SymbolicAWEModels.DYNAMIC
        if haskey(r, :type)
            t = lowercase(String(r.type))
            if t == "static"
                kind = SymbolicAWEModels.STATIC
            elseif t == "dynamic"
                kind = SymbolicAWEModels.DYNAMIC
            elseif t == "knot"
                kind = SymbolicAWEModels.DYNAMIC
            elseif t == "pulley"
                kind = SymbolicAWEModels.STATIC
            end
        end
        push!(points, SymbolicAWEModels.Point(
            Int(r.id),
            [Float64(r.x), Float64(r.y), Float64(r.z)],
            kind
        ))
    end
    # warn on duplicate IDs
    ids = getfield.(points, :idx)
    (length(unique(ids)) == length(ids)) || @warn "Duplicate node ids detected in bridle_nodes."

    # --- Wing Sections Points (LE/TE) ---
    if haskey(data, "wing_sections")
        wing_sections_data = data["wing_sections"]
        headers = wing_sections_data["headers"]
        for row in wing_sections_data["data"]
            section = Dict(zip(headers, row))
            # Add LE point (idx will be set later)
            push!(points, SymbolicAWEModels.Point(
                0,
                [Float64(section["LE_x"]), Float64(section["LE_y"]), Float64(section["LE_z"])],
                SymbolicAWEModels.DYNAMIC
            ))
            # Add TE point (idx will be set later)
            push!(points, SymbolicAWEModels.Point(
                0,
                [Float64(section["TE_x"]), Float64(section["TE_y"]), Float64(section["TE_z"])],
                SymbolicAWEModels.DYNAMIC
            ))
        end
    end

    # Reindex all points so that point.idx == i
    for (i, pt) in enumerate(points)
        points[i] = SymbolicAWEModels.Point(i, pt.pos_cad, pt.type)
    end

    # --- Lines ---
    line_rows = parse_table(data["bridle_lines"])
    # warn on duplicate names, skip rows with missing name
    names = String[]
    line_dict = Dict{String,Any}()
    for r in line_rows
        if !haskey(r, :name) || r.name === nothing
            @warn "Skipping bridle_lines row with missing name: $r"
            continue
        end
        push!(names, String(r.name))
        line_dict[String(r.name)] = r
    end
    (length(unique(names)) == length(names)) || @warn "Duplicate line names detected in bridle_lines."

    # --- Connections ---
    conn_data = []
    if haskey(data["bridle_connections"], "headers")
        conn_rows = parse_table(data["bridle_connections"])
        for r in conn_rows
            # Accept both 3- and 4-column connections
            if haskey(r, :ck) && r.ck !== nothing
                push!(conn_data, (String(r.name), Int(r.ci), Int(r.cj), Int(r.ck)))
            else
                push!(conn_data, (String(r.name), Int(r.ci), Int(r.cj), 0))
            end
        end
    else
        for (k, row) in enumerate(data["bridle_connections"]["data"])
            # skip empty rows
            if isempty(row)
                @warn "Skipping empty bridle_connections row $k"
                continue
            end
            # skip comment rows (first element must be a String)
            if isa(row[1], String) && startswith(row[1], "#")
                @warn "Skipping comment bridle_connections row $k: $row"
                continue
            end
            # Accept both 3- and 4-column connections
            if length(row) == 4
                push!(conn_data, (String(row[1]), Int(row[2]), Int(row[3]), Int(row[4])))
            elseif length(row) == 3
                push!(conn_data, (String(row[1]), Int(row[2]), Int(row[3]), 0))
            else
                @warn "Skipping malformed bridle_connections row $k: $row"
            end
        end
    end

    # # basic index validation
    n_pts = length(points)
    #     for (k, (_, ci, cj, _)) in enumerate(conn_data)
    #         if ci == 0 || cj == 0 || ci < 1 || cj < 1 || ci > n_pts || cj > n_pts
    #             @warn "Skipping bridle_connections row $k: references out-of-range point index ci=$ci cj=$cj (n_pts=$n_pts)"
    #             continue
    #         end
    # end

    # --- Segments ---
    segments = SymbolicAWEModels.Segment[]
    
    for (i, (name, ci, cj, ck)) in enumerate(conn_data)
        if ck != 0 && (ck < 1 || ck > n_pts)
            @warn "bridle_connections row $i: ck index $ck is out of range (n_pts=$n_pts); ignored."
        end
        line = get(line_dict, name, nothing)
        if line === nothing
            @warn "Line $name not found for connection $i"
            continue
        end
        # NOTE: assuming YAML diameter is in meters; convert to mm for Segment
        push!(segments, SymbolicAWEModels.Segment(
            i, set, (ci, cj), SymbolicAWEModels.POWER_LINE;
            l0 = Float64(line.rest_length),
            diameter_mm = Float64(line.diameter) * 1000,
            compression_frac = 0.01
        ))
        # `ck` is intentionally ignored for now (documented behavior)
    end

    # --- Wing Sections Segments (LE-TE) ---
    if haskey(data, "wing_sections")
        wing_sections_data = data["wing_sections"]
        n_sections = length(wing_sections_data["data"])
        # Find the indices of the LE/TE points just added (they are always the last 2*n_sections points)
        first_le_idx = length(points) - 2 * n_sections + 1
        for k in 1:n_sections
            le_idx = first_le_idx + 2 * (k - 1)
            te_idx = first_le_idx + 2 * (k - 1) + 1
            # Defensive: check indices
            if le_idx < 1 || te_idx < 1 || le_idx > length(points) || te_idx > length(points)
                @warn "Skipping wing_sections segment: invalid LE/TE indices le_idx=$le_idx te_idx=$te_idx"
                continue
            end
            l0 = norm(points[le_idx].pos_cad .- points[te_idx].pos_cad)
            # Try to get diameter from YAML if present, else use default
            diameter_mm = DEFAULT_WING_DIAMETER_MM
            if haskey(wing_sections_data, "headers") && "diameter" in wing_sections_data["headers"]
                diam_idx = findfirst(==("diameter"), wing_sections_data["headers"])
                diam_val = wing_sections_data["data"][k][diam_idx]
                if diam_val !== nothing
                    diameter_mm = Float64(diam_val) * 1000
                end
            end
            push!(segments, SymbolicAWEModels.Segment(
                length(segments) + 1, set, (le_idx, te_idx), SymbolicAWEModels.POWER_LINE;
                l0 = l0,
                diameter_mm = diameter_mm,
                compression_frac = 0.01
            ))
        end
    end

    # --- Transforms ---
    # Use neutral_transform to ensure no net rotation/translation
    transforms = [neutral_transform(1; points=points, base_point_idx=1, rot_point_idx=2, base_pos=Vector{Float64}(points[1].pos_cad))]

    return SymbolicAWEModels.SystemStructure(system_name, set; points, segments, transforms)
end

# ------------------- Example usage -------------------

### Change to your model name
# model_name = "pyramid_model"
model_name = "TUDELFT_V3_KITE" 

println("\n\nYAML GEOMETRY LOADER\n", "="^40)
set = SymbolicAWEModels.load_settings(model_name)  # Loads as Dict
geometry_path = joinpath(dirname(@__DIR__), "data", model_name, "wing_geometry.yaml")
sys = load_sys_struct_from_yaml(geometry_path; system_name=model_name, set=set)

function plot3d_saddle(points, segments; title::AbstractString="3D Plotly Plot")
    # Use pos_cad if pos_w is not available
    get_pos(p) = hasproperty(p, :pos_w) ? p.pos_w : p.pos_cad
    x = [get_pos(p)[1] for p in points]
    y = [get_pos(p)[2] for p in points]
    z = [get_pos(p)[3] for p in points]
    plt = Plots.scatter3d(x, y, z; markersize=2, markerstrokewidth=0, title=title, xlabel="X (m)", ylabel="Y (m)", zlabel="Z (m)", legend=false)
    for s in segments
        i, j = s.point_idxs
        if i == 0 || j == 0
            @warn "Segment with zero index detected: $(s)"
            continue
        end
        Plots.plot3d!([x[i], x[j]], [y[i], y[j]], [z[i], z[j]]; alpha=1, linewidth=1, color=:black)
    end
    fixed_idx = [i for (i, pt) in enumerate(points) if pt.type == SymbolicAWEModels.STATIC]
    !isempty(fixed_idx) && Plots.scatter3d!(x[fixed_idx], y[fixed_idx], z[fixed_idx]; markersize=2, markercolor=:red, markerstrokewidth=0)
    display(plt)
end

# Plot the loaded system (optional)
plot3d_saddle(sys.points, sys.segments, title=model_name)

