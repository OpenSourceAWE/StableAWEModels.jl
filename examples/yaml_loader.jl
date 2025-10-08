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

# Helper struct to collect raw points before reindexing
struct RawPoint
    raw_id::Int
    pos::Vector{Float64}
    type::SymbolicAWEModels.DynamicsType
    mass::Float64
    body_damping::Float64
    world_damping::Float64
end

"""
        load_sys_struct_from_yaml(yaml_path::AbstractString; system_name="from_yaml", set=nothing)

Build a `SymbolicAWEModels.SystemStructure` from a TU Delft style structural YAML file.

# Expected top-level blocks
- `bridle_point_node`: `[x,y,z]` location of the KCU/bridle origin (optional).
- `wing_particles`: table with headers `[id,x,y,z]` defining wing node positions.
- `wing_connections`: table with headers `[name,ci,cj]` describing structural members.
- `wing_elements`: table with headers `[name,l0,k,c,m,linktype]` supplying segment properties.
- `bridle_particles`: table with headers `[id,x,y,z]` defining bridle node positions.
- `bridle_connections`: table with headers `[name,ci,cj,ck]`; `ck` is currently logged but unused.
- `bridle_elements`: table with headers `[name,l0,d,material,linktype]` supplying bridle segment data.
- Material property blocks (e.g. `dyneema`) can be specified with `youngs_modulus`, `density`, and `damping_per_stiffness`.

Indices are assumed to be the raw IDs given in the YAML (0-based for the optional KCU node). They are reindexed internally to satisfy `SystemStructure` requirements.
"""
function load_sys_struct_from_yaml(yaml_path::AbstractString; system_name="from_yaml", set=nothing())
    data = YAML.load_file(yaml_path)

    # Use provided settings or fall back to base settings
    set === nothing && (set = SymbolicAWEModels.load_settings("base"))

    # Collect material properties if present (e.g., dyneema block)
    material_props = Dict{String,Dict{String,Float64}}()
    for (name, block) in data
        if isa(block, Dict) && haskey(block, "youngs_modulus")
            props = Dict{String,Float64}()
            for (k, v) in block
                props[String(k)] = Float64(v)
            end
            material_props[String(name)] = props
        end
    end

    # Identify fixed point ids (stored as raw YAML ids, 0-based as documented)
    fixed_ids = Set{Int}()
    if haskey(data, "fixed_point_indices")
        for idx in data["fixed_point_indices"]
            push!(fixed_ids, Int(idx))
        end
    end

    raw_points = RawPoint[]

    # Optional bridle/KCU node (id 0 by convention)
    if haskey(data, "bridle_point_node")
        pos_vec = Float64.(data["bridle_point_node"])
        raw_id = 0
        ptype = raw_id in fixed_ids ? SymbolicAWEModels.STATIC : SymbolicAWEModels.DYNAMIC
        mass = get(data, "kcu_mass", 0.0)
        push!(raw_points, RawPoint(raw_id, pos_vec, ptype, Float64(mass), 0.0, 0.0))
    end

    # Wing particles
    if haskey(data, "wing_particles")
        headers = Symbol.(String.(data["wing_particles"]["headers"]))
        for row in data["wing_particles"]["data"]
            section = Dict(zip(headers, row))
            raw_id = Int(section[:id])
            pos_vec = [Float64(section[:x]), Float64(section[:y]), Float64(section[:z])]
            ptype = raw_id in fixed_ids ? SymbolicAWEModels.STATIC : SymbolicAWEModels.DYNAMIC
            push!(raw_points, RawPoint(raw_id, pos_vec, ptype, 0.0, 0.0, 0.0))
        end
    end

    # Bridle particles
    if haskey(data, "bridle_particles")
        headers = Symbol.(String.(data["bridle_particles"]["headers"]))
        for row in data["bridle_particles"]["data"]
            section = Dict(zip(headers, row))
            raw_id = Int(section[:id])
            pos_vec = [Float64(section[:x]), Float64(section[:y]), Float64(section[:z])]
            ptype = raw_id in fixed_ids ? SymbolicAWEModels.STATIC : SymbolicAWEModels.DYNAMIC
            push!(raw_points, RawPoint(raw_id, pos_vec, ptype, 0.0, 0.0, 0.0))
        end
    end

    # Sort raw points to ensure deterministic ordering before reindexing
    sort!(raw_points, by = rp -> rp.raw_id)

    # Reindex points sequentially while keeping a mapping from raw IDs to internal indices
    points = SymbolicAWEModels.Point[]
    raw_to_idx = Dict{Int,Int}()
    for rp in raw_points
        idx = length(points) + 1
        push!(points, SymbolicAWEModels.Point(
            idx,
            rp.pos,
            rp.type;
            mass = rp.mass,
            body_frame_damping = rp.body_damping,
            world_frame_damping = rp.world_damping
        ))
        raw_to_idx[rp.raw_id] = idx
    end

    isempty(points) && error("No points found in YAML file $(yaml_path).")
    length(points) < 2 && error("At least two points are required to define a transform in $(yaml_path).")

    # --- Segment helpers -----------------------------------------------------
    wing_element_dict = Dict{String,Dict{Symbol,Any}}()
    if haskey(data, "wing_elements")
        headers = Symbol.(String.(data["wing_elements"]["headers"]))
        for row in data["wing_elements"]["data"]
            section = Dict(zip(headers, row))
            wing_element_dict[String(section[:name])] = section
        end
    end

    bridle_element_dict = Dict{String,Dict{Symbol,Any}}()
    if haskey(data, "bridle_elements")
        headers = Symbol.(String.(data["bridle_elements"]["headers"]))
        for row in data["bridle_elements"]["data"]
            section = Dict(zip(headers, row))
            bridle_element_dict[String(section[:name])] = section
        end
    end

    # --- Segments ------------------------------------------------------------
    segments = SymbolicAWEModels.Segment[]

    # Wing structural connections
    if haskey(data, "wing_connections")
        headers = Symbol.(String.(data["wing_connections"]["headers"]))
        for row in data["wing_connections"]["data"]
            section = Dict(zip(headers, row))
            name = String(section[:name])
            ci = Int(section[:ci])
            cj = Int(section[:cj])
            haskey(raw_to_idx, ci) || (@warn "wing segment $name references missing point id $ci"; continue)
            haskey(raw_to_idx, cj) || (@warn "wing segment $name references missing point id $cj"; continue)
            elem = get(wing_element_dict, name, Dict{Symbol,Any}())
            l0 = haskey(elem, :l0) ? Float64(elem[:l0]) : norm(raw_points[findfirst(rp -> rp.raw_id == ci, raw_points)].pos .- raw_points[findfirst(rp -> rp.raw_id == cj, raw_points)].pos)
            axial_stiffness = haskey(elem, :k) ? Float64(elem[:k]) : NaN
            axial_damping   = haskey(elem, :c) ? Float64(elem[:c]) : NaN
            diameter_mm = DEFAULT_WING_DIAMETER_MM
            if haskey(elem, :d)
                diameter_mm = Float64(elem[:d]) * 1000
            end
            push!(segments, SymbolicAWEModels.Segment(
                length(segments) + 1,
                set,
                (raw_to_idx[ci], raw_to_idx[cj]),
                SymbolicAWEModels.POWER_LINE;
                l0 = l0,
                diameter_mm = diameter_mm,
                axial_stiffness = axial_stiffness,
                axial_damping = axial_damping,
                compression_frac = 0.01
            ))
        end
    end

    # Bridle connections and pulleys
    pulleys = SymbolicAWEModels.Pulley[]
    if haskey(data, "bridle_connections")
        headers = Symbol.(String.(data["bridle_connections"]["headers"]))
        for (row_idx, row) in enumerate(data["bridle_connections"]["data"])
            section = Dict(zip(headers, row))
            name = String(section[:name])
            ci = Int(section[:ci])
            cj = Int(section[:cj])
            haskey(raw_to_idx, ci) || (@warn "bridle segment $name references missing point id $ci"; continue)
            haskey(raw_to_idx, cj) || (@warn "bridle segment $name references missing point id $cj"; continue)
            
            # Check if this is a pulley (has ck column)
            ck = haskey(section, :ck) ? Int(section[:ck]) : nothing
            is_pulley = ck !== nothing
            
            if is_pulley
                # This is a pulley connection: ci and cj are anchor points, ck is the pulley node
                haskey(raw_to_idx, ck) || (@warn "pulley $name references missing point id $ck"; continue)
                
                # Create pulley with anchor points (ci, cj)
                push!(pulleys, SymbolicAWEModels.Pulley(
                    length(pulleys) + 1,
                    (raw_to_idx[ci], raw_to_idx[cj]),
                    SymbolicAWEModels.DYNAMIC
                ))
                
                # Get element properties
                elem = get(bridle_element_dict, name, nothing)
                if elem === nothing
                    @warn "No bridle element properties found for pulley $name. Skipping."
                    continue
                end
                
                # Get total rest length from element properties
                total_l0 = Float64(elem[:l0])
                diameter_mm = haskey(elem, :d) ? Float64(elem[:d]) * 1000 : NaN
                diameter_m = haskey(elem, :d) ? Float64(elem[:d]) : NaN
                
                # Find raw points to calculate geometric distances
                idx_ci = findfirst(rp -> rp.raw_id == ci, raw_points)
                idx_cj = findfirst(rp -> rp.raw_id == cj, raw_points)
                idx_ck = findfirst(rp -> rp.raw_id == ck, raw_points)
                
                if idx_ci === nothing || idx_cj === nothing || idx_ck === nothing
                    @warn "Could not find raw points for pulley $name"
                    continue
                end
                
                pos_ci = raw_points[idx_ci].pos
                pos_cj = raw_points[idx_cj].pos
                pos_ck = raw_points[idx_ck].pos
                
                # Calculate straight-line distances (like Python version)
                len_ci_cj = norm(pos_ci - pos_cj)
                len_cj_ck = norm(pos_cj - pos_ck)
                len_total = len_ci_cj + len_cj_ck
                
                # Proportionally divide rest length based on geometric distances
                l0_ci_cj = (len_ci_cj / len_total) * total_l0
                l0_cj_ck = (len_cj_ck / len_total) * total_l0
                
                # Get material properties
                mat_name = haskey(elem, :material) ? String(elem[:material]) : ""
                mat_props = get(material_props, mat_name, nothing)
                
                # Calculate stiffness and damping based on total length (like Python)
                E = mat_props !== nothing ? get(mat_props, "youngs_modulus", NaN) : NaN
                damping_frac = mat_props !== nothing ? get(mat_props, "damping_per_stiffness", 0.0) : 0.0
                
                if !isnan(E) && !isnan(diameter_m) && total_l0 > 0
                    area = π * (diameter_m / 2)^2
                    # Use total_l0 for stiffness calculation (Python approach)
                    k_total = E * area / total_l0
                    c_total = damping_frac * k_total
                else
                    k_total = NaN
                    c_total = NaN
                end
                
                # Add pulley mass to the pulley node
                pulley_mass = get(data, "pulley_mass", 0.0)
                pulley_idx = findfirst(rp -> rp.raw_id == cj, raw_points)
                if pulley_idx !== nothing
                    # Note: In the current structure, mass is set during Point creation
                    # This would need to be handled differently if mass needs to be updated
                end
                
                compression_frac = 0.0
                if haskey(elem, :linktype)
                    lt = lowercase(String(elem[:linktype]))
                    compression_frac = lt == "noncompressive" ? 0.0 : 0.01
                end
                
                # Create segment from ci to cj (first half of pulley)
                push!(segments, SymbolicAWEModels.Segment(
                    length(segments) + 1,
                    set,
                    (raw_to_idx[ci], raw_to_idx[cj]),
                    SymbolicAWEModels.BRIDLE;
                    l0 = l0_ci_cj,
                    diameter_mm = diameter_mm,
                    axial_stiffness = k_total,
                    axial_damping = c_total,
                    compression_frac = compression_frac
                ))
                
                # Create segment from cj to ck (second half of pulley)
                push!(segments, SymbolicAWEModels.Segment(
                    length(segments) + 1,
                    set,
                    (raw_to_idx[cj], raw_to_idx[ck]),
                    SymbolicAWEModels.BRIDLE;
                    l0 = l0_cj_ck,
                    diameter_mm = diameter_mm,
                    axial_stiffness = k_total,
                    axial_damping = c_total,
                    compression_frac = compression_frac
                ))
            else
                # Regular bridle segment (no pulley)
                elem = get(bridle_element_dict, name, nothing)
                if elem === nothing
                    @warn "No bridle element properties found for segment $name (row $row_idx). Using defaults."
                end
                diameter_mm = elem !== nothing && haskey(elem, :d) ? Float64(elem[:d]) * 1000 : NaN
                l0 = elem !== nothing && haskey(elem, :l0) ? Float64(elem[:l0]) : 0.0
                axial_stiffness = NaN
                axial_damping = NaN
                if elem !== nothing
                    mat_name = haskey(elem, :material) ? String(elem[:material]) : ""
                    mat_props = get(material_props, mat_name, nothing)
                    if mat_props !== nothing
                        E = get(mat_props, "youngs_modulus", NaN)
                        damping_frac = get(mat_props, "damping_per_stiffness", 0.0)
                        diameter_m = haskey(elem, :d) ? Float64(elem[:d]) : NaN
                        if !isnan(E) && !isnan(diameter_m) && l0 > 0
                            area = π * (diameter_m / 2)^2
                            axial_stiffness = E * area / l0
                            axial_damping = damping_frac * axial_stiffness
                        end
                    end
                end
                compression_frac = 0.0
                if elem !== nothing && haskey(elem, :linktype)
                    lt = lowercase(String(elem[:linktype]))
                    compression_frac = lt == "noncompressive" ? 0.0 : 0.01
                end
                push!(segments, SymbolicAWEModels.Segment(
                    length(segments) + 1,
                    set,
                    (raw_to_idx[ci], raw_to_idx[cj]),
                    SymbolicAWEModels.BRIDLE;
                    l0 = l0,
                    diameter_mm = diameter_mm,
                    axial_stiffness = axial_stiffness,
                    axial_damping = axial_damping,
                    compression_frac = compression_frac
                ))
            end
        end
    end

    # At least one transform is needed; keep neutral to preserve CAD orientation
    transforms = [neutral_transform(1; points=points, base_point_idx=1, rot_point_idx=2, base_pos=Vector{Float64}(points[1].pos_cad))]

    return SymbolicAWEModels.SystemStructure(system_name, set; points, segments, pulleys, transforms)
end


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

"""
    plot3d_v3(points, segments; title::AbstractString="3D Structure")

Create an interactive 3D plot using PlotlyJS showing the structure with points and segments.
Fixed points are shown in red, dynamic points in blue, and segments as black lines.

# Arguments
- `points`: Vector of Point objects with position information
- `segments`: Vector of Segment objects defining connections between points
- `title`: Plot title (default: "3D Structure")
"""
function plot3d_v3(points, segments; title::AbstractString="3D Structure")
    # Use pos_cad if pos_w is not available
    get_pos(p) = hasproperty(p, :pos_w) ? p.pos_w : p.pos_cad
    
    # Extract point coordinates
    x = [get_pos(p)[1] for p in points]
    y = [get_pos(p)[2] for p in points]
    z = [get_pos(p)[3] for p in points]
    
    # Separate fixed and dynamic points
    fixed_idx = [i for (i, pt) in enumerate(points) if pt.type == SymbolicAWEModels.STATIC]
    dynamic_idx = [i for (i, pt) in enumerate(points) if pt.type != SymbolicAWEModels.STATIC]
    
    # Create traces
    traces = PlotlyJS.GenericTrace[]
    
    # Add segments as lines
    for s in segments
        i, j = s.point_idxs
        if i == 0 || j == 0
            @warn "Segment with zero index detected: $(s)"
            continue
        end
        segment_trace = PlotlyJS.scatter3d(
            x=[x[i], x[j]],
            y=[y[i], y[j]],
            z=[z[i], z[j]],
            mode="lines",
            line=PlotlyJS.attr(color="black", width=2),
            showlegend=false,
            hoverinfo="skip"
        )
        push!(traces, segment_trace)
    end
    
    # Add dynamic points
    if !isempty(dynamic_idx)
        dynamic_trace = PlotlyJS.scatter3d(
            x=x[dynamic_idx],
            y=y[dynamic_idx],
            z=z[dynamic_idx],
            mode="markers",
            marker=PlotlyJS.attr(
                size=4,
                color="blue",
                symbol="circle"
            ),
            name="Dynamic Points",
            hovertemplate="Point %{text}<br>x: %{x:.3f}<br>y: %{y:.3f}<br>z: %{z:.3f}<extra></extra>",
            text=string.(dynamic_idx)
        )
        push!(traces, dynamic_trace)
    end
    
    # Add fixed points
    if !isempty(fixed_idx)
        fixed_trace = PlotlyJS.scatter3d(
            x=x[fixed_idx],
            y=y[fixed_idx],
            z=z[fixed_idx],
            mode="markers",
            marker=PlotlyJS.attr(
                size=6,
                color="red",
                symbol="diamond"
            ),
            name="Fixed Points",
            hovertemplate="Point %{text} (FIXED)<br>x: %{x:.3f}<br>y: %{y:.3f}<br>z: %{z:.3f}<extra></extra>",
            text=string.(fixed_idx)
        )
        push!(traces, fixed_trace)
    end
    
    # Create layout
    layout = PlotlyJS.Layout(
        title=title,
        scene=PlotlyJS.attr(
            xaxis=PlotlyJS.attr(title="X (m)"),
            yaxis=PlotlyJS.attr(title="Y (m)"),
            zaxis=PlotlyJS.attr(title="Z (m)"),
            aspectmode="data"
        ),
        showlegend=true,
        hovermode="closest"
    )
    
    # Create and display plot
    plt = PlotlyJS.plot(traces, layout)
    display(plt)
    return plt
end


# ------------------- Example usage -------------------

### Change to your model name
# # model_name = "pyramid_model"
# model_name = "TUDELFT_V3_KITE" 

# println("\n\nYAML GEOMETRY LOADER\n", "="^40)
# set = SymbolicAWEModels.load_settings(model_name)  # Loads as Dict
# geometry_path = joinpath(dirname(@__DIR__), "data", model_name, "struc_geometry.yaml")
# sys = load_sys_struct_from_yaml(geometry_path; system_name=model_name, set=set)

# Plot the loaded system (optional)
# Use plot3d_saddle for Plots.jl backend:
# plot3d_saddle(sys.points, sys.segments, title=model_name)

# Use plot3d_v3 for interactive PlotlyJS visualization:
# plot3d_v3(sys.points, sys.segments, title=model_name)
