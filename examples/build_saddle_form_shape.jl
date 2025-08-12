"""
Build Saddle Form Shape

Creates the saddle form geometry with diamond/rhombic mesh and plots the initial shape.
No simulation - just geometry creation and visualization.
"""

using SymbolicAWEModels
using LinearAlgebra
using Plots
using KiteUtils
using YAML
using PlotlyJS
plotlyjs()  # Use PlotlyJS backend for interactive 3D plots
# gr()
# ------------------------ Settings ------------------------

function build_settings()
    # Use the data/system.yaml file to avoid permission issues
    data_dir = joinpath(dirname(@__DIR__), "data")
    system_yaml = joinpath(data_dir, "system.yaml")
    set = Settings(system_yaml)
    set.v_wind     = 0.0
    set.g_earth    = 0.0
    set.e_tether   = 1.0e9
    set.rho_tether = 1.0
    return set
end

# ------------------- Transform Helpers -------------------

"""
Create a neutral transform (pure translation; no rotation applied).
- Moves `base_point_idx` to `base_pos`
- Keeps the current orientation (no net rotation) by matching target angles to current angles
"""
function neutral_transform(idx::Int; points::Vector{Point}, base_point_idx::Int=1,
                           rot_point_idx::Int=2, base_pos::Vector{Float64}=[0.0,0.0,0.0])
    # direction in CAD frame used by reinit! after translation
    dir0 = points[rot_point_idx].pos_cad .- points[base_point_idx].pos_cad
    # guard: if degenerate, bump the rot point
    if norm(dir0) < 1e-12
        error("neutral_transform: base_point_idx and rot_point_idx coincide; pick different indices")
    end
    curr_elev   = KiteUtils.calc_elevation(dir0)
    curr_azim   = -KiteUtils.azimuth_east(dir0)
    # target = current  => no net rotation
    return Transform(idx, curr_elev, curr_azim, 0.0;
                     base_pos=base_pos, base_point_idx=base_point_idx,
                     rot_point_idx=rot_point_idx)
end

# ------------------- YAML Data Loading -------------------

"""
    load_saddle_yaml(path::String)

Read a saddle dataset from YAML.
Returns (positions, connections, fixed_indices, params, masses)
"""
function load_saddle_yaml(path::String)
    data = YAML.load_file(path)

    haskey(data, "nodes") || error("YAML missing 'nodes'")
    haskey(data, "connections") || error("YAML missing 'connections'")

    one_based = get(data, "one_based", true)

    # positions & fixed & masses
    nodes = data["nodes"]
    N = length(nodes)
    positions = Vector{Vector{Float64}}(undef, N)
    fixed = Int[]
    masses = fill(0.1, N)  # default dynamic mass

    for (i, n) in enumerate(nodes)
        x = Float64(n["x"]); y = Float64(n["y"]); z = Float64(n["z"])
        positions[i] = [x,y,z]
        if get(n, "fixed", false)
            push!(fixed, i)
        end
        if haskey(n, "mass"); masses[i] = Float64(n["mass"]); end
    end

    # connections
    conns = Tuple{Int,Int}[]
    for c in data["connections"]
        i = Int(c["i"]); j = Int(c["j"])
        if !one_based
            i += 1; j += 1
        end
        push!(conns, (i,j))
    end

    # params (optional)
    params = Dict{String,Any}()
    if haskey(data, "params")
        for (k,v) in data["params"]; params[string(k)] = v; end
    end

    return positions, conns, fixed, params, masses
end

"""
    build_saddle_from_yaml(yaml_file::String, set; prestretch_frac=0.98, diameter_mm=3.0)

Build saddle form points and segments from YAML file.
Returns (points, segments, fixed_nodes::Vector{Int}, positions).
"""
function build_saddle_from_yaml(yaml_file::String, set; prestretch_frac=0.98, diameter_mm=3.0)
    !(0.9 ≤ prestretch_frac < 1.0) && error("prestretch_frac ∈ [0.9,1.0)")
    positions, connections, fixed_idx, params, masses = load_saddle_yaml(yaml_file)
    haskey(params, "n") && (length(positions) == params["n"] || error("nodes ≠ params.n"))

    fixed = Set(fixed_idx)
    points = Point[]
    for i in eachindex(positions)
        pos = positions[i]
        if i ∈ fixed
            push!(points, Point(i, pos, STATIC; transform_idx=1))
        else
            push!(points, Point(i, pos, DYNAMIC; transform_idx=1, mass=masses[i]))
        end
    end

    segments = Segment[]
    for (sid, (i,j)) in enumerate(connections)
        ℓ  = norm(positions[i] .- positions[j])
        l0 = prestretch_frac * ℓ
        push!(segments, Segment(sid, set, (i,j), BRIDLE; l0, compression_frac=0.001, diameter_mm=diameter_mm))
    end
    return points, segments, fixed_idx, positions
end

# ------------------------- 3D Plotting -------------------------

function plot3d_saddle(points, segments; title="Saddle Form")
    x = [p.pos_w[1] for p in points]; y = [p.pos_w[2] for p in points]; z = [p.pos_w[3] for p in points]
    p = scatter3d(x, y, z; markersize=3, markerstrokewidth=0, title, xlabel="X (m)", ylabel="Y (m)", zlabel="Z (m)", legend=false)
    for s in segments
        i,j = s.point_idxs
        plot3d!([x[i],x[j]],[y[i],y[j]],[z[i],z[j]]; alpha=0.6, linewidth=1)
    end
    fixed_idx = [i for (i,p) in enumerate(points) if p.type == STATIC]
    !isempty(fixed_idx) && scatter3d!(x[fixed_idx], y[fixed_idx], z[fixed_idx]; markersize=5, markercolor=:red, markerstrokewidth=0)
    return p
end

# -------------------------- Main --------------------------

"""
    build_and_plot_saddle_shape(; yaml_file="saddle.yaml", 
                                prestretch_frac=0.98, diameter_mm=3.0)

Creates the saddle form geometry from YAML file and plots the initial shape.
Returns the system structure for further analysis if needed.
"""
function build_and_plot_saddle_shape(; yaml_file="saddle.yaml", prestretch_frac=0.98, diameter_mm=3.0)
    println("Building saddle form from YAML...")
    set = build_settings()
    path = joinpath(@__DIR__, yaml_file); isfile(path) || error("Missing $path")

    points, segments, fixed_nodes, positions = build_saddle_from_yaml(path, set; prestretch_frac, diameter_mm)
    println("Created $(length(points)) points and $(length(segments)) segments")
    println("Fixed nodes: $(length(fixed_nodes))")

    # Preserve absolute positions (no rotation, no translation)
    transforms = [neutral_transform(1; points, base_point_idx=1, rot_point_idx=2, base_pos=Vector(points[1].pos_cad))]

    sys = SymbolicAWEModels.SystemStructure("saddle_form", set; points, segments, transforms)
    
    # Plot directly from the system structure points (avoid init! solver issues)
    println("Creating 3D plot...")
    fig = plot3d_saddle(sys.points, sys.segments; title="Saddle Form ($(length(points)) nodes)")
    display(fig)
    println("✓ Saddle form visualization complete!")
    return sys
end

# # Run when file is executed directly or included
# if abspath(PROGRAM_FILE) == @__FILE__
#     build_and_plot_saddle_shape(yaml_file="saddle.yaml")
# end

# When included from menu.jl, run automatically
build_and_plot_saddle_shape(yaml_file="saddle.yaml")
