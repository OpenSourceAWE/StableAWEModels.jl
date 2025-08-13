# SPDX-FileCopyrightText: 2025 Bart van de Lint
# SPDX-License-Identifier: MPL-2.0

"""
Saddle Form Benchmark (compact)

- Diamond/rhombic mesh with alternating full/offset rows
- Boundary nodes fixed with a saddle z-profile
- Self-stress via prestretch (l0 = prestretch_frac * current_length)
- Damped dynamic relaxation to equilibrium
"""

using SymbolicAWEModels, ControlPlots
using LinearAlgebra
using KiteUtils
using YAML

# ------------------------ Settings ------------------------

function build_settings(; sample_freq, abs_tol=1e-6, rel_tol=1e-6,
                       damping=nothing, e_tether=1.0e9, rho_tether=1.0, g_earth=0.0)
    # Use the data/system.yaml file to avoid permission issues (like build_saddle_form_shape.jl)
    data_dir = joinpath(dirname(@__DIR__), "data")
    system_yaml = joinpath(data_dir, "system.yaml")
    set = Settings(system_yaml)
    set.v_wind     = 0.0
    set.g_earth    = g_earth
    set.sample_freq = sample_freq
    set.abs_tol     = abs_tol
    set.rel_tol     = rel_tol
    set.e_tether    = e_tether
    set.rho_tether  = rho_tether
    set.l_tether = 10  # Default length for tether segments in meters   
    if damping !== nothing
        set.damping = damping
        println("Using manual damping: ", damping)
    else
        println("Using base damping: ", set.damping)
    end
    g_earth ≠ 0 && println("Warning: gravity ≠ 0; shape not purely prestress-driven.")
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
function build_saddle_from_yaml(yaml_file::String, set; prestretch_frac=0.98, diameter_mm=3.0,world_frame_damping=0.1)
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
            push!(points, Point(i, pos, DYNAMIC; transform_idx=1, mass=masses[i], world_frame_damping=world_frame_damping))
        end
    end

    segments = Segment[]
    for (sid, (i,j)) in enumerate(connections)
        ℓ  = norm(positions[i] .- positions[j])
        l0 = prestretch_frac * ℓ
        push!(segments, Segment(sid, set, (i,j), BRIDLE; l0, compression_frac=0.001, diameter_mm=diameter_mm))
    end
    
    println("YAML Mesh: nodes=$(length(points)), segments=$(length(segments)), fixed=$(length(fixed_idx))")
    return points, segments, fixed_idx, positions
end

# ------------------------ Solver --------------------------

"""
    relax!(sam; max_steps=10_000, vtol=1e-4, report_every=1000)

Velocity-based convergence (robust across SymbolicAWEModels versions).
Returns (steps, runtime, max_speed).
"""
function relax!(sam; max_steps=100, vtol=1e-4, report_every=1000)
    t0 = time(); maxv = Inf; steps = 0
    for k in 1:max_steps
        next_step!(sam)
        steps = k
        maxv = 0.0
        @inbounds for p in sam.sys_struct.points
            p.type === DYNAMIC || continue
            v = norm(p.vel_w)
            v > maxv && (maxv = v)
        end
        (k % report_every == 0) && println(" step $k: max|v|=$(round(maxv,digits=6))")
        if maxv < vtol
            println("Converged at step $k with max|v|=$(round(maxv,digits=6))")
            break
        end
    end
    return steps, (time()-t0), maxv
end

# ------------------------- 3D Plotting -------------------------

function plot3d_saddle(points, segments; title="Saddle Form")
    x = [p.pos_w[1] for p in points]; y = [p.pos_w[2] for p in points]; z = [p.pos_w[3] for p in points]
    p = Plots.scatter3d(x, y, z; markersize=2, markerstrokewidth=0, title, xlabel="X (m)", ylabel="Y (m)", zlabel="Z (m)", legend=false)
    for s in segments
        i,j = s.point_idxs
        Plots.plot3d!([x[i],x[j]],[y[i],y[j]],[z[i],z[j]]; alpha=1, linewidth=1,color=:black)
    end
    fixed_idx = [i for (i,p) in enumerate(points) if p.type == STATIC]
    !isempty(fixed_idx) && Plots.scatter3d!(x[fixed_idx], y[fixed_idx], z[fixed_idx]; markersize=2, markercolor=:red, markerstrokewidth=0)
    return p
end

# -------------------------- Main --------------------------

"""
    main(; yaml_file="saddle.yaml", prestretch_frac=0.98, diameter_mm=3.0,
          sample_freq=100, abs_tol=1e-6, rel_tol=1e-6,
          damping=nothing, max_steps=10, vtol=1e-4)

Run saddle form simulation using YAML geometry.
Returns (sam, XYZ_initial, XYZ_final).
"""
function main(; yaml_file="saddle_gridsize5.yaml", prestretch_frac=0.98, diameter_mm=3.0,
               sample_freq=10, abs_tol=1e-6, rel_tol=1e-6,
               damping=nothing, max_steps=10, vtol=1e-4)

    println("SADDLE FORM (YAML-based) - SymbolicAWEModels")
    set = build_settings(; sample_freq, abs_tol, rel_tol, damping, e_tether=1e9, rho_tether=1.0, g_earth=0.0)

    # Load geometry from YAML file (same as build_saddle_form_shape.jl)
    yaml_path = joinpath(@__DIR__, yaml_file)
    isfile(yaml_path) || error("Missing YAML file: $yaml_path")
    
    points, segments, fixed_nodes, positions = build_saddle_from_yaml(
        yaml_path, set; prestretch_frac, diameter_mm)

    # Preserve absolute positions (no rotation, no translation) - neutral transform
    transforms = [neutral_transform(1; points, base_point_idx=1, rot_point_idx=2, base_pos=Vector(points[1].pos_cad))]

    sys = SymbolicAWEModels.SystemStructure("saddle_form", set; points, segments, transforms)
    sam = SymbolicAWEModel(set, sys)
    init!(sam; remake=false)

    XYZ0 = [copy(p.pos_w) for p in sys.points]
    fig0 = plot3d_saddle(sam.sys_struct.points, sam.sys_struct.segments; title="Initial State")
    display(fig0)


    # running the simulation
    n_steps = 200
    runtime = n_steps / set.sample_freq
    maxv = 0
    ControlPlots.plot(sam.sys_struct, 0.0; zoom=false)
    for i in 1:n_steps
        current_time = i/set.sample_freq
        ControlPlots.plot(sam, current_time; zoom=false)
        next_step!(sam)
    end


    # steps, runtime, maxv = relax!(sam; max_steps, vtol)
    XYZf = [copy(p.pos_w) for p in sam.sys_struct.points]


    figf = plot3d_saddle(sam.sys_struct.points, sam.sys_struct.segments; title="Final State")
    display(figf)


    # quick stats
    disp = [norm(XYZf[i] .- XYZ0[i]) for i in eachindex(XYZ0) if i ∉ Set(fixed_nodes)]
    println("Max disp=$(round(maximum(disp),digits=4)) m, Avg disp=$(round(sum(disp)/length(disp),digits=4)) m")

    return sam, XYZ0, XYZf
end

if abspath(PROGRAM_FILE) == @__FILE__
    main()
end

# When included from menu.jl, run main() automatically
main(yaml_file="saddle_gridsize5.yaml")
