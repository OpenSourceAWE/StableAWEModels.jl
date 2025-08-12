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

# ------------------------ Settings ------------------------

function build_settings(; sample_freq=100, abs_tol=1e-6, rel_tol=1e-6,
                       damping=nothing, e_tether=1.0e9, rho_tether=1.0, g_earth=0.0)
    set = load_settings("base")
    set.v_wind     = 0.0
    set.g_earth    = g_earth
    set.sample_freq = sample_freq
    set.abs_tol     = abs_tol
    set.rel_tol     = rel_tol
    set.e_tether    = e_tether
    set.rho_tether  = rho_tether
    if damping !== nothing
        set.damping = damping
        println("Using manual damping: ", damping)
    else
        println("Using base damping: ", set.damping)
    end
    g_earth ≠ 0 && println("Warning: gravity ≠ 0; shape not purely prestress-driven.")
    return set
end

# ------------------- Mesh + Connectivity -------------------

"""
    build_saddle_points_segments(grid_size, grid_length, grid_height, set;
                                 prestretch_frac=0.98, diameter_mm=3.0)

Returns (points, segments, fixed_nodes::Vector{Int}, dx).
"""
function build_saddle_points_segments(N::Int, L::Real, H::Real, set;
                                      prestretch_frac=0.98, diameter_mm=3.0)
    N < 3 && error("grid_size must be ≥ 3")
    !(0.9 ≤ prestretch_frac < 1.0) && error("prestretch_frac must be in [0.9, 1.0)")

    dx = L / (N - 1)
    n_total = N^2 + (N - 1)^2
    positions = Vector{Vector{Float64}}(undef, 0)
    node_map = Dict{Tuple{Int,Int}, Int}()       # (row, col) -> idx
    idx_to_rc = Vector{Tuple{Int,Int}}(undef, n_total)

    # Build rows: even rows (full N), odd rows (N-1) offset by dx/2; y advances by dx/2
    idx = 1
    for row in 0:(2N - 2)
        y = (row * dx) / 2
        if iseven(row)
            for col in 0:(N - 1)
                x = col * dx
                push!(positions, [x, y, H/2])
                node_map[(row, col)] = idx
                idx_to_rc[idx] = (row, col)
                idx += 1
            end
        else
            for col in 0:(N - 2)
                x = (col + 0.5) * dx
                push!(positions, [x, y, H/2])
                node_map[(row, col)] = idx
                idx_to_rc[idx] = (row, col)
                idx += 1
            end
        end
    end
    @assert length(positions) == n_total "Position count mismatch"

    # Fixed nodes (boundary): top row, bottom row, left/right edges (non-corner odd rows only)
    top    = [node_map[(0, c)]     for c in 0:(N-1)]
    bottom = [node_map[(2N-2, c)]  for c in 0:(N-1)]
    
    # Left/right edges: only middle odd rows (not adjacent to top/bottom)
    # For diamond mesh, odd rows that connect to interior
    middle_odd_rows = [r for r in 1:2:(2N-3) if r < 2N-3]  # exclude row 2N-3 which connects to bottom
    left   = [node_map[(r, 0)]     for r in middle_odd_rows]
    right  = [node_map[(r, N-2)]   for r in middle_odd_rows]
    
    fixed_nodes = vcat(top, bottom, left, right)
    expected_count = 2*N + 2*(N-2)  # top + bottom + left + right
    @assert length(fixed_nodes) == expected_count "Fixed-node count mismatch: got $(length(fixed_nodes)), expected $expected_count"
    fixed_set = Set(fixed_nodes)

    # Boundary z profile (saddle)
    dl = H / L * dx
    zline = [i*dl for i in 0:(N-1)]
    boundary_z = vcat(zline, reverse(zline), zline[2:end-1], reverse(zline[2:end-1]))
    @assert length(boundary_z) == length(fixed_nodes)

    for (k, idxn) in enumerate(fixed_nodes)
        positions[idxn][3] = boundary_z[k]
    end

    # Points (STATIC boundary, DYNAMIC interior)
    points = Point[]
    for i in 1:n_total
        pos = positions[i]
        if i ∈ fixed_set
            push!(points, Point(i, pos, STATIC; transform_idx=1))
        else
            push!(points, Point(i, pos, DYNAMIC; mass=0.1, transform_idx=1))
        end
    end

    # Connectivity (neighbors between adjacent parity rows), unique pairs
    pairs = Set{Tuple{Int,Int}}()
    for i in 1:n_total
        row, col = idx_to_rc[i]
        neigh = Int[]
        if iseven(row)
            # link to odd rows above/below
            if row > 0
                (col > 0)         && haskey(node_map, (row-1, col-1)) && push!(neigh, node_map[(row-1, col-1)])
                (col < N-1)       && haskey(node_map, (row-1, col))   && push!(neigh, node_map[(row-1, col)])
            end
            if row < 2N - 2
                (col > 0)         && haskey(node_map, (row+1, col-1)) && push!(neigh, node_map[(row+1, col-1)])
                (col < N-1)       && haskey(node_map, (row+1, col))   && push!(neigh, node_map[(row+1, col)])
            end
        else
            # link to even rows above/below
            haskey(node_map, (row-1, col))     && push!(neigh, node_map[(row-1, col)])
            haskey(node_map, (row-1, col+1))   && push!(neigh, node_map[(row-1, col+1)])
            (row < 2N - 2) && haskey(node_map, (row+1, col))   && push!(neigh, node_map[(row+1, col)])
            (row < 2N - 2) && haskey(node_map, (row+1, col+1)) && push!(neigh, node_map[(row+1, col+1)])
        end
        for j in neigh
            i == j && continue
            a,b = i<j ? (i,j) : (j,i)
            push!(pairs, (a,b))
        end
    end

    # Segments (prestretch)
    segments = Segment[]
    for (sid, (i,j)) in enumerate(sort!(collect(pairs)))
        ℓ = norm(positions[i] .- positions[j])
        l0 = prestretch_frac * ℓ
        push!(segments, Segment(sid, set, (i,j), BRIDLE; l0, compression_frac=0.001, diameter_mm=diameter_mm))
    end

    println("Mesh: N=$N, nodes=$(length(points)), segments=$(length(segments)), dx=$(round(dx,digits=3)) m")
    return points, segments, fixed_nodes, dx
end

# ------------------------ Solver --------------------------

"""
    relax!(sam; max_steps=10_000, vtol=1e-4, report_every=1000)

Velocity-based convergence (robust across SymbolicAWEModels versions).
Returns (steps, runtime, max_speed).
"""
function relax!(sam; max_steps=10_000, vtol=1e-4, report_every=1000)
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

# ------------------------- Plots --------------------------

function plot_state(title_str, XYZ, segments, fixed_nodes::Vector{Int})
    x = getindex.(XYZ, 1); y = getindex.(XYZ, 2); z = getindex.(XYZ, 3)
    fig = ControlPlots.plot3d(x, y, z; title=title_str, markersize=4, markerstrokewidth=0)
    for s in segments
        i,j = s.point_idxs
        ControlPlots.plot3d!([x[i],x[j]], [y[i],y[j]], [z[i],z[j]]; color=:gray, alpha=0.7, linewidth=1.0)
    end
    xf = x[fixed_nodes]; yf = y[fixed_nodes]; zf = z[fixed_nodes]
    ControlPlots.plot3d!(xf, yf, zf; color=:red, markersize=6, markerstrokewidth=0)
    display(fig)
    return fig
end

# -------------------------- Main --------------------------

"""
    main(; grid_size=5, grid_length=10.0, grid_height=5.0,
          prestretch_frac=0.98, diameter_mm=3.0,
          sample_freq=100, abs_tol=1e-6, rel_tol=1e-6,
          damping=nothing, max_steps=10_000, vtol=1e-4)

Returns (sam, XYZ_initial, XYZ_final).
"""
function main(; grid_size=5, grid_length=10.0, grid_height=5.0,
               prestretch_frac=0.98, diameter_mm=3.0,
               sample_freq=100, abs_tol=1e-6, rel_tol=1e-6,
               damping=nothing, max_steps=10_000, vtol=1e-4)

    println("SADDLE FORM - SymbolicAWEModels")
    set = build_settings(; sample_freq, abs_tol, rel_tol, damping, e_tether=1e9, rho_tether=1.0, g_earth=0.0)

    points, segments, fixed_nodes, dx = build_saddle_points_segments(
        grid_size, grid_length, grid_height, set; prestretch_frac, diameter_mm)

    # Create a simple transform (no rotation, fixed in place)
    transforms = [Transform(1, 0.0, 0.0, 0.0; base_pos=[0.0, 0.0, 0.0], base_point_idx=1, rot_point_idx=1)]

    sys = SymbolicAWEModels.SystemStructure("saddle_form", set; points, segments, transforms)
    sam = SymbolicAWEModel(set, sys)
    init!(sam; remake=false)

    XYZ0 = [copy(p.pos_w) for p in sys.points]
    fig0 = plot_state("Initial (N=$grid_size, prestretch=$(prestretch_frac), dx=$(round(dx,digits=3)) m)",
                      XYZ0, segments, fixed_nodes)

    steps, runtime, maxv = relax!(sam; max_steps, vtol)

    XYZf = [copy(p.pos_w) for p in sam.sys_struct.points]
    figf = plot_state("Final (steps=$steps, runtime=$(round(runtime,digits=3)) s, max|v|=$(round(maxv,digits=6)))",
                      XYZf, segments, fixed_nodes)

    # quick stats
    disp = [norm(XYZf[i] .- XYZ0[i]) for i in eachindex(XYZ0) if i ∉ Set(fixed_nodes)]
    println("Max disp=$(round(maximum(disp),digits=4)) m, Avg disp=$(round(sum(disp)/length(disp),digits=4)) m")

    return sam, XYZ0, XYZf
end

if abspath(PROGRAM_FILE) == @__FILE__
    main()
end
