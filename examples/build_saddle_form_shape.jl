"""
Build Saddle Form Shape

Creates the saddle form geometry with diamond/rhombic mesh and plots the initial shape.
No simulation - just geometry creation and visualization.
"""

using SymbolicAWEModels
using LinearAlgebra
using Plots
plotlyjs()  # Use PlotlyJS backend for interactive 3D plots

# ------------------------ Settings ------------------------

function build_settings()
    set = load_settings("base")
    set.v_wind     = 0.0
    set.g_earth    = 0.0
    set.e_tether   = 1.0e9
    set.rho_tether = 1.0
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

# ------------------------- 3D Plotting -------------------------

function plot3d_saddle_shape(points, segments, fixed_nodes; title="Saddle Form Shape")
    # Extract coordinates
    x = [p.pos_w[1] for p in points]
    y = [p.pos_w[2] for p in points]
    z = [p.pos_w[3] for p in points]
    
    # Create base plot with all points
    p = scatter3d(x, y, z, 
                  markersize=3, 
                  markercolor=:blue,
                  markerstrokewidth=0,
                  title=title,
                  xlabel="X (m)", 
                  ylabel="Y (m)", 
                  zlabel="Z (m)",
                  legend=false)
    
    # Add segments as lines
    for seg in segments
        i, j = seg.point_idxs
        plot3d!([x[i], x[j]], [y[i], y[j]], [z[i], z[j]], 
                color=:gray, alpha=0.6, linewidth=1)
    end
    
    # Highlight fixed nodes in red
    if !isempty(fixed_nodes)
        x_fixed = x[fixed_nodes]
        y_fixed = y[fixed_nodes]
        z_fixed = z[fixed_nodes]
        scatter3d!(x_fixed, y_fixed, z_fixed, 
                   markersize=5, 
                   markercolor=:red,
                   markerstrokewidth=0)
    end
    
    return p
end

# -------------------------- Main --------------------------

"""
    build_and_plot_saddle_shape(; grid_size=5, grid_length=10.0, grid_height=5.0,
                                prestretch_frac=0.98, diameter_mm=3.0)

Creates the saddle form geometry and plots the initial shape.
Returns the system structure for further analysis if needed.
"""
function build_and_plot_saddle_shape(; grid_size=5, grid_length=10.0, grid_height=5.0,
                                     prestretch_frac=0.98, diameter_mm=3.0, use_3d=true)

    println("SADDLE FORM SHAPE - SymbolicAWEModels")
    set = build_settings()

    points, segments, fixed_nodes, dx = build_saddle_points_segments(
        grid_size, grid_length, grid_height, set; prestretch_frac, diameter_mm)

    # Create a simple transform (no rotation, fixed in place)
    transforms = [Transform(1, 0.0, 0.0, 0.0; base_pos=[0.0, 0.0, 0.0], base_point_idx=1, rot_point_idx=1)]

    sys = SymbolicAWEModels.SystemStructure("saddle_form", set; points, segments, transforms)

    println("Saddle Form Shape (N=$grid_size, prestretch=$(prestretch_frac), dx=$(round(dx,digits=3)) m)")
    
    # Standard SymbolicAWEModels plot
    plot(sys, 0.0; zoom=false)
    
    # Optional 3D plot
    if use_3d
        title_str = "Saddle Form (N=$grid_size, dx=$(round(dx,digits=3))m)"
        fig3d = plot3d_saddle_shape(points, segments, fixed_nodes; title=title_str)
        display(fig3d)
    end

    return sys
end

# Run when file is executed directly or included
if abspath(PROGRAM_FILE) == @__FILE__
    build_and_plot_saddle_shape()
end

# When included from menu.jl, run automatically
build_and_plot_saddle_shape()
