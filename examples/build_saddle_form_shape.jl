"""
Build Saddle Form Shape

Creates the saddle form geometry with diamond/rhombic mesh and plots the initial shape.
No simulation - just geometry creation and visualization.
"""

using SymbolicAWEModels
using LinearAlgebra
using Plots
using KiteUtils
# using PlotlyJS
# plotlyjs()  # Use PlotlyJS backend for interactive 3D plots
gr()
# ------------------------ Settings ------------------------

function build_settings()
    set = SymbolicAWEModels.load_settings("base")
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

    # Fixed nodes (boundary): match Python script exactly
    # top_edge: first row (size N)
    # bottom_edge: last row (size N)  
    # left_edge: middle odd rows only, excluding ones adjacent to top & bottom → exactly N-2 nodes
    # right_edge: same pattern → N-2 nodes
    # Total fixed nodes = 2N + 2(N-2) = 4N - 4
    
    top    = [node_map[(0, c)]     for c in 0:(N-1)]
    bottom = [node_map[(2N-2, c)]  for c in 0:(N-1)]
    
    # Left/right edges: middle odd rows only (exclude the last odd row adjacent to bottom)
    # For N ≥ 3: odd rows are 1, 3, 5, ..., (2N-3)
    # middle odd rows should exclude the one just above the bottom (2N-3)
    left_rows = (N > 3) ? collect(1:2:(2N-5)) : (N == 3 ? [1] : Int[])
    left   = [node_map[(r, 0)]     for r in left_rows]
    right  = [node_map[(r, N-2)]   for r in left_rows]
    
    fixed_nodes = vcat(top, bottom, left, right)
    # Exact Python count: 4N - 4
    expected_count = 4*N - 4  
    println("Debug: N=$N, left_rows=$left_rows, expected_count=$expected_count, actual=$(length(fixed_nodes))")
    @assert length(fixed_nodes) == expected_count "Fixed-node count mismatch: got $(length(fixed_nodes)), expected $expected_count"
    fixed_set = Set(fixed_nodes)

    # Boundary z profile (saddle) - identical to Python
    dl = H / L * dx
    zline = [i*dl for i in 0:(N-1)]
    
    # Build z array in order [top, bottom, left, right] with lengths [N, N, N-2, N-2]
    boundary_z = vcat(zline, reverse(zline), zline[2:end-1], reverse(zline[2:end-1]))
    println("Debug: length(boundary_z)=$(length(boundary_z)), length(fixed_nodes)=$(length(fixed_nodes))")
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
    return points, segments, fixed_nodes, dx, positions  # Also return raw positions
end

# ------------------------- 3D Plotting -------------------------

"""
    check_edges(points, top, bottom, left, right, dx, L)

Debug function to validate that boundary nodes are correctly positioned.
"""
function check_edges(points, top, bottom, left, right, dx, L)
    x = [p.pos_w[1] for p in points]
    y = [p.pos_w[2] for p in points]
    z = [p.pos_w[3] for p in points]

    println("\nEdge Validation:")
    println("Top y ≈ ", unique(round.(y[top], digits=6)))
    println("Bottom y ≈ ", unique(round.(y[bottom], digits=6)))
    println("Left x (odd rows) ≈ ", unique(round.(x[left], digits=6)))
    println("Right x (odd rows) ≈ ", unique(round.(x[right], digits=6)))

    # In this diamond layout, left/right x should be around dx/2 and L - dx/2
    println("Expected: left x ≈ ", round(dx/2, digits=6), 
            ", right x ≈ ", round(L - dx/2, digits=6))
end

"""
    validate_saddle_shape(points, grid_size, grid_length, grid_height)

Validate that the saddle shape has the expected properties.
"""
function validate_saddle_shape(points, grid_size, grid_length, grid_height)
    println("\nSaddle Shape Validation:")
    
    # Get all fixed (STATIC) nodes
    fixed_idx = [i for (i, p) in enumerate(points) if p.type == STATIC]
    
    if isempty(fixed_idx)
        println("ERROR: No fixed nodes found!")
        return false
    end
    
    x_fixed = [points[i].pos_w[1] for i in fixed_idx]
    y_fixed = [points[i].pos_w[2] for i in fixed_idx]
    z_fixed = [points[i].pos_w[3] for i in fixed_idx]
    
    # Check that we have boundary nodes at expected locations
    x_min, x_max = extrema(x_fixed)
    y_min, y_max = extrema(y_fixed)
    z_min, z_max = extrema(z_fixed)
    
    println("Fixed node ranges:")
    println("  X: $(round(x_min, digits=3)) to $(round(x_max, digits=3)) (expected: 0.0 to $(grid_length))")
    println("  Y: $(round(y_min, digits=3)) to $(round(y_max, digits=3)) (expected: 0.0 to $(grid_length))")
    println("  Z: $(round(z_min, digits=3)) to $(round(z_max, digits=3)) (expected: 0.0 to $(grid_height))")
    
    # Check expected count
    expected_fixed = 4 * grid_size - 4
    println("Fixed nodes: $(length(fixed_idx)) (expected: $expected_fixed)")
    
    # Validate saddle shape: z should vary from 0 to grid_height along x-direction
    x_range_ok = abs(x_min - 0.0) < 1e-3 && abs(x_max - grid_length) < 1e-3
    y_range_ok = abs(y_min - 0.0) < 1e-3 && abs(y_max - grid_length) < 1e-3
    z_range_ok = abs(z_min - 0.0) < 1e-3 && abs(z_max - grid_height) < 1e-3
    count_ok = length(fixed_idx) == expected_fixed
    
    success = x_range_ok && y_range_ok && z_range_ok && count_ok
    
    if success
        println("✓ Saddle shape validation PASSED")
    else
        println("✗ Saddle shape validation FAILED")
        !x_range_ok && println("  - X range incorrect")
        !y_range_ok && println("  - Y range incorrect") 
        !z_range_ok && println("  - Z range incorrect")
        !count_ok && println("  - Fixed node count incorrect")
    end
    
    return success
end

function plot3d_saddle_shape(points, segments; title="Saddle Form Shape")
    # Extract coordinates
    x = [p.pos_w[1] for p in points]
    y = [p.pos_w[2] for p in points]
    z = [p.pos_w[3] for p in points]
    
    # Debug: print some info
    println("Number of points: $(length(points))")
    println("Number of segments: $(length(segments))")
    
    # Check for NaNs/issues
    if any(isnan, x) || any(isnan, y) || any(isnan, z)
        println("WARNING: Found NaN values in pos_w coordinates!")
        println("Sample pos_w values:")
        for i in 1:min(5, length(points))
            println("  Point $i: pos_w = $(points[i].pos_w)")
        end
        return nothing
    end
    
    println("X range: $(minimum(x)) to $(maximum(x))")
    println("Y range: $(minimum(y)) to $(maximum(y))")
    println("Z range: $(minimum(z)) to $(maximum(z))")
    
    # Create base plot with all points
    p = Plots.scatter3d(x, y, z, 
                  markersize=3, 
                  markercolor=:blue,
                  markerstrokewidth=0,
                  title=title,
                  xlabel="X (m)", 
                  ylabel="Y (m)", 
                  zlabel="Z (m)",
                  legend=false)
    
    # Robust segment handling - handle both point_idxs and points fields
    getpair(seg) = hasproperty(seg, :point_idxs) ? seg.point_idxs : seg.points
    for seg in segments
        i, j = getpair(seg)
        Plots.plot3d!([x[i], x[j]], [y[i], y[j]], [z[i], z[j]], 
                color=:gray, alpha=0.6, linewidth=1)
    end
    
    # Highlight fixed nodes by type (STATIC), not by pre-saved indices
    fixed_idx = [i for (i, p) in enumerate(points) if p.type == STATIC]
    if !isempty(fixed_idx)
        x_fixed = x[fixed_idx]
        y_fixed = y[fixed_idx]
        z_fixed = z[fixed_idx]
        Plots.scatter3d!(x_fixed, y_fixed, z_fixed, 
                   markersize=5, 
                   markercolor=:red,
                   markerstrokewidth=0)
        println("Number of fixed (STATIC) nodes: $(length(fixed_idx))")
    end
    
    return p
end

# Alternative plotting function using raw positions
function plot3d_from_positions(positions, segments, fixed_nodes; title="Saddle Form Shape")
    x = getindex.(positions, 1)
    y = getindex.(positions, 2)
    z = getindex.(positions, 3)
    
    println("Plotting from raw positions:")
    println("X range: $(minimum(x)) to $(maximum(x))")
    println("Y range: $(minimum(y)) to $(maximum(y))")
    println("Z range: $(minimum(z)) to $(maximum(z))")
    
    # Create base plot with all points
    p = Plots.scatter3d(x, y, z, 
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
        Plots.plot3d!([x[i], x[j]], [y[i], y[j]], [z[i], z[j]], 
                color=:gray, alpha=0.6, linewidth=1)
    end
    
    # Highlight fixed nodes in red
    if !isempty(fixed_nodes)
        x_fixed = x[fixed_nodes]
        y_fixed = y[fixed_nodes]
        z_fixed = z[fixed_nodes]
        Plots.scatter3d!(x_fixed, y_fixed, z_fixed, 
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
function build_and_plot_saddle_shape(; grid_size=3, grid_length=10.0, grid_height=5.0,
                                     prestretch_frac=0.98, diameter_mm=3.0, use_3d=true)

    println("SADDLE FORM SHAPE - SymbolicAWEModels")
    set = build_settings()

    points, segments, fixed_nodes, dx, positions = build_saddle_points_segments(
        grid_size, grid_length, grid_height, set; prestretch_frac, diameter_mm)

    # Create a neutral transform that preserves the current orientation
    # This avoids unwanted rotations while still satisfying the transform requirements
    transforms = [neutral_transform(1; points, base_point_idx=1, rot_point_idx=2, 
                                    base_pos=[0.0, 0.0, 0.0])]
    
    sys = SymbolicAWEModels.SystemStructure("saddle_form", set; points, segments, transforms)
    sam = SymbolicAWEModel(set, sys)
    init!(sam; remake=false)  # Initialize to set pos_w properly

    println("Saddle Form Shape (N=$grid_size, prestretch=$(prestretch_frac), dx=$(round(dx,digits=3)) m)")
    
    # Use the initialized points from the sam object (after init!)
    pts = sam.sys_struct.points
    segs = sam.sys_struct.segments
    
    # Debug: Compare raw positions vs pos_w after transform
    println("\nTransform validation:")
    println("Sample point comparison (raw pos vs pos_w after init!):")
    for i in 1:min(3, length(pts))
        raw_pos = positions[i]
        init_pos = pts[i].pos_w
        diff = norm(init_pos .- raw_pos)
        println("  Point $i: raw=$(round.(raw_pos, digits=3)) → pos_w=$(round.(init_pos, digits=3)), diff=$(round(diff, digits=6))")
    end
    
    # Validate that fixed nodes are still in the right places
    fixed_idx_post = [i for (i, p) in enumerate(pts) if p.type == STATIC]
    println("Fixed nodes count: $(length(fixed_idx_post)) (should be $(4*grid_size - 4))")
    
    if length(fixed_idx_post) != 4*grid_size - 4
        println("WARNING: Fixed node count mismatch!")
    end
    
    # Optional: validate edges (for debugging)
    # Reconstruct edge indices from initialized points based on their positions
    top = [i for (i,p) in enumerate(pts) if p.type == STATIC && abs(p.pos_w[2] - 0.0) < 1e-6]
    bottom = [i for (i,p) in enumerate(pts) if p.type == STATIC && abs(p.pos_w[2] - grid_length) < 1e-6]  
    left = [i for (i,p) in enumerate(pts) if p.type == STATIC && abs(p.pos_w[1] - dx/2) < 1e-6]
    right = [i for (i,p) in enumerate(pts) if p.type == STATIC && abs(p.pos_w[1] - (grid_length - dx/2)) < 1e-6]
    check_edges(pts, top, bottom, left, right, dx, grid_length)
    
    # Validate the overall saddle shape
    validate_saddle_shape(pts, grid_size, grid_length, grid_height)
    
    # 3D plot using Plots.jl with initialized points
    if use_3d
        title_str = "Saddle Form (N=$grid_size, dx=$(round(dx,digits=3))m)"
        
        # Option 1: Plot with pos_w from initialized sam object
        fig3d = plot3d_saddle_shape(pts, segs; title=title_str)
        display(fig3d)
        
        # Option 2: Plot using raw positions (matches Python approach more closely)
        # This avoids any transform effects and shows the pure geometry
        println("\nAlternative plot using raw positions:")
        fig3d_raw = plot3d_from_positions(positions, segments, fixed_nodes; title=title_str * " (Raw)")
        display(fig3d_raw)
    end

    return sys
end

# # Run when file is executed directly or included
# if abspath(PROGRAM_FILE) == @__FILE__
#     build_and_plot_saddle_shape()
# end

# When included from menu.jl, run automatically
build_and_plot_saddle_shape(grid_size=3)
