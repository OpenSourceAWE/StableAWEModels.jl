# SPDX-FileCopyrightText: 2025 Bart van de Lint
#
# SPDX-License-Identifier: MPL-2.0

#TODO: this function should have a:
# - compute_saddle_form function
# - plot_saddle_form function
# - run_saddle_form function that runs when this script is executed

"""
Validation example: self-stressed “saddle” mesh.

A rectangular grid of spring–damper elements, with all boundary nodes fixed,  
and all interior nodes initially at mid-height. We “self-stress” the network  
(i.e. springs have zero rest length), then let it relax under no external loads  
to find the saddle shape. The final positions of the interior nodes are plotted  
against the initial flat layout.
"""

using SymbolicAWEModels, Plots

# ─── Problem Parameters ─────────────────────────────────────────────────────────

const grid_size   = 5       # number of nodes along one side
const grid_length = 10.0    # horizontal span of the grid [m]
const grid_height =  5.0    # maximum height of boundary nodes [m]
const k_t         =  1.0    # spring stiffness for all edges [N/m]
const c_d         =  1.0    # damping coefficient [N·s/m]
const m_node      =  1.0    # mass of each interior node [kg]
const abs_tol     = 1e-8
const rel_tol     = 1e-8
const n_steps     = 500

# ─── Connectivity & Initial Layout ──────────────────────────────────────────────

"""
Builds the list of fixed-node indices (1-based) and the list of
edge‐connections as pairs of node indices.
"""
function build_mesh(grid_size)
    # total node count: square grid + staggered
    n = grid_size^2 + (grid_size - 1)^2

    # identify boundary nodes
    top    = collect(1:grid_size)
    bottom = collect(n-grid_size+1:n)
    left   = [ (2*grid_size-1)*(i-1)+1 for i in 2:grid_size-1 ]
    right  = [ i + grid_size-1 for i in left ]
    fixed_nodes = union(top, bottom, left, right)

    # connect each interior node to its four neighbors
    edges = Tuple{Int,Int}[]
    for i in 1:n
        if i ∉ fixed_nodes
            # neighbor offsets in the unpacked grid
            push!(edges, (i, i - grid_size))
            push!(edges, (i, i - grid_size + 1))
            push!(edges, (i, i + grid_size - 1))
            push!(edges, (i, i + grid_size))
        end
    end

    # remove duplicates (i,j) vs (j,i)
    unique_edges = Set{Tuple{Int,Int}}()
    for (i,j) in edges
        if (j,i) ∉ unique_edges
            push!(unique_edges, (i,j))
        end
    end

    return collect(unique_edges), fixed_nodes
end

"""
Returns the initial 3D coordinates for each of the n nodes,
and a Boolean mask for which nodes are fixed.
"""
function build_initial_positions(grid_size, fixed_nodes, grid_length, grid_height)
    # planar layout, alternating rows
    d  = grid_length/(grid_size-1)
    coords = Vector{SVector{3,Float64}}()
    # build planar X–Y coords
    xy = Vector{SVector{2,Float64}}()
    even = [ i*d for i in 0:grid_size-1 ]
    odd  = [ i*d + d/2 for i in 0:grid_size-2 ]
    for row in 0:grid_size-1
        if row % 2 == 0
            for x in even; push!(xy, SVector(x, row*d)) end
        else
            for x in odd;  push!(xy, SVector(x, row*d)) end
        end
    end
    # two halves: top and bottom
    xy = vcat(xy, xy[1:end-1])  # matches n = grid_size^2 + (grid_size-1)^2

    # build Z: boundary nodes at ±grid_height, interior at 0
    z = Vector{Float64}(undef, length(xy))
    for idx in 1:length(xy)
        if idx in fixed_nodes
            # assign height along the boundary in order: top→bottom→left→right
            # simplest: top and bottom at +grid_height/2 and –grid_height/2, others too
            z[idx] = (xy[idx][2] < grid_length/2 ? +grid_height/2 : -grid_height/2)
        else
            z[idx] = 0.0
        end
    end

    # combine to 3D SVector
    positions = [ SVector(xy[i]..., z[i]) for i in 1:length(xy) ]
    is_fixed  = [ i in fixed_nodes for i in 1:length(xy) ]

    return positions, is_fixed
end

# Build mesh and initial state
edges, fixed_nodes = build_mesh(grid_size)
positions, is_fixed = build_initial_positions(grid_size, fixed_nodes, grid_length, grid_height)
n_nodes = length(positions)

# ─── Build SymbolicAWEModels System ─────────────────────────────────────────────

# 1) Settings
set = Settings()
set.abs_tol  = abs_tol
set.rel_tol  = rel_tol
set.g        = 0.0          # no external gravity
set.v_wind   = SVector(0.,0.,0.)
set.k_tether = k_t
set.c_tether = c_d

# 2) Points
points = Point[]
for i in 1:n_nodes
    status = is_fixed[i] ? STATIC : DYNAMIC
    kwargs = is_fixed[i] ? () : (; mass=m_node)
    push!(points, Point(i, positions[i], status; kwargs...))
end

# 3) Segments
segments = Segment[]
for (idx,(i,j)) in enumerate(edges)
    push!(segments, Segment(idx, set, (i, j), TETHER))
end

# 4) System structure
sys_struct = SymbolicAWEModels.SystemStructure(
    "saddle_form", set;
    points, segments
)

# 5) Model and initialization
sam = SymbolicAWEModel(set, sys_struct)
init!(sam; remake=false)

# ─── Self-stress & Relaxation ────────────────────────────────────────────────────

stress_self!(sam)

for step in 1:n_steps
    next_step!(sam)
end

# ─── Extract & Plot ──────────────────────────────────────────────────────────────

# initial vs final positions of interior nodes
init_pos = [ positions[i]       for i in 1:n_nodes if !is_fixed[i] ]
final_pos = [ sam.x[1:3,i] |> SVector for i in 1:n_nodes if !is_fixed[i] ]

xs_i = [p[1] for p in init_pos];  zs_i = [p[3] for p in init_pos]
xs_f = [p[1] for p in final_pos]; zs_f = [p[3] for p in final_pos]

plot(xs_i, zs_i, seriestype=:scatter, marker=:circle, label="Initial", xlabel="x [m]", ylabel="z [m]")
scatter!(xs_f, zs_f, marker=:star5, label="Relaxed")
title!("Saddle Form: initial vs self-stressed mesh")

println("Saddle form relaxation complete.")
