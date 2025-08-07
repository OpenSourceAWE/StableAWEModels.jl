# SPDX-FileCopyrightText: 2025 Bart van de Lint
#
# SPDX-License-Identifier: MPL-2.0

"""
Validation example: a massless pulley of mass m suspended between two fixed anchors by a 
massless, inextensible rope of total length L.  The pulley (treated here as a point‐mass) 
is released from an initial displacement and settles under gravity to a position that can 
be found analytically:

    2·√((L/2)² + z²) = L  ⇒  z = –√((L/2)² – (d/2)²)

where d is the horizontal separation of the two anchors.

This script compares the simulated equilibrium position to the analytic result.
"""

using SymbolicAWEModels, Plots

# ─── Problem Parameters ─────────────────────────────────────────────────────────

L_rope   = 6.0    # total rope length [m]
d_anchor = 4.0    # distance between the two fixed anchors [m]
m_pulley = 1.5    # mass of the pulley [kg]
g        = 9.807  # gravity [m/s²]

# Initial displacement: pulley starts lower than its resting sag by 20%
z_rest     = -sqrt((L_rope/2)^2 - (d_anchor/2)^2)
z_initial  = 1.2 * z_rest      # 20% “over-sagged”
x_initial  =  0.0              # symmetry ⇒ x = 0
abs_tol    = 1e-7
rel_tol    = 1e-7
damping    = 0.2               # tether damping to speed convergence

# ─── Build Settings ─────────────────────────────────────────────────────────────

set = Settings()
set.v_wind    = 0.0
set.l_tether  = L_rope         # total rope length
set.g         = g
set.abs_tol   = abs_tol
set.rel_tol   = rel_tol
set.c_tether  = damping        # damping in the tether elements

# ─── Define Points ──────────────────────────────────────────────────────────────

points = Point[]
# two fixed anchors, at z = 0
push!(points, Point(1, [-d_anchor/2, 0.0, 0.0], STATIC))
push!(points, Point(2, [ d_anchor/2, 0.0, 0.0], STATIC))
# the pulley node, dynamic, initially below the anchors
push!(points, Point(3, [x_initial, 0.0, z_initial], DYNAMIC; mass=m_pulley))

# ─── Define Rope Segments ───────────────────────────────────────────────────────

segments = Segment[]
# one segment from anchor 1 → pulley
push!(segments, Segment(1, set, (1, 3), TETHER))
# one segment from pulley → anchor 2
push!(segments, Segment(2, set, (3, 2), TETHER))

# ─── Define Pulley ──────────────────────────────────────────────────────────────

pulleys = Pulley[]
# this pulley “device” sits at node 3 and routes the rope between points 1 and 2
push!(pulleys, Pulley(1, (1, 2), DYNAMIC))

# ─── Assemble System ───────────────────────────────────────────────────────────

# no extra transforms needed for a simple ideal pulley
sys_struct = SymbolicAWEModels.SystemStructure(
    "pulley_validation", set;
    points, segments, pulleys
)

sam = SymbolicAWEModel(set, sys_struct)
init!(sam; remake=false)

# ─── Simulate to Equilibrium ────────────────────────────────────────────────────

n_steps = 1000
for i in 1:n_steps
    next_step!(sam)
end

# ─── Extract Final Position ────────────────────────────────────────────────────

sim_pos    = sam.x[1:3, 3]   # [x, y, z] of pulley node
analytic_z = z_rest
analytic_pos = [0.0, 0.0, analytic_z]

println("Simulated final position:   ", sim_pos)
println("Analytical resting position: ", analytic_pos)
println("Error norm: ", norm(sim_pos .- analytic_pos))

# ─── Plot Trajectory (Optional) ────────────────────────────────────────────────

# collect the pulley path during the run
xs = [sam.x_hist[1,3,i] for i in 1:size(sam.x_hist, 3)]
zs = [sam.x_hist[3,3,i] for i in 1:size(sam.x_hist, 3)]

plot(xs, zs, label="Trajectory", xlabel="x [m]", ylabel="z [m]")
scatter!([sim_pos[1]], [sim_pos[3]], label="Sim final", marker=:circle)
scatter!([analytic_pos[1]], [analytic_pos[3]], label="Analytical", marker=:star5)
title!("Pulley under gravity: validation")

# savefig("pulley_validation.png")  # uncomment to save

println("Done. The simulation’s final and analytical points should coincide.")
