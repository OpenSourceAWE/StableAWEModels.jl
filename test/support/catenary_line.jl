# SPDX-FileCopyrightText: 2025 Bart van de Lint
#
# SPDX-License-Identifier: MPL-2.0

"""
Validation example: tether fixed at both ends, under perpendicular gravity,
resulting in a catenary line.

This script sets up a simple system with a single tether, fixed at two points
at equal height, under gravity. The simulation compares the final shape to
the analytical catenary solution.

Requirements:
- SymbolicAWEModels
- ControlPlots

"""

using SymbolicAWEModels, ControlPlots

# --- Settings
set = Settings()
set.v_wind = 0.0               # No wind for pure catenary
set.l_tether = 10.0            # Tether length [m]
set.abs_tol = 1e-6
set.rel_tol = 1e-6
set.g = 9.807
set.n_segments = 12            # Number of tether segments (n points = n_segments+1)
set.rho_tether = 0.1           # Tether mass density [kg/m]
set.tether_diam = 0.01         # Tether diameter [m]
set.c_tether = 0.1             # Damping, can tune for faster settling

# --- Points (nodes)
points = Point[]
push!(points, Point(1, [0.0, 0.0, 0.0], STATIC))           # Anchor left
push!(points, Point(2, [set.l_tether, 0.0, 0.0], STATIC))  # Anchor right

# Add dynamic points along initial straight line
for i in 1:set.n_segments-1
    x = i * set.l_tether / set.n_segments
    push!(points, Point(i+2, [x, 0.0, 0.0], DYNAMIC))
end

# --- Segments (springs/dampers)
segments = Segment[]
for i in 1:set.n_segments
    # Connect consecutive points
    push!(segments, Segment(i, set, (i, i+1), TETHER))
end

# --- System structure
sys_struct = SymbolicAWEModels.SystemStructure("tether_catenary", set; points, segments)

# --- Construct symbolic model
sam = SymbolicAWEModel(set, sys_struct)

# --- Initialize system
init!(sam; remake=false)

# --- Simulate until static equilibrium
n_steps = 1000
for i in 1:n_steps
    next_step!(sam)
    # Optionally plot the system at intervals for debugging/animation
    # if i % 50 == 0
    #     plot(sam, i/set.sample_freq; zoom=false)
    # end
end

# --- Extract simulated point positions (final time)
sim_points = [sam.x[1:3, j] for j in 1:length(points)]

# --- Analytical catenary solution (for plotting)
function analytical_catenary(L::Float64, sag::Float64=undef; n=100)
    # L: horizontal distance between anchors
    # sag: vertical sag at midpoint (if not given, estimate from sim)
    # Returns arrays of x and z positions
    if sag === undef
        sag = maximum([p[3] for p in sim_points])  # Use sim midpoint
    end
    a = (0.25*L^2 - sag^2)/(2sag)
    x = range(0, L; length=n)
    z = a * cosh.((x .- 0.5L) ./ a)
    z .-= z[1]  # Shift so that z=0 at left anchor
    return x, z
end

L = set.l_tether
n = length(sim_points)
sim_x = [p[1] for p in sim_points]
sim_z = [p[3] for p in sim_points]
x_anal, z_anal = analytical_catenary(L, minimum(sim_z); n=200)

# --- Plot simulation vs analytical catenary
using Plots
plot(sim_x, sim_z, lw=3, label="Simulation", xlabel="x [m]", ylabel="z [m]", legend=:bottom)
plot!(x_anal, z_anal, lw=2, ls=:dash, label="Analytical catenary")
title!("Tether under gravity: catenary validation")

# Optionally save the plot
# savefig("catenary_validation.png")

println("Done! The simulated tether shape should match the analytical catenary.")

