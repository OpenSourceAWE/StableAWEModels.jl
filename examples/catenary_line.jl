# SPDX-FileCopyrightText: 2025 Bart van de Lint
#
# SPDX-License-Identifier: MPL-2.0

"""
Validation example: tether fixed at both ends, under gravity,
resulting in a catenary line.

This script sets up a simple system with multiple tether segments, fixed at two points
at equal height, under gravity. The simulation compares the final shape to
the analytical catenary solution.

Requirements:
- SymbolicAWEModels
- ControlPlots
"""

using SymbolicAWEModels, VortexStepMethod, ControlPlots

# --- Settings
set = load_settings("base")  # Loads from data/base/settings.yaml
set.v_wind = 0.0               # No wind for pure catenary
set.abs_tol = 1e-8
set.rel_tol = 1e-8
set.l_tether = 8        # Set tether length for plot size

# Catenary parameters
horizontal_span = 8          # Horizontal distance between anchors [m]
n_segments = 10                # Number of tether segments (n points = n_segments+1)
total_length = 10.0            # Total unstretched length of tether [m]
segment_mass = 0.1             # Mass per segment [kg]

# --- Points (nodes)
points = Point[]
push!(points, Point(1, [0.0, 0.0, 5.0], STATIC))                               # Left anchor
push!(points, Point(2, [horizontal_span, 0.0, 5.0], STATIC))                  # Right anchor

# Add dynamic points along initial straight line
for i in 1:n_segments-1
    x = i * horizontal_span / n_segments
    push!(points, Point(i+2, [x, 0.0, 5.0], DYNAMIC; mass=segment_mass))
end

# --- Segments (springs/dampers)
segments = Segment[]
l0_per_segment = total_length / n_segments  # Rest length per segment

for i in 1:n_segments
    # Connect consecutive points
    push!(segments, Segment(i, set, (i, i+1), POWER_LINE; l0=l0_per_segment, diameter_mm=4.0))
end

# --- Transforms
transforms = [Transform(1, 0.0, 0.0, 0.0; base_pos=[0.0, 0.0, 5.0], base_point_idx=1, rot_point_idx=2)]

# --- System structure
sys_struct = SymbolicAWEModels.SystemStructure("catenary", set; points, segments, transforms)

# Plot initial state
plot(sys_struct, 0.0; zoom=false)

# --- Construct symbolic model
sam = SymbolicAWEModel(set, sys_struct)
init!(sam; remake=false)

# --- Simulate until static equilibrium
println("Simulating catenary formation...")
n_steps = 2000
for i in 1:n_steps
    next_step!(sam)
    # Plot every 100 steps to show evolution
    if i % 200 == 0
        plot(sam, i/set.sample_freq; zoom=false)
        println("Step $i/$n_steps")
    end
end

# --- Final plot
plot(sam, n_steps/set.sample_freq; zoom=false)
println("Simulation complete. The tether should now show a catenary shape under gravity.")

# --- Extract final positions for validation
println("\nFinal point positions:")
for i in 1:length(points)
    pos = sam.sys_state.X[3*(i-1)+1:3*i]
    println("Point $i: x=$(round(pos[1], digits=3)), z=$(round(pos[3], digits=3))")
end