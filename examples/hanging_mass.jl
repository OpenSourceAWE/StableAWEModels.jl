# SPDX-FileCopyrightText: 2025 Bart van de Lint
#
# SPDX-License-Identifier: MPL-2.0

using SymbolicAWEModels, VortexStepMethod, ControlPlots

set = Settings("system.yaml")
set.v_wind = 0  # No wind
set.abs_tol = 1e-6
set.rel_tol = 1e-6

points = Point[]
segments = Segment[]

# Create two points: anchor point (static) and hanging mass (dynamic)
push!(points, Point(1, [2.0, 0.0, 5.0], STATIC))          # Anchor point at height 5m
push!(points, Point(2, [2.0, 0.0, 2.0], DYNAMIC; mass=10.0)) # Hanging mass at height 2m, 1kg

# Create single segment connecting the points
push!(segments, Segment(1, set, (1, 2), BRIDLE; l0=2.0, diameter=1.5))  # 1.5mm diameter, 2m rest length

# Transform to position the system  
transforms = [Transform(1, -deg2rad(90.0), 0.0, 0.0; base_pos=[2.0, 0.0, 5.0], base_point_idx=1, rot_point_idx=2)]

# Create system structure
sys_struct = SymbolicAWEModels.SystemStructure("hanging_mass", set; points, segments, transforms)

# Plot initial state
plot(sys_struct, 0.0; zoom=false)

# Create and initialize the symbolic model
sam = SymbolicAWEModel(set, sys_struct)
init!(sam; remake=false)

# Run simulation for 100 steps
for i in 1:1000
    plot(sam, i/set.sample_freq; zoom=false)
    next_step!(sam)
end