# SPDX-FileCopyrightText: 2025 Bart van de Lint
#
# SPDX-License-Identifier: MPL-2.0

using SymbolicAWEModels, VortexStepMethod, ControlPlots

println("\n\nHanging Mass Example\n", "="^40)
### Loading Settings
# Use load_settings() to temporarily create system.yaml pointing to the desired subdirectory
set = load_settings("base")  # Loads from data/base/settings.yaml
set.v_wind = 0  # No wind
set.sample_freq = 1  # Increase to 100 Hz for better visualization (dt = 0.01s)
set.abs_tol = 1e-6     # Higher precision for better dynamics resolution
set.rel_tol = 1e-6     # Higher precision for better dynamics resolution


# Create two points: anchor point (static) and hanging mass (dynamic)
points = Point[]
push!(points, Point(1, [2.0, 0.0, 5.0], STATIC))          # Anchor point at height 5m
push!(points, Point(2, [2.0, 0.0, 2], DYNAMIC; mass=1.0)) # Hanging mass at height 2m, 1kg

### Create single segment connecting the points
# l0 is the rest length a bit shorter than the distances between the initial points
# compression_frac is set to 0.001, meaning the spring has 0.1% compresive stiffness compared to elongation stiffness
# diameter_mm is the diameter of the bridle segment in millimeters
# As the same E modulus (e_tether) is used, this determines the stiffness:
#     axial_stiffness = set.e_tether * (diameter_m/2)^2 * π
# and the damping:
#     axial_damping = (set.damping / set.c_spring) * axial_stiffness
# where the set. refers to defined values in settings.yaml
segments = Segment[]
push!(segments, Segment(1, set, (1, 2), BRIDLE; l0=4, compression_frac=0.001, diameter_mm=5))  # 5mm diameter, 4m rest length

### Transform to position the system
# The base position is set to [2.0, 0.0, 5.0], which is the anchor point position.
# The rot_point_idx is set to 2, which refers to the hanging mass point.
# The orientation from base to rot point is vertical, meaning the Z-axis points downwards.
# To transfer this back to an x-axis aligned with elev=0 and azimuth=0, we need to rotate the system.
# This is done by the Transform constructor: using -90 degrees around the Z-axis and translating it to the anchor point position.

transforms = [Transform(1, -deg2rad(90.0), 0.0, 0.0; base_pos=[2.0, 0.0, 5.0], base_point_idx=1, rot_point_idx=2)]

### Create system structure
# The system structure consists of:
# - name: "hanging_mass"
# - settings: `set`
# - points: `points`
# - segments: `segments`
# - transforms: `transforms`

sys_struct = SymbolicAWEModels.SystemStructure("hanging_mass", set; points, segments, transforms)

### Analyze damping response
include("damping_analysis.jl")
freq, zeta, recommended_damping = analyze_damping_response(sys_struct, set; verbose=true)
println("\n Setting recommended damping to: ", recommended_damping)
set.damping = recommended_damping  # Update settings with recommended damping

### Plot initial state
# even though the tether is not used here, it defines the size of the plot
# and therefore we must set it to a reasonable length
set.l_tether = 5.0
plot(sys_struct, 0.0; zoom=false)

# Create and initialize the symbolic model
sam = SymbolicAWEModel(set, sys_struct)
init!(sam; remake=false)

# Run simulation for longer time with smaller steps for better convergence visualization
for i in 1:30
    plot(sam, i/set.sample_freq; zoom=false)
    next_step!(sam)
end