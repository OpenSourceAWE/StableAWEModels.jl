# SPDX-FileCopyrightText: 2025 TU Delft V3 Kite Structural Model
#
# SPDX-License-Identifier: MPL-2.0

"""
TU Delft V3 Kite Structural Example

This example loads the TU Delft V3 kite structural geometry from a YAML file
and visualizes it using an interactive 3D PlotlyJS plot.

The geometry includes:
- Wing structural members (spars, ribs, etc.)
- Bridle line connections
- Fixed and dynamic attachment points
"""

using SymbolicAWEModels, VortexStepMethod, ControlPlots
using LinearAlgebra
# using Plots, PlotlyJS
using YAML
using CSV
using DataFrames

# Set PlotlyJS as the backend for interactive 3D visualization
# plotlyjs()

# Include the yaml_loader with all helper functions
include("../yaml_loader.jl")

println("\n\nTU Delft V3 Kite - Structural Model\n", "="^50)

# Model configuration
model_name = "TUDELFT_V3_KITE"

# Load settings for the V3 kite
println("Loading settings from: data/$model_name/")
set = SymbolicAWEModels.load_settings(model_name)

# Path to the structural geometry YAML file
geometry_path = joinpath(dirname(dirname(@__DIR__)), "data", model_name, "struc_geometry.yaml")

if !isfile(geometry_path)
    error("Geometry file not found at: $geometry_path")
end

println("Loading structural geometry from: $geometry_path")

# Load the system structure from YAML
sys = load_sys_struct_from_yaml(geometry_path; system_name=model_name, set=set)

# Fix pulley conflicts: The YAML has overlapping pulley definitions where node 34
# is used as both a pulley point and connection point in different pulleys.
# We need to remove the conflicting pulley from the system.
println("\n⚠️  Fixing pulley geometry conflicts...")

# Find segments that are part of multiple pulleys
segment_pulley_count = Dict{Int, Int}()
for pulley in sys.pulleys
    for seg_idx in pulley.segment_idxs
        segment_pulley_count[seg_idx] = get(segment_pulley_count, seg_idx, 0) + 1
    end
end

# Find segments with conflicts (part of 2+ pulleys)
conflicting_segments = [seg_idx for (seg_idx, count) in segment_pulley_count if count > 1]

if !isempty(conflicting_segments)
    println("  Found $(length(conflicting_segments)) segments in multiple pulleys: $conflicting_segments")
    
    # Remove pulleys that contain conflicting segments
    # Keep only the first pulley that references each segment
    segments_to_keep = Set{Int}()
    pulleys_to_remove = Int[]
    
    for (pulley_idx, pulley) in enumerate(sys.pulleys)
        has_conflict = any(seg_idx in conflicting_segments for seg_idx in pulley.segment_idxs)
        
        if has_conflict
            # Check if any of this pulley's segments are already kept
            already_kept = any(seg_idx in segments_to_keep for seg_idx in pulley.segment_idxs)
            
            if already_kept
                # This pulley conflicts with an earlier one - remove it
                push!(pulleys_to_remove, pulley_idx)
                println("  Removing pulley $pulley_idx (conflicts with earlier pulley)")
            else
                # Keep this pulley and mark its segments as kept
                union!(segments_to_keep, pulley.segment_idxs)
                println("  Keeping pulley $pulley_idx")
            end
        end
    end
    
    # Remove conflicting pulleys
    if !isempty(pulleys_to_remove)
        # Remove in reverse order to maintain correct indices
        for idx in reverse(sort(pulleys_to_remove))
            deleteat!(sys.pulleys, idx)
        end
        println("  Removed $(length(pulleys_to_remove)) conflicting pulleys")
    end
end

println("  Final pulley count: $(length(sys.pulleys))")

# Print summary information
println("\nStructural Model Summary:")
println("  Total points: $(length(sys.points))")
println("  Total segments: $(length(sys.segments))")

fixed_count = count(p -> p.type == SymbolicAWEModels.STATIC, sys.points)
dynamic_count = count(p -> p.type == SymbolicAWEModels.DYNAMIC, sys.points)
println("  Fixed points: $fixed_count")
println("  Dynamic points: $dynamic_count")

# Calculate bounding box
get_pos(p) = hasproperty(p, :pos_w) ? p.pos_w : p.pos_cad
x_coords = [get_pos(p)[1] for p in sys.points]
y_coords = [get_pos(p)[2] for p in sys.points]
z_coords = [get_pos(p)[3] for p in sys.points]

println("\nBounding Box:")
println("  X: [$(minimum(x_coords)), $(maximum(x_coords))] m")
println("  Y: [$(minimum(y_coords)), $(maximum(y_coords))] m")
println("  Z: [$(minimum(z_coords)), $(maximum(z_coords))] m")

##TODO: later should be turned-on again
# # Create interactive 3D visualization using PlotlyJS
# println("\nGenerating initial 3D plot...")
# plot3d_v3(sys.points, sys.segments, title="TU Delft V3 Kite - Structural Model (Initial)")

println("\nVisualization complete!")
println("  - Red diamonds: Fixed attachment points")
println("  - Blue circles: Dynamic structural nodes")
println("  - Black lines: Structural members (spars, ribs, bridle lines)")
println("  - Hover over points to see details")
println("  - Use mouse to rotate, zoom, and pan")

# ===== Load aerodynamic forces from CSV file =====
println("\n" * "="^60)
println("Loading aerodynamic forces from CSV")
println("="^60)

# Path to the loads CSV file
loads_csv_path = joinpath(dirname(dirname(@__DIR__)), "data", model_name, "vsm_computed_loads.csv")

if !isfile(loads_csv_path)
    error("Loads CSV file not found at: $loads_csv_path")
end

# Read the CSV file
loads_df = CSV.read(loads_csv_path, DataFrame)
println("Loaded $(nrow(loads_df)) force entries from CSV")

# Display summary of forces
println("\nForce Summary:")
println("  Total force X: $(round(sum(loads_df.force_x), digits=2)) N")
println("  Total force Y: $(round(sum(loads_df.force_y), digits=2)) N")
println("  Total force Z: $(round(sum(loads_df.force_z), digits=2)) N")
println("  Total magnitude: $(round(sum(loads_df.force_magnitude), digits=2)) N")

# ===== Apply forces to dynamic points =====
println("\nApplying forces to structural nodes...")

# Get dynamic points (excluding fixed points)
dynamic_points = filter(p -> p.type == SymbolicAWEModels.DYNAMIC, sys.points)
println("  Number of dynamic points: $(length(dynamic_points))")
println("  Number of force entries: $(nrow(loads_df))")

# Skip first force entry (node 1 is the fixed point)
# Map forces starting from row 2 to dynamic points
if length(dynamic_points) != (nrow(loads_df) - 1)
    @warn "Mismatch: $(length(dynamic_points)) dynamic points vs $(nrow(loads_df)-1) force entries (excluding fixed point)"
    println("  Will apply forces to first $(min(length(dynamic_points), nrow(loads_df)-1)) points")
end

# Apply forces from CSV to dynamic points (skip first row which is for fixed point)
# Note: disturb field is immutable SVector, so we need to copy the values element by element
let n_forces_applied = 0
    for (i, point) in enumerate(dynamic_points)
        if (i + 1) <= nrow(loads_df)  # +1 to skip first CSV row
            force_row = loads_df[i + 1, :]  # Start from row 2
            
            # Copy force values to the disturb field (element-wise assignment)
            point.disturb[1] = force_row.force_x
            point.disturb[2] = force_row.force_y
            point.disturb[3] = force_row.force_z
            
            n_forces_applied += 1
        end
    end
    
    println("  Successfully applied $n_forces_applied forces to dynamic points")
end

# ===== Run structural simulation =====
println("\n" * "="^60)
println("Running structural simulation with aerodynamic loads")
println("="^60)

# Simulation parameters
set.sample_freq = 10  # 10 Hz sampling (dt = 0.1s)
simulation_time = 5.0  # 5 seconds
n_steps = Int(simulation_time * set.sample_freq)

println("  Simulation time: $simulation_time seconds")
println("  Time step: $(1/set.sample_freq) seconds")
println("  Total steps: $n_steps")

# Analyze system for optimal damping
println("\nAnalyzing system dynamics...")
include("damping_analysis.jl")
freq, zeta, recommended_damping = analyze_damping_response(
    sys, set; verbose=true
)
println("  Setting recommended damping to: $recommended_damping")
set.damping = recommended_damping

# Create symbolic model
println("\nInitializing symbolic model...")
sam = SymbolicAWEModel(set, sys)
init!(sam; remake=false)

# Plot initial state with forces applied
plot(sam, 0.0; zoom=false)

# Run simulation
println("\nRunning simulation...")
for i in 1:n_steps
    current_time = i / set.sample_freq
    next_step!(sam)
    
    # Plot every 5 steps
    if i % 5 == 0
        plot(sam, current_time; zoom=false)
        println("  Step $i/$n_steps (t = $(round(current_time, digits=2))s)")
    end
end

# ===== Final results =====
println("\n" * "="^60)
println("Simulation Complete!")
println("="^60)

# Calculate displacement statistics
println("\nDisplacement Analysis:")
displacements = Float64[]
for point in dynamic_points
    initial_pos = point.pos_cad
    final_pos = point.pos_w
    displacement = norm(final_pos - initial_pos)
    push!(displacements, displacement)
end

println("  Mean displacement: $(round(mean(displacements), digits=4)) m")
println("  Max displacement: $(round(maximum(displacements), digits=4)) m")
println("  Min displacement: $(round(minimum(displacements), digits=4)) m")

# Final plot
println("\nGenerating final 3D plot...")
plot(sam, simulation_time; zoom=false)

nothing
