# SPDX-FileCopyrightText: 2025 TU Delft V3 Kite Structural Model
#
# SPDX-License-Identifier: MPL-2.0

"""
TU Delft V3 Kite Structural Example

This example loads the TU Delft V3 kite structural geometry from a YAML file
and visualizes it using an interactive 3D Makie plot.

The geometry includes:
- Wing structural members (spars, ribs, etc.)
- Bridle line connections
- Fixed and dynamic attachment points
"""

using SymbolicAWEModels, VortexStepMethod
using LinearAlgebra
using YAML
using StaticArrays

# External aerodynamic loads for TUDELFT V3 kite (array rows match node indices)
const F_AERO_WING = [
     0.00000000e+00   0.00000000e+00   0.00000000e+00;
    -2.55684556e+00   1.95970467e+01   2.87074883e+00;
    -3.01386405e-02   2.22824114e-01   3.48008143e-02;
    -1.89232124e+01   9.58494908e+01   6.18829102e+01;
    -2.41222676e+00   1.23422296e+01   7.95456719e+00;
    -4.30387165e+01   1.37202855e+02   1.83362930e+02;
    -6.13500886e+00   1.97691945e+01   2.59619180e+01;
    -5.31750877e+01   9.56172558e+01   2.78414677e+02;
    -6.70024895e+00   1.27649189e+01   3.50790650e+01;
    -4.87888099e+01   3.21247890e+01   3.26042612e+02;
    -5.40699135e+00   3.21478649e+00   3.45192279e+01;
    -4.87888099e+01  -3.21247890e+01   3.26042612e+02;
    -5.40699135e+00  -3.21478649e+00   3.45192279e+01;
    -5.31750877e+01  -9.56172558e+01   2.78414677e+02;
    -6.70024895e+00  -1.27649189e+01   3.50790650e+01;
    -4.30387165e+01  -1.37202855e+02   1.83362930e+02;
    -6.13500886e+00  -1.97691945e+01   2.59619180e+01;
    -1.89232124e+01  -9.58494908e+01   6.18829102e+01;
    -2.41222676e+00  -1.23422296e+01   7.95456719e+00;
    -2.55684556e+00  -1.95970467e+01   2.87074883e+00;
    -3.01386405e-02  -2.22824114e-01   3.48008143e-02;
     0.00000000e+00   0.00000000e+00   0.00000000e+00;
     0.00000000e+00   0.00000000e+00   0.00000000e+00;
     0.00000000e+00   0.00000000e+00   0.00000000e+00;
     0.00000000e+00   0.00000000e+00   0.00000000e+00;
     0.00000000e+00   0.00000000e+00   0.00000000e+00;
     0.00000000e+00   0.00000000e+00   0.00000000e+00;
     0.00000000e+00   0.00000000e+00   0.00000000e+00;
     0.00000000e+00   0.00000000e+00   0.00000000e+00;
     0.00000000e+00   0.00000000e+00   0.00000000e+00;
     0.00000000e+00   0.00000000e+00   0.00000000e+00;
     0.00000000e+00   0.00000000e+00   0.00000000e+00;
     0.00000000e+00   0.00000000e+00   0.00000000e+00;
     0.00000000e+00   0.00000000e+00   0.00000000e+00;
     0.00000000e+00   0.00000000e+00   0.00000000e+00;
     0.00000000e+00   0.00000000e+00   0.00000000e+00;
     0.00000000e+00   0.00000000e+00   0.00000000e+00;
     0.00000000e+00   0.00000000e+00   0.00000000e+00;
     0.00000000e+00   0.00000000e+00   0.00000000e+00
]

# Simulation parameters
const N_PLOTS = 3  # number of static snapshots to render
n_steps = 50
remake_cache = false

# Include the yaml_loader with all helper functions
include("../yaml_loader.jl")
try
    include("../plotly_plots.jl")  # provides plot3d_v3 and refresh_plot3d!
catch
    @warn "plotly_plots.jl not found relative to this script; expecting plot3d_v3 and refresh_plot3d! to be available in scope."
end

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

# Apply provided static loads to system points
println("\nApplying provided static loads to structural nodes...")
for (i, point) in enumerate(sys.points)
    if i <= size(F_AERO_WING, 1)
        point.disturb[1] = F_AERO_WING[i, 1]
        point.disturb[2] = F_AERO_WING[i, 2]
        point.disturb[3] = F_AERO_WING[i, 3]
    end
end
println("  Successfully applied static loads to $(min(length(sys.points), size(F_AERO_WING, 1))) points")

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
# # Create interactive 3D visualization using Makie
# println("\nGenerating initial 3D plot...")
# plot3d_v3(sys.points, sys.segments, title="TU Delft V3 Kite - Structural Model (Initial)")

println("\nVisualization complete!")
println("  - Red diamonds: Fixed attachment points")
println("  - Blue circles: Dynamic structural nodes")
println("  - Black lines: Structural members (spars, ribs, bridle lines)")
println("  - Hover over points to see details")
println("  - Use mouse to rotate, zoom, and pan")

# ===== Run structural simulation =====
println("\n" * "="^60)
println("Running structural simulation with aerodynamic loads")
println("="^60)

# Simulation parameters
set.sample_freq = 10  # 10 Hz sampling (dt = 0.1s)

Δt = 1.0 / set.sample_freq
println("  Simulation time: $(n_steps * Δt) seconds")
println("  Time step: $Δt seconds")
println("  Total steps: $n_steps")

# Analyze system for optimal damping
println("\nAnalyzing system dynamics...")
include("damping_analysis.jl")
freq, zeta, recommended_damping = analyze_damping_response(
    sys, set; verbose=true
)
println("  Setting recommended damping to: $recommended_damping")
set.unit_damping = recommended_damping

# Create symbolic model
println("\nInitializing symbolic model...")
sam = SymbolicAWEModel(set, sys)

# Ensure provided static loads are set on the model structure before initialization
for (i, point) in enumerate(sam.sys_struct.points)
    if i <= size(F_AERO_WING, 1)
        point.disturb[1] = F_AERO_WING[i, 1]
        point.disturb[2] = F_AERO_WING[i, 2]
        point.disturb[3] = F_AERO_WING[i, 3]
    end
end

init!(sam; remake=remake_cache)

# Re-apply provided static loads after initialization
for (i, point) in enumerate(sam.sys_struct.points)
    if i <= size(F_AERO_WING, 1)
        point.disturb[1] = F_AERO_WING[i, 1]
        point.disturb[2] = F_AERO_WING[i, 2]
        point.disturb[3] = F_AERO_WING[i, 3]
    end
end

# Determine which steps to capture for static plots
num_samples = max(N_PLOTS, 2)
snapshot_steps = unique!(sort!(round.(Int, range(0, stop=n_steps, length=num_samples))))
snapshot_steps[1] != 0 && pushfirst!(snapshot_steps, 0)
snapshot_steps[end] != n_steps && push!(snapshot_steps, n_steps)

# Store snapshots
snapshots = Dict{Int, Vector{SymbolicAWEModels.Point}}(0 => deepcopy(sam.sys_struct.points))

# Run simulation
println("\nRunning simulation...")
for step in 1:n_steps
    t = step * Δt
    
    # Re-apply provided static loads at each step
    for (i, point) in enumerate(sam.sys_struct.points)
        if i <= size(F_AERO_WING, 1)
            point.disturb[1] = F_AERO_WING[i, 1]
            point.disturb[2] = F_AERO_WING[i, 2]
            point.disturb[3] = F_AERO_WING[i, 3]
        end
    end
    
    next_step!(sam)
    
    # Store current state if requested
    if step in snapshot_steps
        snapshots[step] = deepcopy(sam.sys_struct.points)
    end
    
    # Print progress
    if step % 20 == 0
        println("  Step $step/$n_steps (t = $(round(t, digits=2))s)")
    end
end

# Ensure final state is captured
snapshots[n_steps] = get(snapshots, n_steps, deepcopy(sam.sys_struct.points))

captured_steps = sort!(collect(keys(snapshots)))

# ===== Final results =====
println("\n" * "="^60)
println("Simulation Complete!")
println("="^60)

# Get dynamic points for displacement analysis
dynamic_points = filter(p -> p.type == SymbolicAWEModels.DYNAMIC, sys.points)

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

# Render static plots for captured snapshots
println("\nRendering $(length(captured_steps)) static plots...")
for (idx, step) in enumerate(captured_steps)
    points_snapshot = snapshots[step]
    t = step * Δt
    plot_title = "TUDELFT V3 Kite – Step $(step) (t=$(round(t, digits=2)) s)"
    plot3d_v3(points_snapshot, sam.sys_struct.segments; title=plot_title)
    println("  Rendered static snapshot $(idx)/$(length(captured_steps)) at step $step")
end

println("\nDone. Created $(length(captured_steps)) static Makie plots covering the simulation window.")

nothing
