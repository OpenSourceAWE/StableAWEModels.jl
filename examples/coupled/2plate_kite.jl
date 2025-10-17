# Copyright (c) 2025 2-Plate Kite Coupled Simulation
# SPDX-License-Identifier: MPL-2.0

"""
2-Plate Kite Coupled Aerodynamic-Structural Simulation

This example loads the 2-plate kite model with aerodynamics and runs
a time-marching simulation with coupled aerodynamic-structural updates.
"""

using SymbolicAWEModels
using VortexStepMethod
using LinearAlgebra
using YAML
using StaticArrays

@info "Loading 2-plate kite model..."

# ------------------------------------------------------------------------------
# Helper to recompute aerodynamic loads with a fresh non-linear VSM solve.
# This bypasses the default linearised update so the coupled model always uses
# forces consistent with the current apparent wind state.
# ------------------------------------------------------------------------------
function refresh_vsm_forces!(sam::SymbolicAWEModel)
    prob = sam.prob
    isempty(sam.sys_struct.wings) && return
    vsm_y = prob.get_vsm_y(sam.integrator)
    groups = sam.sys_struct.groups

    for wing in sam.sys_struct.wings
        y_state = @view vsm_y[wing.idx, :]
        va = SVector{3, Float64}(y_state[1], y_state[2], y_state[3])
        omega = SVector{3, Float64}(y_state[4], y_state[5], y_state[6])

        VortexStepMethod.set_va!(wing.vsm_aero, va, omega)

        n_groups_total = length(groups)
        if n_groups_total > 0 && !isempty(wing.group_idxs)
            twist_angles = Vector{Float64}(undef, length(wing.group_idxs))
            for (local_idx, group_idx) in enumerate(wing.group_idxs)
                twist_angles[local_idx] = y_state[6 + Int(group_idx)]
            end
            VortexStepMethod.group_deform!(wing.vsm_wing, twist_angles, nothing; smooth=false)
            VortexStepMethod.reinit!(wing.vsm_aero; init_aero=false)
        end

        if n_groups_total > 0 && !isempty(wing.group_idxs)
            moment_frac = groups[Int(first(wing.group_idxs))].moment_frac
            VortexStepMethod.solve!(wing.vsm_solver, wing.vsm_aero; log=false, moment_frac=moment_frac)
        else
            VortexStepMethod.solve!(wing.vsm_solver, wing.vsm_aero; log=false)
        end

        sol = wing.vsm_solver.sol
        wing.vsm_x[1:3] .= sol.force
        wing.vsm_x[4:6] .= sol.moment
        if n_groups_total > 0 && !isempty(wing.group_idxs)
            for (local_idx, group_idx) in enumerate(wing.group_idxs)
                wing.vsm_x[6 + local_idx] = sol.group_moment_dist[Int(group_idx)]
            end
        end
        wing.vsm_jac .= 0.0
        wing.vsm_y .= y_state
    end

    prob.set_sys(sam.integrator, sam.sys_struct)
    return nothing
end

# ============= User settings =============
const MODEL_NAME = "2plate_kite"
const GEOM_PATH  = joinpath("data", MODEL_NAME, "struc_geometry.yaml")
const SIM_TIME   = 10.0  # Total simulation time in seconds
const N_PLOTS    = 3      # Number of static snapshots to render
const REMAKE_CACHE = false
# =========================================

# Include helper utilities first
yaml_loader_path = joinpath(@__DIR__, "..", "yaml_loader.jl")
if isfile(yaml_loader_path)
    include(yaml_loader_path)  # provides load_sys_struct_from_yaml
else
    @warn "yaml_loader.jl not found; expecting load_sys_struct_from_yaml to be available"
end

plotly_helpers_path = joinpath(@__DIR__, "..", "plotly_plots.jl")
if isfile(plotly_helpers_path)
    include(plotly_helpers_path)  # provides plot3d_v3
else
    @warn "plotly_plots.jl not found; plotting may not work"
end

# Load settings for the 2-plate kite
set = SymbolicAWEModels.load_settings(MODEL_NAME)

# Load system structure from YAML (not from .obj files)
@info "Building system structure from YAML: $GEOM_PATH"
isfile(GEOM_PATH) || error("Geometry file not found at: $GEOM_PATH")
sys = load_sys_struct_from_yaml(GEOM_PATH; system_name=MODEL_NAME, set=set)

# Create symbolic model with aerodynamics using the YAML-loaded structure
@info "Creating SymbolicAWEModel with aerodynamics..."
sam = SymbolicAWEModel(set, sys)

# Initialize the model
@info "Initializing model..."
SymbolicAWEModels.init!(sam; remake=REMAKE_CACHE)

# Calculate simulation parameters
Δt = 1.0 / max(1, hasproperty(set, :sample_freq) ? set.sample_freq : 100)
n_steps = round(Int, SIM_TIME / Δt)
@info "Running simulation for $(SIM_TIME) seconds ($n_steps steps, Δt = $(round(Δt, digits=4)) s)..."

# Determine which steps to capture for plotting (include start and end)
num_samples = max(N_PLOTS, 2)
snapshot_steps = unique!(sort!(round.(Int, range(0, stop=n_steps, length=num_samples))))
snapshot_steps[1] != 0 && pushfirst!(snapshot_steps, 0)
snapshot_steps[end] != n_steps && push!(snapshot_steps, n_steps)

# Store snapshots
snapshots = Dict{Int, Vector{Point}}(0 => deepcopy(sam.sys_struct.points))

# Plot initial state
@info "Plotting initial state..."
try
    plot3d_v3(sam.sys_struct.points, sam.sys_struct.segments; 
              title="2-Plate Kite - Initial State (t=0)")
catch e
    @warn "Could not plot initial state: $e"
end

# Time-marching loop with coupled aerodynamics
@info "Starting time-marching simulation with coupled aerodynamics..."

for step in 1:n_steps
    t = step * Δt

    refresh_vsm_forces!(sam)

    # Advance simulation one step using freshly computed VSM loads (linear refresh disabled)
    next_step!(sam; dt=Δt, vsm_interval=0)

    # Store current state if requested for plotting
    if step in snapshot_steps
        snapshots[step] = deepcopy(sam.sys_struct.points)
    end

    # Print progress periodically
    if step % 100 == 0 || step == n_steps
        @info "  Step $step/$n_steps (t = $(round(t, digits=2)) s / $(SIM_TIME) s)"
    end
end

# Ensure final state is captured
snapshots[n_steps] = get(snapshots, n_steps, deepcopy(sam.sys_struct.points))

captured_steps = sort!(collect(keys(snapshots)))

@info "Simulation complete. Rendering $(length(captured_steps)) static plots..."

for (idx, step) in enumerate(captured_steps)
    points_snapshot = snapshots[step]
    t = step * Δt
    plot_title = "2-Plate Kite (Coupled Aero) – Step $(step) (t=$(round(t, digits=2)) s)"
    try
        plot3d_v3(points_snapshot, sam.sys_struct.segments; title=plot_title)
        @info "  Rendered static snapshot $(idx)/$(length(captured_steps)) at step $step"
    catch e
        @warn "  Could not render snapshot at step $step: $e"
    end
end

# Print final statistics
println("\n" * "="^60)
println("Final Simulation Results (t = $(SIM_TIME) s)")
println("="^60)

# Calculate wing position and orientation if available
if length(sam.sys_struct.wings) > 0
    wing = sam.sys_struct.wings[1]
    pos = wing.pos_w
    println("  Wing position: [$(round(pos[1], digits=2)), $(round(pos[2], digits=2)), $(round(pos[3], digits=2))] m")
    println("  Elevation: $(round(rad2deg(wing.elevation), digits=2))°")
    println("  Azimuth: $(round(rad2deg(wing.azimuth), digits=2))°")
    println("  Heading: $(round(rad2deg(wing.heading), digits=2))°")
end

# Calculate average position of dynamic points
dynamic_points = filter(p -> p.type == SymbolicAWEModels.DYNAMIC, sam.sys_struct.points)
if !isempty(dynamic_points)
    avg_x = mean([p.pos_w[1] for p in dynamic_points])
    avg_y = mean([p.pos_w[2] for p in dynamic_points])
    avg_z = mean([p.pos_w[3] for p in dynamic_points])
    println("  Average dynamic point position: [$(round(avg_x, digits=2)), $(round(avg_y, digits=2)), $(round(avg_z, digits=2))] m")
end

println("="^60)

@info "Simulation complete! Created $(length(captured_steps)) static plots."

nothing
