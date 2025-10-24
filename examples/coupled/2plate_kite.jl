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
using GLMakie

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
        if any(isnan, y_state[1:6])
            @warn "Skipping VSM force refresh for wing $(wing.idx) due to NaNs in state" y_state
            continue
        end
        va = [y_state[1], y_state[2], y_state[3]]
        omega = [y_state[4], y_state[5], y_state[6]]

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
const SIM_TIME   = 0.5  # Total simulation time in seconds
const N_STEPS    = 10     # Number of time steps (use small number for fast development)
const REMAKE_CACHE = false
# =========================================

# Include helper utilities first
yaml_loader_path = joinpath(@__DIR__, "..", "yaml_loader.jl")
if isfile(yaml_loader_path)
    include(yaml_loader_path)  # provides load_sys_struct_from_yaml
else
    @warn "yaml_loader.jl not found; expecting load_sys_struct_from_yaml to be available"
end

# Load settings for the 2-plate kite
set = SymbolicAWEModels.load_settings(MODEL_NAME)
if hasproperty(set, :c_spring) || hasproperty(set, :damping)
    @info "Legacy tether settings still present" c_spring=(hasproperty(set, :c_spring) ? getproperty(set, :c_spring) : missing) damping=(hasproperty(set, :damping) ? getproperty(set, :damping) : missing)
end
@info "Tether parameters (post-load)" axial_stiffness=set.axial_stiffness axial_damping=set.axial_damping d_tether=set.d_tether cd_tether=set.cd_tether

# The SystemStructure factory will now call create_2plate_sys_struct()
# which loads both structural (struc_geometry.yaml) and aerodynamic (aero_geometry.yaml) data
@info "Creating 2-plate kite system structure..."
sys = SymbolicAWEModels.SystemStructure(set)
segment_props = [(idx=seg.idx, k=seg.axial_stiffness, c=seg.axial_damping, d=seg.diameter) for seg in sys.segments]
@info "Segment mechanical properties" segment_props

@info "System loaded with $(length(sys.points)) points, $(length(sys.segments)) segments, $(length(sys.wings)) wings"

# Create symbolic model with the 2-plate system
@info "Creating SymbolicAWEModel with aerodynamics..."
sam = SymbolicAWEModel(set, sys)

# Test with head-on wind
set.upwind_dir = 0.0  

# Calculate wind vector components in X,Y,Z
upwind_rad = deg2rad(set.upwind_dir)
# Check if wind_elevation exists in settings, default to 0.0
wind_elev = hasfield(typeof(set), :wind_elevation) ? set.wind_elevation : 0.0
wind_elev_rad = deg2rad(wind_elev)
# Wind convention: 0° = North (+X), -90° = East (+Y), 90° = West (-Y), 180° = South (-X)
wind_x = set.v_wind * cos(wind_elev_rad) * cos(upwind_rad)
wind_y = set.v_wind * cos(wind_elev_rad) * sin(upwind_rad)
wind_z = set.v_wind * sin(wind_elev_rad)

# Print wind settings with vector components
@info "Wind Configuration:" wind_speed=set.v_wind upwind_dir=set.upwind_dir wind_elevation=wind_elev
@info "  → Wind direction: $(set.upwind_dir)° (0° = North/head-on, -90° = East, 90° = West, 180° = South)"
@info "  → Wind speed: $(set.v_wind) m/s at reference height $(set.h_ref) m"
@info "  → Wind vector: [$(round(wind_x, digits=3)), $(round(wind_y, digits=3)), $(round(wind_z, digits=3))] m/s (X, Y, Z components)"

# Initialize the model
@info "Initializing model..."
SymbolicAWEModels.init!(sam; remake=REMAKE_CACHE, lin_vsm=false)

# Calculate simulation parameters
n_steps = N_STEPS
Δt = SIM_TIME / n_steps
@info "Running simulation for $(SIM_TIME) seconds ($n_steps steps, Δt = $(round(Δt, digits=4)) s)..."

# Store snapshots for every step
snapshots = Dict{Int, Vector{SymbolicAWEModels.Point}}(0 => deepcopy(sam.sys_struct.points))

# Print initial node coordinates
@info "Initial state (Step 0) node coordinates:"
for (i, point) in enumerate(sam.sys_struct.points)
    pos = point.pos_w
    vel = point.vel_w
    @info "  Node $i: pos=[$(round(pos[1], digits=3)), $(round(pos[2], digits=3)), $(round(pos[3], digits=3))] m, " *
          "vel=[$(round(vel[1], digits=3)), $(round(vel[2], digits=3)), $(round(vel[3], digits=3))] m/s"
end

## Plot initial state
@info "Plotting initial state..."
try
    fig = plot(sam.sys_struct)
    display(fig)
catch e
    @warn "Could not plot initial state: $e"
end

# Time-marching loop with coupled aerodynamics
@info "Starting time-marching simulation with coupled aerodynamics..."
@info "Enforcing Y-symmetry constraint (setting Y-velocities to zero)"
@info "Number of wings in system: $(length(sam.sys_struct.wings))"

for step in 1:n_steps
    t = step * Δt

    refresh_vsm_forces!(sam)

    # Print VSM aerodynamic loads on wings
    if !isempty(sam.sys_struct.wings)
        @info "Step $step VSM aerodynamic loads:"
        for wing in sam.sys_struct.wings
            force = wing.vsm_x[1:3]
            moment = wing.vsm_x[4:6]
            @info "  Wing $(wing.idx): " *
                  "Force=[$(round(force[1], digits=3)), $(round(force[2], digits=3)), $(round(force[3], digits=3))] N, " *
                  "Moment=[$(round(moment[1], digits=3)), $(round(moment[2], digits=3)), $(round(moment[3], digits=3))] N⋅m"
        end
    else
        if step == 1
            @warn "No wings found in system structure - VSM loads not available"
        end
    end

    # Advance simulation one step using freshly computed VSM loads (linear refresh disabled)
    integrator_t = sam.integrator.t
    integrator_dt = sam.integrator.dt
    integrator_norm = norm(sam.integrator.u)
    integrator_max = maximum(abs, sam.integrator.u)
    @info "Integrator state before step" step integrator_t integrator_dt integrator_norm integrator_max
    try
        next_step!(sam; dt=Δt, vsm_interval=0)
    catch err
        @error "next_step! failed" step integrator_t=sam.integrator.t integrator_dt=sam.integrator.dt integrator_norm=norm(sam.integrator.u) integrator_max=maximum(abs, sam.integrator.u) exception=(err, catch_backtrace())
        rethrow(err)
    end

    # SYMMETRY CONSTRAINT: Force Y-velocities to zero to maintain X-Z plane symmetry
    for point in sam.sys_struct.points
        if point.type == SymbolicAWEModels.DYNAMIC
            point.vel_w[2] = 0.0  # Zero out Y-velocity
        end
    end

    # Store current state for every step
    snapshots[step] = deepcopy(sam.sys_struct.points)

    # Print node coordinates for this step
    @info "Step $step node coordinates:"
    for (i, point) in enumerate(sam.sys_struct.points)
        pos = point.pos_w
        vel = point.vel_w
        @info "  Node $i: pos=[$(round(pos[1], digits=3)), $(round(pos[2], digits=3)), $(round(pos[3], digits=3))] m, " *
              "vel=[$(round(vel[1], digits=3)), $(round(vel[2], digits=3)), $(round(vel[3], digits=3))] m/s"
    end

    # Print progress periodically
    if step % max(1, div(n_steps, 10)) == 0 || step == n_steps
        @info "  Step $step/$n_steps (t = $(round(t, digits=2)) s / $(SIM_TIME) s)"
    end
end

captured_steps = sort!(collect(keys(snapshots)))

@info "Simulation complete. Creating interactive animation with $(length(captured_steps)) frames..."

# Create interactive animation with slider controls
# Using lock_limits=false to let the camera auto-fit as the kite moves
bbox = ((-20.0, 20.0), (-20.0, 20.0), (0.0, 40.0))
fig = animate(sam.sys_struct, snapshots; dt=Δt, autoplay=false, loop=true, bbox=bbox)
display(fig)

@info "Animation created! Use the slider to step through time, or press Play to animate."

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

@info "Simulation complete! Interactive animation with $(length(captured_steps)) frames ready."

nothing
