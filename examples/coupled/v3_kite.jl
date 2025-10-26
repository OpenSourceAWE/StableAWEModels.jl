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
using Statistics
using YAML
using GLMakie

@info "Loading 2-plate kite model..."

# ============= User settings =============
const MODEL_NAME = "v3"
const GEOM_PATH  = joinpath("data", MODEL_NAME, "struc_geometry.yaml")
const SIM_TIME   = 10.0
const N_STEPS    = 600
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
set.upwind_dir = -90.0  

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

# Create logger for recording simulation
using KiteUtils
logger = Logger(length(sam.sys_struct.points), n_steps + 1)
sys_state = SysState(sam)
sys_state.time = 0.0
log!(logger, sys_state)

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
fig = plot(sam.sys_struct)
display(fig)

# Time-marching loop with coupled aerodynamics
@info "Starting time-marching simulation with coupled aerodynamics..."
@info "Enforcing Y-symmetry constraint (setting Y-velocities to zero)"
@info "Number of wings in system: $(length(sam.sys_struct.wings))"

for step in 1:n_steps
    t = step * Δt

    try
        next_step!(sam; dt=Δt, vsm_interval=1)
    catch err
        @error "next_step! failed" step integrator_t=sam.integrator.t integrator_dt=sam.integrator.dt integrator_norm=norm(sam.integrator.u) integrator_max=maximum(abs, sam.integrator.u) exception=(err, catch_backtrace())
        rethrow(err)
    end

    # Log current state
    update_sys_state!(sys_state, sam)
    sys_state.time = t
    log!(logger, sys_state)

    # Print progress periodically
    if step % max(1, div(n_steps, 10)) == 0 || step == n_steps
        @info "  Step $step/$n_steps (t = $(round(t, digits=2)) s / $(SIM_TIME) s)"
    end
end

@info "Simulation complete. Creating interactive replay viewer with $(length(logger)) frames..."

# Create interactive replay viewer with slider controls using the logged data
save_log(logger, "tmp_run")
syslog = load_log("tmp_run")
scene = replay(syslog, sam.sys_struct; autoplay=false, loop=true)

@info "Replay viewer created! Use the slider to step through time, or press Play to replay."

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

@info "Simulation complete! Interactive replay viewer with $(length(logger)) frames ready."

nothing
