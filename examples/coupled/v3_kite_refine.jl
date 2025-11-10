# Copyright (c) 2025 Jelle Poland, Bart van de Lint
# SPDX-License-Identifier: MPL-2.0

"""
TU Delft V3 Kite with REFINE Wing Type

This example demonstrates the REFINE wing type, which applies VSM panel forces
directly to structural points instead of using quaternion-based rigid body dynamics.

Key differences from standard v3_kite.jl:
- Wing points are DYNAMIC and deform under loads
- VSM panel forces are lumped to structural points via inverse distance weighting
- No group twist dynamics (REFINE wings cannot have groups)
- Two-way coupling: structure deforms → VSM panels update → forces update
- NO linearization (forces come directly from nonlinear VSM solve)
"""

using SymbolicAWEModels
using VortexStepMethod
using LinearAlgebra
using Statistics
using YAML
using GLMakie

@info "Loading v3 kite model with REFINE wing type..."

# ============= User settings =============
const MODEL_NAME = "v3"
const SIM_TIME   = 4.0
const FPS        = 200
const N_STEPS    = Int(round(FPS * SIM_TIME))
const REMAKE_CACHE = false
const INITIAL_DAMPING = 10.0  # Initial world frame damping [N·s/m]
const DECAY_TIME = 1.0        # Time for damping to decay to zero [s]
# =========================================

# Load settings for the v3 kite
set_data_path("data/v3")
set = Settings("system.yaml")
if hasproperty(set, :c_spring) || hasproperty(set, :damping)
    @info "Legacy tether settings still present" c_spring=(hasproperty(set, :c_spring) ? getproperty(set, :c_spring) : missing) damping=(hasproperty(set, :damping) ? getproperty(set, :damping) : missing)
end
@info "Tether parameters (post-load)" axial_stiffness=set.axial_stiffness axial_damping=set.axial_damping d_tether=set.d_tether cd_tether=set.cd_tether

# Load v3 system structure directly from YAML (automatically creates REFINE wing)
@info "Loading v3 kite system structure from YAML..."
model_name = hasproperty(set, :model_name) ? set.model_name : MODEL_NAME
struc_yaml = hasproperty(set, :struc_geometry_path) ? set.struc_geometry_path :
    joinpath("data", model_name, "struc_geometry.yaml")
sys = load_sys_struct_from_yaml(struc_yaml; system_name=model_name, set=set)

# Initialize damping to starting value
SymbolicAWEModels.set_world_frame_damping(sys, INITIAL_DAMPING)

# Verify REFINE wing setup
@assert length(sys.wings) > 0 "No wings in system"
@assert sys.wings[1].wing_type == SymbolicAWEModels.REFINE "Wing should be REFINE type"
@assert length(sys.groups) == 0 "REFINE wings should have no groups"

wing_points = [p for p in sys.points if p.type == WING]
@info "REFINE wing setup:" n_wing_points=length(wing_points) n_panels=length(sys.wings[1].vsm_aero.panels) n_segments=length(sys.segments)

segment_props = [(idx=seg.idx, k=seg.axial_stiffness, c=seg.axial_damping, d=seg.diameter) for seg in sys.segments]
@info "Segment mechanical properties (first 10)" segment_props[1:min(10, length(segment_props))]

@info "System loaded with $(length(sys.points)) points, $(length(sys.segments)) segments, $(length(sys.wings)) wings"

# Create symbolic model with the v3 REFINE system
@info "Creating SymbolicAWEModel with REFINE aerodynamics..."
@info "Note: Linearization is disabled for REFINE wings (too many DOFs)"
sam = SymbolicAWEModel(set, sys)

# Test with head-on wind
set.upwind_dir = -90.0

# Calculate wind vector components in X,Y,Z
upwind_rad = deg2rad(set.upwind_dir)
# Check if wind_elevation exists in settings, default to 0.0
wind_elev = hasfield(typeof(set), :wind_elevation) ? set.wind_elevation : 0.0
wind_elev_rad = deg2rad(wind_elev)

# Print wind settings with vector components
@info "Wind Configuration:" wind_speed=set.v_wind upwind_dir=set.upwind_dir wind_elevation=wind_elev
@info "  → Wind direction: $(set.upwind_dir)° (0° = North, 90° = East, 180° = South, -90° = West)"
@info "  → Wind speed: $(set.v_wind) m/s at reference height $(set.h_ref) m"

# Initialize the model
# NOTE: REFINE wings do NOT use linearization (too expensive with many structural DOFs)
@info "Initializing model without VSM linearization..."
# sys.points[20].fix_sphere=true
# sys.points[2].fix_sphere=true
# sys.points[10].fix_sphere=true
# sys.points[12].fix_sphere=true
SymbolicAWEModels.init!(sam; remake=REMAKE_CACHE, ignore_l0=false)
wing = sam.sys_struct.wings[1]
vsm_aero = wing.vsm_aero
vsm_solver = wing.vsm_solver
vsm_wing = wing.vsm_wing
results = VortexStepMethod.solve(vsm_solver, vsm_aero; log=true)
body_y_coordinates = [panel.aero_center[2] for panel in vsm_aero.panels]
plot_distribution(
    [body_y_coordinates],
    [results],
    ["VSM"];
    title="CAD_spanwise_distributions",
    data_type=".pdf",
    is_save=false,
    is_show=true,
)

# Calculate simulation parameters
n_steps = N_STEPS
Δt = SIM_TIME / n_steps
@info "Running simulation for $(SIM_TIME) seconds ($n_steps steps, Δt = $(round(Δt, digits=4)) s)..."

[point.fix_static = true for point in sys.points if point.type == WING]
@time next_step!(sam; dt=10.0)
[point.fix_static = false for point in sys.points if point.type == WING]
# [point.fix_sphere = true for point in sys.points if point.type == WING]
# sys.points[1].fix_sphere=true
# @time for i in 1:n_steps
#     next_step!(sam; dt=Δt)
# end
# [point.fix_sphere = false for point in sys.points if point.type == WING]
# sys.points[1].fix_sphere=false

# Create logger for recording simulation
using KiteUtils
logger = Logger(length(sam.sys_struct.points), n_steps + 1)
sys_state = SysState(sam)
sys_state.time = 0.0
log!(logger, sys_state)

## Plot initial state (user must display if desired)
@info "Creating initial plot..."
scene = plot(sam.sys_struct)
display(scene)

# Time-marching loop with coupled aerodynamics
@info "Starting time-marching simulation with coupled REFINE aerodynamics..."
@info "  REFINE wing: Panel forces → lumped to $(length(wing_points)) structural points"
@info "  NO linearization: Forces computed directly from nonlinear VSM solve each timestep"
@info "  Two-way coupling: Structure deforms → VSM panels update → Forces update"

for step in 1:n_steps
    t = step * Δt

    # Update damping: linearly decay from INITIAL_DAMPING to 0 over DECAY_TIME
    if t <= DECAY_TIME
        current_damping = INITIAL_DAMPING * (1.0 - t / DECAY_TIME)
        SymbolicAWEModels.set_world_frame_damping(sam.sys_struct, current_damping)
    else
        SymbolicAWEModels.set_world_frame_damping(sam.sys_struct, 0.0)
    end

    try
        # Update VSM panels from deformed structure (if needed, done internally)
        # Then solve VSM and apply lumped forces to points
        next_step!(sam; dt=Δt, vsm_interval=1)
    catch err
        @error "next_step! failed" step integrator_t=sam.integrator.t integrator_dt=sam.integrator.dt integrator_norm=norm(sam.integrator.u) integrator_max=maximum(abs, sam.integrator.u) exception=(err, catch_backtrace())
        rethrow(err)
    end

    # Log current state
    update_sys_state!(sys_state, sam)
    sys_state.time = t
    log!(logger, sys_state)

    # Print progress periodically with REFINE-specific info
    if step % max(1, div(n_steps, 10)) == 0 || step == n_steps
        # Calculate average wing point position for tracking deformation
        avg_wing_pos = mean([p.pos_w for p in wing_points])
        current_damp = t <= DECAY_TIME ? INITIAL_DAMPING * (1.0 - t / DECAY_TIME) : 0.0
        @info "  Step $step/$n_steps (t = $(round(t, digits=2)) s)" avg_wing_z=round(avg_wing_pos[3], digits=2) damping=round(current_damp, digits=2)
    end
end

@info "Simulation complete. Creating interactive replay viewer with $(length(logger)) frames..."

# Create interactive replay viewer with slider controls using the logged data
save_log(logger, "tmp_run_refine")
syslog = load_log("tmp_run_refine")
scene = replay(syslog, sam.sys_struct; autoplay=false, loop=true)
display(scene)

@info "Replay viewer created! Use the slider to step through time, or press Play to replay."

# Print final statistics
println("\n" * "="^60)
println("Final Simulation Results (t = $(SIM_TIME) s) - REFINE Wing")
println("="^60)

# Calculate wing structural point statistics
if !isempty(wing_points)
    avg_x = mean([p.pos_w[1] for p in wing_points])
    avg_y = mean([p.pos_w[2] for p in wing_points])
    avg_z = mean([p.pos_w[3] for p in wing_points])
    println("  Average wing point position: [$(round(avg_x, digits=2)), $(round(avg_y, digits=2)), $(round(avg_z, digits=2))] m")

    # Calculate span of wing points
    y_coords = [p.pos_w[2] for p in wing_points]
    span = maximum(y_coords) - minimum(y_coords)
    println("  Wing span: $(round(span, digits=2)) m")

    # Calculate deformation from initial positions
    displacements = [norm(p.pos_w - p.pos_cad) for p in wing_points]
    avg_displacement = mean(displacements)
    max_displacement = maximum(displacements)
    println("  Average structural displacement: $(round(avg_displacement, digits=3)) m")
    println("  Maximum structural displacement: $(round(max_displacement, digits=3)) m")
end

# Calculate average position of dynamic (non-wing) points
dynamic_points = filter(p -> p.type == SymbolicAWEModels.DYNAMIC && p.type != SymbolicAWEModels.WING, sam.sys_struct.points)
if !isempty(dynamic_points)
    avg_x = mean([p.pos_w[1] for p in dynamic_points])
    avg_y = mean([p.pos_w[2] for p in dynamic_points])
    avg_z = mean([p.pos_w[3] for p in dynamic_points])
    println("  Average dynamic point position: [$(round(avg_x, digits=2)), $(round(avg_y, digits=2)), $(round(avg_z, digits=2))] m")
end

println("\n  REFINE wing type: Direct panel forces applied to $(length(wing_points)) structural points")
println("  VSM panels: $(length(sys.wings[1].vsm_aero.panels)) panels with forces lumped via inverse distance weighting")
println("  NO linearization: Forces computed directly from nonlinear VSM solve")
println("="^60)

@info "Simulation complete! Interactive replay viewer with $(length(logger)) frames ready."
@info "Note: Compare with standard v3_kite.jl to see the difference between REFINE and QUATERNION wing types"

nothing
