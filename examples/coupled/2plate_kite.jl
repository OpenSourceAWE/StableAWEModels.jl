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
MODEL_NAME = "2plate_kite"
SIM_TIME   = 10.0
N_STEPS    = 600
REMAKE_CACHE = false
RAMP_TIME = 2.0            # Time to ramp inputs from 0 to magnitude [s]
STEERING_MAGNITUDE = 0.1   # Final steering line length offset [m] (differential)
DEPOWER_MAGNITUDE = -0.0    # Final depower line length offset [m] (both lines shorten)
# =========================================

# Set data path to the 2plate_kite project folder
pkg_root = dirname(dirname(@__DIR__))
set_data_path(joinpath(pkg_root, "data", MODEL_NAME))

# Load settings and VSM settings
set = Settings("system.yaml")
vsm_set_path = joinpath(get_data_path(), "vsm_settings.yaml")
vsm_set = VortexStepMethod.VSMSettings(vsm_set_path)

# Load system structure directly from YAML
@info "Creating 2-plate kite system structure..."
struc_yaml = joinpath(get_data_path(), "refine_struc_geometry.yaml")
sys = SymbolicAWEModels.load_sys_struct_from_yaml(struc_yaml; system_name=MODEL_NAME, set=set, vsm_set=vsm_set)
sys.winches[1].brake = false

@info "System loaded with $(length(sys.points)) points, $(length(sys.segments)) segments, $(length(sys.wings)) wings"

# Create symbolic model with the 2-plate system
@info "Creating SymbolicAWEModel with aerodynamics..."
sys.winches[:main_winch].brake = true
sam = SymbolicAWEModel(set, sys)

# Store baseline steering line lengths for ramping
l0_left_base = sam.sys_struct.segments[:kcu_steering_left].l0
l0_right_base = sam.sys_struct.segments[:kcu_steering_right].l0

# Initialize the model
@info "Initializing model..."
SymbolicAWEModels.init!(sam; remake=REMAKE_CACHE, lin_vsm=false)

# Calculate simulation parameters
n_steps = N_STEPS
Δt = SIM_TIME / n_steps
@info "Running simulation for $(SIM_TIME) seconds ($n_steps steps, Δt = $(round(Δt, digits=4)) s)..."

# Create logger for recording simulation
using KiteUtils
logger = Logger(sam, n_steps + 1)
sys_state = SysState(sam)
sys_state.time = 0.0
log!(logger, sys_state)

# Time-marching loop with coupled aerodynamics
@info "Starting time-marching simulation with coupled aerodynamics..."

for step in 1:n_steps
    t = step * Δt

    # Ramp steering and depower from 0 to magnitude over RAMP_TIME
    ramp = clamp(t / RAMP_TIME, 0.0, 1.0)
    steering = STEERING_MAGNITUDE * ramp
    depower = DEPOWER_MAGNITUDE * ramp
    sam.sys_struct.segments[:kcu_steering_left].l0 = l0_left_base - steering + depower
    sam.sys_struct.segments[:kcu_steering_right].l0 = l0_right_base + steering + depower

    try
        next_step!(sam; dt=Δt, vsm_interval=1)
    catch err
        if err isa AssertionError
            @error "next_step! failed"
            break
        else
            rethrow(err)
        end
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
