#!/usr/bin/env julia
# Copyright (c) 2025 Bart van de Lint
# SPDX-License-Identifier: MPL-2.0

"""
Settle REFINE wing with world frame damping

This script loads the V3 kite with REFINE wing type, runs it with constant
world frame damping while maintaining alignment through repositioning, and
updates the YAML point positions to the settled equilibrium state.

The script creates a backup of the original YAML file before modifying it.
"""

using SymbolicAWEModels
using VortexStepMethod
using LinearAlgebra
using KiteUtils
using OrdinaryDiffEqCore
using GLMakie

# Configuration
const BODY_DAMPING = 100.0  # Ns/m
const NUM_STEPS = 200
const DT = 0.05  # seconds
const STRUC_PATH = "data/v3/struc_geometry.yaml"
const AERO_PATH = "data/v3/aero_geometry.yaml"
const STEERING_PERCENTAGE = 0.0  # steering [-100, 100]
const DEPOWER_PERCENTAGE = 20.0   # depower [0, 100]

# V3 Kite steering/depower calibration (from KCU documentation)
const STEERING_L0 = 1.506  # Neutral steering tape length (m)
const STEERING_GAIN = 1.2  # Maximum differential (m) at |u_s| = 1
const DEPOWER_L0 = 0.45
const DEPOWER_GAIN = 5.0

"""
    steering_percentage_to_lengths(percentage)

Convert steering percentage to left/right tape lengths (m).
Percentage convention: negative = left turn, positive = right turn.
"""
function steering_percentage_to_lengths(percentage)
    u_s = percentage / 100.0  # Convert percentage to [-1, 1]
    L_left = STEERING_L0 - (STEERING_GAIN / 2.0) * u_s
    L_right = STEERING_L0 + (STEERING_GAIN / 2.0) * u_s
    return L_left, L_right
end

"""
    depower_percentage_to_length(percentage)

Convert depower percentage to tape length (m).
"""
function depower_percentage_to_length(percentage)
    u_p = percentage / 100.0  # Convert percentage to [0, 1]
    L_depower = DEPOWER_L0 + DEPOWER_GAIN * u_p
    return L_depower
end

@info "Settling REFINE wing with world frame damping..."
@info "Configuration" BODY_DAMPING NUM_STEPS DT total_time=NUM_STEPS*DT

# Load settings
set_data_path("data/v3")
set = Settings("system.yaml")

# Load VSMSettings
vsm_set_path = joinpath(get_data_path(), "vsm_settings.yaml")
vsm_set = VortexStepMethod.VSMSettings(vsm_set_path; data_prefix=false)
vsm_set.wings[1].n_panels = 36

# Load system structure with REFINE wing type
struc_yaml_path = joinpath("data", "v3", "struc_geometry.yaml")
sys = load_sys_struct_from_yaml(struc_yaml_path;
    system_name="v3", set,
    wing_type=SymbolicAWEModels.REFINE, vsm_set)

# Set constant world frame damping
SymbolicAWEModels.set_body_frame_damping(sys, BODY_DAMPING)

wing_points = [p for p in sys.points if p.type == WING]
@info "System setup" n_wing_points=length(wing_points) n_points=length(sys.points) n_segments=length(sys.segments)

# Create symbolic model
sam = SymbolicAWEModel(set, sys)

# Initialize model
SymbolicAWEModels.init!(sam; remake=false, ignore_l0=false, remake_vsm=true)

# Enable winch brakes to lock tether length
for winch in sys.winches
    winch.brake = true
end
@info "Winch brakes enabled"

# Initialize steering and depower from configuration
L_left, L_right = steering_percentage_to_lengths(STEERING_PERCENTAGE)
L_depower = depower_percentage_to_length(DEPOWER_PERCENTAGE)

sys.segments[87].l0 = L_left   # Left steering tape
sys.segments[89].l0 = L_right  # Right steering tape
sys.segments[88].l0 = L_depower  # Depower tape

@info "Steering/depower initialized" steering=STEERING_PERCENTAGE depower=DEPOWER_PERCENTAGE L_left L_right L_depower

# Get initial transform values to maintain
transform = sys.transforms[1]
target_elevation = transform.elevation
target_azimuth = transform.azimuth
initial_heading = sys.wings[1].heading

@info "Initial state" elevation=rad2deg(target_elevation) azimuth=rad2deg(target_azimuth) heading=rad2deg(initial_heading)

# Create logger
logger = Logger(sam, NUM_STEPS + 1)
sys_state = SysState(sam)
sys_state.time = 0.0
log!(logger, sys_state)

# Simulation loop with repositioning every step
@info "Starting settling simulation..."
for step in 1:NUM_STEPS
    t = step * DT

    # Calculate delta heading needed to maintain alignment
    current_heading = sys.wings[1].heading
    delta_heading = initial_heading - current_heading

    # Update transform to maintain orientation
    transform.elevation = target_elevation
    transform.azimuth = target_azimuth
    transform.heading = delta_heading

    # Apply repositioning to maintain alignment
    SymbolicAWEModels.reposition!(sys.transforms, sys)

    # Reinitialize the solver after repositioning
    OrdinaryDiffEqCore.reinit!(sam.integrator)

    # Advance one timestep
    try
        next_step!(sam; dt=DT, vsm_interval=1)
    catch err
        if err isa AssertionError
            @error "Simulation failed" step t
            break
        end
        rethrow(err)
    end

    # Log state
    update_sys_state!(sys_state, sam)
    sys_state.time = t
    log!(logger, sys_state)

    # Progress updates
    if step % 20 == 0 || step == NUM_STEPS
        wing = sys.wings[1]
        @info "Step $step/$NUM_STEPS (t = $(round(t, digits=2)) s)" elevation=round(rad2deg(wing.elevation), digits=2) azimuth=round(rad2deg(wing.azimuth), digits=2) heading=round(rad2deg(wing.heading), digits=2)
    end
end

@info "Simulation completed. Updating YAML with settled positions..."

# Update YAML files with settled positions
update_yaml_from_sys_struct!(sys, STRUC_PATH, AERO_PATH)

# Save and load the log
log_name = "settle_refine_wing"
save_log(logger, log_name)
syslog = load_log(log_name)

@info "Creating plots..."

# Create plots using the Makie extension
fig = plot(sam.sys_struct, syslog;
           plot_tether=true,
           plot_aero_force=true)
display(fig)

@info "Plot created! Settling complete."
