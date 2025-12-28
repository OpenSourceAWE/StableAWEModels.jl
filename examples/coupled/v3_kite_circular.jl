# Copyright (c) 2025 Jelle Poland, Bart van de Lint
# SPDX-License-Identifier: MPL-2.0

"""
TU Delft V3 Kite: REFINE Wing

This example runs a REFINE wing model and plots the results.
"""

using SymbolicAWEModels
using VortexStepMethod
using LinearAlgebra
using Statistics
using GLMakie
using KiteUtils
using DiscretePIDs
using Dates

# Heading PID controller parameters
MAX_HEADING = 10.0  # Maximum heading amplitude [degrees]
PERIOD = 60.0       # Oscillation period [seconds]
HEADING_P = 0.0     # Proportional gain
HEADING_I = 0.1     # Integral gain
HEADING_D = 0.0     # Derivative gain

# Tether winch PID controller parameters
WINCH_P = 1000.0    # Proportional gain [N/m]
WINCH_I = 100.0     # Integral gain [N/(m·s)]
WINCH_D = 50.0      # Derivative gain [N·s/m]

"""
    run_v3_kite(; kwargs...)

Run a v3 kite simulation using the REFINE wing type.

# Keyword Arguments
- `sim_time::Float64=300.0`: Simulation duration [s]
- `fps::Int=1`: Frames per second for logging
- `remake_cache::Bool=false`: Force rebuild of cached model
- `initial_damping::Float64=10.0`: Initial world frame damping [N·s/m]
- `decay_time::Float64=2.0`: Time for damping decay [s]
- `max_steering::Float64=0.2`: Maximum steering line length change [m]
- `show_plots::Bool=false`: Display 3D plots during simulation
- `v_wind::Float64=15.4`: Wind speed [m/s]
- `upwind_dir::Float64=-90.0`: Wind direction [°]
- `max_heading::Float64=50.0`: Maximum heading amplitude for sine wave [°]
- `period::Float64=20.0`: Oscillation period for sine wave [s]
- `heading_p::Float64=0.0`: Proportional gain for heading controller
- `heading_i::Float64=0.1`: Integral gain for heading controller
- `heading_d::Float64=0.0`: Derivative gain for heading controller
- `winch_p::Float64=1000.0`: Proportional gain for winch controller [N/m]
- `winch_i::Float64=100.0`: Integral gain for winch controller [N/(m·s)]
- `winch_d::Float64=50.0`: Derivative gain for winch controller [N·s/m]

# Returns
- `SysLog`: The simulation log containing time history data
"""
function run_v3_kite(;
                     sim_time=300.0,
                     fps=1,
                     remake_cache=false,
                     initial_damping=100.0,
                     decay_time=2.0,
                     up = 0.5858,
                     us=0.1,
                     show_plots=false,
                     v_wind=15.4,
                     upwind_dir=-90.0,
                     steering_ramp_time=25.0,
                     max_heading=50.0,
                     period=20.0,
                     heading_p=0.0,
                     heading_i=0.1,
                     heading_d=0.0,
                     winch_p=1000.0,
                     winch_i=100.0,
                     winch_d=50.0)

    wing_type = SymbolicAWEModels.REFINE
    wing_type_str = "REFINE"
    @info "Running v3 kite simulation with REFINE wing type..."

    # Load settings
    set_data_path("data/v3")
    set = Settings("system.yaml")
    set.v_wind = v_wind
    set.upwind_dir = upwind_dir

    # Load YAML structure path
    model_name = "v3_refine"
    struc_yaml_path = joinpath("data", "v3", "struc_geometry.yaml")

    # Load VSMSettings
    vsm_set_path = joinpath(get_data_path(), "vsm_settings.yaml")
    vsm_set = VortexStepMethod.VSMSettings(vsm_set_path; data_prefix=false)

    # Use 36 panels for both wing types (matches vsm_settings.yaml default)
    vsm_set.wings[1].n_panels = 36
    # Note: n_unrefined_sections is automatically inferred from YAML geometry

    # Load system structure with wing_type and vsm_set parameters
    sys = load_sys_struct_from_yaml(struc_yaml_path;
        system_name=model_name, set, wing_type, vsm_set)


    # Initialize damping
    SymbolicAWEModels.set_world_frame_damping(sys, initial_damping)

    wing_points = [p for p in sys.points if p.type == WING]
    n_unrefined = sys.wings[1].vsm_wing.n_unrefined_sections
    @info "REFINE wing setup:" n_wing_points=length(wing_points) n_groups=length(sys.groups) n_unrefined=n_unrefined n_panels=length(sys.wings[1].vsm_aero.panels) n_segments=length(sys.segments)

    # Create symbolic model
    sam = SymbolicAWEModel(set, sys)

    # Apply steering
    # sys.segments[87].l0 += max_steering
    # sys.segments[89].l0 -= max_steering

    # Initialize model
    n_unrefined = sys.wings[1].vsm_wing.n_unrefined_sections
    SymbolicAWEModels.init!(sam; remake=remake_cache, ignore_l0=false, remake_vsm=true)
    # hold tether length (no dynamic reel-in)
    sam.sys_struct.winches[1].brake = true

    # Create logger
    n_steps = Int(round(fps * sim_time))
    Δt = sim_time / n_steps
    logger = Logger(sam, n_steps + 1)
    sys_state = SysState(sam)
    sys_state.time = 0.0
    log!(logger, sys_state)

    # Store nominal segment lengths for PID control
    nominal_l0_87 = sys.segments[87].l0
    nominal_l0_88 = sys.segments[88].l0
    nominal_l0_89 = sys.segments[89].l0

    # Storage for heading setpoint
    heading_setpoint = Float64[]
    push!(heading_setpoint, 0.0)  # Fixed steering: keep heading setpoint at zero
    # # Create PID controller for heading
    # max_heading_rad = deg2rad(max_heading)
    # angular_freq = 2π / period  # rad/s
    # heading_pid = DiscretePID(;
    #     K = heading_p > 0 ? heading_p : 1.0,
    #     Ti = heading_i > 0 ? 1.0 / heading_i : false,
    #     Td = heading_d > 0 ? heading_d : false,
    #     Ts = Δt,
    #     umin = -abs(max_steering),
    #     umax = abs(max_steering)
    # )
    # @info "Heading PID controller initialized" max_heading period heading_p heading_i heading_d
    # @info "  Sine wave: ±$(max_heading)°, period=$(period)s"

    # Keep winch torque at zero; brake is engaged to fix tether length
    winch = sys.winches[1]
    winch.set_value = 0.0

    # Optional initial plot
    if show_plots
        scene = plot(sam.sys_struct)
        display(scene)
    end

    ## ACTUATION
    #len_power-tape = 200 + 5000 * up
    #len_steering_tape = 1600 + 1400 * us
    #us < 0 shortens right tape, and lengths left tape, causing a right turn
    #RIGHT len_steering_tape (us = -1) = 1600 - 1400 = 200
    #LEFT len_steering_tape (us = 1) = 1600 + 1400 = 3000

    steering_tape_change = 1400 * us / 1000  # Convert mm to m
    power_tape_change = ((200 + 5000 * up) / 1000) - nominal_l0_88  # Convert mm to m


    # Time-marching loop
    @info "Starting simulation: $n_steps steps, Δt = $(round(Δt, digits=4)) s"
    @info " Initial lengths (m): segment 87: $(round(nominal_l0_87, digits=4)), segment 88: $(round(nominal_l0_88, digits=4)), segment 89: $(round(nominal_l0_89, digits=4))"
    @info " Steering tape change (m): $(round(steering_tape_change, digits=4)), Power tape change (m): $(round(power_tape_change, digits=4))"
    sim_start_time = time()
    aoa_log_interval_steps = max(1, Int(round(3.0 / Δt)))  # roughly every 3 seconds

    wings = sam.sys_struct.wings
    wing = wings[1]

    for step in 1:n_steps
        t = step * Δt

        # Update damping
        if t <= decay_time
            current_damping = initial_damping * (1.0 - t / decay_time)
            SymbolicAWEModels.set_world_frame_damping(sam.sys_struct, current_damping)
        else
            SymbolicAWEModels.set_world_frame_damping(sam.sys_struct, 0.0)
        end

        # # PID heading control with sine wave setpoint
        # target_heading_rad = max_heading_rad * sin(angular_freq * t)
        # current_heading = sam.sys_struct.wings[1].heading
        # steering_control = heading_pid(target_heading_rad, current_heading, 0.0)
        # push!(heading_setpoint, target_heading_rad)

        # Fixed tether length: brake engaged; only steering ramp is applied
        ramp_factor = min(t / steering_ramp_time, 1.0)
        steering_control = steering_tape_change * ramp_factor
        power_control = power_tape_change * ramp_factor
        push!(heading_setpoint, 0.0)  # Keep heading setpoint flat for plotting

        # Apply power control
        # sys.segments[88].l0 = nominal_l0_88 + power_control

        # Apply differential steering (opposite signs for turning moment)
        sys.segments[87].l0 = nominal_l0_87 + steering_control
        sys.segments[89].l0 = nominal_l0_89 - steering_control


        # Convert force to torque: τ = -r/G * F + friction
        winch_torque = 0.0
        sys.winches[1].set_value = -winch_torque

        # Advance simulation
        try
            next_step!(sam; set_values=[-winch_torque], dt=Δt, vsm_interval=1)
        catch err
            if err isa AssertionError
                @error "next_step! failed" step
                break
            end
            rethrow(err)
        end

        # Log state
        update_sys_state!(sys_state, sam)
        sys_state.time = t
        log!(logger, sys_state)

        # Log AoA periodically (value stored in sys_state by update_sys_state!)
        if step % aoa_log_interval_steps == 0
            alpha = sys_state.AoA
            alpha = atan(wing.va_b[3], wing.va_b[1])
            vsm_alpha = wing.vsm_solver.sol.alpha_dist[length(wing.vsm_solver.sol.alpha_dist) ÷ 2 + (length(wing.vsm_solver.sol.alpha_dist) % 2)]
            @info "---> Angle of attack" t=round(t, digits=2) sys_state.AoA=round(rad2deg(alpha), digits=2) vsm_alpha=round(rad2deg(vsm_alpha), digits=2)

        end

        # Progress updates
        if step % max(1, div(n_steps, 10)) == 0 || step == n_steps
            elapsed = time() - sim_start_time
            times_realtime = t / elapsed
            @info "  Step $step/$n_steps (t = $(round(t, digits=2)) s)" times_realtime=round(times_realtime, digits=2)
        end
    end

    # Calculate performance
    total_wall_time = time() - sim_start_time
    final_times_realtime = sim_time / total_wall_time
    @info "Simulation completed: $wing_type_str" wall_time=round(total_wall_time, digits=2) times_realtime=round(final_times_realtime, digits=2)

    # Save and load log
    log_name = "tmp_run_$(lowercase(wing_type_str))"
    save_log(logger, log_name)
    syslog = load_log(log_name)

    # Permanent save
    save_dir = joinpath("processed_data", "v3_kite")
    isdir(save_dir) || mkpath(save_dir)
    timestamp = Dates.format(Dates.now(), "yyyy_mm_dd_HH_MM")
    up_tag = Int(round(up*100))
    us_tag = Int(round(us*100))
    log_name = "up_$(up_tag)" * "_" * "us_$(us_tag)" * "_" * "vw_$(v_wind)" * "_" * timestamp 
    save_log(logger, log_name; path=save_dir)


    return syslog, sam, heading_setpoint
end

# ============= Main Execution =============

@info "V3 Kite: Running REFINE simulation..."
us = 0.2 # in URIs kite as a sensor it goes up to about 0.3
up = 0.5858
vw = 15.0
sim_time = 10.0

syslog_refine, sam_refine, heading_setpoint_refine = run_v3_kite(
    sim_time=sim_time, fps=60, show_plots=false, up=up, us=us, v_wind=vw,
    max_heading=MAX_HEADING, period=PERIOD,
    heading_p=HEADING_P, heading_i=HEADING_I, heading_d=HEADING_D,
    winch_p=WINCH_P, winch_i=WINCH_I, winch_d=WINCH_D)

@info "Simulation complete. Creating plots..."

fig = plot(sam_refine.sys_struct, syslog_refine;
           plot_turn_rates=true,
           plot_reelout=false,
        #    plot_tether=true,
        #    plot_aero_force=true,
        #    plot_aero_moment=true,
        #    plot_tether_moment=true,
        #    plot_twist=true,
           plot_aoa=true,
           plot_heading=false,
        #    plot_old_heading=true,
        #    plot_distance=true,
        #    plot_cone_angle=true,
           plot_elevation=true,
           plot_azimuth=true,
           plot_winch_force=false,
           plot_set_values=false)
display(fig)

@info "Plot created!"

nothing
