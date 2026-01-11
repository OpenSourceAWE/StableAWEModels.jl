# Copyright (c) 2025 Jelle Poland, Bart van de Lint
# SPDX-License-Identifier: MPL-2.0

"""
TU Delft V3 Kite: Batch run for zenith hold followed by circular flight

This script runs multiple parameter combinations for the two-phase v3 kite simulation.
"""

using SymbolicAWEModels
using VortexStepMethod
using LinearAlgebra
using Statistics
using GLMakie
using KiteUtils
using DiscretePIDs
using Dates
using StaticArrays

"""
    adjust_tether_length!(sam::SymbolicAWEModel, tether_length_raw; tether_point_idxs=39:44)

Update the winch rest length, reposition tether points in CAD/body frames,
and reapply the main transform so the wing stays at the requested tether radius.
"""
function adjust_tether_length!(sam::SymbolicAWEModel, tether_length_raw; tether_point_idxs=39:44)
    tether_length = float(tether_length_raw)
    sys = sam.sys_struct
    set = sam.set

    if !isempty(set.l_tethers)
        set.l_tethers[1] = tether_length
    end

    n_points = length(tether_point_idxs)
    for (n, p_idx) in enumerate(tether_point_idxs)
        pos = (0.0, 0.0, -n * tether_length / n_points)
        sys.points[p_idx].pos_cad .= pos
        sys.points[p_idx].pos_b .= pos
    end

    if !isempty(sys.transforms)
        transform = sys.transforms[1]
        if !isempty(sys.wings) && norm(sys.wings[1].pos_w) > 0
            target_pos = sys.wings[1].pos_w
            transform.elevation = KiteUtils.calc_elevation(target_pos)
            transform.azimuth = KiteUtils.azimuth_east(target_pos)
        end
        SymbolicAWEModels.reinit!([transform], sys)
    end

    if !isempty(sys.winches)
        winch = sys.winches[1]
        winch.tether_len = tether_length
        winch.tether_vel = 0.0
        winch.brake = true
    end
    return nothing
end

"""
    adjust_elevation!(sam::SymbolicAWEModel, elevation_deg)

Update the transform elevation to the specified value in degrees.
"""
function adjust_elevation!(sam::SymbolicAWEModel, elevation_deg)
    sys = sam.sys_struct
    
    if !isempty(sys.transforms)
        transform = sys.transforms[1]
        transform.elevation = deg2rad(elevation_deg)
        SymbolicAWEModels.reinit!([transform], sys)
    end
    return nothing
end

"""
    run_v3_kite(wing_type::WingType; kwargs...)

Run a two-phase v3 kite simulation with the specified wing type
(zenith azimuth hold, then circular flight). The two phases can have
independent durations/FPS.

# Arguments
- `wing_type::WingType`: Either `REFINE` or `QUATERNION`

# Keyword Arguments
- `sim_time_zenith::Float64=300.0`: Duration of zenith phase [s]
- `fps_zenith::Int=1`: Frames per second for zenith logging
- `sim_time_circles::Float64=300.0`: Duration of circular flight phase [s]
- `fps_circles::Int=1`: Frames per second for circular flight logging
- `remake_cache::Bool=false`: Force rebuild of cached model
- `initial_damping::Float64=10.0`: Initial body frame damping [N·s/m]
- `damping_pattern::Vector{Float64}=[0.0, 1.0, 1.0]`: Per-axis damping pattern [x, y, z]
- `decay_time::Float64=2.0`: Time for damping decay [s] (applied in both phases)
- `min_damping::Float64=0.0`: Minimum damping after decay [N·s/m]
- `up::Float64=0.4`: Power-tape input (0–1) mapped to segment 88 length
- `ramp_time_up::Float64=25.0`: Power-tape ramp time during zenith phase [s]
- `ramp_time_us::Float64=25.0`: Steering/power ramp time during circular phase [s]
- `start_ramp_time::Float64=0.0`: Time offset before starting power-tape ramp [s] (zenith phase)
- `max_us_zenith::Float64=0.1`: Max steering tape change for azimuth PID [m/1.4]
- `us::Float64=0.1`: Steering tape command for the circular phase
- `show_plots::Bool=false`: Display 3D plots during simulation
- `v_wind::Float64=15.4`: Wind speed [m/s]
- `v_wind_base::Float64=15.0`: Baseline wind used for circular phase ramp [m/s]
- `upwind_dir::Float64=-90.0`: Wind direction [°]
- `heading_p::Float64=0.0`: Proportional gain for azimuth controller
- `heading_i::Float64=0.1`: Integral gain for azimuth controller
- `heading_d::Float64=0.0`: Derivative gain for azimuth controller
- `target_azimuth::Float64=0.0`: Azimuth setpoint during zenith phase [rad]
- `winch_p::Float64=1000.0`: Proportional gain for winch controller [N/m]
- `winch_i::Float64=100.0`: Integral gain for winch controller [N/(m·s)]
- `winch_d::Float64=50.0`: Derivative gain for winch controller [N·s/m]
- `tether_length::Float64=150.0`: Tether length [m]
- `elevation::Union{Nothing,Float64}=nothing`: Initial elevation angle [°] (overrides YAML if provided)
- `g_earth::Union{Nothing,Float64}=nothing`: Gravitational acceleration [m/s²] (overrides YAML if provided)
- `save_subdir::AbstractString=""`: Subfolder under `processed_data/v3_kite` for permanent saves
- `run_tag::AbstractString=""`: Extra tag appended to the log name

# Returns
- `SysLog`: The simulation log containing time history data
"""
function run_v3_kite(wing_type::WingType;
                     sim_time_zenith=300.0,
                     fps_zenith=1,
                     sim_time_circles=300.0,
                     fps_circles=1,
                     remake_cache=false,
                     initial_damping=100.0,
                     damping_pattern=[0.0, 1.0, 1.0],
                     decay_time=2.0,
                     min_damping=1.0,
                     up=0.4,
                     ramp_time_up=25.0,
                     ramp_time_us=25.0,
                     start_ramp_time=0.0,
                     max_us_zenith=0.1,
                     show_plots=false,
                     v_wind=15.4,
                     upwind_dir=-90.0,
                     heading_p=0.0,
                     heading_i=0.1,
                     heading_d=0.0,
                     winch_p=1000.0,
                     winch_i=100.0,
                     winch_d=50.0,
                     REFINE = true,
                     target_azimuth=0.0,
                     us=0.1,
                     v_wind_base=15.0,
                     tether_length=150.0,
                     elevation=nothing,
                     g_earth=nothing,
                     save_subdir="",
                     run_tag="")

    wing_type_str = wing_type == SymbolicAWEModels.REFINE ? "REFINE" : "QUATERNION"
    @info "Running v3 kite simulation with $wing_type_str wing type (zenith -> circular)..."

    # Load settings
    set_data_path("data/v3")
    set = Settings("system.yaml")
    set.v_wind = v_wind
    set.upwind_dir = upwind_dir
    if g_earth !== nothing
        set.g_earth = g_earth
    end

    # Load YAML structure path
    model_name = wing_type == QUATERNION ? "v3_quat" : "v3_refine"
    struc_yaml_path = joinpath("data", "v3", "CORRECT_struc_geometry.yaml")

    # Load VSMSettings
    vsm_set_path = joinpath(get_data_path(), "CORRECT_vsm_settings.yaml")
    vsm_set = VortexStepMethod.VSMSettings(vsm_set_path; data_prefix=false)

    # Use 36 panels for both wing types (matches vsm_settings.yaml default)
    vsm_set.wings[1].n_panels = 36
    # Note: n_unrefined_sections is automatically inferred from YAML geometry

    # Load system structure with wing_type and vsm_set parameters
    sys = load_sys_struct_from_yaml(struc_yaml_path;
        system_name=model_name, set, wing_type, vsm_set)

    # Initialize damping with per-axis values [x, y, z]
    SymbolicAWEModels.set_body_frame_damping(sys, damping_pattern * initial_damping)

    wing_points = [p for p in sys.points if p.type == WING]
    n_unrefined = sys.wings[1].vsm_wing.n_unrefined_sections
    @info "$wing_type_str wing setup:" n_wing_points=length(wing_points) n_groups=length(sys.groups) n_unrefined=n_unrefined n_panels=length(sys.wings[1].vsm_aero.panels) n_segments=length(sys.segments)

    # Create symbolic model
    sam = SymbolicAWEModel(set, sys)
    adjust_tether_length!(sam, tether_length)
    
    # Adjust elevation if provided
    if elevation !== nothing
        adjust_elevation!(sam, elevation)
    end

    # Initialize model
    n_unrefined = sys.wings[1].vsm_wing.n_unrefined_sections
    SymbolicAWEModels.init!(sam; remake=remake_cache, ignore_l0=false, remake_vsm=true)

    # Create logger (two phases with independent settings)
    n_steps_zenith = max(1, Int(round(fps_zenith * sim_time_zenith)))
    Δt_zenith = sim_time_zenith / n_steps_zenith
    n_steps_circles = max(1, Int(round(fps_circles * sim_time_circles)))
    Δt_circles = sim_time_circles / n_steps_circles
    total_steps = n_steps_zenith + n_steps_circles
    logger = Logger(sam, total_steps + 1)
    sys_state = SysState(sam)
    sys_state.time = 0.0
    log!(logger, sys_state)

    # Store nominal segment lengths for PID/control
    nominal_l0_87 = sys.segments[87].l0
    nominal_l0_88 = sys.segments[88].l0
    nominal_l0_89 = sys.segments[89].l0
    power_tape_change = ((200 + 5000 * up) / 1000) - nominal_l0_88  # Convert mm to m

    # Storage for azimuth setpoint (targeting azimuth=0)
    azimuth_setpoint = Float64[]
    push!(azimuth_setpoint, target_azimuth)

    # Create PID controller for azimuth (reuse heading gains/limits)
    max_steering = max_us_zenith * 1.4  # Max steering line length change [m]

    azimuth_pid = DiscretePID(;
        K = heading_p > 0 ? heading_p : 1.0,
        Ti = heading_i > 0 ? 1.0 / heading_i : false,
        Td = heading_d > 0 ? heading_d : false,
        Ts = Δt_zenith,
        umin = -abs(max_steering),
        umax = abs(max_steering)
    )
    @info "Azimuth PID controller initialized" target_azimuth target_azimuth_deg=rad2deg(target_azimuth) heading_p heading_i heading_d

    # Store nominal tether length for winch controller
    nominal_tether_length = sys.winches[1].tether_len

    # Initialize winch torque based on initial tether force
    winch = sys.winches[1]
    initial_force = norm(winch.force)
    initial_torque = -winch.drum_radius / winch.gear_ratio * initial_force + winch.friction
    winch.set_value = initial_torque
    @info "Winch initialized" initial_force initial_torque drum_radius=winch.drum_radius gear_ratio=winch.gear_ratio friction=winch.friction

    # Create PID controller for tether winch (maintain constant length)
    max_force = 50000.0  # N
    winch_pid = DiscretePID(;
        K = winch_p,
        Ti = winch_i > 0 ? winch_p / winch_i : false,
        Td = winch_d > 0 ? winch_d / winch_p : false,
        Ts = Δt_zenith,
        umin = -max_force,
        umax = max_force
    )
    @info "Winch PID controller initialized" nominal_tether_length winch_p winch_i winch_d

    # Optional initial plot
    if show_plots
        scene = plot(sam.sys_struct)
        display(scene)
    end

    # Phase 1: Zenith/azimuth hold
    @info "Starting zenith phase: $n_steps_zenith steps, Δt = $(round(Δt_zenith, digits=4)) s"
    @info "Power-tape ramp" up=round(up, digits=3) target_l0_88=round(nominal_l0_88 + power_tape_change, digits=4) ramp_time=ramp_time_up start_ramp_time=start_ramp_time
    sim_start_time = time()

    for step in 1:n_steps_zenith
        t = step * Δt_zenith

        # Update damping: decay to min_damping
        current_damping = max(initial_damping * (1.0 - t / decay_time), min_damping)
        SymbolicAWEModels.set_body_frame_damping(sam.sys_struct, damping_pattern * current_damping)

        # PID azimuth control: drive azimuth to target
        current_azimuth = sam.sys_struct.wings[1].azimuth
        steering_control = azimuth_pid(target_azimuth, current_azimuth, 0.0)

        push!(azimuth_setpoint, target_azimuth)

        # Ramp power-tape change on segment 88
        if t >= start_ramp_time
            elapsed_ramp = t - start_ramp_time
            ramp_factor = min(elapsed_ramp / ramp_time_up, 1.0)
            power_control = power_tape_change * ramp_factor
        else
            power_control = 0.0
        end
        sys.segments[88].l0 = nominal_l0_88 + power_control

        # Apply differential steering (opposite signs for turning moment)
        sys.segments[87].l0 = nominal_l0_87 + steering_control
        sys.segments[89].l0 = nominal_l0_89 - steering_control

        # Winch PID: maintain constant tether length
        current_tether_len = winch.tether_len
        force_correction = winch_pid(nominal_tether_length, current_tether_len, 0.0)
        winch_torque = -winch.drum_radius / winch.gear_ratio * force_correction + winch.friction
        winch.set_value = winch_torque

        # Advance simulation
        try
            next_step!(sam; set_values=[winch_torque], dt=Δt_zenith, vsm_interval=1)
        catch err
            @error "Simulation failed at zenith phase" t step error=err
            rethrow(err)
        end

        # Log state
        update_sys_state!(sys_state, sam)
        sys_state.time = t
        log!(logger, sys_state)
    end

    # Phase 2: Circular flight
    @info "Switching to circular flight phase" phase_time=round(sys_state.time, digits=2)
    SymbolicAWEModels.set_body_frame_damping(sam.sys_struct, damping_pattern * initial_damping)

    # Fix tether length for circular flight
    winch.brake = true
    winch.set_value = 0.0

    steering_tape_change = 1400 * us / 1000  # Convert mm to m
    vw_change = v_wind - v_wind_base

    power_target = nominal_l0_88 + power_tape_change
    steer_target_87 = nominal_l0_87 + steering_tape_change
    steer_target_89 = nominal_l0_89 - steering_tape_change
    power_start = sys.segments[88].l0
    steer_start_87 = sys.segments[87].l0
    steer_start_89 = sys.segments[89].l0

    for step in 1:n_steps_circles
        t = sim_time_zenith + step * Δt_circles
        phase_t = step * Δt_circles

        # Update damping: decay to min_damping
        current_damping = max(initial_damping * (1.0 - phase_t / decay_time), min_damping)
        SymbolicAWEModels.set_body_frame_damping(sam.sys_struct, damping_pattern * current_damping)

        # Ramp steering and power during circular phase
        ramp_factor = min(phase_t / ramp_time_us, 1.0)
        
        # Ramp power
        sys.segments[88].l0 = power_start + (power_target - power_start) * ramp_factor
        
        # Ramp steering
        sys.segments[87].l0 = steer_start_87 + (steer_target_87 - steer_start_87) * ramp_factor
        sys.segments[89].l0 = steer_start_89 + (steer_target_89 - steer_start_89) * ramp_factor

        # Ramp wind speed
        sam.sys_struct.set.v_wind = v_wind_base + vw_change * ramp_factor

        # Advance simulation
        try
            next_step!(sam; set_values=[0.0], dt=Δt_circles, vsm_interval=1)
        catch err
            @error "Simulation failed at circular phase" t step error=err
            rethrow(err)
        end

        # Log state
        update_sys_state!(sys_state, sam)
        sys_state.time = t
        log!(logger, sys_state)
    end

    # Calculate performance
    total_wall_time = time() - sim_start_time
    total_sim_time = sim_time_zenith + sim_time_circles
    final_times_realtime = total_sim_time / total_wall_time
    @info "Simulation completed: $wing_type_str (zenith + circular)" wall_time=round(total_wall_time, digits=2) times_realtime=round(final_times_realtime, digits=2)

    lt_tag = Int(round(tether_length))

    # Save and load log
    log_name = "tmp_run_$(lowercase(wing_type_str))_lt_$(lt_tag)"
    save_log(logger, log_name)
    syslog = load_log(log_name)

    # Permanent save
    save_root = joinpath("processed_data", "v3_kite")
    save_dir = isempty(save_subdir) ? save_root : joinpath(save_root, save_subdir)
    isdir(save_dir) || mkpath(save_dir)
    timestamp = Dates.format(Dates.now(), "yyyy_mm_dd_HH_MM_SS")
    up_tag = Int(round(up*100))
    us_tag = Int(round(us*100))
    v_wind_tag = Int(round(v_wind))
    elev_tag = elevation !== nothing ? Int(round(elevation)) : "yaml"
    g_tag = g_earth !== nothing ? Int(round(g_earth*10)) : "yaml"
    log_name = "zenith_circle__up_$(up_tag)" * "_" * "us_$(us_tag)" * "_" * "vw_$(v_wind_tag)" * "_" * "lt_$(lt_tag)" * "_" * "el_$(elev_tag)" * "_" * "g_$(g_tag)"
    if !isempty(run_tag)
        log_name *= "_" * run_tag
    end
    log_name *= "_date_" * timestamp
    save_log(logger, log_name; path=save_dir)

    return syslog, sam, azimuth_setpoint
end

# ==========================================
# =============== Batch Run ================
# ==========================================

# Batch sweep configuration
elevation_vals  = [20,25,30,35,45,50,55,60,65,70,75,80,85]  # Initial elevation angles [°]
g_earth_vals    = [0.0]  # Gravitational acceleration [m/s²]
us_vals         = [0.0]  # Steering inputs
up_vals         = [0.18]#[0.42]  # Power inputs
vw_vals         = [8.6, 19.8]#[7.6]  # Wind speeds [m/s]
lt_vals         = [268] #[262]  # Tether lengths [m]

batch_tag = "zenith_2019_batch_" * Dates.format(Dates.now(), "yyyy_mm_dd_HH_MM_SS")

# Simulation settings
sim_time_zenith = 500.0
sim_time_circles = 0.0
fps_zenith = 60
fps_circles = 60
decay_time = 10.0
start_ramp_time = 1.0
ramp_time_up = 10.0
ramp_time_us = 5.0
initial_damping = 100.0
damping_pattern = [0.0, 1.0, 1.0]
min_damping = 1.0
max_us_zenith = 0.02

failed_runs = NamedTuple[]

for (run_id, (elev, g, us, up, vw, lt)) in enumerate(Iterators.product(elevation_vals, g_earth_vals, us_vals, up_vals, vw_vals, lt_vals))
    run_tag = "run_" * lpad(string(run_id), 3, '0')
    @info "Starting run" run_id elevation=elev g_earth=g us up vw lt batch_tag
    try
        syslog, sam, azimuth_setpoint = run_v3_kite(SymbolicAWEModels.REFINE;
            # General settings
            v_wind=vw,
            v_wind_base=vw,
            up=up,
            tether_length=lt,
            elevation=elev,
            g_earth=g,
            # Zenith initialization flight
            sim_time_zenith=sim_time_zenith,
            fps_zenith=fps_zenith,
            start_ramp_time=start_ramp_time,
            ramp_time_up=ramp_time_up,
            initial_damping=initial_damping,
            damping_pattern=damping_pattern,
            min_damping=min_damping,
            decay_time=decay_time,
            max_us_zenith=max_us_zenith,
            target_azimuth=0.0,
            # Circular flight
            sim_time_circles=sim_time_circles,
            fps_circles=fps_circles,
            ramp_time_us=ramp_time_us,
            us=us,
            # Batch settings
            save_subdir=batch_tag,
            run_tag=run_tag
        )
        @info "Completed run" run_id elevation=elev g_earth=g us up vw lt
    catch err
        @error "Failed run" run_id elevation=elev g_earth=g us up vw lt error=err
        push!(failed_runs, (run_id=run_id, elevation=elev, g_earth=g, us=us, up=up, vw=vw, lt=lt, error=err))
    end
    GC.gc()
end

if !isempty(failed_runs)
    fail_path = joinpath("processed_data", "v3_kite", batch_tag, "failed_runs.txt")
    open(fail_path, "w") do io
        for (i, fr) in enumerate(failed_runs)
            println(io, "Run $(fr.run_id): elevation=$(fr.elevation), g_earth=$(fr.g_earth), us=$(fr.us), up=$(fr.up), vw=$(fr.vw), lt=$(fr.lt)")
            println(io, "  Error: $(fr.error)")
        end
    end
    @info "Wrote failure list" path=fail_path
end

@info "Batch run completed" total_runs=length(collect(Iterators.product(elevation_vals, g_earth_vals, us_vals, up_vals, vw_vals, lt_vals))) failed=length(failed_runs) batch_tag


# ==========================================
# =============== Batch Run ================
# ==========================================

# Batch sweep configuration
elevation_vals  = [20,25,30,35,45,50,55,60,65,70,75,80,85]  # Initial elevation angles [°]
g_earth_vals    = [0.0]  # Gravitational acceleration [m/s²]
us_vals         = [0.0]  # Steering inputs
up_vals         = [0.42]  # Power inputs
vw_vals         = [7.8 ,19.7]  # Wind speeds [m/s]
lt_vals         = [262]  # Tether lengths [m]

batch_tag = "zenith_2025_batch_" * Dates.format(Dates.now(), "yyyy_mm_dd_HH_MM_SS")

# Simulation settings
sim_time_zenith = 200.0
sim_time_circles = 0.0
fps_zenith = 60
fps_circles = 60
decay_time = 10.0
start_ramp_time = 1.0
ramp_time_up = 10.0
ramp_time_us = 5.0
initial_damping = 100.0
damping_pattern = [0.0, 1.0, 1.0]
min_damping = 1.0
max_us_zenith = 0.02

failed_runs = NamedTuple[]

for (run_id, (elev, g, us, up, vw, lt)) in enumerate(Iterators.product(elevation_vals, g_earth_vals, us_vals, up_vals, vw_vals, lt_vals))
    run_tag = "run_" * lpad(string(run_id), 3, '0')
    @info "Starting run" run_id elevation=elev g_earth=g us up vw lt batch_tag
    try
        syslog, sam, azimuth_setpoint = run_v3_kite(SymbolicAWEModels.REFINE;
            # General settings
            v_wind=vw,
            v_wind_base=vw,
            up=up,
            tether_length=lt,
            elevation=elev,
            g_earth=g,
            # Zenith initialization flight
            sim_time_zenith=sim_time_zenith,
            fps_zenith=fps_zenith,
            start_ramp_time=start_ramp_time,
            ramp_time_up=ramp_time_up,
            initial_damping=initial_damping,
            damping_pattern=damping_pattern,
            min_damping=min_damping,
            decay_time=decay_time,
            max_us_zenith=max_us_zenith,
            target_azimuth=0.0,
            # Circular flight
            sim_time_circles=sim_time_circles,
            fps_circles=fps_circles,
            ramp_time_us=ramp_time_us,
            us=us,
            # Batch settings
            save_subdir=batch_tag,
            run_tag=run_tag
        )
        @info "Completed run" run_id elevation=elev g_earth=g us up vw lt
    catch err
        @error "Failed run" run_id elevation=elev g_earth=g us up vw lt error=err
        push!(failed_runs, (run_id=run_id, elevation=elev, g_earth=g, us=us, up=up, vw=vw, lt=lt, error=err))
    end
    GC.gc()
end

if !isempty(failed_runs)
    fail_path = joinpath("processed_data", "v3_kite", batch_tag, "failed_runs.txt")
    open(fail_path, "w") do io
        for (i, fr) in enumerate(failed_runs)
            println(io, "Run $(fr.run_id): elevation=$(fr.elevation), g_earth=$(fr.g_earth), us=$(fr.us), up=$(fr.up), vw=$(fr.vw), lt=$(fr.lt)")
            println(io, "  Error: $(fr.error)")
        end
    end
    @info "Wrote failure list" path=fail_path
end

@info "Batch run completed" total_runs=length(collect(Iterators.product(elevation_vals, g_earth_vals, us_vals, up_vals, vw_vals, lt_vals))) failed=length(failed_runs) batch_tag
