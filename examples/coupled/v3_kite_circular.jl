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
using CairoMakie
GLMakie.activate!()
using KiteUtils
using DiscretePIDs
using Dates

include("utils.jl")

# Geometry configuration (must match settle_refine_wing.jl output)
TETHER_LENGTH = 248          # Total tether length (m)
TE_FRAC = 0.95               # Factor for TE wires (segments 20-28), 1.0 = no change
TIP_REDUCTION = 0.4          # Tip LE reduction (m), 0.0 = no change
GEOM_SUFFIX = build_geom_suffix(V3_DEPOWER_L0, TIP_REDUCTION, TE_FRAC)
STRUC_YAML_PATH = "data/v3/struc_geometry_$(GEOM_SUFFIX).yaml"
AERO_YAML_PATH = "data/v3/aero_geometry_$(GEOM_SUFFIX).yaml"

"""
    run_v3_kite(; kwargs...)

Run a v3 kite simulation using the REFINE wing type.
Uses exponential decay world frame damping (WORLD_DAMPING, DECAY_STEPS constants).

# Keyword Arguments
- `sim_time::Float64=300.0`: Simulation duration [s]
- `fps::Int=4`: Frames per second for logging
- `remake_cache::Bool=false`: Force rebuild of cached model
- `damping_pattern::Vector{Float64}=[0.0, 30.0, 60.0]`: Body frame damping [x, y, z] (N·s/m)
- `up::Float64=40.0`: Depower percentage [0, 100]
- `us::Float64=10.0`: Steering percentage [-100, 100], positive = right turn
- `show_plots::Bool=false`: Display 3D plots during simulation
- `v_wind::Float64=15.4`: Wind speed [m/s]
- `ramp_time::Float64=25.0`: Time to ramp steering/depower [s]

# Returns
- `SysLog`: The simulation log containing time history data
"""
function run_v3_kite(;
                     sim_time=300.0,
                     fps=4,
                     remake_cache=false,
                     damping_pattern=[0.0, 30.0, 60.0],
                     up=40.0,
                     us=10.0,
                     show_plots=false,
                     v_wind=15.4,
                     upwind_dir=-90.0,
                     ramp_time=25.0,
                     v_wind_base=15)

    wing_type = SymbolicAWEModels.REFINE
    wing_type_str = "REFINE"
    @info "Running v3 kite simulation" n_steps=Int(round(fps * sim_time)) STRUC_YAML_PATH AERO_YAML_PATH

    # Load settings
    set_data_path("data/v3")
    set = Settings("system.yaml")
    set.v_wind = v_wind
    set.upwind_dir = upwind_dir
    set.l_tether = TETHER_LENGTH
    set.v_reel_outs[1] = 0.0

    # Load YAML structure path
    model_name = "v3"

    # Load VSMSettings
    vsm_set_path = joinpath(get_data_path(), "vsm_settings_reduced_for_coupling.yaml")
    vsm_set = VortexStepMethod.VSMSettings(vsm_set_path; data_prefix=false)
    vsm_set.wings[1].geometry_file = AERO_YAML_PATH

    # Use 36 panels for both wing types (matches vsm_settings.yaml default)
    vsm_set.wings[1].n_panels = 36
    # Note: n_unrefined_sections is automatically inferred from YAML geometry

    # Load system structure with wing_type and vsm_set parameters
    sys = load_sys_struct_from_yaml(STRUC_YAML_PATH;
        system_name=model_name, set, wing_type, vsm_set)

    sys.transforms[1].elevation = deg2rad(20.0)
    sys.transforms[1].azimuth = deg2rad(20.0)
    # Initialize damping with per-axis values [x, y, z]
    SymbolicAWEModels.set_body_frame_damping(sys, damping_pattern, 1:38)

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
    # SymbolicAWEModels.reinit!(sam, sam.prob, SymbolicAWEModels.FBDF())

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

    # Optional initial plot
    if show_plots
        scene = plot(sam.sys_struct)
        display(scene)
    end

    # Calculate target tape lengths from percentages using shared functions
    L_left_target, L_right_target = steering_percentage_to_lengths(us)
    L_depower_target = depower_percentage_to_length(up)
    vw_change = v_wind - v_wind_base


    # Time-marching loop
    @info "Starting simulation: $n_steps steps, Δt = $(round(Δt, digits=4)) s"
    @info " Initial lengths (m): segment 87: $(round(nominal_l0_87, digits=4)), segment 88: $(round(nominal_l0_88, digits=4)), segment 89: $(round(nominal_l0_89, digits=4))"
    @info " Target lengths (m): L_left: $(round(L_left_target, digits=4)), L_right: $(round(L_right_target, digits=4)), L_depower: $(round(L_depower_target, digits=4))"
    sim_start_time = time()
    aoa_log_interval_steps = max(1, Int(round(3.0 / Δt)))  # roughly every 3 seconds

    # Storage for segment stretch statistics (after t > 1.0)
    max_stretch_samples = Float64[]
    mean_stretch_samples = Float64[]
    max_idx_samples = Int[]

    # Storage for tape lengths plotting
    tape_times = Float64[]
    tape_steering_pct = Float64[]
    tape_depower_pct = Float64[]

    wings = sam.sys_struct.wings
    wing = wings[1]

    for step in 1:n_steps
        t = step * Δt

        # Body frame damping stays constant at damping_pattern

        # Ramp from initial to target tape lengths
        ramp_factor = min(t / ramp_time, 1.0)
        push!(heading_setpoint, 0.0)  # Keep heading setpoint flat for plotting

        # Apply ramped steering, instant depower
        sys.segments[87].l0 = nominal_l0_87 + ramp_factor * (L_left_target - nominal_l0_87)
        sys.segments[89].l0 = nominal_l0_89 + ramp_factor * (L_right_target - nominal_l0_89)
        sys.segments[88].l0 = L_depower_target  # Instant depower

        # Log tape lengths as percentages
        push!(tape_times, t)
        push!(tape_steering_pct, ramp_factor * us)
        push!(tape_depower_pct, up)  # Constant depower

        # Update wind speed linearly
        sam.sys_struct.set.v_wind = v_wind_base + vw_change

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

        # Plot body frame at t=10.0 (right turn frame)
        if abs(t - 10.0) < Δt / 2
            @info "Plotting body frame at t=10.0 (right turn)..."

            # Load right turn frame points from CSV
            extra_pts_csv = "data/v3/right_turn_reelout_frame_7362.csv"
            extra_pts, extra_groups = load_extra_points(extra_pts_csv, sam.sys_struct)

            # Save PDFs for all three views
            CairoMakie.activate!()
            for dir in (:front, :side, :top)
                scene = plot_body_frame_local([sam.sys_struct];
                    extra_points=extra_pts, extra_groups=extra_groups,
                    dir, point_idxs=1:38, labels=["sim"])
                pdf_filename = "data/v3/circular_body_frame_$(dir)_right_turn_$(GEOM_SUFFIX).pdf"
                save(pdf_filename, scene)
                @info "Plot saved" pdf_filename
            end

            # Display interactive plot
            GLMakie.activate!()
            scene = plot_body_frame_local([sam.sys_struct];
                extra_points=extra_pts, extra_groups=extra_groups,
                dir=:side, point_idxs=1:38, labels=["sim"])
            scr = display(scene)
            @info "Close the plot window to continue..."
            wait(scr)
        end

        # Collect segment stretch statistics after t > 1.0
        if t > 1.0
            max_stretch, mean_stretch, max_idx = segment_stretch_stats(sam.sys_struct)
            push!(max_stretch_samples, max_stretch)
            push!(mean_stretch_samples, mean_stretch)
            push!(max_idx_samples, max_idx)
        end

        # Log AoA periodically (value stored in sys_state by update_sys_state!)
        # if step % aoa_log_interval_steps == 0
        #     alpha = sys_state.AoA
        #     alpha = atan(wing.va_b[3], wing.va_b[1])
        #     vsm_alpha = wing.vsm_solver.sol.alpha_dist[length(wing.vsm_solver.sol.alpha_dist) ÷ 2 + (length(wing.vsm_solver.sol.alpha_dist) % 2)]
        #     @info "---> Angle of attack" t=round(t, digits=2) sys_state.AoA=round(rad2deg(alpha), digits=2) vsm_alpha=round(rad2deg(vsm_alpha), digits=2)

        # end

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

    # Report segment stretch statistics (for t > 1.0)
    if !isempty(max_stretch_samples)
        overall_max_stretch = maximum(max_stretch_samples)
        overall_mean_stretch = mean(mean_stretch_samples)
        max_stretch_idx = max_idx_samples[argmax(max_stretch_samples)]
        @info "Segment stretch statistics (t > 1.0):" max_relative=round(overall_max_stretch, digits=6) max_percentage=round(overall_max_stretch*100, digits=4) mean_relative=round(overall_mean_stretch, digits=6) mean_percentage=round(overall_mean_stretch*100, digits=4) max_segment_idx=max_stretch_idx
    end

    # Save and load log
    log_name = "tmp_run_$(lowercase(wing_type_str))"
    save_log(logger, log_name)
    syslog = load_log(log_name)

    # Permanent save
    save_dir = joinpath("processed_data", "v3_kite")
    isdir(save_dir) || mkpath(save_dir)
    timestamp = Dates.format(Dates.now(), "yyyy_mm_dd_HH_MM")
    up_tag = Int(round(up))
    us_tag = Int(round(us))
    v_wind_tag = Int(round(v_wind))
    log_name = "circle__up_$(up_tag)" * "_" * "us_$(us_tag)" * "_" * "vw_$(v_wind_tag)" * "_date_" * timestamp
    save_log(logger, log_name; path=save_dir)

    # Create tape data for plotting
    tape_data = (
        time = tape_times,
        steering = tape_steering_pct,
        depower = tape_depower_pct
    )

    return syslog, sam, heading_setpoint, tape_data
end

# ==========================================
# ============= Main Execution =============
# ==========================================
us = 20.0   # Steering percentage [-100, 100], positive = right turn
up = 39   # Depower percentage [0, 100]
vw = 8.0    # Wind speed [m/s]

sim_time = 40.0
ramp_time = 2.0
fps = 60
damping_pattern = [0.0, 60.0, 120.0]


syslog, sam, heading_setpoint, tape_data = run_v3_kite(
    sim_time=sim_time, fps=fps,
    up=up, us=us, v_wind=vw,
    ramp_time=ramp_time, damping_pattern=damping_pattern,
)


fig = plot(sam.sys_struct, syslog;
    plot_turn_rates=true, plot_reelout=false, plot_gk=true,
    plot_aoa=true, plot_heading=false, plot_elevation=true,
    plot_azimuth=true, plot_winch_force=false, plot_set_values=false,
    tape_lengths=[tape_data])

scene = replay(syslog, sam.sys_struct)

scr1 = display(fig)
wait(scr1)
scr2 = display(scene)
wait(scr2)



# Report final geometric AoA using hardcoded mid-panel corners (world frame)
last_state = syslog.syslog[end]
X = last_state.X; Y = last_state.Y; Z = last_state.Z
# Mid-panel corners: 10,11,12,13 (11/13 front; 10/12 back)
back = 0.5 .* ([X[10], Y[10], Z[10]] .+ [X[12], Y[12], Z[12]])
front = 0.5 .* ([X[11], Y[11], Z[11]] .+ [X[13], Y[13], Z[13]])

delta_z = front[3] - back[3]
delta_x = front[1] - back[1]
aoa_wrt_horizontal = -rad2deg(atan(delta_z, delta_x))
# @info "alpha wrt horizontal $(round(aoa_wrt_horizontal, digits=2))"

chord_w = front .- back
wing = sam.sys_struct.wings[1]
v_app_w = wing.R_b_w * wing.va_b
@info "v_app" v_app_w=round.(v_app_w, digits=2)
aoa_geom_deg = rad2deg(acos(clamp(dot(chord_w, v_app_w) / (norm(chord_w) * norm(v_app_w) + 1e-12), -1.0, 1.0)))
@info "alpha wrt v_app $(round(aoa_geom_deg, digits=2))"

nothing
