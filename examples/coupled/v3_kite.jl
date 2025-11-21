# Copyright (c) 2025 Jelle Poland, Bart van de Lint
# SPDX-License-Identifier: MPL-2.0

"""
TU Delft V3 Kite Comparison: REFINE vs QUATERNION Wing Types

This example compares REFINE and QUATERNION wing models:
- REFINE: Panel forces applied to structural points, deformable wing
- QUATERNION: Rigid body dynamics with group twist DOFs

The script runs both simulations and plots them for direct comparison.
"""

using SymbolicAWEModels
using VortexStepMethod
using LinearAlgebra
using Statistics
using GLMakie
using KiteUtils

"""
    run_v3_kite(wing_type::WingType; kwargs...)

Run a v3 kite simulation with the specified wing type.

# Arguments
- `wing_type::WingType`: Either `REFINE` or `QUATERNION`

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

# Returns
- `SysLog`: The simulation log containing time history data
"""
function run_v3_kite(wing_type::WingType;
                     sim_time=300.0,
                     fps=1,
                     remake_cache=false,
                     initial_damping=10.0,
                     decay_time=2.0,
                     max_steering=0.1,
                     show_plots=false,
                     v_wind=15.4,
                     upwind_dir=-90.0)

    wing_type_str = wing_type == SymbolicAWEModels.REFINE ? "REFINE" : "QUATERNION"
    @info "Running v3 kite simulation with $wing_type_str wing type..."

    # Load settings
    set_data_path("data/v3")
    set = Settings("system.yaml")
    set.v_wind = v_wind
    set.upwind_dir = upwind_dir

    # Load YAML structure path
    model_name = wing_type == QUATERNION ? "v3_quat" : "v3_refine"
    struc_yaml_path = joinpath("data", "v3", "struc_geometry.yaml")

    # Load VSMSettings based on wing type
    vsm_set_path = joinpath(get_data_path(), "vsm_settings.yaml")
    vsm_set = VortexStepMethod.VSMSettings(vsm_set_path; data_prefix=false)

    # Modify VSM settings based on wing type
    if wing_type == SymbolicAWEModels.REFINE
        # REFINE uses 36 panels and 9 groups
        vsm_set.wings[1].n_panels = 36
        vsm_set.wings[1].n_groups = 9
    else
        # QUATERNION uses 40 panels and 10 groups
        vsm_set.wings[1].n_panels = 40
        vsm_set.wings[1].n_groups = 10
    end

    # Load system structure with wing_type and vsm_set parameters
    sys = load_sys_struct_from_yaml(struc_yaml_path;
        system_name=model_name, set, wing_type, vsm_set)


    # Initialize damping
    SymbolicAWEModels.set_world_frame_damping(sys, initial_damping)

    wing_points = [p for p in sys.points if p.type == WING]
    @info "$wing_type_str wing setup:" n_wing_points=length(wing_points) n_groups=length(sys.groups) n_panels=length(sys.wings[1].vsm_aero.panels) n_segments=length(sys.segments)

    # Create symbolic model
    sam = SymbolicAWEModel(set, sys)

    # Apply steering
    sys.segments[87].l0 += max_steering
    sys.segments[89].l0 -= max_steering

    # Initialize model
    @show sys.wings[1].vsm_wing.n_groups
    SymbolicAWEModels.init!(sam; remake=remake_cache, ignore_l0=false, remake_vsm=true)

    # Stabilization phase
    @info "Stabilizing system..."
    [point.fix_static = true for point in sys.points if point.type == WING]
    if wing_type == QUATERNION
        sys.wings[1].fix_sphere = true
    end
    @time next_step!(sam; dt=10.0)
    [point.fix_static = false for point in sys.points if point.type == WING]
    if wing_type == QUATERNION
        sys.wings[1].fix_sphere = false
    end

    # Create logger
    n_steps = Int(round(fps * sim_time))
    Δt = sim_time / n_steps
    logger = Logger(sam, n_steps + 1)
    sys_state = SysState(sam)
    sys_state.time = 0.0
    log!(logger, sys_state)

    # Optional initial plot
    if show_plots
        scene = plot(sam.sys_struct)
        display(scene)
    end

    # Time-marching loop
    @info "Starting simulation: $n_steps steps, Δt = $(round(Δt, digits=4)) s"
    sim_start_time = time()

    for step in 1:n_steps
        t = step * Δt

        # Update damping
        if t <= decay_time
            current_damping = initial_damping * (1.0 - t / decay_time)
            SymbolicAWEModels.set_world_frame_damping(sam.sys_struct, current_damping)
        else
            SymbolicAWEModels.set_world_frame_damping(sam.sys_struct, 0.0)
        end

        # Advance simulation
        try
            next_step!(sam; dt=Δt, vsm_interval=1)
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

    return syslog, sam
end

# ============= Main Execution =============

@info "V3 Kite Comparison: Running REFINE and QUATERNION simulations..."

# Run both simulations
# syslog_refine, sam_refine = run_v3_kite(SymbolicAWEModels.REFINE; sim_time=60.0, fps=60, show_plots=false)
syslog_quat, sam_quat = run_v3_kite(QUATERNION; sim_time=60.0, fps=24, show_plots=false)

@info "Both simulations complete. Creating comparison plots..."

# # Create comparison plot
# fig = plot(sam_refine.sys_struct, [syslog_refine, syslog_quat];
#            plot_turn_rates=true,
#            plot_azimuth=true,
#            plot_elevation=true,
#            plot_aoa=true,
#            plot_heading=true,
#            plot_default=false,
#            plot_aero_force=true)

# display(fig)

# @info "Comparison plot created!"

nothing
