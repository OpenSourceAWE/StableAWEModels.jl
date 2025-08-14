# Copyright (c) 2025 Bart van de Lint
# SPDX-License-Identifier: MPL-2.0

using Pkg
if ! ("ControlPlots" ∈ keys(Pkg.project().dependencies))
    using TestEnv; TestEnv.activate()
end
using SymbolicAWEModels, KiteUtils, LinearAlgebra, ControlPlots, OrdinaryDiffEqBDF

"""
    sim_reposition!(sam; dt, total_time, reposition_interval_s, target_elevation_deg,
                    target_azimuth_deg, prn)

Run a simulation that periodically resets the kite's elevation and azimuth.

This function simulates the AWE model and, at a specified time interval, calls
`reposition!` to reposition the kite to a target elevation and azimuth. It logs
the entire simulation and returns a `SysLog`.

# Arguments
- `sam::SymbolicAWEModel`: Initialized AWE model.

# Keywords
- `dt::Float64`: Time step [s]. Defaults to `1/sam.set.sample_freq`.
- `total_time::Float64`: Total simulation duration [s]. Defaults to 20.0.
- `reposition_interval_s::Float64`: The interval in seconds at which to reset the pose. Defaults to 5.0.
- `target_elevation_deg::Float64`: The target elevation in degrees for repositioning. Defaults to 45.0.
- `target_azimuth_deg::Float64`: The target azimuth in degrees for repositioning. Defaults to 0.0.
- `prn::Bool`: If true, prints status messages during the simulation. Defaults to true.

# Returns
- `SysLog`: A log of the complete simulation.
"""
function sim_reposition!(
    sam::SymbolicAWEModel;
    dt=1/sam.set.sample_freq,
    total_time=20.0,
    reposition_interval_s=5.0,
    target_elevation_deg=45.0,
    target_azimuth_deg=0.0,
    prn=true
)
    # 1. --- Initialization ---
    sys_struct = sam.sys_struct
    steps = Int(round(total_time / dt))
    reposition_interval_steps = Int(round(reposition_interval_s / dt))
    set_values = zeros(Float64, steps, length(sys_struct.winches))
    vsm_interval = 1 ÷ dt
    
    logger = Logger(length(sys_struct.points), steps)
    sys_state = SysState(sam)

    if prn
        println("--- Starting simulation with periodic repositioning ---")
        println("Total time: $(total_time)s, Reposition interval: $(reposition_interval_s)s")
    end

    # 2. --- Simulation Loop ---
    time = @elapsed for step in 1:steps
        t = (step-1) * dt
        
        # Hold the kite in place by countering the tether forces with winch torques
        set_values[step, :] = -sam.set.drum_radius .* [norm(winch.force) for winch in sys_struct.winches]
        
        # --- Repositioning Logic ---
        if step > 1 && (step - 1) % reposition_interval_steps == 0
            if prn
                println("\n>>> Time: $(round(t, digits=2))s. Repositioning kite...")
                println(">>> Target Elevation: $(target_elevation_deg)°, Target Azimuth: $(target_azimuth_deg)°")
            end
            
            # Update the transform with the new target pose
            sys_struct.transforms[1].elevation = deg2rad(target_elevation_deg)
            sys_struct.transforms[1].azimuth   = deg2rad(target_azimuth_deg)
            sys_struct.transforms[1].heading   = deg2rad(10) - sys_struct.wings[1].heading
            
            # Apply the transformation without changing velocities
            SymbolicAWEModels.reposition!(sys_struct.transforms, sys_struct)
            
            # Reinitialize the solver to handle the state discontinuity
            SymbolicAWEModels.reinit!(sam, sam.prob, FBDF())

            if prn
                # Verify the new pose after one step
                next_step!(sam; dt=dt, set_values=set_values[step, :], vsm_interval)
                updated_elevation_deg = rad2deg(sys_struct.wings[1].elevation)
                println(">>> Pose updated. New Elevation is now: $(round(updated_elevation_deg, digits=2)) degrees.\n")
            end
        else
            # --- Normal simulation step ---
            next_step!(sam; dt=dt, set_values=set_values[step, :], vsm_interval)
        end
        
        # Log the state at the current time step
        update_sys_state!(sys_state, sam)
        sys_state.time = t
        log!(logger, sys_state)
    end

    if prn
        println("--- Simulation Finished ---")
        println("Runtime: $time")
        println("Times realtime: $(dt*steps/time)")
    end

    # Save and return the log
    mkpath(get_data_path())
    save_log(logger, "tmp_reposition_run")
    return load_log("tmp_reposition_run")
end


# Initialize the model
set = Settings("system.yaml")
set.sample_freq = 100
sam = SymbolicAWEModel(set)
init!(sam)
find_steady_state!(sam)

# Run the simulation with repositioning
lg = sim_reposition!(
    sam,
    total_time=10.0,
    reposition_interval_s=0.5,
    target_elevation_deg=50.0,
    target_azimuth_deg=10.0
)

# Plot the results
plot(sam.sys_struct, lg)


