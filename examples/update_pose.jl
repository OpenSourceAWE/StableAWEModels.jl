# Copyright (c) 2025 Bart van de Lint
# SPDX-License-Identifier: MPL-2.0

using Pkg
if ! ("ControlPlots" ∈ keys(Pkg.project().dependencies))
    using TestEnv; TestEnv.activate()
end
using SymbolicAWEModels, LinearAlgebra, ControlPlots, OrdinaryDiffEqBDF

set = Settings("system.yaml")
sam = SymbolicAWEModel(set)

"""
    run_pose_update_example()

Demonstrates how to use `update_pose!` to adjust the elevation of a kite system
during a simulation.
"""
function run_pose_update_example()
    # 1. --- Initialize the Model ---
    # Load default settings and create a symbolic AWE model.
    init!(sam)
    find_steady_state!(sam)
    sys_struct = sam.sys_struct
    
    # Simulation parameters
    dt = 1/set.sample_freq
    total_time = 20.0
    steps = Int(round(total_time / dt))
    set_values = zeros(Float64, steps, length(sys_struct.winches))

    println("--- Starting Simulation ---")

    # 2. --- Run the first part of the simulation ---
    # Simulate for 5 seconds to let the system settle.
    for step in 1:Int(round(1.0 / dt))
        set_values[step, :] = -sam.set.drum_radius .* [norm(winch.force) for winch in sys_struct.winches]
        next_step!(sam; dt=dt, set_values=set_values[step, :])
        plot(sam, dt*(step-1))
    end

    # 3. --- Adjust the Pose ---
    # Get the elevation before the update. Note that the wing's elevation is a derived
    # property and is updated internally during the simulation steps.
    # Here we read it directly from the transform object.
    initial_elevation_deg = rad2deg(sys_struct.transforms[1].elevation)
    println("\n>>> Time: 5.0s. Initial Elevation: $(round(initial_elevation_deg, digits=2)) degrees.")

    # Define the new target elevation (increase by 10 degrees)
    new_elevation_rad = sys_struct.transforms[1].elevation + deg2rad(10.0)
    sys_struct.transforms[1].elevation = new_elevation_rad
    
    println(">>> Applying `update_pose!` to set elevation to $(round(rad2deg(new_elevation_rad), digits=2)) degrees...")
    
    # Apply the transformation. This will rotate the entire system to the new elevation
    # while keeping the velocities of all components intact.
    SymbolicAWEModels.update_pose!(sys_struct.transforms, sys_struct)
    
    # To verify, we can check the wing's derived elevation property after one more step.

    SymbolicAWEModels.reinit!(sam, sam.prob, FBDF())
    next_step!(sam; dt=dt, set_values=set_values[Int(round(5.0 / dt)) + 1, :])
    updated_elevation_deg = rad2deg(sys_struct.wings[1].elevation)
    println(">>> Pose updated. New Elevation at t=5.0s is now: $(round(updated_elevation_deg, digits=2)) degrees.\n")

    # 4. --- Continue the simulation ---
    # Run the simulation for the remaining time from the new pose.
    println("--- Resuming Simulation ---")
    start_step = Int(round(1.0 / dt)) + 2
    for step in start_step:steps
         set_values[step, :] = -sam.set.drum_radius .* [norm(winch.force) for winch in sys_struct.winches]
         next_step!(sam; dt=dt, set_values=set_values[step, :])
         plot(sam, dt*(step-1))
    end
    
    final_elevation_deg = rad2deg(sys_struct.wings[1].elevation)
    println("\n--- Simulation Finished at t=$(total_time)s ---")
    println("Final Elevation: $(round(final_elevation_deg, digits=2)) degrees.")
end

# To run the example:
run_pose_update_example()

