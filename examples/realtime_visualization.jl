# Copyright (c) 2025 Bart van de Lint
# SPDX-License-Identifier: MPL-2.0

"""
Real-Time 3D Visualization Example

This example demonstrates how to create a custom simulation loop with real-time
3D visualization using the time-based plot() API.

Key features:
- Real-time 3D updates using plot(sys_struct, time)
- Automatic time display and background pane updates
- Proper sleep timing to maintain real-time speed
- Interactive camera control during simulation (hover, click-to-zoom)
- Configurable visualization frame rate
- Clean API - just call plot(sys_struct, t) to update!
"""

using GLMakie
using SymbolicAWEModels
using KiteUtils
using LinearAlgebra
using Printf

# ============================================================================
# SIMULATION PARAMETERS
# ============================================================================

dt = 0.05                    # Time step [s]
total_time = 20.0           # Total simulation time [s]
vsm_interval = 3            # VSM update interval
realtime_factor = 1.0       # 1.0 = realtime, 2.0 = 2x speed, 0.5 = half speed
plot_interval = 1           # Update plot every N steps (1 = every step)
vector_scale = 1.0         # Scale for wing orientation arrows

# Steering parameters
steering_freq = 0.5         # Hz - full left-right cycle frequency
steering_magnitude = 10.0   # Magnitude of steering input [Nm]
bias = 0.2                  # Steering bias

# ============================================================================
# INITIALIZE MODEL
# ============================================================================

println("Initializing model...")
set_data_path("data")
set = Settings("system.yaml")
set.profile_law = 3
sam = SymbolicAWEModel(set)
init!(sam)
find_steady_state!(sam)

sys_struct = sam.sys_struct
steps = Int(round(total_time / dt))
num_winches = length(sys_struct.winches)

# ============================================================================
# CREATE INITIAL PLOT
# ============================================================================

println("Creating 3D visualization window...")

# Create initial 3D plot using time-based API (time=0.0)
# This automatically creates observables and sets up the scene
scene = plot(sys_struct, 0.0; vector_scale, size=(1400, 900))

# Add progress text overlay (time display is already added by plot function)
progress_text = Observable("Progress: 0%")
text!(scene, progress_text, position = Point2f(1380, 40), space = :pixel,
      fontsize = 20, color = :black, align = (:right, :top))

# ============================================================================
# PREPARE CONTROL INPUTS
# ============================================================================

println("Preparing control inputs...")
set_values = zeros(Float64, steps, num_winches)

for step in 1:steps
    t = (step-1) * dt
    steering = steering_magnitude * cos(2π * steering_freq * t + bias)
    set_values[step, :] = [0.0, steering, -steering]
end

# ============================================================================
# RUN REAL-TIME SIMULATION
# ============================================================================

println("\nStarting real-time simulation...")
println("  Total time: $(total_time)s")
println("  Time step: $(dt)s")
println("  Realtime factor: $(realtime_factor)x")
println("  Plot update interval: every $(plot_interval) step(s)")
println("\nSimulation running... (you can interact with the 3D view)\n")

# Initialize state
logger = SymbolicAWEModels.Logger(length(sys_struct.points), steps+1)
sys_state = SysState(sam)
sys_state.time = 0.0
SymbolicAWEModels.log!(logger, sys_state)

steady_torque = SymbolicAWEModels.calc_steady_torque(sam)
torque_damp = 0.9
set_torques = similar(set_values)

# Simulation loop with real-time visualization
start_time = time()
for step in 1:steps
    global steady_torque  # Declare that we're modifying the global variable
    t = step * dt
    target_elapsed = t / realtime_factor

    # Calculate control torques
    steady_torque = torque_damp * steady_torque + (1-torque_damp) * SymbolicAWEModels.calc_steady_torque(sam)
    set_torques[step, :] = steady_torque .+ set_values[step, :]

    # Simulation step
    try
        next_step!(sam; set_values=set_torques[step, :], dt, vsm_interval)
    catch e
        if e isa AssertionError
            @warn "Simulation crashed at t=$t"
            break
        else
            rethrow(e)
        end
    end

    # Update system state and log
    SymbolicAWEModels.update_sys_state!(sys_state, sam)
    sys_state.time = t
    SymbolicAWEModels.log!(logger, sys_state)

    # Update visualization
    if step % plot_interval == 0
        # Update plot using time-based API
        # This automatically updates observables, time display, and background panes
        plot(sys_struct, t; vector_scale)

        # Update progress text overlay
        progress_text[] = @sprintf("Progress: %d%%", round(Int, 100 * step / steps))

        # Force Makie to process events and update display
        sleep(0.001)
    end

    # Sleep to maintain real-time pacing
    actual_elapsed = time() - start_time
    sleep_time = max(0.0, target_elapsed - actual_elapsed)
    sleep(sleep_time)

    # Print progress every 10%
    if step % (steps ÷ 10) == 0
        @printf("  %.0f%% complete (t=%.1fs)\n", 100 * step / steps, t)
    end
end

total_elapsed = time() - start_time
println("\nSimulation complete!")
println("  Total runtime: $(round(total_elapsed, digits=2))s")
println("  Speedup: $(round(total_time / total_elapsed, digits=2))x realtime")

# ============================================================================
# SAVE AND PLOT RESULTS
# ============================================================================

println("\nSaving results...")
mkpath(get_data_path())
SymbolicAWEModels.save_log(logger, "tmp_realtime_run")
lg = load_log("tmp_realtime_run")

println("Creating post-simulation plots...")
fig_results = plot(sam.sys_struct, lg; plot_default=false, plot_heading=true, plot_aoa=true)
display(fig_results)

println("\nDone! Close the windows to exit.")
