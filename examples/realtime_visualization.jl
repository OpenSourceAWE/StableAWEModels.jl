# Copyright (c) 2025 Bart van de Lint
# SPDX-License-Identifier: MPL-2.0

"""
Real-Time 3D Visualization Example with Keyboard Control

This example demonstrates how to create a custom simulation loop with real-time
3D visualization and interactive keyboard control using the time-based plot() API.

Key features:
- Real-time 3D updates using plot(sys_struct, time)
- Interactive keyboard control: Arrow keys for steering
  - Left Arrow:  Turn left   [0.0, -mag, mag]
  - Right Arrow: Turn right  [0.0, mag, -mag]
  - Down Arrow:  Power       [0.0, -mag, -mag]
  - Up Arrow:    Depower     [0.0, mag, mag]
  - ESC:         Stop simulation
- Automatic time display and background pane updates
- Proper sleep timing to maintain real-time speed
- Interactive camera control during simulation (hover, click-to-zoom)
- Configurable visualization frame rate
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
total_time = 60.0           # Total simulation time [s]
vsm_interval = 3            # VSM update interval
realtime_factor = 1.0       # 1.0 = realtime, 2.0 = 2x speed, 0.5 = half speed
plot_interval = 1           # Update plot every N steps (1 = every step)
vector_scale = 1.0         # Scale for wing orientation arrows

# Steering parameters
steering_magnitude = 5.0   # Magnitude of steering input [Nm]

# Video recording parameters
record_video = false        # Set to true to record video to MP4
output_filename = "data/kite_simulation.mp4"  # Output video filename
framerate = 20              # Video framerate (frames per second)

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

# Add progress text overlay and keyboard instructions
progress_text = Observable("Progress: 0%")
text!(scene, progress_text, position = Point2f(1380, 40), space = :pixel,
      fontsize = 20, color = :black, align = (:right, :top))

# Add keyboard control instructions
instructions = """
Keyboard Controls:
← Left   → Right
↓ Down   ↑ Up
ESC to Stop
"""
text!(scene, instructions, position = Point2f(20, 130), space = :pixel,
      fontsize = 16, color = :darkblue, align = (:left, :top))

# ============================================================================
# SETUP KEYBOARD CONTROL
# ============================================================================

println("Setting up keyboard controls...")

# Current steering state (will be updated by keyboard events)
current_steering = Observable([0.0, 0.0, 0.0])  # [power, steer_left, steer_right]
stop_simulation = Ref(false)  # Flag to stop simulation on ESC

# Keyboard event handler
on(events(scene).keyboardbutton) do event
    if event.action == Keyboard.press || event.action == Keyboard.repeat
        mag = steering_magnitude
        if event.key == Keyboard.left
            # Left: turn left
            current_steering[] = [0.0, -mag, mag]
        elseif event.key == Keyboard.right
            # Right: turn right
            current_steering[] = [0.0, mag, -mag]
        elseif event.key == Keyboard.down
            # Down: pitch down
            current_steering[] = [0.0, -mag, -mag]
        elseif event.key == Keyboard.up
            # Up: pitch up
            current_steering[] = [0.0, mag, mag]
        elseif event.key == Keyboard.escape
            # ESC: stop simulation
            stop_simulation[] = true
            println("\nESC pressed - stopping simulation...")
        end
    elseif event.action == Keyboard.release
        # Release: return to neutral (but not for ESC)
        if event.key != Keyboard.escape
            current_steering[] = [0.0, 0.0, 0.0]
        end
    end
end

println("Keyboard controls active:")
println("  ← Left:  Turn left")
println("  → Right: Turn right")
println("  ↓ Down:  Pitch down")
println("  ↑ Up:    Pitch up")
println("  ESC:     Stop simulation")

# ============================================================================
# RUN REAL-TIME SIMULATION
# ============================================================================

println("\nStarting real-time simulation...")
println("  Total time: $(total_time)s")
println("  Time step: $(dt)s")
println("  Realtime factor: $(realtime_factor)x")
println("  Plot update interval: every $(plot_interval) step(s)")
if record_video
    println("  Recording video to: $(output_filename)")
    println("  Video framerate: $(framerate) fps")
end
println("\nSimulation running... Use arrow keys to control the kite!\n")

# Initialize state
logger = SymbolicAWEModels.Logger(length(sys_struct.points), steps+1)
sys_state = SysState(sam)
sys_state.time = 0.0
SymbolicAWEModels.log!(logger, sys_state)

steady_torque = SymbolicAWEModels.calc_steady_torque(sam)
torque_damp = 0.9

# Simulation loop with real-time visualization
start_time = time()
simulation_time = 0.0  # Track actual time spent in simulation (next_step!)

# Initialize video recording if enabled
io = record_video ? VideoStream(scene; framerate) : nothing

for step in 1:steps
    # Check if user pressed ESC to stop
    if stop_simulation[]
        break
    end

    global steady_torque, simulation_time  # Declare that we're modifying global variables
    t = step * dt
    target_elapsed = t / realtime_factor

    # Calculate control torques using current keyboard input
    steady_torque = torque_damp * steady_torque + (1-torque_damp) * SymbolicAWEModels.calc_steady_torque(sam)
    control_input = steady_torque .+ current_steering[]

    # Simulation step - measure only this time
    step_start = time()
    try
        next_step!(sam; set_values=control_input, dt, vsm_interval)
    catch e
        if e isa AssertionError
            @warn "Simulation crashed at t=$t"
            break
        else
            rethrow(e)
        end
    end
    simulation_time += time() - step_start

    # Update system state and log
    SymbolicAWEModels.update_sys_state!(sys_state, sam)
    sys_state.time = t
    SymbolicAWEModels.log!(logger, sys_state)

    # Update visualization
    if step % plot_interval == 0
        # Update plot using time-based API
        # This automatically updates observables, time display, and background panes
        plot(sys_struct, t; vector_scale)
        display(scene)

        # Update progress text overlay
        progress_text[] = @sprintf("Progress: %d%%", round(Int, 100 * step / steps))

        # Record frame if video recording is enabled
        if record_video
            recordframe!(io)
        end

        # Force Makie to process events and update display
        sleep(0.001)
    end

    actual_elapsed = time() - start_time
    sleep_time = max(0.0, target_elapsed - actual_elapsed)
    sleep(sleep_time)

    # Print progress every 10%
    if step % (steps ÷ 10) == 0
        @printf("  %.0f%% complete (t=%.1fs)\n", 100 * step / steps, t)
    end
end

# Save video if recording was enabled
if record_video
    save(output_filename, io)
    println("\nVideo saved to: $(output_filename)")
end

total_elapsed = time() - start_time
println("\nSimulation complete!")
println("  Total runtime: $(round(total_elapsed, digits=2))s")
println("  Simulation time (next_step!): $(round(simulation_time, digits=2))s")
println("  Speedup (simulation only): $(round(total_time / simulation_time, digits=2))x realtime")
println("  Overall speedup: $(round(total_time / total_elapsed, digits=2))x realtime")

# ============================================================================
# SAVE AND PLOT RESULTS
# ============================================================================

println("\nSaving results...")
mkpath(get_data_path())
SymbolicAWEModels.save_log(logger, "tmp_realtime_run")
lg = load_log("tmp_realtime_run")

println("Creating post-simulation plots...")
fig_results = plot(sam.sys_struct, lg)
display(fig_results)

println("\nDone! Close the windows to exit.")
