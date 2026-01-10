# Copyright (c) 2025 Jelle Poland, Bart van de Lint
# SPDX-License-Identifier: MPL-2.0

"""
Real-Time 3D Visualization for V3 Kite with Keyboard Control

This example demonstrates real-time visualization of the v3 kite with REFINE wing
type and interactive keyboard control for steering and power/depower.

Key features:
- Real-time 3D updates using plot!(sys_struct)
- Interactive keyboard control:
  - Left Arrow:  Turn left (increase seg 89, decrease seg 87)
  - Right Arrow: Turn right (increase seg 87, decrease seg 89)
  - Down Arrow:  Power (decrease seg 88)
  - Up Arrow:    Depower (increase seg 88)
  - ESC:         Stop simulation
- REFINE wing with VSM panel forces lumped to structural points
- Configurable steering and power/depower magnitudes
"""

using GLMakie
using SymbolicAWEModels
using VortexStepMethod
using KiteUtils
using LinearAlgebra
using Statistics
using Printf

# ============================================================================
# SIMULATION PARAMETERS
# ============================================================================

const MODEL_NAME = "v3"
const SIM_TIME = 60.0
const FPS = 20
const N_STEPS = Int(round(FPS * SIM_TIME))
const REMAKE_CACHE = false

dt = 1.0 / FPS              # Time step [s]
realtime_factor = 1.0       # 1.0 = realtime
plot_interval = 1           # Update plot every N steps
vector_scale = 1.0          # Scale for wing orientation arrows

# Control parameters
max_steering = 0.3         # Max steering line length change [m]
max_power_depower = 0.3     # Max power/depower length change [m]
steering_rate = 0.001        # Steering change per step [m]
power_rate = 0.002           # Power/depower change per step [m]

# Initial damping
initial_damping = 10.0      # Initial world frame damping [N·s/m]
decay_time = 2.0            # Time for damping to decay to zero [s]

# Video recording
record_video = false
output_filename = "data/v3_realtime.mp4"
framerate = 20

# ============================================================================
# INITIALIZE MODEL
# ============================================================================

println("Initializing v3 kite model with REFINE wing...")
set_data_path("data/v3")
set = Settings("system.yaml")

# Load v3 system structure from YAML
model_name = hasproperty(set, :model_name) ? set.model_name : MODEL_NAME
struc_yaml = hasproperty(set, :struc_geometry_path) ?
    set.struc_geometry_path :
    joinpath("data", model_name, "struc_geometry_stable.yaml")
sys = load_sys_struct_from_yaml(struc_yaml; system_name=model_name, set=set)

# Initialize damping
SymbolicAWEModels.set_world_frame_damping(sys, initial_damping, 1:38)

# Verify REFINE wing setup
@assert length(sys.wings) > 0 "No wings in system"
@assert sys.wings[1].wing_type == SymbolicAWEModels.REFINE "Wing should be REFINE type"

wing_points = [p for p in sys.points if p.type == WING]
@info "REFINE wing setup:" n_wing_points=length(wing_points) n_panels=length(sys.wings[1].vsm_aero.panels) n_segments=length(sys.segments)

# Create symbolic model
sam = SymbolicAWEModel(set, sys)

# Initialize without VSM linearization (REFINE wings)
@info "Initializing model..."
SymbolicAWEModels.init!(sam; remake=REMAKE_CACHE, ignore_l0=false)

# Store initial segment lengths for control
seg_87_initial = sys.segments[87].l0
seg_88_initial = sys.segments[88].l0
seg_89_initial = sys.segments[89].l0

@info "Control segments initial lengths:" seg_87=seg_87_initial seg_88=seg_88_initial seg_89=seg_89_initial

# Settle initial state
@info "Settling initial state..."
[point.fix_static = true for point in sys.points if point.type == WING]
next_step!(sam; dt=10.0)
[point.fix_static = false for point in sys.points if point.type == WING]

# ============================================================================
# CREATE INITIAL PLOT
# ============================================================================

println("Creating 3D visualization window...")
scene = plot(sys; vector_scale, size=(1400, 900))
display(scene)

# Add progress text overlay
progress_text = Observable("Progress: 0%")
text!(scene, progress_text, position = Point2f(1380, 40), space = :pixel,
      fontsize = 20, color = :black, align = (:right, :top))

# Add keyboard control instructions
instructions = """
Keyboard Controls:
← Turn Left   → Turn Right
↓ Power       ↑ Depower
ESC to Stop
"""
text!(scene, instructions, position = Point2f(20, 130), space = :pixel,
      fontsize = 16, color = :darkblue, align = (:left, :top))

# Add control state display
control_text = Observable("Steering: 0.00m | Power: 0.00m")
text!(scene, control_text, position = Point2f(20, 60), space = :pixel,
      fontsize = 14, color = :darkgreen, align = (:left, :top))

# ============================================================================
# SETUP KEYBOARD CONTROL
# ============================================================================

println("Setting up keyboard controls...")

# Current control state
current_steering_delta = Ref(0.0)      # Steering offset from neutral
current_power_delta = Ref(0.0)         # Power/depower offset from neutral
stop_simulation = Ref(false)

# Keyboard event handler
on(events(scene).keyboardbutton) do event
    if event.action == Keyboard.press || event.action == Keyboard.repeat
        if event.key == Keyboard.left
            # Turn left: increase seg 89, decrease seg 87
            current_steering_delta[] = clamp(
                current_steering_delta[] - steering_rate,
                -max_steering, max_steering
            )
        elseif event.key == Keyboard.right
            # Turn right: increase seg 87, decrease seg 89
            current_steering_delta[] = clamp(
                current_steering_delta[] + steering_rate,
                -max_steering, max_steering
            )
        elseif event.key == Keyboard.down
            # Power: decrease seg 88
            current_power_delta[] = clamp(
                current_power_delta[] - power_rate,
                -max_power_depower, max_power_depower
            )
        elseif event.key == Keyboard.up
            # Depower: increase seg 88
            current_power_delta[] = clamp(
                current_power_delta[] + power_rate,
                -max_power_depower, max_power_depower
            )
        elseif event.key == Keyboard.escape
            stop_simulation[] = true
            println("\nESC pressed - stopping simulation...")
        end
    elseif event.action == Keyboard.release
        # On release, gradually return to neutral (could implement this)
        # For now, maintain current position
    end
end

println("Keyboard controls active:")
println("  ← Left:   Turn left")
println("  → Right:  Turn right")
println("  ↓ Down:   Power")
println("  ↑ Up:     Depower")
println("  ESC:      Stop simulation")

# ============================================================================
# RUN REAL-TIME SIMULATION
# ============================================================================

println("\nStarting real-time simulation...")
println("  Total time: $(SIM_TIME)s")
println("  Time step: $(dt)s")
println("  Realtime factor: $(realtime_factor)x")
if record_video
    println("  Recording video to: $(output_filename)")
end
println("\nSimulation running... Use arrow keys to control!\n")

# Initialize logger
logger = Logger(sam, N_STEPS + 1)
sys_state = SysState(sam)
sys_state.time = 0.0
log!(logger, sys_state)

# Simulation loop
start_time = time()
simulation_time = 0.0

# Initialize video recording if enabled
io = record_video ? VideoStream(scene; framerate) : nothing

for step in 1:N_STEPS
    if stop_simulation[]
        break
    end

    global simulation_time
    t = step * dt
    target_elapsed = t / realtime_factor

    # Update damping: linearly decay over decay_time
    if t <= decay_time
        current_damping = initial_damping * (1.0 - t / decay_time)
        SymbolicAWEModels.set_world_frame_damping(sys, current_damping, 1:38)
    else
        SymbolicAWEModels.set_world_frame_damping(sys, 0.0, 1:38)
    end

    # Apply control by updating segment lengths
    # Steering: seg 87 (right), seg 89 (left)
    # Positive steering_delta = turn right (increase 87, decrease 89)
    sys.segments[87].l0 = seg_87_initial + current_steering_delta[]
    sys.segments[89].l0 = seg_89_initial - current_steering_delta[]

    # Power/depower: seg 88
    # Positive power_delta = depower (increase 88)
    sys.segments[88].l0 = seg_88_initial + current_power_delta[]

    # Update control display
    control_text[] = @sprintf(
        "Steering: %.3fm | Power: %.3fm",
        current_steering_delta[], current_power_delta[]
    )

    # Simulation step
    step_start = time()
    try
        next_step!(sam; dt, vsm_interval=1)
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
    update_sys_state!(sys_state, sam)
    sys_state.time = t
    log!(logger, sys_state)

    # Update visualization
    if step % plot_interval == 0
        plot!(sys; vector_scale)
        progress_text[] = @sprintf("Progress: %d%%", round(Int, 100 * step / N_STEPS))

        if record_video
            recordframe!(io)
        end

        sleep(0.001)
    end

    actual_elapsed = time() - start_time
    sleep_time = max(0.0, target_elapsed - actual_elapsed)
    sleep(sleep_time)

    # Print progress every 10%
    if step % (N_STEPS ÷ 10) == 0
        avg_wing_pos = mean([p.pos_w for p in wing_points])
        @printf("  %.0f%% complete (t=%.1fs, wing_z=%.2fm)\n",
            100 * step / N_STEPS, t, avg_wing_pos[3])
    end
end

# Save video if recording
if record_video
    save(output_filename, io)
    println("\nVideo saved to: $(output_filename)")
end

total_elapsed = time() - start_time
println("\nSimulation complete!")
println("  Total runtime: $(round(total_elapsed, digits=2))s")
println("  Simulation time (next_step!): $(round(simulation_time, digits=2))s")
println("  Speedup (simulation only): $(round(SIM_TIME / simulation_time, digits=2))x realtime")
println("  Overall speedup: $(round(SIM_TIME / total_elapsed, digits=2))x realtime")

# ============================================================================
# SAVE AND REPLAY
# ============================================================================

println("\nSaving results...")
mkpath(get_data_path())
save_log(logger, "tmp_realtime_v3")
syslog = load_log("tmp_realtime_v3")

println("Creating interactive replay viewer...")
replay_scene = replay(syslog, sys; autoplay=false, loop=true)
display(replay_scene)

println("\nDone! Use the slider to replay the simulation.")

# Print final statistics
println("\n" * "="^60)
println("Final Results - V3 Kite REFINE Wing")
println("="^60)

if !isempty(wing_points)
    avg_x = mean([p.pos_w[1] for p in wing_points])
    avg_y = mean([p.pos_w[2] for p in wing_points])
    avg_z = mean([p.pos_w[3] for p in wing_points])
    println("  Average wing position: [$(round(avg_x, digits=2)), $(round(avg_y, digits=2)), $(round(avg_z, digits=2))] m")

    y_coords = [p.pos_w[2] for p in wing_points]
    span = maximum(y_coords) - minimum(y_coords)
    println("  Wing span: $(round(span, digits=2)) m")

    displacements = [norm(p.pos_w - p.pos_cad) for p in wing_points]
    avg_displacement = mean(displacements)
    max_displacement = maximum(displacements)
    println("  Average structural displacement: $(round(avg_displacement, digits=3)) m")
    println("  Maximum structural displacement: $(round(max_displacement, digits=3)) m")
end

println("\n  Final control state:")
println("    Steering offset: $(round(current_steering_delta[], digits=3)) m")
println("    Power offset: $(round(current_power_delta[], digits=3)) m")
println("="^60)

nothing
