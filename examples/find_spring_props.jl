using SymbolicAWEModels, VortexStepMethod, KiteUtils, WinchModels
using ControlPlots, Statistics, LinearAlgebra

# Assuming 'sam' setup code from your snippet has been run
set = Settings("system.yaml")
dt = 1/set.sample_freq
sam = SymbolicAWEModel(set)
SymbolicAWEModels.init!(sam)
sys = sam.sys

function response(sam, steps, τ_step, τ_0;
                  abs_tol=1e-8,
                  consecutive_steps_needed=5)

    winches = sam.sys_struct.winches
    initial_tether_lens = [winches[j].tether_len for j in 1:3]
    @show initial_tether_lens

    if length(sam.sys_struct.wings) > 0
        sam.integrator.ps[sam.sys.fix_wing] = true
    end

    tether_lens = zeros(3, steps+1)
    tether_lens[:, 1] .= initial_tether_lens # Store the initial lengths
    set_values = τ_0 .+ τ_step               # Apply τ_step
    settled_steps = 0
    @time for step in 1:steps
        next_step!(sam; set_values, vsm_interval=0)
        for j in 1:3
            tether_lens[j, step+1] = winches[j].tether_len
        end

        # Check absolute delta for all tethers
        step_deltas = abs.(tether_lens[:, step+1] .- tether_lens[:, step])
        max_delta = maximum(step_deltas)
        if max_delta < abs_tol
            settled_steps += 1
        else
            settled_steps = 0
        end
        if settled_steps >= consecutive_steps_needed
            println("Stopped at step $step: all tethers within $abs_tol for $consecutive_steps_needed steps.")
            tether_lens[:, step+2:end] .= tether_lens[:, step+1]
            break
        end
    end

    return tether_lens
end

find_steady_state!(sam; t=10.0, dt=3.0)
steps = 1000
τ_step = -1.0 # Newtons
τ_0 = -sam.set.drum_radius .* sam.integrator[sys.winch_force] 
tether_lens = response(sam, steps, τ_step, τ_0)

display(plotx(
    dt .* collect(1:steps+1), 
    tether_lens[1,:], tether_lens[2,:], tether_lens[3,:];
    title="Force step response",
    ylabels=["Power tether [m]", "Left tether [m]", "Right tether [m]"],
))

println("Analysis of Tether Stiffness and Damping:")

k_values = zeros(3)
c_values = zeros(3)

for j in 1:3
    println("\n--- Tether $(j) ---")

    tether_len_series = tether_lens[j, :]
    initial_len = tether_len_series[1]
    final_len = mean(tether_len_series[end-5:end]) # Average last 5 points

    delta_x_ss = final_len - initial_len
    println("Initial Tether Length: $(initial_len) m")
    println("Final (Steady-State) Tether Length: $(final_len) m")
    println("Steady-State Change in Length (Delta_x_ss): $(delta_x_ss) m")

    if abs(delta_x_ss) < 1e-6 # Avoid division by zero or very small numbers
        println("Warning: Steady-state change in length is too small for Tether $(j). Cannot reliably calculate k.")
        k_values[j] = NaN
        c_values[j] = NaN
        continue
    end

    # Calculate Spring Stiffness (k)
    F_step = τ_step / -sam.set.drum_radius
    k = F_step / delta_x_ss * initial_len
    k_values[j] = -k
    println("Spring stiffness constant (k): $(k) N")

    # Calculate Time Constant (tau)
    # Target value for (1 - 1/ℯ) of the change
    target_len = initial_len + (1 - 1/ℯ) * delta_x_ss

    # Find the index where the series crosses this target value
    tau_index = -1
    for i in eachindex(tether_len_series)
        if (delta_x_ss > 0 && tether_len_series[i] >= target_len) ||
           (delta_x_ss < 0 && tether_len_series[i] <= target_len)
            tau_index = i
            break
        end
    end

    if tau_index == -1
        println("Warning: Could not find time constant (tau) for Tether $(j).
            Response might not have settled enough.")
        tau = NaN
    else
        # Interpolate for better accuracy, or just use the found index
        # For simplicity, using the index
        tau = (tau_index - 1) * dt # (index - 1) because time starts from 0 for the change
        println("Time Constant (tau): $(tau) seconds")
    end

    # Calculate Damping Coefficient (c)
    if !isnan(tau)
        c = k * tau
        c_values[j] = c
        println("Damping coefficient constant (c): $(c) Ns")
    else
        c_values[j] = NaN
    end
end

println("\nSummary of Results:")
for j in 1:3
    println("Tether $(j): k = $(k_values[j]) N, c = $(c_values[j]) Ns")
end

set = Settings("system.yaml")
stiffness = k_values ./ set.l_tether
l0=set.l_tether .+ τ_0./-sam.set.drum_radius ./ stiffness
points = [
    Point(1, [0, 0, l0[1]], STATIC)
    Point(2, [0, 0, l0[2]], STATIC)
    Point(3, [0, 0, l0[3]], STATIC)
    Point(4, [0, 0, 0], STATIC)
    Point(5, [0, 0, 0], STATIC)
    Point(6, [0, 0, 0], STATIC)
]
segments = [
    Segment(1, (1,4), k_values[1], c_values[1], √2 * 1e-3 * set.power_tether_diameter)
    Segment(2, (2,5), k_values[2], c_values[2], 1e-3 * set.steering_tether_diameter)
    Segment(3, (3,6), k_values[3], c_values[3], 1e-3 * set.steering_tether_diameter)
]
tethers = [
    Tether(1, [1], 4)
    Tether(2, [2], 5)
    Tether(3, [3], 6)
]
winches = [
    Winch(1, TorqueControlledMachine(set), [1]; tether_len=set.l_tether)
    Winch(2, TorqueControlledMachine(set), [2]; tether_len=set.l_tether)
    Winch(3, TorqueControlledMachine(set), [3]; tether_len=set.l_tether)
]
transforms = [
    Transform(1, deg2rad(90), 0.0, 0.0; base_point_idx=4, base_pos=zeros(3), rot_point_idx=1)
]
sys_struct = SystemStructure("one_seg_tether", set;
    points, segments, tethers, winches, transforms)
ssam = SymbolicAWEModel(set, sys_struct)
init!(ssam)
@show ssam.integrator[ssam.sys.winch_force]

steps = 1000
tether_lens = response(ssam, steps, τ_step, τ_0)

display(plotx(
    dt .* collect(1:steps+1), 
    tether_lens[1,:], tether_lens[2,:], tether_lens[3,:];
    title="Force step response",
    ylabels=["Power tether [m]", "Left tether [m]", "Right tether [m]"],
))
