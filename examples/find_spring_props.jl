using SymbolicAWEModels, VortexStepMethod, KiteUtils, WinchModels
using ControlPlots, Statistics, LinearAlgebra

# Assuming 'sam' setup code from your snippet has been run
set = Settings("system.yaml")
dt = 1/set.sample_freq
sam = SymbolicAWEModel(set)
SymbolicAWEModels.init!(sam)
sys = sam.sys

function response(sam, steps)
    find_steady_state!(sam; t=10.0, dt=3.0)

    winches = sam.sys_struct.winches
    initial_tether_lens = [winches[j].tether_len for j in 1:3]
    @show initial_tether_lens

    winches = sam.sys_struct.winches
    sam.integrator.ps[sys.fix_wing] = true

    delta_F = -1.0 # Newtons
    tether_lens = zeros(3, steps+1)
    tether_lens[:, 1] .= initial_tether_lens # Store the initial lengths
    set_values = -sam.set.drum_radius .* sam.integrator[sys.winch_force] .+ delta_F # Apply delta_F
    last_abs = Inf
    rel_error = Inf
    @time for step in 1:steps
        global last_abs, rel_error
        next_step!(sam; set_values, vsm_interval=0)
        [tether_lens[j, step+1] = winches[j].tether_len for j in 1:3] # Store after step
        abs_error = abs(tether_lens[1, step+1] - tether_lens[1, step])
        rel_error = abs(abs_error - last_abs)
        last_abs = abs_error
        if rel_error < 1e-8
            println("Relative error: $rel_error \t Absolute error: $abs_error")
            println("Stopped at step $step")
            tether_lens[:, step+2:end] .= tether_lens[:, step+1]
            break
        end
        step += 1
    end
    return tether_lens
end

steps = 1000
tether_lens = response(sam, steps)

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
    k = delta_F / delta_x_ss * initial_len
    k_values[j] = k
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

points = [
    Point(1, [0, 0, set.l_tether], STATIC)
    Point(2, [0, 0, 0], STATIC)
]
segments = [
    Segment(1, (1,2), k_values[1], c_values[1], 1e-3*set.power_tether_diameter)
]
tethers = [
    Tether(1, [1])
]
winches = [
    Winch(1, TorqueControlledMachine(set), [1])
]
transforms = [
    Transform(1, deg2rad(90), 0.0, 0.0; base_point_idx=1, base_pos=zeros(3), rot_point_idx=2)
]
sys_struct = SystemStructure("one_seg_tether", set;
    points, segments, tethers, winches, transforms)
ssam = SymbolicAWEModel(set, sys_struct)
init!(ssam)


