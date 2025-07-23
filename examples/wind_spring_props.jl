using SymbolicAWEModels, VortexStepMethod, KiteUtils, WinchModels
using ControlPlots, Statistics, LinearAlgebra
using OrdinaryDiffEqCore

# Assuming 'sam' setup code from your snippet has been run
set = Settings("system.yaml")
set.sample_freq = 800
set.abs_tol = 1e-5
set.rel_tol = 1e-5
dt = 1/set.sample_freq
sam = SymbolicAWEModel(set)
SymbolicAWEModels.init!(sam)
sys = sam.sys

tsys = SymbolicAWEModels.create_tether_sys_struct(set)
tsam = SymbolicAWEModel(set, tsys)
init!(tsam)

find_steady_state!(sam; t=10.0, dt=3.0)
SymbolicAWEModels.copy!(sam.sys_struct, tsam.sys_struct)
OrdinaryDiffEqCore.reinit!(tsam.integrator; reinit_dae=true)
SymbolicAWEModels.update_sys_struct!(tsam, tsam.sys_struct)
F_0 = [-tsys.points[i].force for i in 1:4]

steps = 200
F_step = -0.1

function response(sam, steps, F_step, F_0;
                  abs_tol=1e-6,
                  consecutive_steps_needed=10)

    points = sam.sys_struct.points

    initial_tether_lens = [norm(points[i].pos_w) for i in 1:4]
    [points[i].disturb .= F_0[i] .+ F_step * normalize(points[i].pos_w) for i in 1:4]
    @show initial_tether_lens

    tether_lens = zeros(4, steps+1)
    tether_lens[:, 1] .= initial_tether_lens # Store the initial lengths
    settled_steps = 0
    @time for step in 1:steps
        next_step!(sam; vsm_interval=0)
        for j in 1:4
            tether_lens[j, step+1] = norm(points[j].pos_w)
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

tether_lens = response(tsam, steps, F_step, F_0)

display(plotx(
    dt .* collect(1:steps+1), 
    tether_lens[1,:], tether_lens[2,:], tether_lens[3,:], tether_lens[4,:];
    title="Force step response",
    ylabels=["Left power [m]", "Right power [m]", "Left steering [m]", "Right steering [m]"],
))

k_values = zeros(4)
c_values = zeros(4)

for j in 1:4
    tether_len_series = tether_lens[j, :]
    initial_len = tether_len_series[1]
    final_len = tether_len_series[end]

    delta_x_ss = final_len - initial_len

    if abs(delta_x_ss) < 1e-6 # Avoid division by zero or very small numbers
        println("Warning: Steady-state change in length is too small for Tether $(j). Cannot reliably calculate k.")
        k_values[j] = NaN
        c_values[j] = NaN
        continue
    end

    # Calculate Spring Stiffness (k)
    k = F_step / delta_x_ss * initial_len
    k_values[j] = k

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
    end

    # Calculate Damping Coefficient (c)
    if !isnan(tau)
        c = k * tau
        c_values[j] = c
    else
        c_values[j] = NaN
    end
end

println("Summary of Results:")
for j in 1:4
    println("Tether $(j): k = $(k_values[j]) N, c = $(c_values[j]) Ns")
end

set.segments = 1
ssys = SymbolicAWEModels.create_tether_sys_struct(set; 
                                                  axial_stiffness=k_values, 
                                                  axial_damping=c_values)
ssam = SymbolicAWEModel(set, ssys)
init!(ssam)

forces = [F ⋅ normalize(point.pos_w) for (F, point) in zip(F_0, ssys.points[1:4])]
SymbolicAWEModels.copy!(sam.sys_struct, ssam.sys_struct)
OrdinaryDiffEqCore.reinit!(ssam.integrator; reinit_dae=true)
SymbolicAWEModels.update_sys_struct!(ssam, ssam.sys_struct)

tether_lens = response(ssam, steps, 0.0,    F_0)
tether_lens = response(ssam, steps, F_step, F_0)

display(plotx(
    dt .* collect(1:steps+1), 
    tether_lens[1,:], tether_lens[2,:], tether_lens[3,:], tether_lens[4,:];
    title="Force step response",
    ylabels=["Left power [m]", "Right power [m]", "Left steering [m]", "Right steering [m]"],
))

