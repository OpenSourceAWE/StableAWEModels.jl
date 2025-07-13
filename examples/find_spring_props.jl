using SymbolicAWEModels, VortexStepMethod, KiteUtils, ControlPlots, Statistics

# Assuming 'sam' setup code from your snippet has been run
set = Settings("system.yaml")
sam = SymbolicAWEModel(set)
SymbolicAWEModels.init!(sam)
sys = sam.sys
find_steady_state!(sam; t=10.0, dt=3.0)

winches = sam.sys_struct.winches
initial_tether_lens = [winches[j].tether_len for j in 1:3]
@show initial_tether_lens

steps = 100
winches = sam.sys_struct.winches
sam.integrator.ps[sys.fix_wing] = true
next_step!(sam; dt=1e-6, vsm_interval=0)

delta_F = 1.0 # Newtons
tether_lens = zeros(3, steps+1)
tether_lens[:, 1] .= initial_tether_lens # Store the initial lengths
set_values = -sam.set.drum_radius .* sam.integrator[sys.winch_force] .+ delta_F # Apply delta_F
@time for i in 1:steps
    next_step!(sam; set_values, vsm_interval=0)
    @show winches[1].tether_len
    @show sam.integrator[sys.tether_len[1]]
    [tether_lens[j, i + 1] = winches[j].tether_len for j in 1:3] # Store after step
end

display(plotx(
    dt .* collect(1:steps+1), 
    tether_lens[1,:], tether_lens[2,:], tether_lens[3,:];
    title="Force step response",
    ylabels=["Power tether [m]", "Left tether [m]", "Right tether [m]"],
))

dt = 1/set.sample_freq

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

    if abs(delta_x_ss) < 1e-9 # Avoid division by zero or very small numbers
        println("Warning: Steady-state change in length is too small for Tether $(j). Cannot reliably calculate k.")
        k_values[j] = NaN
        c_values[j] = NaN
        continue
    end

    # Calculate Spring Stiffness (k)
    k = delta_F / delta_x_ss * initial_len
    k_values[j] = k
    println("Spring stiffness constant (k): $(k) N/m")

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
    println("Tether $(j): k = $(k_values[j]) N/m, c = $(c_values[j]) Ns/m")
end

# You can also plot to visually inspect the fits:
# Using ControlPlots, you could try to overlay theoretical exponential response.
# For example, for tether 1:
# t_sim = (0:steps) .* dt_sim
# theoretical_x = initial_tether_lens[1] .+ delta_x_ss_tether1 .* (1 .- exp.(-t_sim ./ tau_tether1))
# plotx(t_sim, all_tether_lens[1,:], theoretical_x) # You'd need to replace with actual values
