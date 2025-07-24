
# Helper: Check if all values after index i are within p% band of steady-state
function in_percent_band(x, steady, delta_x, i, p)
    tol = p/100 * abs(delta_x)
    # All subsequent points must be within steady ± p%
    all(abs.(x[i:end] .- steady) .<= tol)
end

function calc_spring_props(sam::SymbolicAWEModel, tether_lens, F_step; p=5, prn=false)
    @unpack tethers, segments = sam.sys_struct
    set = sam.set
    dt = 1/set.sample_freq

    k_values = zeros(4)
    c_values = zeros(4)

    diameters = [segments[tether.segment_idxs[1]].diameter for tether in tethers]
    mass_per_meter = set.rho_tether * π * ((diameters/2).^2)

    for j in eachindex(tethers)
        tether_len_series = tether_lens[j, :]
        initial_len = tether_len_series[1]
        final_len = tether_len_series[end]
        delta_x_ss = final_len - initial_len

        m = mass_per_meter[j] * 0.5 * tethers[j].stretched_len
        @assert m > 0

        if abs(delta_x_ss) < 1e-6
            @warn "Steady-state change too small for Tether $j; skipping."
            k_values[j] = NaN; c_values[j] = NaN
            continue
        end

        # Spring stiffness
        k = F_step / delta_x_ss
        k_values[j] = k
        ω_n = sqrt(k/m)

        # Find settling time index according to your in_Xpercent_band function:
        T_s_index = -1
        for i in 1:length(tether_len_series)
            if in_percent_band(tether_len_series, final_len, delta_x_ss, i, p)
                T_s_index = i
                break
            end
        end

        if T_s_index == -1
            @warn "Could not find settling time ($p% criterion) for Tether $j; using fallback."
            # Fallback: use time constant method to find tau
            target_len = initial_len + (1 - 1/ℯ) * delta_x_ss
            tau_idx = findfirst(i ->
                (delta_x_ss > 0 && tether_len_series[i] >= target_len) ||
                (delta_x_ss < 0 && tether_len_series[i] <= target_len),
                1:length(tether_len_series))
            if tau_idx === nothing
                @warn "Cannot determine tau for Tether $j"
                c_values[j] = NaN
            else
                tau = (tau_idx-1)*dt
                c = k * tau
                c_values[j] = c
                println("Tether $j fallback c=", c)
            end
            continue
        end

        T_s = (T_s_index - 1) * dt
        # Calculate damping ratio based on variable percentage settling criterion:
        X = -log(p / 100)
        zeta = X / (ω_n * T_s)
        c = 2 * zeta * sqrt(k * m)
        c_values[j] = c
        prn && println("Tether $j: ω_n=$(round(ω_n,digits=3)) rad/s,
                T_s=$(round(T_s,digits=3)) s, 
                ζ=$(round(zeta,digits=4)), c=$(round(c,digits=4)) Ns/m")
    end

    prn && println("Summary of Results:")
    prn && for j in eachindex(tethers)
        println("Tether $(j): k = $(k_values[j]) N/m, c = $(c_values[j]) Ns/m")
    end
    return k_values, c_values
end

function step(sam::SymbolicAWEModel, steps, F_step, F_0;
                  abs_tol=1e-6,
                  consecutive_steps_needed=10,
                  prn=false)

    @unpack points, tethers = sam.sys_struct

    initial_tether_lens = [norm(points[i].pos_w) for i in eachindex(tethers)]
    [points[i].disturb .= F_0[i] .+ F_step * normalize(points[i].pos_w) for i in eachindex(tethers)]

    tether_lens = zeros(length(tethers), steps+1)
    tether_lens[:, 1] .= initial_tether_lens # Store the initial lengths
    settled_steps = 0
    for step in 1:steps
        next_step!(sam; vsm_interval=0)
        for j in eachindex(tethers)
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
            prn && println("Stopped at step $step: all tethers within $abs_tol for $consecutive_steps_needed steps.")
            tether_lens[:, step+2:end] .= tether_lens[:, step+1]
            break
        end
    end
    return tether_lens
end

