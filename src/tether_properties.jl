
# Copies the state from one sam to another sam
function copy_to_simple!(sys::SystemStructure, ssys::SystemStructure)
    (sys.name != "ram") && error("provide a ram sys as the first argument")
    (ssys.name != "simple_ram") && error("provide a simple ram sys as the second argument")

    # copy point pos and vel
    for (tether, stether) in zip(sys.tethers, ssys.tethers)
        (length(stether.segment_idxs) != 1) && 
            error("Provide a simple system structure with 1-segment tethers.")

        # copy ground point of the tether
        point_idx = sys.segments[tether.segment_idxs[end]].point_idxs[2]
        spoint_idx = ssys.segments[stether.segment_idxs[1]].point_idxs[2]
        ssys.points[spoint_idx].pos_w .= sys.points[point_idx].pos_w
        ssys.points[spoint_idx].vel_w .= sys.points[point_idx].vel_w
    end

    # copy wing state
    swing = ssys.wings[1]
    wing = sys.wings[1]
    swing.pos_w .= wing.pos_w
    swing.vel_w .= wing.vel_w
    swing.ω_b .= wing.ω_b
    swing.Q_b_w .= wing.Q_b_w
    # update non-group pos
    ssys.points[1].pos_w .= wing.pos_w + wing.R_b_w * ssys.points[1].pos_b
    ssys.points[2].pos_w .= wing.pos_w + wing.R_b_w * ssys.points[2].pos_b

    # copy twist
    (length(sys.groups) != 4) && error("Sys should have 4 groups.")
    (length(ssys.groups) != 2) && error("Simple sys should have 2 groups.")
    ssys.groups[1].twist = (sys.groups[1].twist + sys.groups[2].twist) / 2
    ssys.groups[2].twist = (sys.groups[3].twist + sys.groups[4].twist) / 2
    ssys.groups[1].twist_ω = (sys.groups[1].twist_ω + sys.groups[2].twist_ω) / 2
    ssys.groups[2].twist_ω = (sys.groups[3].twist_ω + sys.groups[4].twist_ω) / 2

    # match moment by changing moment frac
    # TODO: add aero force
    moment = [group.moment for group in sys.groups]
    moment = [mean(moment[1:2]), mean(moment[3:4])]
    steering_force = [norm(sys.winches[2].force), norm(sys.winches[3].force)]
    for sgroup in ssys.groups
        x_airf = normalize(sgroup.chord)
        init_z_airf = x_airf × sgroup.y_airf
        z_airf = x_airf * sin(sgroup.twist) + init_z_airf * cos(sgroup.twist)
        force = steering_force[sgroup.idx] * normalize(swing.pos_w) ⋅ (swing.R_b_w * z_airf)
        r = moment[sgroup.idx] / force
        @show r steering_force[sgroup.idx] force
        spoint = ssys.points[sgroup.point_idxs[1]]
        spoint.pos_b .= sgroup.le_pos + sgroup.chord * (r / norm(sgroup.chord) + sgroup.moment_frac)

        # update pos_w for correct tether len
        chord_b = spoint.pos_b .- sgroup.le_pos
        normal = chord_b × sgroup.y_airf
        pos_b = sgroup.le_pos + cos(sgroup.twist) * chord_b - 
                sin(sgroup.twist) * normal
        spoint.pos_w .= swing.pos_w + swing.R_b_w * pos_b
    end

    # match winch force by changing tether length
    for (swinch, winch) in zip(ssys.winches, sys.winches)
        swinch.tether_len = 0.0
        for tether_idx in winch.tether_idxs
            stether = ssys.tethers[tether_idx]
            ssegment = ssys.segments[stether.segment_idxs[1]]
            spoint_idxs = ssegment.point_idxs
            slen = norm(ssys.points[spoint_idxs[1]].pos_w .-
                                ssys.points[spoint_idxs[2]].pos_w)
            stiffness = ssegment.axial_stiffness / slen
            nt = length(winch.tether_idxs)
            swinch.tether_len += (slen - norm(winch.force)/stiffness/nt) / nt
        end
        swinch.tether_vel = winch.tether_vel
    end
end

# Helper: Check if all values after index i are within p% band of steady-state
function in_percent_band(x, steady, delta_x, i, p)
    tol = p/100 * abs(delta_x)
    # All subsequent points must be within steady ± p%
    all(abs.(x[i:end] .- steady) .<= tol)
end

function calc_spring_props(sam::SymbolicAWEModel, tsam::SymbolicAWEModel; prn=false)
    find_steady_state!(sam; t=10.0, dt=3.0)
    copy!(sam.sys_struct, tsam.sys_struct)
    OrdinaryDiffEqCore.reinit!(tsam.integrator; reinit_dae=true)
    update_sys_struct!(tsam, tsam.sys_struct)

    F_0 = [-tsam.sys_struct.points[i].force for i in 1:4]
    steps = 200
    F_step = -0.1
    tether_lens = step(tsam, steps, F_step, F_0)
    k_values, c_values = calc_spring_props(sam, tether_lens, F_step; prn)
    return k_values .* tether_lens[:,1], c_values .* tether_lens[:,1]
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

