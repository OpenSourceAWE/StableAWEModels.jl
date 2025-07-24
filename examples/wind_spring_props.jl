using SymbolicAWEModels, VortexStepMethod, KiteUtils, WinchModels
using ControlPlots, Statistics, LinearAlgebra
using OrdinaryDiffEqCore
using UnPack

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

tether_lens = SymbolicAWEModels.step(tsam, steps, F_step, F_0)

display(plotx(
    dt .* collect(1:steps+1), 
    tether_lens[1,:].-tether_lens[1,1], 
    tether_lens[2,:].-tether_lens[2,1], 
    tether_lens[3,:].-tether_lens[3,1], 
    tether_lens[4,:].-tether_lens[4,1];
    title="Force step response",
    ylabels=["Left power [m]", "Right power [m]", "Left steering [m]", "Right steering [m]"],
))

# Helper: Check if all values after index i are within 2% band of steady-state
function in_percent_band(x, steady, delta_x, i, p)
    tol = p/100 * abs(delta_x)
    # All subsequent points must be within steady ± 2%
    all(abs.(x[i:end] .- steady) .<= tol)
end

function calc_spring_props(sam; p=5, prn=false)
    @unpack tethers, segments = sam.sys_struct

    k_values = zeros(4)
    c_values = zeros(4)

    diameters = [segments[tether.segment_idxs[1]].diameter for tether in tethers]
    mass_per_meter = set.rho_tether * π * ((diameters/2).^2)

    for j in eachindex(tethers)
        tether_len_series = tether_lens[j, :]
        initial_len = tether_len_series[1]
        final_len = tether_len_series[end]
        delta_x_ss = final_len - initial_len

        # Effective mass approximation (as before)
        m = mass_per_meter[j] * 0.5 * set.l_tether

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

k_values, c_values = calc_spring_props(sam)

set.segments = 1
ssys = SymbolicAWEModels.create_tether_sys_struct(set; 
                                                  axial_stiffness=k_values.*set.l_tether, 
                                                  axial_damping=c_values.*set.l_tether)
ssam = SymbolicAWEModel(set, ssys)
init!(ssam)

forces = [F ⋅ normalize(point.pos_w) for (F, point) in zip(F_0, ssys.points[1:4])]
SymbolicAWEModels.copy!(sam.sys_struct, ssam.sys_struct)
OrdinaryDiffEqCore.reinit!(ssam.integrator; reinit_dae=true)
SymbolicAWEModels.update_sys_struct!(ssam, ssam.sys_struct)

stether_lens = SymbolicAWEModels.step(ssam, steps, 0.0,    F_0)
stether_lens = SymbolicAWEModels.step(ssam, steps, F_step, F_0)

display(plotx(
    dt .* collect(1:steps+1), 
    stether_lens[1,:].-stether_lens[1,1], 
    stether_lens[2,:].-stether_lens[2,1], 
    stether_lens[3,:].-stether_lens[3,1], 
    stether_lens[4,:].-stether_lens[4,1];
    title="Force step response",
    ylabels=["Left power [m]", "Right power [m]", "Left steering [m]", "Right steering [m]"],
))

@info "Difference at t=0: $(stether_lens[:,1] .- tether_lens[:,1])"

