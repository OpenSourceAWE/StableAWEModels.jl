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

function step(sam, steps, F_step, F_0;
                  abs_tol=1e-6,
                  consecutive_steps_needed=10,
                  stop_at_peak=false)

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

tether_lens = step(tsam, steps, F_step, F_0)

display(plotx(
    dt .* collect(1:steps+1), 
    tether_lens[1,:], tether_lens[2,:], tether_lens[3,:], tether_lens[4,:];
    title="Force step step",
    ylabels=["Left power [m]", "Right power [m]", "Left steering [m]", "Right steering [m]"],
))

# Damping ratio from overshoot (PO = fractional overshoot, i.e. (peak - final)/|final - initial|)
function damping_ratio_from_PO(PO)
    if PO <= 0 || PO >= 1
        error("PO must be between 0 and 1 (exclusive) for this method.")
    end
    ζ = -log(PO) / sqrt(π^2 + log(PO)^2)
    return ζ
end

k_values = zeros(4)
c_values = zeros(4)

diameters = [
    set.power_tether_diameter/1000
    set.power_tether_diameter/1000
    set.steering_tether_diameter/1000
    set.steering_tether_diameter/1000
]
mass_per_meter = set.rho_tether * π * ((diameters/2).^2)

for j in 1:4
    tether_len_series = tether_lens[j, :]
    initial_len = tether_len_series[1]
    final_len = tether_len_series[end]
    delta_x_ss = final_len - initial_len

    # Effective mass approximation (as before)
    m = mass_per_meter[j] * 0.5 * set.l_tether

    if abs(delta_x_ss) < 1e-6
        println("Warning: Steady-state change too small for Tether $j; skipping.")
        k_values[j] = NaN; c_values[j] = NaN
        continue
    end

    # Spring stiffness
    k = F_step / delta_x_ss
    k_values[j] = k

    # Determine peak for overshoot
    peak_val = delta_x_ss > 0 ? maximum(tether_len_series) : minimum(tether_len_series)
    PO = abs((peak_val - final_len) / delta_x_ss)
    ζ = damping_ratio_from_PO(PO)
    c = 2 * ζ * sqrt(k * m)
    c_values[j] = c
end

println("Summary of Results:")
for j in 1:4
    println("Tether $(j): k = $(k_values[j]) N/m, c = $(c_values[j]) Ns/m")
end

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

tether_lens = step(ssam, steps, 0.0,    F_0)
tether_lens = step(ssam, steps, F_step, F_0)

display(plotx(
    dt .* collect(1:steps+1), 
    tether_lens[1,:], tether_lens[2,:], tether_lens[3,:], tether_lens[4,:];
    title="Force step step",
    ylabels=["Left power [m]", "Right power [m]", "Left steering [m]", "Right steering [m]"],
))

