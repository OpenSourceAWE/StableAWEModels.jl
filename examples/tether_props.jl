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
k_values, c_values = SymbolicAWEModels.calc_spring_props(sam, tether_lens, F_step; prn=true)

display(plotx(
    dt .* collect(1:steps+1), 
    tether_lens[1,:].-tether_lens[1,1], 
    tether_lens[2,:].-tether_lens[2,1], 
    tether_lens[3,:].-tether_lens[3,1], 
    tether_lens[4,:].-tether_lens[4,1];
    title="Force step response",
    ylabels=["Left power [m]", "Right power [m]", "Left steering [m]", "Right steering [m]"],
))

set.segments = 1
ssys = SymbolicAWEModels.create_tether_sys_struct(set; 
                                                  axial_stiffness=k_values.*set.l_tether, 
                                                  axial_damping=c_values.*set.l_tether)
ssam = SymbolicAWEModel(set, ssys)
init!(ssam)

forces = [F ⋅ normalize(point.pos_w) for (F, point) in zip(F_0, ssys.points[1:4])]
SymbolicAWEModels.copy_to_simple!(sam.sys_struct, ssam.sys_struct)
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

