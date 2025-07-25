using SymbolicAWEModels, VortexStepMethod, KiteUtils, ControlPlots, Statistics
using OrdinaryDiffEqCore

set = Settings("system.yaml")

sam = SymbolicAWEModel(set, "ram")
init!(sam)
ssam = SymbolicAWEModel(set, "simple_ram")
init!(ssam)
tsam = SymbolicAWEModel(set, "tether")
init!(tsam)

find_steady_state!(sam; t=10.0, dt=3.0)
SymbolicAWEModels.copy!(sam.sys_struct, tsam.sys_struct)
OrdinaryDiffEqCore.reinit!(tsam.integrator; reinit_dae=true)
SymbolicAWEModels.update_sys_struct!(tsam, tsam.sys_struct)

SymbolicAWEModels.copy_to_simple!(sam.sys_struct, ssam.sys_struct)

plot(sim_oscillate(sam))

