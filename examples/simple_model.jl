using SymbolicAWEModels, VortexStepMethod, KiteUtils, ControlPlots, Statistics
using OrdinaryDiffEqCore

set = Settings("system.yaml")

sam = SymbolicAWEModel(set, "ram")
init!(sam)
tsam = SymbolicAWEModel(set, "tether")
init!(tsam)

axial_stiffness, axial_damping = SymbolicAWEModels.calc_spring_props(sam, tsam)

ssam = SymbolicAWEModel(set, "simple_ram"; axial_stiffness, axial_damping)
init!(ssam)

sim_oscillate!(sam; total_time=1.0)

SymbolicAWEModels.copy_to_simple!(sam.sys_struct, ssam.sys_struct)
OrdinaryDiffEqCore.reinit!(ssam.integrator; reinit_dae=true)
SymbolicAWEModels.update_sys_struct!(ssam, ssam.sys_struct)

plot(sim_oscillate!(ssam))

