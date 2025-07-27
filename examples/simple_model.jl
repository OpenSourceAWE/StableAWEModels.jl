# SPDX-FileCopyrightText: 2025 Bart van de Lint <bart@vandelint.net>
#
# SPDX-License-Identifier: MPL-2.0

using SymbolicAWEModels, VortexStepMethod, KiteUtils, ControlPlots, Statistics
using OrdinaryDiffEqCore

set = Settings("system.yaml")

sam = SymbolicAWEModel(set, "ram")
init!(sam)
tsam = SymbolicAWEModel(set, "tether")
init!(tsam)

axial_stiffness, axial_damping = SymbolicAWEModels.calc_spring_props(sam, tsam; prn=true)

ssam = SymbolicAWEModel(set, "simple_ram"; axial_stiffness, axial_damping)
init!(ssam)

sim_oscillate!(sam; total_time=1.0)

SymbolicAWEModels.copy_to_simple!(sam.sys_struct, ssam.sys_struct)
OrdinaryDiffEqCore.reinit!(ssam.integrator; reinit_dae=true)
SymbolicAWEModels.update_sys_struct!(ssam, ssam.sys_struct)

sl = sim_oscillate!(sam; total_time=5.0, prn=true) # TODO: add first frac ram model
display(plot(sam.sys_struct, sl))

sl = sim_oscillate!(ssam; total_time=5.0, prn=true)
display(plot(ssam.sys_struct, sl))

