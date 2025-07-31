# SPDX-FileCopyrightText: 2025 Bart van de Lint <bart@vandelint.net>
#
# SPDX-License-Identifier: MPL-2.0

using SymbolicAWEModels, VortexStepMethod, KiteUtils, ControlPlots, Statistics
using OrdinaryDiffEqCore

set = Settings("system.yaml")
sam = SymbolicAWEModel(set, "ram")
init!(sam)

tset = Settings("system.yaml")
tsam = SymbolicAWEModel(tset, "tether")
init!(tsam)

sset = Settings("system.yaml")
ssam = SymbolicAWEModel(sset, "simple_ram")
init!(ssam)

sim_oscillate!(sam; total_time=1.0)
SymbolicAWEModels.copy_to_simple!(sam, tsam, ssam)

bias = 0.2
sl = sim_oscillate!(sam; total_time=5.0, prn=true, bias) # TODO: add first frac ram model
display(plot(sam.sys_struct, sl))

sl = sim_oscillate!(ssam; total_time=5.0, prn=true, bias)
display(plot(ssam.sys_struct, sl))

