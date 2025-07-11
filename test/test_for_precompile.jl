# SPDX-FileCopyrightText: 2025 Bart van de Lint
# SPDX-License-Identifier: MPL-2.0

using SymbolicAWEModels, KiteUtils, LinearAlgebra, ControlPlots
sam_set = Settings("system_ram.yaml")
sam_set.segments = 3
set_values = [-50, 0.0, 0.0]  # Set values of the torques of the three winches. [Nm]
sam_set.quasi_static = false
sam_set.physical_model = "ram"
s = SymbolicAWEModel(sam_set)

# Initialize at elevation
SymbolicAWEModels.init!(s; prn=false, precompile=true)
find_steady_state!(s)
plot(s, 0.0)
steps = Int(round(10 / 0.05))
logger = Logger(length(s.sys_struct.points), steps)
sys_state = SysState(s)
next_step!(s)
simple_linearize!(s)

@info "Precompile script has completed execution."

