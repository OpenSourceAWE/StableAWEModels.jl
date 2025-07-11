# Copyright (c) 2024, 2025 Bart van de Lint, Uwe Fechner
# SPDX-License-Identifier: MIT

using Timers
tic()
@info "Loading packages "

using SymbolicAWEModels, KiteUtils
toc()

@info "Creating wing, aero, vsm_solver, sys_struct and s:"
set = Settings("system_ram.yaml")
set_values = [-50, 0.0, 0.0]  # Set values of the torques of the three winches. [Nm]
s = SymbolicAWEModel(set)
SymbolicAWEModels.init!(s; remake=false) # doesn't remake the model if it exists
@info "System initialized at:"
toc()

