# Copyright (c) 2024, 2025 Bart van de Lint, Uwe Fechner
# SPDX-License-Identifier: MIT

using Timers
using Pkg
if ! ("LaTeXStrings" ∈ keys(Pkg.project().dependencies))
    using TestEnv; TestEnv.activate()
end
using ControlPlots, LaTeXStrings
using KiteModels, LinearAlgebra, Statistics

# Initialize model
set = load_settings("system_ram.yaml")
set.segments = 3
set_values = [-50, 0.0, 0.0]  # Set values of the torques of the three winches. [Nm]
set.quasi_static = false
set.physical_model = "ram"
set.elevation = 45
set.azimuth = 0.0
set.heading = 0.0

@info "Creating wing, aero, vsm_solver, sys_struct and symbolic_awe_model:"
wing = RamAirWing(set; prn=false)
sys_struct = SystemStructure(set, wing)

init!(sys_struct, set)
init!(sys_struct, set)
plot(sys_struct, 0.0; zoom=false, front=true)
transform = sys_struct.transforms[1]

