# Copyright (c) 2025 Bart van de Lint
# SPDX-License-Identifier: MPL-2.0

using GLMakie
using SymbolicAWEModels, KiteUtils, LinearAlgebra, OrdinaryDiffEqBDF

# Initialize the model
set = Settings("system.yaml")
set.sample_freq = 100
sam = SymbolicAWEModel(set)
init!(sam)
find_steady_state!(sam)

# Run the simulation with repositioning
lg = SymbolicAWEModels.sim_reposition!(
    sam,
    total_time=10.0,
    reposition_interval_s=0.5,
    target_elevation=deg2rad(50.0),
    target_azimuth=deg2rad(10.0),
    target_heading=deg2rad(10.0)
)

# Plot the results
plot(sam.sys_struct, lg)

