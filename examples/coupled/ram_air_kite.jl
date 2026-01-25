# Copyright (c) 2024, 2025 Bart van de Lint, Uwe Fechner
# SPDX-License-Identifier: MPL-2.0

"""
Ram Air Kite Simulation Example

This example uses RamAirKite.jl to run a simulation with sinusoidal steering.

Usage:
    julia --project=examples examples/coupled/ram_air_kite.jl
"""

using Timers
tic()
@info "Loading packages..."
using GLMakie
using RamAirKite

toc()

# Create simulation configuration
config = RamAirSimConfig(
    physical_model = "ram",      # Options: "ram", "simple_ram", "4_attach_ram"
    sim_time = 10.0,
    dt = 0.05,
    vsm_interval = 3,
    steering_freq = 0.5,         # Hz - full left-right cycle frequency
    steering_magnitude = 1.0,    # Nm
)

# Create and initialize model
sam = create_ram_air_model(config)
init!(sam; remake=false)
plot(sam.sys_struct)

# Find steady state
find_steady_state!(sam)

# Adjust bias for 4_attach_ram model
bias = config.steering_bias
if config.physical_model == "4_attach_ram"
    bias = 0.05
end

# Run oscillating simulation
sl, _ = sim_oscillate!(sam;
    dt = config.dt,
    total_time = config.sim_time,
    vsm_interval = config.vsm_interval,
    steering_freq = config.steering_freq,
    steering_magnitude = config.steering_magnitude,
    bias = bias,
    prn = true)

# Plot and replay
fig = plot(sam.sys_struct, sl)
scr = display(fig)
wait(scr)
replay(sl, sam.sys_struct)
