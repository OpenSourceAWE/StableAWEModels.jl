# Copyright (c) 2025 Bart van de Lint, Uwe Fechner
# SPDX-License-Identifier: MPL-2.0

using Timers
tic()
@info "Loading packages "
using GLMakie
using SymbolicAWEModels

PLOT = false
using ControlSystemsBase

if ! @isdefined SIMPLE
    SIMPLE = false
end

toc()

# Simulation parameters
dt = 0.05
total_time = 1.0  # Longer simulation to see oscillations
vsm_interval = 3
steps = Int(round(total_time / dt))

# Steering parameters
steering_freq = 1/2  # Hz - full left-right cycle frequency
steering_magnitude = 10.0      # Magnitude of steering input [Nm]

# Initialize model
set = Settings("ram_air_kite/system.yaml")
set_values = [-50, 0.0, 0.0]  # Set values of the torques of the three winches. [Nm]

@info "Creating wing, aero, vsm_solver, sys_struct and symbolic_awe_model:"
sam = SymbolicAWEModel(set)
sam.set.abs_tol = 1e-3
sam.set.rel_tol = 1e-3
toc()

# Initialize at elevation
set.l_tethers[2] += 0.4
set.l_tethers[3] += 0.4
set.elevation = 70.0
init!(sam; remake=false, reload=false)
sys = sam.sys

@info "System initialized at:"
toc()

# Stabilize system
find_steady_state!(sam; t=10.0, dt=1.0)
u0 = -sam.set.drum_radius .* sam.integrator[sys.winch_force]
sam.set_set_values(sam.integrator, u0)
simple_linearize!(sam; tstab=10.0)
lin_sam = ss(sam.simple_lin_model.A, 
             sam.simple_lin_model.B, 
             sam.simple_lin_model.C, 
             sam.simple_lin_model.D)

logger = Logger(length(sam.sys_struct.points), steps)
sys_state = SysState(sam)
t = 0.0
runtime = 0.0
integ_runtime = 0.0
bias = set.quasi_static ? 0.45 : 0.35
t0 = sam.integrator.t
set_values = zeros(3, steps)

try
    for i in 1:steps
        local steering
        global t, set_values, runtime, integ_runtime
        PLOT && plot(sam, t)
        
        # Calculate steering inputs based on cosine wave
        sign = t > 0.5 ? -1 : 1
        set_values[:,i] = -sam.set.drum_radius .* sam.integrator[sys.winch_force]
        set_values[:,i] .+= sign .* [10.0, steering_magnitude, -steering_magnitude]  # Opposite steering for left/right
        _vsm_interval = vsm_interval
        # Step simulation
        steptime = @elapsed next_step!(sam; set_values=set_values[:,i], 
            dt, vsm_interval=vsm_interval)
        t_new = sam.integrator.t
        integ_steptime = sam.t_step
        t = t_new - t0  # Adjust for initial stabilization time

        # Track performance after initial transient
        if (t > total_time/2)
            runtime += steptime
            integ_runtime += integ_steptime
        end

        # Log state variables
        update_sys_state!(sys_state, sam)
        sys_state.time = t
        log!(logger, sys_state)
    end
catch e
    if isa(e, AssertionError)
        @show t
        println(e)
    else
        rethrow(e)
    end
end
@info "Total time without plotting:"
toc()

# Save results
save_log(logger, "tmp")
lg = load_log("tmp")
sl = lg.syslog

# Linear simulation
lin_res = lsim(lin_sam, set_values .- u0, sl.time)

@info "Simulation completed" steps=length(sl.time) final_heading=rad2deg(sl.heading[end])

@info "Performance:" times_realtime=(total_time/2)/runtime integrator_times_realtime=(total_time/2)/integ_runtime

# 55x realtime (PLOT=false, CPU: Intel i9-9980HK (16) @ 5.000GHz)
# 40-65x realtime (PLOT=false, CPU: Intel i9-9980HK (16) @ 5.000GHz) - commit 6620ed5d0a38e96930615aad9a66e4cd666955f2
# 40x realtime (PLOT=false, CPU: Intel i9-9980HK (16) @ 5.000GHz) - commit 88a78894038d3cbd50fbff83dfbe5c26266b0637
