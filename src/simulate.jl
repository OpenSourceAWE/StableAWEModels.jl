# Copyright (c) 2024, 2025 Bart van de Lint, Uwe Fechner
# SPDX-License-Identifier: MIT

"""
    sim_oscillate(set, sam, dt, total_time; steering_freq=0.5, steering_magnitude=10.0, vsm_interval=3)

Run a simulation with oscillating steering input on the given AWE model.

# Arguments
- `sam::SymbolicAWEModel`: Initialized AWE model.

# Keywords
- `dt::Float64`: Time step [s].
- `total_time::Float64`: Simulation duration [s].
- `steering_freq`: Steering oscillation frequency [Hz] (default 0.5).
- `steering_magnitude`: Steering torque amplitude [Nm] (default 10.0).
- `vsm_interval`: Value state machine interval (default 3).

# Returns
- `SysLog`: Logged simulation data.
"""
function sim_oscillate!(
    sam::SymbolicAWEModel;
    dt=1/sam.set.sample_freq,
    total_time=10.0,
    steering_freq=0.5,
    steering_magnitude=10.0,
    vsm_interval=3,
    bias = sam.set.quasi_static ? 0.45 : 0.35
)
    steps = Int(round(total_time / dt))
    logger = Logger(length(sam.sys_struct.points), steps)
    sys_state = SysState(sam)

    for step in 1:steps
        t = (step-1)*dt
        steering = steering_magnitude * cos(2π * steering_freq * t + bias)
        set_values = -sam.set.drum_radius .* [norm(winch.force) for winch in sam.sys_struct.winches]
        if t > 1.0
            set_values .+= [0.0, steering, -steering]
        end
        try
            next_step!(sam; set_values, dt, vsm_interval=vsm_interval)
        catch e
            if e isa AssertionError
                @warn "Crashed at t=$t"
                break
            else
                rethrow(e)
            end
        end
        update_sys_state!(sys_state, sam)
        sys_state.time = t
        log!(logger, sys_state)
    end

    save_log(logger, "tmp_run")
    lg = load_log("tmp_run")
    return lg
end
