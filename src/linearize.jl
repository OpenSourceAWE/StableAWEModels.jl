# Copyright (c) 2025 Bart van de Lint
# SPDX-License-Identifier: MPL-2.0

function find_steady_state!(s::SymbolicAWEModel, integ=s.integrator; t=1.0, dt=1/s.set.sample_freq)
    old_state = s.get_stabilize(integ)
    s.set_stabilize(integ, true)
    for _ in 1:Int(round(t÷dt))
        next_step!(s; dt, vsm_interval=1)
    end
    s.set_stabilize(integ, old_state)
    update_sys_struct!(s, s.sys_struct)
    return nothing
end

function linearize_vsm!(s::SymbolicAWEModel, integ=s.integrator)
    wings = s.sys_struct.wings
    if length(wings) > 0
        vsm_y = s.get_vsm_y(integ)
        for wing in wings
            wing.vsm_y .= vsm_y[wing.idx,:]
            @show wing.vsm_y
            res = VortexStepMethod.linearize(
                s.vsm_solvers[wing.idx], 
                s.vsm_aeros[wing.idx], 
                wing.vsm_y;
                va_idxs=1:3, 
                omega_idxs=4:6,
                theta_idxs=7:6+length(s.sys_struct.groups),
                moment_frac=s.sys_struct.groups[1].moment_frac
            )
            wing.vsm_jac .= res[1]
            wing.vsm_x .= res[2]
        end
        s.set_psys(integ, s.sys_struct)
    end
    nothing
end

function linearize!(s::SymbolicAWEModel; set_values=s.get_set_values(s.integrator))
    isnothing(s.lin_prob) && error("Run init! with remake=true and lin_outputs=...")
    s.set_lin_vsm(s.lin_prob, s.get_vsm(s.integrator))
    s.set_lin_set_values(s.lin_prob, set_values)
    s.set_lin_unknowns(s.lin_prob, s.get_unknowns(s.integrator))
    s.lin_model = solve(s.lin_prob)[1]
    return s.lin_model
end

function getstate(sys_struct::SystemStructure)
    @unpack wings, winches = sys_struct
    c = copy
    wing = wings[1]
    tether_len = [winch.tether_len for winch in winches]
    tether_vel = [winch.tether_vel for winch in winches]
    return (c(wing.pos_w), c(wing.vel_w), c(wing.wind_disturb), c(wing.R_b_w), c(wing.ω_b),
                tether_len, tether_vel)
end

function setstate!(sys_struct::SystemStructure, state)
    @unpack wings, winches = sys_struct
    wing = wings[1]
    pos_w, vel_w, wind_disturb, R_b_w, ω_b, tether_len, tether_vel = state
    wing.pos_w .= pos_w
    wing.vel_w .= vel_w
    wing.wind_disturb .= wind_disturb
    wing.R_b_w .= R_b_w
    wing.ω_b .= ω_b
    for winch in winches
        winch.tether_len = tether_len[winch.idx]
        winch.tether_vel = tether_vel[winch.idx]
    end
end

function set_measured!(sys_struct::SystemStructure, 
    heading, turn_rate,
    tether_len, tether_vel
)
    @unpack wings, winches = sys_struct
    wing = wings[1]

    # get variables from integrator
    distance = norm(wing.pos_w)
    R_t_w = calc_R_t_w(wing.elevation, wing.azimuth) # rotation of tether to world, similar to view rotation, but always pointing up
    R_v_w = calc_R_v_w(wing.pos_w, wing.R_b_w[:,1])
    
    # get wing_pos, rotate it by elevation and azimuth around the x and z axis
    wing.pos_w .= R_t_w * [0, 0, distance + tether_len[1] - winches[1].tether_len]
    # wing_vel from elevation_vel and azimuth_vel
    # TODO: now I uderstand, vel should not be stabilized, only pos things
    wing.vel_w .= R_t_w * [-wing.elevation_vel, wing.azimuth_vel, 0.0]
    wing.wind_disturb .= R_t_w * [0.0, 0.0, -tether_vel[1]]
    # find quaternion orientation from heading, R_b_w and R_t_w
    R_b_w = zeros(3,3)
    cur_heading = calc_heading(R_t_w, R_v_w)
    d_heading = heading - cur_heading
    for i in 1:3
        R_b_w[:,i] .= R_t_w * rotate_around_z(R_t_w' * wing.R_b_w[:,i], d_heading)
    end
    wing.R_b_w = R_b_w
    # adjust the turn rates for observed turn rate
    wing.ω_b .= wing.R_b_w' * R_t_w * [wing.turn_rate[1], wing.turn_rate[2], turn_rate]
    # directly set tether length
    for winch in winches
        winch.tether_len = tether_len[winch.idx]
        winch.tether_vel = tether_vel[winch.idx]
    end
    return nothing
end

function jacobian(f::Function, x::AbstractVector, ϵ::AbstractVector)
    n = length(x)
    fx = f(x)
    m = length(fx)
    J = zeros(m, n)
    for i in 1:n
        x_perturbed = copy(x)
        x_perturbed[i] += ϵ[i]
        J[:, i] = (f(x_perturbed) - fx) / ϵ[i]
    end
    return J
end

function simple_linearize!(s::SymbolicAWEModel; tstab=10.0)
    integ = s.integrator
    update_sys_struct!(s, s.sys_struct)
    state0 = getstate(s.sys_struct)
    old_stab = s.get_stabilize(integ)
    s.set_stabilize(integ, true)
    lin_x0 = s.get_lin_x(integ)
    u0 = [winch.set_value for winch in s.sys_struct.winches]
    @unpack A, B, C, D = s.simple_lin_model
    A .= 0.0
    B .= 0.0
    C .= 0.0
    D .= 0.0

    # TODO: add sparsity pattern for the known zeros
    function f(x, u)
        heading = x[1]
        turn_rate = x[2]
        tether_len = x[3:5]
        tether_vel = x[6:8]
        set_measured!(s.sys_struct, heading, turn_rate,
                      tether_len, tether_vel)
        s.set_set_values(integ, u)
        OrdinaryDiffEqCore.reinit!(integ)
        OrdinaryDiffEqCore.step!(integ, tstab)
        return s.get_lin_dx(integ)
    end

    # yes it looks weird to step in an output function, but this is a steady state finder rather than output
    function h(x, u)
        heading = x[1]
        turn_rate = x[2]
        tether_len = x[3:5]
        tether_vel = x[6:8]
        set_measured!(s.sys_struct, heading, turn_rate,
                      tether_len, tether_vel)
        s.set_set_values(integ, u)
        OrdinaryDiffEqCore.reinit!(integ)
        OrdinaryDiffEqCore.step!(integ, tstab)
        return s.get_lin_y(integ)
    end

    f_x(x) = f(x, u0)
    f_u(u) = f(lin_x0, u)
    h_x(x) = h(x, u0)
    h_u(u) = h(lin_x0, u)

    @unpack segments, winches, tethers = s.sys_struct
    segment = segments[tethers[1].segment_idxs[1]]
    mass_per_meter = s.set.rho_tether * π * (segment.diameter/2)^2
    mass = winches[1].tether_len * mass_per_meter + s.set.mass

    # calculate jacobian
    ϵ_x = [0.001, 0.1, 0.001, 0.001, 0.001, 0.1, 0.1, 0.1]
    ϵ_u = [1.0, 0.1, 0.1]
    A .= jacobian(f_x, lin_x0, ϵ_x)
    B .= jacobian(f_u, u0, ϵ_u)
    C .= jacobian(h_x, lin_x0, ϵ_x)
    D .= 0.0
    D[4,1] = -mass * B[6,1]
    A[:,1] .= 0.0 # Aero moment due to change in heading cannot be found in steady state
    C[4,1] = 0.0
    s.set_set_values(integ, u0)
    s.set_stabilize(integ, old_stab)
    setstate!(s.sys_struct, state0)
    OrdinaryDiffEqCore.reinit!(integ)
    return s.simple_lin_model
end

