# Copyright (c) 2025 Bart van de Lint and Uwe Fechner
# SPDX-License-Identifier: MPL-2.0

function generate_prob_getters(sys_struct, sys)
    c = collect
    @unpack wings, groups, pulleys, winches, tethers, segments = sys_struct
    get_wing_state, get_vsm_y, get_segment_state, get_group_state, get_pulley_state,
    get_winch_state, get_tether_state, set_set_values, get_set_values = ntuple(i -> nothing, 9)

    if length(wings) > 0
        wing_vars = c.([
            sys.Q_b_w, sys.ω_b, sys.wing_pos, sys.wing_vel, sys.wing_acc, sys.va_wing_b,
            sys.wind_vel_wing, sys.aero_force_b, sys.aero_moment_b, sys.elevation,
            sys.elevation_vel, sys.elevation_acc, sys.azimuth, sys.azimuth_vel,
            sys.azimuth_acc, sys.heading, sys.turn_rate, sys.turn_acc, sys.course,
            sys.angle_of_attack
        ])
        get_wing_state = getu(sys, wing_vars)
        get_vsm_y = getu(sys, sys.y)
    end
    if length(segments) > 0; get_segment_state = getu(sys, c.([sys.spring_force, sys.len])); end
    if length(groups) > 0; get_group_state = getu(sys, c.([sys.twist_angle, sys.twist_ω, sys.group_tether_force, sys.group_tether_moment, sys.group_aero_moment])); end
    if length(pulleys) > 0; get_pulley_state = getu(sys, c.([sys.pulley_len, sys.pulley_vel])); end
    if length(winches) > 0
        get_winch_state = getu(sys, c.([sys.tether_len, sys.tether_vel, sys.set_values, sys.winch_force_vec]))
        set_set_values = setp(sys, sys.set_values)
        get_set_values = getp(sys, sys.set_values)
    end
    if length(tethers) > 0; get_tether_state = getu(sys, c(sys.stretched_len)); end
    set_sys = setp(sys, sys.psys)
    set_set = setp(sys, sys.pset)
    get_struct_state = getu(sys, sys.wind_vec_gnd)
    get_point_state = getu(sys, c.([sys.pos, sys.vel, sys.point_force]))
    return (; get_wing_state, get_vsm_y, get_segment_state, get_group_state,
              get_pulley_state, get_winch_state, get_tether_state, set_set_values,
              get_set_values, set_sys, set_set, get_struct_state, get_point_state)
end

function generate_simple_lin_model(sys_struct, sys, lin_y_vec)
    @unpack wings, winches = sys_struct
    if length(wings) == 1
        lin_x_vec = [
            sys.heading[1], sys.turn_rate[1,3],
            sys.tether_len[1], sys.tether_len[2], sys.tether_len[3],
            sys.tether_vel[1], sys.tether_vel[2], sys.tether_vel[3]
        ]
        lin_dx_vec = [
            sys.turn_rate[1,3], sys.turn_acc[1,3],
            sys.tether_vel[1], sys.tether_vel[2], sys.tether_vel[3],
            sys.tether_acc[1], sys.tether_acc[2], sys.tether_acc[3]
        ]
        get_lin_x = getu(sys, lin_x_vec)
        get_lin_dx = getu(sys, lin_dx_vec)
        get_lin_y = getu(sys, lin_y_vec)
        nx, ny, nu = length(lin_x_vec), length(lin_y_vec), length(winches)
        simple_lin_model = (; A=zeros(nx,nx), B=zeros(nx,nu), C=zeros(ny,nx), D=zeros(ny,nu))
        return (; simple_lin_model, get_lin_x, get_lin_dx, get_lin_y)
    end
    return (simple_lin_model=nothing, get_lin_x=nothing, get_lin_dx=nothing, get_lin_y=nothing)
end

function generate_lin_getters(sys)
    set_lin_set_values = nothing
    if hasproperty(sys, :set_values)
        set_lin_set_values = setp(sys, sys.set_values)
    end
    set_lin_sys = setp(sys, sys.psys)
    set_lin_set = setp(sys, sys.pset)
    return (set_lin_set_values=set_lin_set_values, set_lin_sys=set_lin_sys, set_lin_set=set_lin_set)
end

function generate_control_funcs(model, inputs, outputs)
    (f_ip, f_oop), dvs, psym, io_sys = @suppress_err ModelingToolkit.generate_control_function(
        model, inputs; simplify=false
    )
    nu, nx, ny = length(inputs), length(dvs), length(outputs)
    (h_oop, h_ip) = ModelingToolkit.build_explicit_observed_function(
        io_sys, outputs; inputs, return_inplace=true
    )
    return (f_oop=f_oop, f_ip=f_ip, h_oop=h_oop, h_ip=h_ip, nu=nu, nx=nx, ny=ny,
            dvs=dvs, psym=psym, io_sys=io_sys)
end

function load_serialized_model!(sam, model_path; remake=false, reload=false)
    set_hash = get_set_hash(sam.set)
    sys_struct_hash = get_sys_struct_hash(sam.sys_struct)
    
    if ispath(model_path) && !remake
        if set_hash != sam.serialized_model.set_hash || sys_struct_hash != sam.serialized_model.sys_struct_hash
            sam.serialized_model = SerializedModel(; set_hash, sys_struct_hash)
            return false
        end
        if !isnothing(sam.serialized_model.full_sys) && !reload
            return true
        end
        try
            serialized_model = deserialize(model_path)
            if set_hash == serialized_model.set_hash && sys_struct_hash == serialized_model.sys_struct_hash
                sam.serialized_model = serialized_model
                return true
            end
        catch e
            @warn "Failure to deserialize $model_path: $(typeof(e))"
        end
    end
    sam.serialized_model = SerializedModel(; set_hash, sys_struct_hash)
    return false
end

function maybe_create_prob!(sam; create_prob=true, prn=true)
    if create_prob && isnothing(sam.prob)
        local sys
        time = @elapsed @suppress_err sys = mtkcompile(sam.full_sys; inputs=sam.inputs)
        prn && println("\tSimplified the System for ODEProblem in $time seconds.")

        dt = SimFloat(1/sam.set.sample_freq)
        time = @elapsed prob = ODEProblem(sys, sam.defaults, (0.0, dt); u0_prior=sam.guesses)
        prn && println("\tCreated the ODEProblem in $time seconds.")

        time = @elapsed getters = generate_prob_getters(sam.sys_struct, sys)
        prn && println("\tCreated state getters and setters in $time seconds.")

        sam.prob = ProbWithAttributes(; prob, getters...)
        return true
    end
    return false
end

function maybe_create_simple_lin_model!(sam, lin_outputs; create_simple_lin_model=true, prn=true)
    if create_simple_lin_model && isnothing(sam.simple_lin_model)
        sys = sam.prob.sys
        time = @elapsed slm_attrs = generate_simple_lin_model(sam.sys_struct, sys, lin_outputs)
        if !isnothing(slm_attrs.simple_lin_model)
            sam.simple_lin_model = SimpleLinModelWithAttributes(; slm_attrs...)
        end
        prn && println("\tCreated simplified linear model in $time seconds.")
        return true
    end
    return false
end

function maybe_create_lin_prob!(sam, lin_outputs; create_lin_prob=true, prn=true)
    if create_lin_prob && (isnothing(sam.lin_prob) ||
           length(lin_outputs) != length(sam.lin_prob.lin_outputs) ||
           !all(string.(lin_outputs) .== string.(sam.lin_prob.lin_outputs)))
        time = @elapsed @suppress_err begin
            lin_fun, lin_sys = linearization_function(sam.full_sys, [sam.inputs...], lin_outputs;
                                                      op=sam.defaults, guesses=sam.guesses)
            lin_prob = LinearizationProblem(lin_fun, 0.0)
            getters = generate_lin_getters(lin_sys)
            sam.lin_prob = LinProbWithAttributes(; lin_prob,
                                                 lin_outputs,
                                                 getters...)
        end
        prn && println("\tCreated the LinearizationProblem in $time seconds.")
        return true
    end
    return false
end

function maybe_create_control_functions!(sam, lin_outputs; create_control_func=false, prn=true)
    if create_control_func && (isnothing(sam.control_funcs) ||
           length(lin_outputs) != length(sam.lin_prob.lin_outputs) ||
           !all(string.(lin_outputs) .== string.(sam.lin_prob.lin_outputs)))
        inputs = [sam.inputs...]
        time = @elapsed result = generate_control_funcs(sam.full_sys, inputs, lin_outputs)
        sam.control_funcs = ControlFuncWithAttributes(; result...)
        prn && println("\tCreated the control functions in $time seconds.")
        return true
    end
    return false
end

function init!(sam::SymbolicAWEModel;
    solver=nothing, adaptive=true, prn=true,
    remake=false, reload=false,
    lin_outputs=nothing,
    create_prob::Bool=true,
    create_simple_lin_model::Bool=true,
    create_lin_prob::Bool=true,
    create_control_func::Bool=false,
    lin_vsm::Bool=true
)
    prn && @info "Initializing $(sam.sys_struct.name) model..."
    time = @elapsed begin
        if isnothing(solver)
            solver = if sam.set.solver == "FBDF"
                sam.set.quasi_static ? FBDF(nlsolve=OrdinaryDiffEqNonlinearSolve.NLNewton(relax=sam.set.relaxation)) : FBDF()
            elseif sam.set.solver == "QNDF"
                @warn "This solver is not tested."
                QNDF()
            else
                error("Unavailable solver for SymbolicAWEModel: $(sam.set.solver).")
            end
        end

        if isnothing(lin_outputs)
            @variables begin
                heading(t)[1]
                angle_of_attack(t)[1]
                tether_len(t)[1:3]
                winch_force(t)[1:3]
            end
            lin_outputs = Num[]
            if length(sam.sys_struct.wings) > 0
                push!(lin_outputs, heading[1], angle_of_attack[1])
            end
            if length(sam.sys_struct.winches) > 0
                push!(lin_outputs, tether_len[1], winch_force[1])
            end
        end

        model_path = joinpath(KiteUtils.get_data_path(), get_model_name(sam.set))
        loaded = load_serialized_model!(sam, model_path; remake, reload)
        changed = false
        if !loaded
            sam.inputs = create_sys!(sam, sam.sys_struct; prn)
            changed = true
        end
        
        changed |= maybe_create_prob!(sam; create_prob, prn)
        changed |= maybe_create_simple_lin_model!(sam, lin_outputs; create_simple_lin_model, prn)
        changed |= maybe_create_lin_prob!(sam, lin_outputs; create_lin_prob, prn)
        changed |= maybe_create_control_functions!(sam, lin_outputs; create_control_func, prn)

        if changed
            prn && @info "Serializing model to $model_path..."
            serialize(model_path, sam.serialized_model)
        end

        reinit!(sam.sys_struct, sam.set)
        create_prob && !isnothing(sam.prob) && reinit!(sam, sam.prob, solver; adaptive, reload, lin_vsm)
        create_lin_prob && !isnothing(sam.lin_prob) && reinit!(sam, sam.lin_prob)
    end
    prn && @info "$(sam.sys_struct.name) model initialized in $time seconds."
    return sam.integrator
end

"""
    reinit!(sam::SymbolicAWEModel, lin_prob::ModelingToolkit.LinearizationProblem)

Reinitializes a `LinearizationProblem` with the current system and settings parameters.

This function updates the internal parameter vectors of the linearization problem
with the latest values from the `SymbolicAWEModel`'s `sys_struct` and `set` fields.
"""
function reinit!(sam::SymbolicAWEModel, prob::LinProbWithAttributes)
    prob.set_sys(prob.prob, sam.sys_struct)
    prob.set_set(prob.prob, sam.set)
    nothing
end

"""
    reinit!(s::SymbolicAWEModel, prob::ODEProblem, solver; prn, precompile, reload, lin_outputs) -> (ODEIntegrator, Bool)

Reinitializes an existing kite power system model's ODE integrator.

This function resets the integrator's state with new values from `s.set`,
allowing for the simulation to be restarted from a new initial condition
without needing to rebuild the entire symbolic model.

# Arguments
- `s::SymbolicAWEModel`: The kite power system state object.
- `prob::ODEProblem`: The ODE problem to be solved.
- `solver`: The solver to be used.

# Keyword Arguments
- `adaptive::Bool=true`: Whether to use adaptive time-stepping.
- `reload::Bool=true`: Force reloading the model from disk.
- `lin_vsm::Bool=true`: If `true`, linearizes the VSM model after reinitialization.

# Returns
- `(ODEIntegrator, Bool)`: A tuple containing the reinitialized integrator and a success flag.
"""
function reinit!(
    sam::SymbolicAWEModel,
    prob::ProbWithAttributes,
    solver;
    adaptive=true,
    reload=true, 
    lin_vsm=true
)
    if isnothing(sam.integrator) || !successful_retcode(sam.integrator.sol) || reload
        dt = SimFloat(1/sam.set.sample_freq)
        sam.integrator = OrdinaryDiffEqCore.init(prob.prob, solver; 
            adaptive, dt, tspan=(0.0, dt), abstol=sam.set.abs_tol, reltol=sam.set.rel_tol, 
            save_on=false, save_everystep=false)
    end
    prob.set_sys(sam.integrator, sam.sys_struct)
    prob.set_set(sam.integrator, sam.set)
    OrdinaryDiffEqCore.reinit!(sam.integrator; reinit_dae=true)
    lin_vsm && linearize_vsm!(sam, sam.prob)
    update_sys_struct!(sam.prob, sam.integrator, sam.sys_struct)
    return sam.integrator, true
end

"""
    get_set_hash(set::Settings; fields)

Calculates a SHA1 hash for a subset of fields in the `Settings` object.
This is used to check if a cached model is still valid.
"""
function get_set_hash(set::Settings; 
        fields=[:segments, :model, :foil_file, :physical_model, :quasi_static, :winch_model]
    )
    h = zeros(UInt8, 1)
    for field in fields
        value = getfield(set, field)
        h = sha1(string((value, h)))
    end
    return h
end

"""
    get_sys_struct_hash(sys_struct::SystemStructure)

Calculates a SHA1 hash for the topology of a `SystemStructure`.
This is used to check if a cached model is still valid.
"""
function get_sys_struct_hash(sys_struct::SystemStructure)
    @unpack points, groups, segments, pulleys, tethers, winches, wings, transforms = sys_struct
    data_parts = []
    for point in points
        push!(data_parts, ("point", point.idx, point.wing_idx, Int(point.type)))
    end
    for segment in segments
        push!(data_parts, ("segment", segment.idx, segment.point_idxs))
    end
    for group in groups
        push!(data_parts, ("group", group.idx, group.point_idxs, Int(group.type)))
    end
    for pulley in pulleys
        push!(data_parts, ("pulley", pulley.idx, pulley.segment_idxs, Int(pulley.type)))
    end
    for tether in tethers
        push!(data_parts, ("tether", tether.idx, tether.segment_idxs))
    end
    for winch in winches
        model_type = winch.model isa TorqueControlledMachine
        push!(data_parts, ("winch", winch.idx, model_type, winch.tether_idxs))
    end
    for wing in wings
        push!(data_parts, ("wing", wing.idx, wing.group_idxs))
    end
    for transform in transforms
        push!(data_parts, ("transform", transform.idx, transform.wing_idx, transform.rot_point_idx, 
                 transform.base_point_idx, transform.base_transform_idx))
    end
    content = string(data_parts)
    return sha1(content)
end
