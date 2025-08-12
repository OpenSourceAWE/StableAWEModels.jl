# Copyright (c) 2025 Bart van de Lint and Uwe Fechner
# SPDX-License-Identifier: MPL-2.0

"""
    load_serialized_model!(s, model_path; remake=false, reload=false)

Attempts to load a serialized model from a file and validate it.

The function checks if the model exists at `model_path` and if its hashes
match the current settings and system structure. If a match is found and
`remake` is false, it loads the model.

# Arguments
- `s`: The `SymbolicAWEModel` to load into.
- `model_path`: The path to the serialized model file.
- `remake::Bool=false`: If true, forces the model to be re-created even if a file exists.
- `reload::Bool=false`: If true, forces a reload from disk even if a model is already loaded.

# Returns
- `Bool`: `true` if the model was successfully loaded, `false` otherwise.
"""
function load_serialized_model!(sam, model_path; remake=false, reload=false)
    set_hash = get_set_hash(sam.set)
    sys_struct_hash = get_sys_struct_hash(sam.sys_struct)
    
    # Check if we have a serialized model file and don't want to remake it
    if ispath(model_path) && !remake
        # If the current hashes do not match the stored ones, reset the serialized model
        if set_hash != sam.set_hash || sys_struct_hash != sam.sys_struct_hash
            sam.serialized_model = SerializedModel(; set_hash, sys_struct_hash)
            return false
        end

        # If the full system is already loaded and reload not requested, keep it as is
        if !isnothing(sam.full_sys) && !reload
            return true
        end

        # Attempt to deserialize the stored serialized model
        try
            serialized_model = deserialize(model_path)
            if set_hash == serialized_model.set_hash &&
               sys_struct_hash == serialized_model.sys_struct_hash
                # Success: assign deserialized model
                sam.serialized_model = serialized_model
                return true
            end
            # Fall through to recreate serialized model below
        catch e
            @warn "Failure to deserialize $model_path: $(typeof(e))"
            # Fall through to recreate serialized model below
        end
    end
    # If no file or remake requested, initialize a new SerializedModel
    sam.serialized_model = SerializedModel(; set_hash, sys_struct_hash)
    return false
end

"""
    maybe_create_prob!(sam, lin_outputs; create_prob=true, prn=true)

Conditionally creates and stores the `ODEProblem` and its attributes. (Updated)
"""
function maybe_create_prob!(sam, lin_outputs; create_prob=true, prn=true)
    if create_prob && isnothing(sam.prob)
        prn && @info "Simplifying the System..."
        local sys
        time = @elapsed @suppress_err sys = mtkcompile(sam.full_sys; inputs=sam.inputs)
        prn && @info "Simplified the System in $time seconds."

        prn && @info "Creating the ODEProblem..."
        dt = SimFloat(1/sam.set.sample_freq)
        local prob
        time = @elapsed prob = ODEProblem(sys, sam.defaults, (0.0, dt); u0_prior=sam.guesses)
        prn && @info "Created the ODEProblem in $time seconds."

        prn && @info "Creating getter functions..."
        local getters
        time = @elapsed getters = generate_getters(sam.sys_struct, sys, lin_outputs)
        prn && @info "Created the getter functions in $time seconds."

        sam.prob = ProbWithAttributes(;
            prob=prob,
            getters...
        )
        return true
    end
    return false
end

"""
    maybe_create_lin_prob!(sam, lin_outputs; create_lin_prob=true, prn=true)

Conditionally creates and stores a `LinearizationProblem`. (Fully Refactored)
"""
function maybe_create_lin_prob!(sam, lin_outputs; create_lin_prob=true, prn=true)
    if create_lin_prob && (isnothing(sam.lin_prob) ||
           length(lin_outputs) != length(sam.lin_prob.lin_outputs) ||
           !all(string.(lin_outputs) .== string.(sam.lin_prob.lin_outputs)))

        prn && @info "Creating the LinearizationProblem..."
        local lin_fun, lin_sys, lin_prob_instance
        time = @elapsed @suppress_err begin
            lin_fun, lin_sys = linearization_function(sam.full_sys, [sam.inputs...], lin_outputs;
                                                      op=sam.defaults, guesses=sam.guesses)
            lin_prob_instance = LinearizationProblem(lin_fun, 0.0)
            
            # Call the new pure generator function
            getters = generate_lin_getters(lin_sys)

            # Construct the attribute struct
            sam.lin_prob = LinProbWithAttributes(;
                lin_prob=lin_prob_instance,
                lin_outputs=lin_outputs,
                getters...
            )
        end
        prn && @info "Created the LinearizationProblem in $time seconds."
        return true
    end
    return false
end

"""
    maybe_create_control_functions!(sam, lin_outputs; create_control_func=false, prn=true)

Conditionally generates control functions. (Fully Refactored)
"""
function maybe_create_control_functions!(sam, lin_outputs; create_control_func=false, prn=true)
    if create_control_func && (isnothing(sam.control_funcs) ||
           length(lin_outputs) != length(sam.lin_prob.lin_outputs) ||
           !all(string.(lin_outputs) .== string.(sam.lin_prob.lin_outputs)))

        prn && @info "Creating the control functions..."
        inputs_vec = [sam.inputs...]
        
        time = @elapsed result = generate_control_funcs(sam.full_sys, inputs_vec, lin_outputs)

        sam.control_funcs = ControlFuncWithAttributes(;
            result...
        )
        prn && @info "Created the control functions in $time seconds."
        return true
    end
    return false
end

"""
    generate_lin_getters(sys) -> NamedTuple

Creates setter functions for the linearization problem.
"""
function generate_lin_getters(sys)
    set_lin_set_values = nothing
    # Assuming winches are present if set_values is
    if hasproperty(sys, :set_values)
        set_lin_set_values = setp(sys, sys.set_values)
    end
    set_lin_sys = setp(sys, sys.psys)
    set_lin_set = setp(sys, sys.pset)

    return (
        set_lin_set_values=set_lin_set_values,
        set_lin_sys=set_lin_sys,
        set_lin_set=set_lin_set
    )
end

"""
    generate_control_funcs(model, inputs, outputs) -> NamedTuple

Generates the full suite of control-related functions (f, h, etc.).
"""
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

"""
    generate_getters(sys_struct, sys, lin_y_vec) -> NamedTuple

Creates getter and setter functions based on the system's symbolic components.

This is a pure function that takes symbolic system definitions and returns a
NamedTuple containing all the generated functions and the simple linear model.
It does not modify any state.

# Arguments
- `sys_struct::SystemStructure`: The high-level definition of the system structure.
- `sys::ModelingToolkit.System`: The simplified MTK system.
- `lin_y_vec`: A vector of symbolic variables for the linear model's output `y`.

# Returns
- A `NamedTuple` with keys for each getter/setter and the `simple_lin_model`.
"""
function generate_getters(sys_struct, sys, lin_y_vec)
    c = collect
    @unpack wings, groups, pulleys, winches, tethers, segments = sys_struct

    # Initialize local variables for all potential outputs
    get_lin_x, get_lin_dx, get_lin_y = nothing, nothing, nothing
    get_wing_state, get_vsm_y = nothing, nothing
    get_segment_state, get_group_state, get_pulley_state = nothing, nothing, nothing
    get_winch_state, get_tether_state = nothing, nothing
    set_set_values, get_set_values = nothing, nothing
    simple_lin_model = nothing

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

        nx = length(lin_x_vec)
        ny = length(lin_y_vec)
        nu = length(winches)
        simple_lin_model = (
            A = zeros(nx, nx),
            B = zeros(nx, nu),
            C = zeros(ny, nx),
            D = zeros(ny, nu)
        )
    end

    if length(wings) > 0
        wing_state_vars = c.([
            sys.Q_b_w, sys.ω_b, sys.wing_pos, sys.wing_vel, sys.wing_acc,
            sys.va_wing_b, sys.wind_vel_wing, sys.aero_force_b, sys.aero_moment_b,
            sys.elevation, sys.elevation_vel, sys.elevation_acc,
            sys.azimuth, sys.azimuth_vel, sys.azimuth_acc,
            sys.heading, sys.turn_rate, sys.turn_acc,
            sys.course, sys.angle_of_attack,
        ])
        get_wing_state = getu(sys, wing_state_vars)
        get_vsm_y = getu(sys, sys.y)
    end

    if length(segments) > 0
        get_segment_state = getu(sys, c.([sys.spring_force, sys.len]))
    end

    if length(groups) > 0
        group_state_vars = c.([
            sys.twist_angle, sys.twist_ω, sys.group_tether_force,
            sys.group_tether_moment, sys.group_aero_moment,
        ])
        get_group_state = getu(sys, group_state_vars)
    end
    
    if length(pulleys) > 0
        get_pulley_state = getu(sys, c.([sys.pulley_len, sys.pulley_vel]))
    end

    if length(winches) > 0
        winch_state_vars = c.([
             sys.tether_len, sys.tether_vel, sys.set_values, sys.winch_force_vec,
        ])
        get_winch_state = getu(sys, winch_state_vars)
        set_set_values = setp(sys, sys.set_values)
        get_set_values = getp(sys, sys.set_values)
    end

    if length(tethers) > 0
        get_tether_state = getu(sys, c(sys.stretched_len))
    end

    set_sys = setp(sys, sys.psys)
    set_set = setp(sys, sys.pset)
    get_struct_state = getu(sys, sys.wind_vec_gnd)
    get_point_state = getu(sys, c.([sys.pos, sys.vel, sys.point_force]))
    
    # Package all results into a single NamedTuple and return it
    return (
        get_lin_x=get_lin_x, get_lin_dx=get_lin_dx, get_lin_y=get_lin_y,
        get_wing_state=get_wing_state, get_vsm_y=get_vsm_y,
        get_segment_state=get_segment_state, get_group_state=get_group_state,
        get_pulley_state=get_pulley_state, get_winch_state=get_winch_state,
        get_tether_state=get_tether_state, set_set_values=set_set_values,
        get_set_values=get_set_values, simple_lin_model=simple_lin_model,
        set_sys=set_sys, set_set=set_set, get_struct_state=get_struct_state,
        get_point_state=get_point_state
    )
end

"""
    init!(s::SymbolicAWEModel; kwargs...)

Orchestrates the entire model loading, building, and integrator creation process.

This is the main initialization function. It first checks for a cached model file,
creates the symbolic model if necessary, and then initializes the ODE integrator.
It also handles the creation of linearization and control functions if requested.

# Keyword Arguments
- `solver`: The ODE solver to use (e.g., `FBDF`). Defaults to the solver specified in `s.set`.
- `adaptive::Bool=true`: Whether to use an adaptive time-stepping algorithm.
- `prn::Bool=true`: Whether to print progress information to the console.
- `remake::Bool=false`: Forces the model to be re-created from scratch, ignoring cached files.
- `reload::Bool=false`: Forces a reload from a cached file, even if the model is already in memory.
- `lin_outputs::Vector{Num}=nothing`: Specifies the outputs for the linearized model.
- `create_prob::Bool=true`: If `true`, creates the `ODEProblem`.
- `create_lin_prob::Bool=true`: If `true`, creates the `LinearizationProblem`.
- `create_control_func::Bool=false`: If `true`, generates control functions.
- `lin_vsm::Bool=true`: If `true`, linearizes the VSM model during reinitialization.

# Returns
- `ODEIntegrator`: The initialized ODE integrator.
"""
function init!(sam::SymbolicAWEModel;
    solver=nothing, adaptive=true, prn=true,
    remake=false, reload=false,
    lin_outputs=nothing,
    create_prob::Bool=true,
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
        changed = false # wether or not any changes were made to the serialized model
        if !loaded
            sam.inputs = create_sys!(sam, sam.sys_struct; prn)
            changed = true
        end
        changed = changed | maybe_create_prob!(sam, lin_outputs; create_prob, prn)
        changed = changed | maybe_create_lin_prob!(sam, lin_outputs; create_lin_prob, prn)
        changed = changed | maybe_create_control_functions!(sam, lin_outputs;
                                                           create_control_func, prn)
        if changed
            prn && @info "Serializing model."
            sam.lin_outputs = lin_outputs
            serialize(model_path, sam.serialized_model)
        end

        reinit!(sam.sys_struct, sam.set)
        create_prob && reinit!(sam, sam.prob, solver; adaptive, reload, lin_vsm)
        create_lin_prob && reinit!(sam, sam.lin_prob)
        # create_control_func && reinit!(s, sam.control_functions)
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
function reinit!(sam::SymbolicAWEModel, lin_prob::ModelingToolkit.LinearizationProblem)
    sam.set_lin_sys(lin_prob, sam.sys_struct)
    sam.set_lin_set(lin_prob, sam.set)
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
    prob::ODEProblem,
    solver;
    adaptive=true,
    reload=true, 
    lin_vsm=true
)
    if isnothing(sam.integrator) || !successful_retcode(sam.integrator.sol) || reload
        dt = SimFloat(1/sam.set.sample_freq)
        sam.integrator = OrdinaryDiffEqCore.init(prob, solver; 
            adaptive, dt, tspan=(0.0, dt), abstol=sam.set.abs_tol, reltol=sam.set.rel_tol, 
            save_on=false, save_everystep=false)
    end
    sam.set_sys(sam.integrator, sam.sys_struct)
    sam.set_set(sam.integrator, sam.set)
    OrdinaryDiffEqCore.reinit!(sam.integrator; reinit_dae=true)
    lin_vsm && linearize_vsm!(sam)
    update_sys_struct!(sam, sam.sys_struct)
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
