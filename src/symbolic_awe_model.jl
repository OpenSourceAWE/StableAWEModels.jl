# Copyright (c) 2024, 2025 Bart van de Lint and Uwe Fechner
# SPDX-License-Identifier: MIT

const LinType = @NamedTuple{A::Matrix{SimFloat}, B::Matrix{SimFloat}, C::Matrix{SimFloat}, D::Matrix{SimFloat}}
const GetSetNothing = Union{AbstractIndexer, Nothing}

"""
    @with_kw mutable struct SerializedModel{...}

A type-stable container for the compiled and serialized components of a `SymbolicAWEModel`.

This struct holds the products of the `ModelingToolkit.jl` compilation process,
such as the simplified system, the ODE problem, and various generated functions.
Caching these components allows for fast reloading of a pre-compiled model,
avoiding the time-consuming symbolic processing step on subsequent runs with the
same configuration.

$(TYPEDFIELDS)
"""
@with_kw mutable struct SerializedModel
    set_hash::Vector{UInt8}
    sys_struct_hash::Vector{UInt8}
    "Simplified system of the mtk model"
    sys::Union{ModelingToolkit.System, Nothing} = nothing
    "Unsimplified system of the mtk model"
    full_sys::Union{ModelingToolkit.System, Nothing} = nothing
    "Linearization function of the mtk model"
    lin_prob::Union{ModelingToolkit.LinearizationProblem, Nothing} = nothing
    lin_outputs::Union{Vector{Union{Symbolics.Arr, Symbolics.Num}}, Nothing} = nothing
    "ODE function of the mtk model"
    prob::Union{OrdinaryDiffEqCore.ODEProblem, Nothing} = nothing
    control_functions::Union{NamedTuple, Nothing} = nothing
    inputs::Union{Symbolics.Arr, Vector{Num}, Nothing} = nothing

    defaults::Vector{Pair} = Pair[]
    guesses::Vector{Pair} = Pair[]

    set_sys::GetSetNothing = nothing
    set_set_values::GetSetNothing = nothing
    set_set::GetSetNothing = nothing
    set_lin_set_values::GetSetNothing = nothing
    set_lin_sys::GetSetNothing = nothing
    set_lin_set::GetSetNothing = nothing
    
    get_set_values::GetSetNothing = nothing
    get_wing_state::GetSetNothing = nothing
    get_vsm_y::GetSetNothing = nothing
    get_segment_state::GetSetNothing = nothing
    get_winch_state::GetSetNothing = nothing
    get_tether_state::GetSetNothing = nothing
    get_struct_state::GetSetNothing = nothing
    get_point_state::GetSetNothing = nothing
    get_pulley_state::GetSetNothing = nothing
    get_group_state::GetSetNothing = nothing
    get_lin_x::GetSetNothing = nothing
    get_lin_dx::GetSetNothing = nothing
    get_lin_y::GetSetNothing = nothing

    lin_model::Union{LinType, Nothing} = nothing
    simple_lin_model::Union{LinType, Nothing} = nothing
end

"""
    mutable struct SymbolicAWEModel <: AbstractKiteModel

The main state container for a kite power system model, built using `ModelingToolkit.jl`.

This struct holds the complete state of the simulation, including the physical
structure (`SystemStructure`), the compiled model (`SerializedModel`), the atmospheric
model, and the ODE integrator.

Users typically interact with this model through high-level functions like
[`init!`](@ref) and [`next_step!`](@ref) rather than accessing its fields directly.

# Type Parameters
- `S`: Scalar type, typically `SimFloat`.
- `V`: Vector type, typically `KVec3`.
- `P`: Number of tether points in the system.

$(TYPEDFIELDS)
"""
@with_kw mutable struct SymbolicAWEModel <: AbstractKiteModel
    "Reference to the settings struct"
    set::Settings
    "Reference to the point mass system with points, segments, pulleys and tethers"
    sys_struct::SystemStructure
    "Container for the compiled and serialized model components"
    serialized_model::Any # Initially untyped, becomes a concrete SerializedModel after init
    "Reference to the atmospheric model as implemented in the package AtmosphericModels"
    am::AtmosphericModel = AtmosphericModel(set)
    "The ODE integrator for the full nonlinear model"
    integrator::Union{OrdinaryDiffEqCore.ODEIntegrator, Nothing} = nothing
    "Relative start time of the current time interval"
    t_0::SimFloat = 0.0
    "Number of next_step! calls"
    iter::Int64 = 0
    "Time spent in the VSM linearization step"
    t_vsm::SimFloat  = zero(SimFloat)
    "Time spent in the ODE integration step"
    t_step::SimFloat = zero(SimFloat)
    "Vector of tether length set-points"
    set_tether_len::Vector{SimFloat} = zeros(SimFloat, 3)
end

"""
    Base.getproperty(sam::SymbolicAWEModel, sym::Symbol)

Overloads `getproperty` to allow direct access to fields within the nested `serialized_model`.
This provides a convenient way to access compiled functions and other model
components without explicitly referencing `sam.serialized_model`.
"""
function Base.getproperty(sam::SymbolicAWEModel, sym::Symbol)
    if hasfield(SymbolicAWEModel, sym)
        getfield(sam, sym)
    else
        getproperty(getfield(sam, :serialized_model), sym)
    end
end

"""
    Base.setproperty!(sam::SymbolicAWEModel, sym::Symbol, val)

Overloads `setproperty!` to allow direct setting of fields within the nested `serialized_model`.
This allows you to change properties of the compiled model as if they were
fields of the `SymbolicAWEModel` itself.
"""
function Base.setproperty!(sam::SymbolicAWEModel, sym::Symbol, val)
    if hasfield(SymbolicAWEModel, sym)
        setfield!(sam, sym, val)
    else
        serialized_model = getfield(sam, :serialized_model)
        setproperty!(serialized_model, sym, val)
    end
end

"""
    SymbolicAWEModel(set::Settings, sys_struct::SystemStructure; kwargs...)

Constructs a `SymbolicAWEModel` from an existing `SystemStructure`.

This is the primary inner constructor. It takes a `SystemStructure` that defines the
physical layout of the kite system and prepares it for symbolic model generation.

# Arguments
- `set::Settings`: Configuration parameters.
- `sys_struct::SystemStructure`: The physical system definition.
- `kwargs...`: Further keyword arguments passed to the `SymbolicAWEModel` constructor.

# Returns
- `SymbolicAWEModel`: A model ready for symbolic equation generation via [`init!`](@ref).
"""
function SymbolicAWEModel(
    set::Settings, 
    sys_struct::SystemStructure;
    kwargs...
)
    set_hash = get_set_hash(set)
    sys_struct_hash = get_sys_struct_hash(sys_struct)
    # Initialize with an untyped, empty SerializedModel. It will be replaced by a 
    # fully typed one during init!
    serialized_model = SerializedModel(; set_hash, sys_struct_hash)
    return SymbolicAWEModel(; set, sys_struct, serialized_model, kwargs...)
end

"""
    SymbolicAWEModel(set::Settings; kwargs...)

Constructs a default `SymbolicAWEModel` with automatically generated components.

This convenience constructor automatically creates a complete AWE model:
- It first builds a `SystemStructure` based on the wing geometry and settings.
- Then, it assembles everything into a ready-to-use symbolic model.

# Arguments
- `set::Settings`: Configuration parameters.
- `kwargs...`: Further keyword arguments passed to the `SystemStructure` constructor.

# Returns
- `SymbolicAWEModel`: A model ready for symbolic equation generation via [`init!`](@ref).
"""
function SymbolicAWEModel(set::Settings; kwargs...)
    sys_struct = SystemStructure(set; kwargs...)
    return SymbolicAWEModel(set, sys_struct)
end

"""
    SymbolicAWEModel(set::Settings, name::String; kwargs...)

Constructs a `SymbolicAWEModel` for a specific named physical model.

This convenience constructor sets the `physical_model` field of the `Settings`
struct and then proceeds to create the model.
"""
function SymbolicAWEModel(set::Settings, name::String; kwargs...)
    set.physical_model = name
    sys_struct = SystemStructure(set; kwargs...)
    return SymbolicAWEModel(set, sys_struct)
end

"""
    update_sys_state!(ss::SysState, s::SymbolicAWEModel, zoom=1.0)

Updates a `SysState` object with the current state values from the `SymbolicAWEModel`.

This function takes the raw data from the model's internal integrator and populates
the fields of the user-friendly `SysState` struct, converting units (e.g., radians
to degrees) and calculating derived values like AoA and roll/pitch/yaw angles.

# Arguments
- `ss::SysState`: The state struct to be updated.
- `s::SymbolicAWEModel`: The source model.
- `zoom::SimFloat=1.0`: A scaling factor for the position coordinates.
"""
function update_sys_state!(ss::SysState, sam::SymbolicAWEModel, zoom=1.0)
    ss.time = isnothing(sam.integrator) ? 0.0 : sam.integrator.t # Use integrator time
    @unpack points, groups, segments, pulleys, winches, wings = sam.sys_struct

    # Get the state vectors from the integrator
    if length(winches) > 0
        for winch in winches
            ss.l_tether[winch.idx] = winch.tether_len
            ss.v_reelout[winch.idx] = winch.tether_vel
            ss.force[winch.idx] = norm(winch.force)
            ss.set_torque[winch.idx] = winch.set_value
        end
    end
    if length(groups) > 0
        for group in groups
            ss.twist_angles[group.idx] = group.twist
        end
        ss.depower = rad2deg(mean(ss.twist_angles)) # Average twist for depower
        ss.steering = rad2deg(ss.twist_angles[length(groups)] - ss.twist_angles[1])
    end
    if length(wings) > 0
        wing = wings[1]
        ss.acc = norm(wing.acc_w) # Use the norm of the wing's acceleration vector
        ss.orient .= wing.Q_b_w
        ss.turn_rates .= wing.turn_rate
        ss.elevation = wing.elevation
        ss.azimuth = wing.azimuth
        ss.heading = wing.heading
        ss.course = wing.course
        # Apparent Wind and Aerodynamics
        ss.v_app = norm(wing.va_b)
        ss.v_wind_kite .= wing.v_wind
        # Calculate AoA and Side Slip from apparent wind in body frame
        if ss.v_app > 1e-6 # Avoid division by zero
            ss.AoA = atan(wing.va_b[3], wing.va_b[1])
            ss.side_slip = asin(wing.va_b[2] / norm(wing.va_b))
        else
            ss.AoA = 0.0
            ss.side_slip = 0.0
        end
        ss.aero_force_b .= wing.aero_force_b
        ss.aero_moment_b .= wing.aero_moment_b
        ss.vel_kite .= wing.vel_w
        # Calculate Roll, Pitch, Yaw from Quaternion
        q = wing.Q_b_w
        sinr_cosp = 2 * (q[1] * q[2] + q[3] * q[4])
        cosr_cosp = 1 - 2 * (q[2] * q[2] + q[3] * q[3])
        ss.roll = atan(sinr_cosp, cosr_cosp)
        sinp = 2 * (q[1] * q[3] - q[4] * q[2])
        ss.pitch = abs(sinp) >= 1 ? copysign(pi / 2, sinp) : asin(sinp)
        siny_cosp = 2 * (q[1] * q[4] + q[2] * q[3])
        cosy_cosp = 1 - 2 * (q[3] * q[3] + q[4] * q[4])
        ss.yaw = atan(siny_cosp, cosy_cosp)
    end
    for point in points
        ss.X[point.idx] = point.pos_w[1] * zoom
        ss.Y[point.idx] = point.pos_w[2] * zoom
        ss.Z[point.idx] = point.pos_w[3] * zoom
    end
    ss.v_wind_gnd .= sam.sys_struct.wind_vec_gnd
    nothing
end

"""
    SysState(s::SymbolicAWEModel, zoom=1.0)

Constructs a `SysState` object from a `SymbolicAWEModel`.

This is a convenience constructor that creates a new `SysState` object and populates it
with the current state of the provided model.

# Arguments
- `s::SymbolicAWEModel`: The source model.
- `zoom::SimFloat=1.0`: A scaling factor for the position coordinates.

# Returns
- `SysState`: A new state struct representing the current model state.
"""
function SysState(s::SymbolicAWEModel, zoom=1.0)
    ss = SysState{length(s.sys_struct.points)}()
    update_sys_state!(ss, s, zoom)
    ss
end

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
            # If hashes do not match after deserialization, reset serialized model
            if set_hash != serialized_model.set_hash || sys_struct_hash != serialized_model.sys_struct_hash
                sam.serialized_model = SerializedModel(; set_hash, sys_struct_hash)
                return false
            end
            # Success: assign deserialized model
            sam.serialized_model = serialized_model
            return true
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
    maybe_create_prob!(s, lin_outputs; create_prob=true, prn=true)

Conditionally creates and stores the `ODEProblem` and generated getter/setter functions.

This function simplifies the full symbolic system and creates the `ODEProblem` needed
for simulation. It also generates and stores optimized functions for accessing
and updating model state, which is a key performance optimization.
"""
function maybe_create_prob!(sam, lin_outputs; create_prob=true, prn=true)
    if create_prob && isnothing(sam.prob)
        prn && @info "Simplifying the System..."
        time = @elapsed @suppress_err begin
            sam.sys = mtkcompile(sam.full_sys; inputs = sam.inputs,
                                 additional_passes = [ModelingToolkit.IfLifting])
        end
        prn && @info "Simplified the System in $time seconds."

        prn && @info "Creating the ODEProblem..."
        dt = SimFloat(1/sam.set.sample_freq)
        time = @elapsed sam.prob = ODEProblem(sam.sys, sam.defaults, (0.0, dt); sam.guesses)
        prn && @info "Created the ODEProblem in $time seconds."

        prn && @info "Creating getter functions..."
        generate_getters!(sam, sam.sys, lin_outputs)
        prn && @info "Created the getter functions in $time seconds."
        return true
    end
    return false
end

"""
    maybe_create_lin_prob!(s, lin_outputs; create_lin_prob=true, prn=true)

Conditionally creates and stores a `LinearizationProblem` for the model.

This function is called if a linearization problem is needed and not already present.
It sets up the symbolic functions necessary to linearize the system around an operating point.
"""
function maybe_create_lin_prob!(sam, lin_outputs; create_lin_prob=true, prn=true)
    if create_lin_prob
        if isnothing(sam.lin_prob) ||
           isnothing(sam.lin_outputs) ||
           length(lin_outputs) != length(sam.serialized_model.lin_outputs) ||
           !all(string.(lin_outputs) .== string.(sam.serialized_model.lin_outputs)) 

            prn && @info "Creating the LinearizationProblem..."
            time = @elapsed @suppress_err begin
                lin_fun, lin_sys = linearization_function(sam.full_sys, [sam.inputs...], lin_outputs;
                                                      op=sam.defaults, guesses=sam.guesses)
                sam.lin_prob = LinearizationProblem(lin_fun, 0.0)
                generate_lin_getters!(sam, lin_sys)
            end
            prn && @info "Created the LinearizationProblem in $time seconds."
            return true
        end
    end
    return false
end

"""
    maybe_create_control_functions!(s, lin_outputs; create_control_func=false, prn=true)

Conditionally generates control functions if requested and not yet present.

This function creates the functions required for advanced control, such as state and
output mapping from inputs, by using `ModelingToolkit.jl`'s control function
generation capabilities.
"""
function maybe_create_control_functions!(sam, lin_outputs; create_control_func=false, prn=true)
    if create_control_func
        if isnothing(sam.control_functions) ||
           length(lin_outputs) != length(sam.serialized_model.lin_outputs) ||
           !all(string.(lin_outputs) .== string.(sam.serialized_model.lin_outputs)) 

            function generate_f_h(model, inputs, outputs)
                (f_ip, f_oop), dvs, psym, io_sys = @suppress_err ModelingToolkit.generate_control_function(
                    model, inputs; simplify=false
                )
                nu, nx, ny = length(inputs), length(dvs), length(outputs)
                (h_oop, h_ip) = ModelingToolkit.build_explicit_observed_function(
                    io_sys, outputs; inputs, return_inplace = true
                )
                return (f_oop=f_oop, f_ip=f_ip, h_oop=h_oop, h_ip=h_ip, nu=nu, nx=nx, ny=ny,
                          dvs=dvs, psym=psym, io_sys=io_sys)
            end
            prn && @info "Creating the control functions..."
            inputs = [sam.inputs...]
            time = @elapsed sam.serialized_model.control_functions = generate_f_h(sam.full_sys, inputs, lin_outputs)
            prn && @info "Created the control functions in $time seconds."
            return true
        end
    end
    return false
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

function generate_lin_getters!(sam, sys)
    sam.set_lin_set_values = nothing
    if length(sam.sys_struct.winches) > 0
        sam.set_lin_set_values = setp(sys, sys.set_values)
    end
    sam.set_lin_sys=setp(sys, sys.psys)
    sam.set_lin_set=setp(sys, sys.pset)
    nothing
end

"""
    generate_getters(sys, sys_struct, lin_prob, lin_y_vec) -> NamedTuple

Generates and compiles optimized getter and setter functions for the model.

This internal function uses the symbolic system definition from `ModelingToolkit.jl`
to create fast, non-allocating functions for accessing and modifying the system's
state and parameters directly within the ODE integrator's data structures. This is
a key optimization that avoids symbolic lookups during the simulation loop.

# Returns
- `NamedTuple`: A named tuple containing all the generated functions and a simplified
                linear model representation.
"""
function generate_getters!(sam, sys, lin_y_vec)
    c = collect
    @unpack wings, groups, pulleys, winches, tethers, segments = sam.sys_struct

    # Initialize all potential functions to nothing
    sam.get_lin_x, sam.get_lin_dx, sam.get_lin_y = nothing, nothing, nothing
    sam.get_wing_state, sam.get_vsm_y = nothing, nothing
    sam.get_segment_state, sam.get_group_state, sam.get_pulley_state = nothing, nothing, nothing
    sam.get_winch_state, sam.get_tether_state = nothing, nothing
    sam.set_set_values, sam.get_set_values = nothing, nothing
    sam.simple_lin_model = nothing

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
        sam.get_lin_x = getu(sys, lin_x_vec)
        sam.get_lin_dx = getu(sys, lin_dx_vec)
        sam.get_lin_y = getu(sys, lin_y_vec)

        nx = length(lin_x_vec)
        ny = length(lin_y_vec)
        nu = length(winches)
        sam.simple_lin_model = (
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
        sam.get_wing_state = getu(sys, wing_state_vars)
        sam.get_vsm_y = getu(sys, sys.y)
    end

    if length(segments) > 0
        sam.get_segment_state = getu(sys, c.([sys.spring_force, sys.len]))
    end

    if length(groups) > 0
        group_state_vars = c.([
            sys.twist_angle, sys.twist_ω, sys.group_tether_force,
            sys.group_tether_moment, sys.group_aero_moment,
        ])
        sam.get_group_state = getu(sys, group_state_vars)
    end
    
    if length(pulleys) > 0
        sam.get_pulley_state = getu(sys, c.([sys.pulley_len, sys.pulley_vel]))
    end

    if length(winches) > 0
        winch_state_vars = c.([
             sys.tether_len, sys.tether_vel, sys.set_values, sys.winch_force_vec,
        ])
        sam.get_winch_state = getu(sys, winch_state_vars)
        sam.set_set_values = setp(sys, sys.set_values)
        sam.get_set_values = getp(sys, sys.set_values)
    end

    if length(tethers) > 0
        sam.get_tether_state = getu(sys, c(sys.stretched_len))
    end

    sam.set_sys = setp(sys, sys.psys)
    sam.set_set = setp(sys, sys.pset)
    sam.get_struct_state = getu(sys, sys.wind_vec_gnd)
    sam.get_point_state = getu(sys, c.([sys.pos, sys.vel, sys.point_force]))
    nothing
end

"""
    next_step!(s::SymbolicAWEModel, integrator::ODEIntegrator; set_values, dt, vsm_interval)

Take a simulation step, using the provided integrator.

This is a convenience method that calls the main `next_step!` function.
"""
function next_step!(s::SymbolicAWEModel, integrator::OrdinaryDiffEqCore.ODEIntegrator; set_values=nothing, dt=1/s.set.sample_freq, vsm_interval=1)
    !(s.integrator === integrator) && error("The ODEIntegrator doesn't belong to the SymbolicAWEModel")
    next_step!(s; set_values=set_values, dt=dt, vsm_interval=vsm_interval)
end

"""
    next_step!(s::SymbolicAWEModel; set_values, dt, vsm_interval)

Take a simulation step forward in time.

This function advances the simulation by one time step, optionally updating control
inputs and re-linearizing the VSM model. It then updates the `SystemStructure`
with the new state from the ODE integrator.

# Arguments
- `s::SymbolicAWEModel`: The kite power system state object.

# Keyword Arguments
- `set_values=nothing`: New values for the control inputs. If `nothing`, the current values are used.
- `dt=1/s.set.sample_freq`: Time step size [s].
- `vsm_interval=1`: The interval (in steps) to re-linearize the VSM model. If 0, it is not re-linearized.
"""
function next_step!(sam::SymbolicAWEModel; set_values=nothing, dt=1/sam.set.sample_freq, vsm_interval=1)
    if (!isnothing(set_values)) 
        sam.set_set_values(sam.integrator, set_values)
    end
    if vsm_interval != 0 && sam.iter % vsm_interval == 0
        sam.t_vsm = @elapsed linearize_vsm!(sam)
    end
    
    sam.t_0 = sam.integrator.t
    sam.t_step = @elapsed OrdinaryDiffEqCore.step!(sam.integrator, dt, true)
    if !successful_retcode(sam.integrator.sol)
        @warn "Return code for solution: $(sam.integrator.sol.retcode)"
    end
    @assert successful_retcode(sam.integrator.sol)
    sam.iter += 1
    update_sys_struct!(sam, sam.sys_struct)
    return nothing
end

"""
    update_sys_struct!(s::SymbolicAWEModel, sys_struct::SystemStructure, integ=s.integrator)

Updates the high-level `SystemStructure` from the low-level integrator state vector.

This function reads the raw state vector from the ODE integrator and uses the generated
getter functions to populate the human-readable fields in the `SystemStructure`. This
synchronization step is crucial for making the simulation results accessible.
"""
function update_sys_struct!(sam::SymbolicAWEModel, sys_struct::SystemStructure, integ=sam.integrator)
    @unpack points, groups, segments, pulleys, winches, tethers, wings = sys_struct
    pos, vel, force = sam.get_point_state(integ)
    for point in points
        point.pos_w .= pos[:, point.idx]
        point.vel_w .= vel[:, point.idx]
        point.force .= force[:, point.idx]
    end
    if length(pulleys) > 0
        len, vel = sam.get_pulley_state(integ)
        for pulley in pulleys
            pulley.len = len[pulley.idx]
            pulley.vel = vel[pulley.idx]
        end
    end
    if length(segments) > 0
        spring_force, len = sam.get_segment_state(integ)
        for segment in segments
            segment.force = spring_force[segment.idx]
            segment.len = len[segment.idx]
        end
    end
    if length(groups) > 0
        twist, twist_ω, tether_force, tether_moment, aero_moment = sam.get_group_state(integ)
        for group in groups
            group.twist = twist[group.idx]
            group.twist_ω = twist_ω[group.idx]
            group.tether_force = tether_force[group.idx]
            group.tether_moment = tether_moment[group.idx]
            group.aero_moment = aero_moment[group.idx]
        end
    end
    if length(winches) > 0
        tether_len, tether_vel, set_value, winch_force_vec = sam.get_winch_state(integ)
        for winch in winches
            winch.tether_len = tether_len[winch.idx]
            winch.tether_vel = tether_vel[winch.idx]
            winch.set_value = set_value[winch.idx]
            winch.force .= winch_force_vec[:, winch.idx]
        end
    end
    if length(tethers) > 0
        stretched_len = sam.get_tether_state(integ)
        for tether in tethers
            tether.stretched_len = stretched_len[tether.idx]
        end
    end
    if length(wings) > 0
        wing_state = sam.get_wing_state(integ)
        Q_b_w, ω_b, pos_w, vel_w, acc_w, va_b, v_wind, 
            aero_force_b, aero_moment_b, elevation, elevation_vel,
            elevation_acc, azimuth, azimuth_vel, azimuth_acc,
            heading, turn_rate, turn_acc, course, aoa = wing_state
        for wing in wings
            wing.Q_b_w .= Q_b_w[wing.idx, :]
            wing.ω_b .= ω_b[wing.idx, :]
            wing.pos_w .= pos_w[wing.idx, :]
            wing.vel_w .= vel_w[wing.idx, :]
            wing.acc_w .= acc_w[wing.idx, :]
            wing.va_b .= va_b[wing.idx, :]
            wing.v_wind .= v_wind[wing.idx, :]
            wing.aero_force_b .= aero_force_b[wing.idx, :]
            wing.aero_moment_b .= aero_moment_b[wing.idx, :]
            wing.elevation = elevation[wing.idx]
            wing.elevation_vel = elevation_vel[wing.idx]
            wing.elevation_acc = elevation_acc[wing.idx]
            wing.azimuth = azimuth[wing.idx]
            wing.azimuth_vel = azimuth_vel[wing.idx]
            wing.azimuth_acc = azimuth_acc[wing.idx]
            wing.heading = heading[wing.idx]
            wing.turn_rate .= turn_rate[wing.idx, :]
            wing.turn_acc .= turn_acc[wing.idx, :]
            wing.course = course[wing.idx]
            wing.aoa = aoa[wing.idx]
        end
    end
    sam.sys_struct.wind_vec_gnd .= sam.get_struct_state(integ)
    return nothing
end

"""
    get_model_name(set::Settings; precompile=false)

Constructs a unique filename for the serialized model based on its configuration.
The filename includes the Julia version, physical model, dynamics type, and number of
segments to ensure that the correct cached model is loaded.
"""
function get_model_name(set::Settings; precompile=false)
    suffix = ""
    ver = "$(VERSION.major).$(VERSION.minor)"
    if precompile
        suffix = ".default"
    end
    dynamics_type = ifelse(set.quasi_static, "static", "dynamic")
    return "model_$(ver)_$(set.physical_model)_$(dynamics_type)_$(set.segments)_seg.bin$suffix"
end

"""
    calc_aoa(s::SymbolicAWEModel)

Calculates the mean angle of attack [rad] over the wingspan from the VSM solver.
"""
function calc_aoa(sam::SymbolicAWEModel)
    alpha_array = sam.sys_struct.wings[1].vsm_solver.sol.alpha_array
    middle = length(alpha_array) ÷ 2
    return iseven(length(alpha_array)) ? (0.5 * (alpha_array[middle] + alpha_array[middle+1])) : alpha_array[middle+1]
end

"""
    unstretched_length(s::SymbolicAWEModel)

Returns the unstretched tether length [m] for each winch.
"""
unstretched_length(sam::SymbolicAWEModel) = [winch.tether_len for winch in sam.sys_struct.winches]

"""
    tether_length(s::SymbolicAWEModel)

Returns the current stretched tether length [m] for each tether.
"""
tether_length(sam::SymbolicAWEModel) = [tether.stretched_len for tether in sam.sys_struct.tethers]

"""
    calc_height(s::SymbolicAWEModel)

Returns the height (z-position) [m] of the wing.
"""
calc_height(sam::SymbolicAWEModel) = sam.sys_struct.wings[1].pos_w[3]

"""
    winch_force(s::SymbolicAWEModel)

Returns the winch force [N] for each winch.
"""
winch_force(sam::SymbolicAWEModel) = [norm(winch.force) for winch in sam.sys_struct.winches]

"""
    spring_forces(s::SymbolicAWEModel)

Returns the spring force [N] for each tether segment.
"""
spring_forces(sam::SymbolicAWEModel) = [segment.force for segment in sam.sys_struct.segments]

"""
    pos(s::SymbolicAWEModel)

Returns a vector of the position vectors [m] for each point in the system.
"""
pos(sam::SymbolicAWEModel) = [point.pos_w for point in sam.sys_struct.points]

"""
    min_chord_len(s::SymbolicAWEModel)

Calculates the minimum chord length of the wing at the tip.
"""
function min_chord_len(sam::SymbolicAWEModel)
    min_len = Inf
    for wing in sam.sys_struct.wings
        vsm_wing = wing.vsm_wing
        le_pos = [vsm_wing.le_interp[i](vsm_wing.gamma_tip) for i in 1:3]
        te_pos = [vsm_wing.te_interp[i](vsm_wing.gamma_tip) for i in 1:3]
        min_len = min(norm(le_pos - te_pos), min_len)
    end
    return min_len
end

"""
    set_depower_steering!(s::SymbolicAWEModel, depower, steering)

Sets the kite's depower and steering by adjusting the tether length set-points.
"""
function set_depower_steering!(sam::SymbolicAWEModel, depower, steering)
    len = sam.set_tether_len
    len .= [winch.tether_len for winch in sam.sys_struct.winches]
    depower *= min_chord_len(sam)
    steering *= min_chord_len(sam)
    len[2] = 0.5 * (2*depower + 2*len[1] + steering)
    len[3] = 0.5 * (2*depower + 2*len[1] - steering)
    return nothing
end

"""
    set_v_wind_ground!(s::SymbolicAWEModel, v_wind_gnd=s.set.v_wind, upwind_dir=-π/2)

Sets the ground wind speed [m/s] and upwind direction [rad] in the model.
"""
function set_v_wind_ground!(sam::SymbolicAWEModel, v_wind_gnd=sam.set.v_wind, upwind_dir=-pi/2)
    sam.set.v_wind = v_wind_gnd
    sam.set.upwind_dir = rad2deg(upwind_dir)
    sam.set_set(sam.integrator, sam.set)
    return nothing
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
