# Copyright (c) 2024, 2025 Bart van de Lint and Uwe Fechner
# SPDX-License-Identifier: MIT

const LinType = @NamedTuple{A::Matrix{SimFloat}, B::Matrix{SimFloat}, C::Matrix{SimFloat}, D::Matrix{SimFloat}}

"""
    @with_kw mutable struct SerializedModel

A container for the compiled and serialized components of a `SymbolicAWEModel`.

This struct holds the products of the `ModelingToolkit.jl` compilation process,
such as the simplified system, the ODE problem, and various generated functions.
Storing these allows for fast reloading of a pre-compiled model, avoiding the
time-consuming symbolic processing step on subsequent runs with the same configuration.

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

    defaults::Vector{Pair} = Pair[]
    guesses::Vector{Pair} = Pair[]

    set_psys::Union{Function, Nothing}           = nothing
    set_set_values::Union{Function, Nothing}   = nothing
    set_set::Union{Function, Nothing}           = nothing
    set_lin_set_values::Union{Function, Nothing}= nothing
    
    get_set_values::Union{Function, Nothing}   = nothing
    get_wing_state::Union{Function, Nothing}    = nothing
    get_vsm_y::Union{Function, Nothing}         = nothing
    get_segment_state::Union{Function, Nothing} = nothing
    get_winch_state::Union{Function, Nothing}   = nothing
    get_tether_state::Union{Function, Nothing}  = nothing
    get_struct_state::Union{Function, Nothing}  = nothing
    get_point_state::Union{Function, Nothing}   = nothing
    get_pulley_state::Union{Function, Nothing}  = nothing
    get_group_state::Union{Function, Nothing}   = nothing
    get_spring_force::Union{Function, Nothing}  = nothing
    get_lin_x::Union{Function, Nothing}         = nothing
    get_lin_dx::Union{Function, Nothing}        = nothing
    get_lin_y::Union{Function, Nothing}         = nothing

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
    serialized_model::SerializedModel
    "Reference to the atmospheric model as implemented in the package AtmosphericModels"
    am::AtmosphericModel = AtmosphericModel(set)
    "The ODE integrator for the full nonlinear model"
    integrator::Union{OrdinaryDiffEqCore.ODEIntegrator, Nothing} = nothing
    "The ODE integrator for the linearized model"
    lin_integ::Union{OrdinaryDiffEqCore.ODEIntegrator, Nothing} = nothing
    "Relative start time of the current time interval"
    t_0::SimFloat = 0.0
    "Number of solve! calls"
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

Overload `getproperty` to allow direct access to fields within the nested `serialized_model`.
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

Overload `setproperty!` to allow direct setting of fields within the nested `serialized_model`.
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

Constructs a `SymbolicAWEModel` from a `SystemStructure`.

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
    serialized_model = SerializedModel(; set_hash, sys_struct_hash)
    return SymbolicAWEModel(; set, sys_struct, serialized_model, kwargs...)
end

"""
    SymbolicAWEModel(set::Settings; kwargs...)

Constructs a default `SymbolicAWEModel` with automatically generated components.

This convenience constructor creates a complete AWE model using default configurations:
- Builds a `SystemStructure` based on the wing geometry and settings.
- Assembles everything into a ready-to-use symbolic model.

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
"""
function SymbolicAWEModel(set::Settings, name::String; kwargs...)
    set.physical_model = name
    sys_struct = SystemStructure(set; kwargs...)
    return SymbolicAWEModel(set, sys_struct)
end

"""
    update_sys_state!(ss::SysState, s::SymbolicAWEModel, zoom=1.0)

Update the high-level `SysState` log object from the detailed `SymbolicAWEModel` state.

This function extracts key performance indicators and state variables from the full
simulation model and populates the simpler `SysState` struct, which is used for
logging and visualization. It calculates derived quantities like roll, pitch, yaw,
angle of attack, and side slip.
"""
function update_sys_state!(ss::SysState, s::SymbolicAWEModel, zoom=1.0)
    ss.time = isnothing(s.integrator) ? 0.0 : s.integrator.t # Use integrator time
    @unpack points, groups, segments, pulleys, winches, wings = s.sys_struct

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
        # AoA: angle between v_app projected onto xz-plane and x-axis
        # Side Slip: angle between v_app and the xz-plane
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
        # roll (x-axis rotation)
        sinr_cosp = 2 * (q[1] * q[2] + q[3] * q[4])
        cosr_cosp = 1 - 2 * (q[2] * q[2] + q[3] * q[3])
        ss.roll = atan(sinr_cosp, cosr_cosp)
        # pitch (y-axis rotation)
        sinp = 2 * (q[1] * q[3] - q[4] * q[2])
        if abs(sinp) >= 1
            ss.pitch = copysign(pi / 2, sinp) # use 90 degrees if out of range
        else
            ss.pitch = asin(sinp)
        end
        # yaw (z-axis rotation)
        siny_cosp = 2 * (q[1] * q[4] + q[2] * q[3])
        cosy_cosp = 1 - 2 * (q[3] * q[3] + q[4] * q[4])
        ss.yaw = atan(siny_cosp, cosy_cosp)
    end
    for point in points
        ss.X[point.idx] = point.pos_w[1] * zoom
        ss.Y[point.idx] = point.pos_w[2] * zoom
        ss.Z[point.idx] = point.pos_w[3] * zoom
    end
    ss.v_wind_gnd .= s.sys_struct.wind_vec_gnd
    nothing
end

"""
    SysState(s::SymbolicAWEModel, zoom=1.0)

Create and populate a `SysState` logging object from the current state of a `SymbolicAWEModel`.
"""
function SysState(s::SymbolicAWEModel, zoom=1.0)
    ss = SysState{length(s.sys_struct.points)}()
    update_sys_state!(ss, s, zoom)
    ss
end

"""
    init!(s::SymbolicAWEModel; solver, adaptive, prn, precompile, remake, reload, lin_outputs) -> ODEIntegrator

Initialize a kite power system model, creating or loading a compiled version.

If a serialized model exists for the current configuration, it will load that model
(fast path). Otherwise, it will create a new model from scratch (slow path).

# Fast Path (serialized model exists):
1. Loads existing `ODEProblem` from disk.
2. Calls `reinit!` to update state variables from settings.
3. Sets up the ODE integrator.

# Slow Path (no serialized model):
1. Creates the full symbolic `ModelingToolkit.System` with all equations.
2. Simplifies the system equations.
3. Creates the `ODEProblem` and `LinearizationProblem`.
4. Serializes the compiled model to disk for future use.
5. Proceeds with the fast path.

# Arguments
- `s::SymbolicAWEModel`: The kite system state object.

# Keyword Arguments
- `solver`: Solver algorithm to use. If `nothing`, defaults based on `s.set.solver`.
- `adaptive::Bool=true`: Use adaptive time stepping.
- `prn::Bool=false`: Print progress information.
- `precompile::Bool=false`: Build a generic problem for precompilation.
- `remake::Bool=false`: Force the system to be rebuilt, even if a serialized model exists.
- `reload::Bool=false`: Force the system to reload the serialized model from disk.
- `lin_outputs::Vector{Num}=nothing`: List of symbolic variables for linearization.

# Returns
- `OrdinaryDiffEqCore.ODEIntegrator`: The initialized ODE integrator.
"""
function init!(s::SymbolicAWEModel; 
    solver=nothing, adaptive=true, prn=false, 
    precompile=false, remake=false, reload=false, 
    delta=nothing, stiffness_factor=nothing,
    lin_outputs=nothing
)
    if isnothing(solver)
        solver = if s.set.solver == "FBDF"
            if s.set.quasi_static
                FBDF(nlsolve=OrdinaryDiffEqNonlinearSolve.NLNewton(relax=s.set.relaxation))
            else
                FBDF()
            end
        elseif s.set.solver == "QNDF"
            @warn "This solver is not tested."
            QNDF()
        else
            error("Unavailable solver for SymbolicAWEModel: $(s.set.solver).")
        end
    end
    function init(s)
        reinit!(s.sys_struct, s.set)
        
        inputs = create_sys!(s, s.sys_struct; prn)
        prn && @info "Simplifying the System..."
        time = @elapsed @suppress_err begin
            sys = mtkcompile(s.full_sys; inputs, additional_passes = [ModelingToolkit.IfLifting])
            s.sys = sys
        end
        prn && @info "Simplified the System in $time seconds"

        prn && @info "Creating the ODEProblem..."
        dt = SimFloat(1/s.set.sample_freq)
        time = @elapsed s.prob = ODEProblem(s.sys, s.defaults, (0.0, dt); s.guesses)
        prn && @info "Created the ODEProblem in $time seconds"

        if isnothing(lin_outputs)
            lin_outputs = Num[]
            if length(s.sys_struct.wings) > 0
                push!(lin_outputs, s.sys.heading[1])
                push!(lin_outputs, s.sys.angle_of_attack[1])
            end
            if length(s.sys_struct.winches) > 0
                push!(lin_outputs, s.sys.tether_len[1])
                push!(lin_outputs, s.sys.winch_force[1])
            end
        end
        prn && @info "Creating the LinearizationProblem..."
        time = @elapsed @suppress_err begin
            lin_fun, _ = linearization_function(s.full_sys, [inputs...], lin_outputs; op=s.defaults, guesses=s.guesses)
            s.lin_prob = LinearizationProblem(lin_fun, 0.0)
            s.lin_outputs = lin_outputs
        end
        prn && @info "Created the LinearizationProblem in $time seconds"

        generate_getters!(s, s.lin_outputs)
        s.set_hash = get_set_hash(s.set)
        s.sys_struct_hash = get_sys_struct_hash(s.sys_struct)
        serialize(model_path, s.serialized_model)
        s.integrator = nothing
        return nothing
    end
    model_path = joinpath(KiteUtils.get_data_path(), get_model_name(s.set; precompile))
    if !ispath(model_path) || remake
        init(s)
    end
    _, success = reinit!(s, solver; adaptive, precompile, reload, s.lin_outputs, prn)
    if !success
        rm(model_path)
        @info "Rebuilding the system. This can take some minutes..."
        init(s)
        reinit!(s, solver; adaptive, precompile, s.lin_outputs, prn, reload=true)
    end
    return s.integrator
end


"""
    reinit!(s::SymbolicAWEModel, solver; prn, precompile, reload, lin_outputs) -> (ODEIntegrator, Bool)

Reinitialize an existing kite power system model with new state values from `s.set`.

This function performs the following operations:
1. If no integrator exists yet, or `reload=true`:
   - Loads a serialized `ODEProblem` from disk.
   - Initializes a new ODE integrator.
   - Checks for hash mismatches in settings or structure.
2. Initializes the `SystemStructure` with new positions and orientations.
3. Sets initial values for all state variables in the integrator.
4. Reinitializes the ODE integrator with the new state.
5. Updates the linearized aerodynamic model.

This is more efficient than `init!` as it reuses the existing model structure.

# Arguments
- `s::SymbolicAWEModel`: The kite power system state object.
- `solver`: The solver to be used.
- `prn::Bool=false`: Whether to print progress information.
- `precompile::Bool=false`: Load the precompiled version of the model.
- `reload::Bool=true`: Force reloading the model from disk.
- `lin_outputs::Vector{Num}=Num[]`: Outputs for the linearized model.

# Returns
- `(ODEIntegrator, Bool)`: A tuple containing the reinitialized integrator and a success flag.
"""
function reinit!(
    s::SymbolicAWEModel,
    solver;
    adaptive=true,
    prn=false, 
    reload=true, 
    precompile=false,
    lin_outputs=Num[]
)
    isnothing(s.sys_struct) && error("SystemStructure not defined")

    if isnothing(s.prob) || reload
        model_path = joinpath(KiteUtils.get_data_path(), get_model_name(s.set; precompile))
        if !ispath(model_path)
            error("$model_path not found. Run init!(s::SymbolicAWEModel) first.")
        end
        if prn
            @info "Loading model from $model_path"
        end # model_1.11_ram_dynamic_3_seg.bin
        try
            s.serialized_model = deserialize(model_path)
        catch e
            @warn "Failure to deserialize $model_path !"
            return s.integrator, false
        end
        if isnothing(lin_outputs)
            sys = s.sys
            lin_outputs = Num[]
            if length(s.sys_struct.wings) > 0
                push!(lin_outputs, sys.heading[1])
                push!(lin_outputs, sys.angle_of_attack[1])
            end
            if length(s.sys_struct.winches) > 0
                push!(lin_outputs, sys.tether_len[1])
                push!(lin_outputs, sys.winch_force[1])
            end
        end
        if length(lin_outputs) != length(s.serialized_model.lin_outputs) ||
                !all(string.(lin_outputs) .== string.(s.serialized_model.lin_outputs)) 
            @warn "The linear model outputs have changed."
            return s.integrator, false
        elseif (get_set_hash(s.set) != s.serialized_model.set_hash)
            @warn "The Settings have changed."
            return s.integrator, false
        elseif (get_sys_struct_hash(s.sys_struct) != s.serialized_model.sys_struct_hash)
            @warn "The SystemStructure has changed."
            return s.integrator, false
        s.get_distance = (integ) -> get_distance(integ)
        end
    end
    if isnothing(s.integrator) || !successful_retcode(s.integrator.sol) || reload
        t = @elapsed begin
            dt = SimFloat(1/s.set.sample_freq)
            s.sys = s.prob.f.sys
            s.integrator = OrdinaryDiffEqCore.init(s.prob, solver; 
                adaptive, dt, abstol=s.set.abs_tol, reltol=s.set.rel_tol, 
                save_on=false, save_everystep=false)
            s.lin_integ = OrdinaryDiffEqCore.init(s.prob, solver; 
                adaptive, dt, abstol=s.set.abs_tol, reltol=s.set.rel_tol, 
                save_on=false, save_everystep=false)
        end
        prn && @info "Initialized integrator in $t seconds"
    end

    reinit!(s.sys_struct, s.set)
    s.set_psys(s.integrator, s.sys_struct)
    s.set_set(s.integrator, s.set)
    OrdinaryDiffEqCore.reinit!(s.integrator; reinit_dae=true)
    linearize_vsm!(s)
    update_sys_struct!(s, s.sys_struct)
    return s.integrator, true
end

"""
    generate_getters!(s, lin_y_vec)

Generate and compile optimized getter and setter functions for the model.

This internal function uses the symbolic system definition from `ModelingToolkit.jl`
to create fast, non-allocating functions for accessing and modifying the system's
state and parameters directly within the ODE integrator's data structures. This is
a key optimization that avoids symbolic lookups during the simulation loop.
"""
function generate_getters!(s, lin_y_vec)
    sys = s.sys
    c = collect
    @unpack wings, groups, pulleys, winches, tethers, segments = s.sys_struct

    if length(wings) == 1
        lin_x_vec = [
            sys.heading[1]
            sys.turn_rate[1,3]
            sys.tether_len[1]
            sys.tether_len[2]
            sys.tether_len[3]
            sys.tether_vel[1]
            sys.tether_vel[2]
            sys.tether_vel[3]
        ]
        lin_dx_vec = [
            sys.turn_rate[1,3]
            sys.turn_acc[1,3]
            sys.tether_vel[1]
            sys.tether_vel[2]
            sys.tether_vel[3]
            sys.tether_acc[1]
            sys.tether_acc[2]
            sys.tether_acc[3]
        ]
        nx = length(lin_x_vec)
        ny = length(lin_y_vec)
        nu = length(winches)
        s.simple_lin_model = (
            A = zeros(nx, nx),
            B = zeros(nx, nu),
            C = zeros(ny, nx),
            D = zeros(ny, nu)
        ) 
        get_lin_x = getu(sys, lin_x_vec)
        s.get_lin_x = (integ) -> get_lin_x(integ)
        get_lin_dx = getu(sys, lin_dx_vec)
        s.get_lin_dx = (integ) -> get_lin_dx(integ)
        get_lin_y = getu(sys, lin_y_vec)
        s.get_lin_y = (integ) -> get_lin_y(integ)
    end

    if length(wings) > 0
        get_wing_state = getu(sys, c.([
            sys.Q_b_w,           # Orientation quaternion
            sys.ω_b,             # Angular velocity (body frame)
            sys.wing_pos,         # Position vector (world frame)
            sys.wing_vel,         # Velocity vector (world frame)
            sys.wing_acc,
            sys.va_wing_b,           # Apparent wind vector (body frame)
            sys.wind_vel_wing,         # Wind vector (body frame)
            sys.aero_force_b,   # Aerodynamic force vector (body frame)
            sys.aero_moment_b,  # Aerodynamic moment vector (body frame)
            sys.elevation,      # Elevation angle
            sys.elevation_vel,
            sys.elevation_acc,
            sys.azimuth,       # Azimuth angle
            sys.azimuth_vel,
            sys.azimuth_acc,
            sys.heading,        # Heading angle
            sys.turn_rate,
            sys.turn_acc,
            sys.course,         # Course angle
            sys.angle_of_attack,
        ]))
        s.get_wing_state = (integ) -> get_wing_state(integ)
        get_vsm_y = getu(sys, sys.y)
        s.get_vsm_y = (integ) -> get_vsm_y(integ)
    end

    if length(segments) > 0
        get_segment_state = getu(sys, c.([
            sys.spring_force,
            sys.len,
        ]))
        s.get_segment_state = (integ) -> get_segment_state(integ)
    end

    if length(groups) > 0
        get_group_state = getu(sys, c.([
            sys.twist_angle,     # Twist angle per group
            sys.twist_ω,       # Twist velocity per group
            sys.group_tether_force,
            sys.group_tether_moment,
            sys.group_aero_moment,
        ]))
        s.get_group_state = (integ) -> get_group_state(integ)
    end
    
    if length(pulleys) > 0
        get_pulley_state = getu(sys, c.([
            sys.pulley_len,      # Position vector (world frame)
            sys.pulley_vel,      # Velocity vector (world frame)
        ]))
        s.get_pulley_state = (integ) -> get_pulley_state(integ)
    end

    if length(winches) > 0
        get_winch_state = getu(sys, c.([
             sys.tether_len,   # Unstretched len per winch
             sys.tether_vel,      # Reeling velocity per winch
             sys.set_values,
             sys.winch_force_vec,     # Force at winch connection point per winch
        ]))
        s.get_winch_state = (integ) -> get_winch_state(integ)
        get_set_values = getp(sys, sys.set_values)
        s.get_set_values = (integ) -> get_set_values(integ)

        set_set_values = setp(sys, sys.set_values)
        s.set_set_values = (integ, val) -> set_set_values(integ, val)
        if !isnothing(s.lin_prob)
            set_lin_set_values = setp(s.lin_prob, sys.set_values)
            s.set_lin_set_values = (lin_prob, val) -> set_lin_set_values(lin_prob, val)
        end
    end

    if length(tethers) > 0
        get_tether_state = getu(sys, c(sys.stretched_len))
        s.get_tether_state = (integ) -> get_tether_state(integ)
    end

    set_psys = setp(sys, sys.psys)
    s.set_psys = (integ, val) -> set_psys(integ, val)
    set_set = setp(sys, sys.pset)
    s.set_set = (integ, val) -> set_set(integ, val)
    
    get_point_state = getu(sys, c.([
         sys.pos,             # Particle positions
         sys.vel,             # Kite center acceleration vector (world frame)
         sys.point_force
    ]))
    s.get_point_state = (integ) -> get_point_state(integ)
    get_spring_force = getu(sys, sys.spring_force)
    s.get_spring_force = (integ) -> get_spring_force(integ)
    get_struct_state = getu(sys, sys.wind_vec_gnd)
    s.get_struct_state = (integ) -> get_struct_state(integ)
    nothing
end

"""
    next_step!(s::SymbolicAWEModel, integrator::ODEIntegrator; set_values, dt, vsm_interval)

Take a simulation step, using the provided integrator.
"""
function KiteUtils.next_step!(s::SymbolicAWEModel, integrator::OrdinaryDiffEqCore.ODEIntegrator; set_values=nothing, dt=1/s.set.sample_freq, vsm_interval=1)
    !(s.integrator === integrator) && error("The ODEIntegrator doesn't belong to the SymbolicAWEModel")
    next_step!(s; set_values, upwind_dir, dt, vsm_interval)
end

"""
    next_step!(s::SymbolicAWEModel; set_values, dt, vsm_interval)

Take a simulation step forward in time.

This function performs the following steps:
1. Optionally updates the control inputs (`set_values`).
2. Optionally linearizes the VSM (Vortex Step Method) aerodynamic model.
3. Steps the ODE integrator forward by `dt` seconds.
4. Checks for a successful return code from the integrator.
5. Updates the high-level `SystemStructure` with the new state from the integrator.
6. Increments the iteration counter.

# Arguments
- `s::SymbolicAWEModel`: The kite power system state object.

# Keyword Arguments
- `set_values=nothing`: New values for the control inputs. If `nothing`, the current values are used.
- `dt=1/s.set.sample_freq`: Time step size [s].
- `vsm_interval=1`: Interval (in steps) to re-linearize the VSM model. If 0, it is not re-linearized.

# Returns
- `Nothing`
"""
function KiteUtils.next_step!(s::SymbolicAWEModel; set_values=nothing, dt=1/s.set.sample_freq, vsm_interval=1)
    if (!isnothing(set_values)) 
        s.set_set_values(s.integrator, set_values)
    end
    if vsm_interval != 0 && s.iter % vsm_interval == 0
        s.t_vsm = @elapsed linearize_vsm!(s)
    end
    
    s.t_0 = s.integrator.t
    s.t_step = @elapsed OrdinaryDiffEqCore.step!(s.integrator, dt, true)
    if !successful_retcode(s.integrator.sol)
        @warn "Return code for solution: $(s.integrator.sol.retcode)"
    end
    @assert successful_retcode(s.integrator.sol)
    s.iter += 1
    update_sys_struct!(s, s.sys_struct)
    return nothing
end

"""
    update_sys_struct!(s::SymbolicAWEModel, sys_struct::SystemStructure, integ=s.integrator)

Update the high-level `SystemStructure` from the low-level integrator state vector.

This function is the reverse of `update_sys_state!`. It reads the raw state vector
from the ODE integrator and uses the generated getter functions to populate the
human-readable fields in the `SystemStructure` (e.g., `point.pos_w`, `winch.tether_len`).
This is called after each simulation step to keep the high-level struct synchronized.
"""
function update_sys_struct!(s::SymbolicAWEModel, sys_struct::SystemStructure, integ=s.integrator)
    @unpack points, groups, segments, pulleys, winches, tethers, wings = sys_struct
    pos, vel, force = s.get_point_state(integ)
    for point in points
        point.pos_w .= pos[:, point.idx]
        point.vel_w .= vel[:, point.idx]
        point.force .= force[:, point.idx]
    end
    if length(pulleys) > 0
        len, vel = s.get_pulley_state(integ)
        for pulley in pulleys
            pulley.len = len[pulley.idx]
            pulley.vel = vel[pulley.idx]
        end
    end
    if length(segments) > 0
        spring_force, len = s.get_segment_state(integ)
        for segment in segments
            segment.force = spring_force[segment.idx]
            segment.len = len[segment.idx]
        end
    end
    if length(groups) > 0
        twist, twist_ω, tether_force, tether_moment, aero_moment = s.get_group_state(integ)
        for group in groups
            group.twist = twist[group.idx]
            group.twist_ω = twist_ω[group.idx]
            group.tether_force = tether_force[group.idx]
            group.tether_moment = tether_moment[group.idx]
            group.aero_moment = aero_moment[group.idx]
        end
    end
    if length(winches) > 0
        tether_len, tether_vel, set_value, winch_force_vec = s.get_winch_state(integ)
        for winch in winches
            winch.tether_len = tether_len[winch.idx]
            winch.tether_vel = tether_vel[winch.idx]
            winch.set_value = set_value[winch.idx]
            winch.force .= winch_force_vec[:, winch.idx]
        end
    end
    if length(tethers) > 0
        stretched_len = s.get_tether_state(integ)
        for tether in tethers
            tether.stretched_len = stretched_len[tether.idx]
        end
    end
    if length(wings) > 0
        Q_b_w, ω_b, pos_w, vel_w, acc_w, va_b, v_wind, 
            aero_force_b, aero_moment_b, elevation, elevation_vel,
            elevation_acc, azimuth, azimuth_vel, azimuth_acc,
            heading, turn_rate, turn_acc, course, aoa = s.get_wing_state(integ)
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
    s.sys_struct.wind_vec_gnd .= s.get_struct_state(integ)
    return nothing
end

"""
    get_model_name(set::Settings; precompile=false)

Construct a unique filename for the serialized model based on its configuration.

The filename encodes the Julia version, physical model name, dynamics type, and
number of segments, ensuring that different model configurations are cached separately.
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

Calculate the mean angle of attack [rad] over the wingspan from the VSM solver.
"""
function calc_aoa(s::SymbolicAWEModel)
    alpha_array = s.sys_struct.wings[1].vsm_solver.sol.alpha_array
    middle = length(alpha_array) ÷ 2
    if iseven(length(alpha_array))
        return 0.5alpha_array[middle] + 0.5alpha_array[middle+1]
    else
        return alpha_array[middle+1]
    end
end

"""Returns the unstretched tether length [m] for each winch."""
unstretched_length(s::SymbolicAWEModel) = [winch.tether_len for winch in s.sys_struct.winches]

"""Returns the current tether length of the symbolic AWE model."""
tether_length(s::SymbolicAWEModel) = [winch.tether_len for winch in s.sys_struct.winches]

"""Returns the height (z-position) of the wing in the symbolic AWE model."""
calc_height(s::SymbolicAWEModel) = [wing.pos_w[3] for wing in s.sys_struct.wings]

"""Returns the winch force [N] for each winch."""
winch_force(s::SymbolicAWEModel) = [norm(winch.force) for winch in s.sys_struct.winches]

"""Returns the spring force [N] for each segment."""
spring_forces(s::SymbolicAWEModel) = [segment.force for segment in s.sys_struct.segments]

"""Returns the position vector [m] for each point."""
function pos(s::SymbolicAWEModel)
    return [point.pos_w for point in s.sys_struct.points]
end    

"""
    min_chord_len(s::SymbolicAWEModel)

Calculate the minimum chord length of the wing at the tip.
"""
function min_chord_len(s::SymbolicAWEModel)
    min_len = Inf
    for wing in s.sys_struct.wings
        vsm_wing = wing.vsm_wing
        le_pos = [vsm_wing.le_interp[i](vsm_wing.gamma_tip) for i in 1:3]
        te_pos = [vsm_wing.te_interp[i](vsm_wing.gamma_tip) for i in 1:3]
        min_len = min(norm(le_pos - te_pos), min_len)
    end
    return min_len
end

"""
    set_depower_steering!(s::SymbolicAWEModel, depower, steering)

Set the kite's depower and steering by adjusting the tether length set-points.

Depower and steering values are scaled by the minimum chord length of the wing.
"""
function set_depower_steering!(s::SymbolicAWEModel, depower, steering)
    len = s.set_tether_len
    len .= tether_length(s)
    depower *= min_chord_len(s)
    steering *= min_chord_len(s)
    len[2] = 0.5 * (2*depower + 2*len[1] + steering)
    len[3] = 0.5 * (2*depower + 2*len[1] - steering)
    return nothing
end

"""
    set_v_wind_ground!(s::SymbolicAWEModel, v_wind_gnd=s.set.v_wind, upwind_dir=-π/2)

Set the ground wind speed [m/s] and upwind direction [rad].

Direction convention: 0 for North, π/2 for East, π for South, -π/2 for West (default).
"""
function set_v_wind_ground!(s::SymbolicAWEModel, v_wind_gnd=s.set.v_wind, upwind_dir=-pi/2)
    s.set.v_wind = v_wind_gnd
    s.set.upwind_dir = rad2deg(upwind_dir)
    s.set_set(s.integrator, s.set)
    return nothing
end

"""
    get_set_hash(set::Settings; fields)

Calculate a SHA1 hash for a subset of fields in the `Settings` object.

This is used to detect changes in the configuration that would require
rebuilding the serialized model.
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

Calculate a SHA1 hash for the topology of a `SystemStructure`.

This hash uniquely identifies the number and connectivity of points, segments,
groups, etc. It is used to detect changes in the system's structure that would
require rebuilding the serialized model.
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

