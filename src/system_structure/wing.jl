# Copyright (c) 2025 Bart van de Lint
# SPDX-License-Identifier: MPL-2.0

"""
Wing types and VSM-related wing code.

This file contains:
- AbstractWing, BaseWing, VSMWing structs
- Wing constructors and helper functions
- VSM panel adjustment functions
"""

# ==================== ABSTRACT WING ==================== #

"""
    abstract type AbstractWing

Abstract base type for all wing implementations.

Concrete subtypes must implement rigid body dynamics and provide a reference frame
for attached points and groups.
"""
abstract type AbstractWing end

# ==================== BASE WING ==================== #

"""
    mutable struct BaseWing <: AbstractWing

A rigid wing body that can have multiple groups of points attached to it.

The wing provides a rigid body reference frame for attached points and groups.
Points with `type == WING` move rigidly with the wing body according to the
wing's orientation matrix `R_b_w` and position `pos_w`.

# Special Properties
The wing's orientation can be accessed as a rotation matrix or a quaternion:
```julia
R_matrix = wing.R_b_w
wing.R_b_w = R_matrix

quat = wing.Q_b_w
wing.Q_b_w = quat
```

$(TYPEDFIELDS)
"""
mutable struct BaseWing <: AbstractWing
    const idx::Int64

    # Structural information
    group_idxs::Vector{Int64}
    const transform_idx::Int64
    const R_b_c::Matrix{SimFloat}
    const pos_cad::KVec3
    const inertia_principal::KVec3
    const wing_type::WingType

    # Differential variables in world frame, updated during simulation
    const Q_b_w::Vector{SimFloat}
    const ω_b::KVec3
    const pos_w::KVec3
    const vel_w::KVec3
    const acc_w::KVec3

    # Derived variables and parameters, updated during simulation
    const wind_disturb::KVec3
    drag_frac::SimFloat
    const va_b::KVec3 # apparent wind in body frame
    const v_wind::KVec3 # wind velocity in world frame at the wing
    const aero_force_b::KVec3 # aerodynamic force in body frame
    const aero_moment_b::KVec3 # aerodynamic moment in body frame
    const tether_moment::KVec3 # tether moment in world frame
    const tether_force::KVec3 # tether force in world frame
    elevation::SimFloat
    elevation_vel::SimFloat
    elevation_acc::SimFloat
    azimuth::SimFloat
    azimuth_vel::SimFloat
    azimuth_acc::SimFloat
    heading::SimFloat
    const turn_rate::KVec3
    const turn_acc::KVec3
    course::SimFloat
    aoa::SimFloat
    fix_sphere::Bool
    y_damping::SimFloat
    z_disturb::SimFloat
end

function Base.getproperty(wing::BaseWing, sym::Symbol)
    if sym == :R_b_w
        return quaternion_to_rotation_matrix(wing.Q_b_w)
    else
        return getfield(wing, sym)
    end
end

function Base.setproperty!(wing::BaseWing, sym::Symbol, value)
    if sym == :R_b_w
        wing.Q_b_w .= rotation_matrix_to_quaternion(value)
    else
        setfield!(wing, sym, value)
    end
end

# ==================== VSM WING ==================== #

"""
    mutable struct VSMWing <: AbstractWing

A wing that uses the Vortex Step Method (VSM) for aerodynamic computations.

This struct extends the base wing functionality with VSM-specific aerodynamic
modeling capabilities, including vortex wake computations and aerodynamic loads.

$(TYPEDFIELDS)
"""
mutable struct VSMWing <: AbstractWing
    # Base wing functionality
    base::BaseWing

    # VSM aerodynamics
    vsm_aero::VortexStepMethod.BodyAerodynamics
    vsm_wing::VortexStepMethod.AbstractWing
    vsm_solver::VortexStepMethod.Solver

    # VSM state and linearization
    vsm_y::Vector{SimFloat}
    vsm_x::Vector{SimFloat}
    vsm_jac::Matrix{SimFloat}

    # REFINE-specific fields (Nothing for QUATERNION wings)
    point_to_vsm_point::Union{Nothing, Dict{Int64, Tuple{Int64, Symbol}}}
    wing_segments::Union{Nothing,
        Vector{Tuple{Int64, Int64}}}

    # Orientation reference points for REFINE wings
    # (Nothing for QUATERNION wings)
    # Used to calculate R_b_w from structural deformation
    # Can specify single point or vector of points to average:
    #   (12, 13) - point 12 to point 13
    #   (12, [13, 14]) - point 12 to average of points 13,14
    #   ([11, 12], [13, 14]) - average of 11,12 to average of 13,14
    # Z-axis: Normal to wing plane, Y-axis: Spanwise, X = Y × Z (chord)
    z_ref_points::Union{Nothing, Tuple{Union{Int64, Vector{Int64}}, Union{Int64, Vector{Int64}}}}
    y_ref_points::Union{Nothing, Tuple{Union{Int64, Vector{Int64}}, Union{Int64, Vector{Int64}}}}

    # KCU origin point for REFINE wings
    # (Nothing for QUATERNION wings)
    # Defines wing.pos_w = pos[:, origin_idx] to track structural deformation
    origin_idx::Union{Nothing, Int64}

    # Additional aerodynamic force scale to compensate chord length errors (REFINE)
    aero_scale_chord::SimFloat

    # Body frame z-axis offset for VSM aerodynamics (QUATERNION only)
    # Shifts VSM panel positions in positive z direction (body frame)
    # to adjust moment arm for improved stability
    aero_z_offset::SimFloat

    function VSMWing(base::BaseWing, vsm_aero, vsm_wing, vsm_solver, vsm_y, vsm_x, vsm_jac, point_to_vsm_point, wing_segments, z_ref_points, y_ref_points, origin_idx, aero_scale_chord, aero_z_offset)
        new(base, vsm_aero, vsm_wing, vsm_solver, vsm_y, vsm_x, vsm_jac, point_to_vsm_point, wing_segments, z_ref_points, y_ref_points, origin_idx, aero_scale_chord, aero_z_offset)
    end
end

# Delegate property access to base wing for VSMWing
function Base.getproperty(wing::VSMWing, sym::Symbol)
    if sym in (:base, :vsm_aero, :vsm_wing, :vsm_solver, :vsm_y, :vsm_x, :vsm_jac, :point_to_vsm_point, :wing_segments, :z_ref_points, :y_ref_points, :origin_idx, :aero_scale_chord, :aero_z_offset)
        return getfield(wing, sym)
    elseif sym == :vsm_aoa
        # Compute mean angle of attack from VSM solver solution
        solver = getfield(wing, :vsm_solver)
        return mean(solver.sol.alpha_array)
    else
        return getproperty(getfield(wing, :base), sym)
    end
end

function Base.setproperty!(wing::VSMWing, sym::Symbol, value)
    if sym in (:base, :vsm_aero, :vsm_wing, :vsm_solver, :vsm_y, :vsm_x, :vsm_jac, :point_to_vsm_point, :wing_segments, :z_ref_points, :y_ref_points, :origin_idx, :aero_scale_chord, :aero_z_offset)
        setfield!(wing, sym, value)
    else
        setproperty!(getfield(wing, :base), sym, value)
    end
end

# ==================== CONSTRUCTORS ==================== #

"""
    BaseWing(idx::Int64, group_idxs::Vector{Int64}, R_b_c::Matrix{SimFloat},
             pos_cad::KVec3, inertia_principal::KVec3; transform_idx=1, y_damping=150.0, wing_type=QUATERNION)

Constructs a `BaseWing` object representing a rigid body reference frame.

# Arguments
- `idx::Int64`: Unique identifier for the wing.
- `group_idxs::Vector{Int64}`: Indices of groups attached to this wing.
- `R_b_c::Matrix{SimFloat}`: Rotation matrix from body frame to CAD frame.
- `pos_cad::KVec3`: Position of wing center of mass in CAD frame.
- `inertia_principal::KVec3`: Principal moments of inertia [Ixx, Iyy, Izz] in body frame.

# Keyword Arguments
- `transform_idx::Int64=1`: Transform used for initial positioning and orientation.
- `y_damping::SimFloat=150.0`: Damping coefficient for lateral motion.
- `wing_type::WingType=QUATERNION`: Wing aerodynamic model type.

# Returns
- `BaseWing`: A new base wing object.
"""
function BaseWing(idx, group_idxs::AbstractVector, R_b_c::AbstractMatrix,
                  pos_cad, inertia_principal; transform_idx=1, y_damping=150.0, wing_type::WingType=QUATERNION)
    return BaseWing(idx,
        # Structural information
        group_idxs, transform_idx, R_b_c, pos_cad, inertia_principal, wing_type,
        # Differential variables in world frame, updated during simulation
        zeros(4), zeros(KVec3),
        zeros(KVec3), zeros(KVec3), zeros(KVec3),
        # Derived variables and parameters, updated during simulation
        zeros(KVec3), one(SimFloat),
        zeros(KVec3), zeros(KVec3), zeros(KVec3), zeros(KVec3),
        zeros(KVec3), zeros(KVec3),
        0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
        zeros(KVec3), zeros(KVec3), 0.0, 0.0, false,
        y_damping, 0.0)
end

"""
    VSMWing(idx::Int64, set::Settings, group_idxs::Vector{Int64},
            R_b_c::Matrix{SimFloat}, pos_cad::KVec3;
            transform_idx=1, y_damping=150.0,
            wing_type=QUATERNION, point_to_vsm_point=nothing,
            wing_segments=nothing, x_ref_points=nothing,
            y_ref_points=nothing)

Constructs a `VSMWing` object with Vortex Step Method aerodynamics.
Creates vsm_wing, vsm_aero, and vsm_solver internally.

# Arguments
- `idx::Int64`: Unique identifier for the wing.
- `set::Settings`: Settings object for VSM configuration.
- `group_idxs::Vector{Int64}`: Indices of groups (QUATERNION only).
- `R_b_c::Matrix{SimFloat}`: Rotation matrix body→CAD.
- `pos_cad::KVec3`: Position of wing COM in CAD frame.

# Keyword Arguments
- `transform_idx::Int64=1`: Transform for initial positioning.
- `y_damping::SimFloat=150.0`: Lateral damping coefficient.
- `wing_type::WingType=QUATERNION`: Aerodynamic model type.
- `point_to_vsm_point`: 1:1 structural point to VSM point mapping (REFINE only).
- `wing_segments`: LE/TE pairs (REFINE only).
- `x_ref_points`: Chord direction reference (REFINE only).
- `y_ref_points`: Span direction reference (REFINE only).
- `aero_z_offset::SimFloat=0.0`: Body frame z-offset for VSM panels (QUATERNION only).

# Returns
- `VSMWing`: A new VSM wing object.
"""
function VSMWing(idx::Int, set::Settings,
                 group_idxs::AbstractVector,
                 vsm_set::VortexStepMethod.VSMSettings;
                 R_b_c::Union{Nothing,AbstractMatrix}=nothing,
                 pos_cad::Union{Nothing,AbstractVector}=nothing,
                 transform_idx=1, y_damping=150.0,
                 inertia_diag=nothing,
                 wing_type::WingType=QUATERNION,
                 point_to_vsm_point::Union{Nothing, Dict{Int64, Tuple{Int64, Symbol}}}=nothing,
                 wing_segments::Union{Nothing,
                     Vector{Tuple{Int64, Int64}}}=nothing,
                 z_ref_points::Union{Nothing,
                     Tuple{Union{Int64, Vector{Int64}}, Union{Int64, Vector{Int64}}}}=nothing,
                 y_ref_points::Union{Nothing,
                     Tuple{Union{Int64, Vector{Int64}}, Union{Int64, Vector{Int64}}}}=nothing,
                 origin_idx::Union{Nothing, Int64}=nothing,
                 aero_scale_chord::SimFloat=0.0,
                 aero_z_offset::SimFloat=0.0)

    # Validation
    if wing_type == REFINE
        @assert length(group_idxs) == 0
            "REFINE wings cannot have groups"
        @assert !isnothing(origin_idx)
            "REFINE wings require origin_idx to define KCU position"
    else
        @assert isnothing(point_to_vsm_point)
            "QUATERNION wings: no point_to_vsm_point"
        @assert isnothing(wing_segments)
            "QUATERNION wings: no wing_segments"
        @assert isnothing(z_ref_points)
            "QUATERNION wings: no z_ref_points"
        @assert isnothing(y_ref_points)
            "QUATERNION wings: no y_ref_points"
        @assert isnothing(origin_idx)
            "QUATERNION wings don't use origin_idx"
    end

    # Create VSM wing, aero, and solver
    vsm_wing = VortexStepMethod.Wing(set, vsm_set; prn=false)
    vsm_aero = VortexStepMethod.BodyAerodynamics([vsm_wing])
    vsm_solver = VortexStepMethod.Solver(vsm_aero;
        solver_type=VortexStepMethod.NONLIN,
        atol=2e-8, rtol=2e-8)

    # Set defaults from actual vsm_wing if not provided
    if isnothing(R_b_c) || isnothing(pos_cad)
        isnothing(R_b_c) && (R_b_c = vsm_wing.R_cad_body)
        isnothing(pos_cad) && (pos_cad = vsm_wing.T_cad_body)
    end

    # Compute inertia
    inertia_vec = isnothing(inertia_diag) ?
        wing_inertia_principal(vsm_wing) : inertia_diag

    base = BaseWing(idx, group_idxs, R_b_c, pos_cad,
                    inertia_vec; transform_idx,
                    y_damping, wing_type)

    # Size vsm state vectors based on wing type
    if wing_type == REFINE
        nx = 3 * length(vsm_aero.panels)
        ny = 0
    else
        # QUATERNION: use number of unrefined sections
        n_unrefined = vsm_wing.n_unrefined_sections
        ny = 3 + n_unrefined + 3  # va(3) + twist(n_unrefined) + ω(3)
        nx = 3 + 3 + n_unrefined  # force(3) + moment(3) + unrefined_moments(n_unrefined)
    end

    return VSMWing(base, vsm_aero, vsm_wing, vsm_solver,
                   zeros(SimFloat, ny), zeros(SimFloat, nx),
                   zeros(SimFloat, nx, ny),
                   point_to_vsm_point, wing_segments,
                   z_ref_points, y_ref_points, origin_idx, aero_scale_chord,
                   aero_z_offset)
end

"""
    VSMWing(idx, vsm_aero, vsm_wing, vsm_solver,
            group_idxs, R_b_c, pos_cad)

Legacy constructor accepting pre-created VSM objects directly.
Kept for backward compatibility with predefined structures.

# Arguments
- `idx::Int`: Wing identifier
- `vsm_aero`: Pre-created BodyAerodynamics
- `vsm_wing`: Pre-created VortexStepMethod.Wing
- `vsm_solver`: Pre-created Solver
- `group_idxs`: Group indices
- `R_b_c`: Rotation matrix body→CAD
- `pos_cad`: Position in CAD frame

# Returns
- `VSMWing`: Wing with QUATERNION type
"""
function VSMWing(idx::Int, vsm_aero, vsm_wing, vsm_solver,
                 group_idxs::AbstractVector,
                 R_b_c::AbstractMatrix,
                 pos_cad::AbstractVector)
    inertia_vec = wing_inertia_principal(vsm_wing)
    base = BaseWing(idx, group_idxs, R_b_c, pos_cad, inertia_vec)
    # Use number of unrefined sections
    n_unrefined = vsm_wing.n_unrefined_sections
    ny = 3 + n_unrefined + 3  # va(3) + twist(n_unrefined) + ω(3)
    nx = 3 + 3 + n_unrefined  # force(3) + moment(3) + unrefined_moments(n_unrefined)
    return VSMWing(base, vsm_aero, vsm_wing, vsm_solver,
        zeros(SimFloat, ny), zeros(SimFloat, nx),
        zeros(SimFloat, nx, ny),
        nothing, nothing, nothing, nothing, nothing, 0.0, 0.0)
end

"""
    Wing(idx, vsm_aero, vsm_wing, vsm_solver, group_idxs, R_b_c, pos_cad; transform_idx)

Constructs a `VSMWing` object (backward compatibility constructor).

This is a convenience constructor that creates a VSMWing for backward compatibility
with existing code. New code should use `VSMWing(...)` directly.

# Arguments
- `idx::Int64`: Unique identifier for the wing.
- `vsm_aero`, `vsm_wing`, `vsm_solver`: Vortex Step Method components.
- `group_idxs::Vector{Int64}`: Indices of groups attached to this wing.
- `R_b_c::Matrix{SimFloat}`: Rotation matrix from body frame to CAD frame.
- `pos_cad::KVec3`: Position of wing center of mass in CAD frame.

# Keyword Arguments
- `transform_idx::Int64=1`: Transform used for initial positioning and orientation.
- `y_damping::SimFloat=150.0`: Damping coefficient for lateral motion.

# Returns
- `VSMWing`: A new VSM wing object.
"""
function SymbolicAWEModels.Wing(idx, vsm_aero, vsm_wing, vsm_solver, group_idxs, R_b_c,
                                pos_cad; kwargs...)
    return VSMWing(idx, vsm_aero, vsm_wing, vsm_solver, group_idxs, R_b_c, pos_cad; kwargs...)
end

# ==================== HELPER FUNCTIONS ==================== #

"""
    wing_inertia_principal(vsm_wing)

Extract principal moments of inertia from a VortexStepMethod wing.
Returns diagonal of inertia tensor if available, otherwise returns ones.
"""
function wing_inertia_principal(vsm_wing)
    if hasproperty(vsm_wing, :inertia_tensor) && size(vsm_wing.inertia_tensor) == (3, 3)
        diag_vals = diag(vsm_wing.inertia_tensor)
        return MVector{3, SimFloat}(diag_vals)
    end
    return MVector{3, SimFloat}(ones(SimFloat, 3))
end

"""
    adjust_vsm_panels_to_origin!(vsm_wing, origin_offset)

Adjust VSM panel positions when body frame origin changes.

When QUATERNION wings are loaded from YAML, the panel positions in aero_geometry.yaml
are specified in an absolute body frame. However, the body frame origin is adjusted
to the mean position of all WING points. This function updates all panel positions
to be relative to the new origin by subtracting the offset.

# Arguments
- `vsm_wing`: VortexStepMethod.Wing with sections to adjust
- `origin_offset`: Vector [x, y, z] to subtract from panel positions
"""
function adjust_vsm_panels_to_origin!(vsm_wing, origin_offset)
    for section in vsm_wing.refined_sections
        section.LE_point .-= origin_offset
        section.TE_point .-= origin_offset
    end
    for section in vsm_wing.non_deformed_sections
        section.LE_point .-= origin_offset
        section.TE_point .-= origin_offset
    end
    for section in vsm_wing.unrefined_sections
        section.LE_point .-= origin_offset
        section.TE_point .-= origin_offset
    end
end

"""
    apply_aero_z_offset!(vsm_wing, aero_z_offset)

Apply z-axis offset to VSM panel positions in body frame.

For QUATERNION wings, this shifts the aerodynamic center of pressure
in the positive z-direction (body frame) to adjust the moment arm.
This is applied AFTER the COM adjustment.

# Arguments
- `vsm_wing`: VortexStepMethod.Wing with sections to adjust
- `aero_z_offset`: Distance to shift panels in +z direction [m]
"""
function apply_aero_z_offset!(vsm_wing, aero_z_offset)
    if abs(aero_z_offset) > 1e-10
        offset_vec = [0.0, 0.0, aero_z_offset]
        for section in vsm_wing.refined_sections
            section.LE_point .+= offset_vec
            section.TE_point .+= offset_vec
        end
        for section in vsm_wing.non_deformed_sections
            section.LE_point .+= offset_vec
            section.TE_point .+= offset_vec
        end
        for section in vsm_wing.unrefined_sections
            section.LE_point .+= offset_vec
            section.TE_point .+= offset_vec
        end
    end
end

"""
    VortexStepMethod.Wing(set::Settings; prn=true, kwargs...)

Create a `Wing` geometry object from the settings provided.

This constructor checks for .obj and .dat files in the model directory.
If found, it uses `VortexStepMethod.ObjWing(obj_path, dat_path)` to load the wing.
Otherwise, it falls back to loading from `aero_geometry.yaml`.

This is a constructor helper that reads geometry from the `Settings` object
and initializes the `Wing` object from `VortexStepMethod.jl`.
"""
function VortexStepMethod.Wing(set::Settings, vsm_set::VortexStepMethod.VSMSettings; prn=true, kwargs...)
    # Check for .obj and .dat files in the model directory
    model_dir = get_data_path()
    obj_path = joinpath(model_dir, set.model)
    dat_path = joinpath(model_dir, set.foil_file)

    if isfile(obj_path) && isfile(dat_path)
        # Use ObjWing constructor (default path)
        prn && @info "Loading wing from .obj/.dat files"

        if set.physical_model == "simple_ram"
            n_unrefined_sections = 2
        else
            n_unrefined_sections = 4
        end

        return VortexStepMethod.ObjWing(obj_path, dat_path;
            mass=set.mass, crease_frac=set.crease_frac, n_unrefined_sections,
            align_to_principal=true, prn, kwargs...
        )
    end

    # Fallback: load from aero_geometry.yaml using provided vsm_set
    prn && @info "Using provided VSMSettings for wing creation"
    return VortexStepMethod.Wing(vsm_set; kwargs...)
end

"""
    calc_pos(wing::Wing, gamma, frac)

Calculate a position on the kite based on spanwise (`gamma`) and chordwise (`frac`) parameters.
"""
function calc_pos(wing::VortexStepMethod.Wing, gamma, frac)
    le_pos = [wing.le_interp[i](gamma) for i in 1:3]
    chord = [wing.te_interp[i](gamma) for i in 1:3] .- le_pos
    pos = le_pos .+ chord .* frac
    return pos
end

"""
    cad_to_body_frame(wing::Wing, pos)

Transform a position from the CAD frame to the wing's body frame.
"""
function cad_to_body_frame(wing::VortexStepMethod.Wing, pos)
    return wing.R_cad_body * (pos + wing.T_cad_body)
end
