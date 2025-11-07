# Copyright (c) 2025 Bart van de Lint
# SPDX-License-Identifier: MPL-2.0

"""
    VortexStepMethod.Wing(set::Settings; prn=true, kwargs...)

Create a `Wing` geometry object from the settings provided.

This constructor checks for .obj and .dat files in the model directory.
If found, it uses `VortexStepMethod.ObjWing(obj_path, dat_path)` to load the wing.
Otherwise, it falls back to loading from `aero_geometry.yaml`.

This is a constructor helper that reads geometry from the `Settings` object
and initializes the `Wing` object from `VortexStepMethod.jl`.
"""
function VortexStepMethod.Wing(set::Settings; prn=true, kwargs...)
    # Check for .obj and .dat files in the model directory
    model_dir = get_data_path()
    obj_path = joinpath(model_dir, set.model)
    dat_path = joinpath(model_dir, set.foil_file)

    if isfile(obj_path) && isfile(dat_path)
        # Use ObjWing constructor (default path)
        prn && @info "Loading wing from .obj/.dat files"

        if set.physical_model == "simple_ram"
            n_groups=2
        else
            n_groups=4
        end

        return VortexStepMethod.ObjWing(obj_path, dat_path;
            mass=set.mass, crease_frac=set.crease_frac, n_groups,
            align_to_principal=true, prn, kwargs...
        )
    end

    # Fallback: load from aero_geometry.yaml
    vsm_set_path = joinpath(model_dir, "vsm_settings.yaml")
    prn && @info "No .obj/.dat files found, using aero_geometry.yaml: $vsm_set_path"
    return VortexStepMethod.Wing(VSMSettings(vsm_set_path; data_prefix=false);
                                 kwargs...)
end

"""
    SegmentType `POWER_LINE` `STEERING_LINE` `BRIDLE`

Enumeration for the type of a tether segment.

# Elements
- `POWER_LINE`: A segment belonging to a main power line.
- `STEERING_LINE`: A segment belonging to a steering line.
- `BRIDLE`: A segment belonging to the bridle system.
"""
@enum SegmentType begin
    POWER_LINE
    STEERING_LINE
    BRIDLE
end

"""
    DynamicsType `DYNAMIC` `QUASI_STATIC` `WING` `STATIC`

Enumeration for the dynamic model governing a point's motion.

# Elements
- `DYNAMIC`: The point is a dynamic point mass, moving according to Newton's second law.
- `QUASI_STATIC`: The point's acceleration is constrained to zero, representing a force equilibrium.
- `WING`: The point is rigidly attached to a wing body and moves with it.
- `STATIC`: The point's position is fixed in the world frame.
"""
@enum DynamicsType begin
    DYNAMIC
    QUASI_STATIC
    WING
    STATIC
end

"""
    WingType `QUATERNION` `REFINE`

Enumeration for the aerodynamic model type of a wing.

# Elements
- `QUATERNION`: Wing uses quaternion-based rigid body dynamics with twist groups.
  Aerodynamic forces/moments are applied to the wing center of mass.
- `REFINE`: Wing uses refined per-panel forces directly applied to structural points.
  VSM panel forces are lumped to WING-type points with no rigid body constraint.
"""
@enum WingType begin
    QUATERNION
    REFINE
end

"""
    mutable struct Point

A point mass, representing a node in the mass-spring system.

$(TYPEDFIELDS)
"""
mutable struct Point
    const idx::Int16
    const transform_idx::Int16 # idx of wing used for initial orientation
    const wing_idx::Int16
    const pos_cad::KVec3
    const pos_b::KVec3 # pos relative to wing COM in body frame
    const pos_w::KVec3 # pos in world frame
    const vel_w::KVec3 # vel in world frame
    const disturb::KVec3 # disturbing force
    const force::KVec3
    const aero_force::KVec3 # aerodynamic force in body frame (for REFINE WING points)
    const type::DynamicsType
    mass::SimFloat
    body_frame_damping::SimFloat
    world_frame_damping::SimFloat
    fix_sphere::Bool
    fix_static::Bool
end

"""
    Point(idx, pos_cad, type; wing_idx=1, vel_w=zeros(KVec3), transform_idx=1, mass=0.0)

Constructs a `Point` object, which can be of four different [`DynamicsType`](@ref)s:
- `STATIC`: The point does not move. ``\\ddot{\\mathbf{r}} = \\mathbf{0}``
- `DYNAMIC`: The point moves according to Newton's second law. ``\\ddot{\\mathbf{r}} = \\mathbf{F}/m``
- `QUASI_STATIC`: The acceleration is constrained to be zero by solving a nonlinear problem. ``\\mathbf{F}/m = \\mathbf{0}``
- `WING`: The point has a static position in the rigid body wing frame. ``\\mathbf{r}_w = \\mathbf{r}_{wing} + \\mathbf{R}_{b\\rightarrow w} \\mathbf{r}_b``

where:
- ``\\mathbf{r}`` is the point position vector
- ``\\mathbf{F}`` is the net force acting on the point
- ``m`` is the point mass
- ``\\mathbf{r}_w`` is the position in world frame
- ``\\mathbf{r}_{wing}`` is the wing center position
- ``\\mathbf{R}_{b\\rightarrow w}`` is the rotation matrix from body to world frame
- ``\\mathbf{r}_b`` is the position in body frame

# Arguments
- `idx::Int16`: Unique identifier for the point.
- `pos_cad::KVec3`: Position of the point in the CAD frame.
- `type::DynamicsType`: Dynamics type of the point (`STATIC`, `DYNAMIC`, etc.).

# Keyword Arguments
- `wing_idx::Int16=1`: Index of the wing this point is attached to.
- `vel_w::KVec3=zeros(KVec3)`: Initial velocity of the point in world frame.
- `transform_idx::Int16=1`: Index of the transform used for initial positioning.
- `mass::Float64=0.0`: Mass of the point [kg].
- `body_frame_damping::Float64=0.0`: Damping coefficient for bridle points.
- `world_frame_damping::Float64=0.0`: Damping coefficient for world frame damping.
- `fix_sphere::Bool=false`: If true, constrains the point to a sphere.
- `fix_static::Bool=false`: If true, dynamically freezes the point (behaves like STATIC).

# Returns
- `Point`: A new `Point` object.
"""
function Point(idx, pos_cad, type;
    wing_idx=1, vel_w=zeros(KVec3), transform_idx=1,
    mass=0.0, body_frame_damping=0.0, world_frame_damping=0.0,
    fix_sphere=false, fix_static=false
)
    Point(idx, transform_idx, wing_idx, pos_cad, zeros(KVec3), zeros(KVec3),
        vel_w, zeros(KVec3), zeros(KVec3), zeros(KVec3), type, mass,
        body_frame_damping, world_frame_damping, fix_sphere, fix_static)
end

"""
    mutable struct Group

A set of bridle lines that share the same twist angle and trailing edge angle.

$(TYPEDFIELDS)
"""
mutable struct Group
    const idx::Int16
    const point_idxs::Vector{Int16}
    const gamma::SimFloat  # Spanwise parameter (-1 to 1)
    le_pos::KVec3  # Leading edge position
    chord::KVec3   # Chord vector in body frame
    y_airf::KVec3  # Spanwise vector in local panel frame
    const type::DynamicsType
    moment_frac::SimFloat
    damping::SimFloat
    twist::SimFloat
    twist_ω::SimFloat
    tether_force::SimFloat
    tether_moment::SimFloat
    aero_moment::SimFloat
end

"""
    Group(idx, point_idxs, gamma, type, moment_frac; damping=50.0)

Constructs a `Group` object representing a collection of points on a
kite body that share a common twist deformation.

A `Group` models the local deformation of a kite wing section through
twist dynamics. All points within a group undergo the same twist
rotation about the chord vector.

The group geometry (le_pos, chord, y_airf) is calculated later in the
SystemStructure constructor once the VSM wing is available.

# Arguments
- `idx::Int16`: Unique identifier for the group.
- `point_idxs::Vector{Int16}`: Indices of points that move together.
- `gamma`: Spanwise parameter (-1 to 1) along the wing.
- `type::DynamicsType`: DYNAMIC or QUASI_STATIC.
- `moment_frac::SimFloat`: Chordwise rotation point (0=LE, 1=TE).

# Keyword Arguments
- `damping::SimFloat=50.0`: Damping coefficient for twist dynamics.

# Returns
- `Group`: A new `Group` object (geometry set to zeros initially).
"""
function Group(idx, point_idxs, gamma, type, moment_frac;
               damping=50.0)
    Group(idx, point_idxs, gamma,
          zeros(KVec3), zeros(KVec3), zeros(KVec3),
          type, moment_frac, damping,
          0.0, 0.0, 0.0, 0.0, 0.0)
end

"""
    Group(idx, point_idxs, vsm_wing::Wing, gamma,
          type, moment_frac; damping=50.0)

Legacy constructor that calculates geometry from vsm_wing directly.
Kept for backward compatibility with predefined structures.
"""
function Group(idx, point_idxs, vsm_wing::Wing, gamma,
               type, moment_frac; damping=50.0)
    le_pos = [vsm_wing.le_interp[i](gamma) for i in 1:3]
    chord = [vsm_wing.te_interp[i](gamma) for i in 1:3] .- le_pos
    y_airf = normalize([vsm_wing.le_interp[i](gamma-0.01)
        for i in 1:3] - le_pos)
    Group(idx, point_idxs, gamma, le_pos, chord, y_airf,
          type, moment_frac, damping,
          0.0, 0.0, 0.0, 0.0, 0.0)
end

"""
    mutable struct Segment

A segment representing a spring-damper connection from one point to another.

$(TYPEDFIELDS)
"""
mutable struct Segment
    const idx::Int16
    const point_idxs::Tuple{Int16, Int16}
    axial_stiffness::SimFloat
    axial_damping::SimFloat
    l0::SimFloat
    compression_frac::SimFloat
    diameter::SimFloat
    len::SimFloat # current len of the segment
    force::SimFloat # current force in the segment
end

"""
    Segment(idx, point_idxs, axial_stiffness, axial_damping, diameter; l0, compression_frac)

Inner constructor for a `Segment` object. See [`Segment`](@ref) for details.
"""
function Segment(idx, point_idxs, axial_stiffness, axial_damping, diameter; 
    l0=zero(SimFloat), compression_frac=0.1
)
    Segment(idx, point_idxs, axial_stiffness, axial_damping, l0, compression_frac, 
        diameter, zero(SimFloat), zero(SimFloat))
end

"""
    Segment(idx, set, point_idxs, type; l0, compression_frac, axial_stiffness, axial_damping)

Constructs a `Segment` object representing an elastic spring-damper connection between two points.

The segment follows Hooke's law with damping and aerodynamic drag:

**Spring-Damper Force:**
```math
\\mathbf{F}_{spring} = \\left[k(l - l_0) - c\\dot{l}\\right]\\hat{\\mathbf{u}}
```

**Aerodynamic Drag:**
```math
\\mathbf{F}_{drag} = \\frac{1}{2}\\rho C_d A |\\mathbf{v}_a| \\mathbf{v}_{a,\\perp}
```

**Total Force:**
```math
\\mathbf{F}_{total} = \\mathbf{F}_{spring} + \\mathbf{F}_{drag}
```

where:
- ``k = \\frac{E \\pi d^2/4}{l}`` is the axial stiffness
- ``l`` is current length, ``l_0`` is unstretched length
- ``c = \\frac{\\xi}{c_{spring}} k`` is damping coefficient
- ``\\hat{\\mathbf{u}} = \\frac{\\mathbf{r}_2 - \\mathbf{r}_1}{l}`` is unit vector along segment
- ``\\dot{l} = (\\mathbf{v}_1 - \\mathbf{v}_2) \\cdot \\hat{\\mathbf{u}}`` is extension rate
- ``\\mathbf{v}_{a,\\perp}`` is apparent wind velocity perpendicular to segment

# Arguments
- `idx::Int16`: Unique identifier for the segment.
- `set::Settings`: The settings object containing material properties.
- `point_idxs::Tuple{Int16, Int16}`: Tuple containing the indices of the two points.
- `type::SegmentType`: Type of the segment (`POWER_LINE`, `STEERING_LINE`, `BRIDLE`).

# Keyword Arguments
- `l0::SimFloat=zero(SimFloat)`: Unstretched length [m]. Calculated from point positions if zero.
- `compression_frac::SimFloat=0.0`: Stiffness reduction factor in compression.
- `diameter_mm::Float64=NaN`: Tether diameter [mm]. If `NaN`, uses default from settings.
- `axial_stiffness::Float64=NaN`: Axial stiffness [N]. If `NaN`, it's calculated from diameter and material properties.
- `axial_damping::Float64=NaN`: Axial damping [Ns]. If `NaN`, it's calculated from settings.

# Returns
- `Segment`: A new `Segment` object.
"""
function Segment(idx, set, point_idxs, type;
    l0=zero(SimFloat), compression_frac=0.0, diameter_mm=NaN, axial_stiffness=NaN, axial_damping=NaN
)
    # Set default diameter from settings if not specified
    if isnan(diameter_mm)
        (type == BRIDLE) && (diameter_mm = set.bridle_tether_diameter)
        (type == POWER_LINE) && (diameter_mm = set.power_tether_diameter)
        (type == STEERING_LINE) && (diameter_mm = set.steering_tether_diameter)
    end
    # Convert diameter from mm to m
    diameter_m = 0.001 * diameter_mm

    # Compute axial_stiffness if not provided
    if isnan(axial_stiffness)
        axial_stiffness = set.e_tether * (diameter_m/2)^2 * π
        if type == BRIDLE
            stiffness_frac = 0.01
        else
            stiffness_frac = 1.0
        end
        axial_stiffness *= stiffness_frac
    end

    # Compute axial_damping if not provided
    if isnan(axial_damping)
        # Use rel_damping if available, otherwise compute from axial_damping/axial_stiffness ratio
        if hasproperty(set, :rel_damping) && set.rel_damping != 0.0
            axial_damping = set.rel_damping * axial_stiffness
        elseif hasproperty(set, :axial_damping) && hasproperty(set, :axial_stiffness) &&
                set.axial_damping != 0.0
            axial_damping = (set.axial_damping / set.axial_stiffness) * axial_stiffness
        else
            @warn "Axial damping is zero!"
            axial_damping = 0.0  # fallback if no damping info available
        end
    end

    Segment(idx, point_idxs, axial_stiffness, axial_damping, l0, compression_frac, 
        diameter_m, zero(SimFloat), zero(SimFloat))
end

"""
    mutable struct Pulley

A pulley described by two segments with the common point of the segments being the pulley.

$(TYPEDFIELDS)
"""
mutable struct Pulley
    const idx::Int16
    const segment_idxs::Tuple{Int16, Int16}
    const type::DynamicsType
    sum_len::SimFloat
    len::SimFloat
    vel::SimFloat
end

"""
    Pulley(idx, segment_idxs, type)

Constructs a `Pulley` object that enforces length redistribution between two segments.

The pulley constraint maintains constant total length while allowing force transmission:

**Constraint Equations:**
```math
l_1 + l_2 = l_{total} = \\text{constant}
```

**Force Balance:**
```math
F_{pulley} = F_1 - F_2
```

**Dynamics:**
```math
m\\ddot{l}_1 = F_{pulley} = F_1 - F_2
```

where:
- ``l_1, l_2`` are the lengths of connected segments
- ``F_1, F_2`` are the spring forces in the segments
- ``m = \\rho_{tether} \\pi (d/2)^2 l_{total}`` is the total mass of both segments
- ``\\dot{l}_1 + \\dot{l}_2 = 0`` (velocity constraint)

The pulley can have two [`DynamicsType`](@ref)s:
- `DYNAMIC`: the length redistribution follows Newton's second law: ``m\\ddot{l}_1 = F_1 - F_2``
- `QUASI_STATIC`: the forces are balanced instantaneously: ``F_1 = F_2``

# Arguments
- `idx::Int16`: Unique identifier for the pulley.
- `segment_idxs::Tuple{Int16, Int16}`: Tuple containing the indices of the two segments.
- `type::DynamicsType`: Dynamics type of the pulley (`DYNAMIC` or `QUASI_STATIC`).

# Returns
- `Pulley`: A new `Pulley` object.
"""
function Pulley(idx, segment_idxs, type)
    return Pulley(idx, segment_idxs, type, 0.0, 0.0, 0.0)
end

"""
    mutable struct Tether

A collection of segments that are controlled together by a winch.

$(TYPEDFIELDS)
"""
mutable struct Tether
    const idx::Int16
    const segment_idxs::Vector{Int16}
    const winch_idx::Int16
    stretched_len::SimFloat
end

"""
    Tether(idx, segment_idxs, winch_idx)

Constructs a `Tether` object representing a flexible line composed of multiple segments.

A tether enforces a shared unstretched length constraint across all its constituent segments:

**Length Constraint:**
```math
\\sum_{i \\in \\text{segments}} l_{0,i} = L
```

**Winch Control:**
The unstretched tether length `L` is controlled by a winch.

# Arguments
- `idx::Int16`: Unique identifier for the tether.
- `segment_idxs::Vector{Int16}`: Indices of segments that form this tether.
- `winch_idx::Int16`: Index of the winch controlling this tether.

# Returns
- `Tether`: A new `Tether` object.
"""
function Tether(idx, segment_idxs, winch_idx)
    return Tether(idx, segment_idxs, winch_idx, 0.0)
end

"""
    mutable struct Winch

A set of tethers (or a single tether) connected to a winch mechanism.

$(TYPEDFIELDS)
"""
mutable struct Winch
    const idx::Int16
    const tether_idxs::Vector{Int16}
    tether_len::Union{SimFloat, Nothing}
    tether_vel::SimFloat
    tether_acc::SimFloat
    set_value::SimFloat
    brake::Bool
    const force::KVec3
    gear_ratio::SimFloat
    drum_radius::SimFloat
    f_coulomb::SimFloat
    c_vf::SimFloat
    inertia_total::SimFloat
    friction::SimFloat
end

"""
    Winch(idx, set, tether_idxs; tether_len=0.0, tether_vel=0.0, brake=false)

Constructs a `Winch` object that controls tether length through torque or speed regulation.

The winch acceleration function `α` depends on the winch model type:
- **Torque-controlled**: Direct torque input with motor dynamics.
- **Speed-controlled**: Velocity regulation with internal control loops.

For detailed mathematical models of winch dynamics, motor characteristics, and control algorithms,
see the [WinchModels.jl documentation](https://github.com/aenarete/WinchModels.jl/blob/main/docs/winch.md).

# Arguments
- `idx::Int16`: Unique identifier for the winch.
- `set::Settings`: The main settings object, used to retrieve winch parameters.
- `tether_idxs::Vector{Int16}`: Vector of indices of the tethers connected to this winch.

# Keyword Arguments
- `tether_len::SimFloat=0.0`: Initial tether length [m].
- `tether_vel::SimFloat=0.0`: Initial tether velocity (reel-out rate) [m/s].
- `brake::Bool=false`: If true, the winch brake is engaged.

# Returns
- `Winch`: A new `Winch` object.
"""
function Winch(idx, set::Settings, tether_idxs; tether_len=0.0, tether_vel=0.0, brake=false)
    return Winch(idx, tether_idxs, tether_len, tether_vel, 0.0, 0.0, brake, zeros(KVec3),
                 set.gear_ratio, set.drum_radius, set.f_coulomb, set.c_vf,
                 set.inertia_total, zero(SimFloat))
end

"""
    Winch(idx, tether_idxs, gear_ratio, drum_radius, f_coulomb, c_vf, inertia_total; tether_len=0.0, tether_vel=0.0, brake=false)

Constructs a `Winch` object by directly providing its physical parameters.

This constructor is an alternative to creating a winch from a `Settings` object,
allowing for more modular or programmatic creation of winch components.

# Arguments
- `idx::Int16`: Unique identifier for the winch.
- `tether_idxs::Vector{Int16}`: Vector of indices of the tethers connected to this winch.
- `gear_ratio::SimFloat`: The gear ratio of the winch.
- `drum_radius::SimFloat`: The radius of the winch drum [m].
- `f_coulomb::SimFloat`: Coulomb friction force [N].
- `c_vf::SimFloat`: Viscous friction coefficient [Ns/m].
- `inertia_total::SimFloat`: Total inertia of the motor, gearbox, and drum [kgm²].

# Keyword Arguments
- `tether_len::SimFloat=0.0`: Initial tether length [m].
- `tether_vel::SimFloat=0.0`: Initial tether velocity (reel-out rate) [m/s].
- `brake::Bool=false`: If true, the winch brake is engaged.

# Returns
- `Winch`: A new `Winch` object.
"""
function Winch(idx, tether_idxs, gear_ratio, drum_radius, f_coulomb, c_vf, inertia_total;
               tether_len=0.0, tether_vel=0.0, brake=false)
    return Winch(idx, tether_idxs, tether_len, tether_vel, 0.0, 0.0, brake, zeros(KVec3),
                 gear_ratio, drum_radius, f_coulomb, c_vf, inertia_total, zero(SimFloat))
end

"""
    abstract type AbstractWing

Abstract base type for all wing implementations.

Concrete subtypes must implement rigid body dynamics and provide a reference frame
for attached points and groups.
"""
abstract type AbstractWing end

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
    const idx::Int16

    # Structural information
    const group_idxs::Vector{Int16}
    const transform_idx::Int16
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
    const vsm_y::Vector{SimFloat}
    const vsm_x::Vector{SimFloat}
    const vsm_jac::Matrix{SimFloat}

    # REFINE-specific fields (Nothing for QUATERNION wings)
    point_to_vsm_point::Union{Nothing, Dict{Int16, Tuple{Int16, Symbol}}}
    wing_segments::Union{Nothing,
        Vector{Tuple{Int16, Int16}}}

    # Orientation reference points for REFINE wings
    # (Nothing for QUATERNION wings)
    # Used to calculate R_b_w from structural deformation
    # Can specify single point or vector of points to average:
    #   (12, 13) - point 12 to point 13
    #   (12, [13, 14]) - point 12 to average of points 13,14
    #   ([11, 12], [13, 14]) - average of 11,12 to average of 13,14
    # Z-axis: Normal to wing plane, Y-axis: Spanwise, X = Y × Z (chord)
    z_ref_points::Union{Nothing, Tuple{Union{Int16, Vector{Int16}}, Union{Int16, Vector{Int16}}}}
    y_ref_points::Union{Nothing, Tuple{Union{Int16, Vector{Int16}}, Union{Int16, Vector{Int16}}}}

    # KCU origin point for REFINE wings
    # (Nothing for QUATERNION wings)
    # Defines wing.pos_w = pos[:, origin_idx] to track structural deformation
    origin_idx::Union{Nothing, Int16}

    function VSMWing(base::BaseWing, vsm_aero, vsm_wing, vsm_solver, vsm_y, vsm_x, vsm_jac, point_to_vsm_point, wing_segments, z_ref_points, y_ref_points, origin_idx)
        new(base, vsm_aero, vsm_wing, vsm_solver, vsm_y, vsm_x, vsm_jac, point_to_vsm_point, wing_segments, z_ref_points, y_ref_points, origin_idx)
    end
end

# Delegate property access to base wing for VSMWing
function Base.getproperty(wing::VSMWing, sym::Symbol)
    if sym in (:base, :vsm_aero, :vsm_wing, :vsm_solver, :vsm_y, :vsm_x, :vsm_jac, :point_to_vsm_point, :wing_segments, :z_ref_points, :y_ref_points, :origin_idx)
        return getfield(wing, sym)
    else
        return getproperty(getfield(wing, :base), sym)
    end
end

function Base.setproperty!(wing::VSMWing, sym::Symbol, value)
    if sym in (:base, :vsm_aero, :vsm_wing, :vsm_solver, :vsm_y, :vsm_x, :vsm_jac, :point_to_vsm_point, :wing_segments, :z_ref_points, :y_ref_points, :origin_idx)
        setfield!(wing, sym, value)
    else
        setproperty!(getfield(wing, :base), sym, value)
    end
end

"""
    BaseWing(idx::Int16, group_idxs::Vector{Int16}, R_b_c::Matrix{SimFloat},
             pos_cad::KVec3, inertia_principal::KVec3; transform_idx=1, y_damping=150.0, wing_type=QUATERNION)

Constructs a `BaseWing` object representing a rigid body reference frame.

# Arguments
- `idx::Int16`: Unique identifier for the wing.
- `group_idxs::Vector{Int16}`: Indices of groups attached to this wing.
- `R_b_c::Matrix{SimFloat}`: Rotation matrix from body frame to CAD frame.
- `pos_cad::KVec3`: Position of wing center of mass in CAD frame.
- `inertia_principal::KVec3`: Principal moments of inertia [Ixx, Iyy, Izz] in body frame.

# Keyword Arguments
- `transform_idx::Int16=1`: Transform used for initial positioning and orientation.
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
    VSMWing(idx::Int16, set::Settings, group_idxs::Vector{Int16},
            R_b_c::Matrix{SimFloat}, pos_cad::KVec3;
            transform_idx=1, y_damping=150.0,
            wing_type=QUATERNION, point_to_vsm_point=nothing,
            wing_segments=nothing, x_ref_points=nothing,
            y_ref_points=nothing)

Constructs a `VSMWing` object with Vortex Step Method aerodynamics.
Creates vsm_wing, vsm_aero, and vsm_solver internally.

# Arguments
- `idx::Int16`: Unique identifier for the wing.
- `set::Settings`: Settings object for VSM configuration.
- `group_idxs::Vector{Int16}`: Indices of groups (QUATERNION only).
- `R_b_c::Matrix{SimFloat}`: Rotation matrix body→CAD.
- `pos_cad::KVec3`: Position of wing COM in CAD frame.

# Keyword Arguments
- `transform_idx::Int16=1`: Transform for initial positioning.
- `y_damping::SimFloat=150.0`: Lateral damping coefficient.
- `wing_type::WingType=QUATERNION`: Aerodynamic model type.
- `point_to_vsm_point`: 1:1 structural point to VSM point mapping (REFINE only).
- `wing_segments`: LE/TE pairs (REFINE only).
- `x_ref_points`: Chord direction reference (REFINE only).
- `y_ref_points`: Span direction reference (REFINE only).

# Returns
- `VSMWing`: A new VSM wing object.
"""
function VSMWing(idx::Int, set::Settings,
                 group_idxs::AbstractVector;
                 R_b_c::Union{Nothing,AbstractMatrix}=nothing,
                 pos_cad::Union{Nothing,AbstractVector}=nothing,
                 transform_idx=1, y_damping=150.0,
                 inertia_diag=nothing,
                 wing_type::WingType=QUATERNION,
                 point_to_vsm_point::Union{Nothing, Dict{Int16, Tuple{Int16, Symbol}}}=nothing,
                 wing_segments::Union{Nothing,
                     Vector{Tuple{Int16, Int16}}}=nothing,
                 z_ref_points::Union{Nothing,
                     Tuple{Union{Int16, Vector{Int16}}, Union{Int16, Vector{Int16}}}}=nothing,
                 y_ref_points::Union{Nothing,
                     Tuple{Union{Int16, Vector{Int16}}, Union{Int16, Vector{Int16}}}}=nothing,
                 origin_idx::Union{Nothing, Int16}=nothing)

    # Set defaults from VortexStepMethod if not provided
    if isnothing(R_b_c) || isnothing(pos_cad)
        temp_wing = VortexStepMethod.Wing(set; prn=false)
        isnothing(R_b_c) && (R_b_c = temp_wing.R_cad_body)
        isnothing(pos_cad) &&
            (pos_cad = temp_wing.T_cad_body)
    end

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
    vsm_wing = VortexStepMethod.Wing(set; prn=false)
    vsm_aero = VortexStepMethod.BodyAerodynamics([vsm_wing])
    vsm_solver = VortexStepMethod.Solver(vsm_aero;
        solver_type=VortexStepMethod.NONLIN,
        atol=2e-8, rtol=2e-8)

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
        ny = length(group_idxs) + 3 + 3
        nx = length(group_idxs) + 3 + 3
    end

    return VSMWing(base, vsm_aero, vsm_wing, vsm_solver,
                   zeros(SimFloat, ny), zeros(SimFloat, nx),
                   zeros(SimFloat, nx, ny),
                   point_to_vsm_point, wing_segments,
                   z_ref_points, y_ref_points, origin_idx)
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
    ny = length(group_idxs) + 3 + 3
    nx = length(group_idxs) + 3 + 3
    return VSMWing(base, vsm_aero, vsm_wing, vsm_solver,
        zeros(SimFloat, ny), zeros(SimFloat, nx),
        zeros(SimFloat, nx, ny),
        nothing, nothing, nothing, nothing, nothing)
end

"""
    Wing(idx, vsm_aero, vsm_wing, vsm_solver, group_idxs, R_b_c, pos_cad; transform_idx)

Constructs a `VSMWing` object (backward compatibility constructor).

This is a convenience constructor that creates a VSMWing for backward compatibility
with existing code. New code should use `VSMWing(...)` directly.

# Arguments
- `idx::Int16`: Unique identifier for the wing.
- `vsm_aero`, `vsm_wing`, `vsm_solver`: Vortex Step Method components.
- `group_idxs::Vector{Int16}`: Indices of groups attached to this wing.
- `R_b_c::Matrix{SimFloat}`: Rotation matrix from body frame to CAD frame.
- `pos_cad::KVec3`: Position of wing center of mass in CAD frame.

# Keyword Arguments
- `transform_idx::Int16=1`: Transform used for initial positioning and orientation.
- `y_damping::SimFloat=150.0`: Damping coefficient for lateral motion.

# Returns
- `VSMWing`: A new VSM wing object.
"""
function SymbolicAWEModels.Wing(idx, vsm_aero, vsm_wing, vsm_solver, group_idxs, R_b_c,
                                pos_cad; kwargs...)
    return VSMWing(idx, vsm_aero, vsm_wing, vsm_solver, group_idxs, R_b_c, pos_cad; kwargs...)
end

function wing_inertia_principal(vsm_wing)
    if hasproperty(vsm_wing, :inertia_tensor) && size(vsm_wing.inertia_tensor) == (3, 3)
        diag_vals = diag(vsm_wing.inertia_tensor)
        return MVector{3, SimFloat}(diag_vals)
    end
    return MVector{3, SimFloat}(ones(SimFloat, 3))
end

"""
    mutable struct Transform

Describes the spatial transformation (position and orientation) of system components
relative to a base reference point.

$(TYPEDFIELDS)
"""
mutable struct Transform
    const idx::Int16
    const wing_idx::Union{Int16, Nothing}
    const rot_point_idx::Union{Int16, Nothing}
    const base_point_idx::Union{Int16, Nothing}
    const base_transform_idx::Union{Int16, Nothing}
    elevation::SimFloat # The elevation of the rotating point or kite as seen from the base point
    azimuth::SimFloat # The azimuth of the rotating point or kite as seen from the base point
    heading::SimFloat
    base_pos::Union{KVec3, Nothing}
end

"""
    Transform(idx, elevation, azimuth, heading; base_point_idx, base_pos, base_transform_idx, wing_idx, rot_point_idx)

Constructs a `Transform` object that orients system components using spherical coordinates.

All points and wings with a matching `transform_idx` are transformed together as a rigid body:
1. **Translation**: Translate such that `base_point_idx` is at the specified `base_pos`.
2. **Rotation 1**: Rotate so the target (`wing_idx` or `rot_point_idx`) is at (`elevation`, `azimuth`) relative to the base.
3. **Rotation 2**: Rotate all components by `heading` around the base-target vector.

```math
\\mathbf{r}_{transformed} = \\mathbf{r}_{base} + \\mathbf{R}_{heading} \\circ \\mathbf{R}_{elevation,azimuth}(\\mathbf{r} - \\mathbf{r}_{base})
```

# Arguments
- `idx::Int16`: Unique identifier for the transform.
- `elevation::SimFloat`: Target elevation angle from base [rad].
- `azimuth::SimFloat`: Target azimuth angle from base [rad].
- `heading::SimFloat`: Rotation around base-target vector [rad].

# Keyword Arguments
**Base Reference (choose one method):**
- `base_pos` & `base_point_idx`: Use a fixed position and a reference point.
- `base_transform_idx`: Chain to another transform's position.

**Target Object (choose one):**
- `wing_idx`: The wing to position at (`elevation`, `azimuth`).
- `rot_point_idx`: The point to position at (`elevation`, `azimuth`).

# Returns
- `Transform`: A transform affecting all components with a matching `transform_idx`.
"""
function Transform(idx, elevation, azimuth, heading;
        base_point_idx=nothing, base_pos=nothing, base_transform_idx=nothing,
        wing_idx=nothing, rot_point_idx=nothing)
    (isnothing(wing_idx) == isnothing(rot_point_idx)) && error("Either provide a wing_idx or a rot_point_idx, not both or none.")
    (isnothing(base_pos) == isnothing(base_transform_idx)) && error("Either provide the base_pos or the base_transform_idx, not both or none.")
    (isnothing(base_pos) !== isnothing(base_point_idx)) && error("When providing a base_pos, also provide a base_point_idx.")
    Transform(idx, wing_idx, rot_point_idx, base_point_idx, base_transform_idx, elevation, azimuth, heading, base_pos)
end

"""
    Transform(idx, set, base_point_idx; kwargs...)

Constructor helper to create a `Transform` from a `Settings` object.
"""
function Transform(idx, set, base_point_idx; kwargs...)
    Transform(idx, set.elevations[idx], set.azimuths[idx], set.headings[idx], base_point_idx; kwargs...)
end

"""
    get_rot_pos(transform::Transform, wings, points)

Get the position of the rotating object (wing or point) for a given transform.
"""
function get_rot_pos(transform::Transform, wings, points)
    if !isnothing(transform.wing_idx)
        return wings[transform.wing_idx].pos_w
    elseif !isnothing(transform.rot_point_idx)
        return points[transform.rot_point_idx].pos_w
    end
end

"""
    get_base_pos(transform::Transform, wings, points)

Get the base position for a given transform, resolving chained transforms if necessary.
"""
function get_base_pos(transform::Transform, wings, points)
    curr_base_pos = points[transform.base_point_idx].pos_cad
    if !isnothing(transform.base_pos)
        return transform.base_pos, curr_base_pos
    elseif !isnothing(transform.base_transform_idx)
        base_transform = transforms[transform.base_transform_idx]
        return get_rot_pos(base_transform, wings, points), curr_base_pos
    end
end

"""
    struct SystemStructure

A discrete mass-spring-damper representation of a kite system.

This struct holds all components of the physical model, including points, segments,
winches, and wings, forming a complete description of the kite system's structure.

# Components
- [`Point`](@ref): Point masses.
- [`Group`](@ref): Collections of points for wing deformation.
- [`Segment`](@ref): Spring-damper elements.
- [`Pulley`](@ref): Elements that redistribute line lengths.
- [`Tether`](@ref): Groups of segments controlled by a winch.
- [`Winch`](@ref): Ground-based winches.
- [`Wing`](@ref): Rigid wing bodies.
- [`Transform`](@ref): Spatial transformations for initial positioning.
"""
mutable struct SystemStructure
    const name::String
    const set::Settings
    const points::Vector{Point}
    const groups::Vector{Group}
    const segments::Vector{Segment}
    const pulleys::Vector{Pulley}
    const tethers::Vector{Tether}
    const winches::Vector{Winch}
    const wings::Vector{AbstractWing}
    const transforms::Vector{Transform}
    const y::Array{Float64, 2}
    const x::Array{Float64, 2}
    const jac::Array{Float64, 3}
    const wind_vec_gnd::KVec3
    wind_elevation::SimFloat
    stabilize::Bool
    fix_wing::Bool
end

function Base.getproperty(sys::SystemStructure, sym::Symbol)
    if sym == :diff_vars
        vars = SimFloat[]
        # points
        for point in sys.points
            if point.type == DYNAMIC
                append!(vars, point.pos_w)
                append!(vars, point.vel_w)
            end
        end
        # wings
        for wing in sys.wings
            append!(vars, wing.pos_w)
            append!(vars, wing.vel_w)
            append!(vars, wing.Q_b_w)
            append!(vars, wing.ω_b)
        end
        # groups
        for group in sys.groups
            if group.type == DYNAMIC
                push!(vars, group.twist)
                push!(vars, group.twist_ω)
            end
        end
        # pulleys
        for pulley in sys.pulleys
            if pulley.type == DYNAMIC
                push!(vars, pulley.len)
                push!(vars, pulley.vel)
            end
        end
        # winches
        for winch in sys.winches
            push!(vars, winch.tether_len)
            push!(vars, winch.tether_vel)
        end
        return reshape(vars, :, 1) # Return as a column vector (2D array)
    else
        return getfield(sys, sym)
    end
end

function Base.setproperty!(sys::SystemStructure, sym::Symbol, value)
    if sym == :diff_vars
        flat_value = vec(value) # Ensure value is a flat vector
        offset = 1
        # points
        for point in sys.points
            if point.type == DYNAMIC
                point.pos_w .= @view flat_value[offset:offset+2]
                offset += 3
                point.vel_w .= @view flat_value[offset:offset+2]
                offset += 3
            end
        end
        # wings
        for wing in sys.wings
            wing.pos_w .= @view flat_value[offset:offset+2]
            offset += 3
            wing.vel_w .= @view flat_value[offset:offset+2]
            offset += 3
            wing.Q_b_w .= @view flat_value[offset:offset+3]
            offset += 4
            wing.ω_b .= @view flat_value[offset:offset+2]
            offset += 3
        end
        # groups
        for group in sys.groups
            if group.type == DYNAMIC
                group.twist = flat_value[offset]
                offset += 1
                group.twist_ω = flat_value[offset]
                offset += 1
            end
        end
        # pulleys
        for pulley in sys.pulleys
            if pulley.type == DYNAMIC
                pulley.len = flat_value[offset]
                offset += 1
                pulley.vel = flat_value[offset]
                offset += 1
            end
        end
        # winches
        for winch in sys.winches
            winch.tether_len = flat_value[offset]
            offset += 1
            winch.tether_vel = flat_value[offset]
            offset += 1
        end
        return value
    else
        return setfield!(sys, sym, value)
    end
end

"""
    SystemStructure(name, set; points, groups, segments, pulleys, tethers, winches, wings, transforms)

Constructs a `SystemStructure` object representing a complete kite system.

## Physical Models
- **"ram"**: A model with 4 deformable wing groups and a complex pulley bridle system.
- **"simple_ram"**: A model with 4 deformable wing groups and direct bridle connections.

# Arguments
- `name::String`: Model identifier ("ram", "simple_ram", or a custom name).
- `set::Settings`: Configuration parameters from `KiteUtils.jl`.

# Keyword Arguments
- `points`, `groups`, `segments`, etc.: Vectors of the system components.

# Returns
- `SystemStructure`: A complete system ready for building a `SymbolicAWEModel`.
"""
function SystemStructure(name, set;
        points=Point[],
        groups=Group[],
        segments=Segment[],
        pulleys=Pulley[],
        tethers=Tether[],
        winches=Winch[],
        wings=AbstractWing[],
        transforms=Transform[],
        ignore_l0::Bool=false,
    )
    # If no wings defined, convert WING points to STATIC
    if isempty(wings)
        wing_point_idxs = findall(p -> p.type == WING, points)
        if !isempty(wing_point_idxs)
            @warn "No wings provided but " *
                  "$(length(wing_point_idxs)) WING type " *
                  "points found. Converting to STATIC."
            for idx in wing_point_idxs
                points[idx] = Point(
                    points[idx].idx,
                    points[idx].pos_cad,
                    STATIC;
                    mass = points[idx].mass,
                    body_frame_damping =
                        points[idx].body_frame_damping,
                    world_frame_damping =
                        points[idx].world_frame_damping,
                    transform_idx = points[idx].transform_idx
                )
                points[idx].pos_w .= points[idx].pos_cad
                points[idx].vel_w .= 0.0
            end
        end
    end

    for (i, point) in enumerate(points)
        @assert point.idx == i "Point $(point.idx) != $i"
        # Allow transform_idx=0 (no transform) or valid index
        @assert point.transform_idx == 0 ||
                point.transform_idx <= length(transforms)
    end
    for (i, group) in enumerate(groups)
        @assert group.idx == i
    end
    for (i, segment) in enumerate(segments)
        @assert segment.idx == i
        (segment.l0 ≈ 0) && (segment.l0 = norm(points[segment.point_idxs[1]].pos_cad - points[segment.point_idxs[2]].pos_cad))
    end
    for (i, pulley) in enumerate(pulleys)
        @assert pulley.idx == i
    end
    for (i, tether) in enumerate(tethers)
        @assert tether.idx == i
    end
    for (i, winch) in enumerate(winches)
        @assert winch.idx == i
        if iszero(winch.tether_len)
            for segment_idx in tethers[winch.tether_idxs[1]].segment_idxs
                winch.tether_len += segments[segment_idx].l0
            end
        end
    end
    # Initialize group geometries from VSM wing
    for group in groups
        if iszero(group.le_pos)
            # Find which wing this group belongs to
            for wing in wings
                if group.idx in wing.group_idxs
                    vsm_wing = wing.vsm_wing
                    gamma = group.gamma
                    group.le_pos .= [vsm_wing.le_interp[i](gamma)
                        for i in 1:3]
                    te_pos = [vsm_wing.te_interp[i](gamma)
                        for i in 1:3]
                    group.chord .= te_pos .- group.le_pos
                    le_minus = [vsm_wing.le_interp[i](gamma-0.01)
                        for i in 1:3]
                    group.y_airf .= normalize(
                        le_minus - group.le_pos)
                    break
                end
            end
        end
    end

    for (i, wing) in enumerate(wings)
        @assert wing.idx == i
        # For REFINE wings, set defaults if not provided
        if wing.wing_type == REFINE
            # Build point_to_vsm_point mapping if not provided
            if isnothing(wing.point_to_vsm_point)
                # Get WING-type points for this wing
                wing_point_idxs = findall(
                    p -> p.type == WING && p.wing_idx == wing.idx, points)
                wing_points = [points[idx]
                    for idx in wing_point_idxs]
                wing.point_to_vsm_point =
                    build_point_to_vsm_point_mapping(
                        wing_points, wing.vsm_wing)
            end

            wing_point_idxs = collect(keys(
                wing.point_to_vsm_point))
            wing_points = [points[idx]
                for idx in wing_point_idxs]

            # For REFINE wings, pos_cad should be user-specified (KCU position)
            # or default to vsm_wing.T_cad_body (set in VSMWing constructor)
            # DO NOT calculate as centroid - that would misalign VSM panels

            # Identify wing segments (LE/TE pairs)
            if isnothing(wing.wing_segments)
                wing.wing_segments =
                    identify_wing_segments(wing_points)
            end

            # Set default reference points if not provided
            if isnothing(wing.z_ref_points) ||
               isnothing(wing.y_ref_points)
                segs = wing.wing_segments

                if isnothing(wing.z_ref_points)
                    # Use first segment (center LE-TE) for Z (normal)
                    wing.z_ref_points = segs[1]
                end

                if isnothing(wing.y_ref_points)
                    # Use center LE and mid-span LE for Y (spanwise)
                    mid = length(segs) ÷ 2
                    wing.y_ref_points = (segs[1][1],
                        segs[mid][1])
                end
            end

            # Distribute kite mass to WING points for REFINE wings
            if hasproperty(set, :mass) && set.mass > 0
                n_wing_points = length(wing_points)
                mass_per_point = set.mass / n_wing_points
                for point_idx in wing_point_idxs
                    points[point_idx].mass = mass_per_point
                end
            end
        end
    end
    for (i, transform) in enumerate(transforms)
        @assert transform.idx == i
        set.elevations[i] = rad2deg(transform.elevation)
        set.azimuths[i]   = rad2deg(transform.azimuth)
        set.headings[i]   = rad2deg(transform.heading)
    end
    if length(wings) > 0
        ny = 3+length(wings[1].group_idxs)+3
        nx = 3+3+length(wings[1].group_idxs)
    else
        ny = 0
        nx = 0
    end
    y = zeros(length(wings), ny)
    x = zeros(length(wings), nx)
    jac = zeros(length(wings), nx, ny)
    set.physical_model = name
    sys_struct = SystemStructure(name, set, points, groups, segments, pulleys, tethers,
        winches, wings, transforms, y, x, jac, zeros(KVec3), 0.0, false, false)
    reinit!(sys_struct, set)

    # Recalculate segment rest lengths from current positions if requested
    if ignore_l0
        for segment in sys_struct.segments
            p1 = sys_struct.points[segment.point_idxs[1]]
            p2 = sys_struct.points[segment.point_idxs[2]]
            segment.l0 = norm(p2.pos_w - p1.pos_w)
        end
    end

    return sys_struct
end

"""
    apply_heading(vec, R_t_w, curr_R_t_w, heading)

Apply a heading rotation to a vector.
"""
function apply_heading(vec, R_t_w, curr_R_t_w, heading)
    vec_along_z = rotate_around_z(curr_R_t_w' * vec, heading)
    return R_t_w * vec_along_z
end

"""
    reinit!(transforms::Vector{Transform}, sys_struct::SystemStructure)

Apply the initial spatial transformations to all components in a `SystemStructure`.

This function iterates through all transforms and applies the specified translation
and rotation to position and orient the kite system components correctly in the
world frame at the beginning of a simulation.

If transforms is empty, simply initializes pos_w = pos_cad for all components.
"""
function reinit!(transforms::Vector{Transform}, sys_struct::SystemStructure)
    @unpack points, wings = sys_struct
    
    # Handle the case with no transforms: just copy CAD positions to world positions
    if isempty(transforms)
        for point in points
            point.pos_w .= point.pos_cad
            point.vel_w .= 0.0
        end
        for wing in wings
            wing.pos_w .= wing.pos_cad
            wing.vel_w .= 0.0
            wing.ω_b .= 0.0
            wing.Q_b_w .= rotation_matrix_to_quaternion(wing.R_b_c)
        end
        return  # Early return - no transforms to apply
    end
    
    # Apply transforms
    for transform in transforms
        # ==================== TRANSLATE ==================== #
        base_pos, curr_base_pos = get_base_pos(transform, wings, points)
        T = base_pos - curr_base_pos
        for point in points
            if point.transform_idx == transform.idx
                point.pos_w .= point.pos_cad .+ T
                point.vel_w .= 0.0
            end
        end
        for wing in wings
            if wing.transform_idx == transform.idx
                wing.pos_w .= wing.pos_cad .+ T
                wing.vel_w .= 0.0
            end
        end

        # ==================== ROTATE ==================== #
        curr_rot_pos = get_rot_pos(transform, wings, points)
        rel_pos = curr_rot_pos - base_pos

        # Check if wing and base are aligned (avoid division by zero)
        if norm(rel_pos) < 1e-6
            @warn "Transform #$(transform.idx): Wing and base positions are aligned at $(base_pos). " *
                  "Using identity rotation matrix (no rotation applied)."
            curr_R_t_w = Matrix(1.0I, 3, 3)  # Identity matrix
        else
            curr_R_t_w = calc_R_t_w(rel_pos)
        end

        transform_pos = rotate_around_z(rotate_around_y([1,0,0], -transform.elevation),
                                        -transform.azimuth)
        R_t_w = calc_R_t_w(transform_pos)

        for point in points
            if point.transform_idx == transform.idx
                vec = point.pos_w - base_pos
                point.pos_w .= base_pos +
                    apply_heading(vec, R_t_w, curr_R_t_w, transform.heading)
            end
            if point.type == WING
                wing = wings[point.wing_idx]
                # For QUATERNION wings, calculate pos_b from CAD geometry
                # For REFINE wings, pos_b will be calculated after all transforms are complete
                if wing.wing_type != REFINE
                    point.pos_b .= wing.R_b_c' * (point.pos_cad - wing.pos_cad)
                end
            end
        end
        for wing in wings
            if wing.transform_idx == transform.idx
                vec = wing.pos_w - base_pos
                wing.pos_w .= base_pos + apply_heading(vec, R_t_w, curr_R_t_w, transform.heading)
                R_b_w = zeros(3,3)
                for i in 1:3
                    R_b_w[:, i] .= apply_heading(wing.R_b_c[:, i], R_t_w, curr_R_t_w, transform.heading)
                end
                wing.R_b_w = R_b_w
            end
        end
    end

    # Calculate pos_b for REFINE wing points after all transforms are complete
    # This stores the initial body-frame positions for displacement tracking
    for wing in wings
        if wing.wing_type == REFINE
            # Calculate initial R_b_w and origin from structural geometry
            R_b_w, origin = calc_refine_wing_frame(
                points,
                wing.z_ref_points,
                wing.y_ref_points,
                wing.origin_idx
            )

            # Store R_b_w in wing (used for plotting and VSM updates)
            wing.R_b_w = R_b_w

            # Update wing.pos_w to match origin (KCU position)
            wing.pos_w .= origin

            # Calculate and store pos_b for all WING points of this wing
            for point in points
                if point.type == WING && point.wing_idx == wing.idx
                    # Transform to body frame: pos_b = R_b_w' * (pos_w - origin)
                    point.pos_b .= R_b_w' * (point.pos_w - origin)
                end
            end
        end
    end
end

"""
    reposition!(transforms::Vector{Transform}, sys_struct::SystemStructure)

Update the system's spatial orientation based on its current position, preserving velocities.

This function adjusts the orientation of all components in the `SystemStructure` without
altering their dynamic state. Unlike `reinit!`, it uses the current world positions (`pos_w`)
as the starting point for rotations, rather than resetting from the CAD coordinates.

This function is useful for making real-time adjustments to the system's pose during a simulation.
Crucially, it **preserves the existing velocities (`vel_w`) of all points and wings**.

NOTE: the transform.heading is applied relative to the current heading of the system.

# Arguments
- `sys_struct::SystemStructure`: The system model to update.
"""
function reposition!(transforms::Vector{Transform}, sys_struct::SystemStructure)
    @unpack points, wings = sys_struct
    for transform in transforms
        # Get the current positions of the base and the rotating object
        base_pos = points[transform.base_point_idx].pos_w
        rot_pos = get_rot_pos(transform, wings, points)

        # Calculate the current orientation in spherical coordinates
        curr_rel_pos = rot_pos - base_pos

        # Check if wing and base are aligned (avoid division by zero)
        if norm(curr_rel_pos) < 1e-6
            @warn "Transform #$(transform.idx): Wing and base positions are aligned at $(base_pos). " *
                  "Using identity rotation matrix (no rotation applied)."
            curr_R_t_w = Matrix(1.0I, 3, 3)  # Identity matrix
        else
            curr_R_t_w = calc_R_t_w(curr_rel_pos)
        end

        transform_pos = rotate_around_z(rotate_around_y([1,0,0], -transform.elevation), -transform.azimuth)
        R_t_w = calc_R_t_w(transform_pos)

        # Apply the rotation to all relevant points
        for point in points
            if point.transform_idx == transform.idx
                vec = point.pos_w - base_pos
                point.pos_w .= base_pos + apply_heading(vec, R_t_w, curr_R_t_w, transform.heading)
            end
        end

        # Apply the rotation to all relevant wings
        for wing in wings
            if wing.transform_idx == transform.idx
                # Rotate the wing's position
                vec = wing.pos_w - base_pos
                wing.pos_w .= base_pos + apply_heading(vec, R_t_w, curr_R_t_w, transform.heading)

                # Rotate the wing's orientation matrix
                R_b_w = zeros(3,3)
                current_R_b_w = wing.R_b_w
                for i in 1:3
                    R_b_w[:, i] .= apply_heading(current_R_b_w[:, i], R_t_w, curr_R_t_w, transform.heading)
                end
                wing.R_b_w = R_b_w
            end
        end
    end
end

"""
    get_ref_position_from_points(points::Vector{Point}, ref::Int16)
    get_ref_position_from_points(points::Vector{Point}, refs::Vector{Int16})

Helper to get position (single point or average of multiple).
Used for REFINE wing reference point calculations.
"""
get_ref_position_from_points(points::Vector{Point}, ref::Int16) = points[ref].pos_w
function get_ref_position_from_points(points::Vector{Point}, refs::Vector{Int16})
    n = length(refs)
    return sum(points[idx].pos_w for idx in refs) / n
end

"""
    calc_refine_wing_frame(points::Vector{Point}, z_ref_points, y_ref_points, origin_idx)

Calculate R_b_w rotation matrix and origin position for a REFINE wing from structural point positions.

This function implements the same R_b_w calculation logic as used in `generate_system.jl`
for REFINE wings, ensuring consistency between initialization (`reinit!`) and simulation
(symbolic equations).

# Algorithm
1. Extract reference point positions (with averaging if vectors provided)
2. Calculate Z-axis (normal to wing): normalized vector from z_p1 to z_p2
3. Calculate X-axis (chord direction): Y_temp × Z, where Y_temp is from y_p1 to y_p2
4. Calculate Y-axis (spanwise): Z × X (ensures orthogonality and right-handed system)
5. Extract origin position from origin_idx point

# Arguments
- `points::Vector{Point}`: All structural points (must have pos_w initialized)
- `z_ref_points::Tuple{Union{Int16, Vector{Int16}}, Union{Int16, Vector{Int16}}}`:
  Reference points defining Z-axis (normal direction)
- `y_ref_points::Tuple{Union{Int16, Vector{Int16}}, Union{Int16, Vector{Int16}}}`:
  Reference points defining Y-axis (spanwise direction)
- `origin_idx::Int16`: Point index defining wing origin (KCU position)

# Returns
- `R_b_w::Matrix{SimFloat}`: 3x3 rotation matrix from body frame to world frame
- `origin::KVec3`: Origin position in world frame
"""
function calc_refine_wing_frame(
    points::Vector{Point},
    z_ref_points::Tuple{Union{Int16, Vector{Int16}}, Union{Int16, Vector{Int16}}},
    y_ref_points::Tuple{Union{Int16, Vector{Int16}}, Union{Int16, Vector{Int16}}},
    origin_idx::Int16
)
    # Extract reference point positions (with averaging if vectors provided)
    z_p1, z_p2 = z_ref_points
    y_p1, y_p2 = y_ref_points

    pos_z1 = get_ref_position_from_points(points, z_p1)
    pos_z2 = get_ref_position_from_points(points, z_p2)
    pos_y1 = get_ref_position_from_points(points, y_p1)
    pos_y2 = get_ref_position_from_points(points, y_p2)

    # Build rotation matrix from structural geometry
    # Z direction (normal to wing, normalized)
    z_axis = normalize(pos_z2 - pos_z1)

    # Y temp direction (not necessarily orthogonal yet)
    y_temp = normalize(pos_y2 - pos_y1)

    # X = Y_temp × Z (chord direction, orthogonal to Z)
    x_axis = normalize(y_temp × z_axis)

    # Y = Z × X (ensure orthogonality and right-handed system)
    y_axis = z_axis × x_axis

    # Construct rotation matrix [x y z]
    R_b_w = hcat(x_axis, y_axis, z_axis)

    # Extract origin position
    origin = points[origin_idx].pos_w

    return R_b_w, origin
end

"""
    calc_pos(wing::Wing, gamma, frac)

Calculate a position on the kite based on spanwise (`gamma`) and chordwise (`frac`) parameters.
"""
function calc_pos(wing::Wing, gamma, frac)
    le_pos = [wing.le_interp[i](gamma) for i in 1:3]
    chord = [wing.te_interp[i](gamma) for i in 1:3] .- le_pos
    pos = le_pos .+ chord .* frac
    return pos
end

"""
    create_tether(tether_idx, set, points, segments, tethers, attach_point, type, dynamics_type; z, axial_stiffness, axial_damping)

Procedurally create a multi-segment tether.

This function builds a tether from a specified number of segments, connecting a given
`attach_point` on the kite to a new anchor point on the ground.
"""
function create_tether(tether_idx, set, points, segments, tethers, attach_point, 
                       type, dynamics_type; z=[0,0,1], axial_stiffness=NaN, 
                       axial_damping=NaN, d_pos=zeros(3))
    winch_pos = find_axis_point(attach_point.pos_cad, set.l_tether, z) .+ d_pos
    dir = winch_pos - attach_point.pos_cad
    segment_idxs = Int16[]
    winch_idx = 0
    for i in 1:set.segments
        frac = i / set.segments
        pos = attach_point.pos_cad + frac * dir
        point_idx = length(points)+1 # last point idx
        segment_idx = length(segments)+1 # last segment idx
        if i == 1
            last_idx = attach_point.idx
        else
            last_idx = point_idx-1
        end
        if i == set.segments
            points = [points; Point(point_idx, pos, STATIC)]
            winch_idx = points[end].idx
        else
            points = [points; Point(point_idx, pos, dynamics_type)]
        end
        segments = [segments; Segment(segment_idx, set, (last_idx, point_idx), type;
                                      axial_stiffness, axial_damping)]
        push!(segment_idxs, segment_idx)
    end
    tethers = [tethers; Tether(tether_idx, segment_idxs, winch_idx)]
    return points, segments, tethers, tethers[end].idx
end

"""
    cad_to_body_frame(wing::Wing, pos)

Transform a position from the CAD frame to the wing's body frame.
"""
function cad_to_body_frame(wing::Wing, pos)
    return wing.R_cad_body * (pos + wing.T_cad_body)
end

"""
    find_axis_point(P, l, v=[0,0,1])

Calculate the coordinates of a point `Q` that lies on a line defined by vector `v`
and is at a distance `l` from a given point `P`.
"""
function find_axis_point(P, l, v=[0,0,1])
    # Compute discriminant
    D = (v ⋅ P)^2 - norm(v)^2 * (norm(P)^2 - l^2)
    D < 0 && error("No real solution: l is too small or parameters invalid")
    # Solve quadratic for t, choose solution for negative direction
    t = (v ⋅ P - √D) / norm(v)^2
    # Compute point Q = t * v
    return [t * v[1], t * v[2], t * v[3]]
end

"""
    validate_sys_struct(sys_struct::SystemStructure)

Validate a `SystemStructure` for common configuration errors.

This function checks for issues that can cause initialization failures or
numerical problems during simulation. It emits warnings for suspicious
configurations and throws assertions for definite errors.

# Validations Performed

## Point Validations
- NaN mass (error)
- Negative mass (warning)
- NaN position (error)

## Wing Validations
- NaN position (error)
- Zero or near-zero principal inertia components on QUATERNION wings (error/warning)
- NaN inertia values (error)
- Empty group list for QUATERNION wings (warning)

## Winch Validations
- Zero or negative inertia_total (error)
- Very small inertia_total (warning)
- NaN inertia_total (error)

## Segment Validations
- Unusual diameter outside (0, 1) m range (warning)
- Non-positive rest length (error)
- Zero or negative axial stiffness (warning)
- Negative axial damping (warning)

## Pulley Validations
- Zero total length constraint (error)

## Group Validations
- Inconsistent moment_frac across groups (error)
"""
function validate_sys_struct(sys_struct::SystemStructure)
    @unpack points, groups, segments, pulleys, wings, winches = sys_struct

    # ==================== POINT VALIDATIONS ==================== #
    for point in points
        # Check for NaN mass
        if isnan(point.mass)
            error("Point #$(point.idx) has NaN mass")
        end

        # Warn about negative mass (physically nonsensical but still works)
        if point.mass < 0
            @warn "Point #$(point.idx) has negative mass $(point.mass) kg. " *
                  "This is physically nonsensical."
        end

        # Check for NaN position
        if any(isnan.(point.pos_w))
            error("Point #$(point.idx) has NaN position: pos_w = $(point.pos_w)")
        end
    end

    # ==================== WING VALIDATIONS ==================== #
    for wing in wings
        # Check for NaN position (applies to all wing types)
        if any(isnan.(wing.pos_w))
            error("Wing #$(wing.idx) has NaN position: pos_w = $(wing.pos_w)")
        end

        if wing.wing_type == QUATERNION
            I_b = wing.inertia_principal

            # Check for NaN inertia
            if any(isnan.(I_b))
                error("Wing #$(wing.idx) has NaN inertia: I_b = $I_b")
            end

            # Check for zero or suspiciously small inertia
            for i in 1:3
                if I_b[i] ≈ 0.0
                    error("Wing #$(wing.idx) has zero inertia component " *
                          "I_b[$i] = $(I_b[i]). " *
                          "All principal inertia components must be non-zero.")
                elseif I_b[i] < 1e-6
                    @warn "Wing #$(wing.idx) has very small inertia component " *
                          "I_b[$i] = $(I_b[i]) kg⋅m², may cause numerical issues"
                end
            end

            # Warn if QUATERNION wing has no groups
            if isempty(wing.group_idxs)
                @warn "Wing #$(wing.idx) (QUATERNION) has no groups"
            end
        end
        # REFINE wings don't use rigid body inertia, skip
    end

    # ==================== WINCH VALIDATIONS ==================== #
    for winch in winches
        # Check for NaN inertia
        if isnan(winch.inertia_total)
            error("Winch #$(winch.idx) has NaN inertia_total")
        end

        # Check for zero or negative inertia
        if winch.inertia_total ≈ 0.0
            error("Winch #$(winch.idx) has zero inertia_total. " *
                  "All winches must have non-zero inertia.")
        elseif winch.inertia_total < 0
            error("Winch #$(winch.idx) has negative inertia_total " *
                  "$(winch.inertia_total) kg⋅m². Inertia must be positive.")
        elseif winch.inertia_total < 1e-6
            @warn "Winch #$(winch.idx) has very small inertia_total " *
                  "$(winch.inertia_total) kg⋅m², may cause numerical issues"
        end
    end

    # ==================== SEGMENT VALIDATIONS ==================== #
    for segment in segments
        # Diameter should be in valid range (warn only, not critical)
        if !(0 < segment.diameter < 1)
            @warn "Segment #$(segment.idx) has unusual diameter " *
                  "$(segment.diameter) m (expected range: 0 to 1 m)"
        end

        # Rest length must be positive (after initialization)
        if segment.l0 > 0 && !(segment.l0 > 0)
            error("Segment #$(segment.idx) has non-positive rest length " *
                  "l0 = $(segment.l0) m. All segment rest lengths must be positive.")
        end

        # Warn about zero or negative stiffness/damping
        if segment.axial_stiffness ≈ 0.0
            @warn "Segment #$(segment.idx) has zero axial stiffness"
        elseif segment.axial_stiffness < 0
            @warn "Segment #$(segment.idx) has negative axial stiffness " *
                  "$(segment.axial_stiffness) N"
        end

        if segment.axial_damping < 0
            @warn "Segment #$(segment.idx) has negative axial damping " *
                  "$(segment.axial_damping) N⋅s"
        end
    end

    # ==================== PULLEY VALIDATIONS ==================== #
    for pulley in pulleys
        if pulley.sum_len ≈ 0
            error("Pulley #$(pulley.idx) has zero total length constraint " *
                  "(sum_len = $(pulley.sum_len) m). " *
                  "Pulley constraints must have non-zero total length.")
        end
    end

    # ==================== GROUP VALIDATIONS ==================== #
    if length(groups) > 0
        first_moment_frac = groups[1].moment_frac
        for group in groups
            if !(group.moment_frac ≈ first_moment_frac)
                error("Group #$(group.idx) has moment_frac = " *
                      "$(group.moment_frac), but all groups must have the " *
                      "same moment_frac (first group has $(first_moment_frac))")
            end
        end
    end

    return nothing
end

"""
    reinit!(sys_struct::SystemStructure, set::Settings; ignore_l0=false, remake_vsm=false)

Re-initialize a `SystemStructure` from a `Settings` object.

This function resets various component states (e.g., winch lengths, group twists,
pulley positions) to their initial values as defined in the `Settings` object. It
is typically called before starting a new simulation run.

Pulley lengths are initialized proportionally based on current segment lengths:
`pulley.len = segment1.len / (segment1.len+segment2.len) * pulley.sum_len`

# Keyword Arguments
- `ignore_l0::Bool=false`: If true, recalculate segment rest lengths from current positions
- `remake_vsm::Bool=false`: If true, recreate VSM wing, aerodynamics, and solver from settings.
  This is useful after modifying aero_geometry.yaml or other VSM-related configuration files.
  For REFINE wings, also rebuilds the point_to_vsm_point mapping.
"""
function reinit!(sys_struct::SystemStructure, set::Settings; ignore_l0::Bool=false, remake_vsm::Bool=false)
    @unpack points, groups, segments, pulleys, tethers, winches, wings, transforms = sys_struct

    for winch in winches
        winch.tether_len = set.l_tethers[winch.idx]
        winch.tether_vel    = set.v_reel_outs[winch.idx]
    end

    for group in groups
        group.twist = 0.0
        group.twist_ω = 0.0
    end

    for transform in transforms
        transform.elevation = deg2rad(set.elevations[transform.idx])
        transform.azimuth   = deg2rad(set.azimuths[transform.idx])
        transform.heading   = deg2rad(set.headings[transform.idx])
    end

    for segment in segments
        len = norm(points[segment.point_idxs[1]].pos_cad -
                   points[segment.point_idxs[2]].pos_cad)
        (segment.l0 ≈ 0) && (segment.l0 = len)
        segment.len = len
    end

    for pulley in pulleys
        segment1, segment2 = segments[pulley.segment_idxs[1]],
                             segments[pulley.segment_idxs[2]]
        pulley.sum_len = segment1.l0 + segment2.l0

        # Initialize pulley.len proportional to current segment lengths
        # More accurate for asymmetric bridle configurations
        pulley.len = segment1.len / (segment1.len+segment2.len) *
                     pulley.sum_len

        pulley.vel = 0.0
    end

    reinit!(transforms, sys_struct)

    # Recreate VSM wing and aero if requested
    if remake_vsm
        for wing in wings
            # Recreate VSM wing from settings
            wing.vsm_wing = VortexStepMethod.Wing(set; prn=false)
            wing.vsm_aero = VortexStepMethod.BodyAerodynamics([wing.vsm_wing])
            wing.vsm_solver = VortexStepMethod.Solver(wing.vsm_aero;
                solver_type=VortexStepMethod.NONLIN,
                atol=2e-8, rtol=2e-8)

            # For REFINE wings, rebuild point_to_vsm_point mapping
            if wing.wing_type == REFINE && !isnothing(wing.point_to_vsm_point)
                wing_point_idxs = collect(keys(wing.point_to_vsm_point))
                wing_points = [points[idx] for idx in wing_point_idxs]
                wing.point_to_vsm_point =
                    build_point_to_vsm_point_mapping(wing_points, wing.vsm_wing)
            end
        end
    end

    # Calculate ground-level wind vector with direction rotations
    # Matches symbolic equations in generate_system.jl:1259-1264
    upwind_dir = deg2rad(set.upwind_dir)
    wind_elevation = sys_struct.wind_elevation
    wind_scale_gnd = set.v_wind

    # Base wind vector: [0, -1, 0] points upwind
    wind_vec_base = [0.0, -1.0, 0.0]
    # Rotate by elevation around x-axis (vertical tilt)
    wind_vec_elevated = rotate_around_x(wind_vec_base, wind_elevation)
    # Rotate by upwind direction around z-axis (negative for convention)
    wind_vec_rotated = rotate_around_z(wind_vec_elevated, -upwind_dir)
    # Scale by ground wind speed
    wind_vec_gnd = max(wind_scale_gnd, 1e-6) * wind_vec_rotated

    for wing in wings
        # Calculate wind at wing height using atmospheric model
        # Matches symbolic equations in generate_system.jl:1277-1278
        wing_height = max(wing.pos_w[3], 1.0)  # z-component, minimum 1m
        wind_factor = calc_wind_factor(AtmosphericModel(set), wing_height, set)
        wing.v_wind .= wind_factor * wind_vec_gnd

        if wing.wing_type == REFINE
            # Initialize apparent wind in body frame for REFINE wings
            # va_wing = wind_vel - wing_vel + wind_disturb
            # va_b = R_b_w' * va_wing
            # At initialization: wing_vel typically 0, wind_disturb typically 0
            va_wing_w = wing.v_wind - wing.vel_w + wing.wind_disturb
            @show va_wing_w wing.R_b_w
            wing.va_b .= wing.R_b_w' * va_wing_w
        else
            # Initialize vsm_y for QUATERNION wings (REFINE wings have ny=0)
            if length(wing.vsm_y) >= 3
                wing.vsm_y .= 0.0
                wing.vsm_y[1:3] .= wing.R_b_w' * [set.v_wind, 0., 0.]
            end
        end
    end

    # Validate the system structure after initialization
    validate_sys_struct(sys_struct)

    # Recalculate segment rest lengths from current positions if requested
    if ignore_l0
        for segment in segments
            p1 = points[segment.point_idxs[1]]
            p2 = points[segment.point_idxs[2]]
            segment.l0 = norm(p2.pos_w - p1.pos_w)
        end
    end

    return nothing
end

"""
    copy!(sys1::SystemStructure, sys2::SystemStructure)

Copy the dynamic state from one `SystemStructure` (`sys1`) to another (`sys2`).

This function is designed to transfer the state (positions, velocities, etc.) between
two system models, which can have different levels of fidelity. For example, it can
copy the state from a detailed multi-segment tether model (`sys1`) to a simplified
single-segment model (`sys2`).

The function handles several cases:
- If `sys1` and `sys2` have the same structure, it performs a direct copy of all point states.
- If `sys2` is a simplified (1-segment per tether) version of `sys1`, it copies the
  positions and velocities of the tether endpoints.
- It also copies the state of wings, groups, winches, and pulleys where applicable.
"""
function copy!(sys1::SystemStructure, sys2::SystemStructure)
    simple = false

    # copy point pos and vel
    if length(sys1.points) > 0
        if length(sys1.points) == length(sys2.points)
            for (point1, point2) in zip(sys1.points, sys2.points)
                point2.pos_w .= point1.pos_w
                point2.vel_w .= point1.vel_w
                point2.disturb .= point1.disturb
            end
        # if different number of points, copy only the tether points
        elseif length(sys1.tethers) > 1 && length(sys1.tethers) == length(sys2.tethers)
            for (tether1, tether2) in zip(sys1.tethers, sys2.tethers)
                if length(tether1.segment_idxs) == length(tether2.segment_idxs)
                    # copy the points of the segments of the tethers
                    for (segment_idx1, segment_idx2) in zip(tether1.segment_idxs, tether2.segment_idxs)
                        point_idxs1 = sys1.segments[segment_idx1].point_idxs
                        point_idxs2 = sys2.segments[segment_idx2].point_idxs
                        for (point_idx1, point_idx2) in zip(point_idxs1, point_idxs2)
                            sys2.points[point_idx2].pos_w .= sys1.points[point_idx1].pos_w
                            sys2.points[point_idx2].vel_w .= sys1.points[point_idx1].vel_w
                            sys2.points[point_idx2].disturb .= sys1.points[point_idx1].disturb
                        end
                    end
                elseif length(tether2.segment_idxs) == 1
                    # copy the first and last point of the tether
                    point_idxs1 = [sys1.segments[tether1.segment_idxs[1]].point_idxs[1],
                                   sys1.segments[tether1.segment_idxs[end]].point_idxs[2]]
                    point_idxs2 = sys2.segments[tether2.segment_idxs[1]].point_idxs
                    sys2.points[point_idxs2[1]].pos_w .= sys1.points[point_idxs1[1]].pos_w
                    sys2.points[point_idxs2[2]].pos_w .= sys1.points[point_idxs1[2]].pos_w
                    sys2.points[point_idxs2[1]].vel_w .= sys1.points[point_idxs1[1]].vel_w
                    sys2.points[point_idxs2[2]].vel_w .= sys1.points[point_idxs1[2]].vel_w
                    sys2.points[point_idxs2[1]].disturb .= sys1.points[point_idxs1[1]].disturb
                    sys2.points[point_idxs2[2]].disturb .= sys1.points[point_idxs1[2]].disturb
                    simple = true
                end
            end
        end
    end

    # copy twist and twist_ω of groups
    if length(sys1.groups) > 1 && length(sys1.groups) == length(sys2.groups)
        for (group1, group2) in zip(sys1.groups, sys2.groups)
            group2.twist = group1.twist
            group2.twist_ω = group1.twist_ω
        end
    end

    # copy winch tether lengths and velocities
    if length(sys1.winches) > 1 && length(sys1.winches) == length(sys2.winches)
        for (winch2, winch1) in zip(sys2.winches, sys1.winches)
            if !simple
                winch2.tether_len = winch1.tether_len
                winch2.tether_vel = winch1.tether_vel
            else
                winch2.tether_len = 0.0
                for tether_idx in winch1.tether_idxs
                    tether2 = sys2.tethers[tether_idx]
                    segment2 = sys2.segments[tether2.segment_idxs[1]]
                    point_idxs2 = segment2.point_idxs
                    slen = norm(sys2.points[point_idxs2[1]].pos_w .-
                                        sys2.points[point_idxs2[2]].pos_w)
                    stiffness = segment2.axial_stiffness / slen
                    nt = length(winch1.tether_idxs)
                    winch2.tether_len += (slen - norm(winch1.force)/stiffness/nt) / nt
                end
                winch2.tether_vel = winch1.tether_vel
            end
        end
    end

    # copy pulley lengths and velocities
    if length(sys1.pulleys) > 1 && length(sys1.pulleys) == length(sys2.pulleys)
        for (pulley1, pulley2) in zip(sys1.pulleys, sys2.pulleys)
            pulley2.len = pulley1.len
            pulley2.vel = pulley1.vel
        end
    end

    # copy wing positions and velocities
    if length(sys1.wings) > 1 && length(sys1.wings) == length(sys2.wings)
        for (wing1, wing2) in zip(sys1.wings, sys2.wings)
            wing2.pos_w .= wing1.pos_w
            wing2.vel_w .= wing1.vel_w
            wing2.ω_b .= wing1.ω_b
            wing2.Q_b_w .= wing1.Q_b_w
        end
    end
end

"""
    update_from_sysstate!(sys::SystemStructure, ss::SysState)

Update the dynamic state of a `SystemStructure` from a `SysState` snapshot.

This function copies the state variables that are present in `SysState` (such as point
positions, wing orientations, winch lengths, and twist angles) into an existing `SystemStructure`.
Fields that cannot be populated from `SysState` (such as aerodynamic forces, moments, and
segment forces) are set to `NaN` to prevent them from being plotted.

This is useful for visualizing a `SysLog` by extracting individual `SysState` snapshots
and applying them to a `SystemStructure` for plotting with the Makie extension.

# Arguments
- `sys::SystemStructure`: The system structure to update (must already exist with correct topology).
- `ss::SysState`: The state snapshot to copy from.

# Example
```julia
# Load a system log
lg = load_log(...)

# Create a SystemStructure with the same topology
sys = SystemStructure(se(), "ram")

# Update from a specific time step
update_from_sysstate!(sys, lg.syslog[100])

# Plot the system at that time step
plot(sys)
```

# Notes
- The `SystemStructure` must have been created with the same model configuration as the
  simulation that generated the `SysLog`.
- Aerodynamic and force fields are set to `NaN` and will not be plotted.
- The number of points in `sys` must match the parametric type `P` of `SysState{P}`.
"""
function update_from_sysstate!(sys::SystemStructure, ss::SysState{P}) where P
    @unpack points, groups, winches, wings = sys

    # Verify compatibility
    if length(points) != P
        error("SystemStructure has $(length(points)) points but SysState has $P points")
    end

    # Update point positions (X, Y, Z from SysState)
    for point in points
        point.pos_w[1] = ss.X[point.idx]
        point.pos_w[2] = ss.Y[point.idx]
        point.pos_w[3] = ss.Z[point.idx]
        # Set velocity to zero (not available in basic SysState)
        point.vel_w .= 0.0
        # Set forces to NaN (not available in SysState)
        point.force .= NaN
    end

    # Update wing state if wings exist
    if length(wings) > 0 && length(wings) == 1  # Currently only support single-wing systems
        wing = wings[1]

        # Copy orientation quaternion
        wing.Q_b_w .= ss.orient

        # Copy spherical coordinates
        wing.elevation = Float64(ss.elevation)
        wing.azimuth = Float64(ss.azimuth)
        wing.heading = Float64(ss.heading)

        # Compute wing position from average of points (if wings exist)
        # For a typical system, the wing COM is near the bridle attachment
        wing.pos_w .= [mean(ss.X), mean(ss.Y), mean(ss.Z)]

        # Copy velocity if available in vel_kite
        wing.vel_w .= ss.vel_kite

        # Set angular velocity to NaN (turn_rates in SysState, but need conversion)
        wing.ω_b .= ss.turn_rates

        # Set aerodynamic quantities to NaN (to prevent plotting)
        wing.aero_force_b .= NaN
        wing.aero_moment_b .= NaN
        wing.tether_force .= NaN
        wing.tether_moment .= NaN
        wing.va_b .= NaN
        wing.v_wind .= ss.v_wind_kite
        wing.aoa = Float64(ss.AoA)
        wing.course = Float64(ss.course)
        wing.acc_w .= 0.0
        wing.turn_rate .= ss.turn_rates
        wing.turn_acc .= 0.0
    end

    # Update group twist angles
    n_groups = min(length(groups), 4)  # SysState stores up to 4 twist angles
    for i in 1:n_groups
        if i <= length(groups)
            groups[i].twist = Float64(ss.twist_angles[i])
            groups[i].twist_ω = 0.0  # Not available in SysState
            # Set forces/moments to NaN
            groups[i].tether_force = NaN
            groups[i].tether_moment = NaN
            groups[i].aero_moment = NaN
        end
    end

    # Update winch state
    n_winches = min(length(winches), 4)  # SysState stores up to 4 winches
    for i in 1:n_winches
        if i <= length(winches)
            winches[i].tether_len = Float64(ss.l_tether[i])
            winches[i].tether_vel = Float64(ss.v_reelout[i])
            winches[i].force .= NaN  # Force not directly available
            winches[i].friction = NaN
            winches[i].tether_acc = 0.0
            winches[i].set_value = Float64(ss.set_torque[i])
        end
    end

    # Update global wind vector
    sys.wind_vec_gnd .= ss.v_wind_gnd

    # Calculate segment lengths and forces from current positions and velocities
    # Note: velocities are set to zero, so damping term will be zero
    update_segment_forces!(sys)

    return nothing
end
