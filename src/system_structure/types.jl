# Copyright (c) 2025 Bart van de Lint
# SPDX-License-Identifier: MPL-2.0

"""
Basic type definitions for the system structure components.

This file contains enums and struct definitions for:
- SegmentType, DynamicsType, WingType enums
- Point, Group, Segment, Pulley, Tether, Winch structs
"""

# ==================== ENUMS ==================== #

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

# ==================== POINT ==================== #

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
    const aero_force_b::KVec3 # aerodynamic force in body frame (for REFINE WING points)
    const va_b::KVec3 # apparent velocity in body frame (for VSM per-point va)
    const type::DynamicsType
    extra_mass::SimFloat      # User-provided mass
    total_mass::SimFloat      # extra_mass + segment weights (computed during simulation)
    body_frame_damping::KVec3
    world_frame_damping::KVec3
    area::SimFloat
    drag_coeff::SimFloat
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
- `extra_mass::Float64=0.0`: User-provided mass of the point [kg]. Total mass (including segment weights) is computed during simulation.
- `body_frame_damping::Union{Float64,KVec3}=zeros(KVec3)`: Per-axis damping [x,y,z] for bridle points. Scalar applies to all axes.
- `world_frame_damping::Union{Float64,KVec3}=zeros(KVec3)`: Per-axis damping [x,y,z] for world frame damping. Scalar applies to all axes.
- `fix_sphere::Bool=false`: If true, constrains the point to a sphere.
- `fix_static::Bool=false`: If true, dynamically freezes the point (behaves like STATIC).

# Returns
- `Point`: A new `Point` object.
"""
function Point(idx, pos_cad, type;
    wing_idx=1, vel_w=zeros(KVec3), transform_idx=1,
    extra_mass=0.0, body_frame_damping=zeros(KVec3), world_frame_damping=zeros(KVec3),
    area=0.0, drag_coeff=0.0,
    fix_sphere=false, fix_static=false
)
    # Convert scalar damping to vector (broadcast to all axes)
    bf_damp = body_frame_damping isa Real ? SVector{3,SimFloat}(body_frame_damping, body_frame_damping, body_frame_damping) : SVector{3,SimFloat}(body_frame_damping)
    wf_damp = world_frame_damping isa Real ? SVector{3,SimFloat}(world_frame_damping, world_frame_damping, world_frame_damping) : SVector{3,SimFloat}(world_frame_damping)

    Point(idx, transform_idx, wing_idx, pos_cad, zeros(KVec3), zeros(KVec3),
        vel_w, zeros(KVec3), zeros(KVec3), zeros(KVec3), zeros(KVec3), type, extra_mass, 0.0,
        bf_damp, wf_damp, area, drag_coeff,
        fix_sphere, fix_static)
end

# ==================== GROUP ==================== #

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
    unrefined_section_idxs::Vector{Int16}  # Indices of VSM unrefined sections in this group
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
          0.0, 0.0, 0.0, 0.0, 0.0,
          Int16[])
end

"""
    Group(idx, point_idxs, vsm_wing::Wing, gamma,
          type, moment_frac; damping=50.0)

Legacy constructor that calculates geometry from vsm_wing directly.
Kept for backward compatibility with predefined structures.
"""
function Group(idx, point_idxs, vsm_wing::VortexStepMethod.Wing, gamma,
               type, moment_frac; damping=50.0)
    le_pos = [vsm_wing.le_interp[i](gamma) for i in 1:3]
    chord = [vsm_wing.te_interp[i](gamma) for i in 1:3] .- le_pos
    y_airf = normalize([vsm_wing.le_interp[i](gamma-0.01)
        for i in 1:3] - le_pos)
    Group(idx, point_idxs, gamma, le_pos, chord, y_airf,
          type, moment_frac, damping,
          0.0, 0.0, 0.0, 0.0, 0.0,
          Int16[])
end

# ==================== SEGMENT ==================== #

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

# ==================== PULLEY ==================== #

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

# ==================== TETHER ==================== #

"""
    mutable struct Tether

A collection of segments that are controlled together by a winch.

$(TYPEDFIELDS)
"""
mutable struct Tether
    const idx::Int16
    const segment_idxs::Vector{Int16}
    const winch_point_idx::Int16
    stretched_len::SimFloat
end

"""
    Tether(idx, segment_idxs, winch_point_idx)

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
- `winch_point_idx::Int16`: Index of the ground point where tether attaches to winch.

# Returns
- `Tether`: A new `Tether` object.
"""
function Tether(idx, segment_idxs, winch_point_idx)
    return Tether(idx, segment_idxs, winch_point_idx, 0.0)
end

# ==================== WINCH ==================== #

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

# ==================== TRANSFORM ==================== #

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
    elevation::SimFloat  # [rad]
    azimuth::SimFloat    # [rad]
    heading::SimFloat    # [rad]
    elevation_vel::SimFloat  # [rad/s] angular velocity in elevation direction
    azimuth_vel::SimFloat    # [rad/s] angular velocity in azimuth direction
    turn_rate::SimFloat      # [rad/s] angular velocity around radial axis (not yet implemented)
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
        wing_idx=nothing, rot_point_idx=nothing,
        elevation_vel=0.0, azimuth_vel=0.0, turn_rate=0.0)
    (isnothing(wing_idx) == isnothing(rot_point_idx)) && error("Either provide a wing_idx or a rot_point_idx, not both or none.")
    (isnothing(base_pos) == isnothing(base_transform_idx)) && error("Either provide the base_pos or the base_transform_idx, not both or none.")
    (isnothing(base_pos) !== isnothing(base_point_idx)) && error("When providing a base_pos, also provide a base_point_idx.")
    Transform(idx, wing_idx, rot_point_idx, base_point_idx, base_transform_idx,
              elevation, azimuth, heading, elevation_vel, azimuth_vel, turn_rate, base_pos)
end

"""
    Transform(idx, set, base_point_idx; kwargs...)

Constructor helper to create a `Transform` from a `Settings` object.
"""
function Transform(idx, set, base_point_idx; kwargs...)
    elevation_vel = hasfield(typeof(set), :elevation_vels) ? set.elevation_vels[idx] : 0.0
    azimuth_vel = hasfield(typeof(set), :azimuth_vels) ? set.azimuth_vels[idx] : 0.0
    turn_rate = hasfield(typeof(set), :turn_rates) ? set.turn_rates[idx] : 0.0
    Transform(idx, set.elevations[idx], set.azimuths[idx], set.headings[idx];
              base_point_idx, elevation_vel, azimuth_vel, turn_rate, kwargs...)
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
    get_base_pos(transform::Transform, transforms, wings, points)

Get the base position for a given transform, resolving chained transforms if necessary.
"""
function get_base_pos(transform::Transform, transforms, wings, points)
    curr_base_pos = points[transform.base_point_idx].pos_cad
    if !isnothing(transform.base_pos)
        return transform.base_pos, curr_base_pos
    elseif !isnothing(transform.base_transform_idx)
        base_transform = transforms[transform.base_transform_idx]
        return get_rot_pos(base_transform, wings, points), curr_base_pos
    end
end
