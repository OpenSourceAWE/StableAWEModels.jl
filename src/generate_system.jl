# Copyright (c) 2025 Bart van de Lint
# SPDX-License-Identifier: MPL-2.0

# Implementation of the ram air wing model using ModelingToolkit.jl

"""
    calc_angle_of_attack(va_wing_b)

Calculate the angle of attack [rad] from the apparent wind vector `va_wing_b` in the
body frame.
"""
function calc_angle_of_attack(va_wing_b)
    return atan(va_wing_b[3], va_wing_b[1])
end

"""
    sym_normalize(vec)

Symbolic-safe normalization of a vector. Returns `vec / norm(vec)`.
"""
function sym_normalize(vec)
    return vec / norm(vec)
end

"""
    quaternion_to_rotation_matrix(q)

Convert a quaternion `q` (scalar-first format [w, x, y, z]) to a 3x3 rotation
matrix.
"""
function quaternion_to_rotation_matrix(q)
    w, x, y, z = q[1], q[2], q[3], q[4]

    return [
        1-2*(y*y+z*z) 2*(x*y-z*w) 2*(x*z+y*w)
        2*(x*y+z*w) 1-2*(x*x+z*z) 2*(y*z-x*w)
        2*(x*z-y*w) 2*(y*z+x*w) 1-2*(x*x+y*y)
    ]
end

"""
    rotation_matrix_to_quaternion(R)

Convert a 3x3 rotation matrix `R` to a quaternion (scalar-first format [w, x, y, z]).
This implementation is based on the method that avoids division by zero.
"""
function rotation_matrix_to_quaternion(R)
    tr_ = R[1, 1] + R[2, 2] + R[3, 3]

    if tr_ > 0
        S = sqrt(tr_ + 1.0) * 2
        w = 0.25 * S
        x = (R[3, 2] - R[2, 3]) / S
        y = (R[1, 3] - R[3, 1]) / S
        z = (R[2, 1] - R[1, 2]) / S
    elseif (R[1, 1] > R[2, 2]) && (R[1, 1] > R[3, 3])
        S = sqrt(1.0 + R[1, 1] - R[2, 2] - R[3, 3]) * 2
        w = (R[3, 2] - R[2, 3]) / S
        x = 0.25 * S
        y = (R[1, 2] + R[2, 1]) / S
        z = (R[1, 3] + R[3, 1]) / S
    elseif R[2, 2] > R[3, 3]
        S = sqrt(1.0 + R[2, 2] - R[1, 1] - R[3, 3]) * 2
        w = (R[1, 3] - R[3, 1]) / S
        x = (R[1, 2] + R[2, 1]) / S
        y = 0.25 * S
        z = (R[2, 3] + R[3, 2]) / S
    else
        S = sqrt(1.0 + R[3, 3] - R[1, 1] - R[2, 2]) * 2
        w = (R[2, 1] - R[1, 2]) / S
        x = (R[1, 3] + R[3, 1]) / S
        y = (R[2, 3] + R[3, 2]) / S
        z = 0.25 * S
    end

    return [w, x, y, z]
end

# Component accessors for symbolic registration
rotation_matrix_to_quaternion_w(R) = rotation_matrix_to_quaternion(R)[1]
rotation_matrix_to_quaternion_x(R) = rotation_matrix_to_quaternion(R)[2]
rotation_matrix_to_quaternion_y(R) = rotation_matrix_to_quaternion(R)[3]
rotation_matrix_to_quaternion_z(R) = rotation_matrix_to_quaternion(R)[4]

# Register component functions as symbolic
@register_symbolic rotation_matrix_to_quaternion_w(R::AbstractMatrix)
@register_symbolic rotation_matrix_to_quaternion_x(R::AbstractMatrix)
@register_symbolic rotation_matrix_to_quaternion_y(R::AbstractMatrix)
@register_symbolic rotation_matrix_to_quaternion_z(R::AbstractMatrix)

function calc_wind_factor(am::AtmosphericModel, height, set::Settings)
    if set.profile_law == 0
        return 1.0
    else
        return AtmosphericModels.calc_wind_factor(am, height, set.profile_law)
    end
end
@register_symbolic calc_wind_factor(am::AtmosphericModel, height, set::Settings)

#=
The following `get_*` functions are registered for use within ModelingToolkit.jl's
symbolic context. They act as symbolic "shims" or placeholders that allow the
equation generation logic to access fields from the `SystemStructure` (`psys`) and
`Settings` (`pset`) parameter objects during the symbolic construction of the model.
This avoids hard-coding parameters into the equations, allowing them to be changed
at simulation time.
=#
get_pos_w(sys::SystemStructure, idx::Int16) = sys.points[idx].pos_w
@register_array_symbolic get_pos_w(sys::SystemStructure, idx::Int16) begin
    size = (3,)
    eltype = SimFloat
end
get_vel_w(sys::SystemStructure, idx::Int16) = sys.points[idx].vel_w
@register_array_symbolic get_vel_w(sys::SystemStructure, idx::Int16) begin
    size = (3,)
    eltype = SimFloat
end
get_pos_b(sys::SystemStructure, idx::Int16) = sys.points[idx].pos_b
@register_array_symbolic get_pos_b(sys::SystemStructure, idx::Int16) begin
    size = (3,)
    eltype = SimFloat
end
get_wing_pos_w(sys::SystemStructure, idx::Int16) = sys.wings[idx].pos_w
@register_array_symbolic get_wing_pos_w(sys::SystemStructure, idx::Int16) begin
    size = (3,)
    eltype = SimFloat
end
get_wing_vel_w(sys::SystemStructure, idx::Int16) = sys.wings[idx].vel_w
@register_array_symbolic get_wing_vel_w(sys::SystemStructure, idx::Int16) begin
    size = (3,)
    eltype = SimFloat
end
get_Q_b_w(sys::SystemStructure, idx::Int16) = sys.wings[idx].Q_b_w
@register_array_symbolic get_Q_b_w(sys::SystemStructure, idx::Int16) begin
    size = (4,)
    eltype = SimFloat
end
get_ω_b(sys::SystemStructure, idx::Int16) = sys.wings[idx].ω_b
@register_array_symbolic get_ω_b(sys::SystemStructure, idx::Int16) begin
    size = (3,)
    eltype = SimFloat
end
get_wind_disturb(sys::SystemStructure, idx::Int16) = sys.wings[idx].wind_disturb
@register_array_symbolic get_wind_disturb(sys::SystemStructure, idx::Int16) begin
    size = (3,)
    eltype = SimFloat
end
get_disturb(sys::SystemStructure, idx::Int16) = sys.points[idx].disturb
@register_array_symbolic get_disturb(sys::SystemStructure, idx::Int16) begin
    size = (3,)
    eltype = SimFloat
end
get_le_pos(sys::SystemStructure, idx::Int16) = sys.groups[idx].le_pos
@register_array_symbolic get_le_pos(sys::SystemStructure, idx::Int16) begin
    size = (3,)
    eltype = SimFloat
end
get_vsm_y(sys::SystemStructure, idx::Int16, iy::Int) = sys.wings[idx].vsm_y[iy]
@register_symbolic get_vsm_y(sys::SystemStructure, idx::Int16, iy::Int)
get_vsm_x(sys::SystemStructure, idx::Int16, ix::Int) = sys.wings[idx].vsm_x[ix]
@register_symbolic get_vsm_x(sys::SystemStructure, idx::Int16, ix::Int)
get_vsm_jac(sys::SystemStructure, idx::Int16, ix::Int, iy::Int) =
    sys.wings[idx].vsm_jac[ix, iy]
@register_symbolic get_vsm_jac(sys::SystemStructure, idx::Int16, ix::Int, iy::Int)
get_pulley_len(sys::SystemStructure, idx::Int16) = sys.pulleys[idx].len
@register_symbolic get_pulley_len(sys::SystemStructure, idx::Int16)
get_pulley_vel(sys::SystemStructure, idx::Int16) = sys.pulleys[idx].vel
@register_symbolic get_pulley_vel(sys::SystemStructure, idx::Int16)
get_set_value(sys::SystemStructure, idx::Int16) = sys.winches[idx].set_value
@register_symbolic get_set_value(sys::SystemStructure, idx::Int16)
get_twist(sys::SystemStructure, idx::Int16) = sys.groups[idx].twist
@register_symbolic get_twist(sys::SystemStructure, idx::Int16)
get_group_damping(sys::SystemStructure, idx::Int16) = sys.groups[idx].damping
@register_symbolic get_group_damping(sys::SystemStructure, idx::Int16)
get_twist_ω(sys::SystemStructure, idx::Int16) = sys.groups[idx].twist_ω
@register_symbolic get_twist_ω(sys::SystemStructure, idx::Int16)
get_mass(sys::SystemStructure, idx::Int16) = sys.points[idx].mass
@register_symbolic get_mass(sys::SystemStructure, idx::Int16)
get_l0(sys::SystemStructure, idx::Int16) = sys.segments[idx].l0
@register_symbolic get_l0(sys::SystemStructure, idx::Int16)
get_diameter(sys::SystemStructure, idx::Int16) = sys.segments[idx].diameter
@register_symbolic get_diameter(sys::SystemStructure, idx::Int16)
get_compression_frac(sys::SystemStructure, idx::Int16) =
    sys.segments[idx].compression_frac
@register_symbolic get_compression_frac(sys::SystemStructure, idx::Int16)
get_moment_frac(sys::SystemStructure, idx::Int16) = sys.groups[idx].moment_frac
@register_symbolic get_moment_frac(sys::SystemStructure, idx::Int16)
get_sum_len(sys::SystemStructure, idx::Int16) = sys.pulleys[idx].sum_len
@register_symbolic get_sum_len(sys::SystemStructure, idx::Int16)
get_tether_len(sys::SystemStructure, idx::Int16) = sys.winches[idx].tether_len
@register_symbolic get_tether_len(sys::SystemStructure, idx::Int16)
get_tether_vel(sys::SystemStructure, idx::Int16) = sys.winches[idx].tether_vel
@register_symbolic get_tether_vel(sys::SystemStructure, idx::Int16)
get_axial_stiffness(sys::SystemStructure, idx::Int16) =
    sys.segments[idx].axial_stiffness
@register_symbolic get_axial_stiffness(sys::SystemStructure, idx::Int16)
get_axial_damping(sys::SystemStructure, idx::Int16) =
    sys.segments[idx].axial_damping
@register_symbolic get_axial_damping(sys::SystemStructure, idx::Int16)
get_body_frame_damping(sys::SystemStructure, idx::Int16) = sys.points[idx].body_frame_damping
@register_symbolic get_body_frame_damping(sys::SystemStructure, idx::Int16)
get_world_frame_damping(sys::SystemStructure, idx::Int16) = sys.points[idx].world_frame_damping
@register_symbolic get_world_frame_damping(sys::SystemStructure, idx::Int16)
get_brake(sys::SystemStructure, idx::Int16) = sys.winches[idx].brake
@register_symbolic get_brake(sys::SystemStructure, idx::Int16)
get_fix_point_sphere(sys::SystemStructure, idx::Int16) =
    sys.points[idx].fix_sphere
@register_symbolic get_fix_point_sphere(sys::SystemStructure, idx::Int16)
get_fix_wing_sphere(sys::SystemStructure, idx::Int16) = sys.wings[idx].fix_sphere
@register_symbolic get_fix_wing_sphere(sys::SystemStructure, idx::Int16)
get_drag_frac(sys::SystemStructure, idx::Int16) = sys.wings[idx].drag_frac
@register_symbolic get_drag_frac(sys::SystemStructure, idx::Int16)
get_y_damping(sys::SystemStructure, idx::Int16) = sys.wings[idx].y_damping
@register_symbolic get_y_damping(sys::SystemStructure, idx::Int16)
get_z_disturb(sys::SystemStructure, idx::Int16) = sys.wings[idx].z_disturb
@register_symbolic get_z_disturb(sys::SystemStructure, idx::Int16)
get_winch_gear_ratio(sys::SystemStructure, idx::Int16) =
    sys.winches[idx].gear_ratio
@register_symbolic get_winch_gear_ratio(sys::SystemStructure, idx::Int16)
get_winch_drum_radius(sys::SystemStructure, idx::Int16) =
    sys.winches[idx].drum_radius
@register_symbolic get_winch_drum_radius(sys::SystemStructure, idx::Int16)
get_winch_f_coulomb(sys::SystemStructure, idx::Int16) =
    sys.winches[idx].f_coulomb
@register_symbolic get_winch_f_coulomb(sys::SystemStructure, idx::Int16)
get_winch_c_vf(sys::SystemStructure, idx::Int16) = sys.winches[idx].c_vf
@register_symbolic get_winch_c_vf(sys::SystemStructure, idx::Int16)
get_winch_inertia_total(sys::SystemStructure, idx::Int16) =
    sys.winches[idx].inertia_total
@register_symbolic get_winch_inertia_total(sys::SystemStructure, idx::Int16)
@register_symbolic get_wind_elevation(sys::SystemStructure)
get_wind_elevation(sys::SystemStructure) = sys.wind_elevation
get_set_mass(set::Settings) = set.mass
@register_symbolic get_set_mass(set::Settings)
get_rho_tether(set::Settings) = set.rho_tether
@register_symbolic get_rho_tether(set::Settings)
get_cd_tether(set::Settings) = set.cd_tether
@register_symbolic get_cd_tether(set::Settings)
get_v_wind(set::Settings) = set.v_wind
@register_symbolic get_v_wind(set::Settings)
get_upwind_dir(set::Settings) = set.upwind_dir
@register_symbolic get_upwind_dir(set::Settings)
get_g_earth(set::Settings) = set.g_earth
@register_symbolic get_g_earth(set::Settings)

"""
    force_eqs!(s, system, psys, pset, eqs, defaults, guesses; kwargs...)

Generate the force and constraint equations for the mass-spring-damper components.

This function builds the core dynamic and kinematic equations for all parts of the
tether and bridle system, including points, segments, pulleys, winches, and the
deformable groups on the wing.

# Arguments
- `s::SymbolicAWEModel`: The main model object.
- `system::SystemStructure`: The physical structure definition.
- `psys`, `pset`: Symbolic parameters representing `system` and `s.set`.
- `eqs`, `defaults`, `guesses`: The accumulating vectors for the MTK system.
- `kwargs...`: Symbolic variables for the system's state (e.g., `R_b_w`, `wing_pos`).

# Returns
- A tuple `(eqs, defaults, guesses, tether_wing_force, tether_wing_moment)`
  containing the updated equation lists and the aggregate forces and moments exerted
  by the tethers on the wing.
"""
function force_eqs!(
    s, system, psys, pset, eqs, defaults, guesses;
    R_b_w, wing_pos, wing_vel, wind_vec_gnd, group_aero_moment,
    twist_angle, twist_ω, set_values, fix_wing
)

    @unpack points, groups, segments, pulleys, tethers, winches, wings = system

    # Aggregate forces and moments from tethers onto the wing's center of mass
    tether_wing_force = zeros(Num, length(wings), 3)
    tether_wing_moment = zeros(Num, length(wings), 3)

    # Check if we have any REFINE wings (need aero force per point)
    has_refine_wings = any(wing.wing_type == REFINE for wing in wings)

    @variables begin
        # Point states
        pos(t)[1:3, eachindex(points)]
        vel(t)[1:3, eachindex(points)]
        acc(t)[1:3, eachindex(points)]
        # Point forces and geometry
        point_force(t)[1:3, eachindex(points)]
        disturb_force(t)[1:3, eachindex(points)]
        tether_r(t)[1:3, eachindex(points)]
        point_mass(t)[eachindex(points)]
        chord_b(t)[1:3, eachindex(points)]
        fixed_pos(t)[1:3, eachindex(points)]
        normal(t)[1:3, eachindex(points)]
        pos_b(t)[1:3, eachindex(points)]
        fix_point_sphere(t)[eachindex(points)]
        # Segment forces and rest length
        spring_force_vec(t)[1:3, eachindex(segments)]
        drag_force(t)[1:3, eachindex(segments)]
        l0(t)[eachindex(segments)]
    end

    # REFINE-specific variables: per-point aero forces in body frame
    # (No segment twists - no linearization for REFINE wings)
    if has_refine_wings
        @variables begin
            aero_force_point_b(t)[eachindex(points), 1:3]
        end
    end

    # ==================== POINTS ==================== #
    for point in points
        F::Vector{Num} = zeros(Num, 3)
        mass = get_mass(psys, point.idx)
        for segment in segments
            if point.idx in segment.point_idxs
                mass_per_meter =
                    get_rho_tether(pset) * π *
                    (get_diameter(psys, segment.idx) / 2)^2
                inverted = segment.point_idxs[2] == point.idx
                if inverted
                    F .-= spring_force_vec[:, segment.idx]
                else
                    F .+= spring_force_vec[:, segment.idx]
                end
                mass += mass_per_meter * l0[segment.idx] / 2
                F .+= 0.5 * drag_force[:, segment.idx]
            end
        end

        # The net force on the point. This variable is used by other components.
        eqs = [
            eqs
            point_mass[point.idx] ~ mass
            disturb_force[:, point.idx] ~ get_disturb(psys, point.idx)
        ]

        if point.type == WING
            # Find the wing for this point
            wing = wings[point.wing_idx]

            if wing.wing_type == REFINE
                # REFINE wing: Points are DYNAMIC and receive lumped panel forces
                # Similar to DYNAMIC points but with aero forces included

                # Add aerodynamic forces (calculated in linear_vsm_eqs!)
                aero_force_w = R_b_w[wing.idx, :, :] * aero_force_point_b[point.idx, :]

                eqs = [
                    eqs
                    point_force[:, point.idx] ~
                        F + aero_force_w + [0, 0, -get_g_earth(pset) * mass] + disturb_force[:, point.idx]
                ]

                # Damping terms
                body_frame_damp_vec = get_body_frame_damping(psys, point.idx) * (vel[:, point.idx] - wing_vel[point.wing_idx, :])
                world_frame_damp_vec = get_world_frame_damping(psys, point.idx) * vel[:, point.idx]

                # DYNAMIC point equations
                axis = sym_normalize(pos[:, point.idx])
                eqs = [
                    eqs
                    fix_point_sphere[point.idx] ~ get_fix_point_sphere(psys, point.idx)
                    D(pos[:, point.idx]) ~ ifelse.(fix_point_sphere[point.idx]==true,
                                                    vel[:, point.idx] ⋅ axis * axis,
                                                    vel[:, point.idx]
                                            )
                    D(vel[:, point.idx]) ~ ifelse.(fix_point_sphere[point.idx]==true,
                                                    acc[:, point.idx] ⋅ axis * axis,
                                                    acc[:, point.idx]
                                            )
                    acc[:, point.idx] ~ point_force[:, point.idx] ./ mass - body_frame_damp_vec - world_frame_damp_vec
                ]
                defaults = [
                    defaults
                    [pos[j, point.idx] => get_pos_w(psys, point.idx)[j] for j = 1:3]
                    [vel[j, point.idx] => get_vel_w(psys, point.idx)[j] for j = 1:3]
                ]

            else  # QUATERNION wing (current implementation)
                # Define point_force for QUATERNION wing points
                eqs = [
                    eqs
                    point_force[:, point.idx] ~
                        F + [0, 0, -get_g_earth(pset) * mass] + disturb_force[:, point.idx]
                ]

                found = 0
                group = nothing
                for group_ in groups
                    if point.idx in group_.point_idxs
                        group = group_
                        found += 1
                    end
                end
                !(found in [0, 1]) && error(
                    "Kite point number $(point.idx) is part of $found groups, " *
                    "and should be part of exactly 0 or 1 groups.",
                )

                if found == 1
                    found = 0
                    wing = nothing
                    for wing_ in wings
                        if group.idx in wing_.group_idxs
                            wing = wing_
                            found += 1
                        end
                    end
                    !(found == 1) && error(
                        "Kite group number $(group.idx) is part of $found wings, " *
                        "and should be part of exactly 1 wing.",
                    )

                    eqs = [
                        eqs
                        fixed_pos[:, point.idx] ~ get_le_pos(psys, group.idx)
                        chord_b[:, point.idx] ~
                            get_pos_b(psys, point.idx) .- fixed_pos[:, point.idx]
                        normal[:, point.idx] ~ chord_b[:, point.idx] × group.y_airf
                        pos_b[:, point.idx] ~
                            fixed_pos[:, point.idx] .+
                            cos(twist_angle[group.idx]) * chord_b[:, point.idx] -
                            sin(twist_angle[group.idx]) * normal[:, point.idx]
                    ]
                elseif found == 0
                    eqs = [eqs; pos_b[:, point.idx] ~ get_pos_b(psys, point.idx)]
                end
                eqs = [
                    eqs
                    tether_r[:, point.idx] ~ pos[:, point.idx] -
                                            wing_pos[point.wing_idx, :]
                ]
                tether_wing_moment[point.wing_idx, :] .+=
                    tether_r[:, point.idx] × point_force[:, point.idx]
                tether_wing_force[point.wing_idx, :] .+= point_force[:, point.idx]

                eqs = [
                    eqs
                    pos[:, point.idx] ~
                        wing_pos[point.wing_idx, :] +
                        R_b_w[point.wing_idx, :, :] * pos_b[:, point.idx]
                    vel[:, point.idx] ~ zeros(3)
                    acc[:, point.idx] ~ zeros(3)
                ]
            end
        elseif point.type == STATIC
            # Define point_force for STATIC points
            eqs = [
                eqs
                point_force[:, point.idx] ~
                    F + [0, 0, -get_g_earth(pset) * mass] + disturb_force[:, point.idx]
                pos[:, point.idx] ~ get_pos_w(psys, point.idx)
                vel[:, point.idx] ~ zeros(3)
                acc[:, point.idx] ~ zeros(3)
            ]
        elseif point.type == DYNAMIC
            # Define point_force for DYNAMIC points
            eqs = [
                eqs
                point_force[:, point.idx] ~
                    F + [0, 0, -get_g_earth(pset) * mass] + disturb_force[:, point.idx]
            ]

            if length(wings) > 0
                body_frame_damp_vec = get_body_frame_damping(psys, point.idx) * (vel[:, point.idx] - wing_vel[point.wing_idx, :])
            else
                body_frame_damp_vec = zeros(3)
            end
            world_frame_damp_vec = get_world_frame_damping(psys, point.idx) * vel[:, point.idx]

            axis = sym_normalize(pos[:, point.idx])
            eqs = [
                eqs
                fix_point_sphere[point.idx] ~ get_fix_point_sphere(psys, point.idx)
                D(pos[:, point.idx]) ~  ifelse.(fix_point_sphere[point.idx]==true,
                                                vel[:, point.idx] ⋅ axis * axis,
                                                vel[:, point.idx]
                                        )
                D(vel[:, point.idx]) ~  ifelse.(fix_point_sphere[point.idx]==true,
                                                acc[:, point.idx] ⋅ axis * axis,
                                                acc[:, point.idx]
                                        )
                acc[:, point.idx]    ~ point_force[:, point.idx] ./ mass - body_frame_damp_vec - world_frame_damp_vec
            ]
            defaults = [
                defaults
                [pos[j, point.idx] => get_pos_w(psys, point.idx)[j] for j = 1:3]
                [vel[j, point.idx] => get_vel_w(psys, point.idx)[j] for j = 1:3]
            ]
        elseif point.type == QUASI_STATIC
            # Define point_force for QUASI_STATIC points
            eqs = [
                eqs
                point_force[:, point.idx] ~
                    F + [0, 0, -get_g_earth(pset) * mass] + disturb_force[:, point.idx]
                vel[:, point.idx] ~ zeros(3)
                acc[:, point.idx] ~ zeros(3)
                # For a quasi-static point, the net force must be zero. This becomes
                # the algebraic equation that determines the point's position.
                zeros(3) .~ point_force[:, point.idx]
            ]
            guesses = [
                guesses
                [acc[j, point.idx] => 0 for j = 1:3]
                [pos[j, point.idx] => get_pos_w(psys, point.idx)[j] for j = 1:3]
                [point_force[j, point.idx] => 0 for j = 1:3]
            ]
        else
            error("Unknown point type: $(typeof(point))")
        end
    end

    # REFINE wings do not use segment twist calculations
    # (no linearization, forces come directly from nonlinear VSM)

    # ==================== GROUPS ==================== #
    if length(groups) > 0
        @variables begin
            trailing_edge_angle(t)[eachindex(groups)]
            trailing_edge_ω(t)[eachindex(groups)]
            trailing_edge_α(t)[eachindex(groups)]
            free_twist_angle(t)[eachindex(groups)]
            twist_α(t)[eachindex(groups)]
            group_tether_force(t)[eachindex(groups)]
            group_tether_moment(t)[eachindex(groups)]
            tether_force(t)[eachindex(groups), eachindex(groups[1].point_idxs)]
            tether_moment(t)[eachindex(groups), eachindex(groups[1].point_idxs)]
            r_group(t)[eachindex(groups), eachindex(groups[1].point_idxs)]
            r_vec(t)[eachindex(groups), eachindex(groups[1].point_idxs), 1:3]
        end
    end

    for group in groups
        found = 0
        wing = nothing
        for wing_ in wings
            if group.idx in wing_.group_idxs
                wing = wing_
                found += 1
            end
        end
        !(found == 1) && error(
            "Kite group $(group.idx) is in $found wings; must be in exactly 1.",
        )

        all(iszero.(tether_wing_moment[wing.idx, :])) && error(
            "Tether wing moment is zero. At least one wing connection point " *
            "should not be part of a deforming group.",
        )

        x_airf = normalize(group.chord)
        init_z_airf = x_airf × group.y_airf
        z_airf =
            x_airf * sin(twist_angle[group.idx]) +
            init_z_airf * cos(twist_angle[group.idx])
        for (i, point_idx) in enumerate(group.point_idxs)
            eqs = [
                eqs
                r_vec[group.idx, i, :] ~ (
                    get_pos_b(psys, point_idx) .-
                    (group.le_pos + get_moment_frac(psys, group.idx) * group.chord)
                )
                r_group[group.idx, i] ~
                    r_vec[group.idx, i, :] ⋅ normalize(group.chord)
                tether_force[group.idx, i] ~
                    (point_force[:, point_idx] ⋅ (R_b_w[wing.idx, :, :] * -z_airf))
                tether_moment[group.idx, i] ~
                    r_group[group.idx, i] * tether_force[group.idx, i]
            ]
        end

        # Inertia of a thin rectangular plate rotating around one edge
        inertia =
            1 / 3 * (get_set_mass(pset) / length(groups)) * (norm(group.chord))^2
        @parameters max_twist = deg2rad(90)

        eqs = [
            eqs
            group_tether_force[group.idx] ~ sum(tether_force[group.idx, :])
            group_tether_moment[group.idx] ~ sum(tether_moment[group.idx, :])
            twist_α[group.idx] ~
                (group_aero_moment[group.idx] + group_tether_moment[group.idx]) /
                inertia
            twist_angle[group.idx] ~
                clamp(free_twist_angle[group.idx], -max_twist, max_twist)
        ]
        if group.type == DYNAMIC
            eqs = [
                eqs
                D(free_twist_angle[group.idx]) ~
                    ifelse(fix_wing == true, 0, twist_ω[group.idx])
                D(twist_ω[group.idx]) ~ ifelse(
                    fix_wing == true,
                    0,
                    twist_α[group.idx] -
                    get_group_damping(psys, group.idx) * twist_ω[group.idx],
                )
            ]
            defaults = [
                defaults
                free_twist_angle[group.idx] => get_twist(psys, group.idx)
                twist_ω[group.idx] => get_twist_ω(psys, group.idx)
            ]
        elseif group.type == QUASI_STATIC
            eqs = [eqs; twist_ω[group.idx] ~ 0; twist_α[group.idx] ~ 0]
            guesses = [
                guesses
                free_twist_angle[group.idx] => 0
                twist_angle[group.idx] => 0
            ]
        else
            error("Wrong group type.")
        end
    end

    # ==================== SEGMENTS ==================== #
    @variables begin
        # Spring-damper model
        segment_vec(t)[1:3, eachindex(segments)]
        unit_vec(t)[1:3, eachindex(segments)]
        len(t)[eachindex(segments)]
        rel_vel(t)[1:3, eachindex(segments)]
        spring_vel(t)[eachindex(segments)]
        spring_force(t)[eachindex(segments)]
        stiffness(t)[eachindex(segments)]
        damping(t)[eachindex(segments)]
        # Aerodynamic drag model
        height(t)[eachindex(segments)]
        segment_vel(t)[1:3, eachindex(segments)]
        segment_rho(t)[eachindex(segments)]
        wind_vel(t)[1:3, eachindex(segments)]
        va(t)[1:3, eachindex(segments)]
        area(t)[eachindex(segments)]
        app_perp_vel(t)[1:3, eachindex(segments)]
        # Related component states
        pulley_len(t)[eachindex(pulleys)]
        tether_len(t)[eachindex(winches)]
    end
    for segment in segments
        p1, p2 = segment.point_idxs[1], segment.point_idxs[2]
        guesses = [
            guesses
            [
                segment_vec[i, segment.idx] =>
                    get_pos_w(psys, p2)[i] - get_pos_w(psys, p1)[i] for i = 1:3
            ]
        ]

        # Check if both endpoints are WING points (structural segments within REFINE wings)
        # Only skip aerodynamic drag for segments between two WING points
        # Bridle segments attached to WING points should still have drag
        p1_obj = points[p1]
        p2_obj = points[p2]
        is_wing_structural_segment = (p1_obj.type == WING && p2_obj.type == WING)

        in_pulley = 0
        for pulley in pulleys
            if segment.idx == pulley.segment_idxs[1]
                eqs = [eqs; l0[segment.idx] ~ pulley_len[pulley.idx]]
                in_pulley += 1
            end
            if segment.idx == pulley.segment_idxs[2]
                eqs = [
                    eqs
                    l0[segment.idx] ~
                        get_sum_len(psys, pulley.idx) - pulley_len[pulley.idx]
                ]
                in_pulley += 1
            end
        end
        (in_pulley > 1) && error(
            "Bridle segment $(segment.idx) is in $in_pulley pulleys; " *
            "should be in 0 or 1.",
        )

        if in_pulley == 0
            in_tether = 0
            for tether in tethers
                if segment.idx in tether.segment_idxs
                    in_winch = 0
                    winch_idx = 0
                    for winch in winches
                        if tether.idx in winch.tether_idxs
                            winch_idx = winch.idx
                            in_winch += 1
                        end
                    end
                    (in_winch != 1) && error(
                        "Tether $(tether.idx) is connected to $in_winch " *
                        "winches; should be 1.",
                    )

                    eqs = [
                        eqs
                        l0[segment.idx] ~
                            tether_len[winch_idx] / length(tether.segment_idxs)
                    ]
                    in_tether += 1
                end
            end
            !(in_tether in [0, 1]) && error(
                "Segment $(segment.idx) is in $in_tether tethers; should be 0 or 1.",
            )

            if in_tether == 0
                eqs = [eqs; l0[segment.idx] ~ get_l0(psys, segment.idx)]
            end
        end

        eqs = [
            eqs
            # Spring force equations (Hooke's Law with damping)
            segment_vec[:, segment.idx] ~ pos[:, p2] - pos[:, p1]
            len[segment.idx] ~ norm(segment_vec[:, segment.idx])
            unit_vec[:, segment.idx] ~
                segment_vec[:, segment.idx] / len[segment.idx]
            rel_vel[:, segment.idx] ~ vel[:, p1] - vel[:, p2]
            spring_vel[segment.idx] ~
                rel_vel[:, segment.idx] ⋅ unit_vec[:, segment.idx]
            damping[segment.idx] ~
                get_axial_damping(psys, segment.idx) / len[segment.idx]
            stiffness[segment.idx] ~ ifelse(
                len[segment.idx] > l0[segment.idx],
                get_axial_stiffness(psys, segment.idx) / len[segment.idx],
                get_compression_frac(psys, segment.idx) *
                get_axial_stiffness(psys, segment.idx) / len[segment.idx],
            )
            spring_force[segment.idx] ~ (
                stiffness[segment.idx] * (len[segment.idx] - l0[segment.idx]) -
                damping[segment.idx] * spring_vel[segment.idx]
            )
            spring_force_vec[:, segment.idx] ~
                spring_force[segment.idx] * unit_vec[:, segment.idx]
        ]

        # Skip aerodynamic drag for structural segments within REFINE wings
        # (those segments connect WING points and get forces from VSM panels instead)
        if !is_wing_structural_segment
            eqs = [
                eqs
                # Aerodynamic drag force on the tether segment
                height[segment.idx] ~ max(0.0, 0.5 * (pos[:, p1][3] + pos[:, p2][3]))
                segment_vel[:, segment.idx] ~ 0.5 * (vel[:, p1] + vel[:, p2])
                segment_rho[segment.idx] ~ calc_rho(s.am, height[segment.idx])
                wind_vel[:, segment.idx] ~
                    calc_wind_factor(s.am, max(height[segment.idx], 1.0), pset) *
                    wind_vec_gnd
                va[:, segment.idx] ~
                    wind_vel[:, segment.idx] - segment_vel[:, segment.idx]
                area[segment.idx] ~
                    len[segment.idx] * get_diameter(psys, segment.idx)
                app_perp_vel[:, segment.idx] ~
                    va[:, segment.idx] -
                    (va[:, segment.idx] ⋅ unit_vec[:, segment.idx]) *
                    unit_vec[:, segment.idx]
                drag_force[:, segment.idx] ~
                    (
                        0.5 * segment_rho[segment.idx] * get_cd_tether(pset) *
                        norm(va[:, segment.idx]) * area[segment.idx]
                    ) * app_perp_vel[:, segment.idx]
            ]
        else
            # For wing structural segments, only define drag_force (set to zero)
            # Don't define other aerodynamic variables - they're unused and would create extra equations
            eqs = [
                eqs
                drag_force[:, segment.idx] ~ zeros(3)
            ]
        end
    end

    # ==================== PULLEYS ==================== #
    @variables begin
        pulley_len(t)[eachindex(pulleys)]
        pulley_vel(t)[eachindex(pulleys)]
        pulley_force(t)[eachindex(pulleys)]
        pulley_acc(t)[eachindex(pulleys)]
    end
    @parameters pulley_damp = 5.0
    for pulley in pulleys
        segment = segments[pulley.segment_idxs[1]]
        mass_per_meter =
            get_rho_tether(pset) * π * (get_diameter(psys, segment.idx) / 2)^2
        mass = get_sum_len(psys, pulley.idx) * mass_per_meter
        eqs = [
            eqs
            pulley_force[pulley.idx] ~
                spring_force[pulley.segment_idxs[1]] -
                spring_force[pulley.segment_idxs[2]]
            pulley_acc[pulley.idx] ~ pulley_force[pulley.idx] / mass
        ]
        if pulley.type == DYNAMIC
            eqs = [
                eqs
                D(pulley_len[pulley.idx]) ~ pulley_vel[pulley.idx]
                D(pulley_vel[pulley.idx]) ~
                    pulley_acc[pulley.idx] - pulley_damp * pulley_vel[pulley.idx]
            ]
            defaults = [
                defaults
                pulley_len[pulley.idx] => get_pulley_len(psys, pulley.idx)
                pulley_vel[pulley.idx] => get_pulley_vel(psys, pulley.idx)
            ]
        elseif pulley.type == QUASI_STATIC
            eqs = [eqs; pulley_vel[pulley.idx] ~ 0; pulley_acc[pulley.idx] ~ 0]
            guesses = [
                guesses
                pulley_len[pulley.idx] => get_l0(psys, pulley.segment_idxs[1])
            ]
        else
            error("Wrong pulley type")
        end
    end

    # ==================== WINCHES ==================== #
    @variables begin
        tether_vel(t)[eachindex(winches)]
        tether_acc(t)[eachindex(winches)]
        winch_force(t)[eachindex(winches)]
        winch_force_vec(t)[1:3, eachindex(winches)]
        brake(t)[eachindex(winches)]
        # Winch motor and friction dynamics
        ω_motor(t)[eachindex(winches)]
        tau_friction(t)[eachindex(winches)]
        tau_motor(t)[eachindex(winches)]
        tau_total(t)[eachindex(winches)]
        α_motor(t)[eachindex(winches)]
    end
    for winch in winches
        F = zeros(Num, 3)
        for tether_idx in winch.tether_idxs
            winch_idx = tethers[tether_idx].winch_idx
            (winch_idx > length(points)) &&
                error("Point number $winch_idx does not exist.")
            F .+= point_force[:, winch_idx]
        end

        gear_ratio = get_winch_gear_ratio(psys, winch.idx)
        drum_radius = get_winch_drum_radius(psys, winch.idx)
        f_coulomb = get_winch_f_coulomb(psys, winch.idx)
        c_vf = get_winch_c_vf(psys, winch.idx)
        inertia_total = get_winch_inertia_total(psys, winch.idx)

        # Smooth sign function to avoid discontinuities at zero velocity
        function smooth_sign(x)
            EPSILON = 6
            x / sqrt(x * x + EPSILON * EPSILON)
        end

        eqs = [
            eqs
            brake[winch.idx] ~ get_brake(psys, winch.idx)
            D(tether_len[winch.idx]) ~
                ifelse(brake[winch.idx] == true, 0, tether_vel[winch.idx])
            D(tether_vel[winch.idx]) ~
                ifelse(brake[winch.idx] == true, 0, tether_acc[winch.idx])

            # Winch motor, gear, and friction dynamics
            ω_motor[winch.idx] ~
                gear_ratio / drum_radius * tether_vel[winch.idx]
            tau_friction[winch.idx] ~
                smooth_sign(ω_motor[winch.idx]) * f_coulomb * drum_radius /
                gear_ratio +
                c_vf * ω_motor[winch.idx] * drum_radius^2 / gear_ratio^2
            tau_motor[winch.idx] ~ set_values[winch.idx]
            tau_total[winch.idx] ~
                tau_motor[winch.idx] +
                drum_radius / gear_ratio * winch_force[winch.idx] -
                tau_friction[winch.idx]
            α_motor[winch.idx] ~ tau_total[winch.idx] / inertia_total
            tether_acc[winch.idx] ~
                drum_radius / gear_ratio * α_motor[winch.idx]

            winch_force_vec[:, winch.idx] ~ F
            winch_force[winch.idx] ~ norm(winch_force_vec[:, winch.idx])
        ]
        defaults = [
            defaults
            tether_len[winch.idx] => get_tether_len(psys, winch.idx)
            tether_vel[winch.idx] => get_tether_vel(psys, winch.idx)
        ]
    end

    # ==================== TETHERS ==================== #
    @variables begin
        stretched_len(t)[eachindex(tethers)]
        tether_spring_force(t)[eachindex(tethers)]
    end
    for tether in tethers
        slen = zero(Num)
        tforce = zero(Num)
        for segment_idx in tether.segment_idxs
            slen += len[segment_idx]
            tforce += spring_force[segment_idx]
        end
        tforce /= length(tether.segment_idxs)
        eqs = [
            eqs
            stretched_len[tether.idx] ~ slen
            tether_spring_force[tether.idx] ~ tforce
        ]
    end

    # Return aero_force_point_b if REFINE wings exist, otherwise nothing
    aero_force_point_b_ret = has_refine_wings ? aero_force_point_b : nothing
    return eqs, defaults, guesses, tether_wing_force, tether_wing_moment, aero_force_point_b_ret, pos, vel, acc
end

"""
    wing_eqs!(s, eqs, psys, pset, defaults; kwargs...)

Generate the differential equations for the wing's rigid body dynamics.

This function builds the equations for:
- Quaternion kinematics for the wing's orientation.
- Euler's rotation equations for the angular acceleration.
- Newton's second law for the translational motion of the wing's center of mass.

# Arguments
- `s::SymbolicAWEModel`: The main model object.
- `eqs`, `psys`, `pset`, `defaults`: Accumulating vectors and symbolic parameters.
- `kwargs...`: Symbolic variables for forces, moments, and states.

# Returns
- A tuple `(eqs, defaults)` containing the updated equation and default value lists.
"""
function wing_eqs!(
    s, eqs, psys, pset, defaults;
    tether_wing_force, tether_wing_moment, aero_force_b,
    aero_moment_b, ω_b, α_b, R_b_w, wing_pos, wing_vel, wing_acc, fix_wing,
    pos, vel, acc
)
    wings = s.sys_struct.wings

    # Check if we have any QUATERNION wings (REFINE wings don't need rigid body dynamics)
    has_quaternion = any(w.wing_type == QUATERNION for w in wings)

    # Always declare quaternion variables (REFINE wings will compute Q_b_w from R_b_w)
    @variables begin
        # Intermediate variables for dynamics (only used for QUATERNION wings)
        wing_acc_b(t)[eachindex(wings), 1:3]
        α_b_damped(t)[eachindex(wings), 1:3]
        ω_b_stable(t)[eachindex(wings), 1:3]
        # Orientation states (all wings)
        Q_b_w(t)[eachindex(wings), 1:4]
        Q_vel(t)[eachindex(wings), 1:4]
        # Forces, moments, and other properties
        moment_b(t)[eachindex(wings), 1:3]
        moment_tether_wing(t)[eachindex(wings), 1:3]
        force_tether_wing(t)[eachindex(wings), 1:3]
        wing_mass(t)[eachindex(wings)]
        fix_wing_sphere(t)[eachindex(wings)]
    end

    # Skew-symmetric matrix for quaternion kinematics
    Ω(ω) = [
        0 -ω[1] -ω[2] -ω[3]
        ω[1] 0 ω[3] -ω[2]
        ω[2] -ω[3] 0 ω[1]
        ω[3] ω[2] -ω[1] 0
    ]

    # Helper to get position (single point or average of multiple)
    # Used for REFINE wing reference points
    get_ref_position(pos, ref::Int16) = pos[:, ref]
    function get_ref_position(pos, refs::Vector{Int16})
        n = length(refs)
        return sum(pos[:, idx] for idx in refs) / n
    end

    for wing in wings
        # REFINE wings don't have rigid body dynamics, but we can calculate their
        # orientation and position from structural point positions
        if wing.wing_type == REFINE
            # Calculate R_b_w from reference points defining Y and Z directions
            # Y direction: spanwise (from two points across span)
            # Z direction: normal to wing (from two points defining normal, e.g. LE-TE)
            # X = Y × Z (chord direction, ensures right-handed system)

            z_p1, z_p2 = wing.z_ref_points  # Point indices (or vectors to average)
            y_p1, y_p2 = wing.y_ref_points

            # Get positions (with averaging if vectors provided)
            pos_z1 = get_ref_position(pos, z_p1)
            pos_z2 = get_ref_position(pos, z_p2)
            pos_y1 = get_ref_position(pos, y_p1)
            pos_y2 = get_ref_position(pos, y_p2)

            eqs = [
                eqs
                # Build rotation matrix from structural geometry
                # Z direction (normal to wing, normalized)
                vec(R_b_w[wing.idx, :, 3]) .~ sym_normalize(pos_z2 - pos_z1)
                # Y temp direction (not necessarily orthogonal yet)
                # X = Y_temp × Z (chord direction, orthogonal to Z)
                vec(R_b_w[wing.idx, :, 1]) .~ sym_normalize(
                    sym_normalize(pos_y2 - pos_y1) × R_b_w[wing.idx, :, 3]
                )
                # Y = Z × X (ensure orthogonality and right-handed system)
                vec(R_b_w[wing.idx, :, 2]) .~ R_b_w[wing.idx, :, 3] × R_b_w[wing.idx, :, 1]

                # Define wing position from KCU origin point
                # This ensures wing.pos_w moves with structural deformation
                # and VSM panels (plotted at T_b_w=wing.pos_w) stay aligned
                wing_pos[wing.idx, :] ~ pos[:, wing.origin_idx]
                wing_vel[wing.idx, :] ~ vel[:, wing.origin_idx]
                wing_acc[wing.idx, :] ~ acc[:, wing.origin_idx]
            ]

            # Convert rotation matrix to quaternion for REFINE wings
            R_wing = R_b_w[wing.idx, :, :]
            eqs = [
                eqs
                Q_b_w[wing.idx, 1] ~ rotation_matrix_to_quaternion_w(R_wing)
                Q_b_w[wing.idx, 2] ~ rotation_matrix_to_quaternion_x(R_wing)
                Q_b_w[wing.idx, 3] ~ rotation_matrix_to_quaternion_y(R_wing)
                Q_b_w[wing.idx, 4] ~ rotation_matrix_to_quaternion_z(R_wing)
                Q_vel[wing.idx, :] ~ zeros(4)  # REFINE wings have no quaternion velocity (orientation from geometry)
                # Set variables to zero that are only used for QUATERNION rigid body dynamics
                moment_b[wing.idx, :] ~ zeros(3)
                moment_tether_wing[wing.idx, :] ~ zeros(3)
                force_tether_wing[wing.idx, :] ~ zeros(3)
                wing_mass[wing.idx] ~ 0.0  # Mass is distributed to WING points, not centralized
                fix_wing_sphere[wing.idx] ~ false
                # Set intermediate dynamics variables to zero (unused for REFINE)
                wing_acc_b[wing.idx, :] ~ zeros(3)
                α_b_damped[wing.idx, :] ~ zeros(3)
                ω_b_stable[wing.idx, :] ~ zeros(3)
                # aero_force_b and aero_moment_b will be set in linear_vsm_eqs!
            ]
            continue
        end

        I_b = wing.inertia_principal
        axis = sym_normalize(wing_pos[wing.idx, :])
        axis_b = R_b_w[wing.idx, :, :]' * axis
        eqs = [
            eqs
            fix_wing_sphere[wing.idx] ~ get_fix_wing_sphere(psys, wing.idx)
            # Quaternion kinematics: dQ/dt = 0.5 * Ω(ω) * Q
            [D(Q_b_w[wing.idx, i]) ~ Q_vel[wing.idx, i] for i = 1:4]
            [
                Q_vel[wing.idx, i] ~ 0.5 * sum(
                    Ω(ω_b_stable[wing.idx, :])[i, j] * Q_b_w[wing.idx, j] for
                    j = 1:4
                ) for i = 1:4
            ]
            # Constrain angular velocity for spherical joint
            ω_b_stable[wing.idx, :] ~ ifelse.(
                fix_wing == true,
                zeros(3),
                ifelse.(
                    fix_wing_sphere[wing.idx] == true,
                    ω_b[wing.idx, :] - (ω_b[wing.idx, :] ⋅ axis_b) * axis_b,
                    ω_b[wing.idx, :],
                ),
            )
            # Constrain angular acceleration
            D(ω_b[wing.idx, :]) ~ ifelse.(
                fix_wing == true,
                zeros(3),
                ifelse.(
                    fix_wing_sphere[wing.idx] == true,
                    α_b_damped[wing.idx, :] -
                    (α_b_damped[wing.idx, :] ⋅ axis_b) * axis_b,
                    α_b_damped[wing.idx, :],
                ),
            )
            # Apply damping and disturbances to angular acceleration
            α_b_damped[wing.idx, :] ~ [
                α_b[wing.idx, 1],
                α_b[wing.idx, 2] - get_y_damping(psys, wing.idx) * ω_b[wing.idx, 2],
                α_b[wing.idx, 3] + get_z_disturb(psys, wing.idx),
            ]

            # Convert quaternion to rotation matrix
            [
                R_b_w[wing.idx, :, i] ~
                    quaternion_to_rotation_matrix(Q_b_w[wing.idx, :])[:, i] for
                i = 1:3
            ]

            # Euler's rotation equations
            α_b[wing.idx, 1] ~
                (
                    moment_b[wing.idx, 1] +
                    (I_b[2] - I_b[3]) * ω_b[wing.idx, 2] * ω_b[wing.idx, 3]
                ) / I_b[1]
            α_b[wing.idx, 2] ~
                (
                    moment_b[wing.idx, 2] +
                    (I_b[3] - I_b[1]) * ω_b[wing.idx, 3] * ω_b[wing.idx, 1]
                ) / I_b[2]
            α_b[wing.idx, 3] ~
                (
                    moment_b[wing.idx, 3] +
                    (I_b[1] - I_b[2]) * ω_b[wing.idx, 1] * ω_b[wing.idx, 2]
                ) / I_b[3]

            # Total moment in body frame
            moment_tether_wing[wing.idx, :] ~ tether_wing_moment[wing.idx, :]
            moment_b[wing.idx, :] ~
                aero_moment_b[wing.idx, :] +
                R_b_w[wing.idx, :, :]' * moment_tether_wing[wing.idx, :]

            # Translational dynamics (Newton's second law)
            D(wing_pos[wing.idx, :]) ~ ifelse.(
                fix_wing == true,
                zeros(3),
                ifelse.(
                    fix_wing_sphere[wing.idx] == true,
                    (wing_vel[wing.idx, :] ⋅ axis) * axis,
                    wing_vel[wing.idx, :],
                ),
            )
            D(wing_vel[wing.idx, :]) ~ ifelse.(
                fix_wing == true,
                zeros(3),
                ifelse.(
                    fix_wing_sphere[wing.idx] == true,
                    (wing_acc[wing.idx, :] ⋅ axis) * axis,
                    wing_acc[wing.idx, :],
                ),
            )
            wing_mass[wing.idx] ~ get_set_mass(pset)
            force_tether_wing[wing.idx, :] ~ tether_wing_force[wing.idx, :]
            wing_acc[wing.idx, :] ~
                (
                    force_tether_wing[wing.idx, :] +
                    R_b_w[wing.idx, :, :] * aero_force_b[wing.idx, :]
                ) / wing_mass[wing.idx]
        ]
        defaults = [
            defaults
            [Q_b_w[wing.idx, i] => get_Q_b_w(psys, wing.idx)[i] for i = 1:4]
            [ω_b[wing.idx, i] => get_ω_b(psys, wing.idx)[i] for i = 1:3]
            [wing_pos[wing.idx, i] => get_wing_pos_w(psys, wing.idx)[i] for i = 1:3]
            [wing_vel[wing.idx, i] => get_wing_vel_w(psys, wing.idx)[i] for i = 1:3]
        ]
    end

    return eqs, defaults
end

"""
    rotate_v_around_k(v, k, θ)

Rotate vector `v` around axis `k` by angle `θ` using Rodrigues' rotation formula.
"""
function rotate_v_around_k(v, k, θ)
    k = sym_normalize(k)
    v_rot = v * cos(θ) + (k × v) * sin(θ) + k * (k ⋅ v) * (1 - cos(θ))
    return v_rot
end

"""
    calc_R_v_w(wing_pos, e_x)

Calculate the rotation matrix from the view frame (`_v`) to the world frame (`_w`).

The view frame is defined with its z-axis pointing from the origin to the wing,
and its x-axis aligned with the wing's x-axis projected onto the view plane.
"""
function calc_R_v_w(wing_pos, e_x)
    z = sym_normalize(wing_pos)
    y = sym_normalize(z × e_x)
    x = y × z
    return [x y z]
end

"""
    calc_R_t_w(wing_pos)

Calculate the rotation matrix from the local tether frame (`_t`) to the world
frame (`_w`).

The tether frame is a local spherical coordinate system:
- **z-axis**: Aligned with the tether (radial direction).
- **y-axis**: Azimuthal direction, parallel to the XY plane.
- **x-axis**: Elevation direction, tangent to the sphere (`y × z`).
"""
function calc_R_t_w(wing_pos)
    z = sym_normalize(wing_pos)
    if wing_pos[2] ≈ 0.0 && wing_pos[1] ≈ 0.0
        y = [0, 1, 0]
    else
        y = sym_normalize([-wing_pos[2], wing_pos[1], 0])
    end
    x = y × z
    return [x y z]
end

function sym_calc_R_t_w(wing_pos)
    z = sym_normalize(wing_pos)
    y = sym_normalize([-wing_pos[2], wing_pos[1], 0])
    x = y × z
    return [x y z]
end

"""
    calc_heading(R_t_w, R_v_w)

Calculate the heading angle [rad] of the wing. Heading is the rotation angle
between the tether frame and the view frame around the common z-axis.
"""
function calc_heading(R_t_w::AbstractMatrix, R_v_w::AbstractMatrix)
    heading_vec = R_t_w' * R_v_w[:, 1]
    heading = atan(heading_vec[2], heading_vec[1])
    return heading
end

"""
    scalar_eqs!(s, eqs, psys, pset; kwargs...)

Generate equations for derived scalar kinematic quantities useful for control and
analysis.

This includes elevation, azimuth, heading, course, angle of attack, and their time
derivatives, as well as apparent wind calculations.

# Arguments
- `s::SymbolicAWEModel`: The main model object.
- `eqs`, `psys`, `pset`: Accumulating vectors and symbolic parameters.
- `kwargs...`: Symbolic variables for the system's state.

# Returns
- `eqs`: The updated list of system equations.
"""
function scalar_eqs!(
    s, eqs, psys, pset;
    R_b_w, wind_vec_gnd, va_wing_b, wing_pos,
    wing_vel, wing_acc, twist_angle, ω_b, α_b, R_v_w
)
    @unpack wings = s.sys_struct
    @variables begin
        # Body frame axes and apparent wind
        e_x(t)[eachindex(wings), 1:3]
        e_y(t)[eachindex(wings), 1:3]
        e_z(t)[eachindex(wings), 1:3]
        wind_vel_wing(t)[eachindex(wings), 1:3]
        wind_disturb(t)[eachindex(wings), 1:3]
        va_wing(t)[eachindex(wings), 1:3]
        # Ground wind properties
        upwind_dir(t)
        wind_elevation(t)
        wind_scale_gnd(t)
    end
    eqs = [
        eqs
        upwind_dir ~ deg2rad(get_upwind_dir(pset))
        wind_elevation ~ deg2rad(get_wind_elevation(psys))
        wind_scale_gnd ~ get_v_wind(pset)
        wind_vec_gnd ~
            max(wind_scale_gnd, 1e-6) *
            rotate_around_z(rotate_around_x([0, -1, 0], wind_elevation), -upwind_dir)
    ]
    for wing in wings
        # REFINE wings now have R_b_w, wing_pos, wing_vel calculated from structural geometry
        # But they don't have ω_b, α_b (no rotational dynamics)
        if wing.wing_type == REFINE
            # Extract basis vectors from calculated R_b_w
            # Calculate wind at wing centroid position
            eqs = [
                eqs
                e_x[wing.idx, :] ~ R_b_w[wing.idx, :, 1]
                e_y[wing.idx, :] ~ R_b_w[wing.idx, :, 2]
                e_z[wing.idx, :] ~ R_b_w[wing.idx, :, 3]
                wind_vel_wing[wing.idx, :] ~
                    calc_wind_factor(s.am, max(wing_pos[wing.idx, 3], 1.0), pset) *
                    wind_vec_gnd
                wind_disturb[wing.idx, :] ~ get_wind_disturb(psys, wing.idx)
                va_wing[wing.idx, :] ~
                    wind_vel_wing[wing.idx, :] - wing_vel[wing.idx, :] +
                    wind_disturb[wing.idx, :]
                va_wing_b[wing.idx, :] ~ R_b_w[wing.idx, :, :]' * va_wing[wing.idx, :]
            ]
        else
            eqs = [
                eqs
                e_x[wing.idx, :] ~ R_b_w[wing.idx, :, 1]
                e_y[wing.idx, :] ~ R_b_w[wing.idx, :, 2]
                e_z[wing.idx, :] ~ R_b_w[wing.idx, :, 3]
                wind_vel_wing[wing.idx, :] ~
                    calc_wind_factor(s.am, max(wing_pos[wing.idx, 3], 1.0), pset) *
                    wind_vec_gnd
                wind_disturb[wing.idx, :] ~ get_wind_disturb(psys, wing.idx)
                va_wing[wing.idx, :] ~
                    wind_vel_wing[wing.idx, :] - wing_vel[wing.idx, :] +
                    wind_disturb[wing.idx, :]
                va_wing_b[wing.idx, :] ~ R_b_w[wing.idx, :, :]' * va_wing[wing.idx, :]
            ]
        end
    end
    @variables begin
        # Kinematic quantities
        heading(t)[eachindex(wings)]
        turn_rate(t)[eachindex(wings), 1:3]
        turn_acc(t)[eachindex(wings), 1:3]
        azimuth(t)[eachindex(wings)]
        azimuth_vel(t)[eachindex(wings)]
        azimuth_acc(t)[eachindex(wings)]
        elevation(t)[eachindex(wings)]
        elevation_vel(t)[eachindex(wings)]
        elevation_acc(t)[eachindex(wings)]
        course(t)[eachindex(wings)]
        angle_of_attack(t)[eachindex(wings)]
        R_t_w(t)[eachindex(wings), 1:3, 1:3]
        distance(t)[eachindex(wings)]
        distance_vel(t)[eachindex(wings)]
        distance_acc(t)[eachindex(wings)]
    end

    for wing in wings
        # REFINE wings have wing_pos, wing_vel, wing_acc, R_b_w from structural geometry
        # But they don't have ω_b, α_b (no rotational dynamics)
        if wing.wing_type == REFINE
            x, y, z = wing_pos[wing.idx, :]
            eqs = [
                eqs
                vec(R_v_w[wing.idx, :, :]) .~
                    vec(calc_R_v_w(wing_pos[wing.idx, :], e_x[wing.idx, :]))
                vec(R_t_w[wing.idx, :, :]) .~
                    vec(sym_calc_R_t_w(wing_pos[wing.idx, :]))
                heading[wing.idx] ~
                    calc_heading(R_t_w[wing.idx, :, :], R_v_w[wing.idx, :, :])
                # Rotational quantities are zero (no rigid body rotation)
                ω_b[wing.idx, :] ~ zeros(3)
                α_b[wing.idx, :] ~ zeros(3)
                turn_rate[wing.idx, :] ~ zeros(3)
                turn_acc[wing.idx, :] ~ zeros(3)
                # Translational kinematics use actual centroid motion
                distance[wing.idx] ~ norm(wing_pos[wing.idx, :])
                distance_vel[wing.idx] ~
                    wing_vel[wing.idx, :] ⋅ R_v_w[wing.idx, :, 3]
                distance_acc[wing.idx] ~
                    wing_acc[wing.idx, :] ⋅ R_v_w[wing.idx, :, 3]
                elevation[wing.idx] ~ KiteUtils.calc_elevation(wing_pos[wing.idx, :])
                elevation_vel[wing.idx] ~
                    dot(wing_vel[wing.idx, :], -R_t_w[wing.idx, :, 1]) /
                    distance[wing.idx]
                elevation_acc[wing.idx] ~
                    dot(wing_acc[wing.idx, :], -R_t_w[wing.idx, :, 1]) /
                    distance[wing.idx]
                azimuth[wing.idx] ~ KiteUtils.azimuth_east(wing_pos[wing.idx, :])
                azimuth_vel[wing.idx] ~
                    dot(wing_vel[wing.idx, :], -R_t_w[wing.idx, :, 2]) / norm([x, y])
                azimuth_acc[wing.idx] ~
                    dot(wing_acc[wing.idx, :], -R_t_w[wing.idx, :, 2]) / norm([x, y])
                course[wing.idx] ~
                    atan(-azimuth_vel[wing.idx], elevation_vel[wing.idx])
                # Angle of attack from apparent wind
                angle_of_attack[wing.idx] ~
                    calc_angle_of_attack(va_wing_b[wing.idx, :])
            ]
        else
            x, y, z = wing_pos[wing.idx, :]
            has_groups = !isempty(wing.group_idxs)
            if has_groups
                half_len = wing.group_idxs[1] + length(wing.group_idxs) ÷ 2 - 1
            end

            eqs = [
                eqs
                vec(R_v_w[wing.idx, :, :]) .~
                    vec(calc_R_v_w(wing_pos[wing.idx, :], e_x[wing.idx, :]))
                vec(R_t_w[wing.idx, :, :]) .~
                    vec(sym_calc_R_t_w(wing_pos[wing.idx, :]))
                heading[wing.idx] ~
                    calc_heading(R_t_w[wing.idx, :, :], R_v_w[wing.idx, :, :])
                turn_rate[wing.idx, :] ~
                    R_v_w[wing.idx, :, :]' *
                    (R_b_w[wing.idx, :, :] * ω_b[wing.idx, :])
                turn_acc[wing.idx, :] ~
                    R_v_w[wing.idx, :, :]' *
                    (R_b_w[wing.idx, :, :] * α_b[wing.idx, :])
                distance[wing.idx] ~ norm(wing_pos[wing.idx, :])
                distance_vel[wing.idx] ~
                    wing_vel[wing.idx, :] ⋅ R_v_w[wing.idx, :, 3]
                distance_acc[wing.idx] ~
                    wing_acc[wing.idx, :] ⋅ R_v_w[wing.idx, :, 3]

                elevation[wing.idx] ~ KiteUtils.calc_elevation(wing_pos[wing.idx, :])
                elevation_vel[wing.idx] ~
                    dot(wing_vel[wing.idx, :], -R_t_w[wing.idx, :, 1]) /
                    distance[wing.idx]
                elevation_acc[wing.idx] ~
                    dot(wing_acc[wing.idx, :], -R_t_w[wing.idx, :, 1]) /
                    distance[wing.idx]
                azimuth[wing.idx] ~ KiteUtils.azimuth_east(wing_pos[wing.idx, :])
                azimuth_vel[wing.idx] ~
                    dot(wing_vel[wing.idx, :], -R_t_w[wing.idx, :, 2]) / norm([x, y])
                azimuth_acc[wing.idx] ~
                    dot(wing_acc[wing.idx, :], -R_t_w[wing.idx, :, 2]) / norm([x, y])
                course[wing.idx] ~
                    atan(-azimuth_vel[wing.idx], elevation_vel[wing.idx])

                angle_of_attack[wing.idx] ~
                    calc_angle_of_attack(va_wing_b[wing.idx, :]) +
                    (has_groups ? 0.5 * twist_angle[half_len] + 0.5 * twist_angle[half_len + 1] : 0)
            ]
        end
    end
    return eqs
end

"""
    Base.getindex(x::ModelingToolkit.Symbolics.SymArray, idxs::Vector{Int16})

Extend `Base.getindex` to allow indexing a symbolic array with a vector of
integer indices, which is not natively supported by ModelingToolkit.
"""
function Base.getindex(x::ModelingToolkit.Symbolics.SymArray, idxs::Vector{Int16})
    Num[Base.getindex(x, idx) for idx in idxs]
end

"""
    linear_vsm_eqs!(s, eqs, guesses, psys; kwargs...)

Generate linearized aerodynamic equations using the Vortex Step Method (VSM).

This function approximates the aerodynamic forces and moments using a first-order
Taylor expansion around a pre-calculated operating point. The Jacobian of the
aerodynamic forces w.r.t. the state variables is provided via symbolic parameters.

# Arguments
- `s::SymbolicAWEModel`: The main model object.
- `eqs`, `guesses`, `psys`: Accumulating vectors and symbolic parameters.
- `kwargs...`: Symbolic variables for aerodynamic and state quantities.

# Returns
- A tuple `(eqs, guesses)` containing the updated equation and guess lists.
"""
function linear_vsm_eqs!(
    s, eqs, guesses, psys;
    aero_force_b, aero_moment_b, group_aero_moment,
    twist_angle, va_wing_b, wing_pos, ω_b, R_v_w,
    aero_force_point_b=nothing  # REFINE-specific parameter
)
    @unpack groups, wings, points = s.sys_struct
    if length(wings) == 0
        return eqs, guesses
    end

    # Check for REFINE and QUATERNION wings
    has_refine = any(w.wing_type == REFINE for w in wings)
    has_quaternion = any(w.wing_type == QUATERNION for w in wings)

    # Compute maximum dimensions across all wings
    # ny = dimension of VSM input state for QUATERNION wings: [va(3), ω(3), twist(n_groups)]
    # Only compute if we have QUATERNION wings
    ny_quaternion = has_quaternion ? 3 + length(wings[1].group_idxs) + 3 : 0

    # nx = dimension of VSM output (forces/moments for QUATERNION, panel forces for REFINE)
    nx_values = Int[]
    for wing in wings
        if wing.wing_type == REFINE
            # REFINE: Store per-panel forces [fx_1, fy_1, fz_1, fx_2, ...]
            n_panels = length(wing.vsm_aero.panels)
            push!(nx_values, 3 * n_panels)
        else
            # QUATERNION: Store [aero_force(3), aero_moment(3), group_moments(n_groups)]
            push!(nx_values, 3 + 3 + length(wing.group_idxs))
        end
    end
    nx_max = maximum(nx_values)

    # Declare VSM output forces for all wings (used by both REFINE and QUATERNION)
    @variables begin
        vsm_output_force_prev(t)[eachindex(wings), 1:nx_max]  # Panel forces (REFINE) or [forces, moments] (QUATERNION)
    end

    # Declare linearization variables only for QUATERNION wings
    if has_quaternion
        @variables begin
            # VSM linearization variables (QUATERNION wings only)
            # Linearization: F(state) ≈ F(state₀) + ∂F/∂state|_{state₀} * (state - state₀)
            # where state = [va_wing_b, ω_b, twist_angles] and F = [forces, moments]
            vsm_input_state(t)[eachindex(wings), 1:ny_quaternion]        # Current input state
            vsm_input_state_delta(t)[eachindex(wings), 1:ny_quaternion]  # Δstate for linearization
            vsm_input_state_prev(t)[eachindex(wings), 1:ny_quaternion]   # State at linearization point
            force_jacobian(t)[eachindex(wings), 1:nx_max, 1:ny_quaternion]  # ∂F/∂state Jacobian matrix

            # Aerodynamic quantities (QUATERNION wings only)
            q_inf(t)[eachindex(wings)]                      # Dynamic pressure
            no_scale_aero_force_b(t)[eachindex(wings), 1:3]  # Unscaled force before q_inf scaling
        end
    end

    for wing in wings
        if wing.wing_type == REFINE
            # ==================== REFINE WING: Direct Panel Forces ====================
            # REFINE wings apply VSM panel forces directly to structural points without
            # linearization. Each VSM panel contributes to nearby structural points via
            # inverse distance weighting.
            #
            # Storage: vsm_output_force_prev[wing.idx, 1:nx_refine] stores per-panel forces
            #          as [fx_1, fy_1, fz_1, fx_2, fy_2, fz_2, ...]
            #
            # NO linearization: Forces come from full nonlinear VSM solve each timestep.
            # (Linearization would be prohibitively expensive with many structural DOFs)

            n_panels = length(wing.vsm_aero.panels)
            nx_refine = 3 * n_panels

            # Retrieve current panel forces from VSM solution
            eqs = [
                eqs
                [vsm_output_force_prev[wing.idx, ix] ~ get_vsm_x(psys, wing.idx, ix) for ix = 1:nx_refine]
            ]

            # Lump panel forces to structural WING points
            wing_points = [p for p in points if p.type == WING && p.wing_idx == wing.idx]

            for point in wing_points
                # Weighted sum of forces from nearby panels (weights precomputed in wing.point_to_panels)
                lumped_force_b = zeros(Num, 3)

                for (panel_idx, weight) in wing.point_to_panels[point.idx]
                    panel_force_start = 3*(panel_idx-1) + 1
                    for i in 1:3
                        force_idx = panel_force_start + i - 1
                        lumped_force_b[i] += weight * vsm_output_force_prev[wing.idx, force_idx]
                    end
                end

                eqs = [
                    eqs
                    aero_force_point_b[point.idx, :] ~ lumped_force_b
                ]
            end

            # Set centralized aero variables to zero (REFINE wings distribute forces to points)
            eqs = [
                eqs
                aero_force_b[wing.idx, :] ~ zeros(3)
                aero_moment_b[wing.idx, :] ~ zeros(3)
            ]

        else
            # ==================== QUATERNION WING: Linearized Aerodynamics ====================
            # QUATERNION wings use first-order Taylor expansion to approximate aerodynamic forces:
            #
            #   F(state) ≈ F(state₀) + ∂F/∂state|_{state₀} * (state - state₀)
            #
            # Where:
            #   state = [va_wing_b(3), ω_b(3), twist_angle(n_groups)]  ← VSM inputs
            #   F = [aero_force(3), aero_moment(3), group_moments(n_groups)]  ← VSM outputs
            #   state₀ = previous linearization point (updated periodically by VSM solve)
            #   ∂F/∂state = Jacobian matrix (computed by VSM via finite differences)
            #
            # This allows efficient force updates without re-solving the full VSM system every timestep.

            area = wing.vsm_aero.projected_area
            force_b = no_scale_aero_force_b[wing.idx, :]
            wind_direction_b = sym_normalize(va_wing_b[wing.idx, :])
            drag_force_b = (force_b ⋅ wind_direction_b) * wind_direction_b

            # Dimensions for this wing
            nx_quat = 3 + 3 + length(wing.group_idxs)

            eqs = [
                eqs
                # Dynamic pressure for force scaling
                q_inf[wing.idx] ~
                    0.5 * calc_rho(s.am, wing_pos[wing.idx, 3]) *
                    norm(va_wing_b[wing.idx, :])^2

                # Load linearization data from VSM (state₀, F₀, ∂F/∂state)
                [vsm_input_state_prev[wing.idx, iy] ~ get_vsm_y(psys, wing.idx, iy) for iy = 1:ny_quaternion]
                [vsm_output_force_prev[wing.idx, ix] ~ get_vsm_x(psys, wing.idx, ix) for ix = 1:nx_quat]
                [
                    force_jacobian[wing.idx, ix, iy] ~
                        get_vsm_jac(psys, wing.idx, ix, iy) for ix = 1:nx_quat for
                    iy = 1:ny_quaternion
                ]

                # Current input state for linearization
                vsm_input_state[wing.idx, :] ~ [
                    va_wing_b[wing.idx, :]        # Apparent wind velocity in body frame
                    ω_b[wing.idx, :]              # Angular velocity in body frame
                    twist_angle[wing.group_idxs]  # Wing group twist angles
                ]

                # State deviation from linearization point: Δstate = state - state₀
                vsm_input_state_delta[wing.idx, :] ~ vsm_input_state[wing.idx, :] - vsm_input_state_prev[wing.idx, :]

                # ===== LINEARIZED FORCE CALCULATION =====
                # F ≈ F₀ + ∂F/∂state * Δstate
                # Then scale by dynamic pressure and wing area
                [
                    force_b
                    aero_moment_b[wing.idx, :]
                    group_aero_moment[wing.group_idxs]
                ] ~
                    q_inf[wing.idx] *
                    area *
                    (vsm_output_force_prev[wing.idx, 1:nx_quat] +
                     force_jacobian[wing.idx, 1:nx_quat, :] * vsm_input_state_delta[wing.idx, :])

                # Apply additional drag correction factor
                aero_force_b[wing.idx, :] ~
                    force_b + drag_force_b * (get_drag_frac(psys, wing.idx) - 1)
            ]

            if s.set.quasi_static
                guesses = [
                    guesses
                    [vsm_input_state[wing.idx, iy] => get_vsm_y(psys, wing.idx, iy) for iy = 1:ny_quaternion]
                ]
            end
        end
    end
    return eqs, guesses
end

"""
    create_sys!(s::SymbolicAWEModel, system::SystemStructure; prn=true)

Create the full `ModelingToolkit.ODESystem` for the AWE model.

This is the main top-level function that orchestrates the generation of the entire
set of differential-algebraic equations (DAEs). It calls specialized sub-functions
to build the equations for each part of the system (forces, wing dynamics, scalar
kinematics, linearized aerodynamics) and assembles them into a single `System`.

# Arguments
- `s::SymbolicAWEModel`: The main model object to be populated with the system.
- `system::SystemStructure`: The physical structure definition.
- `prn::Bool=true`: If true, print progress information during system creation.

# Returns
- `set_values`: The symbolic variable representing the control inputs (winch torques).
"""
function create_sys!(s::SymbolicAWEModel, system::SystemStructure; prn = true)
    eqs = []
    defaults = Pair{Num,Any}[]
    guesses = Pair{Num,Any}[]

    @unpack wings, groups, winches, points = system

    # Validation for REFINE wings
    for wing in wings
        if wing.wing_type == REFINE
            # REFINE wings cannot have groups
            @assert length(wing.group_idxs) == 0 "REFINE wing $(wing.idx) cannot have groups"
            @assert !isnothing(wing.point_to_panels) "REFINE wing $(wing.idx) missing point_to_panels mapping"

            # Verify all WING points for this wing are in the mapping
            wing_point_idxs = [p.idx for p in points if p.type == WING && p.wing_idx == wing.idx]
            for point_idx in wing_point_idxs
                @assert haskey(wing.point_to_panels, point_idx) "REFINE wing $(wing.idx) missing mapping for point $(point_idx)"
            end

            prn && println("✓ REFINE wing $(wing.idx) validated: $(length(wing_point_idxs)) points, $(length(wing.vsm_aero.panels)) panels")
        end
    end

    @parameters begin
        psys::SystemStructure = system
        pset::Settings = s.set
        fix_wing = false
    end
    @variables begin
        # Control inputs
        set_values(t)[eachindex(winches)] = zeros(length(winches))
        # Wing rigid body states
        wing_pos(t)[eachindex(wings), 1:3]
        wing_vel(t)[eachindex(wings), 1:3]
        wing_acc(t)[eachindex(wings), 1:3]
        ω_b(t)[eachindex(wings), 1:3]
        α_b(t)[eachindex(wings), 1:3]
        # Wing rotation matrices
        R_b_w(t)[eachindex(wings), 1:3, 1:3]
        R_v_w(t)[eachindex(wings), 1:3, 1:3]
        # Aerodynamic forces and moments
        aero_force_b(t)[eachindex(wings), 1:3]
        aero_moment_b(t)[eachindex(wings), 1:3]
        group_aero_moment(t)[eachindex(groups)]
        # Wing deformation states
        twist_angle(t)[eachindex(groups)]
        twist_ω(t)[eachindex(groups)]
        # Wind and apparent velocity
        wind_vec_gnd(t)[1:3]
        va_wing_b(t)[eachindex(wings), 1:3]
    end

    # Build tether and bridle system equations
    eqs, defaults, guesses, tether_wing_force, tether_wing_moment, aero_force_point_b, pos, vel, acc = force_eqs!(
        s, system, psys, pset, eqs, defaults, guesses;
        R_b_w, wing_pos, wing_vel, wind_vec_gnd, group_aero_moment,
        twist_angle, twist_ω, set_values, fix_wing
    )
    # Build linearized aerodynamic equations
    # Check if we need aero_force_point_b (for REFINE wings)
    has_refine_wings = any(wing.wing_type == REFINE for wing in wings)
    if has_refine_wings
        # aero_force_point_b was already defined in force_eqs!
        eqs, guesses = linear_vsm_eqs!(
            s, eqs, guesses, psys;
            aero_force_b, R_v_w, aero_moment_b, group_aero_moment,
            twist_angle, va_wing_b, wing_pos, ω_b,
            aero_force_point_b=aero_force_point_b
        )
    else
        eqs, guesses = linear_vsm_eqs!(
            s, eqs, guesses, psys;
            aero_force_b, R_v_w, aero_moment_b, group_aero_moment,
            twist_angle, va_wing_b, wing_pos, ω_b
        )
    end
    # Build wing rigid body dynamics equations
    eqs, defaults = wing_eqs!(
        s, eqs, psys, pset, defaults;
        tether_wing_force, tether_wing_moment, aero_force_b, aero_moment_b,
        ω_b, α_b, R_b_w, wing_pos, wing_vel, wing_acc, fix_wing,
        pos, vel, acc
    )
    # Build scalar kinematic and apparent wind equations
    eqs = scalar_eqs!(
        s, eqs, psys, pset;
        R_b_w, wind_vec_gnd, va_wing_b, wing_pos, wing_vel, wing_acc,
        twist_angle, ω_b, α_b, R_v_w
    )

    eqs = Symbolics.scalarize.(reduce(vcat, Symbolics.scalarize.(eqs)))

    time = @elapsed @named sys = System(eqs, t)
    prn && println("\tCreated System in $time seconds.")

    defaults = [
        defaults
        [
            set_values[winch.idx] => get_set_value(psys, winch.idx) for
            winch in winches
        ]
    ]

    s.defaults = defaults
    s.guesses = guesses
    s.full_sys = sys
    return set_values
end
