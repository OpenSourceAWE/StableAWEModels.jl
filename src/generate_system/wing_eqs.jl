# Copyright (c) 2025 Bart van de Lint
# SPDX-License-Identifier: MPL-2.0

# Wing rigid body dynamics equation generation

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
        wing_acc_b(t)[1:3, eachindex(wings)]
        α_b_damped(t)[1:3, eachindex(wings)]
        ω_b_stable(t)[1:3, eachindex(wings)]
        # Orientation states (all wings)
        Q_b_w(t)[1:4, eachindex(wings)]
        Q_vel(t)[1:4, eachindex(wings)]
        # Forces, moments, and other properties
        moment_b(t)[1:3, eachindex(wings)]
        moment_tether_wing(t)[1:3, eachindex(wings)]
        force_tether_wing(t)[1:3, eachindex(wings)]
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
    get_ref_position(pos, ref::Int64) = pos[:, ref]
    function get_ref_position(pos, refs::Vector{Int64})
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
            #
            # NOTE: This symbolic implementation must match calc_refine_wing_frame() in
            # system_structure.jl. Both calculate R_b_w from the same reference points.
            # If you modify this logic, update calc_refine_wing_frame() as well!

            z_p1, z_p2 = wing.z_ref_points  # Point indices (or vectors to average)
            y_p1, y_p2 = wing.y_ref_points

            # Get positions (with averaging if vectors provided)
            # Equivalent to get_ref_position_from_points() in system_structure.jl
            pos_z1 = get_ref_position(pos, z_p1)
            pos_z2 = get_ref_position(pos, z_p2)
            pos_y1 = get_ref_position(pos, y_p1)
            pos_y2 = get_ref_position(pos, y_p2)

            eqs = [
                eqs
                # Build rotation matrix from structural geometry
                # (Same algorithm as calc_refine_wing_frame)
                # Z direction (normal to wing, normalized)
                R_b_w[:, 3, wing.idx] ~ sym_normalize(pos_z2 - pos_z1)
                # Y temp direction (not necessarily orthogonal yet)
                # X = Y_temp × Z (chord direction, orthogonal to Z)
                R_b_w[:, 1, wing.idx] ~ sym_normalize(
                    sym_normalize(pos_y2 - pos_y1) × R_b_w[:, 3, wing.idx]
                )
                # Y = Z × X (ensure orthogonality and right-handed system)
                R_b_w[:, 2, wing.idx] ~ R_b_w[:, 3, wing.idx] × R_b_w[:, 1, wing.idx]

                # Define wing position from KCU origin point
                # This ensures wing.pos_w moves with structural deformation
                # and VSM panels (plotted at T_b_w=wing.pos_w) stay aligned
                wing_pos[:, wing.idx] ~ pos[:, wing.origin_idx]
                wing_vel[:, wing.idx] ~ vel[:, wing.origin_idx]
                wing_acc[:, wing.idx] ~ acc[:, wing.origin_idx]
            ]

            # Convert rotation matrix to quaternion for REFINE wings
            R_wing = R_b_w[:, :, wing.idx]
            eqs = [
                eqs
                Q_b_w[1, wing.idx] ~ rotation_matrix_to_quaternion_w(R_wing)
                Q_b_w[2, wing.idx] ~ rotation_matrix_to_quaternion_x(R_wing)
                Q_b_w[3, wing.idx] ~ rotation_matrix_to_quaternion_y(R_wing)
                Q_b_w[4, wing.idx] ~ rotation_matrix_to_quaternion_z(R_wing)
                Q_vel[:, wing.idx] ~ zeros(4)  # REFINE wings have no quaternion velocity (orientation from geometry)
                # Set variables to zero that are only used for QUATERNION rigid body dynamics
                moment_b[:, wing.idx] ~ zeros(3)
                moment_tether_wing[:, wing.idx] ~ zeros(3)
                force_tether_wing[:, wing.idx] ~ zeros(3)
                wing_mass[wing.idx] ~ 0.0  # Mass is distributed to WING points, not centralized
                fix_wing_sphere[wing.idx] ~ false
                # Set intermediate dynamics variables to zero (unused for REFINE)
                wing_acc_b[:, wing.idx] ~ zeros(3)
                α_b_damped[:, wing.idx] ~ zeros(3)
                ω_b_stable[:, wing.idx] ~ zeros(3)
                # aero_force_b and aero_moment_b will be set in linear_vsm_eqs!
            ]
            continue
        end

        I_b = wing.inertia_principal
        axis = sym_normalize(wing_pos[:, wing.idx])
        axis_b = R_b_w[:, :, wing.idx]' * axis
        eqs = [
            eqs
            fix_wing_sphere[wing.idx] ~ get_fix_wing_sphere(psys, wing.idx)
            # Quaternion kinematics: dQ/dt = 0.5 * Ω(ω) * Q
            [D(Q_b_w[i, wing.idx]) ~ Q_vel[i, wing.idx] for i = 1:4]
            [
                Q_vel[i, wing.idx] ~ 0.5 * sum(
                    Ω(ω_b_stable[:, wing.idx])[i, j] * Q_b_w[j, wing.idx] for
                    j = 1:4
                ) for i = 1:4
            ]
            # Constrain angular velocity for spherical joint
            ω_b_stable[:, wing.idx] ~ ifelse.(
                fix_wing == true,
                zeros(3),
                ifelse.(
                    fix_wing_sphere[wing.idx] == true,
                    ω_b[:, wing.idx] - (ω_b[:, wing.idx] ⋅ axis_b) * axis_b,
                    ω_b[:, wing.idx],
                ),
            )
            # Constrain angular acceleration
            D(ω_b[:, wing.idx]) ~ ifelse.(
                fix_wing == true,
                zeros(3),
                ifelse.(
                    fix_wing_sphere[wing.idx] == true,
                    α_b_damped[:, wing.idx] -
                    (α_b_damped[:, wing.idx] ⋅ axis_b) * axis_b,
                    α_b_damped[:, wing.idx],
                ),
            )
            # Apply damping and disturbances to angular acceleration
            α_b_damped[:, wing.idx] ~ [
                α_b[1, wing.idx],
                α_b[2, wing.idx] - get_y_damping(psys, wing.idx) * ω_b[2, wing.idx],
                α_b[3, wing.idx] + get_z_disturb(psys, wing.idx),
            ]

            # Convert quaternion to rotation matrix
            [
                R_b_w[:, i, wing.idx] ~
                    quaternion_to_rotation_matrix(Q_b_w[:, wing.idx])[:, i] for
                i = 1:3
            ]

            # Euler's rotation equations
            α_b[1, wing.idx] ~
                (
                    moment_b[1, wing.idx] +
                    (I_b[2] - I_b[3]) * ω_b[2, wing.idx] * ω_b[3, wing.idx]
                ) / I_b[1]
            α_b[2, wing.idx] ~
                (
                    moment_b[2, wing.idx] +
                    (I_b[3] - I_b[1]) * ω_b[3, wing.idx] * ω_b[1, wing.idx]
                ) / I_b[2]
            α_b[3, wing.idx] ~
                (
                    moment_b[3, wing.idx] +
                    (I_b[1] - I_b[2]) * ω_b[1, wing.idx] * ω_b[2, wing.idx]
                ) / I_b[3]

            # Total moment in body frame
            moment_tether_wing[:, wing.idx] ~ tether_wing_moment[:, wing.idx]
            moment_b[:, wing.idx] ~
                aero_moment_b[:, wing.idx] +
                R_b_w[:, :, wing.idx]' * moment_tether_wing[:, wing.idx]

            # Translational dynamics (Newton's second law)
            D(wing_pos[:, wing.idx]) ~ ifelse.(
                fix_wing == true,
                zeros(3),
                ifelse.(
                    fix_wing_sphere[wing.idx] == true,
                    (wing_vel[:, wing.idx] ⋅ axis) * axis,
                    wing_vel[:, wing.idx],
                ),
            )
            D(wing_vel[:, wing.idx]) ~ ifelse.(
                fix_wing == true,
                zeros(3),
                ifelse.(
                    fix_wing_sphere[wing.idx] == true,
                    (wing_acc[:, wing.idx] ⋅ axis) * axis,
                    wing_acc[:, wing.idx],
                ),
            )
            wing_mass[wing.idx] ~ get_set_mass(pset)
            force_tether_wing[:, wing.idx] ~ tether_wing_force[:, wing.idx]
            wing_acc[:, wing.idx] ~
                (
                    force_tether_wing[:, wing.idx] +
                    R_b_w[:, :, wing.idx] * aero_force_b[:, wing.idx]
                ) / wing_mass[wing.idx]
        ]
        defaults = [
            defaults
            [Q_b_w[i, wing.idx] => get_Q_b_w(psys, wing.idx)[i] for i = 1:4]
            [ω_b[i, wing.idx] => get_ω_b(psys, wing.idx)[i] for i = 1:3]
            [wing_pos[i, wing.idx] => get_wing_pos_w(psys, wing.idx)[i] for i = 1:3]
            [wing_vel[i, wing.idx] => get_wing_vel_w(psys, wing.idx)[i] for i = 1:3]
        ]
    end

    return eqs, defaults
end
