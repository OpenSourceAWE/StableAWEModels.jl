# Copyright (c) 2025 Bart van de Lint
# SPDX-License-Identifier: MPL-2.0

# Scalar kinematic equation generation

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
                    calc_wind_factor(s.am, wing_pos[wing.idx, 1], wing_pos[wing.idx, 2],
                                     wing_pos[wing.idx, 3], pset) * wind_vec_gnd
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
                    calc_wind_factor(s.am, wing_pos[wing.idx, 1], wing_pos[wing.idx, 2],
                                     wing_pos[wing.idx, 3], pset) * wind_vec_gnd
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
            # Calculate heading using wind-perpendicular frame
            # Normalize wind direction
            wind_norm = sym_normalize(wind_vel_wing[wing.idx, :])
            # Project -e_x onto plane perpendicular to wind
            minus_e_x = -e_x[wing.idx, :]
            proj_on_wind = (minus_e_x ⋅ wind_norm) * wind_norm
            e_x_perp = minus_e_x - proj_on_wind
            # Heading is angle in wind-perpendicular plane
            # x-component: perpendicular to both wind and z
            wind_cross_z = [wind_norm[2], -wind_norm[1], 0]
            heading_x = e_x_perp ⋅ wind_cross_z
            # z-component: world z-axis
            heading_z = e_x_perp[3]
            # Calculate course using same wind-perpendicular projection
            proj_vel_on_wind = (wing_vel[wing.idx, :] ⋅ wind_norm) * wind_norm
            vel_perp = wing_vel[wing.idx, :] - proj_vel_on_wind
            course_x = vel_perp ⋅ wind_cross_z
            course_z = vel_perp[3]
            eqs = [
                eqs
                vec(R_v_w[wing.idx, :, :]) .~
                    vec(calc_R_v_w(wing_pos[wing.idx, :], e_x[wing.idx, :]))
                vec(R_t_w[wing.idx, :, :]) .~
                    vec(sym_calc_R_t_w(wing_pos[wing.idx, :]))
                heading[wing.idx] ~ atan(heading_x, heading_z)
                # Rotational quantities are zero (no rigid body rotation)
                ω_b[wing.idx, :] ~ zeros(3)
                α_b[wing.idx, :] ~ zeros(3)
                turn_rate[wing.idx, :] ~ zeros(3)
                turn_acc[wing.idx, :] ~ zeros(3)
                # Translational kinematics use actual centroid motion
                distance[wing.idx] ~ norm(wing_pos[wing.idx, :])
                distance_vel[wing.idx] ~
                    wing_vel[wing.idx, :] ⋅ R_t_w[wing.idx, :, 3]
                distance_acc[wing.idx] ~
                    wing_acc[wing.idx, :] ⋅ R_t_w[wing.idx, :, 3]
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
                # Course is the direction of velocity in the wind-perpendicular plane
                course[wing.idx] ~ atan(course_x, course_z)
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

            # Calculate heading using wind-perpendicular frame
            # Normalize wind direction
            wind_norm = sym_normalize(wind_vel_wing[wing.idx, :])
            # Project -e_x onto plane perpendicular to wind
            minus_e_x = -e_x[wing.idx, :]
            proj_on_wind = (minus_e_x ⋅ wind_norm) * wind_norm
            e_x_perp = minus_e_x - proj_on_wind
            # Heading is angle in wind-perpendicular plane
            # x-component: perpendicular to both wind and z
            wind_cross_z = [wind_norm[2], -wind_norm[1], 0]
            heading_x = e_x_perp ⋅ wind_cross_z
            # z-component: world z-axis
            heading_z = e_x_perp[3]
            # Calculate course using same wind-perpendicular projection
            proj_vel_on_wind = (wing_vel[wing.idx, :] ⋅ wind_norm) * wind_norm
            vel_perp = wing_vel[wing.idx, :] - proj_vel_on_wind
            course_x = vel_perp ⋅ wind_cross_z
            course_z = vel_perp[3]
            eqs = [
                eqs
                vec(R_v_w[wing.idx, :, :]) .~
                    vec(calc_R_v_w(wing_pos[wing.idx, :], e_x[wing.idx, :]))
                vec(R_t_w[wing.idx, :, :]) .~
                    vec(sym_calc_R_t_w(wing_pos[wing.idx, :]))
                heading[wing.idx] ~ atan(heading_x, heading_z)
                turn_rate[wing.idx, :] ~
                    R_v_w[wing.idx, :, :]' *
                    (R_b_w[wing.idx, :, :] * ω_b[wing.idx, :])
                turn_acc[wing.idx, :] ~
                    R_v_w[wing.idx, :, :]' *
                    (R_b_w[wing.idx, :, :] * α_b[wing.idx, :])
                distance[wing.idx] ~ norm(wing_pos[wing.idx, :])
                distance_vel[wing.idx] ~
                    wing_vel[wing.idx, :] ⋅ R_t_w[wing.idx, :, 3]
                distance_acc[wing.idx] ~
                    wing_acc[wing.idx, :] ⋅ R_t_w[wing.idx, :, 3]

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
                # Course is the direction of velocity in the wind-perpendicular plane
                course[wing.idx] ~ atan(course_x, course_z)

                angle_of_attack[wing.idx] ~
                    calc_angle_of_attack(va_wing_b[wing.idx, :]) +
                    (has_groups ? 0.5 * twist_angle[half_len] + 0.5 * twist_angle[half_len + 1] : 0)
            ]
        end
    end
    return eqs
end
