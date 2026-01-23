# Copyright (c) 2025 Bart van de Lint
# SPDX-License-Identifier: MPL-2.0

# Symbolic accessor functions for ModelingToolkit integration
#
# The following `get_*` functions are registered for use within ModelingToolkit.jl's
# symbolic context. They act as symbolic "shims" or placeholders that allow the
# equation generation logic to access fields from the `SystemStructure` (`psys`) and
# `Settings` (`pset`) parameter objects during the symbolic construction of the model.
# This avoids hard-coding parameters into the equations, allowing them to be changed
# at simulation time.

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
get_va_b(sys::SystemStructure, idx::Int16) = sys.points[idx].va_b
@register_array_symbolic get_va_b(sys::SystemStructure, idx::Int16) begin
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
get_group_y_airf(sys::SystemStructure, idx::Int16) = sys.groups[idx].y_airf
@register_array_symbolic get_group_y_airf(sys::SystemStructure, idx::Int16) begin
    size = (3,)
    eltype = SimFloat
end
get_group_chord(sys::SystemStructure, idx::Int16) = sys.groups[idx].chord
@register_array_symbolic get_group_chord(sys::SystemStructure, idx::Int16) begin
    size = (3,)
    eltype = SimFloat
end
get_group_le_pos(sys::SystemStructure, idx::Int16) = sys.groups[idx].le_pos
@register_array_symbolic get_group_le_pos(sys::SystemStructure, idx::Int16) begin
    size = (3,)
    eltype = SimFloat
end
get_extra_mass(sys::SystemStructure, idx::Int16) = sys.points[idx].extra_mass
@register_symbolic get_extra_mass(sys::SystemStructure, idx::Int16)
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
@register_array_symbolic get_body_frame_damping(sys::SystemStructure, idx::Int16) begin
    size = (3,)
    eltype = SimFloat
end
get_world_frame_damping(sys::SystemStructure, idx::Int16) = sys.points[idx].world_frame_damping
@register_array_symbolic get_world_frame_damping(sys::SystemStructure, idx::Int16) begin
    size = (3,)
    eltype = SimFloat
end
get_point_area(sys::SystemStructure, idx::Int16) = sys.points[idx].area
@register_symbolic get_point_area(sys::SystemStructure, idx::Int16)
get_point_drag_coeff(sys::SystemStructure, idx::Int16) = sys.points[idx].drag_coeff
@register_symbolic get_point_drag_coeff(sys::SystemStructure, idx::Int16)
get_point_aero_force(sys::SystemStructure, idx::Int16, component::Int) = sys.points[idx].aero_force_b[component]
@register_symbolic get_point_aero_force(sys::SystemStructure, idx::Int16, component::Int)
get_brake(sys::SystemStructure, idx::Int16) = sys.winches[idx].brake
@register_symbolic get_brake(sys::SystemStructure, idx::Int16)
get_fix_point_sphere(sys::SystemStructure, idx::Int16) =
    sys.points[idx].fix_sphere
@register_symbolic get_fix_point_sphere(sys::SystemStructure, idx::Int16)
get_fix_static(sys::SystemStructure, idx::Int16) =
    sys.points[idx].fix_static
@register_symbolic get_fix_static(sys::SystemStructure, idx::Int16)
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
