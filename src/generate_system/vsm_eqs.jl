# Copyright (c) 2025 Bart van de Lint
# SPDX-License-Identifier: MPL-2.0

# VSM linearized aerodynamics equation generation

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
    aero_force_point_b=nothing,  # REFINE-specific parameter
    va_point_b=nothing            # REFINE-specific parameter for per-point va
)
    @unpack groups, wings, points = s.sys_struct
    if length(wings) == 0
        return eqs, guesses
    end

    # Check for REFINE and QUATERNION wings
    has_refine = any(w.wing_type == REFINE for w in wings)
    has_quaternion = any(w.wing_type == QUATERNION for w in wings)

    # Compute maximum dimensions across all wings
    # ny = dimension of VSM input state for QUATERNION wings: [va(3), twist(n_unrefined), ω(3)]
    # Only compute if we have QUATERNION wings
    if has_quaternion
        n_unrefined = wings[1].vsm_wing.n_unrefined_sections
        ny_quaternion = 3 + n_unrefined + 3
    else
        ny_quaternion = 0
    end

    # nx = dimension of VSM output (forces/moments for QUATERNION, panel forces for REFINE)
    nx_values = Int[]
    for wing in wings
        if wing.wing_type == REFINE
            # REFINE: Store per-panel forces [fx_1, fy_1, fz_1, fx_2, ...]
            n_panels = length(wing.vsm_aero.panels)
            push!(nx_values, 3 * n_panels)
        else
            # QUATERNION: Store [aero_force(3), aero_moment(3), unrefined_moments(n_unrefined)]
            n_unrefined = wing.vsm_wing.n_unrefined_sections
            push!(nx_values, 3 + 3 + n_unrefined)
        end
    end
    nx_max = maximum(nx_values)

    # Declare VSM output forces for all wings (used by both REFINE and QUATERNION)
    @variables begin
        vsm_output_force_prev(t)[1:nx_max, eachindex(wings)]  # Panel forces (REFINE) or [forces, moments] (QUATERNION)
    end

    # Declare linearization variables only for QUATERNION wings
    if has_quaternion
        @variables begin
            # VSM linearization variables (QUATERNION wings only)
            # Linearization: F(state) ≈ F(state₀) + ∂F/∂state|_{state₀} * (state - state₀)
            # where state = [va_wing_b, ω_b, twist_angles] and F = [forces, moments]
            vsm_input_state(t)[1:ny_quaternion, eachindex(wings)]        # Current input state
            vsm_input_state_delta(t)[1:ny_quaternion, eachindex(wings)]  # Δstate for linearization
            vsm_input_state_prev(t)[1:ny_quaternion, eachindex(wings)]   # State at linearization point
            force_jacobian(t)[1:nx_max, 1:ny_quaternion, eachindex(wings)]  # ∂F/∂state Jacobian matrix

            # Aerodynamic quantities (QUATERNION wings only)
            q_inf(t)[eachindex(wings)]                      # Dynamic pressure
            no_scale_aero_force_b(t)[1:3, eachindex(wings)]  # Unscaled force before q_inf scaling
        end
    end

    for wing in wings
        if wing.wing_type == REFINE
            # ==================== REFINE WING: Direct Section Forces ====================
            # REFINE wings use nonlinear VSM solve each timestep. Unrefined section forces
            # map 1:1 to structural sections (LE+TE points). No linearization.

            n_panels = length(wing.vsm_aero.panels)
            nx_refine = 3 * n_panels

            # Retrieve current panel forces from VSM solution
            eqs = [
                eqs
                [vsm_output_force_prev[ix, wing.idx] ~ get_vsm_x(psys, wing.idx, ix) for ix = 1:nx_refine]
            ]

            # Aero forces computed by distribute_panel_forces_to_points! and stored in point.aero_force_b
            wing_points = [p for p in points if p.type == WING && p.wing_idx == wing.idx]

            for point in wing_points
                eqs = [
                    eqs
                    aero_force_point_b[:, point.idx] ~
                        [get_point_aero_force(psys, point.idx, i) for i in 1:3]
                ]
            end

            # Total wing force is sum of all point forces
            eqs = [
                eqs
                aero_force_b[:, wing.idx] ~ sum([aero_force_point_b[:, p.idx] for p in wing_points])
                aero_moment_b[:, wing.idx] ~ zeros(3)
            ]

        else
            # ==================== QUATERNION WING: Linearized Aerodynamics ====================
            # Linearized forces: F ≈ F₀ + J*(state - state₀)
            # State: [va_b(3), twist(n_unrefined), ω_b(3)]
            # Output: [force(3), moment(3), unrefined_moments(n_unrefined)]

            area = wing.vsm_aero.projected_area
            force_b = no_scale_aero_force_b[:, wing.idx]
            wind_direction_b = sym_normalize(va_wing_b[:, wing.idx])
            drag_force_b = (force_b ⋅ wind_direction_b) * wind_direction_b

            # Dimensions for this wing
            n_unrefined = wing.vsm_wing.n_unrefined_sections
            nx_quat = 3 + 3 + n_unrefined

            eqs = [
                eqs
                # Dynamic pressure for force scaling
                q_inf[wing.idx] ~
                    0.5 * calc_rho(s.am, wing_pos[3, wing.idx]) *
                    norm(collect(va_wing_b[:, wing.idx]))^2

                # Load linearization data from VSM (state₀, F₀, ∂F/∂state)
                [vsm_input_state_prev[iy, wing.idx] ~ get_vsm_y(psys, wing.idx, iy) for iy = 1:ny_quaternion]
                [vsm_output_force_prev[ix, wing.idx] ~ get_vsm_x(psys, wing.idx, ix) for ix = 1:nx_quat]
                [
                    force_jacobian[ix, iy, wing.idx] ~
                        get_vsm_jac(psys, wing.idx, ix, iy) for ix = 1:nx_quat for
                    iy = 1:ny_quaternion
                ]
            ]

            # Build mapping from unrefined section index to group twist
            # Create an array where unrefined_twists[i] = twist of group containing section i
            unrefined_to_group_twist = Vector{Any}(undef, n_unrefined)
            for group_idx in wing.group_idxs
                group = groups[group_idx]
                for unrefined_idx in group.unrefined_section_idxs
                    unrefined_to_group_twist[unrefined_idx] = twist_angle[group.idx]
                end
            end

            eqs = [
                eqs
                # Current input state for linearization
                vsm_input_state[:, wing.idx] ~ [
                    va_wing_b[:, wing.idx]       # Apparent wind velocity in body frame
                    unrefined_to_group_twist       # Twist angles per unrefined section
                    ω_b[:, wing.idx]              # Angular velocity in body frame
                ]

                # State deviation from linearization point: Δstate = state - state₀
                vsm_input_state_delta[:, wing.idx] ~ vsm_input_state[:, wing.idx] - vsm_input_state_prev[:, wing.idx]
            ]

            # Map unrefined moments back to groups
            # Sum moments from all unrefined sections belonging to each group
            group_moment_eqs = []
            for group_idx in wing.group_idxs
                group = groups[group_idx]
                if !isempty(group.unrefined_section_idxs)
                    # Indices in vsm_output_force_prev: [force(3), moment(3), unrefined_moments(...)]
                    # unrefined_moments start at index 7
                    moment_terms = []
                    for unrefined_idx in group.unrefined_section_idxs
                        vsm_output_idx = 6 + unrefined_idx
                        # Linearized moment: moment₀ + jacobian * delta_state
                        linearized_moment = (
                            vsm_output_force_prev[vsm_output_idx, wing.idx] +
                            sum([force_jacobian[vsm_output_idx, iy, wing.idx] *
                                 vsm_input_state_delta[iy, wing.idx]
                                 for iy in 1:ny_quaternion])
                        )
                        push!(moment_terms, linearized_moment)
                    end
                    push!(group_moment_eqs,
                          group_aero_moment[group.idx] ~ sum(moment_terms))
                end
            end

            eqs = [
                eqs
                # ===== LINEARIZED FORCE CALCULATION =====
                # F ≈ F₀ + ∂F/∂state * Δstate
                # Then scale by dynamic pressure and wing area
                [
                    force_b
                    aero_moment_b[:, wing.idx]
                ] ~
                    q_inf[wing.idx] *
                    area *
                    (vsm_output_force_prev[1:6, wing.idx] +
                     force_jacobian[1:6, :, wing.idx] * vsm_input_state_delta[:, wing.idx])

                group_moment_eqs

                # Apply additional drag correction factor
                aero_force_b[:, wing.idx] ~
                    force_b + drag_force_b * (get_drag_frac(psys, wing.idx) - 1)
            ]

            if s.set.quasi_static
                guesses = [
                    guesses
                    [vsm_input_state[iy, wing.idx] => get_vsm_y(psys, wing.idx, iy) for iy = 1:ny_quaternion]
                ]
            end
        end
    end
    return eqs, guesses
end
