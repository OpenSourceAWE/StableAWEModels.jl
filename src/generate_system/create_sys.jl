# Copyright (c) 2025 Bart van de Lint
# SPDX-License-Identifier: MPL-2.0

# Main system creation entry point

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
    eqs = Equation[]
    defaults = Pair{Num,Any}[]
    guesses = Pair{Num,Any}[]

    @unpack points, groups, segments, pulleys, tethers, winches, wings = system

    # Validation for REFINE wings
    for wing in wings
        if wing.wing_type == REFINE
            # REFINE wings cannot have groups
            @assert length(wing.group_idxs) == 0 "REFINE wing $(wing.idx) cannot have groups"
            @assert !isnothing(wing.point_to_vsm_point) "REFINE wing $(wing.idx) missing point_to_vsm_point mapping"

            # Verify all WING points for this wing are in the mapping
            wing_point_idxs = [p.idx for p in points if p.type == WING && p.wing_idx == wing.idx]
            for point_idx in wing_point_idxs
                @assert haskey(wing.point_to_vsm_point, point_idx) "REFINE wing $(wing.idx) missing mapping for point $(point_idx)"
            end

            # Verify 1:1 correspondence: n_structural_points == 2 * n_sections
            n_sections = length(wing.vsm_wing.unrefined_sections)
            @assert length(wing_point_idxs) == 2 * n_sections "REFINE wing $(wing.idx): expected $(2*n_sections) points for $(n_sections) sections, got $(length(wing_point_idxs))"

            prn && println("✓ REFINE wing $(wing.idx) validated: $(length(wing_point_idxs)) points, $(n_sections) sections, $(length(wing.vsm_aero.panels)) panels")
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

    # ==================== INLINED FORCE_EQS! CONTENT ==================== #
    # The following variables and component calls were previously in force_eqs!

    # Declare group geometry symbolic variables
    if length(groups) > 0
        @variables begin
            group_y_airf(t)[eachindex(groups), 1:3]
            group_chord(t)[eachindex(groups), 1:3]
            group_le_pos(t)[eachindex(groups), 1:3]
        end
    else
        group_y_airf = nothing
        group_chord = nothing
        group_le_pos = nothing
    end

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
        point_drag_force(t)[1:3, eachindex(points)]
        spring_sum_force(t)[1:3, eachindex(points)]  # Accumulated spring/drag forces
        disturb_force(t)[1:3, eachindex(points)]
        tether_r(t)[1:3, eachindex(points)]
        point_mass(t)[eachindex(points)]
        chord_b(t)[1:3, eachindex(points)]
        fixed_pos(t)[1:3, eachindex(points)]
        normal(t)[1:3, eachindex(points)]
        pos_b(t)[1:3, eachindex(points)]
        fix_point_sphere(t)[eachindex(points)]
        fix_static(t)[eachindex(points)]
        # Point damping (per-axis)
        body_frame_damping(t)[1:3, eachindex(points)]
        world_frame_damping(t)[1:3, eachindex(points)]
        # Segment forces and rest length
        spring_force_vec(t)[1:3, eachindex(segments)]
        drag_force(t)[1:3, eachindex(segments)]
        l0(t)[eachindex(segments)]
    end

    # Per-point variables for all points
    @variables begin
        va_point_b(t)[1:3, eachindex(points)]
        va_point_w(t)[1:3, eachindex(points)]
        wind_at_point(t)[1:3, eachindex(points)]
        height(t)[eachindex(points)]
    end

    # REFINE-specific variables: per-point aero forces in body frame
    if has_refine_wings
        @variables begin
            aero_force_point_b(t)[1:3, eachindex(points)]
        end
    else
        aero_force_point_b = nothing
    end

    # Pulley and tether length state variables (shared between segments, pulleys, winches)
    @variables begin
        pulley_len(t)[eachindex(pulleys)]
        pulley_vel(t)[eachindex(pulleys)]
        tether_len(t)[eachindex(winches)]
        tether_vel(t)[eachindex(winches)]
    end

    # ==================== CALL COMPONENT FUNCTIONS ==================== #

    # 1. Point equations (generates point dynamics, modifies tether_wing_force/moment in-place)
    eqs, defaults, guesses = point_eqs!(
        s, eqs, defaults, guesses, points, segments, groups, wings, psys, pset;
        R_b_w, wing_pos, wing_vel, wind_vec_gnd, fix_wing, twist_angle,
        pos, vel, acc, point_force, point_mass, spring_force_vec, drag_force, l0,
        spring_sum_force, point_drag_force, disturb_force, tether_r, chord_b, fixed_pos, normal, pos_b,
        fix_point_sphere, fix_static, body_frame_damping, world_frame_damping,
        va_point_b, va_point_w, wind_at_point, height,
        aero_force_point_b, has_refine_wings,
        group_y_airf, tether_wing_force, tether_wing_moment
    )

    # 2. Group equations (deformable wing sections with twist dynamics)
    eqs, defaults, guesses = group_eqs!(
        eqs, defaults, guesses, groups, wings, psys, pset;
        R_b_w, fix_wing, twist_angle, twist_ω, group_aero_moment,
        point_force, tether_wing_moment, group_y_airf, group_chord, group_le_pos
    )

    # 3. Segment equations (spring-damper forces, returns len and spring_force)
    eqs, guesses, len, spring_force = segment_eqs!(
        s, eqs, guesses, points, segments, pulleys, tethers, winches, psys, pset;
        pos, vel, wind_vec_gnd, spring_force_vec, drag_force, l0,
        pulley_len, tether_len
    )

    # 4. Pulley equations (rope distribution)
    eqs, defaults, guesses = pulley_eqs!(
        eqs, defaults, guesses, pulleys, segments, psys, pset;
        spring_force, pulley_len, pulley_vel
    )

    # 5. Winch equations (motor dynamics, tether reeling)
    eqs, defaults = winch_eqs!(
        eqs, defaults, winches, tethers, points, psys, pset;
        point_force, set_values, tether_len, tether_vel
    )

    # 6. Tether equations (stretched length, average force)
    eqs = tether_eqs!(eqs, tethers; len, spring_force)

    # ==================== END INLINED FORCE_EQS! CONTENT ==================== #

    # Build linearized aerodynamic equations
    if has_refine_wings
        eqs, guesses = linear_vsm_eqs!(
            s, eqs, guesses, psys;
            aero_force_b, R_v_w, aero_moment_b, group_aero_moment,
            twist_angle, va_wing_b, wing_pos, ω_b,
            aero_force_point_b=aero_force_point_b,
            va_point_b=va_point_b
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
