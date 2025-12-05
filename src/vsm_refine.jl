# Copyright (c) 2025 Bart van de Lint
# SPDX-License-Identifier: MPL-2.0

"""
Helper functions for REFINE wing type that applies VSM panel forces directly to
structural points.
"""

"""
    identify_wing_segments(wing_points::Vector{Point})

Identify wing segments (LE/TE pairs) from WING-type points.

Assumes points are organized in pairs along the span, with even-numbered points
being leading edge and odd-numbered being trailing edge (or vice versa).

# Arguments
- `wing_points::Vector{Point}`: All WING-type points for a wing (sorted by index)

# Returns
- `Vector{Tuple{Int16, Int16}}`: Vector of (le_point_idx, te_point_idx) pairs defining segments
"""
function identify_wing_segments(wing_points::Vector{Point})
    # Sort points by index to ensure consistent ordering
    sorted_points = sort(wing_points, by=p->p.idx)

    n_points = length(sorted_points)
    @assert n_points % 2 == 0 "Wing must have even number of points (LE/TE pairs)"

    n_segments = n_points ÷ 2
    segments = Tuple{Int16, Int16}[]

    # Group consecutive pairs: (point[1], point[2]), (point[3], point[4]), ...
    for i in 1:n_segments
        le_idx = 2*i - 1
        te_idx = 2*i
        le_point = sorted_points[le_idx]
        te_point = sorted_points[te_idx]

        # Determine which is LE and which is TE by x-coordinate (LE has smaller x)
        if le_point.pos_cad[1] < te_point.pos_cad[1]
            push!(segments, (le_point.idx, te_point.idx))
        else
            push!(segments, (te_point.idx, le_point.idx))
        end
    end

    return segments
end

"""
    build_point_to_vsm_point_mapping(wing_points::Vector{Point}, vsm_wing::VortexStepMethod.AbstractWing)

Build 1:1 mapping from structural WING points to VSM wing section points (LE/TE) using closest-point distance.

For each VSM section point (LE/TE), finds the closest structural point in CAD frame.

# Constraint
Requires: `length(wing_points) == 2 * length(vsm_wing.sections)`

# Arguments
- `wing_points::Vector{Point}`: Structural WING-type points
- `vsm_wing::VortexStepMethod.AbstractWing`: VSM wing with sections

# Returns
- `Dict{Int16, Tuple{Int16, Symbol}}`: Mapping structural_point_idx -> (section_idx, :LE or :TE)

# Algorithm
1. For each section in vsm_wing.sections:
   - Find closest unused structural point to section.LE_point → assign to (section_idx, :LE)
   - Find closest unused structural point to section.TE_point → assign to (section_idx, :TE)
2. Distance measured in CAD/body frame using norm(point.pos_cad - section_point)
"""
function build_point_to_vsm_point_mapping(
    wing_points::Vector{Point},
    vsm_wing::VortexStepMethod.AbstractWing
)
    n_points = length(wing_points)
    n_sections = length(vsm_wing.sections)

    # Validate 1:1 correspondence constraint
    if n_points != 2 * n_sections
        error("REFINE wing requires n_structural_points ($(n_points)) == " *
              "2 * n_vsm_sections ($(n_sections))")
    end

    point_to_vsm_point = Dict{Int16, Tuple{Int16, Symbol}}()
    used_points = Set{Int16}()

    for (section_idx, section) in enumerate(vsm_wing.sections)
        # Map LE_point to closest unused structural point
        le_pos = section.LE_point
        min_dist = Inf
        closest_le_idx = wing_points[1].idx

        for point in wing_points
            if point.idx in used_points
                continue  # Already assigned
            end
            dist = norm(point.pos_cad - le_pos)
            if dist < min_dist
                min_dist = dist
                closest_le_idx = point.idx
            end
        end

        point_to_vsm_point[closest_le_idx] = (Int16(section_idx), :LE)
        push!(used_points, closest_le_idx)

        # Map TE_point to closest unused structural point
        te_pos = section.TE_point
        min_dist = Inf
        closest_te_idx = wing_points[1].idx

        for point in wing_points
            if point.idx in used_points
                continue  # Already assigned
            end
            dist = norm(point.pos_cad - te_pos)
            if dist < min_dist
                min_dist = dist
                closest_te_idx = point.idx
            end
        end

        point_to_vsm_point[closest_te_idx] = (Int16(section_idx), :TE)
        push!(used_points, closest_te_idx)
    end

    return point_to_vsm_point
end

"""
    compute_panel_le_te_forces(panel, cl, cd, cm, alpha_corrected, density, v_a_mag)

Compute leading-edge and trailing-edge forces from aerodynamic coefficients.

Given a panel with aerodynamic coefficients at quarter-chord, this function:
1. Calculates total lift and drag forces using dynamic pressure
2. Uses pitching moment to determine center of pressure location
3. Splits forces between LE and TE to satisfy moment equilibrium

# Arguments
- `panel`: VSM Panel (or group panel) with geometry (x_airf, y_airf, z_airf, chord, width, va)
- `cl`: Lift coefficient [-]
- `cd`: Drag coefficient [-]
- `cm`: Pitching moment coefficient at quarter-chord [-]
- `alpha_corrected`: Corrected angle of attack (including induced velocity effects) [rad]
- `density`: Air density [kg/m³]
- `v_a_mag`: Apparent velocity magnitude [m/s]

# Returns
- `(F_LE, F_TE)`: Tuple of 3D force vectors at leading edge and trailing edge [N]

# Algorithm
1. Dynamic pressure: q = 0.5 × ρ × v_a²
2. Total forces: L = cl × q × chord × width, D = cd × q × chord × width
3. Center of pressure (normalized): x_cp = 0.25 + cm/cl
4. Lift split: L_LE = L × (1 - x_cp), L_TE = L × x_cp
5. Drag split: D_LE = D/2, D_TE = D/2
6. Force directions computed using corrected alpha (matches VortexStepMethod algorithm)

# Note
For REFINE wings, pass the group panel which has averaged geometry and total area.
Uses corrected alpha to properly account for induced velocity effects on force direction.
"""
function compute_panel_le_te_forces(panel, cl, cd, cm, alpha_corrected, density, v_a_mag)
    # Dynamic pressure
    q = 0.5 * density * v_a_mag^2

    # Total lift and drag magnitudes
    L_total = cl * q * panel.chord * panel.width
    D_total = cd * q * panel.chord * panel.width

    # Center of pressure location (normalized by chord)
    # x_cp = 0.25 (quarter-chord) + cm/cl (moment arm)
    x_cp = if abs(cl) > 1e-6
        0.25 + cm / cl
    else
        0.25  # Default to quarter-chord if no lift
    end

    # Clamp to reasonable range [0, 1]
    x_cp = clamp(x_cp, 0.0, 1.0)

    # Split lift between LE and TE using moment equilibrium
    # Moment about LE: L_TE × chord = L_total × (x_cp × chord)
    L_TE = L_total * x_cp
    L_LE = L_total * (1.0 - x_cp)

    # Split drag equally
    D_LE = D_total / 2.0
    D_TE = D_total / 2.0

    # Compute force directions using corrected alpha (matches VortexStepMethod)
    # Induced apparent wind direction in airfoil frame
    dir_induced_va_airfoil = cos(alpha_corrected) * panel.x_airf +
                             sin(alpha_corrected) * panel.z_airf

    # Drag direction: along induced apparent wind
    drag_dir = dir_induced_va_airfoil / (norm(dir_induced_va_airfoil) + 1e-12)

    # Lift direction: perpendicular to drag and spanwise axis
    lift_dir = cross(drag_dir, panel.y_airf)
    lift_dir_mag = norm(lift_dir)
    if lift_dir_mag > 1e-12
        lift_dir = lift_dir / lift_dir_mag
    else
        # Fallback to z_airf if cross product is degenerate
        lift_dir = panel.z_airf
    end

    # Combine lift and drag components at LE and TE
    F_LE = L_LE * lift_dir + D_LE * drag_dir
    F_TE = L_TE * lift_dir + D_TE * drag_dir

    return (F_LE, F_TE)
end

"""
    distribute_panel_forces_to_points!(wing::VSMWing, points::Vector{Point})

Distribute VSM panel forces to structural points using unrefined-level aerodynamic coefficients.

After VSM solve, uses pre-computed unrefined section coefficients (averaged over refined panels)
to compute LE and TE forces for each section, then distributes to structural points.

# Algorithm
1. Initialize all WING point aero_forces to zero
2. For each section (which maps 1:1 to unrefined panels):
   - Get cl, cd, cm from unrefined arrays
   - Get a representative panel from this unrefined section
   - Call compute_panel_le_te_forces() to get LE and TE forces
   - Map forces to structural points via point_to_vsm_point

# Mapping
- Unrefined sections in VortexStepMethod map 1:1 to structural wing segments
- Each unrefined panel connects two adjacent sections
- Each unrefined section aggregates coefficients from all refined panels
- Unrefined coefficients are pre-averaged by the VSM solver
- Forces are distributed 50/50 to the two sections that bound each panel

# Arguments
- `wing::VSMWing`: Wing with REFINE type and solved VSM state
- `points::Vector{Point}`: All structural points (will filter for WING type)
"""
function distribute_panel_forces_to_points!(wing::VSMWing, points::Vector{Point})
    @assert wing.wing_type == REFINE "Can only distribute forces for REFINE wings"

    # Get VSM solution data - use unrefined arrays
    cl_unrefined_array = wing.vsm_solver.sol.cl_unrefined_array
    cd_unrefined_array = wing.vsm_solver.sol.cd_unrefined_array
    cm_unrefined_array = wing.vsm_solver.sol.cm_unrefined_array
    alpha_unrefined_array = wing.vsm_solver.sol.alpha_unrefined_array
    unrefined = wing.vsm_aero.unrefined
    density = wing.vsm_solver.density
    n_sections = length(wing.vsm_wing.sections)

    # Initialize all WING point forces to zero
    for point in points
        if point.type == WING && point.wing_idx == wing.idx
            point.aero_force_b .= 0.0
        end
    end

    # Build inverse mapping: (section_idx, :LE/:TE) -> point_idx
    vsm_point_to_struct = Dict{Tuple{Int16, Symbol}, Int16}()
    for (point_idx, (section_idx, le_or_te)) in wing.point_to_vsm_point
        vsm_point_to_struct[(section_idx, le_or_te)] = point_idx
    end

    # Number of unrefined panels = number of sections - 1
    # Unrefined sections map 1:1 to panels for REFINE wings
    n_unrefined_panels = n_sections - 1

    # For each unrefined panel, distribute forces to adjacent sections
    for panel_idx in 1:n_unrefined_panels
        unrefined_idx = panel_idx  # 1:1 mapping

        # Get unrefined-level coefficients (averaged over all refined panels)
        cl = cl_unrefined_array[unrefined_idx]
        cd = cd_unrefined_array[unrefined_idx]
        cm = cm_unrefined_array[unrefined_idx]
        alpha_corrected = alpha_unrefined_array[unrefined_idx]

        # Get unrefined panel with averaged geometry and total area
        unrefined_panel = unrefined[unrefined_idx]
        v_a_mag = norm(unrefined_panel.va)

        # Panel connects section i and section i+1
        section_i_idx = panel_idx
        section_i_plus_1_idx = panel_idx + 1

        # Compute LE and TE forces using unrefined panel (has correct averaged geometry and total area)
        F_LE, F_TE = compute_panel_le_te_forces(
            unrefined_panel, cl, cd, cm, alpha_corrected, density, v_a_mag
        )

        # Distribute 50% of forces to each adjacent section
        # Section i gets 50%
        le_key_i = (Int16(section_i_idx), :LE)
        te_key_i = (Int16(section_i_idx), :TE)

        if haskey(vsm_point_to_struct, le_key_i)
            point_idx = vsm_point_to_struct[le_key_i]
            points[point_idx].aero_force_b .+= F_LE / 2.0
        end
        if haskey(vsm_point_to_struct, te_key_i)
            point_idx = vsm_point_to_struct[te_key_i]
            points[point_idx].aero_force_b .+= F_TE / 2.0
        end

        # Section i+1 gets 50%
        le_key_i_plus_1 = (Int16(section_i_plus_1_idx), :LE)
        te_key_i_plus_1 = (Int16(section_i_plus_1_idx), :TE)

        if haskey(vsm_point_to_struct, le_key_i_plus_1)
            point_idx = vsm_point_to_struct[le_key_i_plus_1]
            points[point_idx].aero_force_b .+= F_LE / 2.0
        end
        if haskey(vsm_point_to_struct, te_key_i_plus_1)
            point_idx = vsm_point_to_struct[te_key_i_plus_1]
            points[point_idx].aero_force_b .+= F_TE / 2.0
        end
    end

    return nothing
end

"""
    update_vsm_wing_from_structure!(wing::VSMWing, points::Vector{Point})

Update VSM section points (LE/TE) directly from structural point positions using 1:1 mapping.

This creates two-way coupling: structural deformation → VSM sections → aero forces.

# Algorithm
Uses direct 1:1 correspondence between structural points and VSM section points:
1. For each structural WING point:
   - Calculate current position in body frame: pos_b = R_b_w' * (pos_w - origin)
   - Find corresponding VSM section point (LE or TE) via wing.point_to_vsm_point
   - Set VSM section point directly: section.LE_point = pos_b (or TE_point)

# Notes
- Section points are stored in body frame coordinates
- `wing.R_b_w` and `wing.pos_w` are updated each timestep from structural geometry (symbolic equations)
- To get world coordinates: `world_pos = wing.R_b_w * section.LE_point + wing.pos_w`

# Arguments
- `wing::VSMWing`: Wing with REFINE type
- `points::Vector{Point}`: All structural points (will filter for WING type)
"""
function update_vsm_wing_from_structure!(wing::VSMWing, points::Vector{Point})
    @assert wing.wing_type == REFINE "Can only update wing geometry for REFINE wings"

    # Get current R_b_w and origin from wing state
    # (These are updated during simulation from structural geometry)
    R_b_w = wing.R_b_w
    origin = wing.pos_w

    # Update each VSM section point directly from its corresponding structural point
    for (point_idx, (section_idx, le_or_te)) in wing.point_to_vsm_point
        point = points[point_idx]

        # Calculate current position in body frame
        pos_b = R_b_w' * (point.pos_w - origin)

        # Get the section
        section = wing.vsm_wing.sections[section_idx]

        # Set section point directly to body frame position
        if le_or_te == :LE
            section.LE_point .= pos_b
        else  # :TE
            section.TE_point .= pos_b
        end
    end

    # Do NOT call reinit! on wing - only modify sections!
    # body_aero reinit! will update panels from modified sections (called in update_vsm!)
    return nothing
end
