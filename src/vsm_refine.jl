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
    build_point_to_panel_mapping(wing_points::Vector{Point}, vsm_aero::VortexStepMethod.BodyAerodynamics)

Build a mapping from structural wing points to VSM panels for force lumping.

Each structural point receives forces from nearby VSM panels, weighted by inverse distance.
Returns a Dict mapping point_idx => [(panel_idx, weight), ...].

# Arguments
- `wing_points::Vector{Point}`: Structural WING-type points from struc_geometry.yaml
- `vsm_aero::VortexStepMethod.BodyAerodynamics`: VSM aerodynamics with panel geometry

# Returns
- `Dict{Int16, Vector{Tuple{Int16, Float64}}}`: Mapping point_idx -> [(panel_idx, weight), ...]
"""
function build_point_to_panel_mapping(
    wing_points::Vector{Point},
    vsm_aero::VortexStepMethod.BodyAerodynamics
)
    point_to_panels = Dict{Int16, Vector{Tuple{Int16, Float64}}}()

    # Get panels from BodyAerodynamics (panels are stored at the body level, not per wing)
    panels = vsm_aero.panels

    for point in wing_points
        # Find distances from this point to all panel centers
        panel_distances = Tuple{Int16, Float64}[]

        for (panel_idx, panel) in enumerate(panels)
            # Panel center position in body frame (average of 4 corner points)
            panel_center = 0.25 * (panel.LE_point_1 + panel.LE_point_2 +
                                   panel.TE_point_1 + panel.TE_point_2)

            # Distance from structural point to panel center
            dist = norm(point.pos_cad - panel_center)
            push!(panel_distances, (Int16(panel_idx), dist))
        end

        # Sort by distance and keep 3 nearest panels
        sort!(panel_distances, by=x->x[2])
        n_nearest = min(3, length(panel_distances))
        nearest = panel_distances[1:n_nearest]

        # Inverse distance weighting (normalized)
        weights = [1.0/d for (_, d) in nearest]
        total_weight = sum(weights)
        weights ./= total_weight

        # Store mapping
        point_to_panels[point.idx] = [(idx, w) for ((idx, _), w) in zip(nearest, weights)]
    end

    return point_to_panels
end

"""
    extract_panel_forces_to_vsm_state!(wing::VSMWing)

Extract per-panel 3D forces from VSM solution and store in wing.vsm_x for linearization.

The VSM solver provides `f_body_3D::Matrix{Float64}` with shape [3, n_panels], containing
the force vector (in body frame) for each panel. This function flattens this into the
vsm_x vector as: [fx_1, fy_1, fz_1, fx_2, fy_2, fz_2, ...].

# Arguments
- `wing::VSMWing`: Wing with REFINE type and solved VSM state
"""
function extract_panel_forces_to_vsm_state!(wing::VSMWing)
    @assert wing.wing_type == REFINE "Can only extract panel forces for REFINE wings"

    # Get panel forces from VSM solution: f_body_3D is [3, n_panels]
    n_panels = size(wing.vsm_solver.sol.f_body_3D, 2)

    # Flatten into vsm_x: [fx_1, fy_1, fz_1, fx_2, ...]
    for panel_idx in 1:n_panels
        wing.vsm_x[3*(panel_idx-1) + 1] = wing.vsm_solver.sol.f_body_3D[1, panel_idx]
        wing.vsm_x[3*(panel_idx-1) + 2] = wing.vsm_solver.sol.f_body_3D[2, panel_idx]
        wing.vsm_x[3*(panel_idx-1) + 3] = wing.vsm_solver.sol.f_body_3D[3, panel_idx]
    end

    return nothing
end

"""
    update_vsm_wing_from_structure!(wing::VSMWing, points::Vector{Point})

Update VSM wing section LE/TE positions from structural point positions using smooth
inverse distance weighting interpolation.

This creates two-way coupling: structural deformation → VSM wing sections → panels regenerated → aero forces.

Follows VortexStepMethod's deform! pattern: modify sections, then call reinit! to rebuild panels.

# Arguments
- `wing::VSMWing`: Wing with REFINE type
- `points::Vector{Point}`: All structural points (will filter for WING type)
"""
function update_vsm_wing_from_structure!(wing::VSMWing, points::Vector{Point})
    @assert wing.wing_type == REFINE "Can only update wing geometry for REFINE wings"

    # Get structural WING points for this wing
    wing_points = [p for p in points if p.type == WING && p.wing_idx == wing.idx]

    isempty(wing_points) && return  # No structural points to update from

    # Update each wing section from nearby structural points
    # Sections define the wing geometry; panels are derived from sections
    for section in wing.vsm_wing.sections
        # Find section center in CAD frame for distance calculation
        section_center_cad = 0.5 * (section.LE_point + section.TE_point)

        # Find 3 nearest structural points
        distances = [norm(section_center_cad - p.pos_cad) for p in wing_points]
        nearest_indices = sortperm(distances)[1:min(3, length(distances))]

        # Inverse distance weighting
        weights = [1.0/distances[i] for i in nearest_indices]
        total_weight = sum(weights)
        weights ./= total_weight

        # Calculate weighted average displacement from CAD to current world position
        displacement = sum(weights[i] * (wing_points[nearest_indices[i]].pos_w - wing_points[nearest_indices[i]].pos_cad)
                          for i in 1:length(nearest_indices))

        # Update section LE and TE positions (sections are mutable)
        # This follows the deform! pattern: modify sections directly, no reinit! on wing
        section.LE_point .+= displacement
        section.TE_point .+= displacement
    end

    # Do NOT call reinit! on wing - only modify sections!
    # body_aero.reinit! will update panels from modified sections (called in update_vsm!)
    return nothing
end
