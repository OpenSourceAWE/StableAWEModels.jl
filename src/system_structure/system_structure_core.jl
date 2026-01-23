# Copyright (c) 2025 Bart van de Lint
# SPDX-License-Identifier: MPL-2.0

"""
SystemStructure type and main constructor.

This file contains:
- SystemStructure struct definition
- Property accessors (getproperty/setproperty!)
- Main SystemStructure constructor with initialization logic
"""

# ==================== SYSTEM STRUCTURE ==================== #

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
    vsm_set::Union{Nothing, VortexStepMethod.VSMSettings}
end

function Base.getproperty(sys::SystemStructure, sym::Symbol)
    if sym == :total_mass
        # Sum of all point total_mass values (computed during simulation)
        # Falls back to extra_mass if total_mass is 0
        total = 0.0
        for point in getfield(sys, :points)
            if point.total_mass > 0
                total += point.total_mass
            else
                total += point.extra_mass
            end
        end
        return total
    elseif sym == :diff_vars
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
    calc_heading(sys::SystemStructure)

Calculate heading angles for all wings in the system structure.
Returns a vector of heading angles, one per wing.
"""
function calc_heading(sys::SystemStructure)
    wind_norm = normalize(sys.wind_vec_gnd)
    return [calc_heading(wing.R_b_w, wind_norm) for wing in sys.wings]
end

# ==================== CONSTRUCTOR ==================== #

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
        vsm_set=nothing,
    )
    # Load VSMSettings if not provided and wings exist
    if isnothing(vsm_set) && !isempty(wings)
        model_dir = get_data_path()
        vsm_set_path = joinpath(model_dir, "vsm_settings.yaml")
        if isfile(vsm_set_path)
            vsm_set = VortexStepMethod.VSMSettings(vsm_set_path; data_prefix=false)
        end
    end

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
                    extra_mass = points[idx].extra_mass,
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
    # Auto-create groups for QUATERNION wings if needed (before geometry initialization)
    for (i, wing) in enumerate(wings)
        if wing.wing_type == QUATERNION && isempty(wing.group_idxs)
            # Get WING-type points for this wing
            wing_point_idxs = findall(
                p -> p.type == WING && p.wing_idx == wing.idx, points)
            wing_points = [points[idx] for idx in wing_point_idxs]

            # Identify LE/TE pairs
            wing_segments = identify_wing_segments(wing_points)

            # Create a group for each section (LE/TE pair)
            # n_groups = n_unrefined_sections (one group per section)
            vsm_wing = wing.vsm_wing
            new_group_idxs = Int16[]

            # Check if wing has interpolators (from .obj) or not (from YAML)
            has_interpolators = !isnothing(vsm_wing.le_interp)

            # For YAML wings, calculate COM from WING points and update transforms
            if !has_interpolators && !isempty(wing_points)
                calculated_com = mean([p.pos_cad for p in wing_points])
                wing.pos_cad .= calculated_com
                wing.vsm_wing.T_cad_body .= calculated_com
                adjust_vsm_panels_to_origin!(vsm_wing, calculated_com)
                # Apply aero_z_offset after COM adjustment
                apply_aero_z_offset!(vsm_wing, wing.aero_z_offset)
            end

            for (le_idx, te_idx) in wing_segments
                group_idx = length(groups) + 1

                if has_interpolators
                    # For .obj wings, calculate gamma from LE position
                    le_point = points[le_idx]
                    y_le = le_point.pos_cad[2]
                    z_le = le_point.pos_cad[3]
                    # Compute circle_center_z = middle_le_z - radius
                    circle_center_z = vsm_wing.le_interp[3](0.0) - vsm_wing.radius
                    gamma = atan(-y_le, z_le - circle_center_z)

                    # Use constructor with vsm_wing (computes geometry from gamma)
                    new_group = Group(group_idx, [le_idx, te_idx],
                                     vsm_wing, gamma, DYNAMIC, 0.25)
                else
                    # For YAML wings, gamma concept doesn't apply
                    # Use simple constructor (geometry computed from points later)
                    new_group = Group(group_idx, [le_idx, te_idx],
                                     0.0, DYNAMIC, 0.25)
                end

                push!(groups, new_group)
                push!(new_group_idxs, Int16(group_idx))
            end

            # Update wing with new groups and resize vsm arrays
            wing.group_idxs = new_group_idxs

            # Resize vsm arrays based on number of unrefined sections
            n_unrefined = wing.vsm_wing.n_unrefined_sections
            ny = 3 + n_unrefined + 3  # va(3) + twist(n_unrefined) + ω(3)
            nx = 3 + 3 + n_unrefined  # force(3) + moment(3) + unrefined_moments(n_unrefined)
            wing.vsm_y = zeros(SimFloat, ny)
            wing.vsm_x = zeros(SimFloat, nx)
            wing.vsm_jac = zeros(SimFloat, nx, ny)

            @info "Auto-created $(length(new_group_idxs)) groups " *
                  "for QUATERNION wing $(wing.idx)" *
                  (!has_interpolators && !isempty(wing_points) ?
                   " (COM from WING points: " *
                   "[$(round(wing.pos_cad[1], digits=2)), " *
                   "$(round(wing.pos_cad[2], digits=2)), " *
                   "$(round(wing.pos_cad[3], digits=2))])" : "")
        end
    end

    """
        compute_spatial_group_mapping!(the_wing, groups, points)

    Map groups to unrefined sections using spatial proximity.
    Each group is assigned to the closest unrefined section based on distance between centers.
    """
    function compute_spatial_group_mapping!(the_wing::VSMWing, groups::Vector{Group}, points::Vector{Point})
        the_vsm_wing = the_wing.vsm_wing
        n_unrefined = the_vsm_wing.n_unrefined_sections
        n_groups = length(the_wing.base.group_idxs)

        # Compute group centers in body frame
        group_centers = Vector{MVec3}(undef, n_groups)
        for (local_idx, group_idx) in enumerate(the_wing.base.group_idxs)
            group = groups[group_idx]
            le_idx = group.point_idxs[1]
            te_idx = group.point_idxs[2]
            le_pos_b = the_wing.base.R_b_c' * (points[le_idx].pos_cad - the_wing.base.pos_cad)
            te_pos_b = the_wing.base.R_b_c' * (points[te_idx].pos_cad - the_wing.base.pos_cad)
            group_centers[local_idx] = (le_pos_b + te_pos_b) / 2
        end

        # Compute unrefined section centers
        unrefined_centers = Vector{MVec3}(undef, n_unrefined)
        for i in 1:n_unrefined
            le_point = the_vsm_wing.unrefined_sections[i].LE_point
            te_point = the_vsm_wing.unrefined_sections[i].TE_point
            unrefined_centers[i] = (le_point + te_point) / 2
        end

        # Map each group to closest unrefined section
        for (local_idx, group_idx) in enumerate(the_wing.base.group_idxs)
            group = groups[group_idx]
            min_dist = Inf
            closest_idx = 1
            for unrefined_idx in 1:n_unrefined
                dist = norm(group_centers[local_idx] - unrefined_centers[unrefined_idx])
                if dist < min_dist
                    min_dist = dist
                    closest_idx = unrefined_idx
                end
            end
            group.unrefined_section_idxs = Int16[closest_idx]
        end

        # Validate: check all sections are covered
        assigned = Set{Int16}()
        for group_idx in the_wing.base.group_idxs
            union!(assigned, groups[group_idx].unrefined_section_idxs)
        end
        if length(assigned) != n_unrefined
            unassigned = setdiff(1:n_unrefined, assigned)
            @warn "Wing $(the_wing.base.idx): $(length(unassigned)) unrefined sections not assigned to any group: $unassigned"
        end
    end

    # Initialize group-to-unrefined-section mapping for QUATERNION wings
    # Do this BEFORE y_airf calculation so the mapping is available
    for the_wing in wings
        if isa(the_wing, VSMWing) && the_wing.base.wing_type == QUATERNION && !isempty(the_wing.base.group_idxs)
            compute_spatial_group_mapping!(the_wing, groups, points)
        end
    end

    # Initialize group geometries from VSM wing or point positions
    for group in groups
        if iszero(group.chord)
            # Find which wing this group belongs to
            for wing in wings
                if group.idx in wing.group_idxs
                    vsm_wing = wing.vsm_wing

                    if !isnothing(vsm_wing.le_interp)
                        # For .obj wings: use interpolators with gamma
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
                    else
                        # For YAML wings: compute from point positions
                        # group.point_idxs contains [le_idx, te_idx]
                        @assert length(group.point_idxs) >= 2 "Group $(group.idx) needs at least LE and TE points"
                        le_idx = group.point_idxs[1]
                        te_idx = group.point_idxs[2]

                        # Calculate pos_b manually (same as done in reinit!)
                        # pos_b = R_b_c' * (pos_cad - wing.pos_cad)
                        le_point = points[le_idx]
                        te_point = points[te_idx]

                        le_pos_b = wing.R_b_c' * (le_point.pos_cad - wing.pos_cad)
                        te_pos_b = wing.R_b_c' * (te_point.pos_cad - wing.pos_cad)

                        group.le_pos .= le_pos_b
                        group.chord .= te_pos_b .- le_pos_b

                        # y_airf: find the two closest non_deformed_sections
                        group_center = (le_pos_b .+ te_pos_b) ./ 2
                        # Find closest section
                        min_dist1 = Inf
                        closest_idx1 = 1
                        for i in 1:length(vsm_wing.non_deformed_sections)
                            section = vsm_wing.non_deformed_sections[i]
                            section_center = (section.LE_point .+ section.TE_point) ./ 2
                            dist = norm(group_center .- section_center)
                            if dist < min_dist1
                                min_dist1 = dist
                                closest_idx1 = i
                            end
                        end

                        # Use adjacent section to compute local_y
                        section1 = vsm_wing.non_deformed_sections[closest_idx1]
                        if closest_idx1 < length(vsm_wing.non_deformed_sections)
                            section2 = vsm_wing.non_deformed_sections[closest_idx1 + 1]
                            local_y = normalize(section1.LE_point .- section2.LE_point)
                        else
                            section_prev = vsm_wing.non_deformed_sections[closest_idx1 - 1]
                            local_y = normalize(section_prev.LE_point .- section1.LE_point)
                        end
                        group.y_airf .= local_y
                    end
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
                    points[point_idx].extra_mass = mass_per_point
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
        # Use number of unrefined sections
        n_unrefined = wings[1].vsm_wing.n_unrefined_sections
        ny = 3 + n_unrefined + 3
        nx = 3 + 3 + n_unrefined
    else
        ny = 0
        nx = 0
    end
    y = zeros(length(wings), ny)
    x = zeros(length(wings), nx)
    jac = zeros(length(wings), nx, ny)
    set.physical_model = name
    sys_struct = SystemStructure(name, set, points, groups, segments, pulleys, tethers,
        winches, wings, transforms, y, x, jac, zeros(KVec3), 0.0, false, false, vsm_set)
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
