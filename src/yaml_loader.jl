using YAML
using Logging
using LinearAlgebra
using VortexStepMethod, KiteUtils

# Default diameter for wing section segments (in mm)
DEFAULT_WING_DIAMETER_MM = 1.0

# ------------------ helpers ------------------

"""
    resolve_references(row::NamedTuple, property_tables::Dict{String, Dict{String, NamedTuple}})

Resolve string references in a row by looking them up in property tables.
If a field value is a String, check if it exists as a key in any property table.
If found, merge those properties into the current row (current row takes precedence).
"""
function resolve_references(row::NamedTuple, property_tables::Dict{String, Dict{String, NamedTuple}})
    resolved = Dict{Symbol, Any}(pairs(row))

    for (field, value) in pairs(row)
        # Check if this is a string reference (unquoted strings in YAML)
        if value isa String
            # Try to find this reference in any property table
            for (table_name, lookup_dict) in property_tables
                if haskey(lookup_dict, value)
                    # Found! Merge these properties (current row takes precedence)
                    ref_props = lookup_dict[value]
                    for (k, v) in pairs(ref_props)
                        if !haskey(resolved, k) || resolved[k] === nothing
                            resolved[k] = v
                        end
                    end
                    break
                end
            end
        end
    end

    return NamedTuple(resolved)
end

"""
    calculate_derived_properties!(props::Dict{Symbol, Any})

Calculate derived properties like axial_stiffness and axial_damping from material properties.
Modifies props in-place.
"""
function calculate_derived_properties!(props::Dict{Symbol, Any})
    # Calculate axial_stiffness from material properties if missing or if it's a string (material name)
    if haskey(props, :youngs_modulus) && haskey(props, :diameter_mm) && haskey(props, :l0)
        # Check if we need to calculate (missing, nothing, or is a string reference)
        need_calculation = !haskey(props, :axial_stiffness) ||
                          props[:axial_stiffness] === nothing ||
                          props[:axial_stiffness] isa String

        if need_calculation
            d_m = Float64(props[:diameter_mm]) / 1000.0  # mm to m
            A = π * (d_m / 2)^2
            E = Float64(props[:youngs_modulus])
            l0 = Float64(props[:l0])
            props[:axial_stiffness] = E * A / l0
        end
    end

    # Calculate axial_damping from damping coefficient if missing or is a string
    if haskey(props, :damping_per_stiffness) && haskey(props, :axial_stiffness)
        # Only calculate if axial_stiffness is now a number
        if props[:axial_stiffness] isa Number
            need_damping_calc = !haskey(props, :axial_damping) ||
                               props[:axial_damping] === nothing ||
                               props[:axial_damping] isa String

            if need_damping_calc
                props[:axial_damping] = Float64(props[:damping_per_stiffness]) * Float64(props[:axial_stiffness])
            end
        end
    end

    # Set default axial_damping if still missing
    if !haskey(props, :axial_damping) || props[:axial_damping] === nothing || props[:axial_damping] isa String
        props[:axial_damping] = 0.0
    end
end

function parse_table(tbl)::Vector{NamedTuple}
    haskey(tbl, "data") || throw(ArgumentError("table is missing `data`"))

    rows = tbl["data"]
    isempty(rows) && return NamedTuple[]

    # Check format: if first row is a Dict, use dict format; if Array, use header format
    first_row = first(rows)

    if first_row isa AbstractDict
        # Dict format: each row is already a dict with named keys
        # Convert each dict to a NamedTuple
        out = NamedTuple[]
        for row in rows
            nt = NamedTuple{Tuple(Symbol.(keys(row)))}(Tuple(values(row)))
            push!(out, nt)
        end
        return out
    else
        # Array format: requires headers
        haskey(tbl, "headers") || throw(ArgumentError("table with array rows requires `headers`"))
        headers = String.(tbl["headers"])

        out = NamedTuple[]
        for (k, row) in enumerate(rows)
            # skip empty or comment rows
            if isempty(row) || (isa(row[1], String) && startswith(row[1], "#"))
                continue
            end
            # allow missing trailing columns (fill with nothing)
            if length(row) < length(headers)
                row = vcat(row, fill(nothing, length(headers) - length(row)))
            end
            if length(row) > length(headers)
                @warn "Skipping row $k in table: has $(length(row)) values, expected $(length(headers)). Row: $row"
                continue
            end
            nt = NamedTuple{Tuple(Symbol.(headers))}(Tuple(row))
            push!(out, nt)
        end
        return out
    end
end

"""
        load_sys_struct_from_yaml(yaml_path::AbstractString; system_name="from_yaml", set=nothing)

Build a `SystemStructure` from a component-based structural YAML file.

# Expected top-level blocks
- `points`: table with headers `[id,x,y,z,type,mass,body_damping,world_damping]`
  - `type`: STATIC, DYNAMIC, WING, or QUASI_STATIC

- `segments`: table with one of two formats:
  - Direct format: `[id,point_i,point_j,type,l0,diameter_mm,axial_stiffness,axial_damping,compression_frac]`
  - Named format: `[name,point_i,point_j]` (requires `segment_properties` block)

- `segment_properties`: (optional) table with headers `[name,type,l0,diameter_mm,axial_stiffness,axial_damping,compression_frac]`
  - Used with named segment format for shared properties across symmetric segments

- `pulleys`: table with headers `[id,segment_i,segment_j,type]`
  - `type`: DYNAMIC or QUASI_STATIC

- `groups`: (optional) table with headers `[id,point_ids,gamma,type,reference_chord_frac]`
- `tethers`: (optional) table with headers `[id,segment_ids,ground_point_id]`
- `winches`: (optional) table with headers `[id,tether_ids]`
- `wings`: (optional, typically from VSM configuration)
- `transforms`: (optional, typically from settings)
"""
function load_sys_struct_from_yaml(yaml_path::AbstractString; system_name="from_yaml", set=nothing())
    data = YAML.load_file(yaml_path)

    # Use provided settings or fall back to base settings
    set === nothing && (set = load_settings("base"))

    # Parse string types to enums
    function parse_dynamics_type(s::String)
        s_upper = uppercase(s)
        s_upper == "STATIC" && return STATIC
        s_upper == "DYNAMIC" && return DYNAMIC
        s_upper == "WING" && return WING
        s_upper == "QUASI_STATIC" && return QUASI_STATIC
        error("Unknown DynamicsType: $s")
    end

    function parse_segment_type(s::String)
        s_upper = uppercase(s)
        s_upper == "POWER_LINE" && return POWER_LINE
        s_upper == "STEERING_LINE" && return STEERING_LINE
        s_upper == "BRIDLE" && return BRIDLE
        error("Unknown SegmentType: $s")
    end

    # Load points
    points = Point[]
    yaml_id_to_idx = Dict{Int,Int}()  # Map from YAML IDs to 1-based indices

    if haskey(data, "points")
        point_rows = parse_table(data["points"])
        for (i, row) in enumerate(point_rows)
            yaml_id = Int(row.id)
            pos = [Float64(row.x), Float64(row.y), Float64(row.z)]
            ptype = parse_dynamics_type(String(row.type))
            mass = Float64(row.mass)
            body_damping = Float64(row.body_damping)
            world_damping = Float64(row.world_damping)

            push!(points, Point(
                i,  # Use 1-based index
                pos,
                ptype;
                mass = mass,
                body_frame_damping = body_damping,
                world_frame_damping = world_damping,
                transform_idx = Int16(1)
            ))
            points[end].pos_w .= points[end].pos_cad
            points[end].vel_w .= 0.0

            yaml_id_to_idx[yaml_id] = i
        end
    end

    isempty(points) && error("No points found in YAML file $(yaml_path).")

    # Build property tables for reference resolution
    property_tables = Dict{String, Dict{String, NamedTuple}}()

    # Load materials table
    if haskey(data, "materials")
        materials_dict = Dict{String, NamedTuple}()
        material_rows = parse_table(data["materials"])
        for row in material_rows
            name = String(row.name)
            materials_dict[name] = row
        end
        property_tables["materials"] = materials_dict
    end

    # Load elements table (old-style segment properties)
    if haskey(data, "elements")
        elements_dict = Dict{String, NamedTuple}()
        element_rows = parse_table(data["elements"])
        for row in element_rows
            name = String(row.name)
            elements_dict[name] = row
        end
        property_tables["elements"] = elements_dict
    end

    # Load segment_properties table (for backward compatibility)
    if haskey(data, "segment_properties")
        segment_props_dict = Dict{String, NamedTuple}()
        prop_rows = parse_table(data["segment_properties"])
        for row in prop_rows
            name = String(row.name)
            segment_props_dict[name] = row
        end
        property_tables["segment_properties"] = segment_props_dict
    end

    # Load segments
    segments = Segment[]
    if haskey(data, "segments")
        segment_rows = parse_table(data["segments"])

        for (i, row) in enumerate(segment_rows)
            # Resolve any references in the row
            resolved_row = resolve_references(row, property_tables)

            # Convert to mutable dict for derived property calculation
            props = Dict{Symbol, Any}(pairs(resolved_row))

            # Smart detection: check if axial_stiffness is directly provided as a number
            # or if we need to calculate it from material properties
            if haskey(props, :axial_stiffness) && props[:axial_stiffness] isa Number
                # Direct stiffness specified - use it as-is
                axial_stiffness = Float64(props[:axial_stiffness])
                # Check if damping is also provided, otherwise default to 0
                if haskey(props, :axial_damping) && props[:axial_damping] !== nothing
                    axial_damping = Float64(props[:axial_damping])
                else
                    axial_damping = 0.0
                end
            else
                # Material-based or needs calculation from material properties
                calculate_derived_properties!(props)
                if !haskey(props, :axial_stiffness) || props[:axial_stiffness] === nothing
                    error("Segment $i: Unable to determine axial_stiffness. Either provide it directly or specify material properties.")
                end
                axial_stiffness = Float64(props[:axial_stiffness])
                axial_damping = Float64(props[:axial_damping])
            end

            # Extract required fields
            yaml_point_i = Int(props[:point_i])
            yaml_point_j = Int(props[:point_j])
            seg_type = parse_segment_type(String(props[:type]))
            l0 = Float64(props[:l0])
            diameter_mm = Float64(props[:diameter_mm])

            # Handle compression_frac which might be in different positions
            if haskey(props, :compression_frac) && props[:compression_frac] !== nothing
                compression_frac = Float64(props[:compression_frac])
            else
                compression_frac = 0.0
            end

            # Map YAML point IDs to 1-based indices
            if !haskey(yaml_id_to_idx, yaml_point_i)
                error("Segment $i references unknown point ID $yaml_point_i")
            end
            if !haskey(yaml_id_to_idx, yaml_point_j)
                error("Segment $i references unknown point ID $yaml_point_j")
            end
            point_i = yaml_id_to_idx[yaml_point_i]
            point_j = yaml_id_to_idx[yaml_point_j]

            push!(segments, Segment(
                i,  # Use 1-based index
                set,
                (point_i, point_j),
                seg_type;
                l0 = l0,
                diameter_mm = diameter_mm,
                axial_stiffness = axial_stiffness,
                axial_damping = axial_damping,
                compression_frac = compression_frac
            ))
        end
    end

    # Load pulleys
    pulleys = Pulley[]
    if haskey(data, "pulleys")
        pulley_rows = parse_table(data["pulleys"])
        for (i, row) in enumerate(pulley_rows)
            segment_i = Int(row.segment_i)
            segment_j = Int(row.segment_j)
            ptype = parse_dynamics_type(String(row.type))

            # Validate segment indices
            if segment_i < 1 || segment_i > length(segments)
                error("Pulley $i references invalid segment index $segment_i (must be 1-$(length(segments)))")
            end
            if segment_j < 1 || segment_j > length(segments)
                error("Pulley $i references invalid segment index $segment_j (must be 1-$(length(segments)))")
            end

            push!(pulleys, Pulley(
                i,  # Use 1-based index
                (segment_i, segment_j),
                ptype
            ))
        end
    end

    # Load groups (optional, for deformable wings)
    groups = Group[]
    if haskey(data, "groups") && haskey(data["groups"], "data") && data["groups"]["data"] !== nothing && !isempty(data["groups"]["data"])
        # Groups would be loaded here if defined
        # This requires VSM wing instance and gamma values
        @warn "Groups defined in YAML but loading not yet implemented"
    end

    # Load tethers (optional)
    tethers = Tether[]
    if haskey(data, "tethers") && haskey(data["tethers"], "data") && data["tethers"]["data"] !== nothing && !isempty(data["tethers"]["data"])
        # Tethers would be loaded here if defined
        @warn "Tethers defined in YAML but loading not yet implemented"
    end

    # Load winches (optional)
    winches = Winch[]
    if haskey(data, "winches") && haskey(data["winches"], "data") && data["winches"]["data"] !== nothing && !isempty(data["winches"]["data"])
        # Winches would be loaded here if defined
        @warn "Winches defined in YAML but loading not yet implemented"
    end

    # Parse wing type
    function parse_wing_type(s::String)
        s_upper = uppercase(s)
        s_upper == "REFINE" && return REFINE
        s_upper == "QUATERNION" && return QUATERNION
        error("Unknown WingType: $s")
    end

    # Load wings (optional)
    wings = AbstractWing[]
    wings_defined = false
    if haskey(data, "wings") && haskey(data["wings"], "data") && data["wings"]["data"] !== nothing && !isempty(data["wings"]["data"])
        wing_rows = parse_table(data["wings"])
        wings_defined = true

        for row in wing_rows
            wing_id = Int(row.id)
            wing_type = parse_wing_type(String(row.type))

            # Get point IDs for this wing
            yaml_point_ids = Vector{Int}(row.point_ids)
            wing_point_idxs = [yaml_id_to_idx[pid] for pid in yaml_point_ids]
            wing_point_objs = [points[idx] for idx in wing_point_idxs]

            # Validate that the points exist and are WING type
            for idx in wing_point_idxs
                if points[idx].type != WING
                    error("Wing $wing_id references point $idx which is not of type WING")
                end
                # Set wing_idx on the point
                points[idx] = Point(
                    points[idx].idx,
                    points[idx].pos_cad,
                    points[idx].type;
                    wing_idx = Int16(wing_id),
                    mass = points[idx].mass,
                    body_frame_damping = points[idx].body_frame_damping,
                    world_frame_damping = points[idx].world_frame_damping,
                    transform_idx = points[idx].transform_idx
                )
                points[idx].pos_w .= points[idx].pos_cad
                points[idx].vel_w .= 0.0
            end

            # Create VSM wing from settings
            @info "Creating VSM wing $wing_id of type $wing_type..."
            vsm_wing = VortexStepMethod.Wing(set; prn=false)
            vsm_aero = VortexStepMethod.BodyAerodynamics([vsm_wing])
            vsm_solver = VortexStepMethod.Solver(vsm_aero; solver_type=VortexStepMethod.NONLIN, atol=2e-8, rtol=2e-8)

            # Create wing based on type
            if wing_type == REFINE
                # REFINE wing: Direct panel forces to structural points
                # Identify wing segments (LE/TE pairs)
                wing_segments = identify_wing_segments(wing_point_objs)

                # Build panel-to-point force lumping mapping
                point_to_panels = build_point_to_panel_mapping(wing_point_objs, vsm_aero)

                # Identify reference points for orientation tracking
                # Read from YAML if provided, otherwise auto-detect
                if haskey(row, :x_ids) && !isnothing(row.x_ids)
                    # User-specified X reference points (chord direction)
                    yaml_x_ids = Vector{Int}(row.x_ids)
                    x_ref_points = Tuple(Int16(yaml_id_to_idx[xid]) for xid in yaml_x_ids)
                    @info "  Using user-specified X reference points: YAML IDs $yaml_x_ids -> internal $(x_ref_points)"
                else
                    # Auto-detect: Use first wing segment for X direction (chord)
                    x_ref_points = wing_segments[1]  # First LE-TE pair defines forward (X) direction
                    @info "  Auto-detected X reference points: $(x_ref_points)"
                end

                if haskey(row, :y_ids) && !isnothing(row.y_ids)
                    # User-specified Y reference points (span direction)
                    yaml_y_ids = Vector{Int}(row.y_ids)
                    y_ref_points = Tuple(Int16(yaml_id_to_idx[yid]) for yid in yaml_y_ids)
                    @info "  Using user-specified Y reference points: YAML IDs $yaml_y_ids -> internal $(y_ref_points)"
                else
                    # Auto-detect: Use two LE points for span (Y) direction
                    mid_idx = length(wing_segments) ÷ 2
                    y_ref_points = (wing_segments[1][1], wing_segments[mid_idx][1])  # Two LE points for span (Y) direction
                    @info "  Auto-detected Y reference points: $(y_ref_points)"
                end

                # Create REFINE VSMWing (no groups)
                wing = VSMWing(
                    wing_id,
                    vsm_aero,
                    vsm_wing,
                    vsm_solver,
                    Int16[],  # No groups for REFINE
                    vsm_wing.R_cad_body,
                    vsm_wing.T_cad_body;
                    transform_idx=1,
                    y_damping=150.0,
                    wing_type=REFINE,
                    point_to_panels=point_to_panels,
                    wing_segments=wing_segments,
                    x_ref_points=x_ref_points,
                    y_ref_points=y_ref_points
                )
                @info "  ✓ REFINE wing created: $(length(wing_point_objs)) structural points, $(length(vsm_aero.panels)) VSM panels"
            elseif wing_type == QUATERNION
                # QUATERNION wing: Rigid body with group dynamics
                # For now, assume no groups (would need to be specified in YAML)
                group_idxs = Int16[]

                wing = VSMWing(
                    wing_id,
                    vsm_aero,
                    vsm_wing,
                    vsm_solver,
                    group_idxs,
                    vsm_wing.R_cad_body,
                    vsm_wing.T_cad_body;
                    transform_idx=1,
                    y_damping=150.0,
                    wing_type=QUATERNION
                )
                @info "  ✓ QUATERNION wing created with $(length(group_idxs)) groups"
            else
                error("Unsupported wing type: $wing_type")
            end

            push!(wings, wing)
        end
    end

    # Load transforms (optional)
    transforms = Transform[]
    if haskey(data, "transforms") && haskey(data["transforms"], "data") && data["transforms"]["data"] !== nothing && !isempty(data["transforms"]["data"])
        transform_rows = parse_table(data["transforms"])

        for row in transform_rows
            transform_id = Int16(row.id)
            elevation = Float64(row.elevation)
            azimuth = Float64(row.azimuth)
            heading = Float64(row.heading)

            # Parse optional base_pos (can be provided as separate x,y,z or as array)
            base_pos = if haskey(row, :base_pos) && !isnothing(row.base_pos)
                # Array format: base_pos: [x, y, z]
                KVec3(row.base_pos...)
            elseif haskey(row, :base_pos_x) && !isnothing(row.base_pos_x)
                # Separate components: base_pos_x, base_pos_y, base_pos_z
                KVec3(Float64(row.base_pos_x), Float64(row.base_pos_y), Float64(row.base_pos_z))
            else
                nothing
            end

            # Parse optional indices
            base_point_idx = haskey(row, :base_point_idx) && !isnothing(row.base_point_idx) ?
                Int16(yaml_id_to_idx[Int(row.base_point_idx)]) : nothing
            wing_idx = haskey(row, :wing_idx) && !isnothing(row.wing_idx) ?
                Int16(row.wing_idx) : nothing
            rot_point_idx = haskey(row, :rot_point_idx) && !isnothing(row.rot_point_idx) ?
                Int16(yaml_id_to_idx[Int(row.rot_point_idx)]) : nothing
            base_transform_idx = haskey(row, :base_transform_idx) && !isnothing(row.base_transform_idx) ?
                Int16(row.base_transform_idx) : nothing

            # Create Transform
            transform = Transform(
                transform_id,
                elevation,
                azimuth,
                heading;
                base_point_idx = base_point_idx,
                base_pos = base_pos,
                base_transform_idx = base_transform_idx,
                wing_idx = wing_idx,
                rot_point_idx = rot_point_idx
            )

            push!(transforms, transform)
            @info "  ✓ Transform $transform_id created: elevation=$(elevation)°, azimuth=$(azimuth)°, heading=$(heading)°"
        end
    end

    # If no wings are provided, convert WING type points to STATIC with warning
    if !wings_defined
        wing_points = findall(p -> p.type == WING, points)
        if !isempty(wing_points)
            @warn "No wings provided but $(length(wing_points)) WING type points found. Converting to STATIC."
            for idx in wing_points
                points[idx] = Point(
                    points[idx].idx,
                    points[idx].pos_cad,
                    STATIC;
                    mass = points[idx].mass,
                    body_frame_damping = points[idx].body_frame_damping,
                    world_frame_damping = points[idx].world_frame_damping,
                    transform_idx = points[idx].transform_idx
                )
                points[idx].pos_w .= points[idx].pos_cad
                points[idx].vel_w .= 0.0
            end
        end
    end

    return SystemStructure(system_name, set; points, groups, segments, pulleys, tethers, winches, wings, transforms)
end
