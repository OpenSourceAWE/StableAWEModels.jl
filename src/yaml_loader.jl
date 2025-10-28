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
    haskey(tbl, "headers") || throw(ArgumentError("table is missing `headers`"))
    haskey(tbl, "data")    || throw(ArgumentError("table is missing `data`"))
    headers = String.(tbl["headers"])
    rows = tbl["data"]
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
                transform_idx = Int16(0)
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

    # Wings and transforms typically come from VSM and settings
    wings = AbstractWing[]
    transforms = Transform[]

    # If no wings are provided, convert WING type points to DYNAMIC with warning
    if isempty(wings)
        wing_points = findall(p -> p.type == WING, points)
        if !isempty(wing_points)
            @warn "No wings provided but $(length(wing_points)) WING type points found. Converting to DYNAMIC."
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
