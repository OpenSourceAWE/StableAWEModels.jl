using YAML
using Logging
using LinearAlgebra
using VortexStepMethod, KiteUtils

# Default diameter for wing section segments (in mm)
DEFAULT_WING_DIAMETER_MM = 1.0

# ------------------ helpers ------------------

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
    if haskey(data, "points")
        point_rows = parse_table(data["points"])
        for row in point_rows
            idx = Int(row.id)
            pos = [Float64(row.x), Float64(row.y), Float64(row.z)]
            ptype = parse_dynamics_type(String(row.type))
            mass = Float64(row.mass)
            body_damping = Float64(row.body_damping)
            world_damping = Float64(row.world_damping)

            push!(points, Point(
                idx,
                pos,
                ptype;
                mass = mass,
                body_frame_damping = body_damping,
                world_frame_damping = world_damping,
                transform_idx = Int16(0)
            ))
            points[end].pos_w .= points[end].pos_cad
            points[end].vel_w .= 0.0
        end
    end

    isempty(points) && error("No points found in YAML file $(yaml_path).")

    # Load segment properties dictionary (if using named segments)
    segment_props_dict = Dict{String,Dict{Symbol,Any}}()
    if haskey(data, "segment_properties")
        prop_rows = parse_table(data["segment_properties"])
        for row in prop_rows
            name = String(row.name)
            segment_props_dict[name] = Dict{Symbol,Any}(
                :type => String(row.type),
                :l0 => Float64(row.l0),
                :diameter_mm => Float64(row.diameter_mm),
                :axial_stiffness => Float64(row.axial_stiffness),
                :axial_damping => Float64(row.axial_damping),
                :compression_frac => Float64(row.compression_frac)
            )
        end
    end

    # Load segments
    segments = Segment[]
    if haskey(data, "segments")
        segment_rows = parse_table(data["segments"])
        segment_counter = 1

        for row in segment_rows
            # Check if this is direct format (has 'id') or named format (has 'name')
            if haskey(row, :id)
                # Direct format: all properties specified inline
                idx = Int(row.id)
                point_i = Int(row.point_i)
                point_j = Int(row.point_j)
                seg_type = parse_segment_type(String(row.type))
                l0 = Float64(row.l0)
                diameter_mm = Float64(row.diameter_mm)
                axial_stiffness = Float64(row.axial_stiffness)
                axial_damping = Float64(row.axial_damping)
                compression_frac = Float64(row.compression_frac)
            elseif haskey(row, :name)
                # Named format: look up properties from segment_properties
                name = String(row.name)
                point_i = Int(row.point_i)
                point_j = Int(row.point_j)

                if !haskey(segment_props_dict, name)
                    error("Segment named '$name' not found in segment_properties")
                end

                props = segment_props_dict[name]
                idx = segment_counter
                seg_type = parse_segment_type(props[:type])
                l0 = props[:l0]
                diameter_mm = props[:diameter_mm]
                axial_stiffness = props[:axial_stiffness]
                axial_damping = props[:axial_damping]
                compression_frac = props[:compression_frac]
            else
                error("Segment row must have either 'id' or 'name' field")
            end

            push!(segments, Segment(
                idx,
                set,
                (point_i, point_j),
                seg_type;
                l0 = l0,
                diameter_mm = diameter_mm,
                axial_stiffness = axial_stiffness,
                axial_damping = axial_damping,
                compression_frac = compression_frac
            ))
            segment_counter += 1
        end
    end

    # Load pulleys
    pulleys = Pulley[]
    if haskey(data, "pulleys")
        pulley_rows = parse_table(data["pulleys"])
        for row in pulley_rows
            idx = Int(row.id)
            segment_i = Int(row.segment_i)
            segment_j = Int(row.segment_j)
            ptype = parse_dynamics_type(String(row.type))

            push!(pulleys, Pulley(
                idx,
                (segment_i, segment_j),
                ptype
            ))
        end
    end

    # Load groups (optional, for deformable wings)
    groups = Group[]
    if haskey(data, "groups") && haskey(data["groups"], "data") && !isempty(data["groups"]["data"])
        # Groups would be loaded here if defined
        # This requires VSM wing instance and gamma values
        @warn "Groups defined in YAML but loading not yet implemented"
    end

    # Load tethers (optional)
    tethers = Tether[]
    if haskey(data, "tethers") && haskey(data["tethers"], "data") && !isempty(data["tethers"]["data"])
        # Tethers would be loaded here if defined
        @warn "Tethers defined in YAML but loading not yet implemented"
    end

    # Load winches (optional)
    winches = Winch[]
    if haskey(data, "winches") && haskey(data["winches"], "data") && !isempty(data["winches"]["data"])
        # Winches would be loaded here if defined
        @warn "Winches defined in YAML but loading not yet implemented"
    end

    # Wings and transforms typically come from VSM and settings
    wings = AbstractWing[]
    transforms = Transform[]

    return SystemStructure(system_name, set; points, groups, segments, pulleys, tethers, winches, wings, transforms)
end
