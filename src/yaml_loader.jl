using YAML
using Logging
using LinearAlgebra
using VortexStepMethod, KiteUtils

# Default diameter for wing section segments (in mm)
DEFAULT_WING_DIAMETER_MM = 1.0

# ------------------ helpers ------------------

"""
    get_field_or_nothing(::Type{T}, row::NamedTuple,
                         field::Symbol) where T

Convert field to type T if present, otherwise return nothing.

# Examples
```julia
get_field_or_nothing(Int16, row, :idx)  # -> Int16 or nothing
get_field_or_nothing(Tuple{Int16,Int16}, row, :pair)
    # -> (Int16, Int16) or nothing
```
"""
function get_field_or_nothing(::Type{T}, row::NamedTuple,
                               field::Symbol) where T
    if !haskey(row, field) || isnothing(row[field])
        return nothing
    end
    return convert_to_type(T, row[field])
end

"""
    convert_to_type(::Type{T}, value) where T

Convert value to type T. Handles special cases like Tuples.
"""
convert_to_type(::Type{T}, value) where T = T(value)

# Special handling for Tuple types
function convert_to_type(
        ::Type{Tuple{T,T}}, value) where T
    vec = Vector{Int}(value)
    return (T(vec[1]), T(vec[2]))
end

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

    # Check format: if first row is a Dict,
    # use dict format; if Array, use header format
    first_row = first(rows)

    if first_row isa AbstractDict
        # Dict format: each row is already a dict with named keys
        # Convert each dict to a NamedTuple
        out = NamedTuple[]
        for row in rows
            nt = NamedTuple{Tuple(Symbol.(keys(row)))}(
                Tuple(values(row)))
            push!(out, nt)
        end
        return out
    else
        # Array format: requires headers
        haskey(tbl, "headers") ||
            throw(ArgumentError(
                "table with array rows requires `headers`"))
        headers = String.(tbl["headers"])

        out = NamedTuple[]
        for (k, row) in enumerate(rows)
            # skip empty or comment rows
            if isempty(row) ||
               (isa(row[1], String) && startswith(row[1], "#"))
                continue
            end
            # allow missing trailing columns (fill with nothing)
            if length(row) < length(headers)
                row = vcat(row, fill(nothing,
                    length(headers) - length(row)))
            end
            if length(row) > length(headers)
                @warn "Skipping row $k: has $(length(row)) " *
                      "values, expected $(length(headers))."
                continue
            end
            nt = NamedTuple{Tuple(Symbol.(headers))}(Tuple(row))
            push!(out, nt)
        end
        return out
    end
end

"""
    call_yaml_constructor(Constructor, row::NamedTuple,
        args_spec, kwargs_spec; mappings=Dict())

Generic YAML-to-constructor caller. Extracts positional
args and kwargs from YAML row and calls constructor.

# Arguments
- `Constructor`: Constructor function to call
- `row::NamedTuple`: Parsed YAML row
- `args_spec::Vector{Symbol}`: Names for positional args
- `kwargs_spec::Vector{Symbol}`: Names for kwargs

# Keyword Arguments
- `mappings::Dict{Symbol, Function}`: Mapping functions
  that take the row and return the arg value

# Example
```julia
row = (idx=1, x=0.0, y=0.0, z=0.0, type="STATIC")
point = call_yaml_constructor(Point, row,
    [:idx, :pos_cad, :type],  # positional args
    [:mass, :wing_idx];       # kwargs
    mappings=Dict(
        :pos_cad => r -> [Float64(r.x),
            Float64(r.y), Float64(r.z)],
        :type => r -> parse_dynamics_type(
            String(r.type))
    ))
```
"""
function call_yaml_constructor(
        Constructor,
        row::NamedTuple,
        args_spec::Vector{Symbol},
        kwargs_spec::Vector;
        mappings::Dict{Symbol, <:Function}=
            Dict{Symbol, Function}())

    # Extract positional arguments
    args = []
    for arg_name in args_spec
        if haskey(mappings, arg_name)
            push!(args, mappings[arg_name](row))
        elseif haskey(row, arg_name)
            push!(args, row[arg_name])
        else
            error("Missing required arg $arg_name")
        end
    end

    # Extract keyword arguments (only if present)
    kwargs = Dict{Symbol, Any}()
    for kwarg_name in kwargs_spec
        if haskey(mappings, kwarg_name)
            kwargs[kwarg_name] = mappings[kwarg_name](row)
        elseif haskey(row, kwarg_name) &&
               !isnothing(row[kwarg_name])
            kwargs[kwarg_name] = row[kwarg_name]
        end
    end

    return Constructor(args...; kwargs...)
end

"""
        load_sys_struct_from_yaml(yaml_path::AbstractString; system_name="from_yaml", set=nothing)

Build a `SystemStructure` from a component-based structural YAML file.

**IMPORTANT**: All indices (points, segments, etc.) must be sequential
starting from 1 with no gaps.

# Expected top-level blocks
- `points`: table with headers `[id,x,y,z,type,mass,body_damping,world_damping]`
  - `type`: STATIC, DYNAMIC, WING, or QUASI_STATIC
  - `id` must be sequential: 1, 2, 3, ...

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
function load_sys_struct_from_yaml(yaml_path::AbstractString; system_name="from_yaml", set=nothing(), ignore_l0::Bool=false, wing_type::Union{Nothing,WingType}=nothing, vsm_set=nothing)
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
        for (i, row) in enumerate(point_rows)
            # Verify sequential indexing
            @assert Int(row.idx) == i
                "Point indices must be sequential: " *
                "expected $i, got $(row.idx)"

            # Create Point using generic constructor
            point = call_yaml_constructor(Point, row,
                [:idx, :pos_cad, :type],
                [:wing_idx, :transform_idx, :mass,
                 :body_frame_damping, :world_frame_damping];
                mappings=Dict(
                    :pos_cad => r -> KVec3(r.pos_cad...),
                    :type => r -> parse_dynamics_type(
                        String(r.type))
                ))

            point.pos_w .= point.pos_cad
            point.vel_w .= 0.0
            push!(points, point)
        end
    end

    isempty(points) &&
        error("No points found in YAML file $(yaml_path).")

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
            # Resolve references and calculate derived properties
            resolved_row = resolve_references(row, property_tables)
            props = Dict{Symbol, Any}(pairs(resolved_row))
            calculate_derived_properties!(props)

            # Convert back to NamedTuple for constructor
            resolved_row = NamedTuple(props)

            # Create Segment using generic constructor
            segment = call_yaml_constructor(Segment, resolved_row,
                [:idx, :set, :point_idxs, :type],
                [:l0, :diameter_mm, :axial_stiffness,
                 :axial_damping, :compression_frac];
                mappings=Dict(
                    :set => r -> set,
                    :point_idxs => r -> (Int(r.point_i),
                        Int(r.point_j)),
                    :type => r -> parse_segment_type(
                        String(r.type))
                ))

            push!(segments, segment)
        end
    end

    # Load pulleys
    pulleys = Pulley[]
    if haskey(data, "pulleys")
        pulley_rows = parse_table(data["pulleys"])
        for (i, row) in enumerate(pulley_rows)
            pulley = call_yaml_constructor(Pulley, row,
                [:idx, :segment_idxs, :type],
                [];
                mappings=Dict(
                    :segment_idxs => r -> (Int(r.segment_i),
                        Int(r.segment_j)),
                    :type => r -> parse_dynamics_type(String(r.type))
                ))
            push!(pulleys, pulley)
        end
    end

    # Load groups (optional, for deformable wings)
    groups = Group[]
    if haskey(data, "groups") &&
       haskey(data["groups"], "data") &&
       data["groups"]["data"] !== nothing &&
       !isempty(data["groups"]["data"])
        group_rows = parse_table(data["groups"])

        for (i, row) in enumerate(group_rows)
            group = call_yaml_constructor(Group, row,
                [:idx, :point_idxs, :gamma, :type,
                 :moment_frac],
                [:damping];
                mappings=Dict(
                    :point_idxs => r ->
                        Vector{Int}(r.point_ids),
                    :type => r -> parse_dynamics_type(
                        String(r.type))
                ))
            push!(groups, group)
        end
    end

    # Load tethers (optional)
    tethers = Tether[]
    if haskey(data, "tethers") &&
       haskey(data["tethers"], "data") &&
       data["tethers"]["data"] !== nothing &&
       !isempty(data["tethers"]["data"])
        tether_rows = parse_table(data["tethers"])
        for (i, row) in enumerate(tether_rows)
            tether = call_yaml_constructor(Tether, row,
                [:idx, :segment_idxs, :winch_idx],
                [];
                mappings=Dict(
                    :segment_idxs => r -> isnothing(r.segment_idxs) ? Int16[] : Vector{Int16}(r.segment_idxs),
                    :winch_idx => r -> isnothing(r.winch_idx) ? Int16(0) : Int16(r.winch_idx)
                ))
            push!(tethers, tether)
        end
    end

    # Load winches (optional)
    winches = Winch[]
    if haskey(data, "winches") &&
       haskey(data["winches"], "data") &&
       data["winches"]["data"] !== nothing &&
       !isempty(data["winches"]["data"])
        winch_rows = parse_table(data["winches"])
        for (i, row) in enumerate(winch_rows)
            winch = call_yaml_constructor(Winch, row,
                [:idx, :set, :tether_idxs],
                [:tether_len, :tether_vel, :brake];
                mappings=Dict(
                    :set => r -> set,
                    :tether_idxs => r -> Vector{Int16}(r.tether_idxs)
                ))
            push!(winches, winch)
        end
    end

    # Parse wing type
    function parse_wing_type(s::String)
        s_upper = uppercase(s)
        s_upper == "REFINE" && return REFINE
        s_upper == "QUATERNION" && return QUATERNION
        error("Unknown WingType: $s")
    end

    # Parse reference points (can be single point or vector of points to average)
    # Examples: [12, 13], [12, [13, 14]], [[11, 12], [13, 14]]
    function parse_ref_points(row, field)
        !hasfield(typeof(row), field) && return nothing
        val = getfield(row, field)
        val === nothing && return nothing

        # Parse [a, b] or [a, [b, c]]
        @assert length(val) == 2 "ref_points must have 2 elements"

        p1 = val[1] isa Vector ? Int16.(val[1]) : Int16(val[1])
        p2 = val[2] isa Vector ? Int16.(val[2]) : Int16(val[2])

        return (p1, p2)
    end

    # Load wings (optional)
    wings = AbstractWing[]
    if haskey(data, "wings") &&
       haskey(data["wings"], "data") &&
       data["wings"]["data"] !== nothing &&
       !isempty(data["wings"]["data"])
        wing_rows = parse_table(data["wings"])

        for (i, row) in enumerate(wing_rows)
            # Use provided wing_type parameter or parse from YAML
            wt = isnothing(wing_type) ? parse_wing_type(String(row.type)) : wing_type

            # Build kwargs based on wing type
            if wt == REFINE
                # REFINE wings need z_ref_points, y_ref_points, origin_idx
                wing = call_yaml_constructor(VSMWing, row,
                    [:idx, :set, :group_idxs, :vsm_set],
                    [:transform_idx, :y_damping, :wing_type,
                     :z_ref_points, :y_ref_points, :origin_idx];
                    mappings=Dict(
                        :set => r -> set,
                        :group_idxs => r -> Int16[],
                        :vsm_set => r -> vsm_set,
                        :wing_type => r -> wt,
                        :z_ref_points => r ->
                            parse_ref_points(r, :z_ref_points),
                        :y_ref_points => r ->
                            parse_ref_points(r, :y_ref_points),
                        :origin_idx => r ->
                            get_field_or_nothing(Int16, r, :origin_idx)
                    ))
            else  # QUATERNION
                # QUATERNION wings don't use these fields
                wing = call_yaml_constructor(VSMWing, row,
                    [:idx, :set, :group_idxs, :vsm_set],
                    [:transform_idx, :y_damping, :wing_type, :aero_z_offset];
                    mappings=Dict(
                        :set => r -> set,
                        :group_idxs => r -> Int16[],
                        :vsm_set => r -> vsm_set,
                        :wing_type => r -> wt
                    ))
            end
            push!(wings, wing)
        end
    end

    # Load transforms (optional)
    transforms = Transform[]
    if haskey(data, "transforms") &&
       haskey(data["transforms"], "data") &&
       data["transforms"]["data"] !== nothing &&
       !isempty(data["transforms"]["data"])
        transform_rows = parse_table(data["transforms"])

        for row in transform_rows
            transform = call_yaml_constructor(Transform, row,
                [:idx, :elevation, :azimuth, :heading],
                [:base_point_idx, :base_pos,
                 :base_transform_idx, :wing_idx, :rot_point_idx];
                mappings=Dict(
                    :elevation => r -> deg2rad(r.elevation),
                    :azimuth => r -> deg2rad(r.azimuth),
                    :heading => r -> deg2rad(r.heading),
                    :base_pos => r -> KVec3(r.base_pos...),
                    :base_point_idx => r ->
                        Int16(r.base_point_idx),
                    :rot_point_idx => r ->
                        get_field_or_nothing(Int16, r,
                            :rot_point_idx),
                    :wing_idx => r ->
                        get_field_or_nothing(Int16, r, :wing_idx)
                ))
            push!(transforms, transform)
            elev_deg = rad2deg(transform.elevation)
            azim_deg = rad2deg(transform.azimuth)
            head_deg = rad2deg(transform.heading)
            @info "  ✓ Transform $(transform.idx) created: " *
                  "elevation=$(elev_deg)°, " *
                  "azimuth=$(azim_deg)°, heading=$(head_deg)°"
        end
    end

    # SystemStructure constructor now handles WING→STATIC
    # conversion when no wings are defined
    return SystemStructure(system_name, set; points, groups,
        segments, pulleys, tethers, winches, wings, transforms, ignore_l0, vsm_set)
end
