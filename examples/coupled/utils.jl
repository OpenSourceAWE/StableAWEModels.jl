"""
Shared utility functions for coupled examples.
"""

using CSV, DataFrames
using LinearAlgebra
using UnPack
using Rotations

"""
    load_extra_points(csv_path, sys_struct; body_offset)

Load extra points from CSV and transform from camera frame to simulation frame.
Returns (transformed_points, groups) where groups is Vector of (group_name, indices).

CSV has columns: group, idx_in_group, x, y, z.
Alignment: CSV strut3/strut4 LE centers align with sim points 10, 12.
"""
function load_extra_points(csv_path::String, sys_struct; body_offset=[0.3, 0.0, 0.2])
    df = CSV.read(csv_path, DataFrame)

    # CSV strut centers: strut3[1]/strut4[1] are at TE, [end] are at LE
    strut3 = [[r.x, r.y, r.z] for r in eachrow(df) if r.group == "strut3"]
    strut4 = [[r.x, r.y, r.z] for r in eachrow(df) if r.group == "strut4"]
    csv_le_center = (strut3[end] + strut4[end]) / 2
    csv_te_center = (strut3[1] + strut4[1]) / 2

    # Sim reference: points 10, 12 (center LE)
    sim_p10 = collect(sys_struct.points[10].pos_w)
    sim_p12 = collect(sys_struct.points[12].pos_w)
    sim_le_center = (sim_p10 + sim_p12) / 2

    # Direction vectors
    csv_span = normalize(strut4[end] - strut3[end])

    # CSV basis: y=spanwise, z from wing center geometry, x from cross
    csv_y = csv_span
    csv_wing_center = (csv_le_center + csv_te_center) / 2
    csv_z = normalize(csv_wing_center - csv_y * 0.84/2)
    csv_x = cross(csv_y, csv_z)

    # Sim basis: directly from wing rotation matrix
    R_b_w = sys_struct.wings[1].R_b_w
    sim_x = R_b_w[:, 1]
    sim_y = R_b_w[:, 2]
    sim_z = R_b_w[:, 3]

    # Rotation: R * csv_basis = sim_basis
    csv_basis = hcat(csv_x, csv_y, csv_z)
    sim_basis = hcat(sim_x, sim_y, sim_z)
    R = sim_basis * csv_basis'

    # Translation: align LE centers
    T = sim_le_center - R * csv_le_center + R_b_w * body_offset

    # Transform all points (including camera origin marker at zeros)
    all_pts = [[row.x, row.y, row.z] for row in eachrow(df)]
    push!(all_pts, zeros(3))
    transformed = [Tuple(R * p + T) for p in all_pts]

    # Build group indices (1-based)
    groups = Vector{Tuple{String, Vector{Int}}}()
    current_group = ""
    current_indices = Int[]
    for (i, row) in enumerate(eachrow(df))
        if row.group != current_group
            if !isempty(current_indices)
                push!(groups, (current_group, copy(current_indices)))
            end
            current_group = row.group
            current_indices = [i]
        else
            push!(current_indices, i)
        end
    end
    if !isempty(current_indices)
        push!(groups, (current_group, current_indices))
    end

    return transformed, groups
end

"""
    plot_body_frame_local(sys_struct; extra_points, extra_groups, dir)

Plot wing points in 2D body frame. Extra points connected per-strut and LE.
No depth coloring, circle markers, thicker grey segments.

# Arguments
- `sys_struct`: System structure to plot
- `extra_points`: Optional vector of (x,y,z) tuples from load_extra_points
- `extra_groups`: Optional groups from load_extra_points
- `dir::Symbol`: Viewing direction (:side, :front, or :top)
"""
function plot_body_frame_local(sys_struct;
                               extra_points=nothing,
                               extra_groups=nothing,
                               dir::Symbol=:front,
                               point_size=10,
                               extra_point_size=8,
                               figsize=(800, 600))
    @unpack points, wings, segments = sys_struct

    # Update pos_b for REFINE wing points
    for wing in wings
        if wing.wing_type == SymbolicAWEModels.REFINE
            R_w_b = wing.R_b_w'
            for point in points
                if point.wing_idx == wing.idx
                    point.pos_b .= R_w_b * (point.pos_w - wing.pos_w)
                end
            end
        end
    end

    wing_points = [p for p in points if p.type == SymbolicAWEModels.WING]

    # Extract 2D coords based on viewing direction
    if dir == :top
        coords = [(p.pos_b[1], p.pos_b[2]) for p in wing_points]
        xlabel, ylabel = "x [m]", "y [m]"
    elseif dir == :side
        coords = [(p.pos_b[1], p.pos_b[3]) for p in wing_points]
        xlabel, ylabel = "x [m]", "z [m]"
    else  # :front
        coords = [(p.pos_b[2], p.pos_b[3]) for p in wing_points]
        xlabel, ylabel = "y [m]", "z [m]"
    end

    x_vals = [c[1] for c in coords]
    y_vals = [c[2] for c in coords]

    fig = Figure(size=figsize)
    ax = Axis(fig[1, 1]; xlabel, ylabel,
              title="Wing Points (Body Frame)", aspect=DataAspect())

    function get_2d(pos_b)
        if dir == :top
            return (pos_b[1], pos_b[2])
        elseif dir == :side
            return (pos_b[1], pos_b[3])
        else
            return (pos_b[2], pos_b[3])
        end
    end

    # Plot segments (only LE, struts, TE - skip diagonals 29-46) as thick grey lines
    wing_point_idxs = Set(p.idx for p in wing_points)
    for seg in segments
        # Skip diagonal spring segments (indices 29-46)
        if 29 <= seg.idx <= 46
            continue
        end
        from_idx, to_idx = seg.point_idxs
        if from_idx in wing_point_idxs && to_idx in wing_point_idxs
            p1 = points[from_idx]
            p2 = points[to_idx]
            c1 = get_2d(p1.pos_b)
            c2 = get_2d(p2.pos_b)
            lines!(ax, [c1[1], c2[1]], [c1[2], c2[2]];
                   color=(:gray, 0.5), linewidth=3)
        end
    end

    # Plot wing points (single color, circles)
    scatter!(ax, x_vals, y_vals;
             markersize=point_size, color=:blue, marker=:circle)

    # Add point labels
    for (i, p) in enumerate(wing_points)
        px, py = coords[i]
        away_x, away_y = 0.0, 0.0
        for (j, _) in enumerate(wing_points)
            if i != j
                ox, oy = coords[j]
                dx, dy = px - ox, py - oy
                dist = sqrt(dx^2 + dy^2) + 0.01
                away_x += dx / dist^2
                away_y += dy / dist^2
            end
        end
        away_len = sqrt(away_x^2 + away_y^2)
        if away_len > 0
            away_x /= away_len
            away_y /= away_len
        else
            away_x, away_y = 1.0, 1.0
        end
        offset = (12 * sign(away_x), 12 * sign(away_y))
        align_x = away_x >= 0 ? :left : :right
        align_y = away_y >= 0 ? :bottom : :top
        text!(ax, px, py; text=string(p.idx), fontsize=12,
              align=(align_x, align_y), offset=offset)
    end

    # Plot extra points with connections
    if !isnothing(extra_points) && !isnothing(extra_groups)
        wing = wings[1]
        R_w_b = wing.R_b_w'
        extra_body = [R_w_b * (collect(p) - wing.pos_w) for p in extra_points]

        if dir == :top
            extra_coords = [(p[1], p[2]) for p in extra_body]
        elseif dir == :side
            extra_coords = [(p[1], p[3]) for p in extra_body]
        else
            extra_coords = [(p[2], p[3]) for p in extra_body]
        end

        # Draw lines connecting points within each group
        for (gname, indices) in extra_groups
            for i in 1:(length(indices)-1)
                c1 = extra_coords[indices[i]]
                c2 = extra_coords[indices[i+1]]
                lines!(ax, [c1[1], c2[1]], [c1[2], c2[2]];
                       color=(:red, 0.6), linewidth=2)
            end
        end

        # Plot all extra points as circles
        ex_x = [c[1] for c in extra_coords]
        ex_y = [c[2] for c in extra_coords]
        scatter!(ax, ex_x, ex_y;
                 markersize=extra_point_size, color=:red, marker=:circle)
    end

    # Auto-zoom with margin
    if !isempty(x_vals)
        x_min, x_max = extrema(x_vals)
        y_min, y_max = extrema(y_vals)
        margin_x = 0.15 * (x_max - x_min) + 0.3
        margin_y = 0.15 * (y_max - y_min) + 0.3
        limits!(ax, x_min - margin_x, x_max + margin_x,
                    y_min - margin_y, y_max + margin_y)
    end

    # Legend
    legend_elements = [
        MarkerElement(color=:blue, marker=:circle, markersize=10),
    ]
    legend_labels = ["sim"]
    if !isnothing(extra_points)
        push!(legend_elements,
              MarkerElement(color=:red, marker=:circle, markersize=10))
        push!(legend_labels, "photogrammetry")
    end
    Legend(fig[1, 2], legend_elements, legend_labels)

    return fig
end

"""
    wrap_to_pi(angle)

Wrap angle to [-π, π] range.
"""
function wrap_to_pi(angle)
    return mod(angle + π, 2π) - π
end

"""
    euler_to_quaternion(roll_deg, pitch_deg, yaw_deg)

Convert Euler angles (in degrees) to quaternion.
Converts from NED to ENU frame:
  X_ENU = Y_NED (East)
  Y_ENU = X_NED (North)
  Z_ENU = -Z_NED (Up = -Down)
"""
function euler_to_quaternion(roll_deg, pitch_deg, yaw_deg)
    roll_rad = deg2rad(roll_deg)
    pitch_rad = deg2rad(pitch_deg)
    yaw_rad = deg2rad(yaw_deg)
    rot_ned = RotZYX(yaw_rad, pitch_rad, roll_rad)
    R_ned_to_enu = [0.0 1.0 0.0;
                    1.0 0.0 0.0;
                    0.0 0.0 -1.0]
    rot_enu = R_ned_to_enu * Matrix(rot_ned)
    q = SymbolicAWEModels.rotation_matrix_to_quaternion(rot_enu)
    return q
end

"""
    calc_heading(sys_struct::SystemStructure, R_b_w)

Calculate heading angle from rotation matrix, wrapped to [-π, π].
"""
function calc_heading(sys_struct::SymbolicAWEModels.SystemStructure, R_b_w)
    e_x = R_b_w[:, 1]
    wind_norm = [1,0,0]
    minus_e_x = -e_x
    proj_on_wind = dot(minus_e_x, wind_norm) * wind_norm
    e_x_perp = minus_e_x - proj_on_wind
    wind_cross_z = [wind_norm[2], -wind_norm[1], 0]
    heading_x = dot(e_x_perp, wind_cross_z)
    heading_z = e_x_perp[3]
    heading = atan(heading_x, heading_z)
    return wrap_to_pi(heading)
end

"""
    calc_csv_heading(roll_deg, pitch_deg, yaw_deg, sys_struct)

Calculate heading from CSV Euler angles, wrapped to [-π, π].
"""
function calc_csv_heading(roll_deg, pitch_deg, yaw_deg, sys_struct)
    quat = euler_to_quaternion(roll_deg, pitch_deg, yaw_deg)
    R = SymbolicAWEModels.quaternion_to_rotation_matrix(quat)
    heading = calc_heading(sys_struct, R)
    return wrap_to_pi(heading + π)
end
