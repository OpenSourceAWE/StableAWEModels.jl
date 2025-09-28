# Copyright (c) 2025 Bart van de Lint
# SPDX-License-Identifier: MPL-2.0

module SymbolicAWEModelsMakieExt

using Makie
using UnPack
using LinearAlgebra
using StaticArrays
using Statistics
using SymbolicAWEModels

function Makie.plot!(ax, sys::SystemStructure;
                     point_color = :darkred, segment_color = :black,
                     wing_colors = Makie.wong_colors(), vector_scale = 3.0,
                     show_points = true, show_segments = true, show_orient = true)

    plots = Dict{Symbol, Any}()

    # === Plot Segments ===
    if show_segments
        lineseg_points = Point3f[]
        for seg in sys.segments
            p1 = sys.points[seg.point_idxs[1]].pos_w
            p2 = sys.points[seg.point_idxs[2]].pos_w
            push!(lineseg_points, Point3f(p1))
            push!(lineseg_points, Point3f(p2))
        end
        
        num_segments = length(sys.segments)
        seg_colors = Observable(fill(to_color(segment_color), num_segments))

        plots[:segments] = linesegments!(ax, lineseg_points, color=seg_colors,
                                         linewidth=3, label="Segments")
    end

    # === Plot Points ===
    if show_points
        point_positions = [Point3f(p.pos_w) for p in sys.points]
        plots[:points] = scatter!(ax, point_positions, color=point_color, label="Points")
    end

    # === Plot Wings ===
    if show_orient
        plots[:wings] = []
        for (i, wing) in enumerate(sys.wings)
            wing_pos = Point3f(wing.pos_w)
            color = wing_colors[mod1(i, length(wing_colors))]
            p = scatter!(ax, wing_pos, color=color, markersize=4, strokewidth=1, strokecolor=:black, label="Wing $i")
            push!(plots[:wings], p)

            R = wing.R_b_w
            scale = vector_scale
            origins = [wing_pos, wing_pos, wing_pos]
            directions = [Vec3f(R[:, 1]) * scale, Vec3f(R[:, 2]) * scale, Vec3f(R[:, 3]) * scale]

            axis_colors = [:red, :green, :blue]
            p = arrows3d!(ax, origins, directions, color=axis_colors, label="Wing $i Axes")
            push!(plots[:wings], p)
        end
    end

    # === Plot Global Axes ===
    begin
        scale = vector_scale * 1.5 # Make global axes slightly larger
        origins = [Point3f(0, 0, 0), Point3f(0, 0, 0), Point3f(0, 0, 0)]
        directions = [Vec3f(1, 0, 0) * scale, Vec3f(0, 1, 0) * scale, Vec3f(0, 0, 1) * scale]
        axis_colors = [:red, :green, :blue]
        plots[:global_axes] = arrows3d!(ax, origins, directions, color=axis_colors, label="Global Axes")
    end
    
    return plots
end

"""
    point_line_segment_distance(p, a, b)

Calculate the minimum 2D distance from a point `p` to a line segment defined
by points `a` and `b`.
"""
function point_line_segment_distance(p, a, b)
    # Vector from a to b
    ab = b - a
    # Vector from a to p
    ap = p - a
    # Length squared of the segment.
    len_sq = dot(ab, ab)
    # If the segment is just a point, return distance from p to a.
    if len_sq ≈ 0.0
        return norm(ap)
    end
    # Project ap onto ab to find the closest point on the infinite line.
    # t is the normalized position of the projection.
    t = dot(ap, ab) / len_sq
    # Clamp t to the [0, 1] range to stay on the segment.
    t_clamped = clamp(t, 0.0, 1.0)
    # Calculate the closest point on the segment.
    closest_point = a + t_clamped * ab
    # Return the distance from p to that closest point.
    return norm(p - closest_point)
end

function Makie.plot(sys::SystemStructure; 
                    size = (1200, 800), 
                    margin = 10.0, 
                    segment_color = :black, 
                    highlight_color = :red,
                    kwargs...)
    # Use LScene for advanced camera controls
    scene = Scene(; camera=cam3d!, show_axis=false, size, zoommode = :free)
    plots = plot!(scene, sys; segment_color, kwargs...)

    # --- Event Handling for Segments ---
    if haskey(plots, :segments)
        lineseg_plot = plots[:segments]
        seg_colors_obs = lineseg_plot.color
        last_hovered_idx = Ref(-1)
        original_cam_state = Ref{Any}(nothing)

        # --- Event Handler for Robust Hover Highlighting ---
        on(events(scene).mouseposition, priority = 2) do mp
            # This approach is more robust than `pick` for thin lines.
            # It finds the closest segment in 2D screen space.
            min_dist = Inf
            closest_seg_idx = -1
            mouse_pos_2d = Point2f(mp)
            margin_px = 30.0 # pixel margin for hover detection
            if haskey(plots, :segments)
                seg_points_3d = plots[:segments][1][]
                for i in 1:length(sys.segments)
                    p1_3d = seg_points_3d[2 * i - 1]
                    p2_3d = seg_points_3d[2 * i]
                    p1_2d = Makie.project(scene, p1_3d)
                    p2_2d = Makie.project(scene, p2_3d)

                    dist = point_line_segment_distance(mouse_pos_2d, p1_2d, p2_2d)

                    if dist < min_dist
                        min_dist = dist
                        closest_seg_idx = i
                    end
                end
            end

            hover_idx = (min_dist < margin_px) ? closest_seg_idx : -1
            if hover_idx != last_hovered_idx[]
                num_segments = length(sys.segments)
                new_colors = fill(to_color(segment_color), num_segments)
                if hover_idx != -1
                    point1 = sys.points[sys.segments[hover_idx].point_idxs[1]].idx
                    point2 = sys.points[sys.segments[hover_idx].point_idxs[2]].idx
                    println("Highlighting segment: $hover_idx," *
                            " connecting point $point1 to point $point2.")
                    new_colors[hover_idx] = to_color(highlight_color)
                end
                seg_colors_obs[] = new_colors
                last_hovered_idx[] = hover_idx
            end
        end

        # --- Event Handler for Click-to-Zoom ---
        on(events(scene).mousebutton, priority = 2) do event
            if event.button == Mouse.left && event.action == Mouse.press
                if original_cam_state[] === nothing
                    # --- ZOOM IN ---
                    original_cam_state[] = (copy(scene.camera.eyeposition[]),
                                            copy(scene.camera.view_direction[]))
                    hover_idx = last_hovered_idx[]
                    
                    if hover_idx != -1
                        # Zoom to highlighted segment
                        seg = sys.segments[hover_idx]
                        p1_w = sys.points[seg.point_idxs[1]].pos_w
                        p2_w = sys.points[seg.point_idxs[2]].pos_w
                        center = (p1_w + p2_w) / 2.0f0
                        
                        segment_len = norm(p2_w - p1_w)
                        dist_heuristic = segment_len * 1.5 + 2.0
                        
                        cam = scene.camera
                        cam_dir_vec = normalize(cam.eyeposition[] - cam.view_direction[])
                        new_eyepos = center + dist_heuristic * cam_dir_vec
                        update_cam!(scene, new_eyepos, center)
                    else
                        # Zoom to kite (all wings)
                        if !isempty(sys.wings)
                            cam = scene.camera
                            cam_dir_vec = normalize(cam.eyeposition[] - cam.view_direction[])
                            wing = sys.wings[1]
                            len = norm(wing.vsm_aero.panels[1].LE_point_1 -
                                       wing.vsm_aero.panels[end].LE_point_2)
                            dist_heuristic = len * 1.5 + 2.0
                            new_eyepos = wing.pos_w + dist_heuristic * cam_dir_vec
                            update_cam!(scene, new_eyepos, wing.pos_w)
                        end
                    end
                else
                    # --- ZOOM OUT ---
                    eyepos, view_direction = original_cam_state[]
                    update_cam!(scene, eyepos, view_direction)
                    bounding_box = data_limits(scene)
                    update_cam!(scene, bounding_box)
                    original_cam_state[] = nothing
                end
                return Consume(true) # Consume the event to prevent other interactions
            end
            return Consume(false)
        end
    end

    # --- Calculate limits and draw background panes ---
    xlims, ylims, zlims = (-10, 10), (-10, 10), (-10, 10) # Default limits
    if !isempty(sys.points)
        all_x = [p.pos_w[1] for p in sys.points]
        all_y = [p.pos_w[2] for p in sys.points]
        all_z = [p.pos_w[3] for p in sys.points]

        xlims_data = extrema(all_x)
        ylims_data = extrema(all_y)
        zlims_data = extrema(all_z)

        xlims = (xlims_data[1] - margin, xlims_data[2] + margin)
        ylims = (ylims_data[1] - margin, ylims_data[2] + margin)
        zlims = (zlims_data[1] - margin, zlims_data[2] + margin)
    end

    # Manually create background panes
    pane_color = RGBAf(0.95, 0.95, 0.95, 0.8)
    # XZ plane at y_max (since camera is at negative y)
    mesh!(scene, Rect3(Vec3f(xlims[1], ylims[2], zlims[1]), Vec3f(xlims[2]-xlims[1], 0.01, zlims[2]-zlims[1])), color=pane_color)
    # YZ plane at x_max
    mesh!(scene, Rect3(Vec3f(xlims[2], ylims[1], zlims[1]), Vec3f(0.01, ylims[2]-ylims[1], zlims[2]-zlims[1])), color=pane_color)
    # XY plane at z_min
    mesh!(scene, Rect3(Vec3f(xlims[1], ylims[1], zlims[1]), Vec3f(xlims[2]-xlims[1], ylims[2]-ylims[1], 0.01)), color=pane_color)

    # Set initial camera position
    update_cam!(scene, Vec3f(-100, -100, 100), Vec3f(0, 0, 0))
    bounding_box = data_limits(scene)
    update_cam!(scene, bounding_box)
    scene
end

end
