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
        seg_colors = Observable(fill(to_color(segment_color), 2 * num_segments))

        plots[:segments] = linesegments!(ax, lineseg_points, color=seg_colors, linewidth=3, label="Segments")
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

function Makie.plot(sys::SystemStructure; 
                    size = (1200, 800), 
                    margin = 10.0, 
                    segment_color = :black, 
                    highlight_color = :red,
                    kwargs...)
    fig = Figure(; size)
    
    # Use LScene for advanced camera controls
    ax = LScene(fig[1, 1]; show_axis=false)
    plots = plot!(ax, sys; segment_color, kwargs...)

    # --- Event Handling for Segments ---
    if haskey(plots, :segments)
        lineseg_plot = plots[:segments]
        seg_colors_obs = lineseg_plot.color
        
        last_hovered_idx = Ref(-1)

        on(events(ax.scene).mouseposition, priority = 2) do mp
            plot, idx = pick(ax.scene)
            
            if plot === lineseg_plot
                seg_idx = ceil(Int, idx / 2)
                if seg_idx != last_hovered_idx[]
                    if last_hovered_idx[] != -1
                        seg_colors_obs[][2 * last_hovered_idx[] - 1] = to_color(segment_color)
                        seg_colors_obs[][2 * last_hovered_idx[]] = to_color(segment_color)
                    end
                    
                    seg_colors_obs[][2 * seg_idx - 1] = to_color(highlight_color)
                    seg_colors_obs[][2 * seg_idx] = to_color(highlight_color)
                    
                    notify(seg_colors_obs)
                    last_hovered_idx[] = seg_idx
                end
            else
                if last_hovered_idx[] != -1
                    seg_colors_obs[][2 * last_hovered_idx[] - 1] = to_color(segment_color)
                    seg_colors_obs[][2 * last_hovered_idx[]] = to_color(segment_color)
                    notify(seg_colors_obs)
                    last_hovered_idx[] = -1
                end
            end
        end

        on(events(fig).mousebutton) do event
            if event.button == Mouse.left && event.action == Mouse.press
                plot, idx = pick(ax.scene)
                if plot === lineseg_plot
                    seg_idx = ceil(Int, idx / 2)
                    
                    p1 = lineseg_plot[1][][2 * seg_idx - 1]
                    p2 = lineseg_plot[1][][2 * seg_idx]
                    
                    center = (p1 + p2) / 2
                    seg_length = norm(p2 - p1)
                    
                    cam_controls = ax.scene.camera_controls
                    up_vec = cam_controls.upvector[]
                    view_dir = normalize(cam_controls.lookat[] - cam_controls.eyeposition[])
                    side_vec = cross(view_dir, up_vec)

                    new_eyepos = center - view_dir * (seg_length * 2.5) + up_vec * (seg_length * 0.5) + side_vec * (seg_length * 0.5)
                    
                    update_cam!(ax.scene, new_eyepos, center)
                end
            end
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
    mesh!(ax, Rect3(Vec3f(xlims[1], ylims[2], zlims[1]), Vec3f(xlims[2]-xlims[1], 0.01, zlims[2]-zlims[1])), color=pane_color)
    # YZ plane at x_max
    mesh!(ax, Rect3(Vec3f(xlims[2], ylims[1], zlims[1]), Vec3f(0.01, ylims[2]-ylims[1], zlims[2]-zlims[1])), color=pane_color)
    # XY plane at z_min
    mesh!(ax, Rect3(Vec3f(xlims[1], ylims[1], zlims[1]), Vec3f(xlims[2]-xlims[1], ylims[2]-ylims[1], 0.01)), color=pane_color)

    # Set initial camera position
    update_cam!(ax.scene, Vec3f(-10, -10, 10), Vec3f(0, 0, 0))

    return fig
end

end
