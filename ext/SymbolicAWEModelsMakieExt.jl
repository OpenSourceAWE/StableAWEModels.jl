# Copyright (c) 2025 Bart van de Lint
# SPDX-License-Identifier: MPL-2.0

module SymbolicAWEModelsMakieExt

using Makie
using UnPack
using LinearAlgebra
using StaticArrays
using Statistics
using Printf
using KiteUtils
using KiteUtils: SysLog
using SymbolicAWEModels
using VortexStepMethod

function Makie.plot!(ax, sys::SystemStructure;
                     point_color = :darkred, segment_color = :black,
                     wing_colors = Makie.wong_colors(), vector_scale = 0.2,
                     show_points = true, show_segments = true, show_orient = true,
                     # Optional observables for real-time updates
                     segment_points_obs = nothing,
                     point_positions_obs = nothing,
                     wing_origins_obs = nothing,
                     wing_directions_obs = nothing)

    plots = Dict{Symbol, Any}()

    # === Plot Segments ===
    if show_segments
        if isnothing(segment_points_obs)
            # Static plotting: build segment points once
            lineseg_points = Point3f[]
            for seg in sys.segments
                p1 = sys.points[seg.point_idxs[1]].pos_w
                p2 = sys.points[seg.point_idxs[2]].pos_w
                push!(lineseg_points, Point3f(p1))
                push!(lineseg_points, Point3f(p2))
            end
        else
            # Dynamic plotting: use provided observable
            lineseg_points = segment_points_obs
        end

        num_segments = length(sys.segments)
        seg_colors = Observable(fill(to_color(segment_color), num_segments))

        plots[:segments] = linesegments!(ax, lineseg_points, color=seg_colors,
                                         linewidth=2, label="Segments", transparency=true)
    end

    # === Plot Points ===
    if show_points
        if isnothing(point_positions_obs)
            # Static plotting
            point_positions = [Point3f(p.pos_w) for p in sys.points]
        else
            # Dynamic plotting: use provided observable
            point_positions = point_positions_obs
        end
        plots[:points] = scatter!(ax, point_positions, color=point_color, label="Points",
                                  transparency=true)
    end

    # === Plot Wings ===
    if show_orient
        if isnothing(wing_origins_obs) || isnothing(wing_directions_obs)
            # Static plotting: create separate arrows for each wing
            plots[:wings] = []
            plots[:vsm] = []
            for (i, wing) in enumerate(sys.wings)
                wing_pos = Point3f(wing.pos_w)
                R = wing.R_b_w
                scale = vector_scale
                origins = [wing_pos, wing_pos, wing_pos]
                directions = [Vec3f(R[:, 1]) * scale, Vec3f(R[:, 2]) * scale, Vec3f(R[:, 3]) * scale]

                axis_colors = [:red, :green, :blue]
                p = arrows3d!(ax, origins, directions, color=axis_colors, label="Wing $i Axes")
                push!(plots[:wings], p)
                p = plot!(ax, wing.vsm_aero; R_b_w=wing.R_b_w, T_b_w=wing.pos_w)
                push!(plots[:vsm], p)
            end
        else
            # Dynamic plotting: single arrows plot with observables
            axis_colors = repeat([:red, :green, :blue], length(sys.wings))
            plots[:wings] = arrows3d!(ax, wing_origins_obs, wing_directions_obs,
                                     color=axis_colors)
        end
    end

    # === Plot Global Axes ===
    begin
        scale = vector_scale * 10 # Make global axes slightly larger
        origins = [Point3f(0, 0, 0), Point3f(0, 0, 0), Point3f(0, 0, 0)]
        directions = [Vec3f(10, 0, 0), Vec3f(0, 10, 0), Vec3f(0, 0, 10)]
        axis_colors = [:red, :green, :blue]
        plots[:global_axes] = arrows3d!(ax, origins, directions;
                                        shaftradius=0.02,
                                        tipradius=0.06,
                                        tiplength=0.1,
                                        color=axis_colors,
                                        label="Global Axes")
    end

    return plots
end

"""
    SymbolicAWEModels.update_plot_observables!(segment_points_obs, point_positions_obs,
                                               wing_origins_obs, wing_directions_obs,
                                               sys::SystemStructure; vector_scale=0.2)

Update Observable objects for real-time 3D visualization from a SystemStructure.

This function extracts the current state from `sys` (point positions, segment endpoints,
wing orientations) and updates the provided Observable objects. The observables will
trigger automatic updates in any Makie plots that use them.

# Arguments
- `segment_points_obs::Observable`: Observable containing segment line endpoints
- `point_positions_obs::Observable`: Observable containing point positions
- `wing_origins_obs::Observable`: Observable containing wing arrow origins
- `wing_directions_obs::Observable`: Observable containing wing arrow directions
- `sys::SystemStructure`: The system structure to extract data from

# Keyword Arguments
- `vector_scale::Real=0.2`: Scale factor for wing orientation arrows

# Example
```julia
# Create observables
seg_obs = Observable(Point3f[])
pts_obs = Observable(Point3f[])
orig_obs = Observable(Point3f[])
dir_obs = Observable(Vec3f[])

# Create plot with observables
scene = plot(sys_struct; segment_points_obs=seg_obs, point_positions_obs=pts_obs,
             wing_origins_obs=orig_obs, wing_directions_obs=dir_obs)

# In simulation loop:
for step in 1:steps
    next_step!(sam; ...)
    update_plot_observables!(seg_obs, pts_obs, orig_obs, dir_obs, sam.sys_struct)
    sleep(0.001)  # Allow Makie to process updates
end
```
"""
function SymbolicAWEModels.update_plot_observables!(segment_points_obs, point_positions_obs,
                                                     wing_origins_obs, wing_directions_obs,
                                                     sys::SystemStructure; vector_scale=0.2)
    # Update point positions
    if !isnothing(point_positions_obs)
        point_positions_obs[] = [Point3f(p.pos_w) for p in sys.points]
    end

    # Update segment endpoints
    if !isnothing(segment_points_obs)
        seg_points = Point3f[]
        for seg in sys.segments
            p1 = sys.points[seg.point_idxs[1]].pos_w
            p2 = sys.points[seg.point_idxs[2]].pos_w
            push!(seg_points, Point3f(p1))
            push!(seg_points, Point3f(p2))
        end
        segment_points_obs[] = seg_points
    end

    # Update wing orientations
    if !isnothing(wing_origins_obs) && !isnothing(wing_directions_obs)
        origins = Point3f[]
        directions = Vec3f[]
        for wing in sys.wings
            wing_pos = Point3f(wing.pos_w)
            R = wing.R_b_w
            # Add three arrow vectors for each axis (x, y, z in body frame)
            for i in 1:3
                push!(origins, wing_pos)
                push!(directions, Vec3f(R[:, i]) * vector_scale)
            end
        end
        wing_origins_obs[] = origins
        wing_directions_obs[] = directions
    end

    return nothing
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

"""
    Makie.plot(sys::SystemStructure, lg::SysLog; kwargs...)

Create a multi-panel plot of key simulation results from a `SysLog`.

This function visualizes various aspects of the kite's performance and state,
such as turn rates, reel-out speeds, aerodynamic forces, and wing deformation.
Each panel can be individually enabled or disabled via keyword arguments.

# Arguments
- `sys::SystemStructure`: The system structure, used to get component counts (e.g., number of groups).
- `lg::SysLog`: The simulation log data to be plotted.

# Keyword Arguments
- `plot_default::Bool=true`: Defaults to true, enabling all plot panels. If false, all panels are disabled.
- `plot_turn_rates::Bool=false`: Show the panel with the wing's angular velocities (ω_x, ω_y, ω_z).
- `plot_reelout::Bool=plot_default`: Show the panel with the reel-out velocities of the steering winches.
- `plot_aero_force::Bool=plot_default`: Show the panel with the z-component of aerodynamic force.
- `plot_aero_moment::Bool=false`: Show the panel with the y-component of aerodynamic moment.
- `plot_tether_moment::Bool=false`: Show the panel with the y-component of tether-induced moment.
- `plot_twist::Bool=plot_default`: Show the panel with the twist angles for each wing group.
- `plot_aoa::Bool=plot_default`: Show the panel with the angle of attack.
- `plot_heading::Bool=plot_default`: Show the panel with the kite's heading angle.
- `plot_elevation::Bool=false`: Show the panel with the kite's elevation angle.
- `plot_azimuth::Bool=false`: Show the panel with the kite's azimuth angle.
- `plot_winch_force::Bool=plot_default`: Show the panel with the winch forces.
- `plot_set_values::Bool=false`: Show the panel with the set torque values.
- `suffix::String=" - " * sys.name`: Suffix to append to plot labels.
- `size::Tuple=(1200, 800)`: Figure size in pixels.

# Example
```julia
# Plot only the angle of attack and heading
plot(model.sys_struct, log, plot_reelout=false, plot_aero_force=false, plot_twist=false, plot_winch_force=false)
```
"""
function Makie.plot(sys::SystemStructure, lg::SysLog;
                    plot_default=true,
                    plot_reelout=plot_default,
                    plot_aero_force=plot_default,
                    plot_twist=plot_default,
                    plot_aoa=plot_default,
                    plot_heading=plot_default,
                    plot_aero_moment=false,
                    plot_turn_rates=false,
                    plot_elevation=false,
                    plot_azimuth=false,
                    plot_tether_moment=false,
                    plot_winch_force=plot_default,
                    plot_set_values=false,
                    suffix=" - " * sys.name,
                    size=(1200, 800))

    sl = lg.syslog

    # Build list of panels to plot
    panels = []

    if plot_turn_rates
        turn_rates_deg = rad2deg.(hcat(sl.turn_rates...))
        push!(panels, (
            data = [turn_rates_deg[1,:], turn_rates_deg[2,:], turn_rates_deg[3,:]],
            labels = ["ω_x" * suffix, "ω_y" * suffix, "ω_z" * suffix],
            ylabel = "turn rates [°/s]"
        ))
    end

    if plot_reelout
        v_reelout_2 = [sl.v_reelout[i][2] for i in eachindex(sl.v_reelout)]
        v_reelout_3 = [sl.v_reelout[i][3] for i in eachindex(sl.v_reelout)]
        push!(panels, (
            data = [v_reelout_2, v_reelout_3],
            labels = ["v_ro[2]" * suffix, "v_ro[3]" * suffix],
            ylabel = "v_ro [m/s]"
        ))
    end

    if plot_aero_force
        aero_force_z = [sl.aero_force_b[i][3] for i in eachindex(sl.aero_force_b)]
        push!(panels, (
            data = [aero_force_z],
            labels = ["F_aero,z" * suffix],
            ylabel = "aero F [N]"
        ))
    end

    if plot_aero_moment
        moment_y = [sl.aero_moment_b[i][2] for i in eachindex(sl.aero_moment_b)]
        push!(panels, (
            data = [moment_y],
            labels = ["M_aero,y" * suffix],
            ylabel = "aero M [Nm]"
        ))
    end

    if plot_tether_moment
        moment_y = [sl.tether_induced_moment[i][2] for i in eachindex(sl.tether_moment)]
        push!(panels, (
            data = [moment_y],
            labels = ["M_tether,y" * suffix],
            ylabel = "tether M [Nm]"
        ))
    end

    if plot_twist && !isempty(sys.groups)
        twist_angles_deg = rad2deg.(hcat(sl.twist_angles...))[eachindex(sys.groups),:]
        twist_data = [twist_angles_deg[i,:] for i in eachindex(sys.groups)]
        twist_labels = ["twist[$i]" * suffix for i in eachindex(sys.groups)]
        push!(panels, (
            data = twist_data,
            labels = twist_labels,
            ylabel = "twist [°]"
        ))
    end

    if plot_aoa
        AoA_deg = rad2deg.(sl.AoA)
        push!(panels, (
            data = [AoA_deg],
            labels = ["AoA" * suffix],
            ylabel = "AoA [°]"
        ))
    end

    if plot_heading
        heading_deg = rad2deg.(sl.heading)
        push!(panels, (
            data = [heading_deg],
            labels = ["heading" * suffix],
            ylabel = "heading [°]"
        ))
    end

    if plot_elevation
        elevation_deg = rad2deg.(sl.elevation)
        push!(panels, (
            data = [elevation_deg],
            labels = ["elevation" * suffix],
            ylabel = "elevation [°]"
        ))
    end

    if plot_azimuth
        azimuth_deg = rad2deg.(sl.azimuth)
        push!(panels, (
            data = [azimuth_deg],
            labels = ["azimuth" * suffix],
            ylabel = "azimuth [°]"
        ))
    end

    if plot_winch_force
        winch_force = [[sl.winch_force[i][j] for i in eachindex(sl.winch_force)] for j in 1:3]
        push!(panels, (
            data = winch_force,
            labels = ["F_winch,1" * suffix, "F_winch,2" * suffix, "F_winch,3" * suffix],
            ylabel = "Winch force [N]"
        ))
    end

    if plot_set_values
        set_values = [[sl.set_torque[i][j] for i in eachindex(sl.set_torque)] for j in 1:3]
        push!(panels, (
            data = set_values,
            labels = ["T_winch,1" * suffix, "T_winch,2" * suffix, "T_winch,3" * suffix],
            ylabel = "Set torque [Nm]"
        ))
    end

    # Check if there's anything to plot
    if isempty(panels)
        error("No plot sections enabled. Enable at least one plot panel.")
    end

    # Create figure with subplots
    n_panels = length(panels)
    fig = Figure(size=size)

    axes = []
    for (i, panel) in enumerate(panels)
        # Share x-axis with first subplot
        if i == 1
            ax = Axis(fig[i, 1], ylabel=panel.ylabel)
        else
            ax = Axis(fig[i, 1], ylabel=panel.ylabel, xticklabelsvisible=false)
            linkxaxes!(axes[1], ax)
        end

        # Plot each data series in this panel
        for (j, (data_series, label)) in enumerate(zip(panel.data, panel.labels))
            lines!(ax, sl.time, data_series, label=label)
        end

        # Add legend if multiple traces
        if length(panel.data) > 1
            axislegend(ax, position=:rt)
        end

        push!(axes, ax)
    end

    # Add x-label to bottom subplot
    axes[end].xlabel = "time [s]"
    axes[end].xticklabelsvisible = true

    return fig
end

function zoom_out!(scene, cam, plots; relmargin=0.2)
    # --- ROBUST ZOOM OUT ---
    # 1. Get the current camera viewing direction vector
    inv_view_matrix = inv(cam.view[])
    cam_dir_vec = normalize(Vec3f(inv_view_matrix[1, 3],
                                  inv_view_matrix[2, 3],
                                  inv_view_matrix[3, 3]))
    # 2. Get the scene's bounding box and its center (the new target)
    bbox = data_limits(plots)
    center = bbox.origin .+ bbox.widths ./ 2
    # 3. Calculate the distance needed to see the whole box
    radius = norm(bbox.widths) / 2.0
    fov_rad = 2 * atan(1 / cam.projection[][2, 2])
    distance = radius / tan(fov_rad / 2.0)
    # 4. Calculate the new camera position
    new_eyepos = center + cam_dir_vec * (distance * (1+relmargin))
    # 5. Update the camera to the new "fit-all" view
    update_cam!(scene, new_eyepos, center)
end

function Makie.plot(sys::SystemStructure; 
                    size = (1200, 800), 
                    margin = 10.0,
                    relmargin = 0.2,
                    segment_color = :black, 
                    highlight_color = :red,
                    kwargs...)
    # Use LScene for advanced camera controls
    scene = Scene(; camera=cam3d!, show_axis = false, size, zoommode = :free, samples = 16)
    plots = plot!(scene, sys; segment_color, kwargs...)
    
    relevant_plots = AbstractPlot[]
    if haskey(plots, :segments)
        push!(relevant_plots, plots[:segments])
    end
    if haskey(plots, :points)
        push!(relevant_plots, plots[:points])
    end
    if haskey(plots, :wings)
        # plots[:wings] can be either an array (static) or a single plot (observables)
        if plots[:wings] isa AbstractArray
            append!(relevant_plots, plots[:wings])
        else
            push!(relevant_plots, plots[:wings])
        end
    end

    # --- Event Handling for Segments ---
    if haskey(plots, :segments)
        lineseg_plot = plots[:segments]
        seg_colors_obs = lineseg_plot.color
        last_hovered_idx = Ref(-1)
        zoomed_in = Ref(false)

        # --- Hover Labels ---
        # Segment index label at middle of segment
        segment_label = Observable("")
        segment_label_pos = Observable(Point2f(0, 0))
        segment_label_visible = Observable(false)
        text!(scene, segment_label, position = segment_label_pos, space = :pixel,
              fontsize = 14, color = :white, strokecolor = :black, strokewidth = 1,
              align = (:center, :center), visible = segment_label_visible, transparency = true)

        # Point index labels at segment endpoints
        point1_label = Observable("")
        point1_label_pos = Observable(Point2f(0, 0))
        point1_label_visible = Observable(false)
        text!(scene, point1_label, position = point1_label_pos, space = :pixel,
              fontsize = 14, color = :white, strokecolor = :black, strokewidth = 1,
              align = (:center, :center), visible = point1_label_visible, transparency = true)

        point2_label = Observable("")
        point2_label_pos = Observable(Point2f(0, 0))
        point2_label_visible = Observable(false)
        text!(scene, point2_label, position = point2_label_pos, space = :pixel,
              fontsize = 14, color = :white, strokecolor = :black, strokewidth = 1,
              align = (:center, :center), visible = point2_label_visible, transparency = true)

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
                    seg = sys.segments[hover_idx]
                    p1_3d = sys.points[seg.point_idxs[1]].pos_w
                    p2_3d = sys.points[seg.point_idxs[2]].pos_w

                    # Show segment index 20px to the right of middle of segment
                    mid_point_3d = (p1_3d + p2_3d) / 2
                    mid_point_2d = Makie.project(scene, mid_point_3d)
                    segment_label[] = string(hover_idx)
                    segment_label_pos[] = mid_point_2d + Point2f(20, 0)
                    segment_label_visible[] = true

                    # Show point indices 20px to the right of segment endpoints
                    p1_2d = Makie.project(scene, p1_3d)
                    p2_2d = Makie.project(scene, p2_3d)
                    point1 = sys.points[seg.point_idxs[1]].idx
                    point2 = sys.points[seg.point_idxs[2]].idx
                    point1_label[] = string(point1)
                    point1_label_pos[] = p1_2d + Point2f(20, 0)
                    point1_label_visible[] = true
                    point2_label[] = string(point2)
                    point2_label_pos[] = p2_2d + Point2f(20, 0)
                    point2_label_visible[] = true

                    new_colors[hover_idx] = to_color(highlight_color)
                else
                    segment_label_visible[] = false
                    point1_label_visible[] = false
                    point2_label_visible[] = false
                end
                seg_colors_obs[] = new_colors
                last_hovered_idx[] = hover_idx
            end
        end

        # --- Event Handler for Camera Movement ---
        on(scene.camera.view, priority = 1) do _
            # Update label positions when camera moves
            hover_idx = last_hovered_idx[]
            if hover_idx != -1
                seg = sys.segments[hover_idx]
                p1_3d = sys.points[seg.point_idxs[1]].pos_w
                p2_3d = sys.points[seg.point_idxs[2]].pos_w

                # Update segment label position
                mid_point_3d = (p1_3d + p2_3d) / 2
                mid_point_2d = Makie.project(scene, mid_point_3d)
                segment_label_pos[] = mid_point_2d + Point2f(20, 0)

                # Update point label positions
                p1_2d = Makie.project(scene, p1_3d)
                p2_2d = Makie.project(scene, p2_3d)
                point1_label_pos[] = p1_2d + Point2f(20, 0)
                point2_label_pos[] = p2_2d + Point2f(20, 0)
            end
        end

        # --- Update label positions after zoom ---
        on(scene.camera.view, priority = 0) do _
            # Small delay to ensure camera update is complete
            # This helps with label positioning after zoom operations
            sleep(0.01)
            hover_idx = last_hovered_idx[]
            if hover_idx != -1
                seg = sys.segments[hover_idx]
                p1_3d = sys.points[seg.point_idxs[1]].pos_w
                p2_3d = sys.points[seg.point_idxs[2]].pos_w

                # Update segment label position
                mid_point_3d = (p1_3d + p2_3d) / 2
                mid_point_2d = Makie.project(scene, mid_point_3d)
                segment_label_pos[] = mid_point_2d + Point2f(20, 0)

                # Update point label positions
                p1_2d = Makie.project(scene, p1_3d)
                p2_2d = Makie.project(scene, p2_3d)
                point1_label_pos[] = p1_2d + Point2f(20, 0)
                point2_label_pos[] = p2_2d + Point2f(20, 0)
            end
        end

        # --- Event Handler for Click-to-Zoom ---
        zoomed_in = Ref(false)
        on(events(scene).mousebutton, priority = 2) do event
            if event.button == Mouse.left && event.action == Mouse.press
                cam = scene.camera
                if !zoomed_in[] || last_hovered_idx[] != -1
                    # --- ZOOM IN --- (This part remains the same)
                    hover_idx = last_hovered_idx[]
                    if hover_idx != -1
                        seg = sys.segments[hover_idx]
                        p1_w = sys.points[seg.point_idxs[1]].pos_w
                        p2_w = sys.points[seg.point_idxs[2]].pos_w
                        
                        center = (p1_w + p2_w) / 2.0f0
                        segment_len = norm(p2_w - p1_w)
                        dist_heuristic = segment_len * 1.5 + 2.0
                        
                        inv_view_matrix = inv(cam.view[])
                        cam_dir_vec = normalize(Vec3f(inv_view_matrix[1, 3], inv_view_matrix[2, 3], inv_view_matrix[3, 3]))
                        new_eyepos = center + dist_heuristic * cam_dir_vec
                        
                        update_cam!(scene, new_eyepos, center)
                        zoomed_in[] = true
                        
                        # Update label positions after zoom
                        # Small delay to ensure camera update is complete
                        sleep(0.01)
                        p1_3d = sys.points[seg.point_idxs[1]].pos_w
                        p2_3d = sys.points[seg.point_idxs[2]].pos_w
                        
                        # Update segment label position
                        mid_point_3d = (p1_3d + p2_3d) / 2
                        mid_point_2d = Makie.project(scene, mid_point_3d)
                        segment_label_pos[] = mid_point_2d + Point2f(20, 0)
                        
                        # Update point label positions
                        p1_2d = Makie.project(scene, p1_3d)
                        p2_2d = Makie.project(scene, p2_3d)
                        point1_label_pos[] = p1_2d + Point2f(20, 0)
                        point2_label_pos[] = p2_2d + Point2f(20, 0)
                    end
                else
                    zoom_out!(scene, cam, relevant_plots; relmargin)
                    zoomed_in[] = false
                    
                    # Update label positions after zoom out
                    # Small delay to ensure camera update is complete
                    sleep(0.01)
                    hover_idx = last_hovered_idx[]
                    if hover_idx != -1
                        seg = sys.segments[hover_idx]
                        p1_3d = sys.points[seg.point_idxs[1]].pos_w
                        p2_3d = sys.points[seg.point_idxs[2]].pos_w
                        
                        # Update segment label position
                        mid_point_3d = (p1_3d + p2_3d) / 2
                        mid_point_2d = Makie.project(scene, mid_point_3d)
                        segment_label_pos[] = mid_point_2d + Point2f(20, 0)
                        
                        # Update point label positions
                        p1_2d = Makie.project(scene, p1_3d)
                        p2_2d = Makie.project(scene, p2_3d)
                        point1_label_pos[] = p1_2d + Point2f(20, 0)
                        point2_label_pos[] = p2_2d + Point2f(20, 0)
                    end
                end
                return Consume(true) # Consume the event
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
    pane_color = RGBAf(0.95, 0.95, 0.95, 0.3)
    # XZ plane at y_max (since camera is at negative y)
    mesh!(scene, Rect3(Vec3f(xlims[1], ylims[2], zlims[1]), Vec3f(xlims[2]-xlims[1], 0.01, zlims[2]-zlims[1])), color=pane_color)
    # YZ plane at x_max
    mesh!(scene, Rect3(Vec3f(xlims[2], ylims[1], zlims[1]), Vec3f(0.01, ylims[2]-ylims[1], zlims[2]-zlims[1])), color=pane_color)
    # XY plane at z_min
    mesh!(scene, Rect3(Vec3f(xlims[1], ylims[1], zlims[1]), Vec3f(xlims[2]-xlims[1], ylims[2]-ylims[1], 0.01)), color=pane_color)

    # Set initial camera position
    update_cam!(scene, Vec3f(-100, -100, 100), Vec3f(0, 0, 0))
    zoom_out!(scene, scene.camera, relevant_plots; relmargin)
    scene
end

"""
    Makie.plot(sam::SymbolicAWEModel, reltime::Real; kwargs...)

Plot a SymbolicAWEModel at a specific simulation time.

This is a convenience wrapper that extracts wing positions from the SAM
and calls the SystemStructure plotting function.

# Arguments
- `sam::SymbolicAWEModel`: The symbolic AWE model to plot
- `reltime::Real`: Relative time (not used for static plot, but kept for API compatibility)

# Keyword Arguments
All keyword arguments are passed through to `plot(::SystemStructure)`.
Common options include `size`, `margin`, `segment_color`, `highlight_color`, etc.
"""
function Makie.plot(sam::SymbolicAWEModel, reltime::Real=0.0; kwargs...)
    # Extract wing positions if available
    wings = sam.sys_struct.wings
    if length(wings) > 0
        wing_pos = [wing.pos_w for wing in wings]
    else
        wing_pos = nothing
    end

    # Call the SystemStructure plot function
    plot(sam.sys_struct; kwargs...)
end

"""
    replay(lg::SysLog, sys::SystemStructure; replay_speed=1.0, kwargs...)
    
Replay a SysLog with real-time 3D visualization.
    
This function replays a recorded simulation log with 3D visualization, updating the 
system structure at each time step according to the logged states. The replay speed 
can be adjusted to control the visualization speed.
    
# Arguments
- `lg::SysLog`: The simulation log to replay
- `sys::SystemStructure`: The system structure matching the log's topology
    
# Keyword Arguments
- `replay_speed::Real=1.0`: Replay speed factor (1.0 = real-time, 2.0 = 2x speed, etc.)
- All other keyword arguments are passed through to the SystemStructure plot function
    
# Example
```julia
# Replay a log at 2x speed
replay(log, sys_struct, replay_speed=2.0)

# Replay with custom visualization settings
replay(log, sys_struct, replay_speed=0.5, vector_scale=0.3)
```
"""
function SymbolicAWEModels.replay(lg::SysLog, sys::SystemStructure; 
                      replay_speed=1.0, 
                      segment_color = :black,
                      point_color = :darkred,
                      vector_scale = 0.2,
                      size = (1200, 800),
                      kwargs...)
    
    # Create observables for dynamic visualization
    segment_points_obs = Observable(Point3f[])
    point_positions_obs = Observable(Point3f[])
    wing_origins_obs = Observable(Point3f[])
    wing_directions_obs = Observable(Vec3f[])
    
    # Initialize with first state
    if !isempty(lg.syslog)
        update_from_sysstate!(sys, lg.syslog[1])
        update_plot_observables!(
            segment_points_obs, point_positions_obs,
            wing_origins_obs, wing_directions_obs,
            sys; vector_scale
        )
    end
    
    # Create 3D scene using existing plot function with observables
    scene = plot(sys;
                 segment_points_obs,
                 point_positions_obs,
                 wing_origins_obs,
                 wing_directions_obs,
                 segment_color,
                 point_color,
                 vector_scale,
                 size,
                 kwargs...)
    
    # Add text overlay for time display
    sim_time_text = Observable("Time: 0.00 s")
    text!(scene, sim_time_text, position = Point2f(20, 50), space = :pixel,
          fontsize = 24, color = :black, align = (:left, :top))
    
    # Display the scene
    display(scene)
    
    # Replay loop
    start_time = time()
    for (i, ss) in enumerate(lg.syslog)
        # Update system state from log
        update_from_sysstate!(sys, ss)
        
        # Update observables
        update_plot_observables!(
            segment_points_obs, point_positions_obs,
            wing_origins_obs, wing_directions_obs,
            sys; vector_scale
        )
        
        # Update time display
        sim_time_text[] = @sprintf("Time: %.2f s", ss.time)
        
        # Update label positions if labels were visible in the static plot
        # Note: In replay mode, hover labels are not currently supported
        # This is a placeholder for potential future implementation
        
        # Maintain replay speed timing
        if replay_speed > 0 && i < length(lg.syslog)
            next_time = lg.syslog[i+1].time
            current_time = ss.time
            dt = next_time - current_time
            
            target_elapsed = (current_time - lg.syslog[1].time) / replay_speed
            actual_elapsed = time() - start_time
            sleep_time = max(0.0, target_elapsed - actual_elapsed)
            sleep(sleep_time)
        end
        
        # Allow UI updates
        sleep(0.001)
    end
    
    return scene
end

end
