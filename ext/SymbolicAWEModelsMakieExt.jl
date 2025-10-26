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

# Global storage for plot observables (for time-based plotting)
const PLOT_OBSERVABLES = Ref{Union{Nothing, NamedTuple}}(nothing)
const PLOT_SCENE = Ref{Union{Nothing, Scene}}(nothing)
const PLOT_TIME_TEXT = Ref{Union{Nothing, Observable}}(nothing)
const PLOT_BACKGROUND_PANES = Ref{Union{Nothing, Vector}}(nothing)
const PLOT_MARGIN = Ref{Float64}(10.0)

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

    Makie.resize_to_layout!(fig)
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

function _plot_with_panes(sys::SystemStructure;
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
    function calculate_limits(sys, margin)
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
        return xlims, ylims, zlims
    end

    xlims, ylims, zlims = calculate_limits(sys, margin)

    # Create background panes with observables for dynamic updates
    pane_color = RGBAf(0.95, 0.95, 0.95, 0.3)
    pane_extent = 10000.0f0  # Large value for "infinite" extent

    # XZ plane at y_max - extends far in -X and +Z directions
    xz_pane_obs = Observable(Rect3(Vec3f(xlims[1] - pane_extent, ylims[2], zlims[1]),
                                    Vec3f(xlims[2] - xlims[1] + pane_extent, 0.01, pane_extent)))
    xz_pane = mesh!(scene, xz_pane_obs, color=pane_color)

    # YZ plane at x_max - extends far in -Y and +Z directions
    yz_pane_obs = Observable(Rect3(Vec3f(xlims[2], ylims[1] - pane_extent, zlims[1]),
                                    Vec3f(0.01, ylims[2] - ylims[1] + pane_extent, pane_extent)))
    yz_pane = mesh!(scene, yz_pane_obs, color=pane_color)

    # XY plane at z_min - extends far in -X and -Y directions
    xy_pane_obs = Observable(Rect3(Vec3f(xlims[1] - pane_extent, ylims[1] - pane_extent, zlims[1]),
                                    Vec3f(xlims[2] - xlims[1] + pane_extent, ylims[2] - ylims[1] + pane_extent, 0.01)))
    xy_pane = mesh!(scene, xy_pane_obs, color=pane_color)

    # Store pane observables
    pane_observables = [xz_pane_obs, yz_pane_obs, xy_pane_obs]

    # Set initial camera position
    update_cam!(scene, Vec3f(-100, -100, 100), Vec3f(0, 0, 0))
    zoom_out!(scene, scene.camera, relevant_plots; relmargin)

    # Return scene along with pane_observables and margin
    # These will be used by time-based plotting
    return scene, pane_observables, margin
end

# Public API function - returns just the scene for backward compatibility
function Makie.plot(sys::SystemStructure; kwargs...)
    scene, _, _ = _plot_with_panes(sys; kwargs...)
    return scene
end

"""
    Makie.plot(sam::SymbolicAWEModel, reltime::Real; kwargs...)

Plot a SymbolicAWEModel at a specific simulation time.

This is a convenience wrapper that calls `plot(sam.sys_struct, reltime)`.

# Arguments
- `sam::SymbolicAWEModel`: The symbolic AWE model to plot
- `reltime::Real`: Simulation time. Use 0.0 to create a new plot, any other value to update existing plot.

# Keyword Arguments
All keyword arguments are passed through to `plot(::SystemStructure, ::Real)`.
Common options include `size`, `margin`, `segment_color`, `highlight_color`, `vector_scale`, etc.
"""
function Makie.plot(sam::SymbolicAWEModel, reltime::Real=0.0; kwargs...)
    # Delegate to the SystemStructure time-based plot function
    plot(sam.sys_struct, reltime; kwargs...)
end

"""
    Makie.plot(sys::SystemStructure, time::Real; kwargs...)

Plot a SystemStructure at a specific simulation time, with automatic observable management.

When `time == 0.0`, this function creates a new 3D scene with observables for dynamic updates.
The observables are stored globally and reused for subsequent calls.

When `time != 0.0`, this function updates the existing observables from the current
SystemStructure state without creating a new scene.

# Arguments
- `sys::SystemStructure`: The system structure to plot
- `time::Real`: Simulation time. Use 0.0 to create a new plot, any other value to update existing plot.

# Keyword Arguments
All keyword arguments are passed through to `plot(::SystemStructure)`.
Common options include:
- `size::Tuple=(1200, 800)`: Figure size in pixels
- `margin::Real=10.0`: Margin around plot limits
- `relmargin::Real=0.2`: Relative margin for zoom operations
- `segment_color=:black`: Color for tether segments
- `highlight_color=:red`: Color for highlighted segments
- `vector_scale::Real=0.2`: Scale factor for wing orientation arrows
- `point_color=:darkred`: Color for point markers
- `show_points::Bool=true`: Whether to show points
- `show_segments::Bool=true`: Whether to show segments
- `show_orient::Bool=true`: Whether to show wing orientations

# Returns
- When `time == 0.0`: Returns the new Scene object
- When `time != 0.0`: Returns nothing (updates existing scene)

# Example
```julia
# Create initial plot
scene = plot(sys_struct, 0.0)

# In simulation loop, update the plot
for i in 1:100
    next_step!(sam)
    plot(sys_struct, i/sample_freq)
    sleep(0.01)
end
```
"""
function Makie.plot(sys::SystemStructure, time::Real;
                    vector_scale=0.2,
                    kwargs...)
    # Helper function to create new plot
    function create_new_plot()
        # Create new plot with observables
        segment_points_obs = Observable(Point3f[])
        point_positions_obs = Observable(Point3f[])
        wing_origins_obs = Observable(Point3f[])
        wing_directions_obs = Observable(Vec3f[])

        # Initialize observables from current state
        update_plot_observables!(
            segment_points_obs, point_positions_obs,
            wing_origins_obs, wing_directions_obs,
            sys; vector_scale
        )

        # Store observables globally for reuse
        PLOT_OBSERVABLES[] = (
            segment_points_obs = segment_points_obs,
            point_positions_obs = point_positions_obs,
            wing_origins_obs = wing_origins_obs,
            wing_directions_obs = wing_directions_obs
        )

        # Create scene with observables using internal function
        scene, pane_observables, margin = _plot_with_panes(sys;
                    segment_points_obs,
                    point_positions_obs,
                    wing_origins_obs,
                    wing_directions_obs,
                    vector_scale,
                    kwargs...)

        # Add time display overlay
        time_text = Observable(@sprintf("Time: %.2f s", time))
        text!(scene, time_text, position = Point2f(20, 50), space = :pixel,
              fontsize = 24, color = :black, align = (:left, :top))

        # Store scene, time text, and pane observables globally
        PLOT_SCENE[] = scene
        PLOT_TIME_TEXT[] = time_text
        PLOT_BACKGROUND_PANES[] = pane_observables
        PLOT_MARGIN[] = margin

        # Display the scene
        display(scene)

        return scene
    end

    # Check if we need to create a new plot
    if time == 0.0
        # User explicitly requested new plot
        return create_new_plot()
    else
        # Try to update existing plot
        if isnothing(PLOT_OBSERVABLES[]) || isnothing(PLOT_SCENE[])
            # No plot exists, create new one
            @warn "No plot exists. Creating new plot (call with time=0.0 to avoid this warning)."
            return create_new_plot()
        else
            # Check if the scene still has an active display
            scene = PLOT_SCENE[]
            scene_has_display = false

            # Check if scene is in any current screen
            try
                # In Makie, we can check the events.window_open observable
                # If the scene has events and window_open exists, check its value
                if hasfield(typeof(scene), :events) &&
                   hasfield(typeof(scene.events), :window_open)
                    scene_has_display = scene.events.window_open[]
                else
                    # Fallback: assume display exists (we'll create new one if update fails)
                    scene_has_display = true
                end
            catch
                # If checking fails, assume no display
                scene_has_display = false
            end

            if !scene_has_display
                # Display was closed, create new one
                return create_new_plot()
            else
                # Update existing observables
                obs = PLOT_OBSERVABLES[]
                update_plot_observables!(
                    obs.segment_points_obs, obs.point_positions_obs,
                    obs.wing_origins_obs, obs.wing_directions_obs,
                    sys; vector_scale
                )

                # Update time display
                if !isnothing(PLOT_TIME_TEXT[])
                    PLOT_TIME_TEXT[][] = @sprintf("Time: %.2f s", time)
                end

                # Update background panes
                if !isnothing(PLOT_BACKGROUND_PANES[])
                    panes = PLOT_BACKGROUND_PANES[]
                    margin = PLOT_MARGIN[]

                    # Recalculate limits based on current point positions
                    xlims, ylims, zlims = (-10, 10), (-10, 10), (-10, 10)
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

                    # Update pane observables with infinite extent
                    pane_extent = 10000.0f0

                    # XZ plane at y_max - extends far in -X and +Z directions
                    panes[1][] = Rect3(Vec3f(xlims[1] - pane_extent, ylims[2], zlims[1]),
                                       Vec3f(xlims[2] - xlims[1] + pane_extent, 0.01, pane_extent))
                    # YZ plane at x_max - extends far in -Y and +Z directions
                    panes[2][] = Rect3(Vec3f(xlims[2], ylims[1] - pane_extent, zlims[1]),
                                       Vec3f(0.01, ylims[2] - ylims[1] + pane_extent, pane_extent))
                    # XY plane at z_min - extends far in -X and -Y directions
                    panes[3][] = Rect3(Vec3f(xlims[1] - pane_extent, ylims[1] - pane_extent, zlims[1]),
                                       Vec3f(xlims[2] - xlims[1] + pane_extent, ylims[2] - ylims[1] + pane_extent, 0.01))
                end

                return nothing
            end
        end
    end
end

"""
    replay(lg::SysLog, sys::SystemStructure; replay_speed=1.0, autoplay=false, loop=false, kwargs...)

Replay a SysLog with interactive 3D visualization and playback controls.

This function creates an interactive viewer for a recorded simulation log with 3D visualization.
The viewer includes a slider for scrubbing through time, play/pause button, and frame stepping controls.

# Arguments
- `lg::SysLog`: The simulation log to replay
- `sys::SystemStructure`: The system structure matching the log's topology

# Keyword Arguments
- `replay_speed::Real=1.0`: Replay speed factor (1.0 = real-time, 2.0 = 2x speed, etc.)
- `autoplay::Bool=false`: Start playing automatically when opened
- `loop::Bool=false`: Loop playback continuously
- `vector_scale::Real=0.2`: Scale factor for wing orientation arrows
- All other keyword arguments are passed through to the SystemStructure plot function

# Returns
- A Figure with interactive controls and 3D visualization

# Example
```julia
# Create interactive replay viewer
fig = replay(log, sys_struct)

# Auto-play at 2x speed with looping
fig = replay(log, sys_struct, replay_speed=2.0, autoplay=true, loop=true)

# Replay with custom visualization settings
fig = replay(log, sys_struct, replay_speed=0.5, vector_scale=0.3)
```
"""
function SymbolicAWEModels.replay(lg::SysLog, sys::SystemStructure;
                      replay_speed=1.0,
                      autoplay=false,
                      loop=false,
                      vector_scale=0.2,
                      kwargs...)

    n_frames = length(lg.syslog)
    n_frames == 0 && error("Empty SysLog provided for replay")

    # Initialize with first state
    update_from_sysstate!(sys, lg.syslog[1])

    # Create initial plot using plot(sys, 0.0) which sets up observables, scene, and time display
    scene = plot(sys, 0.0; vector_scale, kwargs...)

    # Create figure for the complete viewer
    fig = Figure(size=(1200, 900))

    # Add the scene to the figure
    fig[1, 1] = scene

    # Control panel layout
    control_grid = fig[2, 1] = GridLayout()

    # Create observable for current frame index
    frame_idx = Observable(1)

    # Function to update to a specific frame
    function update_frame!(idx)
        ss = lg.syslog[idx]
        update_from_sysstate!(sys, ss)
        plot(sys, ss.time; vector_scale)
    end

    # Create slider for frame selection
    sl = Slider(control_grid[1, 1:3], range=1:n_frames, startvalue=1)

    # Play/Pause button
    is_playing = Observable(autoplay)
    play_button = Button(control_grid[2, 1], label=@lift($is_playing ? "Pause" : "Play"))

    # Step forward/backward buttons
    step_back_button = Button(control_grid[2, 2], label="<")
    step_forward_button = Button(control_grid[2, 3], label=">")

    # Frame counter and time label
    frame_label = Label(control_grid[3, 1:3],
                       text=@lift("Frame: $($(frame_idx))/$n_frames | Time: $(@sprintf("%.2f", lg.syslog[$(frame_idx)].time)) s | Speed: $(replay_speed)x"),
                       halign=:center)

    # Connect slider to frame updates
    on(sl.value) do val
        frame_idx[] = val
        update_frame!(val)
    end

    # Play button functionality
    on(play_button.clicks) do _
        is_playing[] = !is_playing[]
    end

    # Step buttons
    on(step_back_button.clicks) do _
        sl.value[] = max(1, frame_idx[] - 1)
    end

    on(step_forward_button.clicks) do _
        sl.value[] = min(n_frames, frame_idx[] + 1)
    end

    # Animation loop with replay speed
    @async begin
        while true
            if is_playing[]
                if frame_idx[] < n_frames
                    sl.value[] = frame_idx[] + 1
                    # Calculate sleep time based on actual time difference and replay speed
                    if frame_idx[] > 1
                        dt = lg.syslog[frame_idx[]].time - lg.syslog[frame_idx[] - 1].time
                        sleep(max(0.01, dt / replay_speed))
                    else
                        sleep(0.05)
                    end
                elseif loop
                    sl.value[] = 1  # Loop back to start
                else
                    is_playing[] = false  # Stop at end
                end
            end
            sleep(0.02)  # Check state frequently
        end
    end

    return fig
end

end
