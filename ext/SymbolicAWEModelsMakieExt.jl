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
const PLOT_BACKGROUND_PANES = Ref{Union{Nothing, Vector}}(nothing)
const PLOT_MARGIN = Ref{Float64}(10.0)

"""
    calculate_segment_force_colors(segments, segment_color)

Calculate segment colors based on their force values.
Maps forces from green (low) to red (high) using linear interpolation.

# Arguments
- `segments`: Collection of segments with force field
- `segment_color`: Default color to use when all forces are equal

# Returns
- Vector of RGBAf colors, one per segment
"""
function calculate_segment_force_colors(segments, segment_color)
    forces = [seg.force for seg in segments]
    max_force = maximum(forces)
    min_force = minimum(forces)
    force_range = max_force - min_force

    return [begin
        if force_range > 0
            normalized_force = (seg.force - min_force) / force_range
            # Interpolate from green (0,1,0) to red (1,0,0)
            RGBAf(normalized_force, 1.0 - normalized_force, 0.0, 1.0)
        else
            to_color(segment_color)
        end
    end for seg in segments]
end

function Makie.plot!(ax, sys::SystemStructure;
                     point_color = :darkred, segment_color = :black,
                     wing_colors = Makie.wong_colors(), vector_scale = 0.2,
                     show_points = true, show_segments = true, show_orient = true,
                     show_panes = true, margin = 10.0, force_color = false,
                     # Optional observables for real-time updates
                     segment_points_obs = nothing,
                     point_positions_obs = nothing,
                     wing_origins_obs = nothing,
                     wing_directions_obs = nothing,
                     pane_observables = nothing)

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

        # Calculate segment colors based on force if requested
        if force_color
            seg_colors = Observable(calculate_segment_force_colors(sys.segments, segment_color))
        else
            seg_colors = Observable(fill(to_color(segment_color), num_segments))
        end

        plots[:segments] = linesegments!(ax, lineseg_points, color=seg_colors,
                                         linewidth=2, label="Segments", transparency=true)
        plots[:segment_colors_obs] = seg_colors
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

    # === Calculate system scale for axes and panes ===
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

    # Calculate characteristic length for scaling arrows
    char_length = max(xlims[2] - xlims[1], ylims[2] - ylims[1], zlims[2] - zlims[1])
    axis_length = char_length * 0.15

    # === Plot Global Axes ===
    begin
        origins = [Point3f(0, 0, 0), Point3f(0, 0, 0), Point3f(0, 0, 0)]
        directions = [Vec3f(axis_length, 0, 0), Vec3f(0, axis_length, 0), Vec3f(0, 0, axis_length)]
        axis_colors = [:red, :green, :blue]
        plots[:global_axes] = arrows3d!(ax, origins, directions;
                                        shaftradius=axis_length*0.002,
                                        tipradius=axis_length*0.006,
                                        tiplength=axis_length*0.01,
                                        color=axis_colors,
                                        label="Global Axes")
    end

    # === Plot Background Panes ===
    if show_panes
        pane_color = RGBAf(0.95, 0.95, 0.95, 0.3)
        pane_extent = 10000.0f0  # Large value for "infinite" extent

        if isnothing(pane_observables)
            # Static plotting: create new observables
            xz_pane_obs = Observable(Rect3(Vec3f(xlims[1] - pane_extent, ylims[2], zlims[1]),
                                            Vec3f(xlims[2] - xlims[1] + pane_extent, 0.01, pane_extent)))
            yz_pane_obs = Observable(Rect3(Vec3f(xlims[2], ylims[1] - pane_extent, zlims[1]),
                                            Vec3f(0.01, ylims[2] - ylims[1] + pane_extent, pane_extent)))
            xy_pane_obs = Observable(Rect3(Vec3f(xlims[1] - pane_extent, ylims[1] - pane_extent, zlims[1]),
                                            Vec3f(xlims[2] - xlims[1] + pane_extent, ylims[2] - ylims[1] + pane_extent, 0.01)))
            plots[:pane_observables] = [xz_pane_obs, yz_pane_obs, xy_pane_obs]
        else
            # Dynamic plotting: use provided observables and update them
            xz_pane_obs, yz_pane_obs, xy_pane_obs = pane_observables
            xz_pane_obs[] = Rect3(Vec3f(xlims[1] - pane_extent, ylims[2], zlims[1]),
                                  Vec3f(xlims[2] - xlims[1] + pane_extent, 0.01, pane_extent))
            yz_pane_obs[] = Rect3(Vec3f(xlims[2], ylims[1] - pane_extent, zlims[1]),
                                  Vec3f(0.01, ylims[2] - ylims[1] + pane_extent, pane_extent))
            xy_pane_obs[] = Rect3(Vec3f(xlims[1] - pane_extent, ylims[1] - pane_extent, zlims[1]),
                                  Vec3f(xlims[2] - xlims[1] + pane_extent, ylims[2] - ylims[1] + pane_extent, 0.01))
            plots[:pane_observables] = pane_observables
        end

        # Create mesh plots for the panes
        xz_pane = mesh!(ax, plots[:pane_observables][1], color=pane_color)
        yz_pane = mesh!(ax, plots[:pane_observables][2], color=pane_color)
        xy_pane = mesh!(ax, plots[:pane_observables][3], color=pane_color)
        plots[:panes] = [xz_pane, yz_pane, xy_pane]
    end

    return plots
end

"""
    SymbolicAWEModels.update_plot_observables!(segment_points_obs, point_positions_obs,
                                               wing_origins_obs, wing_directions_obs,
                                               sys::SystemStructure; vector_scale=0.2,
                                               segment_colors_obs=nothing, force_color=false,
                                               segment_color=:black)

Update Observable objects for real-time 3D visualization from a SystemStructure.

# Arguments
- `segment_points_obs::Observable`: Observable containing segment line endpoints
- `point_positions_obs::Observable`: Observable containing point positions
- `wing_origins_obs::Observable`: Observable containing wing arrow origins
- `wing_directions_obs::Observable`: Observable containing wing arrow directions
- `sys::SystemStructure`: The system structure to extract data from

# Keyword Arguments
- `vector_scale::Real=0.2`: Scale factor for wing orientation arrows
- `segment_colors_obs::Union{Observable,Nothing}=nothing`: Observable for segment colors
- `force_color::Bool=false`: Whether to color segments by force
- `segment_color=:black`: Default color for segments when force_color is false

# Example
```julia
# Create observables
seg_obs = Observable(Point3f[])
pts_obs = Observable(Point3f[])
orig_obs = Observable(Point3f[])
dir_obs = Observable(Vec3f[])
seg_colors_obs = Observable(RGBAf[])

# Create plot with observables
scene = plot(sys_struct; segment_points_obs=seg_obs, point_positions_obs=pts_obs,
             wing_origins_obs=orig_obs, wing_directions_obs=dir_obs, force_color=true)

# In simulation loop:
for step in 1:steps
    next_step!(sam; ...)
    update_plot_observables!(seg_obs, pts_obs, orig_obs, dir_obs, sam.sys_struct;
                            segment_colors_obs=seg_colors_obs, force_color=true)
    sleep(0.001)  # Allow Makie to process updates
end
```
"""
function SymbolicAWEModels.update_plot_observables!(segment_points_obs, point_positions_obs,
                                                     wing_origins_obs, wing_directions_obs,
                                                     sys::SystemStructure; vector_scale=0.2,
                                                     segment_colors_obs=nothing, force_color=false,
                                                     segment_color=:black)
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

    # Update segment colors if force_color is enabled
    if !isnothing(segment_colors_obs) && force_color
        segment_colors_obs[] = calculate_segment_force_colors(sys.segments, segment_color)
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
                    force_color = false,
                    kwargs...)
    # Use LScene for advanced camera controls
    scene = Scene(; camera=cam3d!, show_axis = false, size, zoommode = :free, samples = 16)
    plots = plot!(scene, sys; segment_color, margin, force_color, kwargs...)
    
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

                # Calculate base colors (either force-based or uniform)
                if force_color
                    new_colors = calculate_segment_force_colors(sys.segments, segment_color)
                else
                    new_colors = fill(to_color(segment_color), num_segments)
                end

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

        # --- Event Handler for Click-to-Zoom ---
        zoomed_in = Ref(false)
        on(events(scene).mousebutton, priority = 1) do event
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
                        
                        return Consume(true) # Consume the event
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
                    
                    return Consume(true) # Consume the event
                end
            end
            return Consume(false)
        end
    end

    # Extract pane observables from plots (created by plot!())
    pane_observables = haskey(plots, :pane_observables) ? plots[:pane_observables] : nothing

    # Set initial camera position
    update_cam!(scene, Vec3f(-100, -100, 100), Vec3f(0, 0, 0))
    zoom_out!(scene, scene.camera, relevant_plots; relmargin)

    # Return scene along with pane_observables, margin, and plots dict
    # These will be used by time-based plotting
    return scene, pane_observables, margin, plots
end

# Public API function - creates scene with observables for dynamic updates
function Makie.plot(sys::SystemStructure;
                    vector_scale=0.2,
                    force_color=false,
                    segment_color=:black,
                    kwargs...)
    # Create new plot with observables for dynamic updates
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

    # Create scene with observables using internal function
    scene, pane_observables, margin, plots = _plot_with_panes(sys;
                segment_points_obs,
                point_positions_obs,
                wing_origins_obs,
                wing_directions_obs,
                vector_scale,
                force_color,
                segment_color,
                kwargs...)

    # Get the segment colors observable from the plots if available
    segment_colors_obs = nothing
    if haskey(plots, :segment_colors_obs)
        segment_colors_obs = plots[:segment_colors_obs]
    end

    # Store observables globally for reuse by plot!()
    PLOT_OBSERVABLES[] = (
        segment_points_obs = segment_points_obs,
        point_positions_obs = point_positions_obs,
        wing_origins_obs = wing_origins_obs,
        wing_directions_obs = wing_directions_obs,
        segment_colors_obs = segment_colors_obs,
        force_color = force_color,
        segment_color = segment_color
    )

    # Store scene and pane observables globally
    PLOT_SCENE[] = scene
    PLOT_BACKGROUND_PANES[] = pane_observables
    PLOT_MARGIN[] = margin

    return scene
end

"""
    Makie.plot!(sys::SystemStructure; vector_scale=0.2)

Update the currently displayed SystemStructure plot with new data from `sys`.

This function follows standard Makie conventions: `plot!` with `!` mutates the existing
scene by updating its observables. Must be called after an initial `plot(sys)` has created
the scene and observables.

# Arguments
- `sys::SystemStructure`: The system structure with updated state to display

# Keyword Arguments
- `vector_scale::Real=0.2`: Scale factor for wing orientation arrows

# Returns
- `nothing` (mutates existing scene via observables)

# Example
```julia
# Create initial plot
scene = plot(sys_struct)

# In simulation loop, update the plot
for i in 1:100
    next_step!(sam)
    plot!(sys_struct)  # Updates observables
    sleep(0.01)
end
```

# See Also
- `plot(::SystemStructure)`: Create a new scene (non-mutating)
- `update_plot_observables!`: Lower-level observable update function
"""
function Makie.plot!(sys::SystemStructure; vector_scale=0.2)
    # Check if observables exist
    if isnothing(PLOT_OBSERVABLES[]) || isnothing(PLOT_SCENE[])
        error("No existing plot to update. Call plot(sys) first to create a scene.")
    end

    # Check if the scene is still active
    scene = PLOT_SCENE[]
    scene_has_display = false
    try
        if hasfield(typeof(scene), :events) &&
           hasfield(typeof(scene.events), :window_open)
            scene_has_display = scene.events.window_open[]
        else
            scene_has_display = true
        end
    catch
        scene_has_display = false
    end

    if !scene_has_display
        error("Plot window has been closed. Call plot(sys) to create a new scene.")
    end

    # Update existing observables
    obs = PLOT_OBSERVABLES[]
    update_plot_observables!(
        obs.segment_points_obs, obs.point_positions_obs,
        obs.wing_origins_obs, obs.wing_directions_obs,
        sys; vector_scale,
        segment_colors_obs=obs.segment_colors_obs,
        force_color=obs.force_color,
        segment_color=obs.segment_color
    )

    # Update background panes based on current point positions
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

"""
    record(lg::SysLog, sys::SystemStructure, filename::String; framerate=30, kwargs...)

Record a SysLog animation to an MP4 video file.

# Arguments
- `lg::SysLog`: The simulation log to record
- `sys::SystemStructure`: The system structure matching the log's topology
- `filename::String`: Output video filename (e.g., "simulation.mp4")

# Keyword Arguments
- `framerate::Int=30`: Video framerate (frames per second)
- `vector_scale::Real=0.2`: Scale factor for wing orientation arrows
- All other keyword arguments are passed through to the SystemStructure plot function

# Returns
- The Scene object used for recording

# Example
```julia
# Record simulation to MP4
record(log, sys_struct, "output.mp4")

# Record with custom framerate and settings
record(log, sys_struct, "simulation.mp4", framerate=60, vector_scale=0.3)
```
"""
function SymbolicAWEModels.record(lg::SysLog, sys::SystemStructure, filename::String;
                                   framerate::Int=30,
                                   vector_scale::Real=0.2,
                                   kwargs...)
    n_frames = length(lg.syslog)
    n_frames == 0 && error("Empty SysLog provided for recording")

    println("Recording video to: $filename")
    println("Framerate: $framerate fps")
    println("Total frames: $n_frames")

    # Initialize with first state
    update_from_sysstate!(sys, lg.syslog[1])

    # Create initial plot with observables (following standard Makie pattern)
    scene = plot(sys; vector_scale, kwargs...)

    # Record video by stepping through all frames
    Makie.record(scene, filename, 1:n_frames; framerate=framerate) do frame_num
        # Update system state
        update_from_sysstate!(sys, lg.syslog[frame_num])

        # Update plot observables (mutating, following standard Makie pattern)
        plot!(sys; vector_scale)

        # Critical: Yield to GLMakie's render thread to process observable updates
        # Without this, the record function captures the same (first) frame repeatedly
        sleep(0.001)

        # Print progress every 10% or at the end
        if frame_num % max(1, div(n_frames, 10)) == 0 || frame_num == n_frames
            progress_pct = round(100 * frame_num / n_frames, digits=1)
            println("  Progress: $progress_pct% ($frame_num/$n_frames frames)")
        end
    end

    println("Video saved successfully!")
    return scene
end

"""
    replay(lg::SysLog, sys::SystemStructure; replay_speed=1.0, autoplay=false, loop=false, kwargs...)

Replay a SysLog with interactive 3D visualization and playback controls.

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
- A Scene with interactive controls overlaid on the 3D visualization

# UI Controls
- **Slider**: Drag to scrub through time
- **Play/Pause button**: Toggle playback (green when stopped, red when playing)
- **< button**: Step backward one frame
- **> button**: Step forward one frame
- **Info display**: Shows current frame, time, and playback speed

# Example
```julia
# Create interactive replay viewer
scene = replay(log, sys_struct)

# Auto-play at 2x speed with looping
scene = replay(log, sys_struct, replay_speed=2.0, autoplay=true, loop=true)

# Replay with custom visualization settings
scene = replay(log, sys_struct, replay_speed=0.5, vector_scale=0.3)
```

# See Also
- `record`: For saving replay to MP4 video file
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

    # Create initial plot which sets up observables and scene (following standard Makie pattern)
    scene = plot(sys; vector_scale, kwargs...)

    # Create pixel-space subscene for UI controls overlay
    ui_scene = Scene(scene, viewport=scene.viewport, clear=false, camera=campixel!)

    # Create observable for current frame index
    frame_idx = Observable(1)

    # Function to update to a specific frame (using plot! to mutate observables)
    function update_frame!(idx)
        ss = lg.syslog[idx]
        update_from_sysstate!(sys, ss)
        plot!(sys; vector_scale)
    end

    # UI Layout constants (pixel coordinates from bottom-left)
    ui_height = 120
    ui_margin = 20
    button_height = 30
    button_width = 80
    slider_height = 20

    # Get scene size
    scene_width = Observable(scene.viewport[].widths[1])
    scene_height = Observable(scene.viewport[].widths[2])

    # Update scene size on resize
    on(scene.viewport) do area
        scene_width[] = area.widths[1]
        scene_height[] = area.widths[2]
    end

    # --- Slider Implementation ---
    slider_y = ui_margin + button_height + 10
    slider_x_start = ui_margin
    slider_width_obs = @lift($(scene_width) - 2 * ui_margin)

    # Slider track
    slider_track_rect = @lift(Rect2f(slider_x_start, slider_y, $(slider_width_obs), slider_height))
    poly!(ui_scene, slider_track_rect, color=RGBAf(0.3, 0.3, 0.3, 0.8))

    # Slider thumb position
    slider_thumb_x = @lift(slider_x_start + ($(frame_idx) - 1) / (n_frames - 1) * $(slider_width_obs))
    slider_thumb_pos = @lift(Point2f($(slider_thumb_x), slider_y + slider_height / 2))
    scatter!(ui_scene, slider_thumb_pos, color=:white, markersize=15)

    # --- Button Implementation ---
    button_y = ui_margin

    # Play/Pause button
    is_playing = Observable(autoplay)
    play_button_x = ui_margin
    play_button_rect = Rect2f(play_button_x, button_y, button_width, button_height)
    play_button_color = @lift($(is_playing) ? RGBAf(0.8, 0.3, 0.3, 0.8) : RGBAf(0.3, 0.8, 0.3, 0.8))
    poly!(ui_scene, play_button_rect, color=play_button_color)
    play_button_label = @lift($(is_playing) ? "Pause" : "Play")
    text!(ui_scene, play_button_label, position=Point2f(play_button_x + button_width/2, button_y + button_height/2),
          align=(:center, :center), fontsize=14, color=:white)

    # Step backward button
    step_back_x = play_button_x + button_width + 10
    step_back_rect = Rect2f(step_back_x, button_y, button_width, button_height)
    poly!(ui_scene, step_back_rect, color=RGBAf(0.4, 0.4, 0.4, 0.8))
    text!(ui_scene, "<", position=Point2f(step_back_x + button_width/2, button_y + button_height/2),
          align=(:center, :center), fontsize=14, color=:white)

    # Step forward button
    step_forward_x = step_back_x + button_width + 10
    step_forward_rect = Rect2f(step_forward_x, button_y, button_width, button_height)
    poly!(ui_scene, step_forward_rect, color=RGBAf(0.4, 0.4, 0.4, 0.8))
    text!(ui_scene, ">", position=Point2f(step_forward_x + button_width/2, button_y + button_height/2),
          align=(:center, :center), fontsize=14, color=:white)

    # Info label
    info_label_x = step_forward_x + button_width + 20
    info_text = @lift("Frame: $($(frame_idx))/$n_frames | Time: $(@sprintf("%.2f", lg.syslog[$(frame_idx)].time)) s | Speed: $(replay_speed)x")
    text!(ui_scene, info_text, position=Point2f(info_label_x, button_y + button_height/2),
          align=(:left, :center), fontsize=14, color=:white, strokecolor=:black, strokewidth=1)

    # Combined mouse event handling for slider and buttons
    slider_dragging = Ref(false)

    on(events(ui_scene).mousebutton, priority = 2) do event
        if event.button == Mouse.left
            mp = events(ui_scene).mouseposition[]

            if event.action == Mouse.press
                # Check if click is on slider
                track_rect = slider_track_rect[]
                if mp[2] >= track_rect.origin[2] && mp[2] <= track_rect.origin[2] + track_rect.widths[2] &&
                   mp[1] >= track_rect.origin[1] && mp[1] <= track_rect.origin[1] + track_rect.widths[1]
                    slider_dragging[] = true
                    # Update frame immediately
                    rel_pos = clamp((mp[1] - slider_x_start) / slider_width_obs[], 0.0, 1.0)
                    new_idx = round(Int, 1 + rel_pos * (n_frames - 1))
                    frame_idx[] = clamp(new_idx, 1, n_frames)
                    update_frame!(frame_idx[])
                    return Consume(true)
                end
                
                # Check play button
                if mp[1] >= play_button_rect.origin[1] && mp[1] <= play_button_rect.origin[1] + play_button_rect.widths[1] &&
                   mp[2] >= play_button_rect.origin[2] && mp[2] <= play_button_rect.origin[2] + play_button_rect.widths[2]
                    is_playing[] = !is_playing[]
                    return Consume(true)
                end

                # Check step back button
                if mp[1] >= step_back_rect.origin[1] && mp[1] <= step_back_rect.origin[1] + step_back_rect.widths[1] &&
                   mp[2] >= step_back_rect.origin[2] && mp[2] <= step_back_rect.origin[2] + step_back_rect.widths[2]
                    new_idx = max(1, frame_idx[] - 1)
                    frame_idx[] = new_idx
                    update_frame!(new_idx)
                    return Consume(true)
                end

                # Check step forward button
                if mp[1] >= step_forward_rect.origin[1] && mp[1] <= step_forward_rect.origin[1] + step_forward_rect.widths[1] &&
                   mp[2] >= step_forward_rect.origin[2] && mp[2] <= step_forward_rect.origin[2] + step_forward_rect.widths[2]
                    new_idx = min(n_frames, frame_idx[] + 1)
                    frame_idx[] = new_idx
                    update_frame!(new_idx)
                    return Consume(true)
                end
            elseif event.action == Mouse.release
                if slider_dragging[]
                    slider_dragging[] = false
                    return Consume(true)
                end
            end
        end
        return Consume(false)
    end

    on(events(ui_scene).mouseposition, priority=2) do mp
        if slider_dragging[]
            rel_pos = clamp((mp[1] - slider_x_start) / slider_width_obs[], 0.0, 1.0)
            new_idx = round(Int, 1 + rel_pos * (n_frames - 1))
            frame_idx[] = clamp(new_idx, 1, n_frames)
            update_frame!(frame_idx[])
        end
    end

    # Animation loop with replay speed
    @async begin
        while true
            if is_playing[]
                if frame_idx[] < n_frames
                    new_idx = frame_idx[] + 1
                    frame_idx[] = new_idx
                    update_frame!(new_idx)
                    # Calculate sleep time based on actual time difference and replay speed
                    if frame_idx[] > 1
                        dt = lg.syslog[frame_idx[]].time - lg.syslog[frame_idx[] - 1].time
                        sleep(max(0.01, dt / replay_speed))
                    else
                        sleep(0.05)
                    end
                elseif loop
                    frame_idx[] = 1
                    update_frame!(1)
                else
                    is_playing[] = false  # Stop at end
                end
            end
            sleep(0.02)  # Check state frequently
        end
    end

    return scene
end

end
