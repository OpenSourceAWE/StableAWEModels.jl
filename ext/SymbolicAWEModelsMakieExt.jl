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
const PLOT_GEOMETRY_OBS = Ref{Union{Nothing, Observable}}(nothing)  # Single trigger observable
const PLOT_SEGMENT_COLORS_OBS = Ref{Union{Nothing, Observable}}(nothing)  # Separate for force coloring
const PLOT_SCENE = Ref{Union{Nothing, Scene}}(nothing)
const PLOT_BACKGROUND_PANES = Ref{Union{Nothing, Vector}}(nothing)
const PLOT_MARGIN = Ref{Float64}(10.0)
const PLOT_RELEVANT_PLOTS = Ref{Union{Nothing, Vector}}(nothing)
const PLOT_SYSTEM_STRUCTURE = Ref{Union{Nothing, SystemStructure}}(nothing)
const PLOT_VECTOR_SCALE = Ref{Float64}(1.0)
const PLOT_FORCE_COLOR = Ref{Bool}(false)
const PLOT_SEGMENT_COLOR = Ref{Symbol}(:black)
const PLOT_ZOOMED_IN = Ref{Bool}(false)
const PLOT_ZOOM_RELMARGIN = Ref{Float64}(0.2)
const PLOT_ZOOM_SEGMENT_IDX = Ref{Int}(-1)  # Which segment we're zoomed into (-1 = none)
const PLOT_BODY_FRAME = Ref{Bool}(false)  # Whether body frame tracking is active
const PLOT_CAMERA_DISTANCE = Ref{Union{Nothing, Float64}}(nothing)  # Stored camera distance
const PLOT_PREV_BODY_FRAME = Ref{Bool}(false)  # Previous body frame state
const PLOT_PREV_ZOOMED_IN = Ref{Bool}(false)  # Previous zoomed state
const PLOT_PREV_SEGMENT_IDX = Ref{Int}(-1)  # Previous segment index

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
                     wing_colors = Makie.wong_colors(), vector_scale = 1.0,
                     show_points = true, show_segments = true, show_orient = true,
                     show_panes = true, margin = 10.0, force_color = false,
                     plot_vsm = true, plot_aero = true,
                     # Optional observable for real-time updates
                     geometry_obs = nothing,
                     pane_observables = nothing)

    plots = Dict{Symbol, Any}()

    # === Plot Segments ===
    if show_segments
        if isnothing(geometry_obs)
            # Static plotting: build segment points once
            lineseg_points = Point3f[]
            for seg in sys.segments
                p1 = sys.points[seg.point_idxs[1]].pos_w
                p2 = sys.points[seg.point_idxs[2]].pos_w
                push!(lineseg_points, Point3f(p1))
                push!(lineseg_points, Point3f(p2))
            end
        else
            # Dynamic plotting: compute from PLOT_SYSTEM_STRUCTURE when triggered
            lineseg_points = @lift begin
                $geometry_obs  # Trigger dependency
                sys_ref = PLOT_SYSTEM_STRUCTURE[]
                points = Point3f[]
                for seg in sys_ref.segments
                    p1 = sys_ref.points[seg.point_idxs[1]].pos_w
                    p2 = sys_ref.points[seg.point_idxs[2]].pos_w
                    push!(points, Point3f(p1))
                    push!(points, Point3f(p2))
                end
                points
            end
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
        if isnothing(geometry_obs)
            # Static plotting
            point_positions = [Point3f(p.pos_w) for p in sys.points]
        else
            # Dynamic plotting: compute from PLOT_SYSTEM_STRUCTURE when triggered
            point_positions = @lift begin
                $geometry_obs  # Trigger dependency
                [Point3f(p.pos_w) for p in PLOT_SYSTEM_STRUCTURE[].points]
            end
        end
        plots[:points] = scatter!(ax, point_positions, color=point_color, label="Points",
                                  transparency=true)
    end

    # === Plot Wings ===
    if show_orient
        if isnothing(geometry_obs)
            # Static plotting: create separate arrows for each wing
            plots[:wings] = []
            for (i, wing) in enumerate(sys.wings)
                wing_pos = Point3f(wing.pos_w)
                R = wing.R_b_w
                scale = vector_scale
                origins = [wing_pos, wing_pos, wing_pos]
                directions = [Vec3f(R[:, 1]) * scale, Vec3f(R[:, 2]) * scale, Vec3f(R[:, 3]) * scale]

                axis_colors = [:red, :green, :blue]
                p = arrows3d!(ax, origins, directions, color=axis_colors, label="Wing $i Axes")
                push!(plots[:wings], p)
            end
        else
            # Dynamic plotting: compute from PLOT_SYSTEM_STRUCTURE when triggered
            wing_origins_dirs = @lift begin
                $geometry_obs  # Trigger dependency
                sys_ref = PLOT_SYSTEM_STRUCTURE[]
                scale = PLOT_VECTOR_SCALE[]
                origins = Point3f[]
                directions = Vec3f[]
                for wing in sys_ref.wings
                    wing_pos = Point3f(wing.pos_w)
                    R = wing.R_b_w
                    # Add three arrow vectors for each axis (x, y, z in body frame)
                    for i in 1:3
                        push!(origins, wing_pos)
                        push!(directions, Vec3f(R[:, i]) * scale)
                    end
                end
                (origins, directions)
            end
            axis_colors = repeat([:red, :green, :blue], length(sys.wings))
            plots[:wings] = arrows3d!(ax, @lift($wing_origins_dirs[1]), @lift($wing_origins_dirs[2]),
                                     color=axis_colors)
        end
    end

    # === Plot VSM Aerodynamics ===
    # VSM panels - use observables if we're in dynamic mode
    if plot_vsm && !isempty(sys.wings)
        plots[:vsm] = []
        use_obs = !isnothing(geometry_obs)  # If geometry observable exists, use it for VSM too
        for (i, wing) in enumerate(sys.wings)
            p = plot!(ax, wing.vsm_aero; R_b_w=wing.R_b_w, T_b_w=wing.pos_w, use_observables=use_obs)
            push!(plots[:vsm], p)
        end
    end

    # === Plot Aero Forces ===
    if plot_aero
        if isnothing(geometry_obs)
            # Static plotting: build aero force arrows
            aero_origins = Point3f[]
            aero_forces_raw = Vec3f[]

            for wing in sys.wings
                if wing.wing_type == QUATERNION
                    # For QUATERNION wings, use wing.aero_force_b
                    if !iszero(wing.aero_force_b)
                        aero_force_w = wing.R_b_w * wing.aero_force_b
                        push!(aero_origins, Point3f(wing.pos_w))
                        push!(aero_forces_raw, Vec3f(aero_force_w))
                    end
                elseif wing.wing_type == SymbolicAWEModels.REFINE
                    # For REFINE wings, plot both point forces and total wing force
                    # Plot individual point forces
                    for point in sys.points
                        if point.type == WING && point.wing_idx == wing.idx
                            if !iszero(point.aero_force_b)
                                aero_force_w = wing.R_b_w * point.aero_force_b
                                push!(aero_origins, Point3f(point.pos_w))
                                push!(aero_forces_raw, Vec3f(aero_force_w))
                            end
                        end
                    end
                    # Also plot total wing aero force at wing center when vector_scale > 0
                    if vector_scale > 0 && !iszero(wing.aero_force_b)
                        aero_force_w = wing.R_b_w * wing.aero_force_b
                        push!(aero_origins, Point3f(wing.pos_w))
                        push!(aero_forces_raw, Vec3f(aero_force_w))
                    end
                end
            end

            if !isempty(aero_origins)
                # Calculate adaptive force scale
                max_force = maximum(norm.(aero_forces_raw))
                # Scale forces to be similar size as vector_scale
                force_scale = vector_scale / max_force
                aero_directions = [f * force_scale for f in aero_forces_raw]

                plots[:aero_forces] = arrows3d!(ax, aero_origins, aero_directions,
                                               color=:magenta,
                                               label="Aero Forces")
            end
        else
            # Dynamic plotting: compute from PLOT_SYSTEM_STRUCTURE when triggered
            aero_origins_dirs = @lift begin
                $geometry_obs  # Trigger dependency
                sys_ref = PLOT_SYSTEM_STRUCTURE[]
                scale = PLOT_VECTOR_SCALE[]
                origins = Point3f[]
                forces_raw = Vec3f[]

                for wing in sys_ref.wings
                    if wing.wing_type == QUATERNION
                        if !iszero(wing.aero_force_b)
                            aero_force_w = wing.R_b_w * wing.aero_force_b
                            push!(origins, Point3f(wing.pos_w))
                            push!(forces_raw, Vec3f(aero_force_w))
                        end
                    elseif wing.wing_type == SymbolicAWEModels.REFINE
                        for point in sys_ref.points
                            if point.type == WING && point.wing_idx == wing.idx
                                if !iszero(point.aero_force_b)
                                    aero_force_w = wing.R_b_w * point.aero_force_b
                                    push!(origins, Point3f(point.pos_w))
                                    push!(forces_raw, Vec3f(aero_force_w))
                                end
                            end
                        end
                        if scale > 0 && !iszero(wing.aero_force_b)
                            aero_force_w = wing.R_b_w * wing.aero_force_b
                            push!(origins, Point3f(wing.pos_w))
                            push!(forces_raw, Vec3f(aero_force_w))
                        end
                    end
                end

                # Calculate adaptive force scale
                directions = Vec3f[]
                if !isempty(forces_raw)
                    max_force = maximum(norm.(forces_raw))
                    force_scale = scale / max_force
                    directions = [f * force_scale for f in forces_raw]
                end
                (origins, directions)
            end

            plots[:aero_forces] = arrows3d!(ax, @lift($aero_origins_dirs[1]), @lift($aero_origins_dirs[2]),
                                           color=:magenta)
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
        # Use fixed radii so arrows don't get fatter with longer tethers
        plots[:global_axes] = arrows3d!(ax, origins, directions;
                                        shaftradius=0.02,
                                        tipradius=0.08,
                                        tiplength=0.2,
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
    SymbolicAWEModels.update_plot_observables!(sys::SystemStructure)

Trigger plot updates by updating the geometry observable.

The SystemStructure should already be updated via `update_from_sysstate!`.
This function simply triggers the observable which causes Makie to recompute
all geometry from `PLOT_SYSTEM_STRUCTURE[]` via `@lift` expressions.

# Example
```julia
# Create initial plot
scene = plot(sys_struct)

# In simulation loop:
for step in 1:steps
    next_step!(sam; ...)
    update_plot_observables!(sam.sys_struct)
    sleep(0.001)  # Allow Makie to process updates
end
```
"""
function SymbolicAWEModels.update_plot_observables!(sys::SystemStructure)
    # Trigger geometry observable - this causes all @lift expressions to recompute
    if !isnothing(PLOT_GEOMETRY_OBS[])
        PLOT_GEOMETRY_OBS[][] = time()  # Use timestamp as trigger value
    end

    # Update segment colors if force coloring is enabled
    if !isnothing(PLOT_SEGMENT_COLORS_OBS[]) && PLOT_FORCE_COLOR[]
        PLOT_SEGMENT_COLORS_OBS[][] = calculate_segment_force_colors(sys.segments, PLOT_SEGMENT_COLOR[])
    end

    # Update VSM panel meshes
    if !isnothing(PLOT_GEOMETRY_OBS[])
        for wing in sys.wings
            plot!(wing.vsm_aero; R_b_w=wing.R_b_w, T_b_w=wing.pos_w)
        end
    end

    # Auto-update camera to keep geometry centered (runs every frame during replay)
    # Re-runs the appropriate zoom function to track moving geometry
    if !isnothing(PLOT_SCENE[]) && !isnothing(PLOT_RELEVANT_PLOTS[]) && !isnothing(PLOT_SYSTEM_STRUCTURE[])
        scene = PLOT_SCENE[]
        relevant_plots = PLOT_RELEVANT_PLOTS[]
        stored_sys = PLOT_SYSTEM_STRUCTURE[]

        # Detect mode change
        mode_changed = (PLOT_PREV_BODY_FRAME[] != PLOT_BODY_FRAME[] ||
                        PLOT_PREV_ZOOMED_IN[] != PLOT_ZOOMED_IN[] ||
                        PLOT_PREV_SEGMENT_IDX[] != PLOT_ZOOM_SEGMENT_IDX[])

        if mode_changed
            PLOT_CAMERA_DISTANCE[] = nothing  # Force recalculation
            # Update prev state
            PLOT_PREV_BODY_FRAME[] = PLOT_BODY_FRAME[]
            PLOT_PREV_ZOOMED_IN[] = PLOT_ZOOMED_IN[]
            PLOT_PREV_SEGMENT_IDX[] = PLOT_ZOOM_SEGMENT_IDX[]
        end

        # Call zoom functions with stored distance (preserves manual zoom)
        if PLOT_BODY_FRAME[]
            # Body frame mode: continuously track wing orientation
            dist = zoom_body_frame!(scene, scene.camera, stored_sys, PLOT_CAMERA_DISTANCE[])
            PLOT_CAMERA_DISTANCE[] = dist
        elseif PLOT_ZOOMED_IN[] && PLOT_ZOOM_SEGMENT_IDX[] > 0
            # When zoomed in, keep segment centered as it moves
            dist = zoom_in!(scene, scene.camera, stored_sys, PLOT_ZOOM_SEGMENT_IDX[], PLOT_CAMERA_DISTANCE[])
            PLOT_CAMERA_DISTANCE[] = dist
        elseif !PLOT_ZOOMED_IN[]
            # When not zoomed in, keep full view centered on geometry
            dist = zoom_out!(scene, scene.camera, relevant_plots, PLOT_CAMERA_DISTANCE[]; relmargin=PLOT_ZOOM_RELMARGIN[])
            PLOT_CAMERA_DISTANCE[] = dist
        end
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
- `plot_twist::Bool=false`: Show the panel with the twist angles for each wing group.
- `plot_aoa::Bool=plot_default`: Show the panel with the angle of attack.
- `plot_heading::Bool=plot_default`: Show the panel with the kite's heading and course angles.
- `plot_kiteutils_course::Bool=false`: Also plot course calculated using KiteUtils.calc_course.
- `plot_elevation::Bool=false`: Show the panel with the kite's elevation angle.
- `plot_azimuth::Bool=false`: Show the panel with the kite's azimuth angle.
- `plot_distance::Bool=false`: Show the panel with the kite distance from origin (norm of position).
- `plot_cone_angle::Bool=false`: Show the panel with the cone angle (angle between wind vector and normalized kite position).
- `plot_old_heading::Bool=false`: Show the old heading calculated from orientation quaternion (angle between -R_b_w[:,1] and -R_v_w[:,1]).
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
function Makie.plot(sys::SystemStructure, lg::SysLog{N}; kwargs...) where N
    # Wrapper that uses the vector version
    return Makie.plot([sys], [lg]; kwargs...)
end

"""
    Makie.plot(sys::SystemStructure, logs::Vector{SysLog}; kwargs...)

Create a multi-panel plot comparing multiple simulation logs on the same figure.

This method allows plotting multiple syslogs (e.g., from REFINE and QUATERNION models)
on the same panels for direct comparison. Each log's traces are labeled with its
corresponding system name.

# Arguments
- `sys::SystemStructure`: The system structure (can be from any of the models).
- `logs::Vector{SysLog}`: Vector of simulation logs to compare.

# Keyword Arguments
Same as the single-syslog version. See `Makie.plot(sys::SystemStructure, lg::SysLog)` for details.

# Example
```julia
# Compare REFINE vs QUATERNION models
plot(sys_struct, [syslog_refine, syslog_quat];
     plot_turn_rates=true, plot_azimuth=true,
     plot_heading=true, plot_aoa=true,
     plot_default=false, plot_aero_force=true)
```
"""
function Makie.plot(sys::SystemStructure, logs::Vector{<:SysLog}; kwargs...)
    # Wrapper that creates a vector of SystemStructure with the same sys for each log
    syss = [sys for _ in logs]
    return Makie.plot(syss, logs; kwargs...)
end

function Makie.plot(syss::Vector{SystemStructure}, logs::Vector{<:SysLog};
                   plot_default=true,
                   plot_reelout=plot_default,
                   plot_aero_force=plot_default,
                   plot_twist=false,
                   plot_aoa=plot_default,
                   plot_heading=plot_default,
                   plot_kiteutils_course=false,
                   plot_aero_moment=false,
                   plot_turn_rates=false,
                   plot_elevation=false,
                   plot_azimuth=false,
                   plot_tether_moment=false,
                   plot_winch_force=plot_default,
                   plot_set_values=false,
                   plot_distance=false,
                   plot_cone_angle=false,
                   plot_old_heading=false,
                   plot_tether=false,
                   heading_setpoint=nothing,
                   tether_len_setpoint=nothing,
                   size=(1200, 800))

    # Build list of panels to plot by combining data from all logs
    panels = []

    if plot_turn_rates
        all_data = []
        all_labels = []
        all_times = []
        for (i, lg) in enumerate(logs)
            sl = lg.syslog
            suffix = " - " * syss[i].name

            # Calculate heading rate from diff for quaternion wings
            heading_unwrapped = copy(sl.heading)
            for j in 2:length(heading_unwrapped)
                while heading_unwrapped[j] - heading_unwrapped[j-1] > π
                    heading_unwrapped[j] -= 2π
                end
                while heading_unwrapped[j] - heading_unwrapped[j-1] < -π
                    heading_unwrapped[j] += 2π
                end
            end
            heading_rate = diff(rad2deg.(heading_unwrapped)) ./ diff(sl.time)
            push!(all_data, heading_rate)
            push!(all_labels, "dψ/dt" * suffix)
            push!(all_times, sl.time[1:end-1])

            # Also plot z-component (yaw rate ω_z) if available
            turn_rates_rad = hcat(sl.turn_rates...)
            if !all(iszero, turn_rates_rad)
                omega_z_deg = rad2deg.(turn_rates_rad[3, :])
                push!(all_data, omega_z_deg)
                push!(all_labels, "ω_z" * suffix)
                push!(all_times, sl.time)
            end
        end
        push!(panels, (
            data = all_data,
            labels = all_labels,
            times = all_times,
            ylabel = "turn rate [°/s]"
        ))
    end

    if plot_reelout
        all_data = []
        all_labels = []
        all_times = []
        for (i, lg) in enumerate(logs)
            sl = lg.syslog
            suffix = " - " * syss[i].name
            v_ro = [[sl.v_reelout[i][j] for i in eachindex(sl.v_reelout)] for j in 1:3]
            for j in 1:3
                # Only plot if non-zero or if it's index 1
                if j == 1 || !all(iszero, v_ro[j])
                    push!(all_data, v_ro[j])
                    push!(all_labels, "v_ro,$j" * suffix)
                    push!(all_times, sl.time)
                end
            end
        end
        if !isempty(all_data)
            push!(panels, (
                data = all_data,
                labels = all_labels,
                times = all_times,
                ylabel = "reel-out speed [m/s]"
            ))
        end
    end

    if plot_tether
        all_data = []
        all_labels = []
        all_times = []
        for (i, lg) in enumerate(logs)
            sl = lg.syslog
            suffix = " - " * syss[i].name
            # Extract tether length (first element if vector, otherwise scalar)
            l_tether = [length(sl.l_tether[j]) > 0 ? sl.l_tether[j][1] : 0.0
                        for j in eachindex(sl.l_tether)]
            push!(all_data, l_tether)
            push!(all_labels, "l_tether" * suffix)
            push!(all_times, sl.time)

            # Add setpoint if provided
            if !isnothing(tether_len_setpoint)
                is_multi_setpoint = (tether_len_setpoint isa Vector &&
                                     length(tether_len_setpoint) > 0 &&
                                     tether_len_setpoint[1] isa AbstractVector)

                if is_multi_setpoint
                    if i <= length(tether_len_setpoint) && !isnothing(tether_len_setpoint[i])
                        push!(all_data, tether_len_setpoint[i])
                        push!(all_labels, "l_sp" * suffix)
                        push!(all_times, sl.time)
                    end
                else
                    push!(all_data, tether_len_setpoint)
                    push!(all_labels, "l_sp" * suffix)
                    push!(all_times, sl.time)
                end
            end
        end
        push!(panels, (
            data = all_data,
            labels = all_labels,
            times = all_times,
            ylabel = "tether length [m]"
        ))
    end

    if plot_aero_force
        all_data = []
        all_labels = []
        all_times = []
        for (i, lg) in enumerate(logs)
            sl = lg.syslog
            suffix = " - " * syss[i].name
            aero_force_z = [sl.aero_force_b[i][3] for i in eachindex(sl.aero_force_b)]
            push!(all_data, aero_force_z)
            push!(all_labels, "F_aero,z" * suffix)
            push!(all_times, sl.time)
        end
        push!(panels, (
            data = all_data,
            labels = all_labels,
            times = all_times,
            ylabel = "F_aero,z [N]"
        ))
    end

    if plot_aero_moment
        all_data = []
        all_labels = []
        all_times = []
        for (i, lg) in enumerate(logs)
            sl = lg.syslog
            suffix = " - " * syss[i].name
            aero_moment_z = [sl.aero_moment_b[i][3] for i in eachindex(sl.aero_moment_b)]
            push!(all_data, aero_moment_z)
            push!(all_labels, "M_aero,z" * suffix)
            push!(all_times, sl.time)
        end
        push!(panels, (
            data = all_data,
            labels = all_labels,
            times = all_times,
            ylabel = "M_aero,z [Nm]"
        ))
    end

    if plot_tether_moment
        all_data = []
        all_labels = []
        all_times = []
        for (i, lg) in enumerate(logs)
            sl = lg.syslog
            suffix = " - " * syss[i].name
            tether_moment_y = [sl.tether_induced_moment[i][2] for i in eachindex(sl.tether_induced_moment)]
            push!(all_data, tether_moment_y)
            push!(all_labels, "M_tether,y" * suffix)
            push!(all_times, sl.time)
        end
        push!(panels, (
            data = all_data,
            labels = all_labels,
            times = all_times,
            ylabel = "M_tether,y [Nm]"
        ))
    end

    if plot_twist
        all_data = []
        all_labels = []
        all_times = []
        for (i, lg) in enumerate(logs)
            sl = lg.syslog
            suffix = " - " * syss[i].name
            n_groups = length(sl.twist_angles[1])
            if n_groups > 0
                twist_deg = rad2deg.(hcat(sl.twist_angles...))
                for j in 1:n_groups
                    push!(all_data, twist_deg[j, :])
                    push!(all_labels, "β_$j" * suffix)
                    push!(all_times, sl.time)
                end
            end
        end
        if !isempty(all_data)
            push!(panels, (
                data = all_data,
                labels = all_labels,
                times = all_times,
                ylabel = "twist [°]"
            ))
        end
    end

    if plot_aoa
        all_data = []
        all_labels = []
        all_times = []
        for (i, lg) in enumerate(logs)
            sl = lg.syslog
            suffix = " - " * syss[i].name
            aoa_deg = rad2deg.(sl.AoA)
            push!(all_data, aoa_deg)
            push!(all_labels, "α" * suffix)
            push!(all_times, sl.time)
        end
        push!(panels, (
            data = all_data,
            labels = all_labels,
            times = all_times,
            ylabel = "angle of attack [°]"
        ))
    end

    if plot_heading
        all_data = []
        all_labels = []
        all_times = []
        for (i, lg) in enumerate(logs)
            sl = lg.syslog
            suffix = " - " * syss[i].name
            heading_deg = rad2deg.(sl.heading)
            course_deg = rad2deg.(sl.course)
            push!(all_data, heading_deg)
            push!(all_labels, "ψ" * suffix)
            push!(all_times, sl.time)
            push!(all_data, course_deg)
            push!(all_labels, "χ" * suffix)
            push!(all_times, sl.time)
            if plot_kiteutils_course
                course_kiteutils_deg = [rad2deg(KiteUtils.calc_course(sl.orient[i])) for i in eachindex(sl.orient)]
                push!(all_data, course_kiteutils_deg)
                push!(all_labels, "χ_KU" * suffix)
                push!(all_times, sl.time)
            end
            # Add setpoint if provided
            if !isnothing(heading_setpoint)
                # Check if heading_setpoint is a vector of vectors (multiple logs)
                # or a vector of numbers (single log)
                is_multi_setpoint = (heading_setpoint isa Vector &&
                                     length(heading_setpoint) > 0 &&
                                     heading_setpoint[1] isa AbstractVector)

                if is_multi_setpoint
                    # Multiple setpoints (one per log)
                    if i <= length(heading_setpoint) && !isnothing(heading_setpoint[i])
                        setpoint_deg = rad2deg.(heading_setpoint[i])
                        push!(all_data, setpoint_deg)
                        push!(all_labels, "ψ_sp" * suffix)
                        push!(all_times, sl.time)
                    end
                else
                    # Single setpoint for all logs (or single log)
                    setpoint_deg = rad2deg.(heading_setpoint)
                    push!(all_data, setpoint_deg)
                    push!(all_labels, "ψ_sp" * suffix)
                    push!(all_times, sl.time)
                end
            end
        end
        push!(panels, (
            data = all_data,
            labels = all_labels,
            times = all_times,
            ylabel = "heading/course [°]"
        ))
    end

    if plot_old_heading
        all_data = []
        all_labels = []
        all_times = []
        for (i, lg) in enumerate(logs)
            sl = lg.syslog
            suffix = " - " * syss[i].name
            # Calculate old heading from orientation quaternion
            old_heading_rad = similar(sl.heading)
            for i in eachindex(sl.orient)
                R_b_w = KiteUtils.rotation_matrix(sl.orient[i])
                R_v_w = [1.0 0.0 0.0; 0.0 -1.0 0.0; 0.0 0.0 -1.0]
                v1 = -R_b_w[:, 1]
                v2 = -R_v_w[:, 1]
                old_heading_rad[i] = atan(v1[2] - v2[2], v1[1] - v2[1])
            end
            old_heading_deg = rad2deg.(old_heading_rad)
            push!(all_data, old_heading_deg)
            push!(all_labels, "ψ_old" * suffix)
            push!(all_times, sl.time)
        end
        push!(panels, (
            data = all_data,
            labels = all_labels,
            times = all_times,
            ylabel = "old heading [°]"
        ))
    end

    if plot_distance
        all_data = []
        all_labels = []
        all_times = []
        for (i, lg) in enumerate(logs)
            sl = lg.syslog
            suffix = " - " * syss[i].name
            # Get wing origin index
            kite_idx = syss[i].wings[1].origin_idx
            distance = [norm([sl.X[j][kite_idx], sl.Y[j][kite_idx], sl.Z[j][kite_idx]]) for j in eachindex(sl.X)]
            push!(all_data, distance)
            push!(all_labels, "r" * suffix)
            push!(all_times, sl.time)
        end
        push!(panels, (
            data = all_data,
            labels = all_labels,
            times = all_times,
            ylabel = "distance [m]"
        ))
    end

    if plot_cone_angle
        all_data = []
        all_labels = []
        all_times = []
        for (i, lg) in enumerate(logs)
            sl = lg.syslog
            suffix = " - " * syss[i].name
            # Get wing origin index
            kite_idx = syss[i].wings[1].origin_idx
            # Assuming wind_vec_gnd is available in syslog
            cone_angle_rad = similar(sl.heading)
            for j in eachindex(sl.X)
                pos = [sl.X[j][kite_idx], sl.Y[j][kite_idx], sl.Z[j][kite_idx]]
                pos_norm = normalize(pos)
                wind_norm = normalize([sl.v_wind_gnd[1], sl.v_wind_gnd[2], sl.v_wind_gnd[3]])
                cone_angle_rad[j] = acos(clamp(dot(pos_norm, wind_norm), -1.0, 1.0))
            end
            cone_angle_deg = rad2deg.(cone_angle_rad)
            push!(all_data, cone_angle_deg)
            push!(all_labels, "θ_cone" * suffix)
            push!(all_times, sl.time)
        end
        push!(panels, (
            data = all_data,
            labels = all_labels,
            times = all_times,
            ylabel = "cone angle [°]"
        ))
    end

    if plot_elevation
        all_data = []
        all_labels = []
        all_times = []
        for (i, lg) in enumerate(logs)
            sl = lg.syslog
            suffix = " - " * syss[i].name
            elevation_deg = rad2deg.(sl.elevation)
            push!(all_data, elevation_deg)
            push!(all_labels, "elevation" * suffix)
            push!(all_times, sl.time)
        end
        push!(panels, (
            data = all_data,
            labels = all_labels,
            times = all_times,
            ylabel = "elevation [°]"
        ))
    end

    if plot_azimuth
        all_data = []
        all_labels = []
        all_times = []
        for (i, lg) in enumerate(logs)
            sl = lg.syslog
            suffix = " - " * syss[i].name
            azimuth_deg = rad2deg.(sl.azimuth)
            push!(all_data, azimuth_deg)
            push!(all_labels, "azimuth" * suffix)
            push!(all_times, sl.time)
        end
        push!(panels, (
            data = all_data,
            labels = all_labels,
            times = all_times,
            ylabel = "azimuth [°]"
        ))
    end

    if plot_winch_force
        all_data = []
        all_labels = []
        all_times = []
        for (i, lg) in enumerate(logs)
            sl = lg.syslog
            suffix = " - " * syss[i].name
            winch_force = [[sl.winch_force[i][j] for i in eachindex(sl.winch_force)] for j in 1:3]
            for j in 1:3
                # Only plot if non-zero or if it's index 1
                if j == 1 || !all(iszero, winch_force[j])
                    push!(all_data, winch_force[j])
                    push!(all_labels, "F_winch,$j" * suffix)
                    push!(all_times, sl.time)
                end
            end
        end
        if !isempty(all_data)
            push!(panels, (
                data = all_data,
                labels = all_labels,
                times = all_times,
                ylabel = "Winch force [N]"
            ))
        end
    end

    if plot_set_values
        all_data = []
        all_labels = []
        all_times = []
        for (i, lg) in enumerate(logs)
            sl = lg.syslog
            suffix = " - " * syss[i].name
            set_values = [[sl.set_torque[i][j] for i in eachindex(sl.set_torque)] for j in 1:3]
            for j in 1:3
                # Only plot if non-zero or if it's index 1
                if j == 1 || !all(iszero, set_values[j])
                    push!(all_data, set_values[j])
                    push!(all_labels, "T_winch,$j" * suffix)
                    push!(all_times, sl.time)
                end
            end
        end
        if !isempty(all_data)
            push!(panels, (
                data = all_data,
                labels = all_labels,
                times = all_times,
                ylabel = "Set torque [Nm]"
            ))
        end
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
        for (j, (data_series, label, time_vec)) in enumerate(zip(panel.data, panel.labels, panel.times))
            if length(data_series) != length(time_vec)
                @warn "Skipping trace '$label': data length $(length(data_series)) != time length $(length(time_vec))"
                continue
            end
            lines!(ax, time_vec, data_series, label=label)
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

function zoom_out!(scene, cam, plots, distance=nothing; relmargin=0.2)
    # --- ROBUST ZOOM OUT ---
    # 1. Get the current camera viewing direction vector
    inv_view_matrix = inv(cam.view[])
    cam_dir_vec = normalize(Vec3f(inv_view_matrix[1, 3],
                                  inv_view_matrix[2, 3],
                                  inv_view_matrix[3, 3]))
    # 2. Get the scene's bounding box and its center (the new target)
    bbox = data_limits(plots)
    center = bbox.origin .+ bbox.widths ./ 2

    # 3. Calculate distance only if not provided
    if isnothing(distance)
        # Calculate the distance needed to see the whole box
        radius = norm(bbox.widths) / 2.0
        fov_rad = 2 * atan(1 / cam.projection[][2, 2])
        distance = radius / tan(fov_rad / 2.0) * (1 + relmargin)
    end

    # 4. Calculate the new camera position
    new_eyepos = center + cam_dir_vec * distance
    # 5. Update the camera to the new "fit-all" view
    update_cam!(scene, new_eyepos, center)

    return distance
end

function zoom_in!(scene, cam, sys, segment_idx, distance=nothing)
    # --- ZOOM IN ON SEGMENT ---
    # Get current segment endpoints
    seg = sys.segments[segment_idx]
    p1_w = sys.points[seg.point_idxs[1]].pos_w
    p2_w = sys.points[seg.point_idxs[2]].pos_w

    # Calculate segment center
    center = (p1_w + p2_w) / 2.0f0

    # Calculate distance only if not provided
    if isnothing(distance)
        segment_len = norm(p2_w - p1_w)
        distance = segment_len * 1.5 + 2.0
    end

    # Get camera direction
    inv_view_matrix = inv(cam.view[])
    cam_dir_vec = normalize(Vec3f(inv_view_matrix[1, 3],
                                  inv_view_matrix[2, 3],
                                  inv_view_matrix[3, 3]))

    # Calculate new camera position
    new_eyepos = center + distance * cam_dir_vec

    # Update camera
    update_cam!(scene, new_eyepos, center)

    return distance
end

"""
    zoom_body_frame!(scene, cam, sys, distance=nothing)

Set camera to body-frame view tracking the wing orientation.
Camera positioned behind the kite looking forward along body x-axis.

# Arguments
- `scene`: The Makie scene
- `cam`: The camera object
- `sys`: SystemStructure with wing to track
- `distance`: Optional fixed distance to use (if nothing, calculates from geometry)

# Returns
- Camera distance used (for storage/reuse)
"""
function zoom_body_frame!(scene, cam, sys, distance=nothing)
    if isempty(sys.wings)
        @warn "No wings in system, cannot use body frame view"
        return nothing
    end

    wing = sys.wings[1]
    kite_pos = wing.pos_w
    R_b_w = wing.R_b_w

    # Calculate distance only if not provided
    if isnothing(distance)
        # Calculate characteristic system length
        if !isempty(sys.points)
            all_x = [p.pos_w[1] for p in sys.points]
            all_y = [p.pos_w[2] for p in sys.points]
            all_z = [p.pos_w[3] for p in sys.points]

            xlims = extrema(all_x)
            ylims = extrema(all_y)
            zlims = extrema(all_z)

            char_length = max(xlims[2] - xlims[1], ylims[2] - ylims[1], zlims[2] - zlims[1])
        else
            char_length = 10.0
        end
        distance = char_length * 0.5
    end

    # Camera position: kite_pos - R_b_w * [distance, 0, 0]
    # Places camera in front of kite (negative x in body frame), looking along +x axis
    cam_offset_body = [-distance, 0.0, 0.0]
    cam_offset_world = R_b_w * cam_offset_body
    cam_pos = kite_pos + cam_offset_world

    # Set up vector to align with body z-axis BEFORE updating camera
    cam.upvector[] = Vec3f(R_b_w[:, 3])

    # Set camera position and lookat
    update_cam!(scene, Vec3f(cam_pos), Vec3f(kite_pos))

    return distance
end

function _plot_with_panes(sys::SystemStructure;
                    size = (1200, 800),
                    margin = 10.0,
                    relmargin = 0.2,
                    segment_color = :black,
                    highlight_color = :red,
                    force_color = false,
                    body_frame = false,
                    perspective = false,
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
                        PLOT_ZOOMED_IN[] = true  # Track global zoom state
                        PLOT_ZOOM_SEGMENT_IDX[] = hover_idx  # Track which segment we're zoomed into

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
                    zoom_out!(scene, cam, relevant_plots, nothing; relmargin)
                    zoomed_in[] = false
                    PLOT_ZOOMED_IN[] = false  # Track global zoom state
                    PLOT_ZOOM_SEGMENT_IDX[] = -1  # Clear zoomed segment

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

    # Set camera projection type
    cam = scene.camera
    if !perspective
        # Use orthographic projection
        # Get current projection view bounds
        widths = scene.viewport[].widths
        w_half = Float32(widths[1] / 2)
        h_half = Float32(widths[2] / 2)
        cam.projection[] = Makie.orthographicprojection(
            -w_half, w_half,
            -h_half, h_half,
            -10_000f0, 10_000f0
        )
    end

    # Set initial camera position
    if body_frame
        zoom_body_frame!(scene, cam, sys)
    else
        update_cam!(scene, Vec3f(-100, -100, 100), Vec3f(0, 0, 0))
        zoom_out!(scene, cam, relevant_plots, nothing; relmargin)
    end

    # Return scene along with pane_observables, margin, plots dict, and relevant_plots
    # These will be used by time-based plotting
    return scene, pane_observables, margin, plots, relevant_plots
end

# Public API function - creates scene with observables for dynamic updates
function Makie.plot(sys::SystemStructure;
                    vector_scale=1.0,
                    force_color=false,
                    segment_color=:black,
                    plot_aero=true,
                    relmargin=0.2,
                    body_frame=false,
                    perspective=false,
                    kwargs...)
    # Store SystemStructure globally FIRST so @lift expressions can access it
    PLOT_SYSTEM_STRUCTURE[] = sys
    PLOT_VECTOR_SCALE[] = vector_scale
    PLOT_FORCE_COLOR[] = force_color
    PLOT_SEGMENT_COLOR[] = segment_color
    PLOT_BODY_FRAME[] = body_frame

    # Create single geometry trigger observable
    geometry_obs = Observable(0.0)

    # Create scene with observables using internal function
    scene, pane_observables, margin, plots, relevant_plots = _plot_with_panes(sys;
                geometry_obs,
                vector_scale,
                force_color,
                segment_color,
                plot_aero,
                relmargin,
                body_frame,
                perspective,
                kwargs...)

    # Get the segment colors observable from the plots if available
    segment_colors_obs = nothing
    if haskey(plots, :segment_colors_obs)
        segment_colors_obs = plots[:segment_colors_obs]
    end

    # Store observables and settings globally
    PLOT_GEOMETRY_OBS[] = geometry_obs
    PLOT_SEGMENT_COLORS_OBS[] = segment_colors_obs
    PLOT_SCENE[] = scene
    PLOT_BACKGROUND_PANES[] = pane_observables
    PLOT_MARGIN[] = margin
    PLOT_RELEVANT_PLOTS[] = relevant_plots
    # SystemStructure and settings already stored above before plot creation
    PLOT_ZOOMED_IN[] = false  # Initialize zoom state (not zoomed in)
    PLOT_ZOOM_RELMARGIN[] = relmargin  # Store relmargin for auto-updates
    PLOT_ZOOM_SEGMENT_IDX[] = -1  # No segment zoomed initially

    return scene
end

"""
    Makie.plot!(sys::SystemStructure; vector_scale=1.0)

Update the currently displayed SystemStructure plot with new data from `sys`.

This function follows standard Makie conventions: `plot!` with `!` mutates the existing
scene by updating its observables. Must be called after an initial `plot(sys)` has created
the scene and observables.

# Arguments
- `sys::SystemStructure`: The system structure with updated state to display

# Keyword Arguments
- `vector_scale::Real=1.0`: Scale factor for wing orientation arrows

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
function Makie.plot!(sys::SystemStructure; vector_scale=1.0)
    # Check if geometry observable exists
    if isnothing(PLOT_GEOMETRY_OBS[]) || isnothing(PLOT_SCENE[])
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

    # Update vector scale if changed
    if vector_scale != PLOT_VECTOR_SCALE[]
        PLOT_VECTOR_SCALE[] = vector_scale
    end

    # Trigger observable update
    update_plot_observables!(sys)

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
- `vector_scale::Real=1.0`: Scale factor for wing orientation arrows
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
- `vector_scale::Real=1.0`: Scale factor for wing orientation arrows
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
                      vector_scale=1.0,
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

    # Body frame toggle button
    body_frame_button_x = step_forward_x + button_width + 10
    body_frame_button_rect = Rect2f(body_frame_button_x, button_y, button_width, button_height)
    body_frame_obs = Observable(PLOT_BODY_FRAME[])
    body_frame_button_color = @lift($(body_frame_obs) ? RGBAf(0.3, 0.6, 0.8, 0.8) : RGBAf(0.5, 0.5, 0.5, 0.8))
    poly!(ui_scene, body_frame_button_rect, color=body_frame_button_color)
    body_frame_button_label = @lift($(body_frame_obs) ? "Body" : "World")
    text!(ui_scene, body_frame_button_label, position=Point2f(body_frame_button_x + button_width/2, button_y + button_height/2),
          align=(:center, :center), fontsize=14, color=:white)

    # Info label
    info_label_x = body_frame_button_x + button_width + 20
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

                # Check body frame toggle button
                if mp[1] >= body_frame_button_rect.origin[1] && mp[1] <= body_frame_button_rect.origin[1] + body_frame_button_rect.widths[1] &&
                   mp[2] >= body_frame_button_rect.origin[2] && mp[2] <= body_frame_button_rect.origin[2] + body_frame_button_rect.widths[2]
                    PLOT_BODY_FRAME[] = !PLOT_BODY_FRAME[]
                    body_frame_obs[] = PLOT_BODY_FRAME[]
                    # Disable zoom-in mode when switching to body frame
                    if PLOT_BODY_FRAME[]
                        PLOT_ZOOMED_IN[] = false
                        PLOT_ZOOM_SEGMENT_IDX[] = -1
                    end
                    # Force camera update immediately (not just during playback)
                    PLOT_CAMERA_DISTANCE[] = nothing  # Force recalculation on mode change
                    if !isnothing(PLOT_SCENE[]) && !isnothing(PLOT_RELEVANT_PLOTS[]) && !isnothing(PLOT_SYSTEM_STRUCTURE[])
                        scene = PLOT_SCENE[]
                        relevant_plots = PLOT_RELEVANT_PLOTS[]
                        stored_sys = PLOT_SYSTEM_STRUCTURE[]

                        if PLOT_BODY_FRAME[]
                            dist = zoom_body_frame!(scene, scene.camera, stored_sys, PLOT_CAMERA_DISTANCE[])
                            PLOT_CAMERA_DISTANCE[] = dist
                        else
                            dist = zoom_out!(scene, scene.camera, relevant_plots, PLOT_CAMERA_DISTANCE[]; relmargin=PLOT_ZOOM_RELMARGIN[])
                            PLOT_CAMERA_DISTANCE[] = dist
                        end
                    end
                    # Update tracking variables to prevent mode change detection
                    PLOT_PREV_BODY_FRAME[] = PLOT_BODY_FRAME[]
                    PLOT_PREV_ZOOMED_IN[] = PLOT_ZOOMED_IN[]
                    PLOT_PREV_SEGMENT_IDX[] = PLOT_ZOOM_SEGMENT_IDX[]
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
