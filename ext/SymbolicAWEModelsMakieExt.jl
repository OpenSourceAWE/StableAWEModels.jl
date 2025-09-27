# Copyright (c) 2025 Bart van de Lint
# SPDX-License-Identifier: MPL-2.0

module SymbolicAWEModelsMakieExt

using Makie
using LinearAlgebra
using StaticArrays
using SymbolicAWEModels

# This is a type recipe for SystemStructure. It defines how to plot the
# `SystemStructure` type by combining several of Makie's basic plotting functions.
# See: https://docs.makie.org/stable/explanations/recipes/#type_recipes

function Makie.plot!(plot::Plot{<:Any, Tuple{SystemStructure}})
    # Get the SystemStructure object as an Observable.
    # This allows the plot to update automatically if the SystemStructure changes.
    sys_obs = plot[1]

    # Manually provide default values for attributes if they are not passed by the user.
    point_color = haskey(plot, :point_color) ? plot.point_color : :gray
    segment_color = haskey(plot, :segment_color) ? plot.segment_color : :black
    wing_colors = haskey(plot, :wing_colors) ? plot.wing_colors : Makie.wong_colors()
    vector_scale = haskey(plot, :vector_scale) ? plot.vector_scale : 1.0
    show_points = haskey(plot, :show_points) ? plot.show_points : true
    show_segments = haskey(plot, :show_segments) ? plot.show_segments : true

    # === Plot Segments ===
    lineseg_points_obs = lift(sys_obs) do sys
        points = Point3f[]
        for seg in sys.segments
            p1 = sys.points[seg.point_idxs[1]].pos_w
            p2 = sys.points[seg.point_idxs[2]].pos_w
            push!(points, Point3f(p1))
            push!(points, Point3f(p2))
        end
        points
    end
    linesegments!(plot, lineseg_points_obs, color=segment_color, label="Segments", visible=show_segments)

    # === Plot Points ===
    point_positions_obs = lift(sys -> [Point3f(p.pos_w) for p in sys.points], sys_obs)
    scatter!(plot, point_positions_obs, color=point_color, label="Points", visible=show_points)

    # === Plot Wings ===
    sys_val = to_value(sys_obs)
    wc_val = to_value(wing_colors)

    for i in 1:length(sys_val.wings)
        wing_obs = lift(s -> s.wings[i], sys_obs)
        
        wing_pos_obs = lift(w -> Point3f(w.pos_w), wing_obs)
        color = wc_val[mod1(i, length(wc_val))]
        scatter!(plot, wing_pos_obs, color=color, markersize=20, strokewidth=2, strokecolor=:black, label="Wing $i")

        origins_obs = lift(p -> [p, p, p], wing_pos_obs)
        
        directions_obs = lift(wing_obs, vector_scale) do w, scale
            R = w.R_b_w
            [Vec3f(R[:, 1]) * scale, Vec3f(R[:, 2]) * scale, Vec3f(R[:, 3]) * scale]
        end
        
        axis_colors = [:red, :green, :blue]
        arrows3d!(plot, origins_obs, directions_obs, color=axis_colors, label="Wing $i Axes")
    end

    return plot
end

end
