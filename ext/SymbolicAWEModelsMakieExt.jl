# Copyright (c) 2025 Bart van de Lint
# SPDX-License-Identifier: MPL-2.0

module SymbolicAWEModelsMakieExt

using Makie
using LinearAlgebra
using StaticArrays
using Statistics
using SymbolicAWEModels

# This is the core implementation. It plots a `SystemStructure` into a
# pre-existing `Axis3` or `LScene`. This makes the plotting logic composable.
function Makie.plot!(ax::Union{Axis3, LScene}, sys::SystemStructure;
                     point_color = :darkred, segment_color = :black,
                     wing_colors = Makie.wong_colors(), vector_scale = 1.0,
                     show_points = true, show_segments = true, show_orient = true)

    # === Plot Segments ===
    if show_segments
        lineseg_points = Point3f[]
        for seg in sys.segments
            p1 = sys.points[seg.point_idxs[1]].pos_w
            p2 = sys.points[seg.point_idxs[2]].pos_w
            push!(lineseg_points, Point3f(p1))
            push!(lineseg_points, Point3f(p2))
        end
        linesegments!(ax, lineseg_points, color=segment_color, label="Segments")
    end

    # === Plot Points ===
    if show_points
        point_positions = [Point3f(p.pos_w) for p in sys.points]
        scatter!(ax, point_positions, color=point_color, label="Points")
    end

    # === Plot Wings ===
    if show_orient
        for (i, wing) in enumerate(sys.wings)
            wing_pos = Point3f(wing.pos_w)
            color = wing_colors[mod1(i, length(wing_colors))]
            scatter!(ax, wing_pos, color=color, markersize=4, strokewidth=1,
                     strokecolor=:black, label="Wing $i")

            R = wing.R_b_w
            scale = vector_scale
            origins = [wing_pos, wing_pos, wing_pos]
            directions = [
                Vec3f(R[:, 1]) * scale,
                Vec3f(R[:, 2]) * scale,
                Vec3f(R[:, 3]) * scale
            ]

            axis_colors = [:red, :green, :blue]
            arrows3d!(ax, origins, directions, color=axis_colors, label="Wing $i Axes")
        end
    end

    return ax
end

# This is the top-level function that gets called when a user types `plot(sys)`.
# It creates a new Figure and Axis3, and then calls the `plot!` method above.
function Makie.plot(sys::SystemStructure; size = (1200, 800), zoom=false,
                    zoomsize=10, paddingsize=10, kwargs...)
    fig = Figure(; size)

    wing_pos = sys.wings[1].pos_w
    tether_len = mean([winch.tether_len for winch in sys.winches])
    if zoom
        limits = ((wing_pos[1]-zoomsize/2, wing_pos[1]+zoomsize/2),
                  (wing_pos[2]-zoomsize/2, wing_pos[2]+zoomsize/2),
                  (wing_pos[3]-zoomsize/2, wing_pos[3]+zoomsize/2))
    else
        limits = (zeros(2), zeros(2), zeros(2))
        for i in 1:3
            limits[i][1] = ifelse(wing_pos[i] > 0, -paddingsize, -tether_len-paddingsize)
            limits[i][2] = ifelse(wing_pos[i] > 0, tether_len+paddingsize, paddingsize)
        end
    end
    ax = Axis3(fig[1, 1]; aspect = :data,
               xlabel = "X", ylabel = "Y", zlabel = "Z",
               azimuth = 9/8*π, limits, zoommode = :cursor, viewmode = :fit)

    plot!(ax, sys; kwargs...)
    return fig
end

end
