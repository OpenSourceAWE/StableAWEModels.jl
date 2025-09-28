# Copyright (c) 2025 Bart van de Lint
# SPDX-License-Identifier: MPL-2.0

module SymbolicAWEModelsMakieExt

using Makie
using UnPack
using LinearAlgebra
using StaticArrays
using Statistics
using SymbolicAWEModels

function Makie.plot!(ax::Union{Axis3, LScene}, sys::SystemStructure;
                     point_color = :darkred, segment_color = :black,
                     wing_colors = Makie.wong_colors(), vector_scale = 1.0,
                     plot_points = true, plot_segments = true, plot_orient = true,
                     plot_vsm = true, kwargs...)

    # === Plot Segments ===
    if plot_segments
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
    if plot_points
        point_positions = [Point3f(p.pos_w) for p in sys.points]
        scatter!(ax, point_positions, color=point_color, label="Points")
    end

    # === Plot Wings ===
    if plot_orient
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

            if plot_vsm
                plot!(ax, wing.vsm_aero; R_b_w=wing.R_b_w, T_b_w=wing.pos_w, kwargs...)
            end
        end
    end

    return ax
end

# This is the top-level function that gets called when a user types `plot(sys)`.
# It creates a new Figure and Axis3, and then calls the `plot!` method above.
function Makie.plot(sys::SystemStructure; size = (1200, 800), zoom=false,
                    zoommargin=6, margin=[20, 20, 5], kwargs...)
    fig = Figure(; size)

    @unpack wings, winches, tethers, points = sys
    wing_pos = wings[1].pos_w
    if zoom
        limits = ((wing_pos[1]-zoommargin/2, wing_pos[1]+zoommargin/2),
                  (wing_pos[2]-zoommargin/2, wing_pos[2]+zoommargin/2),
                  (wing_pos[3]-zoommargin/2, wing_pos[3]+zoommargin/2))
    else
        limits = ([Inf, -Inf], [Inf, -Inf], [Inf, -Inf])
        for point in points
            for i in 1:3
                if point.pos_w[i] - margin[i] < limits[i][1]
                    limits[i][1] = point.pos_w[i] - margin[i]
                end
                if point.pos_w[i] + margin[i] > limits[i][2]
                    limits[i][2] = point.pos_w[i] + margin[i]
                end
            end
        end
    end
    ax = Axis3(fig[1, 1]; aspect = :data,
               xlabel = "X", ylabel = "Y", zlabel = "Z",
               azimuth = 9/8*π, zoommode = :cursor, viewmode = :fit,
               limits
           )

    plot!(ax, sys; kwargs...)
    return fig
end

end
