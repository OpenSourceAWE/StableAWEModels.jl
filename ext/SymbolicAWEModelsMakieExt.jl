# Copyright (c) 2025 Bart van de Lint
# SPDX-License-Identifier: MPL-2.0

module SymbolicAWEModelsMakieExt

using Makie
using LinearAlgebra
using StaticArrays
using SymbolicAWEModels

# This is the core implementation. It plots a `SystemStructure` into a
# pre-existing `Axis3` or `LScene`. This makes the plotting logic composable.
function Makie.plot!(ax::Union{Axis3, LScene}, sys::SystemStructure;
                     point_color = :darkred, segment_color = :black,
                     wing_colors = Makie.wong_colors(), vector_scale = 3.0,
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
            scatter!(ax, wing_pos, color=color, markersize=4, strokewidth=1, strokecolor=:black, label="Wing $i")

            R = wing.R_b_w
            scale = vector_scale
            origins = [wing_pos, wing_pos, wing_pos]
            directions = [Vec3f(R[:, 1]) * scale, Vec3f(R[:, 2]) * scale, Vec3f(R[:, 3]) * scale]

            axis_colors = [:red, :green, :blue]
            arrows3d!(ax, origins, directions, color=axis_colors, label="Wing $i Axes")
        end
    end
    
    return ax
end

# This is the top-level function that gets called when a user types `plot(sys)`.
# It creates a new Figure and LScene, and then calls the `plot!` method above.
function Makie.plot(sys::SystemStructure; size = (1200, 800), kwargs...)
    fig = Figure(; size)
    Label(fig[0, 1], "System Structure", fontsize = 24)
    ax = LScene(fig[1, 1])
    plot!(ax, sys; kwargs...)
    update_cam!(ax.scene, Vec3f(-10, -10, 10), Vec3f(0, 0, 0))
    return fig
end

end
