# Copyright (c) 2025 Bart van de Lint
# SPDX-License-Identifier: MPL-2.0

module SymbolicAWEModelsMakieExt

using Makie
using LinearAlgebra
using StaticArrays
using SymbolicAWEModels

# This is the core implementation. It plots a `SystemStructure` into a
# pre-existing `Axis3`. This makes the plotting logic composable.
function Makie.plot!(ax::Axis3, sys::SystemStructure;
                     point_color = :gray, segment_color = :black,
                     wing_colors = Makie.wong_colors(), vector_scale = 1.0,
                     show_points = true, show_segments = true)

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
    for (i, wing) in enumerate(sys.wings)
        wing_pos = Point3f(wing.pos_w)
        color = wing_colors[mod1(i, length(wing_colors))]
        scatter!(ax, wing_pos, color=color, markersize=20, strokewidth=2, strokecolor=:black, label="Wing $i")

        R = wing.R_b_w
        scale = vector_scale
        origins = [wing_pos, wing_pos, wing_pos]
        directions = [Vec3f(R[:, 1]) * scale, Vec3f(R[:, 2]) * scale, Vec3f(R[:, 3]) * scale]
        
        axis_colors = [:red, :green, :blue]
        arrows3d!(ax, origins, directions, color=axis_colors, label="Wing $i Axes")
    end
    
    return ax
end

# This is the top-level function that gets called when a user types `plot(sys)`.
# It creates a new Figure and Axis3, and then calls the `plot!` method above.
function Makie.plot(sys::SystemStructure; kwargs...)
    fig = Figure(size = (1200, 800))
    
    # Create the Axis3, passing any user-provided keywords for the axis itself
    ax = Axis3(fig[1, 1], title = "System Structure", aspect = :data)

    # Plot into the axis by calling the `plot!` method we defined.
    # All other keywords are passed to the implementation.
    plot!(ax, sys; kwargs...)

    # Add the camera controls legend
    controls_text = """
    Camera Controls:
    - Left-click + drag: Rotate
    - Right-click + drag: Pan
    - Scroll wheel: Zoom
    """
    Label(fig[1, 2], controls_text, tellwidth=false, justification=:left, halign=:left)

    return fig
end

end
