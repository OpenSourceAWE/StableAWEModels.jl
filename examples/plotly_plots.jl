# SPDX-FileCopyrightText: 2025 SymbolicAWEModels contributors
# SPDX-License-Identifier: MPL-2.0
#
# Makie-based helper utilities for visualising SymbolicAWEModels examples.

using GLMakie
using SymbolicAWEModels
using SymbolicAWEModels: DYNAMIC, QUASI_STATIC, STATIC, WING

const Point3List = Vector{Point3f}

"""
    plot3d_v3(points, segments; title="3D Structure")

Render a static 3D view of the kite structure using Makie.
"""
function plot3d_v3(points, segments; title::AbstractString="3D Structure")
    fig, ax = _setup_axis(title)
    _draw_segments!(ax, points, segments)
    _draw_point_sets!(ax, points)
    display(fig)
    return fig
end

"""
    make_plot3d(points, segments; title="3D Structure")

Create (but do not display) a Makie Figure containing a static 3D plot.
"""
function make_plot3d(points, segments; title::AbstractString="3D Structure")
    fig, ax = _setup_axis(title)
    _draw_segments!(ax, points, segments)
    _draw_point_sets!(ax, points)
    return fig
end

"""
    make_animated_plot3d(states, segments; title="3D Structure", dt=0.02)

Create an interactive 3D plot with a slider that scrubs through a sequence of
system states `states` (each entry should be a vector of `SymbolicAWEModels.Point`).
"""
function make_animated_plot3d(all_states, segments;
                              title::AbstractString="3D Structure",
                              dt::Float64=0.02)
    isempty(all_states) && error("make_animated_plot3d requires at least one state")
    n_frames = length(all_states)

    fig, ax = _setup_axis(title)

    seg_node = Observable(_segment_vertices(all_states[1], segments))
    linesegments!(ax, seg_node; color=:black, linewidth=2)

    dynamic_node = Observable(_points_of_type(all_states[1], DYNAMIC))
    wing_node    = Observable(_points_of_type(all_states[1], WING))
    static_node  = Observable(_points_of_type(all_states[1], STATIC))
    qstatic_node = Observable(_points_of_type(all_states[1], QUASI_STATIC))

    _scatter_if_not_empty!(ax, dynamic_node; color=:dodgerblue, markersize=14)
    _scatter_if_not_empty!(ax, wing_node; color=:firebrick, markersize=16)
    _scatter_if_not_empty!(ax, static_node; color=:grey55, markersize=12)
    _scatter_if_not_empty!(ax, qstatic_node; color=:darkorange, markersize=12)

    slider = Slider(fig[2, 1], range=1:n_frames, startvalue=1,
                    format=i -> "Frame $(i) / $(n_frames)")

    on(slider.value) do idx
        seg_node[] = _segment_vertices(all_states[idx], segments)
        dynamic_node[] = _points_of_type(all_states[idx], DYNAMIC)
        wing_node[]    = _points_of_type(all_states[idx], WING)
        static_node[]  = _points_of_type(all_states[idx], STATIC)
        qstatic_node[] = _points_of_type(all_states[idx], QUASI_STATIC)
    end

    return fig
end

# --------------------------------------------------------------------------- #
# Utility helpers
# --------------------------------------------------------------------------- #

function _setup_axis(title)
    fig = Figure(resolution=(900, 700))
    ax = Axis3(fig[1, 1];
               title,
               xlabel="X (m)",
               ylabel="Y (m)",
               zlabel="Z (m)",
               aspect=:data,
               azimuth=0.8,
               elevation=0.6)
    hidespines!(ax)
    return fig, ax
end

function _draw_segments!(ax, points, segments)
    seg_vertices = _segment_vertices(points, segments)
    isempty(seg_vertices) && return nothing
    linesegments!(ax, seg_vertices; color=:black, linewidth=2)
    return nothing
end

function _draw_point_sets!(ax, points)
    _scatter_if_not_empty!(ax, _points_of_type(points, DYNAMIC);
                           color=:dodgerblue, markersize=14)
    _scatter_if_not_empty!(ax, _points_of_type(points, WING);
                           color=:firebrick, markersize=16)
    _scatter_if_not_empty!(ax, _points_of_type(points, STATIC);
                           color=:grey55, markersize=12)
    _scatter_if_not_empty!(ax, _points_of_type(points, QUASI_STATIC);
                           color=:darkorange, markersize=12)
    return nothing
end

function _points_of_type(points, target_type)
    coords = Point3List()
    for p in points
        p.type == target_type || continue
        push!(coords, Point3f(p.pos_w[1], p.pos_w[2], p.pos_w[3]))
    end
    return coords
end

function _segment_vertices(points, segments)
    verts = Point3List()
    for seg in segments
        i, j = seg.point_idxs
        p1, p2 = points[i], points[j]
        push!(verts, Point3f(p1.pos_w[1], p1.pos_w[2], p1.pos_w[3]))
        push!(verts, Point3f(p2.pos_w[1], p2.pos_w[2], p2.pos_w[3]))
        push!(verts, Point3f(NaN32, NaN32, NaN32)) # break between segments
    end
    return verts
end

function _scatter_if_not_empty!(ax, node::Observable; kwargs...)
    isempty(node[]) && return nothing
    scatter!(ax, node; kwargs...)
    return nothing
end

function _scatter_if_not_empty!(ax, pts::Point3List; kwargs...)
    isempty(pts) && return nothing
    scatter!(ax, pts; kwargs...)
    return nothing
end
