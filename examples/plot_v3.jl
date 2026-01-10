"""
Plot the V3 kite OBJ mesh file.

Usage:
    julia --project=examples examples/plot_v3.jl
"""

using GLMakie
using FileIO
using GeometryBasics
using LinearAlgebra

# Load the OBJ file
obj_path = joinpath(@__DIR__, "..", "data", "v3", "V3_25.obj")
mesh = load(obj_path)

# Rotation matrices
function rot_x(angle)
    c, s = cos(angle), sin(angle)
    [1 0 0; 0 c -s; 0 s c]
end

function rot_z(angle)
    c, s = cos(angle), sin(angle)
    [c -s 0; s c 0; 0 0 1]
end

# Combined rotation: 90deg around x, then -90deg around z
R = rot_z(-π/2) * rot_x(π/2)

# Transform vertices
old_vertices = coordinates(mesh)
new_vertices = [Point3f(R * (Vector(v)) + [0,0,7.3]) for v in old_vertices]

# Create new mesh with transformed vertices
mesh = GeometryBasics.Mesh(new_vertices, faces(mesh))

# Create figure and axis
fig = Figure(size=(1000, 800))
ax = Axis3(fig[1, 1];
    aspect=:data,
    title="V3 Kite Mesh",
    xlabel="x [m]",
    ylabel="y [m]",
    zlabel="z [m]")

# Plot the rotated mesh
mesh!(ax, mesh; color=:lightblue, shading=FastShading)
wireframe!(ax, mesh; color=:gray, linewidth=0.5)

display(fig)
