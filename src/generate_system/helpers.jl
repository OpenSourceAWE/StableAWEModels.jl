# Copyright (c) 2025 Bart van de Lint
# SPDX-License-Identifier: MPL-2.0

# Helper functions for symbolic equation generation

"""
    calc_angle_of_attack(va_wing_b)

Calculate the angle of attack [rad] from the apparent wind vector `va_wing_b` in the
body frame.
"""
function calc_angle_of_attack(va_wing_b)
    return atan(va_wing_b[3], va_wing_b[1])
end

"""
    sym_normalize(vec)

Symbolic-safe normalization of a vector. Returns `vec / norm(vec)`.
"""
function sym_normalize(vec)
    return vec / norm(vec)
end

"""
    quaternion_to_rotation_matrix(q)

Convert a quaternion `q` (scalar-first format [w, x, y, z]) to a 3x3 rotation
matrix.
"""
function quaternion_to_rotation_matrix(q)
    w, x, y, z = q[1], q[2], q[3], q[4]

    return [
        1-2*(y*y+z*z) 2*(x*y-z*w) 2*(x*z+y*w)
        2*(x*y+z*w) 1-2*(x*x+z*z) 2*(y*z-x*w)
        2*(x*z-y*w) 2*(y*z+x*w) 1-2*(x*x+y*y)
    ]
end

"""
    rotation_matrix_to_quaternion(R)

Convert a 3x3 rotation matrix `R` to a quaternion (scalar-first format [w, x, y, z]).
This implementation is based on the method that avoids division by zero.
"""
function rotation_matrix_to_quaternion(R)
    tr_ = R[1, 1] + R[2, 2] + R[3, 3]

    if tr_ > 0
        S = sqrt(tr_ + 1.0) * 2
        w = 0.25 * S
        x = (R[3, 2] - R[2, 3]) / S
        y = (R[1, 3] - R[3, 1]) / S
        z = (R[2, 1] - R[1, 2]) / S
    elseif (R[1, 1] > R[2, 2]) && (R[1, 1] > R[3, 3])
        S = sqrt(1.0 + R[1, 1] - R[2, 2] - R[3, 3]) * 2
        w = (R[3, 2] - R[2, 3]) / S
        x = 0.25 * S
        y = (R[1, 2] + R[2, 1]) / S
        z = (R[1, 3] + R[3, 1]) / S
    elseif R[2, 2] > R[3, 3]
        S = sqrt(1.0 + R[2, 2] - R[1, 1] - R[3, 3]) * 2
        w = (R[1, 3] - R[3, 1]) / S
        x = (R[1, 2] + R[2, 1]) / S
        y = 0.25 * S
        z = (R[2, 3] + R[3, 2]) / S
    else
        S = sqrt(1.0 + R[3, 3] - R[1, 1] - R[2, 2]) * 2
        w = (R[2, 1] - R[1, 2]) / S
        x = (R[1, 3] + R[3, 1]) / S
        y = (R[2, 3] + R[3, 2]) / S
        z = 0.25 * S
    end

    return [w, x, y, z]
end

# Component accessors for symbolic registration
rotation_matrix_to_quaternion_w(R) = rotation_matrix_to_quaternion(R)[1]
rotation_matrix_to_quaternion_x(R) = rotation_matrix_to_quaternion(R)[2]
rotation_matrix_to_quaternion_y(R) = rotation_matrix_to_quaternion(R)[3]
rotation_matrix_to_quaternion_z(R) = rotation_matrix_to_quaternion(R)[4]

# Register component functions as symbolic
@register_symbolic rotation_matrix_to_quaternion_w(R::AbstractMatrix)
@register_symbolic rotation_matrix_to_quaternion_x(R::AbstractMatrix)
@register_symbolic rotation_matrix_to_quaternion_y(R::AbstractMatrix)
@register_symbolic rotation_matrix_to_quaternion_z(R::AbstractMatrix)

function calc_wind_factor(am::AtmosphericModel, pos_x, pos_y, pos_z, set::Settings)
    if set.profile_law == 0
        return 1.0
    elseif set.profile_law == 4
        # Linear scaling: 1.0 at kite position, 0.0 at origin
        return sqrt(pos_x^2 + pos_y^2 + pos_z^2) / set.l_tether
    else
        return AtmosphericModels.calc_wind_factor(am, max(1.0, pos_z), set.profile_law)
    end
end
@register_symbolic calc_wind_factor(am::AtmosphericModel, pos_x, pos_y, pos_z,
                                    set::Settings)

"""
    rotate_v_around_k(v, k, θ)

Rotate vector `v` around axis `k` by angle `θ` using Rodrigues' rotation formula.
"""
function rotate_v_around_k(v, k, θ)
    k = sym_normalize(k)
    v_rot = v * cos(θ) + (k × v) * sin(θ) + k * (k ⋅ v) * (1 - cos(θ))
    return v_rot
end

"""
    calc_R_v_w(wing_pos, e_x)

Calculate the rotation matrix from the view frame (`_v`) to the world frame (`_w`).

The view frame is defined with its z-axis pointing from the origin to the wing,
and its x-axis aligned with the wing's x-axis projected onto the view plane.
"""
function calc_R_v_w(wing_pos, e_x)
    z = sym_normalize(wing_pos)
    y = sym_normalize(z × e_x)
    x = y × z
    # Explicit matrix construction for symbolic compatibility
    return [x[1] y[1] z[1]; x[2] y[2] z[2]; x[3] y[3] z[3]]
end

"""
    calc_R_t_w(wing_pos)

Calculate the rotation matrix from the local tether frame (`_t`) to the world
frame (`_w`).

The tether frame is a local spherical coordinate system:
- **z-axis**: Aligned with the tether (radial direction).
- **y-axis**: Azimuthal direction, parallel to the XY plane.
- **x-axis**: Elevation direction, tangent to the sphere (`y × z`).
"""
function calc_R_t_w(wing_pos)
    z = sym_normalize(wing_pos)
    if wing_pos[2] ≈ 0.0 && wing_pos[1] ≈ 0.0
        y = [0, 1, 0]
    else
        y = sym_normalize([-wing_pos[2], wing_pos[1], 0])
    end
    x = y × z
    # Explicit matrix construction for symbolic compatibility
    return [x[1] y[1] z[1]; x[2] y[2] z[2]; x[3] y[3] z[3]]
end

function sym_calc_R_t_w(wing_pos)
    z = sym_normalize(wing_pos)
    y = sym_normalize([-wing_pos[2], wing_pos[1], 0])
    x = y × z
    # Explicit matrix construction for symbolic compatibility
    return [x[1] y[1] z[1]; x[2] y[2] z[2]; x[3] y[3] z[3]]
end

"""
    Base.getindex(x::ModelingToolkit.Symbolics.Arr, idxs::Vector{Int64})

Extend `Base.getindex` to allow indexing a symbolic array with a vector of
integer indices, which is not natively supported by ModelingToolkit.
"""
function Base.getindex(x::ModelingToolkit.Symbolics.Arr, idxs::Vector{Int64})
    Num[Base.getindex(x, idx) for idx in idxs]
end
