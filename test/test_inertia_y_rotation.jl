# SPDX-FileCopyrightText: 2026 Uwe Fechner
# SPDX-License-Identifier: LGPL-3.0-only

# test_inertia_y_rotation.jl
#
# `calc_inertia_y_rotation` diagonalizes the XZ block of a wing's inertia
# tensor with a closed-form Y-axis rotation, unique for a body symmetric
# about the XZ-plane. `setup_wing_frame!` uses it instead of the generic
# `principal_frame` (full eigendecomposition + permutation search), which
# picks the axis assignment ambiguously whenever two principal moments are
# close and can flip the body frame ~90° relative to the CAD/aero-panel frame.

using Test
using SymbolicAWEModels
import SymbolicAWEModels: calc_inertia_y_rotation, principal_frame
using LinearAlgebra

@testset "calc_inertia_y_rotation" begin
    @testset "diagonalizes a synthetic XZ-symmetric tensor" begin
        I_tensor = [12.0 0.0 3.5;
                    0.0  30.0 0.0;
                    3.5  0.0  10.0]
        I_diag, Ry = calc_inertia_y_rotation(I_tensor)

        @test isapprox(I_diag[1, 3], 0.0; atol=1e-10)
        @test isapprox(I_diag[3, 1], 0.0; atol=1e-10)
        @test Ry * I_tensor * Ry' ≈ I_diag atol=1e-10

        # Proper rotation about Y only
        @test Ry[1, 2] == 0 && Ry[2, 1] == 0 &&
              Ry[2, 3] == 0 && Ry[3, 2] == 0 && Ry[2, 2] == 1
        @test Ry * Ry' ≈ I(3) atol=1e-10
        @test isapprox(det(Ry), 1.0; atol=1e-12)
    end

    @testset "matches the closed-form axis assignment for the A1-15 wing" begin
        # Raw (pre-diagonalization) inertia tensor of a hybrid wing,
        # as returned by `normalized_inertia`. I_xx (271.8) and I_zz (281.2)
        # are close enough that `principal_frame`'s eigendecomposition +
        # permutation search picks a qualitatively different axis assignment
        # (rotated ~90° from this one) than the closed-form solution.
        I_cad = [271.84319281707144    0.01785249304951421  4.080525782000695;
                   0.01785249304951421 57.3152905862639     -7.327690738960784e-5;
                   4.080525782000695  -7.327690738960784e-5  281.23766987051357]

        I_diag, Ry = calc_inertia_y_rotation(I_cad)
        theta = rad2deg(atan(Ry[1, 3], Ry[1, 1]))
        @test isapprox(theta, 69.5095; atol=1e-3)

        # The generic permutation-search rotation is ~90° away from the
        # closed-form one on this near-degenerate tensor.
        _, R_gen = principal_frame(I_cad)
        theta_gen = rad2deg(atan(R_gen[1, 3], R_gen[1, 1]))
        angle_diff = abs(mod(theta - theta_gen + 180, 360) - 180)
        @test isapprox(angle_diff, 90.0; atol=1.0)
    end
end
nothing
