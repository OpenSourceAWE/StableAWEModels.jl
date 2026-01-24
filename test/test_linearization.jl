# SPDX-FileCopyrightText: 2025 Bart van de Lint
# SPDX-License-Identifier: MPL-2.0

using Test, ControlSystemsBase
using SymbolicAWEModels
using LinearAlgebra

# Set up tmpdir if not already done by runtests.jl
if !startswith(get_data_path(), tempdir())
    src_data_path = joinpath(dirname(dirname(pathof(SymbolicAWEModels))), "data", "ram_air_kite")
    tmpdir = mktempdir()
    set_data_path(joinpath(tmpdir, "ram_air_kite"))
    cp(src_data_path, get_data_path(); force=true)
end

set = Settings("system.yaml")
sam = SymbolicAWEModel(set, "ram")

tether_set = Settings("system.yaml")
tether_sam = SymbolicAWEModel(tether_set, "tether")
init!(tether_sam)

simple_set = Settings("system.yaml")
simple_sam = SymbolicAWEModel(simple_set, "simple_ram")
init!(simple_sam)

original_set = Settings("system.yaml")

function reset!(set::Settings)
    for field in fieldnames(Settings)
        setfield!(set, field, getfield(original_set, field))
    end
    return set
end

@testset verbose=true "Linearization" begin
    @testset "Linearize" begin
        old_abs = set.abs_tol
        old_rel = set.rel_tol
        set.abs_tol = 1e-4
        set.rel_tol = 1e-4
        init!(sam)
        init!(simple_sam)

        (; A, B, C, D) = SymbolicAWEModels.linearize!(simple_sam)
        sys = ss(A,B,C,D)
        norm_A = norm(A)
        res = lsim(sys, repeat([-1.0 0.0 -1.0], 2)', [0.0, 0.5])
        println(res.y[:,2])
        @test isapprox(res.y[:,2],
            [-0.0008037289321365251, 0.0004562826732837309, -0.020711457720341487,
                        -0.0017333135190197818], rtol=0.1)

        find_steady_state!(sam)
        (; A, B, C, D) = SymbolicAWEModels.simple_linearize!(sam; tstab=1.0)
        sys = ss(A,B,C,D)
        res = lsim(sys, repeat([-1.0 0.0 -1.0], 2)', [0.0, 0.5])
        println(res.y[:,2])
        @test isapprox(res.y[:,2],
            [0.014234402954620558, -0.0005674058560722778, -0.0186760660540293,
                5.933033873737758], rtol=0.1)

        # test that linearization is state-dependent
        next_step!(simple_sam; dt=1.0)
        (; A, B, C, D) = SymbolicAWEModels.linearize!(simple_sam)
        @test !isapprox(norm(A), norm_A; atol=1e-3)

        set.abs_tol = old_abs
        set.rel_tol = old_rel
    end
end
