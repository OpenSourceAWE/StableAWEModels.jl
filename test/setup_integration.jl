# SPDX-FileCopyrightText: 2026 Bart van de Lint
#
# SPDX-License-Identifier: MPL-2.0

# Integration test for end-user setup.
# Run from the repo root: julia test/setup_integration.jl
#
# This runs everything in a single Julia session to avoid repeated
# startup/compilation overhead.

using Test

const REPO_ROOT = dirname(@__DIR__)

# Create a temporary directory simulating a fresh user project
const USER_DIR = mktempdir()
println("  tmpdir: $USER_DIR")
cd(USER_DIR)

# Activate the temp project and dev-install SymbolicAWEModels
using Pkg
Pkg.activate(".")
Pkg.develop(path=REPO_ROOT)

using SymbolicAWEModels

@testset "End-user setup" begin
    @testset "copy_data()" begin
        SymbolicAWEModels.copy_data()
        for d in ["data/2plate_kite", "data/base",
                  "data/saddle_form"]
            @test isdir(d)
        end
    end

    @testset "copy_examples()" begin
        SymbolicAWEModels.copy_examples()
        for f in ["menu.jl", "hanging_mass.jl",
                   "catenary_line.jl", "pulley.jl",
                   "saddle_form.jl", "coupled_2plate_kite.jl",
                   "coupled_2plate_kite_linear_vsm.jl",
                   "coupled_tether_deflection.jl",
                   "coupled_realtime_visualization.jl",
                   "coupled_linearize.jl",
                   "static_load_2plate_kite.jl",
                   "sam_tutorial.jl"]
            @test isfile(joinpath("examples", f))
        end
    end

    @testset "Examples use GLMakie" begin
        for f in ["menu.jl", "hanging_mass.jl",
                  "catenary_line.jl", "pulley.jl",
                  "saddle_form.jl", "coupled_2plate_kite.jl"]
            content = read(joinpath("examples", f), String)
            @test occursin("using GLMakie", content)
        end
    end

    @testset "Run examples" begin
        Pkg.activate(joinpath(USER_DIR, "examples"))
        Pkg.instantiate()
        for f in ["hanging_mass.jl", "catenary_line.jl",
                  "pulley.jl", "saddle_form.jl",
                  "sam_tutorial.jl",
                  "coupled_tether_deflection.jl",
                  "coupled_2plate_kite.jl",
                  "coupled_2plate_kite_linear_vsm.jl",
                  "coupled_linearize.jl",
                  "static_load_2plate_kite.jl"]
            @testset "run $f" begin
                println("  Running $f...")
                include(joinpath(USER_DIR, "examples", f))
            end
        end
    end

    @testset "README pendulum example" begin
        println("  Running README pendulum example...")
        set_data_path("data/base")

        set = Settings("system.yaml")
        set.v_wind = 0.0

        points = [
            Point(:anchor, [0, 0, 0], STATIC),
            Point(:mass, [0, 0, -50], DYNAMIC; extra_mass=1.0),
        ]
        segments = [
            Segment(:spring, set, :anchor, :mass, BRIDLE)
        ]
        transforms = [
            Transform(:tf, deg2rad(-80), 0.0, 0.0;
                base_pos=[0, 0, 50], base_point=:anchor,
                rot_point=:mass)
        ]

        sys = SystemStructure("pendulum", set;
            points, segments, transforms)
        sam = SymbolicAWEModel(set, sys)
        init!(sam)

        for _ in 1:100
            next_step!(sam)
        end
    end

    @testset "README 2plate kite example" begin
        println("  Running README 2plate kite example...")
        using VortexStepMethod

        set_data_path("data/2plate_kite")

        struc_yaml = joinpath(
            get_data_path(), "quat_struc_geometry.yaml")

        set = Settings("system.yaml")
        vsm_set = VortexStepMethod.VSMSettings(
            joinpath(get_data_path(), "vsm_settings.yaml");
            data_prefix=false)

        sys = load_sys_struct_from_yaml(struc_yaml;
            system_name="2plate_kite", set, vsm_set)

        sam = SymbolicAWEModel(set, sys)
        init!(sam)

        l0_left = sam.sys_struct.segments[
            :kcu_steering_left].l0
        l0_right = sam.sys_struct.segments[
            :kcu_steering_right].l0

        for step in 1:600
            t = step * (10.0 / 600)
            ramp = clamp(t / 2.0, 0.0, 1.0)
            sam.sys_struct.segments[
                :kcu_steering_left].l0 =
                l0_left - 0.1 * ramp
            sam.sys_struct.segments[
                :kcu_steering_right].l0 =
                l0_right + 0.1 * ramp
            next_step!(sam; dt=10.0/600, vsm_interval=1)
        end
    end
end
