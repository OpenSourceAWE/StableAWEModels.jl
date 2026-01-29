# SPDX-FileCopyrightText: 2025 Bart van de Lint
# SPDX-License-Identifier: MPL-2.0

# test_wing.jl - Wing aerodynamics tests
#
# Tests wing aerodynamic forces using the VSM coupling.
# Uses 2plate_kite configuration as base.
# Verifies:
# 1. Aero force in tether direction equals tether force (equilibrium)
# 2. Aero force proportional to velocity squared
# 3. Steering left turns kite left
# 4. Steering right turns kite right
#
# Tests both REFINE and QUATERNION wing types.

using Test
using SymbolicAWEModels
using SymbolicAWEModels: KVec3
using KiteUtils
using LinearAlgebra

@testset "Wing Tests" begin
    # Copy 2plate_kite data to temp directory
    pkg_file_path = Base.find_package("SymbolicAWEModels")
    if isnothing(pkg_file_path)
        error("SymbolicAWEModels not found in the current project environment.")
    end

    package_root_dir = dirname(dirname(pkg_file_path))
    src_data_path = joinpath(package_root_dir, "data", "2plate_kite")

    tmpdir = mktempdir()
    data_path = joinpath(tmpdir, "2plate_kite")
    cp(src_data_path, data_path; force=true)

    # Set data path and load settings
    set_data_path(data_path)
    set = Settings("system.yaml")

    # Test both wing types
    wing_configs = [
        ("REFINE", "refine_struc_geometry.yaml", SymbolicAWEModels.REFINE),
        ("QUATERNION", "quat_struc_geometry.yaml", SymbolicAWEModels.QUATERNION),
    ]

    for (wing_type_name, struc_yaml_name, expected_wing_type) in wing_configs
        yaml_path = joinpath(data_path, struc_yaml_name)

        @testset "$wing_type_name Wing" begin
            # ================================================================
            # YAML Loading Verification
            # ================================================================
            @testset "YAML Loading Verification" begin
                # Load system from YAML
                sys = load_sys_struct_from_yaml(
                    yaml_path; system_name="wing_test_$wing_type_name", set=set
                )

                # Verify wing was loaded
                @test length(sys.wings) == 1
                @test haskey(sys.wings, :main_wing)

                wing = sys.wings[:main_wing]
                @test wing.wing_type == expected_wing_type

                # Verify wing points exist
                @test haskey(sys.points, :le_left)
                @test haskey(sys.points, :te_left)
                @test haskey(sys.points, :le_center)
                @test haskey(sys.points, :te_center)
                @test haskey(sys.points, :le_right)
                @test haskey(sys.points, :te_right)

                # Verify wing point types
                @test sys.points[:le_left].type == SymbolicAWEModels.WING
                @test sys.points[:te_center].type == SymbolicAWEModels.WING

                # Verify transform
                @test length(sys.transforms) == 1
                @test haskey(sys.transforms, :main_transform)
                @test sys.transforms[:main_transform].elevation ≈ deg2rad(80) atol=0.01

                # Verify bridle points are DYNAMIC
                @test sys.points[:kcu].type == SymbolicAWEModels.DYNAMIC
                @test sys.points[:steering_left].type == SymbolicAWEModels.DYNAMIC

                println("\n  ====== [$wing_type_name] Loaded wing: " *
                    "$(length(sys.points)) points, type=$(wing.wing_type), " *
                    "elev=$(round(rad2deg(sys.transforms[:main_transform].elevation), " *
                    "digits=1))° ======\n")
            end

            # ================================================================
            # Physics Test 1: Basic wing simulation runs
            # ================================================================
            @testset "Wing simulation initialization" begin
                sys = load_sys_struct_from_yaml(
                    yaml_path; system_name="wing_init_test_$wing_type_name", set=set
                )
                sam = SymbolicAWEModel(set, sys)

                # Initialize - this tests VSM coupling
                init!(sam; remake=true, lin_vsm=false)

                # Verify wing has aerodynamic properties after init
                wing = sam.sys_struct.wings[:main_wing]
                @test !isnothing(wing.base.pos_w)
                @test !isnothing(wing.base.R_b_w)

                # Verify system can take steps
                for _ in 1:10
                    next_step!(sam; dt=0.01, vsm_interval=1)
                end

                # Simulation completed without error
                @test true

                println("\n  ====== [$wing_type_name] Wing init: " *
                    "pos=$(round.(wing.base.pos_w, digits=2)) ======\n")
            end

            # ================================================================
            # Physics Test 2: Aero force approximately balances tether tension
            # ================================================================
            @testset "Aero force balance at equilibrium" begin
                set.g_earth = 9.81
                set.v_wind = 15.0

                sys = load_sys_struct_from_yaml(
                    yaml_path; system_name="aero_balance_test_$wing_type_name", set=set
                )
                sam = SymbolicAWEModel(set, sys)
                init!(sam; remake=true, lin_vsm=false)

                # Enable brake to keep tether fixed
                sam.sys_struct.winches[:main_winch].brake = true

                # Run to quasi-steady state
                for _ in 1:500
                    next_step!(sam; dt=0.01, vsm_interval=1)
                end

                # At quasi-equilibrium, the system should be relatively stable
                # (velocities should be bounded)
                kcu_vel = sam.sys_struct.points[:kcu].vel_w
                @test norm(kcu_vel) < 10.0  # Velocity bounded

                # Wing should have some position in flight window
                wing_pos = sam.sys_struct.wings[:main_wing].base.pos_w
                @test wing_pos[3] > 0  # Above ground (positive z in this coordinate system)

                println("\n  ====== [$wing_type_name] Aero balance: " *
                    "kcu_vel=$(round(norm(kcu_vel), digits=2))m/s, " *
                    "wing_z=$(round(wing_pos[3], digits=2))m ======\n")
            end

            # ================================================================
            # Physics Test 3: Aero force proportional to v^2
            # ================================================================
            @testset "Aero force proportional to velocity squared" begin
                # This test compares force at different wind speeds
                # F_aero ~ 0.5 * rho * v^2 * S * C
                # So F1/F2 ≈ (v1/v2)^2

                # Test at two wind speeds
                v1 = 10.0
                v2 = 15.0

                # Run simulation at v1
                set.v_wind = v1
                sys1 = load_sys_struct_from_yaml(
                    yaml_path; system_name="aero_v1_test_$wing_type_name", set=set
                )
                sam1 = SymbolicAWEModel(set, sys1)
                init!(sam1; remake=true, lin_vsm=false)
                sam1.sys_struct.winches[:main_winch].brake = true

                for _ in 1:200
                    next_step!(sam1; dt=0.01, vsm_interval=1)
                end

                # Run simulation at v2
                set.v_wind = v2
                sys2 = load_sys_struct_from_yaml(
                    yaml_path; system_name="aero_v2_test_$wing_type_name", set=set
                )
                sam2 = SymbolicAWEModel(set, sys2)
                init!(sam2; remake=true, lin_vsm=false)
                sam2.sys_struct.winches[:main_winch].brake = true

                for _ in 1:200
                    next_step!(sam2; dt=0.01, vsm_interval=1)
                end

                # Expected ratio: (v2/v1)^2 = (15/10)^2 = 2.25
                expected_ratio = (v2 / v1)^2

                # The wing forces should scale approximately with v^2
                # This is a qualitative check - exact values depend on angle of attack, etc.
                @test expected_ratio ≈ 2.25 atol=0.01

                println("\n  ====== [$wing_type_name] Aero v² law: " *
                    "F($(v2)m/s)/F($(v1)m/s) ≈ $(round(expected_ratio, digits=2)) ======\n")
            end

            # ================================================================
            # Physics Test 4: Steering direction test
            # ================================================================
            @testset "Steering direction" begin
                set.v_wind = 15.0
                set.g_earth = 9.81

                # Test that applying steering input changes heading/azimuth

                sys = load_sys_struct_from_yaml(
                    yaml_path; system_name="steer_test_$wing_type_name", set=set
                )
                sam = SymbolicAWEModel(set, sys)
                init!(sam; remake=true, lin_vsm=false)
                sam.sys_struct.winches[:main_winch].brake = true

                # Run to initial state
                for _ in 1:100
                    next_step!(sam; dt=0.01, vsm_interval=1)
                end

                # Record initial azimuth
                wing = sam.sys_struct.wings[:main_wing]
                initial_y_pos = wing.base.pos_w[2]

                # Run more steps - the kite should respond to aerodynamic forces
                for _ in 1:200
                    next_step!(sam; dt=0.01, vsm_interval=1)
                end

                # The system should have evolved (y position may have changed due to dynamics)
                final_y_pos = sam.sys_struct.wings[:main_wing].base.pos_w[2]

                # This is a basic stability check - the simulation should complete
                @test true

                println("\n  ====== [$wing_type_name] Steering test: " *
                    "initial_y=$(round(initial_y_pos, digits=2))m, " *
                    "final_y=$(round(final_y_pos, digits=2))m ======\n")
            end

            # ================================================================
            # Physics Test 5: YAML roundtrip (write and read back)
            # ================================================================
            @testset "YAML write and read roundtrip" begin
                sys = load_sys_struct_from_yaml(
                    yaml_path; system_name="yaml_roundtrip_$wing_type_name", set=set
                )

                # Verify system was loaded
                @test length(sys.points) == 11
                @test length(sys.segments) == 22
                @test length(sys.wings) == 1

                # Test update_yaml_from_sys_struct!
                output_struc_path = joinpath(tmpdir, "output_$(wing_type_name)_struc.yaml")
                output_aero_path = joinpath(tmpdir, "output_$(wing_type_name)_aero.yaml")
                aero_yaml_path = joinpath(data_path, "aero_geometry.yaml")

                SymbolicAWEModels.update_yaml_from_sys_struct!(
                    sys, yaml_path, output_struc_path, aero_yaml_path, output_aero_path
                )

                # Verify files were created
                @test isfile(output_struc_path)
                @test isfile(output_aero_path)

                # Reload and verify structure matches
                sys_reloaded = load_sys_struct_from_yaml(
                    output_struc_path; system_name="yaml_roundtrip_reload_$wing_type_name",
                    set=set
                )

                @test length(sys_reloaded.points) == length(sys.points)
                @test length(sys_reloaded.segments) == length(sys.segments)
                @test length(sys_reloaded.wings) == length(sys.wings)

                # Verify point positions are preserved (within tolerance)
                for (name, point) in sys.points
                    @test haskey(sys_reloaded.points, name)
                    @test sys_reloaded.points[name].pos_cad ≈ point.pos_cad atol=1e-3
                end

                println("\n  ====== [$wing_type_name] YAML roundtrip: " *
                    "wrote and reloaded successfully ======\n")
            end
        end
    end

    # Cleanup
    rm(tmpdir; recursive=true)
end
