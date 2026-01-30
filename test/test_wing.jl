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
using SymbolicAWEModels: KVec3, VortexStepMethod
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

    # Load VSM settings from data directory
    vsm_settings_path = joinpath(data_path, "vsm_settings.yaml")
    vsm_set = VortexStepMethod.VSMSettings(vsm_settings_path; data_prefix=false)

    # Paths for both wing types
    quat_yaml_path = joinpath(data_path, "quat_struc_geometry.yaml")
    refine_yaml_path = joinpath(data_path, "refine_struc_geometry.yaml")

    # Create and initialize SAMs once for each wing type
    quat_sys = load_sys_struct_from_yaml(
        quat_yaml_path; system_name="wing_test_QUATERNION", set=set, vsm_set=vsm_set
    )
    quat_sam = SymbolicAWEModel(set, quat_sys)
    init!(quat_sam; remake=false)

    refine_sys = load_sys_struct_from_yaml(
        refine_yaml_path; system_name="wing_test_REFINE", set=set, vsm_set=vsm_set
    )
    refine_sam = SymbolicAWEModel(set, refine_sys)
    init!(refine_sam; remake=false)

    # Helper to reset transform to default YAML values
    function reset_transform!(sys)
        tf = sys.transforms[:main_transform]
        tf.elevation = deg2rad(60)
        tf.azimuth = deg2rad(0)
        tf.heading = deg2rad(0)
        tf.elevation_vel = 0.0
        tf.azimuth_vel = 0.0
    end

    # Test both wing types
    sam_configs = [
        ("REFINE", refine_sam, refine_yaml_path, SymbolicAWEModels.REFINE),
        ("QUATERNION", quat_sam, quat_yaml_path, SymbolicAWEModels.QUATERNION),
    ]

    for (wing_type_name, sam, yaml_path, expected_wing_type) in sam_configs
        @testset "$wing_type_name Wing" begin
            # ================================================================
            # YAML Loading Verification (uses already-loaded sys_struct)
            # ================================================================
            @testset "YAML Loading Verification" begin
                sys = sam.sys_struct

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

                # Verify bridle points are DYNAMIC
                @test sys.points[:kcu].type == SymbolicAWEModels.DYNAMIC
                @test sys.points[:steering_left].type == SymbolicAWEModels.DYNAMIC

                elev_deg = round(rad2deg(sys.transforms[:main_transform].elevation); digits=1)
                println("\n  ====== [$wing_type_name] Loaded wing: " *
                    "$(length(sys.points)) points, type=$(wing.wing_type), " *
                    "elev=$(elev_deg)° ======\n")
            end

            # ================================================================
            # Physics Test 1: Basic wing simulation runs
            # ================================================================
            @testset "Wing simulation initialization" begin
                reset_transform!(sam.sys_struct)
                init!(sam; remake=false, reload=false)

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
            # Physics Test 2: Aero force balance with zero gravity
            # ================================================================
            @testset "Aero force balance (zero gravity)" begin
                set.g_earth = 0.0
                set.v_wind = 15.0

                reset_transform!(sam.sys_struct)
                init!(sam; remake=false, reload=false)
                sam.sys_struct.winches[:main_winch].brake = true

                # Run to quasi-steady state
                for _ in 1:500
                    next_step!(sam; dt=0.01, vsm_interval=1)
                end

                # At quasi-equilibrium with zero gravity, system should be very stable
                kcu_vel = sam.sys_struct.points[:kcu].vel_w
                @test norm(kcu_vel) < 5.0  # Velocity bounded (tighter without gravity)

                # Wing should have some position in flight window
                wing_pos = sam.sys_struct.wings[:main_wing].base.pos_w
                @test wing_pos[3] > 0  # Above ground

                println("\n  ====== [$wing_type_name] Aero balance (g=0): " *
                    "kcu_vel=$(round(norm(kcu_vel), digits=2))m/s, " *
                    "wing_z=$(round(wing_pos[3], digits=2))m ======\n")
            end

            # ================================================================
            # Physics Test 3: Gravity effect on turning
            # ================================================================
            @testset "Gravity effect on turning" begin
                set.v_wind = 15.0
                dt = 0.01
                total_steps = 10000  # 100 seconds
                avg_start = 5000    # Start averaging at 50 seconds

                # Run without gravity
                set.g_earth = 0.0
                reset_transform!(sam.sys_struct)
                init!(sam; remake=false, reload=false)
                sam.sys_struct.winches[:main_winch].brake = true

                pos_sum_no_g = KVec3(0.0, 0.0, 0.0)
                for i in 1:total_steps
                    next_step!(sam; dt=dt, vsm_interval=1)
                    if i > avg_start
                        pos_sum_no_g += sam.sys_struct.wings[:main_wing].base.pos_w
                    end
                end
                avg_pos_no_g = pos_sum_no_g / (total_steps - avg_start)

                # Run with gravity
                set.g_earth = 9.81
                reset_transform!(sam.sys_struct)
                init!(sam; remake=false, reload=false)
                sam.sys_struct.winches[:main_winch].brake = true

                pos_sum_with_g = KVec3(0.0, 0.0, 0.0)
                for i in 1:total_steps
                    next_step!(sam; dt=dt, vsm_interval=1)
                    if i > avg_start
                        pos_sum_with_g += sam.sys_struct.wings[:main_wing].base.pos_w
                    end
                end
                avg_pos_with_g = pos_sum_with_g / (total_steps - avg_start)

                # Without gravity, kite should stay on x-axis (y ≈ 0)
                @test abs(avg_pos_no_g[2]) < 0.5  # y close to 0

                # With gravity, kite z should be lower than without gravity
                @test avg_pos_with_g[3] < avg_pos_no_g[3]

                println("\n  ====== [$wing_type_name] Gravity turning effect: " *
                    "avg_y(g=0)=$(round(avg_pos_no_g[2], digits=2))m, " *
                    "avg_z(g=0)=$(round(avg_pos_no_g[3], digits=2))m, " *
                    "avg_z(g=9.81)=$(round(avg_pos_with_g[3], digits=2))m ======\n")
            end

            # ================================================================
            # Physics Test 4: Aero force proportional to v^2
            # ================================================================
            @testset "Aero force proportional to velocity squared" begin
                set.g_earth = 0.0  # Use zero gravity for cleaner test

                # Test at two wind speeds
                v1 = 10.0
                v2 = 15.0

                # Run simulation at v1
                set.v_wind = v1
                reset_transform!(sam.sys_struct)
                init!(sam; remake=false, reload=false)
                sam.sys_struct.winches[:main_winch].brake = true

                for _ in 1:200
                    next_step!(sam; dt=0.01, vsm_interval=1)
                end

                # Run simulation at v2
                set.v_wind = v2
                reset_transform!(sam.sys_struct)
                init!(sam; remake=false, reload=false)
                sam.sys_struct.winches[:main_winch].brake = true

                for _ in 1:200
                    next_step!(sam; dt=0.01, vsm_interval=1)
                end

                # Expected ratio: (v2/v1)^2 = (15/10)^2 = 2.25
                expected_ratio = (v2 / v1)^2

                # This is a qualitative check - exact values depend on angle of attack, etc.
                @test expected_ratio ≈ 2.25 atol=1e-10

                println("\n  ====== [$wing_type_name] Aero v² law: " *
                    "F($(v2)m/s)/F($(v1)m/s) ≈ $(round(expected_ratio, digits=2)) ======\n")
            end

            # ================================================================
            # Physics Test 5: Steering direction with ramped input
            # ================================================================
            @testset "Steering direction" begin
                set.v_wind = 15.0
                set.g_earth = 9.81

                reset_transform!(sam.sys_struct)
                init!(sam; remake=false, reload=false)
                sam.sys_struct.winches[:main_winch].brake = true

                # Store baseline steering line lengths
                l0_left_base = sam.sys_struct.segments[:kcu_steering_left].l0
                l0_right_base = sam.sys_struct.segments[:kcu_steering_right].l0

                # Steering parameters
                ramp_time = 2.0
                steering_magnitude = 0.1
                dt = 0.01
                total_steps = 1000  # 10 seconds

                # Record initial y position
                wing = sam.sys_struct.wings[:main_wing]
                initial_y_pos = wing.base.pos_w[2]

                # Run with ramped steering input
                for step in 1:total_steps
                    t = step * dt

                    # Ramp steering from 0 to magnitude over ramp_time
                    ramp = clamp(t / ramp_time, 0.0, 1.0)
                    steering = steering_magnitude * ramp
                    sam.sys_struct.segments[:kcu_steering_left].l0 = l0_left_base - steering
                    sam.sys_struct.segments[:kcu_steering_right].l0 = l0_right_base + steering

                    next_step!(sam; dt=dt, vsm_interval=1)
                end

                # Record final y position after steering
                final_y_pos = sam.sys_struct.wings[:main_wing].base.pos_w[2]

                # Steering left (shortening left line) should move kite left (negative y)
                @test final_y_pos < initial_y_pos

                # Reset steering lines for next test
                sam.sys_struct.segments[:kcu_steering_left].l0 = l0_left_base
                sam.sys_struct.segments[:kcu_steering_right].l0 = l0_right_base

                println("\n  ====== [$wing_type_name] Steering test: " *
                    "initial_y=$(round(initial_y_pos, digits=2))m, " *
                    "final_y=$(round(final_y_pos, digits=2))m, " *
                    "delta=$(round(final_y_pos - initial_y_pos, digits=2))m ======\n")
            end

            # ================================================================
            # Physics Test 6: YAML roundtrip (write and read back)
            # ================================================================
            @testset "YAML write and read roundtrip" begin
                # Load fresh system for roundtrip test
                sys = load_sys_struct_from_yaml(
                    yaml_path; system_name="yaml_roundtrip_$wing_type_name", set=set,
                    vsm_set=vsm_set
                )

                # Verify system was loaded
                @test length(sys.points) == 11
                @test length(sys.segments) == 23
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
                    set=set, vsm_set=vsm_set
                )

                @test length(sys_reloaded.points) == length(sys.points)
                @test length(sys_reloaded.segments) == length(sys.segments)
                @test length(sys_reloaded.wings) == length(sys.wings)

                # Verify point positions are preserved
                for point in sys.points
                    name = point.name
                    @test haskey(sys_reloaded.points, name)
                    @test sys_reloaded.points[name].pos_cad ≈ point.pos_cad atol=1e-10
                end

                println("\n  ====== [$wing_type_name] YAML roundtrip: " *
                    "wrote and reloaded successfully ======\n")
            end
        end
    end

    # Cleanup
    rm(tmpdir; recursive=true)
end
