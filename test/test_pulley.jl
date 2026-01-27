# SPDX-FileCopyrightText: 2025 Bart van de Lint
# SPDX-License-Identifier: MPL-2.0

# test_pulley.jl - Pulley constraint tests
#
# Tests pulley equilibrium and constraint enforcement with stiff tethers.
# Verifies:
# 1. Length constraint: l_left + l_right = constant
# 2. Equilibrium finding when not initialized at equilibrium
# 3. Analytical geometric equilibrium solution

using Test
using SymbolicAWEModels
using KiteUtils
using LinearAlgebra

# ============================================================================
# YAML Configuration - V-shaped bridle with pulley at apex
# ============================================================================
const PULLEY_TEST_YAML = """
##############################
## Pulley Test System ########
##############################

###########################
## Materials ##############
###########################
materials:
  headers: [name, youngs_modulus, density, damping_per_stiffness]
  data:
    - [dyneema, 55000000000.0, 724, 0.00077]

###########################
## Points #################
###########################
# V-shaped bridle: two attachment points at top, pulley point in middle
points:
  headers: [name, pos_cad, type, wing_idx, transform_idx, extra_mass, body_frame_damping, world_frame_damping, area, drag_coeff]
  data:
    # Two attachment points (like wing tips) - symmetric about x=0
    - [attach_left, [-2.0, 0.0, 10.0], STATIC, nothing, nothing, 0.0, 0.0, 0.0, 0.0, 0.0]
    - [attach_right, [2.0, 0.0, 10.0], STATIC, nothing, nothing, 0.0, 0.0, 0.0, 0.0, 0.0]
    # Pulley point - starts OFF-CENTER at x=0.5 to test equilibrium finding
    - [pulley_point, [0.5, 0.0, 5.0], DYNAMIC, nothing, nothing, 1.0, 50.0, 0.0, 0.0, 0.0]
    # Ground anchor
    - [ground, [0.0, 0.0, 0.0], STATIC, nothing, nothing, 0.0, 0.0, 0.0, 0.0, 0.0]

###########################
## Segments ###############
###########################
segments:
  headers: [name, point_i, point_j, type, l0, diameter_mm, unit_stiffness, unit_damping, compression_frac]
  data:
    # Left bridle leg - use dyneema for stiff tether
    - [left_leg, attach_left, pulley_point, BRIDLE, 6.0, 5.0, dyneema, nothing, 0.01]
    # Right bridle leg
    - [right_leg, attach_right, pulley_point, BRIDLE, 6.0, 5.0, dyneema, nothing, 0.01]
    # Main tether to ground
    - [main_tether, pulley_point, ground, BRIDLE, 5.0, 5.0, dyneema, nothing, 0.01]

###########################
## Pulleys ################
###########################
pulleys:
  headers: [name, segment_i, segment_j, type]
  data:
    - [main_pulley, left_leg, right_leg, DYNAMIC]
"""

# YAML for symmetric starting position
const PULLEY_SYMMETRIC_YAML = """
##############################
## Pulley Test - Symmetric ###
##############################

###########################
## Materials ##############
###########################
materials:
  headers: [name, youngs_modulus, density, damping_per_stiffness]
  data:
    - [dyneema, 55000000000.0, 724, 0.00077]

###########################
## Points #################
###########################
points:
  headers: [name, pos_cad, type, wing_idx, transform_idx, extra_mass, body_frame_damping, world_frame_damping, area, drag_coeff]
  data:
    # Two attachment points symmetric about x=0
    - [attach_left, [-2.0, 0.0, 10.0], STATIC, nothing, nothing, 0.0, 0.0, 0.0, 0.0, 0.0]
    - [attach_right, [2.0, 0.0, 10.0], STATIC, nothing, nothing, 0.0, 0.0, 0.0, 0.0, 0.0]
    # Pulley point - starts at x=0 (symmetric)
    - [pulley_point, [0.0, 0.0, 4.34], DYNAMIC, nothing, nothing, 1.0, 50.0, 0.0, 0.0, 0.0]
    # Ground anchor
    - [ground, [0.0, 0.0, 0.0], STATIC, nothing, nothing, 0.0, 0.0, 0.0, 0.0, 0.0]

###########################
## Segments ###############
###########################
segments:
  headers: [name, point_i, point_j, type, l0, diameter_mm, unit_stiffness, unit_damping, compression_frac]
  data:
    - [left_leg, attach_left, pulley_point, BRIDLE, 6.0, 5.0, dyneema, nothing, 0.01]
    - [right_leg, attach_right, pulley_point, BRIDLE, 6.0, 5.0, dyneema, nothing, 0.01]
    - [main_tether, pulley_point, ground, BRIDLE, 4.34, 5.0, dyneema, nothing, 0.01]

###########################
## Pulleys ################
###########################
pulleys:
  headers: [name, segment_i, segment_j, type]
  data:
    - [main_pulley, left_leg, right_leg, DYNAMIC]
"""

@testset "Pulley Tests" begin
    # Write YAML to temp file
    tmpdir = mktempdir()
    yaml_path = joinpath(tmpdir, "test_pulley_geometry.yaml")
    write(yaml_path, PULLEY_TEST_YAML)

    yaml_symmetric_path = joinpath(tmpdir, "test_pulley_symmetric.yaml")
    write(yaml_symmetric_path, PULLEY_SYMMETRIC_YAML)

    # Create minimal settings file
    settings_yaml = """
system:
  sim_time: 10.0
  segments: 1
  sample_freq: 50

solver:
  solver: "FBDF"
  abs_tol: 0.0001
  rel_tol: 0.0001

kite:
  physical_model: "from_yaml"

tether:
  cd_tether: 0.958
  unit_damping: 350.0
  unit_stiffness: 120000.0
  rho_tether: 724.0
  e_tether: 55000000000.0
  rel_damping: 0.00077
  d_tether: 5.0

winch:
  winch_model: "TorqueControlledMachine"
  max_force: 4000
  v_ro_max: 8.0
  drum_radius: 0.110
  gear_ratio: 1.0
  inertia_total: 0.024
  f_coulomb: 10.0
  c_vf: 5.0

environment:
  v_wind: 0.0
  upwind_dir: -90.0
  h_ref: 6.0
  rho_0: 1.225
"""
    settings_path = joinpath(tmpdir, "settings.yaml")
    write(settings_path, settings_yaml)

    system_yaml = """
system:
  settings: settings.yaml
"""
    system_path = joinpath(tmpdir, "system.yaml")
    write(system_path, system_yaml)

    # Set data path and load settings
    set_data_path(tmpdir)
    set = load_settings("system.yaml")

    # Load system structure from YAML
    sys = load_sys_struct_from_yaml(yaml_path; system_name="pulley_test", set=set)

    # ========================================================================
    # YAML Loading Verification
    # ========================================================================
    @testset "YAML Loading Verification" begin
        # Verify points were loaded correctly
        @test length(sys.points) == 4
        @test haskey(sys.points, :attach_left)
        @test haskey(sys.points, :attach_right)
        @test haskey(sys.points, :pulley_point)
        @test haskey(sys.points, :ground)

        # Verify attachment points are STATIC
        @test sys.points[:attach_left].type == SymbolicAWEModels.STATIC
        @test sys.points[:attach_right].type == SymbolicAWEModels.STATIC

        # Verify pulley point is DYNAMIC and off-center
        @test sys.points[:pulley_point].type == SymbolicAWEModels.DYNAMIC
        @test sys.points[:pulley_point].pos_cad[1] == 0.5  # Off-center at x=0.5
        @test sys.points[:pulley_point].extra_mass == 1.0

        # Verify segments
        @test length(sys.segments) == 3
        @test haskey(sys.segments, :left_leg)
        @test haskey(sys.segments, :right_leg)
        @test haskey(sys.segments, :main_tether)

        # Verify segment rest lengths
        @test sys.segments[:left_leg].l0 == 6.0
        @test sys.segments[:right_leg].l0 == 6.0

        # Verify pulley was loaded
        @test length(sys.pulleys) == 1
        @test haskey(sys.pulleys, :main_pulley)

        pulley = sys.pulleys[:main_pulley]
        @test pulley.type == SymbolicAWEModels.DYNAMIC
    end

    # ========================================================================
    # Physics Test 1: Pulley length constraint verification
    # ========================================================================
    @testset "Pulley length constraint" begin
        set.g_earth = 9.81
        set.v_wind = 0.0

        sys = load_sys_struct_from_yaml(yaml_path; system_name="pulley_len_test", set=set)
        sam = SymbolicAWEModel(set, sys)
        init!(sam; remake=true)

        # Get initial total length of pulley segments
        pulley = sam.sys_struct.pulleys[:main_pulley]
        initial_sum_len = pulley.sum_len

        # This should be l0_left + l0_right = 6.0 + 6.0 = 12.0
        @test initial_sum_len ≈ 12.0 atol=0.01

        # Run simulation
        dt = 0.001
        n_steps = 1000  # 1 second

        sum_len_history = Float64[]

        for _ in 1:n_steps
            next_step!(sam; dt=dt, vsm_interval=0)

            # Calculate current lengths
            left_leg = sam.sys_struct.segments[:left_leg]
            right_leg = sam.sys_struct.segments[:right_leg]

            # Get current segment lengths from point positions
            attach_left = sam.sys_struct.points[:attach_left].pos_w
            attach_right = sam.sys_struct.points[:attach_right].pos_w
            pulley_pos = sam.sys_struct.points[:pulley_point].pos_w

            len_left = norm(pulley_pos - attach_left)
            len_right = norm(pulley_pos - attach_right)
            sum_len = len_left + len_right

            push!(sum_len_history, sum_len)
        end

        # Verify sum of lengths remains approximately constant
        # Note: With stiff tethers (dyneema), actual length may vary due to elasticity
        max_deviation = maximum(abs.(sum_len_history .- initial_sum_len))
        @test max_deviation < 0.5  # Less than 0.5m variation (accounting for elastic effects)

        # The final sum should still be close to initial
        @test sum_len_history[end] ≈ initial_sum_len atol=0.5
    end

    # ========================================================================
    # Physics Test 2: Equilibrium finding from off-center start
    # ========================================================================
    @testset "Equilibrium finding" begin
        set.g_earth = 9.81
        set.v_wind = 0.0

        sys = load_sys_struct_from_yaml(yaml_path; system_name="pulley_eq_test", set=set)
        sam = SymbolicAWEModel(set, sys)
        init!(sam; remake=true)

        # Initial x position is off-center at 0.5
        initial_x = sam.sys_struct.points[:pulley_point].pos_w[1]
        @test abs(initial_x) > 0.4  # Verify we start off-center

        # Run simulation until equilibrium
        dt = 0.001
        n_steps = 5000  # 5 seconds

        for _ in 1:n_steps
            next_step!(sam; dt=dt, vsm_interval=0)
        end

        # At equilibrium, pulley should find symmetric position (x ≈ 0)
        # due to equal tensions from symmetric attachment points
        final_x = sam.sys_struct.points[:pulley_point].pos_w[1]

        # Should be close to x=0 (symmetric equilibrium)
        @test abs(final_x) < 0.3  # Within 30cm of center (allowing for numerical damping)
    end

    # ========================================================================
    # Physics Test 3: Analytical equilibrium position
    # ========================================================================
    @testset "Analytical equilibrium position" begin
        # For a symmetric V-bridle:
        # - Attachment points at (-2, 0, 10) and (2, 0, 10)
        # - Total rope length = 12.0 (each leg l0 = 6.0)
        # - At equilibrium with equal tensions: pulley at x=0
        #
        # Geometry at equilibrium (x=0):
        # - Distance from each attachment to pulley: 6.0
        # - sqrt((2-0)^2 + (10-z)^2) = 6
        # - 4 + (10-z)^2 = 36
        # - (10-z)^2 = 32
        # - z = 10 - sqrt(32) ≈ 4.34

        z_eq_analytical = 10.0 - sqrt(32.0)
        @test z_eq_analytical ≈ 4.343 atol=0.001

        set.g_earth = 9.81
        set.v_wind = 0.0

        # Load symmetric configuration (starting near equilibrium)
        sys = load_sys_struct_from_yaml(yaml_symmetric_path; system_name="pulley_analytic", set=set)
        sam = SymbolicAWEModel(set, sys)
        init!(sam; remake=true)

        # Run to settle at equilibrium
        dt = 0.001
        n_steps = 3000

        for _ in 1:n_steps
            next_step!(sam; dt=dt, vsm_interval=0)
        end

        final_pos = sam.sys_struct.points[:pulley_point].pos_w
        final_x = final_pos[1]
        final_z = final_pos[3]

        # Verify equilibrium position
        @test abs(final_x) < 0.2  # Close to x=0
        @test final_z ≈ z_eq_analytical atol=0.5  # Close to analytical z

        # Verify both legs have equal length at equilibrium
        attach_left = sam.sys_struct.points[:attach_left].pos_w
        attach_right = sam.sys_struct.points[:attach_right].pos_w

        len_left = norm(final_pos - attach_left)
        len_right = norm(final_pos - attach_right)

        @test len_left ≈ len_right atol=0.1  # Equal lengths = equal tensions
    end

    # ========================================================================
    # Physics Test 4: Tension balance verification
    # ========================================================================
    @testset "Tension balance at equilibrium" begin
        set.g_earth = 9.81
        set.v_wind = 0.0

        sys = load_sys_struct_from_yaml(yaml_symmetric_path; system_name="pulley_tension", set=set)
        sam = SymbolicAWEModel(set, sys)
        init!(sam; remake=true)

        # Run to equilibrium
        dt = 0.001
        for _ in 1:3000
            next_step!(sam; dt=dt, vsm_interval=0)
        end

        # At equilibrium, the pulley point should have near-zero velocity
        final_vel = sam.sys_struct.points[:pulley_point].vel_w
        @test norm(final_vel) < 0.1  # Velocity near zero

        # The geometry should be symmetric
        final_pos = sam.sys_struct.points[:pulley_point].pos_w
        attach_left = sam.sys_struct.points[:attach_left].pos_w
        attach_right = sam.sys_struct.points[:attach_right].pos_w

        # Unit vectors from pulley to attachments
        vec_left = attach_left - final_pos
        vec_right = attach_right - final_pos

        len_left = norm(vec_left)
        len_right = norm(vec_right)

        # For equal tensions with symmetric geometry:
        # The horizontal components should cancel
        # T_left * (vec_left/len_left) + T_right * (vec_right/len_right) = [0, 0, -W - T_main]

        # Horizontal components
        unit_left = vec_left / len_left
        unit_right = vec_right / len_right

        # If tensions are equal: T * (unit_left + unit_right) should have x ≈ 0
        combined_horizontal = unit_left[1] + unit_right[1]
        @test abs(combined_horizontal) < 0.1  # Horizontal forces balance
    end

    # Cleanup
    rm(tmpdir; recursive=true)
end
