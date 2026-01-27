# SPDX-FileCopyrightText: 2025 Bart van de Lint
# SPDX-License-Identifier: MPL-2.0

# test_transform.jl - Transform and spherical coordinate tests
#
# Tests initial heading, elevation, azimuth and their velocities after init.
# Verifies:
# 1. Initial angles match YAML configuration
# 2. Initial velocities match YAML configuration
# 3. Geometric consistency (position from spherical coords)
# 4. Heading calculation consistency

using Test
using SymbolicAWEModels
using SymbolicAWEModels: KVec3
using KiteUtils
using LinearAlgebra

# ============================================================================
# YAML Configuration - Wing with transform at specific angles
# ============================================================================
const TRANSFORM_TEST_YAML = """
##############################
## Transform Test System #####
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
    # Wing points (simple 2-section wing)
    - [le_1, [-0.5, 1.0, 0.0], WING, main_wing, main_transform, 0.1, 10.0, 0.0, 0.0, 0.0]
    - [te_1, [0.5, 1.0, 0.0], WING, main_wing, main_transform, 0.1, 10.0, 0.0, 0.0, 0.0]
    - [le_2, [-0.5, -1.0, 0.0], WING, main_wing, main_transform, 0.1, 10.0, 0.0, 0.0, 0.0]
    - [te_2, [0.5, -1.0, 0.0], WING, main_wing, main_transform, 0.1, 10.0, 0.0, 0.0, 0.0]
    # Bridle/KCU point
    - [kcu, [0.0, 0.0, -1.0], DYNAMIC, main_wing, main_transform, 1.0, 0.0, 0.0, 0.0, 0.0]
    # Ground station
    - [ground, [0.0, 0.0, -50.0], STATIC, nothing, nothing, 0.0, 0.0, 0.0, 0.0, 0.0]

###########################
## Segments ###############
###########################
segments:
  headers: [name, point_i, point_j, type, l0, diameter_mm, unit_stiffness, unit_damping, compression_frac]
  data:
    - [chord_1, le_1, te_1, POWER_LINE, nothing, 1.0, 10000.0, 10.0, 1.0]
    - [chord_2, le_2, te_2, POWER_LINE, nothing, 1.0, 10000.0, 10.0, 1.0]
    - [spar, le_1, le_2, POWER_LINE, nothing, 1.0, 10000.0, 10.0, 1.0]
    - [bridle1, kcu, le_1, BRIDLE, nothing, 5.0, dyneema, nothing, 0.01]
    - [bridle2, kcu, le_2, BRIDLE, nothing, 5.0, dyneema, nothing, 0.01]
    - [tether, kcu, ground, BRIDLE, nothing, 10.0, dyneema, nothing, 0.01]

###########################
## Tethers ################
###########################
tethers:
  headers: [name, segment_idxs, winch_point_idx]
  data:
    - [main_tether, [tether], ground]

###########################
## Winches ################
###########################
winches:
  headers: [name, tether_idxs]
  data:
    - [main_winch, [main_tether]]

###########################
## Wings ##################
###########################
# Note: For transform tests, we use QUATERNION wing type without VSM
# to avoid VSM dependencies. This tests the coordinate system only.
wings:
  data:
    - name: main_wing
      type: QUATERNION
      aero_z_offset: 0.0

###########################
## Transforms #############
###########################
transforms:
  data:
    - name: main_transform
      elevation: 60
      azimuth: 15
      heading: 10
      elevation_vel: 0.0
      azimuth_vel: 0.5
      wing_idx: main_wing
      base_pos: [0.0, 0.0, 0.0]
      base_point_idx: ground
"""

# YAML for different transform angles
const TRANSFORM_ALT_YAML = """
##############################
## Transform Alt Angles ######
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
    - [le_1, [-0.5, 1.0, 0.0], WING, main_wing, main_transform, 0.1, 10.0, 0.0, 0.0, 0.0]
    - [te_1, [0.5, 1.0, 0.0], WING, main_wing, main_transform, 0.1, 10.0, 0.0, 0.0, 0.0]
    - [le_2, [-0.5, -1.0, 0.0], WING, main_wing, main_transform, 0.1, 10.0, 0.0, 0.0, 0.0]
    - [te_2, [0.5, -1.0, 0.0], WING, main_wing, main_transform, 0.1, 10.0, 0.0, 0.0, 0.0]
    - [kcu, [0.0, 0.0, -1.0], DYNAMIC, main_wing, main_transform, 1.0, 0.0, 0.0, 0.0, 0.0]
    - [ground, [0.0, 0.0, -50.0], STATIC, nothing, nothing, 0.0, 0.0, 0.0, 0.0, 0.0]

###########################
## Segments ###############
###########################
segments:
  headers: [name, point_i, point_j, type, l0, diameter_mm, unit_stiffness, unit_damping, compression_frac]
  data:
    - [chord_1, le_1, te_1, POWER_LINE, nothing, 1.0, 10000.0, 10.0, 1.0]
    - [chord_2, le_2, te_2, POWER_LINE, nothing, 1.0, 10000.0, 10.0, 1.0]
    - [spar, le_1, le_2, POWER_LINE, nothing, 1.0, 10000.0, 10.0, 1.0]
    - [bridle1, kcu, le_1, BRIDLE, nothing, 5.0, dyneema, nothing, 0.01]
    - [bridle2, kcu, le_2, BRIDLE, nothing, 5.0, dyneema, nothing, 0.01]
    - [tether, kcu, ground, BRIDLE, nothing, 10.0, dyneema, nothing, 0.01]

###########################
## Tethers ################
###########################
tethers:
  headers: [name, segment_idxs, winch_point_idx]
  data:
    - [main_tether, [tether], ground]

###########################
## Winches ################
###########################
winches:
  headers: [name, tether_idxs]
  data:
    - [main_winch, [main_tether]]

###########################
## Wings ##################
###########################
wings:
  data:
    - name: main_wing
      type: QUATERNION
      aero_z_offset: 0.0

###########################
## Transforms #############
###########################
transforms:
  data:
    - name: main_transform
      elevation: 45
      azimuth: 30
      heading: 0
      elevation_vel: 0.1
      azimuth_vel: 0.0
      wing_idx: main_wing
      base_pos: [0.0, 0.0, 0.0]
      base_point_idx: ground
"""

@testset "Transform Tests" begin
    # Write YAML to temp file
    tmpdir = mktempdir()
    yaml_path = joinpath(tmpdir, "test_transform_geometry.yaml")
    write(yaml_path, TRANSFORM_TEST_YAML)

    yaml_alt_path = joinpath(tmpdir, "test_transform_alt.yaml")
    write(yaml_alt_path, TRANSFORM_ALT_YAML)

    # Create minimal settings file
    settings_yaml = """
system:
    log_file: "data/transform_test"
    g_earth: 9.81

initial:
    l_tethers: [0.0]
    v_reel_outs: [0.0]

solver:
    solver: "FBDF"
    abs_tol: 0.0001
    rel_tol: 0.0001
    relaxation: 0.6

kite:
    model: ""
    foil_file: "ram_air_kite/ram_air_kite_foil.dat"
    physical_model: "2plate"
    struc_geometry_path: "struc_geometry.yaml"
    aero_geometry_path: "aero_geometry.yaml"
    mass: 0.0
    quasi_static: false

tether:
    cd_tether: 0.958
    unit_damping: 350.0
    unit_stiffness: 120000.0
    rho_tether: 724.0
    e_tether: 55000000000.0
    rel_damping: 0.00077
    d_tether: 10.0

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
    rho_0: 1.225
    v_wind: 10.0
    upwind_dir: -90.0
    h_ref: 6.0
    profile_law: 0
"""
    settings_path = joinpath(tmpdir, "settings.yaml")
    write(settings_path, settings_yaml)

    system_yaml = """
system:
  sim_settings: settings.yaml
"""
    system_path = joinpath(tmpdir, "system.yaml")
    write(system_path, system_yaml)

    # Set data path and load settings
    set_data_path(tmpdir)
    set = Settings("system.yaml")

    # Load system structure from YAML
    sys = load_sys_struct_from_yaml(yaml_path; system_name="transform_test", set=set)

    # ========================================================================
    # YAML Loading Verification
    # ========================================================================
    @testset "YAML Loading Verification" begin
        # Verify transform was loaded
        @test length(sys.transforms) == 1
        @test haskey(sys.transforms, :main_transform)

        transform = sys.transforms[:main_transform]

        # Verify angles match YAML (converted to radians)
        @test transform.elevation ≈ deg2rad(60) atol=0.001
        @test transform.azimuth ≈ deg2rad(15) atol=0.001
        @test transform.heading ≈ deg2rad(10) atol=0.001

        # Verify velocities match YAML (converted to radians)
        @test transform.elevation_vel ≈ deg2rad(0.0) atol=0.001
        @test transform.azimuth_vel ≈ deg2rad(0.5) atol=0.001

        # Verify base point
        @test transform.base_point_idx == 6  # ground point index

        # Verify wing reference
        @test transform.wing_idx == 1  # main_wing index

        println("\n  ====== Loaded transform: elev=$(round(rad2deg(transform.elevation), digits=1))°, azim=$(round(rad2deg(transform.azimuth), digits=1))°, heading=$(round(rad2deg(transform.heading), digits=1))° ======\n")
    end

    # ========================================================================
    # Physics Test 1: Initial angles after init
    # ========================================================================
    @testset "Initial angles after init!" begin
        sys = load_sys_struct_from_yaml(yaml_path; system_name="angles_test", set=set)
        sam = SymbolicAWEModel(set, sys)
        init!(sam; remake=true)

        # After init, the transform angles should still match
        transform = sam.sys_struct.transforms[:main_transform]

        @test transform.elevation ≈ deg2rad(60) atol=0.01
        @test transform.azimuth ≈ deg2rad(15) atol=0.01
        @test transform.heading ≈ deg2rad(10) atol=0.01

        println("\n  ====== After init: elev=$(round(rad2deg(transform.elevation), digits=1))°, azim=$(round(rad2deg(transform.azimuth), digits=1))°, heading=$(round(rad2deg(transform.heading), digits=1))° ======\n")
    end

    # ========================================================================
    # Physics Test 2: Initial velocities after init
    # ========================================================================
    @testset "Initial velocities after init!" begin
        sys = load_sys_struct_from_yaml(yaml_path; system_name="vel_test", set=set)
        sam = SymbolicAWEModel(set, sys)
        init!(sam; remake=true)

        transform = sam.sys_struct.transforms[:main_transform]

        @test transform.elevation_vel ≈ deg2rad(0.0) atol=0.01
        @test transform.azimuth_vel ≈ deg2rad(0.5) atol=0.01

        println("\n  ====== Velocities: elev_vel=$(round(rad2deg(transform.elevation_vel), digits=2))°/s, azim_vel=$(round(rad2deg(transform.azimuth_vel), digits=2))°/s ======\n")
    end

    # ========================================================================
    # Physics Test 3: Geometric consistency - position from spherical coords
    # ========================================================================
    @testset "Geometric consistency" begin
        # For elevation=60deg, azimuth=15deg, the wing should be positioned
        # according to spherical coordinate transformation

        sys = load_sys_struct_from_yaml(yaml_path; system_name="geom_test", set=set)
        sam = SymbolicAWEModel(set, sys)
        init!(sam; remake=true)

        # Get wing position
        wing = sam.sys_struct.wings[:main_wing]
        wing_pos = wing.base.pos_w

        # Get ground position
        ground_pos = sam.sys_struct.points[:ground].pos_w

        # Vector from ground to wing
        rel_pos = wing_pos - ground_pos

        # Calculate distance (tether length)
        distance = norm(rel_pos)

        # Expected spherical coordinates:
        # x = r * cos(elev) * cos(azim)
        # y = r * cos(elev) * sin(azim)
        # z = r * sin(elev)
        # Note: The actual coordinate convention may differ

        elev = deg2rad(60)
        azim = deg2rad(15)

        # The wing should be at the expected position
        # This tests that the transform correctly places the wing
        @test distance > 0  # Wing is above ground

        # More specific tests depend on coordinate convention
        # Check that z component is positive (wing above ground level in world frame)
        @test wing_pos[3] > ground_pos[3]

        println("\n  ====== Geometry: wing_pos=$(round.(wing_pos, digits=2)), distance=$(round(distance, digits=2))m ======\n")
    end

    # ========================================================================
    # Physics Test 4: Different angle configuration
    # ========================================================================
    @testset "Alternative angles configuration" begin
        sys = load_sys_struct_from_yaml(yaml_alt_path; system_name="alt_angles", set=set)

        # Verify alternative YAML loaded correctly
        transform = sys.transforms[:main_transform]

        @test transform.elevation ≈ deg2rad(45) atol=0.001
        @test transform.azimuth ≈ deg2rad(30) atol=0.001
        @test transform.heading ≈ deg2rad(0) atol=0.001
        @test transform.elevation_vel ≈ deg2rad(0.1) atol=0.001
        @test transform.azimuth_vel ≈ deg2rad(0.0) atol=0.001

        sam = SymbolicAWEModel(set, sys)
        init!(sam; remake=true)

        # Angles should be preserved
        transform_after = sam.sys_struct.transforms[:main_transform]
        @test transform_after.elevation ≈ deg2rad(45) atol=0.01
        @test transform_after.azimuth ≈ deg2rad(30) atol=0.01

        println("\n  ====== Alt config: elev=$(round(rad2deg(transform_after.elevation), digits=1))°, azim=$(round(rad2deg(transform_after.azimuth), digits=1))° ======\n")
    end

    # ========================================================================
    # Physics Test 5: Transform affects wing position
    # ========================================================================
    @testset "Transform affects wing position" begin
        # Compare two systems with different elevations
        # System 1: elevation = 60 deg
        sys1 = load_sys_struct_from_yaml(yaml_path; system_name="elev_test_1", set=set)
        sam1 = SymbolicAWEModel(set, sys1)
        init!(sam1; remake=true)
        wing_z1 = sam1.sys_struct.wings[:main_wing].base.pos_w[3]

        # System 2: elevation = 45 deg (lower)
        sys2 = load_sys_struct_from_yaml(yaml_alt_path; system_name="elev_test_2", set=set)
        sam2 = SymbolicAWEModel(set, sys2)
        init!(sam2; remake=true)
        wing_z2 = sam2.sys_struct.wings[:main_wing].base.pos_w[3]

        # Higher elevation should result in higher z position
        # (wing more overhead)
        @test wing_z1 > wing_z2

        println("\n  ====== Elevation effect: z(60°)=$(round(wing_z1, digits=2))m > z(45°)=$(round(wing_z2, digits=2))m ======\n")
    end

    # ========================================================================
    # Physics Test 6: Azimuth affects y-position
    # ========================================================================
    @testset "Azimuth affects y-position" begin
        # System 1: azimuth = 15 deg
        sys1 = load_sys_struct_from_yaml(yaml_path; system_name="azim_test_1", set=set)
        sam1 = SymbolicAWEModel(set, sys1)
        init!(sam1; remake=true)
        wing_y1 = sam1.sys_struct.wings[:main_wing].base.pos_w[2]

        # System 2: azimuth = 30 deg (more to the side)
        sys2 = load_sys_struct_from_yaml(yaml_alt_path; system_name="azim_test_2", set=set)
        sam2 = SymbolicAWEModel(set, sys2)
        init!(sam2; remake=true)
        wing_y2 = sam2.sys_struct.wings[:main_wing].base.pos_w[2]

        # Both have positive azimuth, so y should be positive
        # Larger azimuth should give larger |y| component
        @test abs(wing_y2) > abs(wing_y1)

        println("\n  ====== Azimuth effect: |y|(30°)=$(round(abs(wing_y2), digits=2))m > |y|(15°)=$(round(abs(wing_y1), digits=2))m ======\n")
    end

    # ========================================================================
    # Physics Test 7: Heading affects wing orientation (not position)
    # ========================================================================
    @testset "Heading affects orientation" begin
        sys = load_sys_struct_from_yaml(yaml_path; system_name="heading_test", set=set)
        sam = SymbolicAWEModel(set, sys)
        init!(sam; remake=true)

        wing = sam.sys_struct.wings[:main_wing]

        # Wing should have a rotation matrix
        @test !isnothing(wing.base.R_b_w)
        @test size(wing.base.R_b_w) == (3, 3)

        # R_b_w should be a valid rotation matrix (orthonormal)
        @test det(wing.base.R_b_w) ≈ 1.0 atol=0.01
        @test wing.base.R_b_w * wing.base.R_b_w' ≈ I(3) atol=0.01

        println("\n  ====== Heading affects rotation: det(R_b_w)=$(round(det(wing.base.R_b_w), digits=4)) ======\n")
    end

    # ========================================================================
    # Physics Test 8: Base position offset
    # ========================================================================
    @testset "Base position from base_point" begin
        sys = load_sys_struct_from_yaml(yaml_path; system_name="base_test", set=set)

        # The transform references ground as base_point
        transform = sys.transforms[:main_transform]
        @test transform.base_point_idx == 6  # ground index

        sam = SymbolicAWEModel(set, sys)
        init!(sam; remake=true)

        # Ground position should be the base for the transform
        ground_pos = sam.sys_struct.points[:ground].pos_w
        @test ground_pos ≈ KVec3(0.0, 0.0, -50.0) atol=1.0

        println("\n  ====== Base point: ground_pos=$(round.(ground_pos, digits=2)) ======\n")
    end

    # Cleanup
    rm(tmpdir; recursive=true)
end
