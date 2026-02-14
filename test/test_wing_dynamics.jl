# SPDX-FileCopyrightText: 2025 Bart van de Lint
# SPDX-License-Identifier: MPL-2.0

# test_wing_dynamics.jl - Wing rigid body dynamics tests
#
# Verifies QUATERNION wing dynamics with AERO_NONE.
# A wing-only model (no tethers, no bridle) in free fall under
# gravity should have translational acceleration = g.

using Test
using SymbolicAWEModels
using SymbolicAWEModels: KVec3, VortexStepMethod
using KiteUtils
using LinearAlgebra

# ==================== YAML DEFINITIONS ==================== #

# Minimal YAML: just wing + ground point, no segments/tethers/winches
const WING_FREEFALL_YAML = """
points:
  headers: [name, pos_cad, type, wing_idx, transform_idx,
            extra_mass, body_frame_damping,
            world_frame_damping, area, drag_coeff]
  data:
    - [le_left,   [-0.5, 1.0, 2.0], WING, main_wing,
       ~, 0.5, 0.0, 0.0, 0.0, 0.0]
    - [te_left,   [0.5,  1.0, 2.2], WING, main_wing,
       ~, 0.5, 0.0, 0.0, 0.0, 0.0]
    - [le_center, [-0.5, 0.0, 2.5], WING, main_wing,
       ~, 0.5, 0.0, 0.0, 0.0, 0.0]
    - [te_center, [0.5,  0.0, 2.7], WING, main_wing,
       ~, 0.5, 0.0, 0.0, 0.0, 0.0]
    - [le_right,  [-0.5,-1.0, 2.0], WING, main_wing,
       ~, 0.5, 0.0, 0.0, 0.0, 0.0]
    - [te_right,  [0.5, -1.0, 2.2], WING, main_wing,
       ~, 0.5, 0.0, 0.0, 0.0, 0.0]
    - [ground,    [0.0, 0.0, 0.0],  STATIC, ~,
       ~, 0.0, 0.0, 0.0, 0.0, 0.0]

groups:
  headers: [name, point_idxs, type, moment_frac, damping]
  data:
    - [left,   [le_left, te_left],     DYNAMIC, 0.25, 0.0]
    - [center, [le_center, te_center], DYNAMIC, 0.25, 0.0]
    - [right,  [le_right, te_right],   DYNAMIC, 0.25, 0.0]

wings:
  data:
    - name: main_wing
      type: QUATERNION
      aero_mode: AERO_NONE
      transform_idx: 0
      point_idxs: [le_left, te_left, le_center, te_center,
                    le_right, te_right]
      groups: [left, center, right]
      y_damping: 0.0
      aero_z_offset: 0.0
"""

# Wing + segment pendulum: wing hangs from static anchor via segment.
# attach point is 1.0 above wing CoM; anchor is 1.0 to the side.
# Under gravity + angular damping, settles to vertical alignment:
#   anchor (top) → attach → wing CoM (bottom)
const WING_PENDULUM_YAML = """
points:
  headers: [name, pos_cad, type, wing_idx, transform_idx,
            extra_mass, body_frame_damping,
            world_frame_damping, area, drag_coeff]
  data:
    - [le_left,   [-0.5,  1.0, 0.0], WING, wing,
       ~, 0.5, 0.0, 0.0, 0.0, 0.0]
    - [te_left,   [ 0.5,  1.0, 0.0], WING, wing,
       ~, 0.5, 0.0, 0.0, 0.0, 0.0]
    - [le_center, [-0.5,  0.0, 0.0], WING, wing,
       ~, 0.5, 0.0, 0.0, 0.0, 0.0]
    - [te_center, [ 0.5,  0.0, 0.0], WING, wing,
       ~, 0.5, 0.0, 0.0, 0.0, 0.0]
    - [le_right,  [-0.5, -1.0, 0.0], WING, wing,
       ~, 0.5, 0.0, 0.0, 0.0, 0.0]
    - [te_right,  [ 0.5, -1.0, 0.0], WING, wing,
       ~, 0.5, 0.0, 0.0, 0.0, 0.0]
    - [attach,    [ 0.0,  0.0, 1.0], WING, wing,
       ~, 0.5, 0.0, 0.0, 0.0, 0.0]
    - [anchor,    [ 1.0,  0.0, 0.0], STATIC, ~,
       ~, 0.0, 0.0, 0.0, 0.0, 0.0]

segments:
  headers: [name, point_i, point_j, type, l0,
            diameter_mm, unit_stiffness, unit_damping,
            compression_frac]
  data:
    - [tether, attach, anchor, POWER_LINE, ~,
       1.0, 5000.0, 0.0, 0.01]

groups:
  headers: [name, point_idxs, type, moment_frac, damping]
  data:
    - [left,   [le_left, te_left],     DYNAMIC, 0.25, 0.0]
    - [center, [le_center, te_center], DYNAMIC, 0.25, 0.0]
    - [right,  [le_right, te_right],   DYNAMIC, 0.25, 0.0]

wings:
  data:
    - name: wing
      type: QUATERNION
      aero_mode: AERO_NONE
      transform_idx: 0
      groups: [left, center, right]
      y_damping: 0.0
      angular_damping: 10.0
      aero_z_offset: 0.0
"""

const SETTINGS_YAML = """
system:
    log_file: "data/wing_test"
    g_earth: 9.81

initial:
    l_tethers: [0.0]
    v_reel_outs: [0.0]

solver:
    solver: "FBDF"
    abs_tol: 0.0001
    rel_tol: 0.0001

kite:
    model: ""
    foil_file: "ram_air_kite/ram_air_kite_foil.dat"
    physical_model: "wing_test"
    mass: 0.0
    quasi_static: false

tether:
    cd_tether: 0.958
    unit_damping: 0.0
    unit_stiffness: 0.0
    rho_tether: 724.0
    e_tether: 5.5e10

winch:
    winch_model: "TorqueControlledMachine"
    drum_radius: 0.110
    gear_ratio: 1.0
    inertia_total: 0.024
    f_coulomb: 122.0
    c_vf: 30.6

environment:
    rho_0: 1.225
    v_wind: 0.0
    upwind_dir: -90.0
    profile_law: 0
"""

@testset "Wing Dynamics" begin
    # Copy 2plate_kite data for VSM settings and airfoil files
    pkg_file_path = Base.find_package("SymbolicAWEModels")
    isnothing(pkg_file_path) && error("SymbolicAWEModels not found")
    package_root_dir = dirname(dirname(pkg_file_path))
    src_data_path = joinpath(package_root_dir, "data", "2plate_kite")

    tmpdir = mktempdir()
    data_path = joinpath(tmpdir, "2plate_kite")
    cp(src_data_path, data_path; force=true)

    # Write settings
    settings_path = joinpath(data_path, "settings.yaml")
    write(settings_path, SETTINGS_YAML)
    system_path = joinpath(data_path, "system.yaml")
    write(system_path, "system:\n  sim_settings: settings.yaml\n")

    # Write wing-only YAML
    yaml_path = joinpath(data_path, "wing_freefall.yaml")
    write(yaml_path, WING_FREEFALL_YAML)

    set_data_path(data_path)
    set = Settings("system.yaml")

    vsm_set = VortexStepMethod.VSMSettings(
        joinpath(data_path, "vsm_settings.yaml");
        data_prefix=false
    )

    sys = load_sys_struct_from_yaml(
        yaml_path;
        system_name="wing_freefall",
        set, vsm_set,
        aero_mode=AERO_NONE
    )

    @testset "Model setup" begin
        @test length(sys.wings) == 1
        wing = sys.wings[:main_wing]
        @test wing.wing_type == SymbolicAWEModels.QUATERNION
        @test wing.aero_mode == AERO_NONE
        @test wing.mass ≈ 3.0  # 6 points × 0.5 kg
        @test length(sys.segments) == 0
        @test length(sys.tethers) == 0
        @test length(sys.winches) == 0
        println("  Wing mass: $(wing.mass) kg")
        println("  Inertia: $(wing.inertia_principal)")
    end

    sam = SymbolicAWEModel(set, sys)
    init!(sam; remake=true, prn=true)

    @testset "Free fall acceleration = g" begin
        wing = sam.sys_struct.wings[:main_wing]

        # Let solver settle for a few steps
        dt = 0.01
        for _ in 1:5
            next_step!(sam; dt, vsm_interval=0,
                error_on_unstable=false)
        end

        # Measure velocity before
        vel_before = copy(wing.vel_w)
        t_before = sam.integrator.t

        # Run more steps
        n_steps = 10
        for _ in 1:n_steps
            next_step!(sam; dt, vsm_interval=0,
                error_on_unstable=false)
        end

        vel_after = copy(wing.vel_w)
        t_after = sam.integrator.t
        elapsed = t_after - t_before

        # Compute measured acceleration
        measured_acc = (vel_after - vel_before) / elapsed

        println("  vel_before: $(round.(vel_before; digits=4))")
        println("  vel_after:  $(round.(vel_after; digits=4))")
        println("  elapsed:    $(round(elapsed; digits=4)) s")
        println("  measured_acc: $(round.(measured_acc; digits=4))")
        println("  expected:     [0, 0, -9.81]")

        # z-acceleration should be -g
        @test measured_acc[3] ≈ -9.81 atol=0.1

        # x,y acceleration should be near zero
        @test abs(measured_acc[1]) < 0.1
        @test abs(measured_acc[2]) < 0.1

        # Magnitude check
        @test norm(measured_acc) ≈ 9.81 atol=0.15
    end

    # ==================== Pendulum settling ==================== #
    # Wing + segment: attach point connected to static anchor.
    # Wing should rotate and translate until vertical alignment.

    pendulum_yaml_path = joinpath(data_path, "wing_pendulum.yaml")
    write(pendulum_yaml_path, WING_PENDULUM_YAML)

    pend_sys = load_sys_struct_from_yaml(
        pendulum_yaml_path;
        system_name="wing_pendulum",
        set, vsm_set,
        aero_mode=AERO_NONE
    )

    @testset "Pendulum setup" begin
        @test length(pend_sys.wings) == 1
        @test length(pend_sys.segments) == 1
        @test pend_sys.wings[1].aero_mode == AERO_NONE
        # 7 WING points × 0.5 kg
        @test pend_sys.wings[1].mass ≈ 3.5
    end

    pend_sam = SymbolicAWEModel(set, pend_sys)
    init!(pend_sam; remake=true, prn=true)

    @testset "Pendulum vertical alignment" begin
        wing = pend_sam.sys_struct.wings[1]
        anchor = pend_sam.sys_struct.points[:anchor]
        attach = pend_sam.sys_struct.points[:attach]

        # Simulate until settled (near-critical angular_damping=10)
        dt = 0.05
        n_steps = 2000  # 100 seconds
        for i in 1:n_steps
            next_step!(pend_sam; dt, vsm_interval=0,
                error_on_unstable=false)
        end

        # Get final positions
        anchor_pos = copy(anchor.pos_w)
        attach_pos = copy(attach.pos_w)
        com_pos = copy(wing.pos_w)

        println("  anchor: $(round.(anchor_pos; digits=3))")
        println("  attach: $(round.(attach_pos; digits=3))")
        println("  CoM:    $(round.(com_pos; digits=3))")
        println("  ω_b:    $(round.(wing.ω_b; digits=4))")

        tol = 0.3

        # x-alignment: all at anchor.x
        @test attach_pos[1] ≈ anchor_pos[1] atol=tol
        @test com_pos[1] ≈ anchor_pos[1] atol=tol

        # y-alignment: all near 0
        @test abs(attach_pos[2]) < tol
        @test abs(com_pos[2]) < tol

        # z-ordering: anchor above attach above CoM
        # (gravity is -z, so higher z = higher up)
        @test anchor_pos[3] > attach_pos[3]
        @test attach_pos[3] > com_pos[3]

        # Angular velocity should be near zero (settled)
        @test norm(wing.ω_b) < 0.5
    end

    rm(tmpdir; recursive=true)
end
