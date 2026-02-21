# SPDX-FileCopyrightText: 2025 Bart van de Lint
# SPDX-License-Identifier: MPL-2.0

# test_tether_winch.jl - Winch dynamics tests
#
# Tests winch motor dynamics in isolation using zero-stiffness
# and low-stiffness tethers. Each test verifies one physical
# behavior with tight tolerances and swept parameters.

using Test
using SymbolicAWEModels
using SymbolicAWEModels: KVec3
using KiteUtils
using LinearAlgebra

# ================================================================
# Minimal YAML: weight at z=-50 connected to ground via tether
# ================================================================
const WINCH_TEST_YAML = """
materials:
  headers: [name, youngs_modulus, density, damping_per_stiffness]
  data:
    - [test_mat, 1000.0, 724, 0.001]

points:
  headers: [name, pos_cad, type, wing_idx, transform_idx,
            extra_mass, body_frame_damping, world_frame_damping,
            area, drag_coeff]
  data:
    - [weight, [0.0, 0.0, -50.0], DYNAMIC, nothing, nothing,
       10.0, 0.0, 0.0, 0.0, 0.0]
    - [ground, [0.0, 0.0, 0.0], STATIC, nothing, nothing,
       0.0, 0.0, 0.0, 0.0, 0.0]

segments:
  headers: [name, point_i, point_j, type, l0, diameter_mm,
            unit_stiffness, unit_damping, compression_frac]
  data:
    - [tether_seg, weight, ground, BRIDLE, 50.0, 1.0,
       0.0, 0.0, 0.0]

tethers:
  headers: [name, segment_idxs, winch_point_idx]
  data:
    - [main_tether, [tether_seg], ground]

winches:
  headers: [name, tether_idxs]
  data:
    - [main_winch, [main_tether]]
"""

const WINCH_TEST_SETTINGS = """
system:
    log_file: "data/winch_test"
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
    cd_tether: 0.0
    unit_damping: 350.0
    unit_stiffness: 120000.0
    rho_tether: 724.0
    e_tether: 1000.0
    rel_damping: 0.001
    d_tether: 1.0

winch:
    winch_model: "TorqueControlledMachine"
    max_force: 4000
    v_ro_max: 8.0
    drum_radius: 0.1
    gear_ratio: 1.0
    inertia_total: 0.1
    f_coulomb: 0.0
    c_vf: 0.0

environment:
    rho_0: 1.225
    v_wind: 0.0
    upwind_dir: -90.0
    h_ref: 6.0
    profile_law: 0
"""

@testset "Winch Tests" begin
    # --- Setup: write YAML files, build model ONCE ---
    tmpdir = mktempdir()
    yaml_path = joinpath(tmpdir, "winch_test_geometry.yaml")
    write(yaml_path, WINCH_TEST_YAML)

    settings_path = joinpath(tmpdir, "settings.yaml")
    write(settings_path, WINCH_TEST_SETTINGS)

    system_yaml = "system:\n  sim_settings: settings.yaml\n"
    write(joinpath(tmpdir, "system.yaml"), system_yaml)

    set_data_path(tmpdir)
    set = Settings("system.yaml")

    sys = load_sys_struct_from_yaml(
        yaml_path; system_name="winch_test", set=set)
    sam = SymbolicAWEModel(set, sys)
    init!(sam; remake=true, prn=false)

    winch = sam.sys_struct.winches[:main_winch]
    seg = sam.sys_struct.segments[:tether_seg]
    r = winch.drum_radius   # 0.1 m
    n = winch.gear_ratio     # 1.0

    # ============================================================
    # Test 1: Brake holds tether velocity at zero
    # ============================================================
    @testset "Brake holds" begin
        winch.brake = true
        winch.f_coulomb = 0.0
        winch.c_vf = 0.0
        winch.inertia_total = 0.1
        seg.unit_stiffness = 0.0
        seg.unit_damping = 0.0
        init!(sam; remake=true, prn=false)

        tau_motor = 5.0
        for _ in 1:100
            next_step!(sam; set_values=[tau_motor],
                       dt=0.001, vsm_interval=0)
        end
        @test abs(winch.tether_vel) < 1e-10
    end

    # ============================================================
    # Test 2: Acceleration vs inertia
    # Zero stiffness, zero friction. Only motor torque and inertia.
    # Expected: a = (r/n) * tau_motor / I
    # ============================================================
    @testset "Acceleration vs inertia" begin
        winch.brake = false
        winch.f_coulomb = 0.0
        winch.c_vf = 0.0
        seg.unit_stiffness = 0.0
        seg.unit_damping = 0.0
        tau_motor = 1.0

        for I_test in [0.1, 0.5, 1.0]
            winch.inertia_total = I_test
            init!(sam; remake=true, prn=false)

            dt = 0.001
            n_steps = 20
            for _ in 1:n_steps
                next_step!(sam; set_values=[tau_motor],
                           dt=dt, vsm_interval=0)
            end

            # v(t) = a*t, so a = v / t
            a_measured = winch.tether_vel / (n_steps * dt)
            a_expected = (r / n) * tau_motor / I_test

            @test a_measured ≈ a_expected rtol=0.05
        end
    end

    # ============================================================
    # Test 3: Acceleration vs Coulomb friction
    # Zero stiffness, zero viscous friction, small epsilon.
    # After spin-up (so smooth_sign ≈ 1):
    #   a = (r/n)/I * (tau_motor - f_coulomb * r/n)
    # ============================================================
    @testset "Acceleration vs Coulomb friction" begin
        winch.brake = false
        winch.c_vf = 0.0
        winch.friction_epsilon = 0.01
        winch.inertia_total = 0.1
        seg.unit_stiffness = 0.0
        seg.unit_damping = 0.0
        tau_motor = 1.0
        I = winch.inertia_total

        for f_c in [0.5, 2.0, 5.0]
            winch.f_coulomb = f_c
            init!(sam; remake=true, prn=false)

            dt = 0.001
            # Spin up so smooth_sign(ω, 0.01) ≈ 1
            for _ in 1:50
                next_step!(sam; set_values=[tau_motor],
                           dt=dt, vsm_interval=0)
            end
            v0 = winch.tether_vel

            # Measure acceleration over next steps
            n_meas = 20
            for _ in 1:n_meas
                next_step!(sam; set_values=[tau_motor],
                           dt=dt, vsm_interval=0)
            end
            v1 = winch.tether_vel

            a_measured = (v1 - v0) / (n_meas * dt)
            a_expected = (r / n) / I *
                (tau_motor - f_c * r / n)

            @test a_measured ≈ a_expected rtol=0.05
        end
    end

    # ============================================================
    # Test 4: Terminal velocity vs viscous friction
    # Zero stiffness, zero Coulomb. Motor torque balanced by
    # viscous friction at terminal velocity.
    # Expected: v_term = tau_motor * n / (c_vf * r)
    # ============================================================
    @testset "Terminal velocity vs viscous friction" begin
        winch.brake = false
        winch.f_coulomb = 0.0
        winch.friction_epsilon = 6.0
        winch.inertia_total = 0.01  # Small I for fast settling
        seg.unit_stiffness = 0.0
        seg.unit_damping = 0.0
        tau_motor = 1.0
        I = winch.inertia_total

        for c_vf_test in [50.0, 100.0, 200.0]
            winch.c_vf = c_vf_test
            init!(sam; remake=true, prn=false)

            # Time constant: I / (c_vf * (r/n)^2)
            tau = I / (c_vf_test * (r / n)^2)
            t_settle = 10 * tau
            dt = 0.001
            n_steps = Int(ceil(t_settle / dt))

            for _ in 1:n_steps
                next_step!(sam; set_values=[tau_motor],
                           dt=dt, vsm_interval=0)
            end

            v_expected = tau_motor * n / (c_vf_test * r)
            @test winch.tether_vel ≈ v_expected rtol=0.05
        end
    end

    # ============================================================
    # Test 5: Steady-state tether stretch
    # Brake on (l0 constant), weight sags under gravity.
    # At equilibrium: (unit_stiffness / len) * (len - l0) = m*g
    # Exact: extension = m*g*l0 / (unit_stiffness - m*g)
    # ============================================================
    @testset "Steady-state tether stretch" begin
        winch.brake = true
        winch.f_coulomb = 1.0
        winch.c_vf = 0.5
        winch.friction_epsilon = 6.0
        winch.inertia_total = 0.1

        unit_k = 50000.0  # Low but not too low [N]
        unit_d = unit_k * 0.01  # Damping for settling
        seg.unit_stiffness = unit_k
        seg.unit_damping = unit_d
        init!(sam; remake=true, prn=false)

        # Compute total mass at weight point:
        # m = extra_mass + rho * pi * (d/2)^2 * l0 / 2
        l0 = 50.0
        d = seg.diameter  # meters
        rho = set.rho_tether
        extra_m = sam.sys_struct.points[:weight].extra_mass
        tether_m = rho * π * (d / 2)^2 * l0 / 2
        m = extra_m + tether_m
        g = set.g_earth

        # Run to steady state
        dt = 0.001
        for _ in 1:5000
            next_step!(sam; dt=dt, vsm_interval=0)
        end

        # Exact equilibrium: k * ext = m*g where k = EA/len
        extension = m * g * l0 / (unit_k - m * g)
        expected_z = -(l0 + extension)

        weight_z = sam.sys_struct.points[:weight].pos_w[3]
        @test weight_z ≈ expected_z rtol=0.01
    end

    # ============================================================
    # Test 6: calc_steady_torque at terminal velocity
    # Stiff tether with heavy weight, viscous friction,
    # small motor torque. At terminal velocity the system
    # is in steady state so calc_steady_torque must equal
    # the motor torque we are applying.
    # ============================================================
    @testset "calc_steady_torque at terminal velocity" begin
        winch.brake = false
        winch.f_coulomb = 0.0
        winch.friction_epsilon = 6.0
        winch.inertia_total = 0.01
        winch.c_vf = 100.0

        # Stiff tether to transmit gravity force
        seg.unit_stiffness = 120000.0
        seg.unit_damping = 500.0
        init!(sam; remake=true, prn=false)

        tau_motor = 0.5

        # Time constant: I / (c_vf * (r/n)^2)
        tau_tc = winch.inertia_total /
            (winch.c_vf * (r / n)^2)
        t_settle = 100 * tau_tc
        dt = 0.1
        n_steps = Int(ceil(t_settle / dt))

        for _ in 1:n_steps
            next_step!(sam; set_values=[tau_motor],
                       dt=dt, vsm_interval=0)
        end

        steady = calc_steady_torque(sam)
        @test sam.sys_struct.winches[1].tether_acc ≈ 0.0 atol=0.01
        @test steady[1] ≈ tau_motor rtol=0.01

        println(
            "  calc_steady_torque: " *
            "applied=$(tau_motor), " *
            "computed=$(round(steady[1]; digits=4))"
        )
    end

    rm(tmpdir; recursive=true)
end
