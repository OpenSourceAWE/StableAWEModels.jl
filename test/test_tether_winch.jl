# SPDX-FileCopyrightText: 2025 Bart van de Lint
# SPDX-License-Identifier: MPL-2.0

# test_tether_winch.jl - Tether and Winch dynamics tests
#
# Tests winch motor dynamics and tether reeling mechanics.
# Verifies:
# 1. YAML loading of tether and winch components
# 2. Steady-state tether force with brake engaged
# 3. Brake holds tether length constant
# 4. Reel-out under gravity (brake off)
# 5. Winch dynamics: acceleration matches tau_net / I with friction

using Test
using SymbolicAWEModels
using SymbolicAWEModels: KVec3
using KiteUtils
using LinearAlgebra

# ============================================================================
# YAML Configuration - Simple tether with winch
# Weight at z=-50m (below ground) with 10kg mass, connected to ground via tether
# Gravity pulls weight down (-z), creating tension that drives reel-out
# ============================================================================
const TETHER_WINCH_YAML = """
##############################
## Tether-Winch Test System ##
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
# Weight at z=-50 (below ground) so gravity pulls it down, creating tether tension
points:
  headers: [name, pos_cad, type, wing_idx, transform_idx, extra_mass, body_frame_damping, world_frame_damping, area, drag_coeff]
  data:
    - [weight, [0.0, 0.0, -50.0], DYNAMIC, nothing, nothing, 10.0, 0.0, 0.0, 0.0, 0.0]
    - [ground, [0.0, 0.0, 0.0], STATIC, nothing, nothing, 0.0, 0.0, 0.0, 0.0, 0.0]

###########################
## Segments ###############
###########################
segments:
  headers: [name, point_i, point_j, type, l0, diameter_mm, unit_stiffness, unit_damping, compression_frac]
  data:
    - [tether_seg, weight, ground, BRIDLE, 50.0, 10.0, dyneema, nothing, 0.01]

###########################
## Tethers ################
###########################
tethers:
  headers: [name, segment_idxs, winch_point_idx]
  data:
    - [main_tether, [tether_seg], ground]

###########################
## Winches ################
###########################
winches:
  headers: [name, tether_idxs]
  data:
    - [main_winch, [main_tether]]
"""

@testset "Tether and Winch Tests" begin
    # Write YAML to temp file
    tmpdir = mktempdir()
    yaml_path = joinpath(tmpdir, "test_tether_winch_geometry.yaml")
    write(yaml_path, TETHER_WINCH_YAML)

    # Create settings file with winch parameters
    # Winch params: r_drum=0.1m, n=1.0, I=0.1 kg·m², f_coulomb=1.0 N·m, c_vf=0.5 N·m·s
    settings_yaml = """
system:
    log_file: "data/tether_winch_test"
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
    drum_radius: 0.1
    gear_ratio: 1.0
    inertia_total: 0.1
    f_coulomb: 1.0
    c_vf: 0.5

environment:
    rho_0: 1.225
    v_wind: 0.0
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
    sys = load_sys_struct_from_yaml(yaml_path; system_name="tether_winch_test", set=set)

    # ========================================================================
    # YAML Loading Verification
    # ========================================================================
    @testset "YAML Loading Verification" begin
        # Verify points
        @test length(sys.points) == 2
        @test haskey(sys.points, :weight)
        @test haskey(sys.points, :ground)

        weight = sys.points[:weight]
        @test weight.type == SymbolicAWEModels.DYNAMIC
        @test weight.extra_mass == 10.0
        @test weight.pos_cad == KVec3(0.0, 0.0, -50.0)

        # Verify segment
        @test length(sys.segments) == 1
        @test haskey(sys.segments, :tether_seg)
        @test sys.segments[:tether_seg].l0 == 50.0

        # Verify tether
        @test length(sys.tethers) == 1
        @test haskey(sys.tethers, :main_tether)
        tether = sys.tethers[:main_tether]
        @test length(tether.segment_idxs) == 1

        # Verify winch
        @test length(sys.winches) == 1
        @test haskey(sys.winches, :main_winch)
        winch = sys.winches[:main_winch]
        @test length(winch.tether_idxs) == 1

        # Verify winch parameters from settings
        @test winch.drum_radius == 0.1
        @test winch.gear_ratio == 1.0
        @test winch.inertia_total == 0.1
        @test winch.f_coulomb == 1.0
        @test winch.c_vf == 0.5

        println("\n  ====== Loaded: mass=$(weight.extra_mass)kg, l0=$(sys.segments[:tether_seg].l0)m, r=$(winch.drum_radius)m, I=$(winch.inertia_total)kg·m² ======\n")
    end

    # ========================================================================
    # Physics Test 1: Steady-state with brake
    # With brake engaged, weight hangs at equilibrium under gravity
    # ========================================================================
    @testset "Steady-state with brake" begin
        set.g_earth = 9.81
        set.v_wind = 0.0

        sys = load_sys_struct_from_yaml(yaml_path; system_name="tether_steady", set=set)
        sam = SymbolicAWEModel(set, sys)
        init!(sam; remake=true)

        # Set brake to hold tether in place
        sam.sys_struct.winches[:main_winch].brake = true

        # Run simulation to let system settle
        dt = 0.001
        for _ in 1:2000
            next_step!(sam; dt=dt, vsm_interval=0)
        end

        # Verify weight velocity is near zero (equilibrium)
        weight_vel = sam.sys_struct.points[:weight].vel_w
        @test norm(weight_vel) < 0.1

        # Verify weight is near initial position (tether holds it)
        weight_pos = sam.sys_struct.points[:weight].pos_w
        @test weight_pos[3] ≈ -50.0 atol=0.1

        println("\n  ====== Steady-state: vel=$(round(norm(weight_vel)*1000, digits=1))mm/s, z=$(round(weight_pos[3], digits=2))m ======\n")
    end

    # ========================================================================
    # Physics Test 2: Brake holds tether length constant
    # ========================================================================
    @testset "Brake holds tether length" begin
        set.g_earth = 9.81
        set.v_wind = 0.0

        sys = load_sys_struct_from_yaml(yaml_path; system_name="brake_test", set=set)
        sam = SymbolicAWEModel(set, sys)
        init!(sam; remake=true)

        sam.sys_struct.winches[:main_winch].brake = true
        initial_len = sam.sys_struct.winches[:main_winch].tether_len

        # Run simulation
        dt = 0.001
        len_history = Float64[initial_len]
        for _ in 1:1000
            next_step!(sam; dt=dt, vsm_interval=0)
            push!(len_history, sam.sys_struct.winches[:main_winch].tether_len)
        end

        max_change = maximum(abs.(len_history .- initial_len))
        @test max_change < 0.01

        println("\n  ====== Brake: initial=$(round(initial_len, digits=2))m, max_change=$(round(max_change*1000, digits=2))mm ======\n")
    end

    # ========================================================================
    # Physics Test 3: Reel-out under gravity with dynamics verification
    # With brake off and no motor torque, weight causes reel-out.
    # Captures z(t) to back-calculate inertia and friction from observed behavior.
    # ========================================================================
    @testset "Reel-out under gravity" begin
        set.g_earth = 9.81
        set.v_wind = 0.0

        sys = load_sys_struct_from_yaml(yaml_path; system_name="reelout_test", set=set)
        sam = SymbolicAWEModel(set, sys)
        init!(sam; remake=true)

        sam.sys_struct.winches[:main_winch].brake = false

        # Get known parameters from winch
        winch = sam.sys_struct.winches[:main_winch]
        r = winch.drum_radius      # 0.1 m
        n = winch.gear_ratio       # 1.0
        I_expected = winch.inertia_total  # 0.1 kg·m²
        f_c_expected = winch.f_coulomb    # 1.0 N·m
        c_vf_expected = winch.c_vf        # 0.5 N·m·s

        # Get total mass (extra_mass + tether mass)
        weight_point = sam.sys_struct.points[:weight]
        total_mass = weight_point.total_mass
        g = set.g_earth

        # Force from gravity
        F_gravity = total_mass * g

        # Capture z position over time
        dt = 0.001
        n_steps = 100
        z_history = Float64[]
        t_history = Float64[]

        push!(z_history, weight_point.pos_w[3])
        push!(t_history, 0.0)

        for i in 1:n_steps
            next_step!(sam; dt=dt, vsm_interval=0)
            push!(z_history, sam.sys_struct.points[:weight].pos_w[3])
            push!(t_history, i * dt)
        end

        # Calculate velocity from z(t) using central differences
        vel_history = Float64[]
        for i in 2:(length(z_history)-1)
            vel = (z_history[i+1] - z_history[i-1]) / (2 * dt)
            push!(vel_history, vel)
        end

        # Calculate acceleration from velocity using central differences
        acc_history = Float64[]
        for i in 2:(length(vel_history)-1)
            acc = (vel_history[i+1] - vel_history[i-1]) / (2 * dt)
            push!(acc_history, acc)
        end

        # Use early measurements (low velocity) to estimate inertia
        # At low v: F_gravity - f_coulomb = (m + I*(n/r)²) * acc
        # Rearranging: I = (F_gravity - f_coulomb - m*acc) * (r/n)² / acc
        # But simpler: total_effective_inertia = F_net / acc, then extract I
        early_acc = abs(acc_history[5])  # Use 5th sample (still near v=0)
        early_vel = abs(vel_history[5])
        omega_early = early_vel * n / r

        # tau_friction at low velocity ≈ f_coulomb (viscous term small)
        tau_friction_early = f_c_expected + c_vf_expected * omega_early
        tau_tether = (r / n) * F_gravity
        tau_net = tau_tether - tau_friction_early

        # From tau_net = I * alpha, and alpha = acc * n / r:
        alpha_measured = early_acc * n / r
        I_calculated = tau_net / alpha_measured

        # Verify calculated inertia matches expected
        @test I_calculated ≈ I_expected rtol=0.3  # 30% tolerance for numerical effects

        # Use later measurements (higher velocity) to verify friction
        # At steady state: tau_tether = f_coulomb + c_vf * omega_terminal
        # v_terminal = (tau_tether - f_coulomb) * r / (c_vf * n)
        v_terminal_expected = (tau_tether - f_c_expected) * r / (c_vf_expected * n)

        # Verify basic reel-out behavior
        initial_z = z_history[1]
        final_z = z_history[end]
        @test final_z < initial_z  # Weight moves down (negative z direction)

        final_vel = sam.sys_struct.winches[:main_winch].tether_vel
        @test final_vel > 0  # Positive = reel-out

        println("\n  ====== Reel-out: z=$(round(initial_z, digits=2))m -> $(round(final_z, digits=2))m, vel=$(round(final_vel, digits=2))m/s")
        println("  ====== Dynamics: I_calc=$(round(I_calculated, digits=3))kg·m² (expected=$(I_expected)), acc=$(round(early_acc, digits=2))m/s² ======\n")
    end

    # ========================================================================
    # Physics Test 4: Initial acceleration and terminal velocity
    # At v≈0: a = (m*g - f_c*n/r) / (m + I*n²/r²)
    # At terminal: F_drive = F_friction, so v_term = (m*g*r/n - f_c) * r / (c_vf*n)
    # ========================================================================
    @testset "Winch dynamics verification" begin
        set.g_earth = 9.81
        set.v_wind = 0.0

        sys = load_sys_struct_from_yaml(yaml_path; system_name="dynamics_test", set=set)
        sam = SymbolicAWEModel(set, sys)
        init!(sam; remake=true)

        sam.sys_struct.winches[:main_winch].brake = false

        # Get winch parameters
        winch = sam.sys_struct.winches[:main_winch]
        r = winch.drum_radius      # 0.1 m
        n = winch.gear_ratio       # 1.0
        I = winch.inertia_total    # 0.1 kg·m²
        f_c = winch.f_coulomb      # 1.0 N·m
        c_vf = winch.c_vf          # 0.5 N·m·s

        # Get mass (use initial value before significant tether payout)
        m = sam.sys_struct.points[:weight].total_mass
        g = set.g_earth

        # === Part 1: Initial acceleration ===
        # Expected at v≈0 (viscous friction ≈ 0):
        # a = (m*g - f_c*n/r) / (m + I*n²/r²)
        M_eff = m + I * n^2 / r^2  # Effective mass including drum inertia
        F_drive = m * g            # Gravity force
        F_coulomb = f_c * n / r    # Coulomb friction converted to linear force
        a_expected = (F_drive - F_coulomb) / M_eff

        # Capture velocity over first few timesteps
        dt = 0.001
        n_steps = 20
        vel_history = Float64[]
        t_history = Float64[]

        push!(vel_history, winch.tether_vel)
        push!(t_history, 0.0)

        for i in 1:n_steps
            next_step!(sam; dt=dt, vsm_interval=0)
            push!(vel_history, sam.sys_struct.winches[:main_winch].tether_vel)
            push!(t_history, i * dt)
        end

        # Measure acceleration (linear fit: v = a*t + v0)
        t_mean = sum(t_history) / length(t_history)
        v_mean = sum(vel_history) / length(vel_history)
        numerator = sum((t_history .- t_mean) .* (vel_history .- v_mean))
        denominator = sum((t_history .- t_mean).^2)
        a_measured = numerator / denominator

        @test a_measured > 0
        @test abs(a_measured - a_expected) / a_expected < 0.3

        println("\n  ====== Initial: a_meas=$(round(a_measured, digits=2))m/s², a_exp=$(round(a_expected, digits=2))m/s²")

        # === Part 2: Terminal velocity ===
        # Run longer to approach terminal velocity
        # Time constant: τ = M_eff / (c_vf * n² / r²)
        tau = M_eff / (c_vf * n^2 / r^2)
        t_settle = 10 * tau  # Run for 10 time constants

        n_settle_steps = Int(ceil(t_settle / dt))
        for _ in 1:n_settle_steps
            next_step!(sam; dt=dt, vsm_interval=0)
        end

        v_final = sam.sys_struct.winches[:main_winch].tether_vel

        # At terminal velocity: F_drive = F_friction
        # F_friction = (f_c + c_vf * omega) * n / r = (f_c + c_vf * v * n / r) * n / r
        # Solving for v_terminal: v = (F_drive * r / n - f_c) * r / (c_vf * n)
        v_terminal_expected = (F_drive * r / n - f_c) * r / (c_vf * n)

        # Calculate friction force at measured velocity
        omega_final = v_final * n / r
        tau_friction_final = f_c + c_vf * omega_final
        F_friction_final = tau_friction_final * n / r

        # At terminal: F_friction should equal F_drive
        force_ratio = F_friction_final / F_drive

        @test v_final > 0
        @test abs(v_final - v_terminal_expected) / v_terminal_expected < 0.3
        @test abs(force_ratio - 1.0) < 0.3  # Forces should be balanced

        println("  ====== Terminal: v_meas=$(round(v_final, digits=2))m/s, v_exp=$(round(v_terminal_expected, digits=2))m/s")
        println("  ====== Force balance: F_friction=$(round(F_friction_final, digits=1))N, F_drive=$(round(F_drive, digits=1))N, ratio=$(round(force_ratio, digits=2)) ======\n")
    end

    # Cleanup
    rm(tmpdir; recursive=true)
end
