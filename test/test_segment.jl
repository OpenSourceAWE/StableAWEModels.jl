# SPDX-FileCopyrightText: 2025 Bart van de Lint
# SPDX-License-Identifier: MPL-2.0

# test_segment.jl - Spring-damper segment dynamics tests
#
# Tests a minimal 2-point system (1 STATIC, 1 DYNAMIC) with 1 segment.
# Verifies:
# 1. No gravity, no wind: point stays still
# 2. With gravity: oscillation frequency and equilibrium position
# 3. Damping ratio from decay

using Test
using SymbolicAWEModels
using KiteUtils
using LinearAlgebra

# ============================================================================
# YAML Configuration - Minimal 2-point system with 1 segment
# ============================================================================
const SEGMENT_TEST_YAML = """
##############################
## Segment Test System #######
##############################

###########################
## Materials ##############
###########################
materials:
  headers: [name, youngs_modulus, density, damping_per_stiffness]
  data:
    - [test_material, 55000000000.0, 724, 0.00077]

###########################
## Points #################
###########################
points:
  headers: [name, pos_cad, type, wing_idx, transform_idx, extra_mass, body_frame_damping, world_frame_damping, area, drag_coeff]
  data:
    - [ground, [0.0, 0.0, 0.0], STATIC, nothing, nothing, 0.0, 0.0, 0.0, 0.0, 0.0]
    - [mass_point, [0.0, 0.0, 10.0], DYNAMIC, nothing, nothing, 1.0, 0.0, 0.0, 0.0, 0.0]

###########################
## Segments ###############
###########################
segments:
  headers: [name, point_i, point_j, type, l0, diameter_mm, axial_stiffness, axial_damping, compression_frac]
  data:
    - [test_segment, ground, mass_point, BRIDLE, 10.0, 5.0, 1000.0, 10.0, 0.1]
"""

@testset "Segment Tests" begin
    # Write YAML to temp file
    tmpdir = mktempdir()
    yaml_path = joinpath(tmpdir, "test_segment_geometry.yaml")
    write(yaml_path, SEGMENT_TEST_YAML)

    # Create minimal settings file for loading
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
  axial_damping: 350.0
  axial_stiffness: 120000.0
  rho_tether: 724.0
  e_tether: 55000000000.0
  rel_damping: 0.00077

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

    # Also create system.yaml that points to settings
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
    sys = load_sys_struct_from_yaml(yaml_path; system_name="segment_test", set=set)

    # ========================================================================
    # YAML Loading Verification
    # ========================================================================
    @testset "YAML Loading Verification" begin
        # Verify points were loaded correctly
        @test length(sys.points) == 2
        @test haskey(sys.points, :ground)
        @test haskey(sys.points, :mass_point)

        # Verify point properties
        ground = sys.points[:ground]
        @test ground.type == SymbolicAWEModels.STATIC
        @test ground.pos_cad == KVec3(0.0, 0.0, 0.0)
        @test ground.extra_mass == 0.0

        mass_point = sys.points[:mass_point]
        @test mass_point.type == SymbolicAWEModels.DYNAMIC
        @test mass_point.pos_cad == KVec3(0.0, 0.0, 10.0)
        @test mass_point.extra_mass == 1.0
        @test mass_point.area == 0.0
        @test mass_point.drag_coeff == 0.0

        # Verify segment was loaded correctly
        @test length(sys.segments) == 1
        @test haskey(sys.segments, :test_segment)

        segment = sys.segments[:test_segment]
        @test segment.l0 == 10.0
        @test segment.axial_stiffness == 1000.0
        @test segment.axial_damping == 10.0
        @test segment.diameter == 0.005  # 5mm in meters
        @test segment.compression_frac == 0.1
    end

    # ========================================================================
    # Physics Test 1: No gravity, no wind - point stays still
    # ========================================================================
    @testset "No gravity, no wind - stationary" begin
        # Reload with zero gravity
        set.g_earth = 0.0
        set.v_wind = 0.0

        sys = load_sys_struct_from_yaml(yaml_path; system_name="segment_test_nograv", set=set)
        sam = SymbolicAWEModel(set, sys)
        init!(sam; remake=true)

        # Record initial position
        initial_z = sam.sys_struct.points[:mass_point].pos_w[3]

        # Run simulation for 1 second
        dt = 0.01
        n_steps = 100
        for _ in 1:n_steps
            next_step!(sam; dt=dt, vsm_interval=0)
        end

        # Position should be unchanged
        final_z = sam.sys_struct.points[:mass_point].pos_w[3]
        @test abs(final_z - initial_z) < 0.001  # Should not move more than 1mm
    end

    # ========================================================================
    # Physics Test 2: With gravity - oscillation and equilibrium
    # ========================================================================
    @testset "With gravity - oscillation dynamics" begin
        # Reload with gravity
        set.g_earth = 9.81
        set.v_wind = 0.0

        sys = load_sys_struct_from_yaml(yaml_path; system_name="segment_test_grav", set=set)
        sam = SymbolicAWEModel(set, sys)
        init!(sam; remake=true)

        # Physics parameters
        m = 1.0  # mass [kg]
        k = 1000.0 / 10.0  # stiffness [N/m] = axial_stiffness / l0
        c = 10.0  # damping [N*s/m]
        g = 9.81  # gravity [m/s^2]
        l0 = 10.0  # rest length [m]

        # Expected equilibrium: z_eq = l0 - m*g/k
        # Spring force at equilibrium: k * delta_l = m * g
        # delta_l = m*g/k = 1.0 * 9.81 / 100 = 0.0981 m
        z_eq_expected = l0 - m * g / k
        @test z_eq_expected ≈ 9.902 atol=0.001

        # Expected natural frequency: omega_n = sqrt(k/m)
        omega_n = sqrt(k / m)
        @test omega_n ≈ 10.0 atol=0.001  # 10 rad/s

        # Expected damping ratio: zeta = c / (2 * sqrt(k*m))
        zeta = c / (2 * sqrt(k * m))
        @test zeta ≈ 0.5 atol=0.001  # Underdamped

        # Damped frequency: omega_d = omega_n * sqrt(1 - zeta^2)
        omega_d = omega_n * sqrt(1 - zeta^2)

        # Run simulation for several periods
        dt = 0.001  # Small timestep for accuracy
        total_time = 3.0  # Several oscillation periods
        n_steps = Int(ceil(total_time / dt))

        z_history = Float64[]
        t_history = Float64[]

        for i in 1:n_steps
            next_step!(sam; dt=dt, vsm_interval=0)
            push!(z_history, sam.sys_struct.points[:mass_point].pos_w[3])
            push!(t_history, i * dt)
        end

        # Check equilibrium (final position should converge)
        z_final = z_history[end]
        @test z_final ≈ z_eq_expected atol=0.05  # Within 5cm of equilibrium

        # Check oscillation occurred (should have crossed equilibrium multiple times)
        crossings = 0
        for i in 2:length(z_history)
            if (z_history[i-1] - z_eq_expected) * (z_history[i] - z_eq_expected) < 0
                crossings += 1
            end
        end
        @test crossings >= 4  # At least 2 full oscillations (4 crossings)

        # Verify damping - amplitude should decrease
        # Find first peak and a later peak
        peaks = Int[]
        for i in 2:(length(z_history)-1)
            if z_history[i] > z_history[i-1] && z_history[i] > z_history[i+1]
                push!(peaks, i)
            end
        end

        if length(peaks) >= 2
            amp1 = abs(z_history[peaks[1]] - z_eq_expected)
            amp2 = abs(z_history[peaks[end]] - z_eq_expected)
            @test amp2 < amp1  # Later amplitude should be smaller
        end
    end

    # ========================================================================
    # Physics Test 3: Recalculate stiffness/damping from oscillation
    # ========================================================================
    @testset "Extract parameters from oscillation" begin
        # Reload with gravity and small damping for clearer oscillation
        set.g_earth = 9.81
        set.v_wind = 0.0

        # Create YAML with lower damping for clearer oscillation
        low_damp_yaml = """
materials:
  headers: [name, youngs_modulus, density, damping_per_stiffness]
  data:
    - [test_material, 55000000000.0, 724, 0.00077]

points:
  headers: [name, pos_cad, type, wing_idx, transform_idx, extra_mass, body_frame_damping, world_frame_damping, area, drag_coeff]
  data:
    - [ground, [0.0, 0.0, 0.0], STATIC, nothing, nothing, 0.0, 0.0, 0.0, 0.0, 0.0]
    - [mass_point, [0.0, 0.0, 10.0], DYNAMIC, nothing, nothing, 1.0, 0.0, 0.0, 0.0, 0.0]

segments:
  headers: [name, point_i, point_j, type, l0, diameter_mm, axial_stiffness, axial_damping, compression_frac]
  data:
    - [test_segment, ground, mass_point, BRIDLE, 10.0, 5.0, 1000.0, 2.0, 0.1]
"""
        # Use lower damping: c = 2.0 → zeta = 2/(2*sqrt(100*1)) = 0.1

        low_damp_path = joinpath(tmpdir, "low_damp_geometry.yaml")
        write(low_damp_path, low_damp_yaml)

        sys = load_sys_struct_from_yaml(low_damp_path; system_name="segment_test_lowdamp", set=set)

        # Verify the damping was loaded
        @test sys.segments[:test_segment].axial_damping == 2.0

        sam = SymbolicAWEModel(set, sys)
        init!(sam; remake=true)

        # Physics parameters
        m = 1.0
        k = 1000.0 / 10.0  # 100 N/m
        c = 2.0
        zeta_expected = c / (2 * sqrt(k * m))  # 0.1
        omega_n_expected = sqrt(k / m)  # 10 rad/s

        # Run simulation and record
        dt = 0.001
        total_time = 5.0
        n_steps = Int(ceil(total_time / dt))

        z_history = Float64[]
        t_history = Float64[]

        for i in 1:n_steps
            next_step!(sam; dt=dt, vsm_interval=0)
            push!(z_history, sam.sys_struct.points[:mass_point].pos_w[3])
            push!(t_history, i * dt)
        end

        # Find peaks for frequency estimation
        peaks_t = Float64[]
        peaks_z = Float64[]
        for i in 2:(length(z_history)-1)
            if z_history[i] > z_history[i-1] && z_history[i] > z_history[i+1]
                push!(peaks_t, t_history[i])
                push!(peaks_z, z_history[i])
            end
        end

        if length(peaks_t) >= 3
            # Estimate period from peak spacing
            periods = diff(peaks_t)
            avg_period = mean(periods)
            omega_d_measured = 2π / avg_period

            # For underdamped: omega_d = omega_n * sqrt(1 - zeta^2)
            omega_d_expected = omega_n_expected * sqrt(1 - zeta_expected^2)

            @test omega_d_measured ≈ omega_d_expected rtol=0.1  # Within 10%

            # Estimate damping ratio from logarithmic decrement
            # delta = ln(A1/A2) = 2*pi*zeta/sqrt(1-zeta^2)
            z_eq = 10.0 - 1.0 * 9.81 / 100.0  # equilibrium position
            amps = abs.(peaks_z .- z_eq)

            if length(amps) >= 2 && amps[1] > 0.001 && amps[2] > 0.001
                log_decrement = log(amps[1] / amps[2])
                # From log decrement: zeta = delta / sqrt(4*pi^2 + delta^2)
                zeta_measured = log_decrement / sqrt(4π^2 + log_decrement^2)

                @test zeta_measured ≈ zeta_expected rtol=0.3  # Within 30% (numerical errors)
            end
        end
    end

    # Cleanup
    rm(tmpdir; recursive=true)
end

# Helper function for mean
function mean(x)
    return sum(x) / length(x)
end
