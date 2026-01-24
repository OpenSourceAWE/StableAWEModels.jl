# SPDX-FileCopyrightText: 2025 Bart van de Lint
# SPDX-License-Identifier: MPL-2.0

using Test
using SymbolicAWEModels
using Statistics, LinearAlgebra

# Set up tmpdir if not already done by runtests.jl
if !startswith(get_data_path(), tempdir())
    src_data_path = joinpath(dirname(dirname(pathof(SymbolicAWEModels))), "data", "ram_air_kite")
    tmpdir = mktempdir()
    set_data_path(joinpath(tmpdir, "ram_air_kite"))
    cp(src_data_path, get_data_path(); force=true)
end

@testset verbose=true "SysState ↔ SystemStructure conversion" begin
    set = Settings("system.yaml")
    original_set = Settings("system.yaml")

    function reset!(set::Settings)
        for field in fieldnames(Settings)
            setfield!(set, field, getfield(original_set, field))
        end
        return set
    end

    @testset "Basic point position updates" begin
        reset!(set)
        sys = create_simple_ram_sys_struct(set)
        P = length(sys.points)

        # Create a SysState with known positions
        ss = SysState{P}()
        for i in 1:P
            ss.X[i] = 10.0 * i
            ss.Y[i] = 5.0 * i
            ss.Z[i] = 20.0 + i
        end

        # Update the SystemStructure
        update_from_sysstate!(sys, ss)

        # Verify point positions were updated
        for point in sys.points
            @test point.pos_w[1] ≈ ss.X[point.idx]
            @test point.pos_w[2] ≈ ss.Y[point.idx]
            @test point.pos_w[3] ≈ ss.Z[point.idx]
            # Velocities should be zero
            @test all(point.vel_w .≈ 0.0)
            # Forces should be NaN
            @test all(isnan.(point.force))
        end
    end

    @testset "Wing state updates" begin
        reset!(set)
        sys = create_simple_ram_sys_struct(set)
        P = length(sys.points)

        ss = SysState{P}()
        # Set wing orientation and angles
        ss.orient .= [0.9239, 0.3827, 0.0, 0.0]  # Example quaternion
        ss.elevation = 0.5
        ss.azimuth = 0.2
        ss.heading = 0.1
        ss.vel_kite .= [1.0, 2.0, 3.0]
        ss.turn_rates .= [0.1, 0.2, 0.3]
        ss.AoA = 0.15
        ss.course = 0.25
        ss.v_wind_kite .= [10.0, 0.0, 0.0]

        # Populate point positions to avoid NaN in wing position calculation
        for i in 1:P
            ss.X[i] = Float32(i)
            ss.Y[i] = Float32(i)
            ss.Z[i] = Float32(10 + i)
        end

        update_from_sysstate!(sys, ss)

        # Verify wing state
        @test length(sys.wings) > 0
        wing = sys.wings[1]
        @test wing.Q_b_w ≈ ss.orient
        @test wing.elevation ≈ ss.elevation
        @test wing.azimuth ≈ ss.azimuth
        @test wing.heading ≈ ss.heading
        @test wing.vel_w ≈ ss.vel_kite
        @test wing.ω_b ≈ ss.turn_rates
        @test wing.aoa ≈ ss.AoA
        @test wing.course ≈ ss.course
        @test wing.v_wind ≈ ss.v_wind_kite

        # Aerodynamic quantities should be NaN
        @test all(isnan.(wing.aero_force_b))
        @test all(isnan.(wing.aero_moment_b))
        @test all(isnan.(wing.tether_force))
        @test all(isnan.(wing.tether_moment))
        @test all(isnan.(wing.va_b))
    end

    @testset "Winch state updates" begin
        reset!(set)
        sys = create_simple_ram_sys_struct(set)
        P = length(sys.points)

        ss = SysState{P}()
        # Set winch data (up to 4 winches supported)
        ss.l_tether .= [50.0, 51.0, 52.0, 53.0]
        ss.v_reelout .= [0.5, 0.6, 0.7, 0.8]
        ss.set_torque .= [100.0, 101.0, 102.0, 103.0]

        update_from_sysstate!(sys, ss)

        # Verify winch state
        n_winches = min(length(sys.winches), 4)
        for i in 1:n_winches
            @test sys.winches[i].tether_len ≈ ss.l_tether[i]
            @test sys.winches[i].tether_vel ≈ ss.v_reelout[i]
            @test sys.winches[i].set_value ≈ ss.set_torque[i]
            @test all(isnan.(sys.winches[i].force))
            @test isnan(sys.winches[i].friction)
        end
    end

    @testset "Group twist updates" begin
        reset!(set)
        sys = create_simple_ram_sys_struct(set)
        P = length(sys.points)

        ss = SysState{P}()
        # Set twist angles (up to 4 groups supported)
        ss.twist_angles .= [0.1, 0.2, 0.3, 0.4]

        update_from_sysstate!(sys, ss)

        # Verify group state
        n_groups = min(length(sys.groups), 4)
        for i in 1:n_groups
            @test sys.groups[i].twist ≈ ss.twist_angles[i]
            @test sys.groups[i].twist_ω ≈ 0.0
            @test isnan(sys.groups[i].tether_force)
            @test isnan(sys.groups[i].tether_moment)
            @test isnan(sys.groups[i].aero_moment)
        end
    end

    @testset "Segment length calculation" begin
        reset!(set)
        sys = create_simple_ram_sys_struct(set)
        P = length(sys.points)

        ss = SysState{P}()
        # Create a specific geometry
        for i in 1:P
            ss.X[i] = Float32(i * 2.0)
            ss.Y[i] = 0.0
            ss.Z[i] = 0.0
        end

        update_from_sysstate!(sys, ss)

        # Check segments have computed lengths and forces
        for segment in sys.segments
            @test segment.len > 0.0
            # Force should be computed (not NaN)
            @test !isnan(segment.force)
            @test isfinite(segment.force)
        end
    end

    @testset "Round-trip consistency" begin
        reset!(set)
        sam = SymbolicAWEModel(set, "simple_ram")
        init!(sam)

        # Extract state
        ss1 = SysState(sam)

        # Create a new SystemStructure and update it from the SysState
        sys2 = create_simple_ram_sys_struct(set)
        update_from_sysstate!(sys2, ss1)

        # Extract state again (through the existing update_sys_state!)
        # We need to use the sam since update_sys_state! expects a SymbolicAWEModel
        # So we test the positions directly instead
        for point in sys2.points
            @test point.pos_w[1] ≈ ss1.X[point.idx] atol=1e-4
            @test point.pos_w[2] ≈ ss1.Y[point.idx] atol=1e-4
            @test point.pos_w[3] ≈ ss1.Z[point.idx] atol=1e-4
        end

        # Test wing state
        if length(sys2.wings) > 0 && length(sam.sys_struct.wings) > 0
            wing_orig = sam.sys_struct.wings[1]
            wing_new = sys2.wings[1]
            @test wing_new.Q_b_w ≈ wing_orig.Q_b_w atol=1e-4
            @test wing_new.elevation ≈ wing_orig.elevation atol=1e-4
            @test wing_new.azimuth ≈ wing_orig.azimuth atol=1e-4
        end
    end

    @testset "Point count validation" begin
        reset!(set)
        sys = create_simple_ram_sys_struct(set)
        P = length(sys.points)

        # Create SysState with wrong number of points
        ss_wrong = SysState{P+5}()

        # Should throw an error
        @test_throws ErrorException update_from_sysstate!(sys, ss_wrong)
    end

    @testset "Global wind vector" begin
        reset!(set)
        sys = create_simple_ram_sys_struct(set)
        P = length(sys.points)

        ss = SysState{P}()
        ss.v_wind_gnd .= [12.5, 1.5, 0.5]

        update_from_sysstate!(sys, ss)

        @test sys.wind_vec_gnd ≈ ss.v_wind_gnd
    end

    @testset "Integration with simulation" begin
        reset!(set)
        sam = SymbolicAWEModel(set, "simple_ram")
        init!(sam)
        find_steady_state!(sam)

        # Run simulation for a few steps
        for _ in 1:10
            next_step!(sam)
        end

        # Extract state at current time
        ss = SysState(sam)

        # Create a fresh SystemStructure
        sys_fresh = create_simple_ram_sys_struct(set)

        # Update it from the simulation state
        update_from_sysstate!(sys_fresh, ss)

        # Verify positions match
        for i in 1:length(sys_fresh.points)
            @test sys_fresh.points[i].pos_w[1] ≈ sam.sys_struct.points[i].pos_w[1] atol=1e-3
            @test sys_fresh.points[i].pos_w[2] ≈ sam.sys_struct.points[i].pos_w[2] atol=1e-3
            @test sys_fresh.points[i].pos_w[3] ≈ sam.sys_struct.points[i].pos_w[3] atol=1e-3
        end
    end

    @testset "Edge case: System with minimal components" begin
        reset!(set)
        set.segments = 1

        # Create a minimal system
        points = [
            Point(1, zeros(3), STATIC; wing_idx=0, transform_idx=1)
            Point(2, [1.0, 0.0, 0.0], DYNAMIC; wing_idx=0, transform_idx=1)
        ]
        segments = [Segment(1, set, (1, 2), BRIDLE)]
        transforms = [Transform(1, 0.0, 0.0, 0.0;
            base_pos=zeros(3), base_point_idx=1, rot_point_idx=2)]
        sys = SystemStructure("minimal", set; points, segments, transforms)

        P = length(points)
        ss = SysState{P}()
        ss.X .= [0.0, 5.0]
        ss.Y .= [0.0, 0.0]
        ss.Z .= [0.0, 10.0]

        # Should not crash
        update_from_sysstate!(sys, ss)

        @test sys.points[2].pos_w[1] ≈ 5.0
        @test sys.points[2].pos_w[3] ≈ 10.0
    end
end
