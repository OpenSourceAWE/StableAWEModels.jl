# SPDX-FileCopyrightText: 2025 Bart van de Lint
# SPDX-License-Identifier: MIT

using Pkg
if ! ("ControlPlots" ∈ keys(Pkg.project().dependencies))
    using TestEnv; TestEnv.activate()
end
using Test, ControlSystemsBase, Printf
using SymbolicAWEModels, ControlPlots
using Statistics, LinearAlgebra

set = Settings("system.yaml")
sam = SymbolicAWEModel(set, "ram")
init!(sam)

tether_set = Settings("system.yaml")
tether_sam = SymbolicAWEModel(tether_set, "tether")
init!(tether_sam)

simple_set = Settings("system.yaml")
simple_sam = SymbolicAWEModel(simple_set, "simple_ram")
init!(simple_sam)

original_set = Settings("system.yaml")

function reset!(set::Settings)
    for field in fieldnames(Settings)
        setfield!(set, field, getfield(original_set, field))
    end
    return set
end

@testset verbose=true "SymbolicAWEModels Tests" begin

    function test_plot(sam)
        @testset "Plotting of SymbolicAWEModel" begin
            function plot_(zoom, front)
                plt.figure("Kite")
                lines, sc, txt = plot(sam, 0.0; zoom, front)
                plt.show(block=false)
                sleep(1)
                @test !isnothing(lines)
                @test length(lines) ≥ 1  # Should have at least one line
                @test !isnothing(sc)     # Should have scatter points
                @test !isnothing(txt)    # Should have time text
            end
            plot_(false, false)
            plot_(false, true)
            plot_(true, false)
            plot_(true, true)
        end
    end

    function init(elevation, azimuth, heading)
        set.elevation = elevation
        set.azimuth = azimuth
        set.heading = heading
        init!(sam)
        ss = SysState(sam)
        @test sam.sys_struct.wings[1].elevation ≈ deg2rad(set.elevation) atol=1e-2
        @test sam.sys_struct.wings[1].azimuth ≈ deg2rad(set.azimuth) atol=1e-2
        @test sam.sys_struct.wings[1].heading ≈ deg2rad(set.heading) atol=1e-2
        @test ss.elevation ≈ deg2rad(set.elevation) atol=1e-2
        @test ss.azimuth ≈ deg2rad(set.azimuth) atol=1e-2
        @test ss.heading ≈ deg2rad(set.heading) atol=1e-2
    end
    
    @testset "Initialization" begin
        @test sam isa SymbolicAWEModel
        @test sam isa AbstractKiteModel
        @test simple_sam isa SymbolicAWEModel
        @test simple_sam isa AbstractKiteModel
        @test tether_sam isa SymbolicAWEModel
        @test tether_sam isa AbstractKiteModel
        
        init!(sam; prn=true)
        init_time = @elapsed init!(sam; prn=true)
        @test init_time < 0.3

        init(ones(3)...)
        init(zeros(3)...)
        test_plot(sam)
    end

    @testset "Tether properties" begin
        reset!(set)
        set.sample_freq = 600
        set.abs_tol = 1e-6
        set.rel_tol = 1e-6
        set.segments = 1
        one_seg_sam = SymbolicAWEModel(set, "ram")
        init!(one_seg_sam)
        one_seg_tether_sam = SymbolicAWEModel(set, "tether")
        init!(one_seg_tether_sam)

        axial_stiffness, axial_damping = 
            SymbolicAWEModels.calc_spring_props(one_seg_sam, one_seg_tether_sam)
        segments = one_seg_sam.sys_struct.segments
        tethers = one_seg_sam.sys_struct.tethers
        segments = [segments[tether.segment_idxs[1]] for tether in tethers]
        real_axial_stiffness = [segment.axial_stiffness for segment in segments]
        real_axial_damping = [segment.axial_damping for segment in segments]
        @test isapprox(real_axial_stiffness, axial_stiffness; rtol=0.02)
        @test isapprox(real_axial_damping, axial_damping; rtol=0.2)

        println("\n--- Tether Spring Properties ---")
        # Print table headers
        @printf "%-8s | %-15s %-15s %-10s | %-15s %-15s %-10s\n" "Tether" "Calc. Stiffness" "Real Stiffness" "Error (%)" "Calc. Damping" "Real Damping" "Error (%)"
        # Print separator line
        println(repeat("-", 100))
        for i in 1:4
            # Calculate relative errors in percent
            stiffness_err = 100 * abs(axial_stiffness[i] - real_axial_stiffness[i]) / real_axial_stiffness[i]
            damping_err   = 100 * abs(axial_damping[i] - real_axial_damping[i]) / real_axial_damping[i]
            # Print data rows
            @printf "%-8d | %-15.2f %-15.2f %-10.2f | %-15.2f %-15.2f %-10.2f\n" i axial_stiffness[i] real_axial_stiffness[i] stiffness_err axial_damping[i] real_axial_damping[i] damping_err
        end
        println()
    end

    @testset "Oscillating simulation" begin
        function test_for_peak_at_steering_freq(sam, steering_freq)
            dt = 0.01
            sl, _ = sim_oscillate!(sam; total_time=5.0, steering_freq, dt)
            @test sl.syslog.elevation[begin] ≈ deg2rad(set.elevation) atol=1e-2
            @test sl.syslog.azimuth[begin] ≈ deg2rad(set.azimuth) atol=1e-2
            @test sl.syslog.heading[begin] ≈ deg2rad(set.heading) atol=1e-2
            @test isapprox(sl.syslog.time, collect(0.0:dt:5.0-dt))
            ControlPlots.plt.close_figs()
            plt = plot(sam.sys_struct, sl)
            display(plt)
            @test plt isa ControlPlots.PlotX
            savefig(joinpath(
                get_data_path(), "oscillate_$(sam.sys_struct.name)_$steering_freq.png"
            ))

            # --- Cross-Correlation Analysis ---
            # 1. Extract the detrended signal and time vector
            heading_signal = sl.syslog.heading
            signal_detrended = heading_signal .- mean(heading_signal)
            t = sl.syslog.time

            # 2. Generate a perfect reference sine and cosine wave at the steering frequency
            ref_sin = sin.(2 * π * steering_freq .* t)
            ref_cos = cos.(2 * π * steering_freq .* t)

            # 3. Calculate the correlation of your signal with the reference waves
            # This is equivalent to finding the magnitude of the Fourier series coefficient
            # for this specific frequency.
            corr_sin = dot(signal_detrended, ref_sin)
            corr_cos = dot(signal_detrended, ref_cos)

            # The magnitude of the correlation is proportional to the amplitude of the
            # oscillation at the steering frequency.
            magnitude_at_freq = sqrt(corr_sin^2 + corr_cos^2)

            # 4. To test, we can compare this magnitude to the magnitude at other
            #    nearby frequencies to ensure it's a local peak.
            #    Let's check against a frequency 20% lower and 20% higher.
            freq_lower = steering_freq * 0.5
            ref_sin_lower = sin.(2 * π * freq_lower .* t)
            ref_cos_lower = cos.(2 * π * freq_lower .* t)
            mag_lower = sqrt(dot(signal_detrended, ref_sin_lower)^2 + dot(signal_detrended, ref_cos_lower)^2)

            freq_higher = steering_freq * 1.5
            ref_sin_higher = sin.(2 * π * freq_higher .* t)
            ref_cos_higher = cos.(2 * π * freq_higher .* t)
            mag_higher = sqrt(dot(signal_detrended, ref_sin_higher)^2 + dot(signal_detrended, ref_cos_higher)^2)

            @show magnitude_at_freq, mag_lower, mag_higher
            @test magnitude_at_freq > mag_lower && magnitude_at_freq > mag_higher
        end

        reset!(set)

        init!(sam)
        find_steady_state!(sam)
        test_for_peak_at_steering_freq(sam, 0.5)

        init!(sam)
        find_steady_state!(sam)
        SymbolicAWEModels.copy_to_simple!(sam, tether_sam, simple_sam)
        test_for_peak_at_steering_freq(simple_sam, 0.5)
    end

    @testset "Turning simulation" begin
        function unwrap!(v::AbstractVector, period::Real=2π)
            offset = 0.0
            for i in 2:length(v)
                diff = v[i] - v[i-1]
                if diff > period / 2
                    offset -= period
                elseif diff < -period / 2
                    offset += period
                end
                v[i] += offset
            end
            return v
        end
        function calc_heading(steering_time, steering_magnitude)
            reset!(set)
            init!(sam)
            find_steady_state!(sam)
            dt = 0.05
            sl, _ = sim_turn!(sam; total_time=10.0, steering_time, steering_magnitude, dt)
            unwrap!(sl.syslog.heading)
            @test sl.syslog.heading[begin] ≈ 0.0 atol=1e-1
            return sl.syslog.heading[end]
        end
        default_heading = calc_heading(1.0, 10.0)
        @test default_heading ≈ 1035 atol=10.0
        short_steer_heading = calc_heading(0.5, 10.0)
        soft_steer_heading = calc_heading(1.0, 5.0)
        # make sure less steering results in less final heading
        @test default_heading - short_steer_heading ≈ 93 atol=10.0
        @test default_heading - soft_steer_heading ≈ 150 atol=10.0
        @show default_heading, short_steer_heading, soft_steer_heading
    end

    @testset "Linearize" begin
        old_abs = set.abs_tol
        old_rel = set.rel_tol
        set.abs_tol = 1e-4
        set.rel_tol = 1e-4
        init!(sam)
        init!(simple_sam)

        (; A, B, C, D) = SymbolicAWEModels.linearize!(simple_sam)
        sys = ss(A,B,C,D)
        norm_A = norm(A)
        res = lsim(sys, repeat([-1.0 0.0 -1.0], 2)', [0.0, 0.5])
        println(res.y[:,2])
        @test isapprox(res.y[:,2], 
            [-0.0008037289321365251, 0.0004562826732837309, -0.020711457720341487, 
                       -0.0017333135190197818], rtol=0.1)

        find_steady_state!(sam; dt=3.0, t=10.0)
        (; A, B, C, D) = SymbolicAWEModels.simple_linearize!(sam; tstab=1.0)
        sys = ss(A,B,C,D)
        res = lsim(sys, repeat([-1.0 0.0 -1.0], 2)', [0.0, 0.5])
        println(res.y[:,2])
        @test isapprox(res.y[:,2],
            [0.015575316961016356, -0.0001989661253600774, -0.017933805715950355, 
                       6.679990358160092], rtol=0.1)

        # test that linearization is state-dependent
        next_step!(simple_sam; dt=1.0)
        (; A, B, C, D) = SymbolicAWEModels.linearize!(simple_sam)
        @test !isapprox(norm(A), norm_A; atol=1e-3)

        set.abs_tol = old_abs
        set.rel_tol = old_rel
    end

    @testset "Just a tether, without winch or kite" begin
        set.segments = 20
        dynamics_type = DYNAMIC

        points = Point[]
        segments = Segment[]

        points = push!(points, Point(1, zeros(3), STATIC; wing_idx=0))

        segment_idxs = Int[]
        for i in 1:set.segments
            point_idx = i+1
            pos = [0.0, 0.0, i * set.l_tether / set.segments]
            push!(points, Point(point_idx, pos, dynamics_type; wing_idx=0))
            segment_idx = i
            push!(segments, Segment(segment_idx, set, (point_idx-1, point_idx), BRIDLE))
            push!(segment_idxs, segment_idx)
        end

        transforms = [Transform(1, deg2rad(-80), 0.0, 0.0; 
            base_pos=[0.0, 0.0, 50.0], base_point_idx=points[1].idx, rot_point_idx=points[end].idx)]
        sys_struct = SymbolicAWEModels.SystemStructure("tether", set; points, segments, transforms)

        sam = SymbolicAWEModel(set, sys_struct)
        sys = sam.sys
        init!(sam; remake=false)
        @test isapprox(sam.integrator[sam.sys.pos[:, end]], [8.682408883346524, 0.0, 0.7596123493895988], atol=1e-2)
        for i in 1:100
            next_step!(sam)
        end
        @test sam.integrator[sam.sys.pos[1, end]] > 0.8set.l_tether
        @test isapprox(sam.integrator[sam.sys.pos[2, end]], 0.0, atol=1.0)
        test_plot(sam)
    end
end
nothing
