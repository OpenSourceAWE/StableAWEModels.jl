# SPDX-FileCopyrightText: 2025 Bart van de Lint
# SPDX-License-Identifier: MPL-2.0

using Test
using SymbolicAWEModels
using Statistics, LinearAlgebra

tmpdir=mktempdir()
mkpath(joinpath(tmpdir, "data"))
old_data_path = get_data_path()
set_data_path(joinpath(tmpdir, "data"))
cp(old_data_path, get_data_path(); force=true)

set = Settings("system.yaml")
sam = SymbolicAWEModel(set, "ram")

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

@testset verbose=true "Simulation" begin
    @testset "Oscillating simulation" begin
        function test_for_peak_at_steering_freq(sam, steering_freq)
            dt = 0.01
            sl, _ = sim_oscillate!(sam; total_time=5.0, steering_freq, dt)
            @test sl.syslog.elevation[begin] ≈ deg2rad(set.elevation) atol=1e-2
            @test sl.syslog.azimuth[begin] ≈ deg2rad(set.azimuth) atol=1e-2
            @test sl.syslog.heading[begin] ≈ deg2rad(set.heading) atol=1e-2
            @test isapprox(sl.syslog.time, collect(dt:dt:5.0))

            # --- Cross-Correlation Analysis with Linear Offset Removal (first/last only) ---
            heading_signal = sl.syslog.heading
            t = sl.syslog.time
            # Compute linear trend using endpoints
            trend = range(heading_signal[1], heading_signal[end], length=length(t))
            signal_detrended = heading_signal .- trend
            # Reference sine and cosine at steering frequency
            ref_sin = sin.(2 * π * steering_freq .* t)
            ref_cos = cos.(2 * π * steering_freq .* t)
            # Correlation at steering frequency
            corr_sin = dot(signal_detrended, ref_sin)
            corr_cos = dot(signal_detrended, ref_cos)
            magnitude_at_freq = sqrt(corr_sin^2 + corr_cos^2)
            # Compare to frequencies 50% lower and 50% higher
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
        @test default_heading ≈ 850 atol=50.0
        short_steer_heading = calc_heading(0.5, 10.0)
        soft_steer_heading = calc_heading(1.0, 5.0)
        # make sure less steering results in less final heading
        @test default_heading - short_steer_heading > 100
        @test default_heading - soft_steer_heading > 100
        @show default_heading, short_steer_heading, soft_steer_heading
    end

    @testset "Reposition simulation" begin
        reset!(set)
        init!(sam)
        find_steady_state!(sam)
        target_elevation = deg2rad(50)
        target_azimuth = deg2rad(10)
        target_heading = deg2rad(10)
        # Run the simulation with repositioning
        lg = SymbolicAWEModels.sim_reposition!(
            sam;
            total_time=5.0,
            reposition_interval_s=1/set.sample_freq,
            target_elevation,
            target_azimuth,
            target_heading,
            prn=false
        )
        @test lg.syslog.heading[end] ≈ target_heading atol=0.05
        @test lg.syslog.elevation[end] ≈ target_elevation atol=0.05
        @test lg.syslog.azimuth[end] ≈ target_azimuth atol=0.05
    end
end

set_data_path(old_data_path)
