#!/usr/bin/env julia
"""
Test heading calculation for a kite flying in a circular downwind path.

This test verifies that the heading angle is calculated correctly as the kite
completes a full 2π rotation along a circular path.

Usage:
  julia --project=test test/test_heading_calculation.jl           # Run all tests
  julia --project=test test/test_heading_calculation.jl circular  # Run only circular path test
  julia --project=test test/test_heading_calculation.jl special   # Run only special cases
"""

using Test
using LinearAlgebra
using SymbolicAWEModels: sym_normalize, sym_calc_R_t_w, calc_R_v_w, rotate_around_z, calc_heading

# Support selective test execution via command-line args
const test_patterns = isempty(ARGS) ? String[] : ARGS

println("Running heading calculation tests...")
if !isempty(test_patterns)
    println("Filtering tests matching: ", test_patterns)
end

# Helper to check if a test should run
function should_run_test(test_name::String)
    isempty(test_patterns) && return true
    for pattern in test_patterns
        if occursin(lowercase(pattern), lowercase(test_name))
            return true
        end
    end
    return false
end

# Define test-specific helper functions
function calc_heading_new(R_t_w, e_x)
    """New simplified heading calculation"""
    e_x_t = R_t_w' * e_x
    return atan(e_x_t[2], e_x_t[1])
end

function calc_heading_old(R_t_w, R_v_w)
    """Old heading calculation"""
    heading_vec = R_t_w' * R_v_w[:, 1]
    heading = atan(heading_vec[2], heading_vec[1])
    return heading
end

if should_run_test("circular")
@testset "Heading Calculation - Circular Path" begin
    # Kite flying in a circle downwind
    # Circle parameters
    elevation = deg2rad(45)  # 45° elevation
    radius = 100.0           # 100m from origin
    circle_radius = 20.0     # 20m circle radius

    # Center of the circle (downwind, at given elevation)
    center_distance = radius
    center = [center_distance * cos(elevation), 0.0,
              center_distance * sin(elevation)]

    n_points = 36  # Test at 36 points around the circle
    angles = range(0, 2π, length=n_points+1)[1:end-1]

    headings_new = Float64[]
    headings_old = Float64[]
    expected_headings = Float64[]

    for (i, angle) in enumerate(angles)
        # Position along circular path in a vertical plane
        # Circle normal is along y-axis (azimuthal direction)
        wing_pos = center + circle_radius * [cos(angle), 0.0, sin(angle)]

        # Calculate R_t_w (tether frame)
        R_t_w = sym_calc_R_t_w(wing_pos)

        # Kite's body x-axis should point tangent to the circle
        # Tangent direction in world frame
        tangent = [-sin(angle), 0.0, cos(angle)]
        e_x = sym_normalize(tangent)

        # Expected heading: angle from tether x-axis to body x-axis
        # in tether frame xy-plane
        e_x_t = R_t_w' * e_x
        expected_heading = atan(e_x_t[2], e_x_t[1])

        # Calculate heading using new method
        heading_new = calc_heading_new(R_t_w, e_x)

        # Calculate heading using old method for comparison
        R_v_w = calc_R_v_w(wing_pos, e_x)
        heading_old = calc_heading_old(R_t_w, R_v_w)

        push!(headings_new, heading_new)
        push!(headings_old, heading_old)
        push!(expected_headings, expected_heading)

        # Test that new method matches expected
        @test heading_new ≈ expected_heading atol=1e-10

        # Test that new method matches old method
        @test heading_new ≈ heading_old atol=1e-10
    end

    println("\n=== Circular Path Heading Test ===")
    println("Number of test points: $n_points")
    println("Heading range (new): [$(rad2deg(minimum(headings_new)))°,
            $(rad2deg(maximum(headings_new)))°]")
    println("Heading range (old): [$(rad2deg(minimum(headings_old)))°,
            $(rad2deg(maximum(headings_old)))°]")

    # Test that heading varies over a significant range (at least π/2)
    heading_range = maximum(headings_new) - minimum(headings_new)
    @test heading_range > π/2

    # Test continuity: check that large jumps are due to atan2 discontinuity
    large_jumps = 0
    for i in 1:(n_points-1)
        diff = abs(headings_new[i+1] - headings_new[i])
        # Normalize diff to account for 2π periodicity
        diff_normalized = min(diff, 2π - diff)

        if diff_normalized > π/4
            large_jumps += 1
            # If there's a large jump, it should be approximately π or 2π
            # (atan2 discontinuity)
            @test (abs(diff - π) < 0.1) || (abs(diff - 2π) < 0.1)
        end
    end
    println("Number of discontinuities detected: $large_jumps")
    # For a full circle, we expect at most 2-3 discontinuities
    @test large_jumps <= 3

    # Test that absolute d_heading is approximately constant
    # for uniform circular motion (ignoring discontinuities)
    println("\n=== Testing Constant |d_heading| ===")
    abs_d_headings = Float64[]
    for i in 1:(n_points-1)
        diff = headings_new[i+1] - headings_new[i]
        diff_normalized = min(abs(diff), 2π - abs(diff))

        # Only consider smooth segments (no discontinuities)
        if diff_normalized < π/4
            push!(abs_d_headings, abs(diff))
        end
    end

    if !isempty(abs_d_headings)
        mean_abs_d_heading = sum(abs_d_headings) / length(abs_d_headings)
        std_abs_d_heading = sqrt(sum((abs_d_headings .- mean_abs_d_heading).^2) /
                                 length(abs_d_headings))

        println("Mean |d_heading|: $(rad2deg(mean_abs_d_heading))°")
        println("Std |d_heading|: $(rad2deg(std_abs_d_heading))°")

        if mean_abs_d_heading > 1e-6
            coeff_var = std_abs_d_heading / mean_abs_d_heading
            println("Coefficient of variation: $coeff_var")
            # For uniform circular motion, |d_heading| should be very consistent
            @test coeff_var < 0.05  # Less than 5% variation
        else
            println("Mean |d_heading| is negligible - heading nearly constant")
            @test std_abs_d_heading < 0.01  # Very small variations
        end
    end
end
end  # if should_run_test("circular")

if should_run_test("horizontal")
@testset "Heading Calculation - Horizontal Circle" begin
    # Kite flying in a horizontal circle (constant elevation)
    elevation = deg2rad(45)
    radius = 100.0

    n_points = 24
    azimuths = range(0, 2π, length=n_points+1)[1:end-1]

    for azimuth in azimuths
        # Position at given azimuth and elevation
        wing_pos = radius * [
            cos(elevation) * cos(azimuth),
            cos(elevation) * sin(azimuth),
            sin(elevation)
        ]

        R_t_w = sym_calc_R_t_w(wing_pos)

        # Kite pointing tangent to azimuthal circle (in flight direction)
        tangent_azimuth = [-sin(azimuth), cos(azimuth), 0.0]
        e_x = sym_normalize(tangent_azimuth)

        # Calculate headings
        heading_new = calc_heading_new(R_t_w, e_x)
        R_v_w = calc_R_v_w(wing_pos, e_x)
        heading_old = calc_heading_old(R_t_w, R_v_w)

        # New and old methods should agree
        @test heading_new ≈ heading_old atol=1e-10

        # For horizontal circular motion, heading should be close to π/2
        # (pointing in azimuthal direction)
        @test abs(heading_new) > π/4
    end
end
end  # if should_run_test("horizontal")

if should_run_test("special")
@testset "Heading Calculation - Special Cases" begin
    # Test 1: Kite pointing radially outward (heading = 0)
    println("\n=== Special Case: Radial Alignment ===")
    wing_pos = [70.0, 70.0, 70.0]
    R_t_w = sym_calc_R_t_w(wing_pos)
    e_x = R_t_w[:, 1]  # Aligned with tether x-axis (elevation dir)

    heading_new = calc_heading_new(R_t_w, e_x)
    R_v_w = calc_R_v_w(wing_pos, e_x)
    heading_old = calc_heading_old(R_t_w, R_v_w)

    println("Heading (new): $(rad2deg(heading_new))°")
    println("Heading (old): $(rad2deg(heading_old))°")

    @test heading_new ≈ 0.0 atol=1e-10
    @test heading_new ≈ heading_old atol=1e-10

    # Test 2: Kite pointing in azimuthal direction (heading = π/2)
    println("\n=== Special Case: Azimuthal Alignment ===")
    e_x = R_t_w[:, 2]  # Aligned with tether y-axis (azimuthal)

    heading_new = calc_heading_new(R_t_w, e_x)
    R_v_w = calc_R_v_w(wing_pos, e_x)
    heading_old = calc_heading_old(R_t_w, R_v_w)

    println("Heading (new): $(rad2deg(heading_new))°")
    println("Heading (old): $(rad2deg(heading_old))°")

    @test abs(heading_new - π/2) < 1e-10
    @test heading_new ≈ heading_old atol=1e-10

    # Test 3: Kite at various headings by rotating body x-axis
    println("\n=== Special Case: Various Heading Angles ===")
    # Directly construct body x-axes at different headings in the tether frame
    for rot_angle in [0, π/6, π/4, π/3, π/2, 2π/3, π, -π/6, -π/4]
        # Create e_x in tether frame at specified heading
        e_x_tether = [cos(rot_angle), sin(rot_angle), 0.0]
        # Transform to world frame
        e_x_world = R_t_w * e_x_tether

        heading_new = calc_heading_new(R_t_w, e_x_world)
        R_v_w = calc_R_v_w(wing_pos, e_x_world)
        heading_old = calc_heading_old(R_t_w, R_v_w)

        # Check heading matches expected angle (within atan2 periodicity)
        @test heading_new ≈ rot_angle atol=1e-10
        @test heading_new ≈ heading_old atol=1e-10
    end
    println("✓ All heading angles calculated correctly")
end
end  # if should_run_test("special")

println("\n=== All Heading Tests Passed ===\n")

nothing
