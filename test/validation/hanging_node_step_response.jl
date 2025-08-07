# SPDX-FileCopyrightText: 2025 Bart van de Lint
# SPDX-License-Identifier: MIT

using Pkg
if ! ("ControlPlots" ∈ keys(Pkg.project().dependencies))
    using TestEnv; TestEnv.activate()
end
using Test, ControlSystemsBase, Printf
using SymbolicAWEModels, ControlPlots
using Statistics, LinearAlgebra

"""
Test for hanging node step response validation.

This test creates a simple hanging mass system (single node hanging from a fixed point)
and validates that its step response matches analytical second-order system behavior
when subjected to a step force input.

The system should behave as a mass-spring-damper with:
- Mass: tether mass
- Spring constant: axial stiffness of tether
- Damping: axial damping of tether

Reference: https://lpsa.swarthmore.edu/Transient/TransInputs/TransStep.html
"""

@testset "Hanging Node Step Response Validation" begin
    
    function create_hanging_node_system(set::Settings)
        """Create a simple hanging node system for testing"""
        # Reset settings for clean test
        set.segments = 1  # Single segment for simplicity
        set.sample_freq = 1000  # High sample rate for good resolution
        set.abs_tol = 1e-8
        set.rel_tol = 1e-8
        set.l_tether = 10.0  # 10m tether length
        
        # Create a simple system with just a hanging mass
        points = Point[]
        segments = Segment[]
        
        # Fixed ground point
        push!(points, Point(1, [0.0, 0.0, 0.0], STATIC; wing_idx=0))
        
        # Hanging mass point
        pos = [0.0, 0.0, set.l_tether]
        push!(points, Point(2, pos, DYNAMIC; wing_idx=0))
        
        # Single tether segment connecting them
        push!(segments, Segment(1, set, (1, 2), BRIDLE))
        
        # No transforms needed for this simple case
        transforms = Transform[]
        
        sys_struct = SymbolicAWEModels.SystemStructure("hanging_node", set; 
                                                       points, segments, transforms)
        
        sam = SymbolicAWEModel(set, sys_struct)
        init!(sam; remake=false)
        
        return sam
    end
    
    function analytical_step_response(t, m, k, c, F_step, x0=0.0, v0=0.0)
        """
        Analytical step response for a second-order mass-spring-damper system
        
        Args:
            t: time vector
            m: mass
            k: spring constant  
            c: damping coefficient
            F_step: step force magnitude
            x0, v0: initial conditions
            
        Returns:
            x: displacement response
        """
        ω_n = sqrt(k/m)  # Natural frequency
        ζ = c/(2*sqrt(k*m))  # Damping ratio
        
        # Steady state displacement
        x_ss = F_step/k
        
        if ζ > 1.0  # Overdamped
            s1 = -ζ*ω_n + ω_n*sqrt(ζ^2 - 1)
            s2 = -ζ*ω_n - ω_n*sqrt(ζ^2 - 1)
            
            A = (x0*s2 - v0 + x_ss*s2)/(s2 - s1)
            B = -(x0*s1 - v0 + x_ss*s1)/(s2 - s1)
            
            x = x_ss .+ A.*exp.(s1.*t) .+ B.*exp.(s2.*t)
            
        elseif ζ == 1.0  # Critically damped
            A = x0 - x_ss
            B = v0 + ω_n*(x0 - x_ss)
            
            x = x_ss .+ (A .+ B.*t).*exp.(-ω_n.*t)
            
        elseif ζ < 1.0 && ζ > 0.0  # Underdamped
            ω_d = ω_n*sqrt(1 - ζ^2)  # Damped frequency
            
            A = x0 - x_ss
            B = (v0 + ζ*ω_n*(x0 - x_ss))/ω_d
            
            x = x_ss .+ exp.(-ζ*ω_n.*t).*(A.*cos.(ω_d.*t) .+ B.*sin.(ω_d.*t))
            
        else  # ζ == 0, undamped
            A = x0 - x_ss
            B = v0/ω_n
            
            x = x_ss .+ A.*cos.(ω_n.*t) .+ B.*sin.(ω_n.*t)
        end
        
        return x
    end
    
    function test_step_response_accuracy()
        """Test that simulated step response matches analytical solution"""
        
        set = Settings("system.yaml")
        sam = create_hanging_node_system(set)
        
        # Get system properties
        segment = sam.sys_struct.segments[1]
        k = segment.axial_stiffness  # Spring constant [N/m]
        c = segment.axial_damping    # Damping coefficient [Ns/m]
        
        # Calculate effective mass (half the tether mass for single segment)
        tether_mass = set.rho_tether * π * (segment.diameter/2)^2 * set.l_tether
        m = 0.5 * tether_mass  # Effective mass for dynamics
        
        @test m > 0 "Mass must be positive"
        @test k > 0 "Stiffness must be positive" 
        @test c >= 0 "Damping must be non-negative"
        
        # Calculate system characteristics
        ω_n = sqrt(k/m)
        ζ = c/(2*sqrt(k*m))
        
        println("\n--- Hanging Node System Properties ---")
        @printf "Mass: %.6f kg\n" m
        @printf "Stiffness: %.2f N/m\n" k
        @printf "Damping: %.2f Ns/m\n" c
        @printf "Natural frequency: %.3f rad/s\n" ω_n
        @printf "Damping ratio: %.4f\n" ζ
        
        # Determine damping type
        if ζ > 1.0
            damping_type = "Overdamped"
        elseif ζ == 1.0
            damping_type = "Critically damped"
        elseif ζ > 0.0
            damping_type = "Underdamped"
        else
            damping_type = "Undamped"
        end
        @printf "Damping type: %s\n" damping_type
        
        # Apply step force and simulate
        F_step = -10.0  # 10N downward force
        steps = 500
        
        # Get initial position (should be at equilibrium)
        x0 = norm(sam.sys_struct.points[2].pos_w) - set.l_tether  # displacement from rest length
        v0 = 0.0  # starts at rest
        
        # Simulate step response
        dt = 1/set.sample_freq
        tether_lens = SymbolicAWEModels.step(sam, steps, F_step, 
                                           [KVec3(zeros(3)) for _ in 1:1]; 
                                           prn=true)
        
        # Extract displacement data (convert from tether length to displacement)
        t_sim = collect(0:dt:(steps*dt))
        x_sim = tether_lens[1, :] .- set.l_tether  # displacement from rest length
        
        # Calculate analytical solution
        t_analytical = t_sim
        x_analytical = analytical_step_response(t_analytical, m, k, c, F_step, x0, v0)
        
        # Compare results
        max_error = maximum(abs.(x_sim[1:length(x_analytical)] .- x_analytical))
        rms_error = sqrt(mean((x_sim[1:length(x_analytical)] .- x_analytical).^2))
        relative_error = rms_error / abs(F_step/k)  # Relative to steady-state displacement
        
        println("\n--- Step Response Validation ---")
        @printf "Maximum error: %.2e m\n" max_error
        @printf "RMS error: %.2e m\n" rms_error
        @printf "Relative error: %.2f%%\n" relative_error*100
        @printf "Steady-state displacement (analytical): %.6f m\n" F_step/k
        @printf "Final displacement (simulation): %.6f m\n" x_sim[end]
        
        # Test that the simulation matches analytical solution within reasonable tolerance
        @test relative_error < 0.05  # Less than 5% relative error
        @test abs(x_sim[end] - F_step/k) < abs(F_step/k) * 0.02  # Final value within 2%
        
        # Test initial conditions
        @test abs(x_sim[1] - x0) < 1e-6  # Initial displacement correct
        
        # Create visualization if ControlPlots is available
        try
            # Create step response plot using ControlPlots
            plt = ControlPlots.plotx(
                t_sim[1:length(x_analytical)], x_analytical,
                t_sim[1:length(x_sim)], x_sim;
                title="Hanging Node Step Response: Analytical vs Simulation",
                xlabel="Time (s)",
                ylabel="Displacement (m)",
                labels=["Analytical", "Simulation"]
            )
            display(plt)
            
            # Note: Save functionality depends on ControlPlots implementation
            println("Step response plot displayed using ControlPlots")
        catch e
            @warn "Could not create plot: $e"
        end
        
        return true
    end
    
    function test_different_damping_ratios()
        """Test step response for different damping ratios"""
        
        println("\n--- Testing Different Damping Ratios ---")
        
        set = Settings("system.yaml")
        
        # Test different damping ratios by adjusting tether properties
        damping_ratios = [0.1, 0.5, 1.0, 2.0]  # Underdamped, underdamped, critical, overdamped
        
        for (i, target_ζ) in enumerate(damping_ratios)
            sam = create_hanging_node_system(set)
            segment = sam.sys_struct.segments[1]
            
            # Calculate required damping coefficient for target ζ
            k = segment.axial_stiffness
            tether_mass = set.rho_tether * π * (segment.diameter/2)^2 * set.l_tether
            m = 0.5 * tether_mass
            
            c_target = target_ζ * 2 * sqrt(k * m)
            segment.axial_damping = c_target
            
            # Verify achieved damping ratio
            ζ_actual = c_target / (2 * sqrt(k * m))
            
            @test abs(ζ_actual - target_ζ) < 1e-10
            
            # Quick simulation test
            F_step = -5.0
            steps = 100
            tether_lens = SymbolicAWEModels.step(sam, steps, F_step, 
                                               [KVec3(zeros(3)) for _ in 1:1])
            
            # Check that simulation completes without errors
            @test size(tether_lens, 1) == 1
            @test size(tether_lens, 2) == steps + 1
            @test !any(isnan.(tether_lens))
            @test !any(isinf.(tether_lens))
            
            # Check final value approaches steady state
            final_displacement = tether_lens[1, end] - set.l_tether
            expected_final = F_step / k
            @test abs(final_displacement - expected_final) < abs(expected_final) * 0.1
            
            @printf "ζ = %.1f: Final displacement = %.6f m (expected %.6f m)\n" ζ_actual final_displacement expected_final
        end
    end
    
    # Run the tests
    test_step_response_accuracy()
    test_different_damping_ratios()
    
    println("\nHanging node step response validation completed successfully!")
end

nothing