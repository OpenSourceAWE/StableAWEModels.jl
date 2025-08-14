# SPDX-FileCopyrightText: 2025 Bart van de Lint
# SPDX-License-Identifier: MPL-2.0

"""
Simplified Damping Analysis for SymbolicAWEModels

Quick analysis to determine optimal damping parameters:
1. Measure current system response
2. Recommend improved damping values

Usage: 
  analyze_damping_response(sys_struct, set; verbose=true)
"""

using SymbolicAWEModels, Printf, Statistics

"""
    analyze_damping_response(sys_struct, set; verbose=false, target_zeta=0.05)

Analyze system damping and recommend improvements.

Returns: (frequency_hz, damping_ratio, recommended_damping)
"""
function analyze_damping_response(sys_struct::SystemStructure, set::Settings; 
                                  verbose=false, target_zeta=0.05, perturbation_dir=[1.0, 0.0, 0.0])
    verbose && println("🔧 Damping Analysis: $(sys_struct.name)")
    
    # Current settings
    current_damping = set.damping
    current_c_spring = set.c_spring
    current_a1 = current_damping / current_c_spring
    
    verbose && @printf "📊 Current: damping=%.0f, c_spring=%.0f, a₁=%.4f [s]\n" current_damping current_c_spring current_a1
    
    try
        # Create and initialize model
        sam = SymbolicAWEModel(set, sys_struct)
        init!(sam)
        find_steady_state!(sam; t=25.0, dt=0.1)
        
        # Find dynamic points
        dynamic_indices = [i for (i, p) in enumerate(sys_struct.points) if p.type == DYNAMIC]
        isempty(dynamic_indices) && return 0.5, 0.1, current_damping  # Fallback
        
        # Apply perturbation and simulate
        target_idx = dynamic_indices[1]
        initial_pos = copy(sam.sys_struct.points[target_idx].pos_w)
        
        # Calculate dynamic perturbation based on system scale
        # Use 1% of the characteristic length (distance to nearest connected point)
        target_point = sam.sys_struct.points[target_idx]
        min_distance = Inf
        
        # Find connected segments to determine characteristic length
        for segment in sam.sys_struct.segments
            p1_idx, p2_idx = segment.point_idxs
            if p1_idx == target_idx || p2_idx == target_idx
                other_idx = p1_idx == target_idx ? p2_idx : p1_idx
                other_pos = sam.sys_struct.points[other_idx].pos_w
                distance = norm(target_point.pos_w - other_pos)
                min_distance = min(min_distance, distance)
            end
        end
        
        # Fallback to overall system size if no segments found
        if min_distance == Inf
            all_positions = [p.pos_w for p in sam.sys_struct.points]
            bbox_size = maximum(maximum(pos) - minimum(pos) for pos in zip(all_positions...))
            min_distance = bbox_size * 0.1  # 10% of bounding box
        end
        
        # Use 1% of characteristic length, with reasonable bounds
        perturbation_magnitude = clamp(0.01 * min_distance, 0.005, 0.1)  # 5mm to 10cm
        
        # Define perturbation direction (could be made configurable)
        perturbation_vector = perturbation_dir .* perturbation_magnitude
        perturbation_unit = perturbation_vector ./ norm(perturbation_vector)
        
        # Clear velocities before perturbation
        for p in sam.sys_struct.points
            if p.type == DYNAMIC
                p.vel_w .= (0.0, 0.0, 0.0)
            end
        end
        
        sam.sys_struct.points[target_idx].pos_w .+= perturbation_vector
        
        verbose && @printf "  Applied %.1fmm perturbation in [%.1f,%.1f,%.1f] direction (%.1f%% of %.2fm characteristic length)\n" (perturbation_magnitude*1000) perturbation_dir... (perturbation_magnitude/min_distance*100) min_distance
        
        # Record response
        dt = 1.0 / set.sample_freq
        n_steps = Int(8.0 / dt)  # 8 seconds
        response = zeros(n_steps)
        response_3d = zeros(3, n_steps)  # Track all 3 components for fallback
        
        for step in 1:n_steps
            next_step!(sam)
            current_pos = sam.sys_struct.points[target_idx].pos_w
            displacement = current_pos - initial_pos
            response[step] = dot(displacement, perturbation_unit)  # Project onto kick direction
            response_3d[:, step] = displacement  # Store full 3D for fallback
        end
        
        # Check if kick-direction response is meaningful, otherwise use most energetic component
        if var(response) < 1e-10  # Very low variance suggests poor alignment
            variances = [var(response_3d[i, :]) for i in 1:3]
            best_component = argmax(variances)
            response = response_3d[best_component, :]
            verbose && @printf "  Switched to component %d (highest variance: %.2e)\n" best_component variances[best_component]
        end
        
        # Analyze response
        freq, zeta = analyze_response(response, dt)
        
        # Recommend new damping
        if zeta > 0.001  # Valid measurement
            delta_a1 = 2*(target_zeta - zeta) / (2π * freq)
            new_a1 = max(current_a1 + delta_a1, 0.001)  # Minimum damping
            recommended_damping = new_a1 * current_c_spring
            
            # Clamp the jump to avoid absurd changes
            scale = clamp(recommended_damping / current_damping, 0.25, 4.0)
            recommended_damping = current_damping * scale
        else
            recommended_damping = current_damping  # Keep current if measurement failed
        end
        
        if verbose
            @printf "🎯 Measured: freq=%.3f Hz, ζ=%.3f (%.1f%%)\n" freq zeta (zeta*100)
            @printf "💡 Recommended damping: %.0f (%.2fx current)\n" recommended_damping (recommended_damping/current_damping)
        end
        
        return freq, zeta, recommended_damping
        
    catch e
        verbose && println("❌ Analysis failed, using current settings")
        return 0.5, 0.1, current_damping
    end
end

"""
    analyze_response(signal, dt)

Extract frequency and damping ratio from response signal.
Returns: (frequency_hz, damping_ratio)
"""
function analyze_response(signal, dt)
    # Remove DC offset
    signal = signal .- mean(signal)
    
    # Find peaks for frequency estimation
    peaks_idx = Int[]
    peaks_val = Float64[]
    
    # Use relative threshold instead of absolute
    thr = 0.01 * maximum(abs, signal)  # 1% of max signal
    
    for i in 2:length(signal)-1
        if signal[i] > signal[i-1] && signal[i] > signal[i+1] && signal[i] > thr
            if isempty(peaks_idx) || (i - peaks_idx[end]) >= 5  # Min distance
                push!(peaks_idx, i)
                push!(peaks_val, abs(signal[i]))
            end
        end
    end
    
    if length(peaks_idx) >= 3
        # Estimate frequency from peak spacing
        period_samples = mean(diff(peaks_idx))  # Full cycle (not 2.0*)
        freq = 1.0 / (period_samples * dt)
        
        # Estimate damping from peak decay
        if length(peaks_val) >= 2
            log_decrement = mean(log.(peaks_val[1:end-1] ./ peaks_val[2:end]))
            zeta = log_decrement / sqrt((2π)^2 + log_decrement^2)
        else
            println("\n --> Warning: Not enough peaks for damping estimation, using fallback; default value zeta=0.05")
            zeta = 0.05  # Default
        end
    else
        # Heavily damped - estimate from settling
        envelope = abs.(signal)
        settling_idx = findlast(envelope .> 0.01 * maximum(envelope))
        
        if settling_idx !== nothing && settling_idx > 10
            settling_time = settling_idx * dt
            freq = 1.0 / (4 * settling_time)  # Rough estimate
            zeta = 0.7  # Overdamped estimate
        else
            println("\n --> Warning: Settling time analysis failed, using fallback; default value freq=0.5, zeta=0.5")
            freq, zeta = 0.5, 0.5  # Fallback
        end
    end
    
    return freq, max(zeta, 0.001)  # Ensure positive damping
end

"""
    recommend_rayleigh_damping(freq, target_zeta_low=0.02, target_zeta_high=0.08)

Simple Rayleigh damping recommendation for two-frequency targeting.
Returns: (a0, a1) coefficients
"""
function recommend_rayleigh_damping(freq, target_zeta_low=0.02, target_zeta_high=0.08)
    ω1 = 2π * freq
    ω2 = 2π * min(10 * freq, 20.0)  # Higher frequency target
    
    # Solve: ζ = 0.5*(a₀/ω + a₁*ω)
    A = [1/(2ω1) ω1/2; 1/(2ω2) ω2/2]
    a0, a1 = A \ [target_zeta_low; target_zeta_high]
    
    return max(a0, 0.0), max(a1, 0.0)  # Ensure non-negative
end
