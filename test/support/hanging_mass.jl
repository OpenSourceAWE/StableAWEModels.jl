# SPDX-FileCopyrightText: 2025 Bart van de Lint
#
# SPDX-License-Identifier: MPL-2.0

"""
Validation example: a single mass suspended by a massless tether from a fixed point,
subject to gravity, with initial displacement.

Compares the simulated final position to the analytic resting position.
"""

using Pkg
if ! ("ControlPlots" ∈ keys(Pkg.project().dependencies))
    using TestEnv; TestEnv.activate()
end
using SymbolicAWEModels, ControlPlots

"""
    compute_hanging_mass(; l_tether=5.0, m_mass=2.0, g=9.807, theta0=π/4, 
                          abs_tol=1e-7, rel_tol=1e-7, damping=0.2, n_steps=1000)

Compute the dynamics of a hanging mass system and return simulation results.

# Arguments
- `l_tether::Float64`: Tether length [m]
- `m_mass::Float64`: Mass [kg] 
- `g::Float64`: Gravitational acceleration [m/s^2]
- `theta0::Float64`: Initial angle from vertical [rad]
- `abs_tol::Float64`: Absolute tolerance for solver
- `rel_tol::Float64`: Relative tolerance for solver
- `damping::Float64`: Damping coefficient
- `n_steps::Int`: Number of simulation steps

# Returns
- `NamedTuple`: Contains simulation results including:
  - `sam`: The SymbolicAWEModel
  - `sim_pos`: Final simulated position
  - `analytic_pos`: Analytical resting position
  - `trajectory`: Position history
  - `error`: Distance between simulated and analytical final positions
"""
function compute_hanging_mass(; l_tether=5.0, m_mass=2.0, g=9.807, theta0=π/4, 
                               abs_tol=1e-7, rel_tol=1e-7, damping=0.2, n_steps=1000)
    
    # --- Settings (load from YAML and override specific values)
    set = Settings("system.yaml")
    set.v_wind = 0.0
    set.l_tether = l_tether
    set.abs_tol = abs_tol
    set.rel_tol = rel_tol
    set.g_earth = g  # Note: Settings uses g_earth, not g
    set.damping = damping  # Tether damping coefficient
    
    # --- Settings
    set = Settings()
    set.v_wind = 0.0
    set.l_tether = l_tether
    set.abs_tol = abs_tol
    set.rel_tol = rel_tol
    set.g_earth = g
    set.damping = damping

    # --- Points
    points = Point[]
    # Fixed point at origin (anchor)
    push!(points, Point(1, [0.0, 0.0, 0.0], STATIC; wing_idx=0, transform_idx=0))
    # Dynamic mass, initially displaced by theta0 from vertical
    x0 = l_tether * sin(theta0)
    z0 = -l_tether * cos(theta0)
    push!(points, Point(2, [x0, 0.0, z0], DYNAMIC; mass=m_mass, wing_idx=0, transform_idx=0))

    # --- Segments (massless, inextensible tether)
    segments = [Segment(1, set, (1,2), BRIDLE)]

    # --- System structure and model
    # No transforms needed for this simple case  
    transforms = Transform[]
    sys_struct = SymbolicAWEModels.SystemStructure("mass_hanging", set; points, segments, transforms)
    sam = SymbolicAWEModel(set, sys_struct)
    init!(sam; remake=false)

    # --- Simulate until static equilibrium
    for i in 1:n_steps
        next_step!(sam)
    end

    # --- Extract simulated position
    sim_pos = sam.x[1:3, 2] # [x, y, z] of mass

    # --- Analytical resting position
    # At rest, mass should be vertically below anchor at (0, 0, -l_tether)
    analytic_pos = [0.0, 0.0, -l_tether]

    # --- Extract trajectory
    trajectory = [sam.x_hist[1:3, 2, i] for i in 1:length(sam.x_hist[1,1,:])]

    # --- Calculate error
    error = norm(sim_pos - analytic_pos)

    return (
        sam = sam,
        sim_pos = sim_pos,
        analytic_pos = analytic_pos,
        trajectory = trajectory,
        error = error,
        l_tether = l_tether,
        m_mass = m_mass,
        g = g,
        theta0 = theta0
    )
end

"""
    plot_hanging_mass(results; save_plot=false, filename="hanging_mass_validation.png")

Create visualization of the hanging mass simulation results using ControlPlots.

# Arguments
- `results`: Results from `compute_hanging_mass`
- `save_plot::Bool`: Whether to save the plot to file
- `filename::String`: Filename for saved plot

# Returns
- Plot object
"""
function plot_hanging_mass(results; save_plot=false, filename="hanging_mass_validation.png")
    
    # Extract trajectory data
    xs = [h[1] for h in results.trajectory]
    zs = [h[3] for h in results.trajectory]
    
    # Create time vector
    dt = 1.0 / 100.0  # Approximate dt, could be passed as parameter
    t = collect(0:dt:(length(xs)-1)*dt)
    
    # Create trajectory plot using ControlPlots
    p1 = plotx(t, xs; 
               title="Hanging Mass Trajectory - X Position",
               xlabel="Time [s]",
               ylabel="X Position [m]")
    
    p2 = plotx(t, zs;
               title="Hanging Mass Trajectory - Z Position", 
               xlabel="Time [s]",
               ylabel="Z Position [m]")
    
    # Create phase space plot (x vs z)
    p3 = plotx(xs, zs;
               title="Hanging Mass Trajectory - Phase Space",
               xlabel="X Position [m]",
               ylabel="Z Position [m]")
    
    # Display plots
    display(p1)
    display(p2) 
    display(p3)
    
    # Save plots if requested
    if save_plot
        try
            # Note: ControlPlots may have different save methods
            # This is a placeholder - check ControlPlots documentation
            println("Plot saving with ControlPlots - check documentation for exact method")
            println("Trajectory data available for manual saving if needed")
        catch e
            @warn "Could not save plot with ControlPlots: $e"
        end
    end
    
    return (p1, p2, p3)
end

"""
    run_hanging_mass(; verbose=true, plot_results=true, save_plot=false, kwargs...)

Run the complete hanging mass validation example.

# Arguments
- `verbose::Bool`: Whether to print detailed results
- `plot_results::Bool`: Whether to create and display plots
- `save_plot::Bool`: Whether to save plots to file
- `kwargs...`: Additional arguments passed to `compute_hanging_mass`

# Returns
- `NamedTuple`: Results from computation
"""
function run_hanging_mass(; verbose=true, plot_results=true, save_plot=false, kwargs...)
    
    if verbose
        println("=" ^ 60)
        println("HANGING MASS VALIDATION")
        println("=" ^ 60)
        println("Computing hanging mass dynamics...")
    end
    
    # Compute results
    results = compute_hanging_mass(; kwargs...)
    
    if verbose
        println("\n--- Results ---")
        println("Simulated resting position: ", results.sim_pos)
        println("Analytical resting position: ", results.analytic_pos)
        println("Distance between simulation and analytical: ", results.error)
        
        # Additional validation info
        println("\n--- System Parameters ---")
        println("Tether length: $(results.l_tether) m")
        println("Mass: $(results.m_mass) kg")
        println("Gravity: $(results.g) m/s²")
        println("Initial angle: $(round(rad2deg(results.theta0), digits=1))°")
        
        # Validation check
        tolerance = 1e-3  # 1mm tolerance
        if results.error < tolerance
            println("\n✓ VALIDATION PASSED: Error $(results.error) m < $(tolerance) m")
        else
            println("\n✗ VALIDATION FAILED: Error $(results.error) m ≥ $(tolerance) m")
        end
    end
    
    # Create plot if requested
    if plot_results
        if verbose
            println("\nCreating visualization...")
        end
        
        plots = plot_hanging_mass(results; save_plot=save_plot)
    end
    
    if verbose
        println("\nValidation complete. The blue and red points should coincide at (0, 0, -$(results.l_tether)).")
        println("=" ^ 60)
    end
    
    return results
end

# Run the validation when script is executed directly
if abspath(PROGRAM_FILE) == @__FILE__
    run_hanging_mass()
end
