using GLMakie
using SymbolicAWEModels, VortexStepMethod
using YAML


# --- Analytical solution functions for simple pulley ---
"""
    pulley_equilibrium_position(anchor1, anchor2, point_mass, rest_length_per_segment, line_diameter_mm, set)

Analytical equilibrium position for a mass hanging from a pulley system.
Uses proper pulley physics with straight segments and force balance.
Each side is a straight segment; line is linear-elastic.
"""
function pulley_equilibrium_position(anchor1, anchor2, point_mass, rest_length_per_segment, line_diameter_mm, set)
    # Calculate system parameters
    cross_section_area = π * ((1e-3)*line_diameter_mm/2)^2
    span = sqrt(sum((anchor2 .- anchor1).^2))
    k_spring = set.e_tether * cross_section_area / rest_length_per_segment  # CORRECTED: k = EA/L0
    
    # Handle rope mass (crude approximation: add full rope mass to point load)
    total_rope_length = 2 * rest_length_per_segment
    mass_line = set.rho_tether * cross_section_area * total_rope_length
    mass_system = point_mass + mass_line
    total_weight = mass_system * set.g_earth
    
    # Geometric parameters
    midpoint_x = (anchor1[1] + anchor2[1]) / 2
    midpoint_y = (anchor1[2] + anchor2[2]) / 2
    midpoint_z = (anchor1[3] + anchor2[3]) / 2
    
    println("Pulley System Parameters:")
    println("   Point mass: ", point_mass, " kg")
    println("   Rest length per segment: ", rest_length_per_segment, " m")
    println("   Total rope length: ", total_rope_length, " m")
    println("   Line diameter: ", line_diameter_mm, " mm")
    println("   E modulus: ", set.e_tether, " Pa")
    println("   Density: ", set.rho_tether, " kg/m^3")
    println("   Gravity: ", set.g_earth, " m/s^2")
    println("   Cross section area: ", round(cross_section_area, digits=6), " m^2")
    println("   Spring constant per segment: ", round(k_spring, digits=2), " N/m")
    println("   Mass system: ", round(mass_system, digits=4), " kg")
    println("   Anchor span: ", round(span, digits=2), " m")
    println("   Total weight: ", round(total_weight, digits=2), " N")
    
    # Solve for sag s >= 0 using proper pulley physics
    # Equilibrium equation: 2 * T(s) * sin(θ) = W
    # where T(s) = k * (ℓ(s) - L0), sin(θ) = s/ℓ(s), ℓ(s) = sqrt((L/2)^2 + s^2)
    function equilibrium_equation(s)
        if s < 0; return Inf; end
        segment_length = hypot(span/2, s)  # ℓ(s) = sqrt((L/2)^2 + s^2)
        extension = segment_length - rest_length_per_segment
        tension = k_spring * extension
        sin_theta = s / segment_length
        vertical_force = 2 * tension * sin_theta
        return vertical_force - total_weight
    end
    
    # Bracket the solution
    s_lo = 0.0
    s_hi = max(rest_length_per_segment, span) + total_weight / max(k_spring, 1e-9)
    
    # Ensure we have a sign change
    f_lo = equilibrium_equation(s_lo)
    f_hi = equilibrium_equation(s_hi)
    
    if f_lo > 0
        error("Rest length per segment appears too short for equilibrium (try larger rest_length_per_segment or smaller load).")
    end
    
    while f_hi < 0
        s_hi *= 2
        f_hi = equilibrium_equation(s_hi)
        if s_hi > 1e6 * max(span, rest_length_per_segment)
            error("Could not bracket a root - system may be unstable.")
        end
    end
    
    # Bisection method to solve for sag
    for _ in 1:200
        s_mid = 0.5 * (s_lo + s_hi)
        f_mid = equilibrium_equation(s_mid)
        
        if abs(f_mid) < 1e-10 || (s_hi - s_lo) < 1e-10
            equilibrium_z = midpoint_z - s_mid
            println("   Equilibrium sag: ", round(s_mid, digits=4), " m")
            return [midpoint_x, midpoint_y, equilibrium_z]
        end
        
        if sign(f_mid) == sign(f_lo)
            s_lo, f_lo = s_mid, f_mid
        else
            s_hi, f_hi = s_mid, f_mid
        end
    end
    
    error("Bisection did not converge - check system parameters.")
end

println("\n\nSimple Pulley Example\n", "="^40)
set = Settings("base/system.yaml")  # Correct way to load settings
set.v_wind = 0
set.l_tether = 5.0 # set l_tether as it affects the plot size
dynamics_type = DYNAMIC
set.sample_freq = 5  # Increase to 100 Hz for better visualization (dt = 0.01s)

# System parameters
point_mass = 2.0  # Mass of the hanging point in kg
rest_length_per_segment = 3.5  # Rest length per segment in meters
line_diameter_mm = 5.0  # Diameter of the segments in mm

# pulley point was placed in the middle, as the side-to-side movement converges very slowly
points = SymbolicAWEModels.Point[]
push!(points, SymbolicAWEModels.Point(1, [0.0, 0.0, 5.0], STATIC))
push!(points, SymbolicAWEModels.Point(2, [5.0, 0.0, 5.0], STATIC))
push!(points, SymbolicAWEModels.Point(3, [2.5, 0.0, 1], DYNAMIC; mass=point_mass))

segments = Segment[]
push!(segments, Segment(1, set, (3,1), BRIDLE,; l0=rest_length_per_segment, compression_frac=0.01, diameter_mm=line_diameter_mm))  
push!(segments, Segment(2, set, (3,2), BRIDLE,; l0=rest_length_per_segment, compression_frac=0.01, diameter_mm=line_diameter_mm))

pulleys = Pulley[]
push!(pulleys, Pulley(1, (1,2), DYNAMIC))

transforms = [Transform(1, -deg2rad(0.0), 0.0, 0.0; base_pos=[0.0, 0.0, 5.0], base_point_idx=1, rot_point_idx=2)]
sys_struct = SymbolicAWEModels.SystemStructure("pulley", set; points, segments, pulleys, transforms)

# Plot initial state using custom 2D Makie plot
fig_initial = Figure(size=(800, 600))
ax_initial = Axis(fig_initial[1, 1], 
                  xlabel="x [m]", ylabel="z [m]",
                  title="Simple Pulley - Initial State",
                  aspect=DataAspect())

# Plot segments
for seg in sys_struct.segments
    p1 = sys_struct.points[seg.point_idxs[1]].pos_w
    p2 = sys_struct.points[seg.point_idxs[2]].pos_w
    lines!(ax_initial, [p1[1], p2[1]], [p1[3], p2[3]], color=:black, linewidth=2)
end

# Plot points
for p in sys_struct.points
    if p.type == STATIC
        scatter!(ax_initial, [p.pos_w[1]], [p.pos_w[3]], color=:blue, markersize=15, label="Static")
    else
        scatter!(ax_initial, [p.pos_w[1]], [p.pos_w[3]], color=:red, markersize=20, label="Dynamic")
    end
end

display(fig_initial)

# Analyze the system, to find optimal damping
include("damping_analysis.jl")
freq, zeta, recommended_damping = analyze_damping_response(sys_struct, set; verbose=true)
println("\n Setting recommended damping to: ", recommended_damping)
set.axial_damping = recommended_damping  # Update settings with recommended damping

# Create the symbolic model
sam = SymbolicAWEModel(set, sys_struct)
init!(sam; remake=false)

# Store states for animation
states = []
push!(states, deepcopy(sam.sys_struct.points))

# Run simulation for 100 time steps
for i in 1:100
    next_step!(sam)
    push!(states, deepcopy(sam.sys_struct.points))
end

# Create animated 2D plot
fig_anim = Figure(size=(800, 600))
ax_anim = Axis(fig_anim[1, 1], 
               xlabel="x [m]", ylabel="z [m]",
               title="Simple Pulley - Animated",
               aspect=DataAspect())

# Get bounds for axis limits
all_x = [p.pos_w[1] for state in states for p in state]
all_z = [p.pos_w[3] for state in states for p in state]
xlims!(ax_anim, minimum(all_x) - 0.5, maximum(all_x) + 0.5)
ylims!(ax_anim, minimum(all_z) - 0.5, maximum(all_z) + 0.5)

# Observable for animation frame
frame = Observable(1)

# Create observables for segments and points
segment_lines = Observable(Point2f[])
static_points = Observable(Point2f[])
dynamic_points = Observable(Point2f[])
trajectory = Observable(Point2f[])

# Update function
function update_frame(f)
    current_state = states[f]
    
    # Update segments
    seg_pts = Point2f[]
    for seg in sam.sys_struct.segments
        p1 = current_state[seg.point_idxs[1]].pos_w
        p2 = current_state[seg.point_idxs[2]].pos_w
        push!(seg_pts, Point2f(p1[1], p1[3]))
        push!(seg_pts, Point2f(p2[1], p2[3]))
    end
    segment_lines[] = seg_pts
    
    # Update points
    static_pts = Point2f[]
    dynamic_pts = Point2f[]
    for p in current_state
        if p.type == STATIC
            push!(static_pts, Point2f(p.pos_w[1], p.pos_w[3]))
        else
            push!(dynamic_pts, Point2f(p.pos_w[1], p.pos_w[3]))
        end
    end
    static_points[] = static_pts
    dynamic_points[] = dynamic_pts
    
    # Update trajectory (trace of dynamic point)
    traj = [Point2f(states[i][3].pos_w[1], states[i][3].pos_w[3]) for i in 1:f]
    trajectory[] = traj
end

# Plot elements
linesegments!(ax_anim, segment_lines, color=:black, linewidth=2)
scatter!(ax_anim, static_points, color=:blue, markersize=15, label="Static")
scatter!(ax_anim, dynamic_points, color=:red, markersize=20, label="Dynamic")
lines!(ax_anim, trajectory, color=(:green, 0.5), linewidth=1, label="Trajectory")
axislegend(ax_anim, position=:rt)

# Initial update
update_frame(1)

# Animation controls
framerate = 30  # fps
dt = 1.0 / set.sample_freq
slowdown = 2  # Make animation slower by this factor
total_frames = length(states)

# Animate
record(fig_anim, "pulley_animation.mp4", 1:total_frames; framerate=framerate) do f
    update_frame(f)
end
println("Animation saved to pulley_animation.mp4")

# Also display the final frame
update_frame(total_frames)
display(fig_anim)

# --- Final comparison with analytical solution ---
println("\n\nSimulation vs Analytical Comparison\n", "="^45)

### Calculate analytical solution
anchor1 = [0.0, 0.0, 5.0]  # First anchor position
anchor2 = [5.0, 0.0, 5.0]  # Second anchor position
equilibrium_pos = pulley_equilibrium_position(anchor1, anchor2, point_mass, rest_length_per_segment, line_diameter_mm, set)

final_pos = sam.sys_struct.points[3].pos_w
final_x, final_y, final_z = final_pos[1], final_pos[2], final_pos[3]

println("Final Results:")
println("  Simulation final z-position: $(round(final_z, digits=4)) m")
println("  Analytical equilibrium z:    $(round(equilibrium_pos[3], digits=4)) m")
println("  Position error (Δz):         $(round(final_z - equilibrium_pos[3], digits=4)) m")

# Calculate relative error
z_error_percent = abs(final_z - equilibrium_pos[3]) / abs(equilibrium_pos[3]) * 100
println("  Relative error:              $(round(z_error_percent, digits=2))%")
