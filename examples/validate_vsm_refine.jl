using LinearAlgebra
using VortexStepMethod
using StaticArrays

# Copy of compute_section_le_te_forces from vsm_refine.jl for standalone validation
function compute_section_le_te_forces(x_airf, y_airf, z_airf, chord, width,
                                       cl, cd, cm, alpha_corrected, density, v_a_mag)
    q = 0.5 * density * v_a_mag^2
    L_total = cl * q * chord * width
    D_total = cd * q * chord * width

    x_cp = abs(cl) > 1e-6 ? clamp(0.25 + cm / cl, 0.0, 1.0) : 0.25
    L_TE = L_total * x_cp
    L_LE = L_total * (1.0 - x_cp)
    D_LE = D_total / 2.0
    D_TE = D_total / 2.0

    dir_induced_va = cos(alpha_corrected) * x_airf + sin(alpha_corrected) * z_airf
    dir_induced_va = dir_induced_va / (norm(dir_induced_va) + 1e-12)

    lift_dir = cross(dir_induced_va, y_airf)
    lift_dir = norm(lift_dir) > 1e-12 ? lift_dir / norm(lift_dir) : z_airf

    spanwise = SVector(0.0, 1.0, 0.0)
    drag_dir = cross(spanwise, lift_dir)
    drag_dir = norm(drag_dir) > 1e-12 ? drag_dir / norm(drag_dir) : dir_induced_va

    F_LE = L_LE * lift_dir + D_LE * drag_dir
    F_TE = L_TE * lift_dir + D_TE * drag_dir
    return (F_LE, F_TE)
end

function validate_forces(solver, wing; label="")
    sol = solver.sol
    n_sections = length(wing.unrefined_sections)
    density = solver.density

    # Compute forces using vsm_refine approach
    total_force_refine = zeros(3)
    total_moment_refine = zeros(3)

    for i in 1:n_sections
        section = wing.unrefined_sections[i]
        cl = sol.cl_unrefined_dist[i]
        cd = sol.cd_unrefined_dist[i]
        cm = sol.cm_unrefined_dist[i]
        alpha = sol.alpha_unrefined_dist[i]
        x_airf = SVector{3}(sol.x_airf_unrefined_dist[i])
        y_airf = SVector{3}(sol.y_airf_unrefined_dist[i])
        z_airf = SVector{3}(sol.z_airf_unrefined_dist[i])
        va = SVector{3}(sol.va_unrefined_dist[i])
        chord = sol.chord_unrefined_dist[i]
        width = sol.width_unrefined_dist[i]
        v_a_mag = norm(va)

        F_LE, F_TE = compute_section_le_te_forces(
            x_airf, y_airf, z_airf, chord, width,
            cl, cd, cm, alpha, density, v_a_mag
        )

        total_force_refine .+= F_LE .+ F_TE

        # Moment from LE/TE positions
        pos_LE = section.LE_point
        pos_TE = section.TE_point
        total_moment_refine .+= cross(pos_LE, F_LE) .+ cross(pos_TE, F_TE)
    end

    # VSM forces/moments
    vsm_force = Vector(sol.force)
    vsm_moment = Vector(sol.moment)

    # Print comparison
    println("\n=== $label ===")
    println("Force  - VSM:    ", round.(vsm_force, digits=2))
    println("Force  - Refine: ", round.(total_force_refine, digits=2))
    println("Force  - Diff:   ", round.(vsm_force - total_force_refine, digits=2),
            " (", round.(100 .* (vsm_force - total_force_refine) ./
                         (abs.(vsm_force) .+ 1e-6), digits=1), "%)")
    println("Moment - VSM:    ", round.(vsm_moment, digits=2))
    println("Moment - Refine: ", round.(total_moment_refine, digits=2))
    println("Moment - Diff:   ", round.(vsm_moment - total_moment_refine, digits=2))

    return (vsm_force, total_force_refine, vsm_moment, total_moment_refine)
end

# Main validation
println("Loading V3 kite settings...")
settings = VSMSettings("v3/vsm_settings_reduced_for_coupling.yaml")
settings.wings[1].geometry_file = "data/v3/aero_geometry.yaml"
wing = Wing(settings)
refine!(wing)
body_aero = BodyAerodynamics([wing])
VortexStepMethod.reinit!(body_aero)
solver = Solver(body_aero, settings)
set_va!(body_aero, settings)

# Test 1: No deformation
println("\nSolving with no deformation...")
solve!(solver, body_aero)
validate_forces(solver, wing; label="No deformation")

# Test 2: Small twist deformation
println("\nApplying small twist deformation...")
VortexStepMethod.unrefined_deform!(wing,
    deg2rad.(range(-5, 5, length=wing.n_unrefined_sections)),
    zeros(wing.n_unrefined_sections); smooth=true)
VortexStepMethod.reinit!(body_aero; init_aero=false)
solve!(solver, body_aero)
validate_forces(solver, wing; label="Twist -5° to +5°")

# Test 3: Larger twist
println("\nApplying larger twist deformation...")
VortexStepMethod.unrefined_deform!(wing,
    deg2rad.(range(-10, 10, length=wing.n_unrefined_sections)),
    zeros(wing.n_unrefined_sections); smooth=true)
VortexStepMethod.reinit!(body_aero; init_aero=false)
solve!(solver, body_aero)
validate_forces(solver, wing; label="Twist -10° to +10°")

println("\nValidation complete.")
