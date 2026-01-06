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
    println("Moment - Diff:   ", round.(vsm_moment - total_moment_refine, digits=2),
            " (", round.(100 .* (vsm_moment - total_moment_refine) ./
                         (abs.(vsm_moment) .+ 1e-6), digits=1), "%)")

    return (vsm_force, total_force_refine, vsm_moment, total_moment_refine)
end

function validate_forces_aerostruc_mapping(solver, wing, body_aero; label="", reference_point=zeros(3))
    sol = solver.sol
    panels = body_aero.panels
    f_body = sol.f_body_3D
    m_body = sol.m_body_3D
    r_ref = SVector{3}(reference_point)

    start_idx = 1
    for w in body_aero.wings
        w === wing && break
        start_idx += w.n_panels
    end

    idxs = start_idx:(start_idx + wing.n_panels - 1)
    Fsum = zero(SVector{3,Float64})
    Msum = zero(SVector{3,Float64})
    for j in idxs
        Fsum += SVector{3}(f_body[:, j])
        Msum += SVector{3}(m_body[:, j])
    end

    total_force = zeros(3)
    total_moment = zeros(3)
    res_norm_max = 0.0
    res_norm_sumsq = 0.0
    mp_norm_max = 0.0
    mp_norm_sumsq = 0.0
    n_panels = 0

    for local_panel_idx in 1:wing.n_panels
        panel_idx = start_idx + local_panel_idx - 1
        panel = panels[panel_idx]
        Fp = SVector{3}(f_body[:, panel_idx])
        Mp = SVector{3}(m_body[:, panel_idx])  # about reference_point

        F1, F2, F3, F4, nodes = compute_aerostruc_loads(panel, Fp, Mp; reference_point=r_ref)
        Fs = (F1, F2, F3, F4)
        total_force .+= sum(Fs)
        total_moment .+= cross(nodes[1] - r_ref, F1) + cross(nodes[2] - r_ref, F2) +
                         cross(nodes[3] - r_ref, F3) + cross(nodes[4] - r_ref, F4)

        # Per-panel residual moment (diagnostic)
        mapped_moment = cross(nodes[1] - r_ref, F1) + cross(nodes[2] - r_ref, F2) +
                        cross(nodes[3] - r_ref, F3) + cross(nodes[4] - r_ref, F4)
        res = Mp - mapped_moment
        res_norm = norm(res)
        res_norm_max = max(res_norm_max, res_norm)
        res_norm_sumsq += res_norm^2
        mp_norm = norm(Mp)
        mp_norm_max = max(mp_norm_max, mp_norm)
        mp_norm_sumsq += mp_norm^2
        n_panels += 1
    end

    vsm_force = Vector(sol.force)
    vsm_moment = Vector(sol.moment)
    res_norm_rms = n_panels > 0 ? sqrt(res_norm_sumsq / n_panels) : 0.0
    mp_norm_rms = n_panels > 0 ? sqrt(mp_norm_sumsq / n_panels) : 0.0
    res_max_pct = mp_norm_max > 1e-12 ? 100 * res_norm_max / mp_norm_max : 0.0
    res_rms_pct = mp_norm_rms > 1e-12 ? 100 * res_norm_rms / mp_norm_rms : 0.0

    println("\n=== $label (aero→struct mapping) ===")
    println("Force  - Panel sum vs VSM: ", round.(Fsum, digits=2), " vs ", round.(vsm_force, digits=2))
    println("Moment - Panel sum vs VSM: ", round.(Msum, digits=2), " vs ", round.(vsm_moment, digits=2))
    println("Force  - VSM:    ", round.(vsm_force, digits=2))
    println("Force  - Aero→Struct: ", round.(total_force, digits=2))
    println("Force  - Diff:   ", round.(vsm_force - total_force, digits=2),
            " (", round.(100 .* (vsm_force - total_force) ./
                         (abs.(vsm_force) .+ 1e-6), digits=1), "%)")
    println("Moment - VSM:    ", round.(vsm_moment, digits=2))
    println("Moment - Aero→Struct: ", round.(total_moment, digits=2))
    println("Moment - Diff:   ", round.(vsm_moment - total_moment, digits=2),
            " (", round.(100 .* (vsm_moment - total_moment) ./
                         (abs.(vsm_moment) .+ 1e-6), digits=1), "%)")
    println("Per-panel mapped moment residuals (norm): max=", round(res_norm_max, digits=4),
            " (", round(res_max_pct, digits=2), "% of max |Mp|), rms=",
            round(res_norm_rms, digits=4), " (", round(res_rms_pct, digits=2), "% of rms |Mp|)")

    return (total_force, total_moment)
end


function compute_aerostruc_loads(panel, F_panel::SVector{3}, M_panel::SVector{3};
    reference_point::SVector{3}=SVector(0.0, 0.0, 0.0))
    
    # Approximate quad corners from aero center, chord, width, and local axes
    # build panel corner nodes
    c_vec = panel.x_airf
    s_vec = panel.y_airf
    chord = panel.chord
    width = panel.width
    r_ac = SVector{3}(panel.aero_center)
    # Assume aero center at quarter-chord midspan
    r_le_mid = r_ac - 0.25 * chord * c_vec
    r_te_mid = r_ac + 0.75 * chord * c_vec
    half_span = 0.5 * width * s_vec
    le_left = r_le_mid - half_span
    le_right = r_le_mid + half_span
    te_right = r_te_mid + half_span
    te_left = r_te_mid - half_span
    
    nodes = (le_left, le_right, te_right, te_left)
    r_ref = reference_point
    r_ac = SVector{3}(panel.aero_center)

    # Midpoints of LE/TE edges
    r_le_mid = 0.5 * (nodes[1] + nodes[2])
    r_te_mid = 0.5 * (nodes[3] + nodes[4])

    # Relative positions
    r_le_rel = r_le_mid - r_ac
    r_te_rel = r_te_mid - r_ac
    d = r_le_rel - r_te_rel  # chord direction (LE→TE)
    d_norm_sq = dot(d, d)
    if d_norm_sq < 1e-12
        # Degenerate chord: just split equally and spanwise-preserve the torque
        F_le = 0.5 * F_panel
        F_te = 0.5 * F_panel
    else
        # Minimum-norm split that preserves moment about r_ref
        w_le = clamp(-dot(r_te_rel, d) / d_norm_sq, 0.0, 1.0)
        r_weighted = w_le * r_le_rel + (1 - w_le) * r_te_rel
        M_cp = M_panel - cross(r_ac - r_ref, F_panel)
        M_target = M_cp - cross(r_weighted, F_panel)
        ΔF = cross(M_target, d) / d_norm_sq
        F_le = w_le * F_panel + ΔF
        F_te = (1 - w_le) * F_panel - ΔF
    end

    # Spanwise split preserving moment about r_ref
    span_split = function (F::SVector{3}, r_mid::SVector{3}, r_left::SVector{3}, r_right::SVector{3})
        a = cross(r_left - r_ref, F)
        b = cross(r_right - r_ref, F)
        m_target = cross(r_mid - r_ref, F)
        ab = a - b
        denom = dot(ab, ab)
        w = denom < 1e-14 ? 0.5 : clamp(dot(ab, m_target - b) / denom, 0.0, 1.0)
        return (w * F, (1 - w) * F)
    end

    F_le_left, F_le_right = span_split(F_le, r_le_mid, nodes[1], nodes[2])
    # nodes[3] = TE right, nodes[4] = TE left → pass left/right accordingly
    F_te_left, F_te_right = span_split(F_te, r_te_mid, nodes[4], nodes[3])

    return (F_le_left, F_le_right, F_te_right, F_te_left, nodes)
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
# validate_forces(solver, wing; label="No deformation")
validate_forces_aerostruc_mapping(solver, wing, body_aero; label="No deformation")

# Test 2: Small twist deformation
println("\nApplying small twist deformation...")
VortexStepMethod.unrefined_deform!(wing,
    deg2rad.(range(-5, 5, length=wing.n_unrefined_sections)),
    zeros(wing.n_unrefined_sections); smooth=true)
VortexStepMethod.reinit!(body_aero; init_aero=false)
solve!(solver, body_aero)
# validate_forces(solver, wing; label="Twist -5° to +5°")
validate_forces_aerostruc_mapping(solver, wing, body_aero; label="Twist -5° to +5°")

# Test 3: Larger twist
println("\nApplying larger twist deformation...")
VortexStepMethod.unrefined_deform!(wing,
    deg2rad.(range(-10, 10, length=wing.n_unrefined_sections)),
    zeros(wing.n_unrefined_sections); smooth=true)
VortexStepMethod.reinit!(body_aero; init_aero=false)
solve!(solver, body_aero)
# validate_forces(solver, wing; label="Twist -10° to +10°")
validate_forces_aerostruc_mapping(solver, wing, body_aero; label="Twist -10° to +10°")

println("\nValidation complete.")
