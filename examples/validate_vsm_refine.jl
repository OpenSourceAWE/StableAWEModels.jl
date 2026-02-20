# Copyright (c) 2025 Bart van de Lint, Jelle Poland
# SPDX-License-Identifier: MPL-2.0

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

# Map refined panel forces/moments directly to LE/TE structural points
function validate_forces_panel_mapping(solver, wing, body_aero; label="", reference_point=zeros(3))
    sol = solver.sol
    panels = body_aero.panels
    f_body = sol.f_body_3D
    m_body = sol.m_body_3D
    r_ref = SVector{3}(reference_point)

    # Find panel index offset for this wing
    start_idx = 1
    for w in body_aero.wings
        w === wing && break
        start_idx += w.n_panels
    end

    idxs = start_idx:(start_idx + wing.n_panels - 1)
    # Diagnostic: check whether summed panel loads equal the reported totals (about r_ref)
    Fsum = zero(SVector{3,Float64})
    Msum = zero(SVector{3,Float64})
    for j in idxs
        Fsum += SVector{3}(f_body[:, j])
        Msum += SVector{3}(m_body[:, j])
    end

    total_force = zeros(3)
    total_moment = zeros(3)

    for local_panel_idx in 1:wing.n_panels
        panel_idx = start_idx + local_panel_idx - 1
        panel = panels[panel_idx]
        Fp = SVector{3}(f_body[:, panel_idx])
        Mp = SVector{3}(m_body[:, panel_idx])  # moment about reference_point

        section_idx = wing.refined_panel_mapping[local_panel_idx]
        section = wing.unrefined_sections[section_idx]
        r_le = SVector{3}(section.LE_point)
        r_te = SVector{3}(section.TE_point)

        F_le, F_te = compute_panel_le_te_forces(panel, Fp, Mp;
            le_point=r_le, te_point=r_te, reference_point=r_ref)

        total_force .+= F_le .+ F_te
        total_moment .+= cross(r_le - r_ref, F_le) .+ cross(r_te - r_ref, F_te)
    end

    vsm_force = Vector(sol.force)
    vsm_moment = Vector(sol.moment)

    println("\n=== $label (panel→LE/TE mapping) ===")
    println("Force  - VSM:    ", round.(vsm_force, digits=2))
    println("Force  - Panel map: ", round.(total_force, digits=2))
    println("Force  - Diff:   ", round.(vsm_force - total_force, digits=2),
            " (", round.(100 .* (vsm_force - total_force) ./
                         (abs.(vsm_force) .+ 1e-6), digits=1), "%)")
    println("Moment - VSM:    ", round.(vsm_moment, digits=2))
    println("Moment - Panel map: ", round.(total_moment, digits=2))
    println("Moment - Diff:   ", round.(vsm_moment - total_moment, digits=2),
            " (", round.(100 .* (vsm_moment - total_moment) ./
                         (abs.(vsm_moment) .+ 1e-6), digits=1), "%)")

    return (total_force, total_moment)
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
    end

    vsm_force = Vector(sol.force)
    vsm_moment = Vector(sol.moment)

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

    return (total_force, total_moment)
end


# Per-panel LE/TE force split that preserves the panel moment about its aero center.
function compute_panel_le_te_forces(panel, F_panel::SVector{3}, M_panel::SVector{3};
    le_point::SVector{3}, te_point::SVector{3}, reference_point::SVector{3}=SVector(0.0, 0.0, 0.0))
    r_le = le_point
    r_te = te_point
    r_ac = SVector{3}(panel.aero_center)
    r_ref = reference_point

    # Convert moment to the panel aero center to avoid reference mismatches
    M_cp = M_panel - cross(r_ac - r_ref, F_panel)
    r_le_rel = r_le - r_ac
    r_te_rel = r_te - r_ac

    d = r_le_rel - r_te_rel  # chord direction (LE→TE)
    d_norm_sq = dot(d, d)
    if d_norm_sq < 1e-12  # degenerate: points coincident
        return F_panel / 2, F_panel / 2
    end

    # Weight that places the base resultant at the aero center (0.75/0.25 for quarter-chord)
    w_le = clamp(-dot(r_te_rel, d) / d_norm_sq, 0.0, 1.0)
    r_weighted = w_le * r_le_rel + (1 - w_le) * r_te_rel

    M_target = M_cp - cross(r_weighted, F_panel)
    ΔF = cross(M_target, d) / d_norm_sq

    F_le = w_le * F_panel + ΔF
    F_te = (1 - w_le) * F_panel - ΔF
    return F_le, F_te
end

# --- Alternative mapping: distribute panel force/moment to four panel corner nodes ---
function line_intersect_clamped(p1::SVector{3}, p2::SVector{3}, p3::SVector{3}, p4::SVector{3})
    p13 = p1 - p3
    p43 = p4 - p3
    p21 = p2 - p1
    d1343 = dot(p13, p43)
    d4321 = dot(p43, p21)
    d1321 = dot(p13, p21)
    d4343 = dot(p43, p43)
    d2121 = dot(p21, p21)
    denom = d2121 * d4343 - d4321 * d4321
    if abs(denom) < 1e-12
        return 0.5 * (p1 + p2)  # fallback: midpoint if nearly parallel
    end
    numer = d1343 * d4321 - d1321 * d4343
    mua = numer / denom
    mua = clamp(mua, 0.0, 1.0)  # clamp to segment p1-p2
    return p1 + mua * p21
end

function force2nodes(F::SVector{3}, Fpoint::SVector{3}, nodes::NTuple{4,SVector{3}}, tangential::SVector{3})
    P1 = line_intersect_clamped(nodes[1], nodes[2], Fpoint, Fpoint + tangential)
    d1 = Fpoint - P1
    M1 = cross(d1, F)

    P2 = line_intersect_clamped(nodes[4], nodes[3], Fpoint, Fpoint + tangential)
    d2 = P2 - P1
    denom = dot(d2, d2)
    denom < 1e-12 && (denom = 1e-12)
    Fp2 = cross(M1, d2) / denom
    Fp1 = F - Fp2

    M3 = cross(P1 - nodes[1], Fp1)
    d3 = nodes[2] - nodes[1]
    denom3 = dot(d3, d3)
    denom3 < 1e-12 && (denom3 = 1e-12)
    F3 = cross(M3, d3) / denom3
    node1 = Fp1 - F3
    node2 = F3

    M4 = cross(P2 - nodes[3], Fp2)
    d4 = nodes[4] - nodes[3]
    denom4 = dot(d4, d4)
    denom4 < 1e-12 && (denom4 = 1e-12)
    F4 = cross(M4, d4) / denom4
    node4 = F4
    node3 = Fp2 - F4

    return (node1, node2, node3, node4)
end

function moment2nodes(M::SVector{3}, Mpoint::SVector{3}, nodes::NTuple{4,SVector{3}}, tangential::SVector{3}, lever::Float64)
    d = lever * tangential
    if norm(d) < 1e-12 || norm(M) < 1e-12
        return (zero(M), zero(M), zero(M), zero(M))
    end
    dF = cross(M, d)
    norm_dF = norm(dF)
    if norm_dF < 1e-12
        return (zero(M), zero(M), zero(M), zero(M))
    end
    dF /= norm_dF
    Fmag = norm(M) / norm(cross(dF, d))
    F = Fmag * dF

    P1 = Mpoint + d
    Fnode1 = force2nodes(F, P1, nodes, tangential)
    Fnode2 = force2nodes(-F, Mpoint, nodes, tangential)

    return (
        Fnode1[1] + Fnode2[1],
        Fnode1[2] + Fnode2[2],
        Fnode1[3] + Fnode2[3],
        Fnode1[4] + Fnode2[4],
    )
end

function build_panel_nodes(panel)
    # Approximate quad corners from aero center, chord, width, and local axes
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
    return (le_left, le_right, te_right, te_left)
end

function split_force_spanwise(F::SVector{3}, r_mid::SVector{3}, r_left::SVector{3},
                              r_right::SVector{3}, r_ref::SVector{3})
    # Choose weights so that applying F at left/right reproduces the moment of F at r_mid (about r_ref)
    a = cross(r_left - r_ref, F)
    b = cross(r_right - r_ref, F)
    m_target = cross(r_mid - r_ref, F)
    ab = a - b
    denom = dot(ab, ab)
    w = if denom < 1e-14
        0.5
    else
        clamp(dot(ab, m_target - b) / denom, 0.0, 1.0)
    end
    F_left = w * F
    F_right = (1 - w) * F
    return F_left, F_right
end

function compute_aerostruc_loads(panel, F_panel::SVector{3}, M_panel::SVector{3};
    reference_point::SVector{3}=SVector(0.0, 0.0, 0.0))
    nodes = build_panel_nodes(panel)
    r_ref = reference_point
    r_le_mid = 0.5 * (nodes[1] + nodes[2])
    r_te_mid = 0.5 * (nodes[3] + nodes[4])

    F_le, F_te = compute_panel_le_te_forces(panel, F_panel, M_panel;
        le_point=r_le_mid, te_point=r_te_mid, reference_point=r_ref)

    # Distribute LE/TE forces along the span while preserving the moment about r_ref
    F_le_left, F_le_right = split_force_spanwise(F_le, r_le_mid, nodes[1], nodes[2], r_ref)
    # nodes[3] = TE right, nodes[4] = TE left → pass left/right accordingly
    F_te_left, F_te_right = split_force_spanwise(F_te, r_te_mid, nodes[4], nodes[3], r_ref)

    return (
        F_le_left,
        F_le_right,
        F_te_right,
        F_te_left,
        nodes,
    )
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
# validate_forces_panel_mapping(solver, wing, body_aero; label="No deformation")
validate_forces_aerostruc_mapping(solver, wing, body_aero; label="No deformation")

# Test 2: Small twist deformation
println("\nApplying small twist deformation...")
VortexStepMethod.unrefined_deform!(wing,
    deg2rad.(range(-5, 5, length=wing.n_unrefined_sections)),
    zeros(wing.n_unrefined_sections); smooth=true)
VortexStepMethod.reinit!(body_aero; init_aero=false)
solve!(solver, body_aero)
# validate_forces(solver, wing; label="Twist -5° to +5°")
# validate_forces_panel_mapping(solver, wing, body_aero; label="Twist -5° to +5°")
validate_forces_aerostruc_mapping(solver, wing, body_aero; label="Twist -5° to +5°")

# Test 3: Larger twist
println("\nApplying larger twist deformation...")
VortexStepMethod.unrefined_deform!(wing,
    deg2rad.(range(-10, 10, length=wing.n_unrefined_sections)),
    zeros(wing.n_unrefined_sections); smooth=true)
VortexStepMethod.reinit!(body_aero; init_aero=false)
solve!(solver, body_aero)
# validate_forces(solver, wing; label="Twist -10° to +10°")
# validate_forces_panel_mapping(solver, wing, body_aero; label="Twist -10° to +10°")
validate_forces_aerostruc_mapping(solver, wing, body_aero; label="Twist -10° to +10°")

println("\nValidation complete.")
