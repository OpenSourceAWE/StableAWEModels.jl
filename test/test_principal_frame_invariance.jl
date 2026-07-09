# SPDX-FileCopyrightText: 2026 Bart van de Lint
# SPDX-License-Identifier: LGPL-3.0-only

# test_principal_frame_invariance.jl
#
# The principal frame is internal ODE state: ANY orthonormal frame that
# diagonalizes the wing's inertia tensor must produce the exact same
# simulation, because every physical quantity goes through
# R_b_to_w = R_p_to_w * R_b_to_p and R_b_to_p absorbs the frame choice.
# A user hit VSM non-convergence when `principal_frame` picked a
# ~90°-flipped axis assignment for the A1-15 wing (the flip leaked into
# the body frame via the old no-ref-points fallback) — that must NOT
# happen, so this test hunts for leaks: it rebuilds the principal frame
# with alternative valid axis choices (including the reported 90°
# Y-flip) and asserts world-frame invariance of the aero force/moment,
# the tether/bridle forces at the WING attachment points, segment
# tensions, winch force, and a short dynamic run — on the quaternion
# (RIGID_DYNAMICS) 2plate kite with twist surfaces (groups).
#
# The compiled model is reused across frame variants: `R_b_to_p`,
# `inertia_principal` are flat params synced from the live struct, and
# `init!` re-derives the principal ODE state via `init_principal_state!`.

using Pkg
if abspath(PROGRAM_FILE) == abspath(@__FILE__)
    Pkg.activate(@__DIR__)
end

using Test
using SymbolicAWEModels
using SymbolicAWEModels: VortexStepMethod, RIGID_DYNAMICS, WING,
    quaternion_to_rotation_matrix
using KiteUtils
using LinearAlgebra

"""
    apply_principal_frame!(wing, S)

Re-express the wing's principal frame through the axis change `S`
(new-principal → old-principal, a proper signed permutation so the
inertia stays diagonal): permutes `inertia_principal` and updates
`R_p_to_c` and `R_b_to_p` consistently.
"""
function apply_principal_frame!(wing, S)
    @assert S * S' ≈ I(3) atol = 1e-12
    @assert det(S) ≈ 1.0 atol = 1e-12
    inertia_new = S' * Diagonal(collect(wing.inertia_principal)) * S
    offdiag = norm(inertia_new - Diagonal(diag(inertia_new)))
    @assert offdiag < 1e-9 * norm(diag(inertia_new))
    wing.inertia_principal .= diag(inertia_new)
    # Per-axis angular damping also lives in the principal frame.
    wing.damping .= diag(S' * Diagonal(collect(wing.damping)) * S)
    wing.R_p_to_c .= wing.R_p_to_c * S
    wing.R_b_to_p .= S' * wing.R_b_to_p
    return nothing
end

"""
    frame_snapshot(sys_struct)

World-frame observables of the current state: wing pose, aero
force/moment (world), per-point positions and forces (tether/bridle
loads at WING points), segment tensions, twist angles, winch force, and
the principal attitude `R_p_to_w` (expected to CHANGE across variants —
used to prove the alternative frame took effect).
"""
function frame_snapshot(sys_struct)
    wing = sys_struct.wings[1]
    R_b_to_w = quaternion_to_rotation_matrix(Vector(wing.Q_b_to_w))
    R_p_to_w = quaternion_to_rotation_matrix(Vector(wing.Q_p_to_w))
    return (;
        R_b_to_w, R_p_to_w,
        wing_pos_w = copy(wing.pos_w),
        omega_b = copy(wing.ω_b),
        aero_force_w = R_b_to_w * Vector(wing.aero_force_b),
        aero_moment_w = R_b_to_w * Vector(wing.aero_moment_b),
        point_pos_w = [copy(point.pos_w) for point in sys_struct.points],
        point_force_w = [copy(point.force) for point in sys_struct.points],
        segment_force = [segment.force for segment in sys_struct.segments],
        twist = [surface.twist for surface in sys_struct.twist_surfaces],
        winch_force = norm(sys_struct.winches[1].force))
end

"""
    angle_between_deg(a, b)

Angle in degrees between two vectors (0 for zero-length input).
"""
function angle_between_deg(a, b)
    scale = norm(a) * norm(b)
    scale ≈ 0 && return 0.0
    return rad2deg(acos(clamp(dot(a, b) / scale, -1.0, 1.0)))
end

"""
    compare_snapshots(snap, ref; rtol, atol_pos, max_angle_deg)

Assert that a variant snapshot matches the baseline in the world frame.
Forces compare relative to the largest force in the baseline (direction
checked separately via `max_angle_deg`), positions with `atol_pos`.
"""
function compare_snapshots(snap, ref; rtol, atol_pos, max_angle_deg)
    @test norm(snap.R_b_to_w - ref.R_b_to_w) < rtol
    @test snap.wing_pos_w ≈ ref.wing_pos_w atol = atol_pos
    @test norm(snap.omega_b - ref.omega_b) <
        rtol * max(norm(ref.omega_b), 1.0)

    force_scale = max(norm(ref.aero_force_w),
        maximum(norm, ref.point_force_w), 1.0)
    @test norm(snap.aero_force_w - ref.aero_force_w) < rtol * force_scale
    @test norm(snap.aero_moment_w - ref.aero_moment_w) <
        rtol * max(norm(ref.aero_moment_w), force_scale)
    @test angle_between_deg(snap.aero_force_w, ref.aero_force_w) <
        max_angle_deg

    for idx in eachindex(ref.point_pos_w)
        @test norm(snap.point_pos_w[idx] - ref.point_pos_w[idx]) < atol_pos
        @test norm(snap.point_force_w[idx] - ref.point_force_w[idx]) <
            rtol * force_scale
        if norm(ref.point_force_w[idx]) > 1e-2 * force_scale
            @test angle_between_deg(
                snap.point_force_w[idx], ref.point_force_w[idx]) <
                max_angle_deg
        end
    end

    seg_scale = max(maximum(abs, ref.segment_force), 1.0)
    for idx in eachindex(ref.segment_force)
        @test abs(snap.segment_force[idx] - ref.segment_force[idx]) <
            rtol * seg_scale
    end
    for idx in eachindex(ref.twist)
        @test abs(snap.twist[idx] - ref.twist[idx]) < rtol
    end
    @test abs(snap.winch_force - ref.winch_force) <
        rtol * max(ref.winch_force, 1.0)
    return nothing
end

"""
    run_pose(sam, set, pose)

Init at a transform pose + wind speed, take one tiny step so all force
outputs are realised, and return the world-frame snapshot.
"""
function run_pose(sam, set, pose)
    elevation, azimuth, heading, v_wind = pose
    transform = sam.sys_struct.transforms[:main_transform]
    transform.elevation = elevation
    transform.azimuth = azimuth
    transform.heading = heading
    transform.elevation_vel = 0.0
    transform.azimuth_vel = 0.0
    set.v_wind = v_wind
    init!(sam; prn=false)
    next_step!(sam; dt=1e-5, vsm_interval=1)
    return frame_snapshot(sam.sys_struct)
end

"""
    run_dynamic(sam, set, poses; steps=20, dt=0.05)

Init at the first pose and step the model, snapshotting every step.
Catches integration-level divergence (the reported failure grew over
time and only blew up at t≈0.9).
"""
function run_dynamic(sam, set, poses; steps=20, dt=0.05)
    run_pose(sam, set, poses[1])
    return [begin
                next_step!(sam; dt, vsm_interval=1)
                frame_snapshot(sam.sys_struct)
            end
            for _ in 1:steps]
end

@testset "Principal frame invariance" begin
    pkg_root = dirname(@__DIR__)
    src_data_path = joinpath(pkg_root, "data", "2plate_kite")
    tmpdir = mktempdir()
    data_path = joinpath(tmpdir, "2plate_kite")
    cp(src_data_path, data_path; force=true)
    set_data_path(data_path)
    set = Settings("system.yaml")
    # Tight tolerances so frame-choice leaks can't hide in integrator noise.
    set.abs_tol = 1e-7
    set.rel_tol = 1e-7

    vsm_set = VortexStepMethod.VSMSettings(
        joinpath(data_path, "vsm_settings.yaml"); data_prefix=false)
    rigid_yaml = joinpath(data_path, "rigid_structural_geometry.yaml")

    sys = load_sys_struct_from_yaml(rigid_yaml;
        system_name="frame_invariance", set, vsm_set,
        aero_mode=AeroLinearized())
    wing = sys.wings[1]
    @test wing.dynamics_type == RIGID_DYNAMICS
    @test length(sys.twist_surfaces) == 3

    sam = SymbolicAWEModel(set, sys)
    init!(sam; prn=false)

    baseline = (R_p_to_c = copy(wing.R_p_to_c),
                R_b_to_p = copy(wing.R_b_to_p),
                inertia_principal = copy(wing.inertia_principal),
                damping = copy(wing.damping))
    function restore_baseline!()
        wing.R_p_to_c .= baseline.R_p_to_c
        wing.R_b_to_p .= baseline.R_b_to_p
        wing.inertia_principal .= baseline.inertia_principal
        wing.damping .= baseline.damping
        return nothing
    end

    poses = [
        (deg2rad(60), 0.0,          0.0,          15.0),
        (deg2rad(70), deg2rad(12),  0.0,          16.0),
        (deg2rad(52), deg2rad(-15), deg2rad(8),   18.0),
    ]

    ref_poses = [run_pose(sam, set, pose) for pose in poses]
    for snap in ref_poses
        @test norm(snap.aero_force_w) > 1.0
        @test maximum(norm, snap.point_force_w) > 1.0
    end
    ref_dynamic = run_dynamic(sam, set, poses)

    # Signed permutations; the first is the flip from the A1-15 failure.
    variants = [
        ("90° about y (reported flip)", [0.0 0 1; 0 1 0; -1 0 0]),
        ("180° about z", [-1.0 0 0; 0 -1 0; 0 0 1]),
        ("cyclic x→y→z→x", [0.0 0 1; 1 0 0; 0 1 0]),
    ]

    for (variant_name, S) in variants
        @testset "$variant_name" begin
            restore_baseline!()
            apply_principal_frame!(wing, S)
            # R_b_to_p must still map body to a frame that is principal.
            @test wing.R_b_to_p ≈ wing.R_p_to_c' * wing.R_b_to_c atol=1e-12

            snap_poses = [run_pose(sam, set, pose) for pose in poses]
            # Prove the variant took effect: principal attitude rotated by S.
            @test norm(snap_poses[1].R_p_to_w -
                ref_poses[1].R_p_to_w * S) < 1e-6
            @test norm(snap_poses[1].R_p_to_w - ref_poses[1].R_p_to_w) > 0.5

            for (snap, ref) in zip(snap_poses, ref_poses)
                compare_snapshots(snap, ref;
                    rtol=1e-4, atol_pos=1e-6, max_angle_deg=0.01)
            end

            snap_dynamic = run_dynamic(sam, set, poses)
            for (snap, ref) in zip(snap_dynamic, ref_dynamic)
                compare_snapshots(snap, ref;
                    rtol=1e-2, atol_pos=1e-3, max_angle_deg=0.5)
            end
        end
    end

    restore_baseline!()
    rm(tmpdir; recursive=true)
end
nothing
