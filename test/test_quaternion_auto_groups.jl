# Copyright (c) 2025 Bart van de Lint, Jelle Poland
# SPDX-License-Identifier: MPL-2.0

# Test auto-creation of groups for QUATERNION wings
using SymbolicAWEModels
using Test
using LinearAlgebra

@testset "QUATERNION wing auto-group creation" begin
    # Copy 2plate_kite data to temp directory
    pkg_root = dirname(@__DIR__)
    src_data_path = joinpath(
        pkg_root, "data", "2plate_kite"
    )
    tmpdir = mktempdir()
    data_path = joinpath(tmpdir, "2plate_kite")
    cp(src_data_path, data_path; force=true)
    set_data_path(data_path)
    set = Settings("system.yaml")

    # Load with REFINE (should have 0 groups)
    println("\n=== Testing REFINE wing (no auto-groups) ===")
    refine_yaml = joinpath(
        data_path, "refine_struc_geometry.yaml"
    )
    vsm_set_path = joinpath(
        data_path, "vsm_settings.yaml"
    )
    vsm_set = SymbolicAWEModels.VortexStepMethod.VSMSettings(
        vsm_set_path; data_prefix=false
    )
    sys_refine = load_sys_struct_from_yaml(
        refine_yaml;
        system_name="2plate_refine", set, vsm_set
    )

    @test length(sys_refine.wings) == 1
    @test sys_refine.wings[1].wing_type == SymbolicAWEModels.REFINE
    @test length(sys_refine.groups) == 0
    @test length(sys_refine.wings[1].group_idxs) == 0
    println("✓ REFINE wing: $(length(sys_refine.groups)) groups")

    # Now test manually creating a QUATERNION wing with WING points
    println("\n=== Testing QUATERNION wing auto-group creation ===")

    # Get WING points from REFINE system
    wing_points = [p for p in sys_refine.points if p.type == SymbolicAWEModels.WING]
    println("Found $(length(wing_points)) WING points")
    @test length(wing_points) == 6  # 3 LE/TE pairs

    # Create a QUATERNION wing with these points
    vsm_wing = SymbolicAWEModels.Wing(
        set, vsm_set; prn=false
    )
    vsm_aero = SymbolicAWEModels.BodyAerodynamics([vsm_wing])
    vsm_solver = SymbolicAWEModels.Solver(vsm_aero;
                                          solver_type=SymbolicAWEModels.NONLIN,
                                          atol=2e-8, rtol=2e-8)

    # Create wing with QUATERNION type and empty group_idxs
    quat_wing = SymbolicAWEModels.VSMWing(
        SymbolicAWEModels.BaseWing(
            :main_wing, Int16[], I(3),
            zeros(3), ones(3);
            wing_type=SymbolicAWEModels.QUATERNION
        ),
        vsm_aero, vsm_wing, vsm_solver,
        Float64[], Float64[], zeros(0, 0),
        nothing, nothing,
        nothing, nothing, nothing, nothing,
        nothing, nothing, 0.0, 0.0
    )

    # Create SystemStructure (should auto-create groups)
    # collect() converts NamedCollection to Vector
    sys_quat = SymbolicAWEModels.SystemStructure(
        "2plate_quat", set;
        points=collect(sys_refine.points),
        segments=collect(sys_refine.segments),
        pulleys=collect(sys_refine.pulleys),
        tethers=collect(sys_refine.tethers),
        winches=collect(sys_refine.winches),
        wings=[quat_wing],
        transforms=collect(sys_refine.transforms)
    )

    # Verify groups were auto-created
    @test length(sys_quat.groups) == 3  # One group per LE/TE pair
    @test length(sys_quat.wings[1].group_idxs) == 3
    @test sys_quat.wings[1].wing_type == SymbolicAWEModels.QUATERNION

    println("✓ QUATERNION wing: $(length(sys_quat.groups)) groups auto-created")

    # Check that geometry was computed from closest VSM panel
    for (i, group) in enumerate(sys_quat.groups)
        println("  Group $i: le_pos = $(group.le_pos), chord_norm = $(norm(group.chord)), points = $(group.point_idxs)")
        @test !iszero(group.chord)  # Chord should be computed from panel
        @test !iszero(group.y_airf)  # y_airf should be computed from panel
    end
end
