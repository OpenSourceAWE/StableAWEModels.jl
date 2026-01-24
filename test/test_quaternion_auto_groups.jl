# Test auto-creation of groups for QUATERNION wings
using SymbolicAWEModels
using Test

@testset "QUATERNION wing auto-group creation" begin
    set_data_path("data/v3")
    set = Settings("system.yaml")

    # Load with REFINE (should have 0 groups)
    println("\n=== Testing REFINE wing (no auto-groups) ===")
    yaml_path = "data/v3/struc_geometry.yaml"
    sys_refine = load_sys_struct_from_yaml(yaml_path; system_name="v3_refine", set=set)

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
    @test length(wing_points) == 20  # 10 LE/TE pairs

    # Create a QUATERNION wing with these points
    vsm_set_path = joinpath(SymbolicAWEModels.get_data_path(), "vsm_settings.yaml")
    vsm_set = SymbolicAWEModels.VortexStepMethod.VSMSettings(vsm_set_path; data_prefix=false)
    vsm_wing = SymbolicAWEModels.Wing(set, vsm_set; prn=false)
    vsm_aero = SymbolicAWEModels.BodyAerodynamics([vsm_wing])
    vsm_solver = SymbolicAWEModels.Solver(vsm_aero;
                                          solver_type=SymbolicAWEModels.NONLIN,
                                          atol=2e-8, rtol=2e-8)

    # Create wing with QUATERNION type and empty group_idxs
    quat_wing = SymbolicAWEModels.VSMWing(
        SymbolicAWEModels.BaseWing(
            1, Int16[], I(3), zeros(3), ones(3);
            wing_type=SymbolicAWEModels.QUATERNION
        ),
        vsm_aero, vsm_wing, vsm_solver,
        Float64[], Float64[], zeros(0, 0),
        nothing, nothing, nothing, nothing, nothing
    )

    # Create SystemStructure (should auto-create groups)
    sys_quat = SymbolicAWEModels.SystemStructure(
        "v3_quat", set;
        points=sys_refine.points,
        segments=sys_refine.segments,
        pulleys=sys_refine.pulleys,
        tethers=sys_refine.tethers,
        winches=sys_refine.winches,
        wings=[quat_wing],
        transforms=sys_refine.transforms
    )

    # Verify groups were auto-created
    @test length(sys_quat.groups) == 10  # One group per LE/TE pair
    @test length(sys_quat.wings[1].group_idxs) == 10
    @test sys_quat.wings[1].wing_type == SymbolicAWEModels.QUATERNION

    println("✓ QUATERNION wing: $(length(sys_quat.groups)) groups auto-created")

    # Check that gamma values were calculated
    for (i, group) in enumerate(sys_quat.groups)
        println("  Group $i: gamma = $(group.gamma), points = $(group.point_idxs)")
        @test !iszero(group.gamma) || i == 5 || i == 6  # Center groups may have gamma ≈ 0
    end
end
