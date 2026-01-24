# SPDX-FileCopyrightText: 2025 Bart van de Lint
# SPDX-License-Identifier: MPL-2.0

using Test, Printf
using SymbolicAWEModels

# Set up tmpdir if not already done by runtests.jl
if !startswith(get_data_path(), tempdir())
    src_data_path = joinpath(dirname(dirname(pathof(SymbolicAWEModels))), "data", "ram_air_kite")
    tmpdir = mktempdir()
    set_data_path(joinpath(tmpdir, "ram_air_kite"))
    cp(src_data_path, get_data_path(); force=true)
end

set = Settings("system.yaml")
sam = SymbolicAWEModel(set, "ram")

tether_set = Settings("system.yaml")
tether_sam = SymbolicAWEModel(tether_set, "tether")
init!(tether_sam)

simple_set = Settings("system.yaml")
simple_sam = SymbolicAWEModel(simple_set, "simple_ram")
init!(simple_sam)

original_set = Settings("system.yaml")

function reset!(set::Settings)
    for field in fieldnames(Settings)
        setfield!(set, field, getfield(original_set, field))
    end
    return set
end

@testset verbose=true "Tether properties" begin
    @testset "Tether spring properties" begin
        reset!(set)
        set.sample_freq = 600
        set.abs_tol = 1e-6
        set.rel_tol = 1e-6
        set.segments = 1
        one_seg_sam = SymbolicAWEModel(set, "ram")
        init!(one_seg_sam)
        one_seg_tether_sam = SymbolicAWEModel(set, "tether")
        init!(one_seg_tether_sam)

        # run twice to make sure state is reset properly
        axial_stiffness, axial_damping =
            SymbolicAWEModels.calc_spring_props(one_seg_sam, one_seg_tether_sam)
        next_step!(one_seg_sam; dt=1.0)
        axial_stiffness, axial_damping =
            SymbolicAWEModels.calc_spring_props(one_seg_sam, one_seg_tether_sam)
        segments = one_seg_sam.sys_struct.segments
        tethers = one_seg_sam.sys_struct.tethers
        segments = [segments[tether.segment_idxs[1]] for tether in tethers]
        real_axial_stiffness = [segment.axial_stiffness for segment in segments]
        real_axial_damping = [segment.axial_damping for segment in segments]
        @test isapprox(real_axial_stiffness, axial_stiffness; rtol=0.02)
        @test isapprox(real_axial_damping, axial_damping; rtol=0.2)

        println("\n--- Tether Spring Properties ---")
        # Print table headers
        @printf "%-8s | %-15s %-15s %-10s | %-15s %-15s %-10s\n" "Tether" "Calc. Stiffness" "Real Stiffness" "Error (%)" "Calc. Damping" "Real Damping" "Error (%)"
        # Print separator line
        println(repeat("-", 100))
        for i in 1:4
            # Calculate relative errors in percent
            stiffness_err = 100 * abs(axial_stiffness[i] - real_axial_stiffness[i]) / real_axial_stiffness[i]
            damping_err   = 100 * abs(axial_damping[i] - real_axial_damping[i]) / real_axial_damping[i]
            # Print data rows
            @printf "%-8d | %-15.2f %-15.2f %-10.2f | %-15.2f %-15.2f %-10.2f\n" i axial_stiffness[i] real_axial_stiffness[i] stiffness_err axial_damping[i] real_axial_damping[i] damping_err
        end
        println()
    end

    @testset "Test calc winch force" begin
        reset!(set)
        init!(sam)
        tether_vel = [winch.tether_vel for winch in sam.sys_struct.winches]
        tether_acc = [winch.tether_acc for winch in sam.sys_struct.winches]
        set_values = [winch.set_value for winch in sam.sys_struct.winches]
        winch_force = SymbolicAWEModels.calc_winch_force(sam.sys_struct, tether_vel, tether_acc, set_values)
        ss = SysState(sam)
        @test all(isapprox(ss.winch_force[1:3], winch_force))
    end

    @testset "Just a tether, without winch or kite" begin
        set.segments = 20
        dynamics_type = DYNAMIC

        points = Point[]
        segments = Segment[]

        points = push!(points, Point(1, zeros(3), STATIC; wing_idx=0))

        segment_idxs = Int[]
        for i in 1:set.segments
            point_idx = i+1
            pos = [0.0, 0.0, i * set.l_tether / set.segments]
            push!(points, Point(point_idx, pos, dynamics_type; wing_idx=0))
            segment_idx = i
            push!(segments, Segment(segment_idx, set, (point_idx-1, point_idx), BRIDLE))
            push!(segment_idxs, segment_idx)
        end

        transforms = [Transform(1, deg2rad(-80), 0.0, 0.0;
            base_pos=[0.0, 0.0, 50.0], base_point_idx=points[1].idx, rot_point_idx=points[end].idx)]
        sys_struct = SymbolicAWEModels.SystemStructure("tether", set; points, segments, transforms)

        sam = SymbolicAWEModel(set, sys_struct)
        init!(sam; remake=false)
        sys = sam.prob.sys
        @test isapprox(sam.integrator[sys.pos[:, end]], [8.682408883346524, 0.0, 0.7596123493895988], atol=1e-2)
        for i in 1:100
            next_step!(sam)
        end
        @test sam.integrator[sys.pos[1, end]] > 0.8set.l_tether
        @test isapprox(sam.integrator[sys.pos[2, end]], 0.0, atol=1.0)
    end
end
