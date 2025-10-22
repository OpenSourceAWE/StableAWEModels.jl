# SPDX-FileCopyrightText: 2025 Bart van de Lint
# SPDX-License-Identifier: MPL-2.0

using Test
using SymbolicAWEModels
using Serialization
using ModelingToolkit: @variables

tmpdir=mktempdir()
mkpath(joinpath(tmpdir, "data"))
old_data_path = get_data_path()
set_data_path(joinpath(tmpdir, "data"))
cp(old_data_path, get_data_path(); force=true)

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

@testset verbose=true "Serialization" begin
    @testset "Serialization and Deserialization" begin
        points = [
            Point(1, zeros(3), DYNAMIC; wing_idx=0, transform_idx=1)
            Point(2, ones(3), DYNAMIC; wing_idx=0, transform_idx=1)
        ]
        segments = [
            Segment(1, set, (1,2), BRIDLE)
        ]
        transforms = [Transform(1, zeros(3)...;
            base_pos=zeros(3), base_point_idx=1, rot_point_idx=1)]
        sys_struct = SystemStructure("one_point", set; points, segments, transforms)
        sam = SymbolicAWEModel(set, sys_struct)
        model_path = joinpath(get_data_path(), SymbolicAWEModels.get_model_name(set))

        function test_init_with_reset(create_prob, create_lin_prob, create_control_func)
            println("Create prob: $create_prob \t"*
                    "lin_prob: $create_lin_prob \t"*
                    "control_func: $create_control_func")
            rm(model_path; force=true)
            sam = SymbolicAWEModel(set, sys_struct)
            init!(sam; create_prob, create_lin_prob, create_control_func, prn=false)
            @test isnothing(sam.prob) == !create_prob
            @test isnothing(sam.lin_prob) == !create_lin_prob
            @test isnothing(sam.control_funcs) == !create_control_func
        end
        test_init_with_reset(false, false, false)
        test_init_with_reset(true, false, false)
        test_init_with_reset(false, true, false)
        test_init_with_reset(false, false, true)

        init!(sam; create_prob=true, create_lin_prob=true, create_control_func=true, prn=false)
        # check if not removed
        init!(sam; create_prob=false, create_lin_prob=false, create_control_func=false, prn=false)
        @test !isnothing(sam.prob)
        @test !isnothing(sam.lin_prob)
        @test !isnothing(sam.control_funcs)

        # same name, check if hash works
        push!(points, Point(3, [0, 0, 2], DYNAMIC; wing_idx=0, transform_idx=1))
        sys_struct2 = SystemStructure("one_point", set; points, segments, transforms)
        sam.sys_struct = sys_struct2
        model_path2 = joinpath(get_data_path(), SymbolicAWEModels.get_model_name(set))
        @test model_path == model_path2
        # should create nothing because hashes are broken
        init!(sam; create_prob=false, create_lin_prob=false, create_control_func=false, prn=false)
        @test isnothing(sam.prob)
        @test isnothing(sam.lin_prob)
        @test isnothing(sam.control_funcs)
        init!(sam; create_prob=true, create_lin_prob=true, create_control_func=true, prn=false)
        @test !isnothing(sam.prob)
        @test !isnothing(sam.lin_prob)
        @test !isnothing(sam.control_funcs)

        # changing from 0 to 1 output
        old_ny = length(sam.outputs)
        outputs = [sam.prob.sys.pos[1,1]]
        init!(sam; outputs)
        @test old_ny == 0
        @test length(sam.outputs) == 1
        lin_model = linearize!(sam)
        @test size(lin_model.C)[1] == 1

        # changing from 1 to 2 outputs
        old_ny = length(sam.outputs)
        outputs = [sam.prob.sys.pos[1,1], sam.prob.sys.pos[2,1]]
        init!(sam; outputs)
        @test old_ny == 1
        @test length(sam.outputs) == 2
        lin_model = linearize!(sam)
        @test size(lin_model.C)[1] == 2
    end
end

set_data_path(old_data_path)
