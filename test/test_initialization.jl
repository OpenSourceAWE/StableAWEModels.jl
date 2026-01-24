# SPDX-FileCopyrightText: 2025 Bart van de Lint
# SPDX-License-Identifier: MPL-2.0

using Test
using SymbolicAWEModels

tmpdir=mktempdir()
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

@testset verbose=true "Initialization" begin
    function init(elevation, azimuth, heading)
        set.elevation = elevation
        set.azimuth = azimuth
        set.heading = heading
        init!(sam)
        ss = SysState(sam)
        @test sam.sys_struct.wings[1].elevation ≈ deg2rad(set.elevation) atol=1e-2
        @test sam.sys_struct.wings[1].azimuth ≈ deg2rad(set.azimuth) atol=1e-2
        @test sam.sys_struct.wings[1].heading ≈ deg2rad(set.heading) atol=1e-2
        @test ss.elevation ≈ deg2rad(set.elevation) atol=1e-2
        @test ss.azimuth ≈ deg2rad(set.azimuth) atol=1e-2
        @test ss.heading ≈ deg2rad(set.heading) atol=1e-2
    end

    @testset "Model types" begin
        @test sam isa SymbolicAWEModel
        @test sam isa AbstractKiteModel
        @test simple_sam isa SymbolicAWEModel
        @test simple_sam isa AbstractKiteModel
        @test tether_sam isa SymbolicAWEModel
        @test tether_sam isa AbstractKiteModel
    end

    @testset "Initialization timing" begin
        init_time = @elapsed init!(sam; prn=true)
        @show init_time
        @test init_time < 700
        init!(sam; prn=true)
        init_time = @elapsed init!(sam; prn=true)
        @test init_time < 0.3
    end

    @testset "Initialization with different angles" begin
        init(ones(3)...)
        init(zeros(3)...)
    end
end

set_data_path(old_data_path)
