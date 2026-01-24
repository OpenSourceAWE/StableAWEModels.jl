# SPDX-FileCopyrightText: 2025 Bart van de Lint
# SPDX-License-Identifier: MPL-2.0

using Test
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

@testset verbose=true "Initialization" begin
    function init(elevation, azimuth, heading)
        transform = sam.sys_struct.transforms[1]
        transform.elevation = deg2rad(elevation)
        transform.azimuth = deg2rad(azimuth)
        transform.heading = deg2rad(heading)
        init!(sam)
        ss = SysState(sam)
        @test sam.sys_struct.wings[1].elevation ≈ transform.elevation atol=1e-2
        @test sam.sys_struct.wings[1].azimuth ≈ transform.azimuth atol=1e-2
        @test sam.sys_struct.wings[1].heading ≈ transform.heading atol=1e-2
        @test ss.elevation ≈ transform.elevation atol=1e-2
        @test ss.azimuth ≈ transform.azimuth atol=1e-2
        @test ss.heading ≈ transform.heading atol=1e-2
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
