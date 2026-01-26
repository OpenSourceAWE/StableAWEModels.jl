# SPDX-FileCopyrightText: 2025 Bart van de Lint
# SPDX-License-Identifier: MPL-2.0

ENV["MPLBACKEND"] = "Agg"
using KiteUtils
using Test
using SymbolicAWEModels

# Set up data path for 2plate_kite tests
pkg_file_path = Base.find_package("SymbolicAWEModels")
if isnothing(pkg_file_path)
    error("SymbolicAWEModels not found in the current project environment.")
else
    package_root_dir = dirname(dirname(pkg_file_path))
    src_data_path = joinpath(package_root_dir, "data", "2plate_kite")
    tmpdir = mktempdir()
    data_path = joinpath(tmpdir, "2plate_kite")
    cp(src_data_path, data_path; force=true)
    @show data_path
    set_data_path(data_path)
end

@testset verbose = true "Testing SymbolicAWEModels..." begin
    # Component tests
    println("--> Point dynamics")
    include("test_point.jl")
    println("--> Segment dynamics")
    include("test_segment.jl")
    println("--> Pulley constraints")
    include("test_pulley.jl")
    println("--> Tether and winch")
    include("test_tether_winch.jl")
    println("--> Transform coordinates")
    include("test_transform.jl")
    println("--> Wing aerodynamics")
    include("test_wing.jl")

    println("--> Quaternion conversions")
    include("test_quaternion_conversions.jl")
    println("--> Helpers")
    include("test_helpers.jl")
    println("--> Code quality")
    include("aqua.jl")
end
