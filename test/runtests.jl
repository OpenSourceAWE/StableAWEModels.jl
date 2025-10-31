# SPDX-FileCopyrightText: 2025 Bart van de Lint
# SPDX-License-Identifier: MPL-2.0

ENV["MPLBACKEND"] = "Agg"
using KiteUtils
pkg_file_path = Base.find_package("SymbolicAWEModels")
if isnothing(pkg_file_path)
    error("SymbolicAWEModels not found in the current project environment.")
else
    package_root_dir = dirname(dirname(pkg_file_path))
    data_path = joinpath(package_root_dir, "data")
    @show data_path
    set_data_path(data_path)
end
using Test
using SymbolicAWEModels

@testset verbose = true "Testing SymbolicAWEModels..." begin
    println("--> Quaternion conversions")
    include("test_quaternion_conversions.jl")
    println("--> Initialization")
    include("test_initialization.jl")
    println("--> Simulation")
    include("test_simulation.jl")
    println("--> Tether properties")
    include("test_tether_properties.jl")
    println("--> Linearization")
    include("test_linearization.jl")
    println("--> Serialization")
    include("test_serialization.jl")
    println("--> State conversion")
    include("test_state_conversion.jl")
    println("--> Helpers")
    include("test_helpers.jl")
    println("--> Code quality")
    include("aqua.jl")
end
