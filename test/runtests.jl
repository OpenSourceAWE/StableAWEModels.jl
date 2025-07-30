# SPDX-FileCopyrightText: 2022, 2024, 2025 Uwe Fechner
# SPDX-License-Identifier: MIT

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
    println("--> 1")
    include("test_sam.jl")
    println("--> 2")
    include("test_helpers.jl")
    println("--> 3")
    include("aqua.jl")
end
