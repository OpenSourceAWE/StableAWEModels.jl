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

exclude = ["test_for_precompile.jl"]
test_files = filter(readdir(@__DIR__)) do f
    startswith(f, "test_") && endswith(f, ".jl") &&
        f ∉ exclude
end
sort!(test_files)

@testset verbose = true "Testing SymbolicAWEModels..." begin
    for f in test_files
        println("--> $f")
        include(f)
    end
end
