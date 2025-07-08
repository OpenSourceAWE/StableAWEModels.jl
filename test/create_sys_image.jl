# SPDX-FileCopyrightText: 2025 Bart van de Lint
# SPDX-License-Identifier: MPL-2.0

using Pkg
if ! ("PackageCompiler" ∈ keys(Pkg.project().dependencies))
    using TestEnv; TestEnv.activate()
    Pkg.update()
end
using PackageCompiler, BenchmarkTools, Documenter

@info "Creating sysimage ..."
PackageCompiler.create_sysimage(
    [:SymbolicAWEModels, :ControlPlots];
    sysimage_path="kps-image_tmp.so",
    precompile_execution_file=joinpath("test", "test_for_precompile.jl")
)
