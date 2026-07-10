# SPDX-FileCopyrightText: 2024 Uwe Fechner
# SPDX-License-Identifier: MIT

using Pkg
if abspath(PROGRAM_FILE) == abspath(@__FILE__)
    Pkg.activate(@__DIR__)
end

using Aqua, StableAWEModels, Test

@testset "Aqua.jl" begin
    Aqua.test_all(
        StableAWEModels;
        stale_deps=(ignore=[:CodecXz, :REPL],), # CodecXz is used during precompilation only
        piracies = false,                        # the norm function is doing piracy for performance reasons
        persistent_tasks = false,
    )
end
nothing
