# SPDX-FileCopyrightText: 2024 Uwe Fechner
# SPDX-License-Identifier: MIT

using Aqua
@testset "Aqua.jl" begin
    Aqua.test_all(
        SymbolicAWEModels;
        stale_deps=(ignore=[:CodecXz, :REPL],), # CodecXz is used during precompilation only
        piracies = false,                        # the norm function is doing piracy for performance reasons
        persistent_tasks = false,
    )
end
