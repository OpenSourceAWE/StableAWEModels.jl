# Copyright (c) 2024, 2025 Bart van de Lint, Uwe Fechner
# SPDX-License-Identifier: MIT

@info "Updating"
using Pkg
Pkg.update()
m1 = "Manifest-v1.$(VERSION.minor).toml"
if !isfile(m1)
    mv("Manifest.toml", m1)
end
if isfile(m1 * ".default")
    mv(m1 * ".default", m1 * ".default.bak"; force=true)
end
cp(m1, m1 * ".default")

using Timers
@info "Loading packages"
tic()
using SymbolicAWEModels 
toc()

prn = true
@info "Creating default models"
models = SymbolicAWEModels.create_default_models(; prn)
@info "Created all models"
toc()

for model in models
    name = SymbolicAWEModels.get_model_name(model.set)
    input_path = joinpath(get_data_path(), name)
    output_path = input_path * ".xz"
    SymbolicAWEModels.compress_binary(input_path, output_path)
    println("Compressed $input_path => $output_path")
end
