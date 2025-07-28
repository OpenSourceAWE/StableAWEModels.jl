# Copyright (c) 2024, 2025 Bart van de Lint, Uwe Fechner
# SPDX-License-Identifier: MIT

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
