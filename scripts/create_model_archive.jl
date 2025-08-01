# Copyright (c) 2024, 2025 Bart van de Lint, Uwe Fechner
# SPDX-License-Identifier: MIT

version = VERSION.minor
cp("Manifest-v1.$version.toml.default", "Manifest-v1.$version.toml"; force=true)
using Pkg
Pkg.instantiate()

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

SymbolicAWEModels.create_model_archive(get_data_path(), 
                                       joinpath(get_data_path(), "models_v1.$version.tar.gz"))

