# Copyright (c) 2025 Bart van de Lint, Jelle Poland
# SPDX-License-Identifier: MPL-2.0

# build and display the html documentation locally
# run with: julia --project=docs scripts/build_docu.jl

using Pkg
Pkg.develop(path=dirname(@__DIR__))
Pkg.instantiate()
using LiveServer; servedocs(launch_browser=true)
