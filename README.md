<!--
SPDX-FileCopyrightText: 2025 Uwe Fechner, Bart van de Lint
SPDX-License-Identifier: MPL-2.0
-->

# SymbolicAWEModels

[![Stable](https://img.shields.io/badge/docs-stable-blue.svg)](https://OpenSourceAWE.github.io/SymbolicAWEModels.jl/stable)
[![Dev](https://img.shields.io/badge/docs-dev-blue.svg)](https://OpenSourceAWE.github.io/SymbolicAWEModels.jl/dev)
[![CI](https://github.com/OpenSourceAWE/SymbolicAWEModels.jl/actions/workflows/CI.yml/badge.svg)](https://github.com/OpenSourceAWE/SymbolicAWEModels.jl/actions/workflows/CI.yml)
[![Coverage](https://codecov.io/gh/OpenSourceAWE/SymbolicAWEModels.jl/branch/main/graph/badge.svg)](https://codecov.io/gh/OpenSourceAWE/SymbolicAWEModels.jl)
[![Aqua QA](https://raw.githubusercontent.com/JuliaTesting/Aqua.jl/master/badge.svg)](https://github.com/JuliaTesting/Aqua.jl)

## Overview

**SymbolicAWEModels.jl** provides modular symbolic models for simulating Airborne Wind Energy (AWE) systems, including:

- One or more wings (kites)
- Tethers (with or without pulleys)
- Winches
- Bridle systems

The kite is modeled as a deforming rigid body with quaternion dynamics for orientation. Aerodynamic forces and moments are computed using the [Vortex Step Method](https://github.com/Albatross-Kite-Transport/VortexStepMethod.jl). Tethers are modeled as point masses connected by spring-damper elements with realistic drag. Winches are modeled as motors/generators that can reel tethers in/out.

### Modular Subcomponents

- **AtmosphericModel** from [AtmosphericModels.jl](https://github.com/aenarete/AtmosphericModels.jl)
- **WinchModel** from [WinchModels.jl](https://github.com/aenarete/WinchModels.jl)
- **Aerodynamics** via [VortexStepMethod.jl](https://github.com/Albatross-Kite-Transport/VortexStepMethod.jl)

This package is part of the Julia Kite Power Tools ecosystem:

![Julia Kite Power Tools](docs/src/kite_power_tools.png)

---

## Installation

Install Julia using [juliaup](https://github.com/JuliaLang/juliaup):

```bash
curl -fsSL https://install.julialang.org | sh
juliaup add release
juliaup default release
```

**Quick Start:**

```bash
mkdir my_kite_project
cd my_kite_project
julia --project="."
```

Then add the package and copy examples:

```julia
using Pkg
pkg"add SymbolicAWEModels"

using SymbolicAWEModels
SymbolicAWEModels.init_module()  # Copies examples and installs dependencies
```

Run the interactive example menu:

```julia
include("examples/menu.jl")
```

> **Note:** The first run will be slow (several minutes) due to compilation. Run a second time for a significant speedup - subsequent runs will be much faster.

See the [Getting Started Guide](https://OpenSourceAWE.github.io/SymbolicAWEModels.jl/dev/getting_started/) for detailed instructions for registry users, cloned package users, and developers.

---

## Kite Models

SymbolicAWEModels provides the building blocks for assembling kite models from
YAML or Julia constructors. Ready-to-use kite models live in dedicated packages:

- **[RamAirKite.jl](https://github.com/OpenSourceAWE/RamAirKite.jl)** — Ram
  air kite with bridle system, 4-tether steering, and deformable wing groups
- **[V3Kite.jl](https://github.com/OpenSourceAWE/V3Kite.jl)** — TU Delft V3
  leading-edge-inflatable kite, YAML-based configuration

### 2-Plate Kite Example

A minimal coupled aero-structural model included in `data/2plate_kite/`:

```julia
using SymbolicAWEModels, VortexStepMethod

set_data_path("data/2plate_kite")
struc_yaml = joinpath(get_data_path(), "quat_struc_geometry.yaml")
set = Settings("system.yaml")
vsm_set = VortexStepMethod.VSMSettings(
    joinpath(get_data_path(), "vsm_settings.yaml"))

sys = load_sys_struct_from_yaml(struc_yaml;
    system_name="2plate_kite", set, vsm_set)
sam = SymbolicAWEModel(set, sys)
init!(sam)
```

For visualization with Makie, see the [Examples](https://OpenSourceAWE.github.io/SymbolicAWEModels.jl/dev/examples/) page.

![2-plate kite structure](docs/src/assets/2plate_kite_structure.png)

---

## See Also

- **Kite models:** [RamAirKite.jl](https://github.com/OpenSourceAWE/RamAirKite.jl), [V3Kite.jl](https://github.com/OpenSourceAWE/V3Kite.jl)
- [Research Fechner](https://research.tudelft.nl/en/publications/?search=Fechner+wind&pageSize=50&ordering=rating&descending=true) – scientific background for winches and tethers
- More kite models: [KiteModels.jl](https://github.com/ufechner7/KiteModels.jl)
- Meta-package: [KiteSimulators.jl](https://github.com/aenarete/KiteSimulators.jl)
- Utilities: [KiteUtils.jl](https://github.com/OpenSourceAWE/KiteUtils.jl)
- Component models: [WinchModels.jl](https://github.com/aenarete/WinchModels.jl), [KitePodModels.jl](https://github.com/aenarete/KitePodModels.jl), [AtmosphericModels.jl](https://github.com/aenarete/AtmosphericModels.jl)
- Controllers and viewers: [KiteControllers.jl](https://github.com/aenarete/KiteControllers.jl), [KiteViewers.jl](https://github.com/aenarete/KiteViewers.jl)
- Aerodynamics: [VortexStepMethod.jl](https://github.com/Albatross-Kite-Transport/VortexStepMethod.jl)

---

## Questions?

- Submit an [issue](https://github.com/OpenSourceAWE/SymbolicAWEModels.jl/issues/new)
- Start a [discussion](https://github.com/OpenSourceAWE/SymbolicAWEModels.jl/discussions/new/choose)
- Ask on [Julia Discourse](https://discourse.julialang.org/)
- Email Bart van de Lint: bart@vandelint.net

**Authors:**  
Bart van de Lint (bart@vandelint.net)  
Uwe Fechner (uwe.fechner.msc@gmail.com)

---

## License

This project is licensed under the [MPL-2.0 License](LICENSE).

---

## Citing SymbolicAWEModels

If you use SymbolicAWEModels in your research, please cite this repository:

```bibtex
@misc{SymbolicAWEModels,
  author = {Bart van de Lint, Uwe Fechner, Jelle Poland},
  title = {{SymbolicAWEModels}: Symbolic airborne wind energy system models},
  year = {2025},
  publisher = {GitHub},
  journal = {GitHub repository},
  howpublished = {\url{[https://github.com/OpenSourceAWE/SymbolicAWEModels.jl]}},
}
```

## Copyright Notice

Technische Universiteit Delft hereby disclaims all copyright interest in the package “SymbolicAWEModels.jl” (symbolic models for airborne wind energy systems) written by the Author(s).

Prof.dr. H.G.C. (Henri) Werij, Dean of Aerospace Engineering, Technische Universiteit Delft.

See copyright notices in the source files and the list of authors in [AUTHORS.md](AUTHORS.md).

**Documentation** [Stable Version](https://OpenSourceAWE.github.io/SymbolicAWEModels.jl/stable) --- [Development Version](https://OpenSourceAWE.github.io/SymbolicAWEModels.jl/dev)
