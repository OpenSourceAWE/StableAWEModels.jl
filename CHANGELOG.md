<!--
SPDX-FileCopyrightText: 2025 Uwe Fechner, Bart van de Lint
SPDX-License-Identifier: MPL-2.0
-->

# v0.2.0 01-08-2025
## Added
- Adds simple model and tether model
- Adds `copy_to_simple!` function, which copies the ram model state to the simple model state, uses the tether model to find the equivalent 1-segment spring properties of the tether
- Adds open-loop sim functions `sim!`, `sim_oscillate!`, `sim_turn!`
- Adds plotting function `plot(sys_struct::SystemStructure, sys_log::SysLog)`
- Adds documentation
- Adds new updated tests: test/test_sam.jl
## Fixed
- Fixes documentation
- Fixes the bug where the kite could not have negative position
## Changed
- Improved precompilation
- Breaking: `Segment` constructor has different arguments
## Removed
- Removed `.bin` files from git, will be added as release artifacts

# v0.1.3 18-07-2025
## Changed
- Add interface keyword arguments to `init!`

# v0.1.2 13-07-2025
## Changed
- Update VortexStepMethod.jl

# v0.1.1 13-07-2025
## Added
- Added a simple linearized model
## Changed
- Improved the reinitialization using scalar settings values
- Update KiteUtils and AtmosphericModels

# v0.1.0
- Moved the SymbolicAWEModel from KiteModels.jl to SymbolicAWEModels.jl
