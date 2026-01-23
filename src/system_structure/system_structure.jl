# Copyright (c) 2025 Bart van de Lint
# SPDX-License-Identifier: MPL-2.0

"""
System Structure Module

This module defines the physical structure of the kite system, including:
- Basic types: Point, Group, Segment, Pulley, Tether, Winch
- Wing types: BaseWing, VSMWing with VSM aerodynamics coupling
- Transform type for spatial positioning
- SystemStructure container and initialization

Files are organized as:
- types.jl: Enums and basic struct definitions (Point, Segment, etc.)
- wing.jl: Wing types and VSM-related code
- transforms.jl: Transform type and heading/rotation functions
- system_structure_core.jl: SystemStructure type and constructor
- utilities.jl: Helper functions and state management
"""

# Include submodules in dependency order
include("types.jl")         # Enums, Point, Group, Segment, Pulley, Tether, Winch
include("wing.jl")          # AbstractWing, BaseWing, VSMWing
include("transforms.jl")    # Transform, heading calculations, reinit!/reposition!
include("system_structure_core.jl")  # SystemStructure type and constructor
include("utilities.jl")     # Helpers, validation, state management
