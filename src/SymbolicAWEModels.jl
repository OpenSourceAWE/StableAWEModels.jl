# Copyright (c) 2025 Bart van de Lint and Uwe Fechner
# SPDX-License-Identifier: MPL-2.0

module SymbolicAWEModels

#======================================================================#
#                         DEPENDENCIES
#======================================================================#

# --- Julia Standard Library & General Utilities ---
using Pkg
using PrecompileTools: @setup_workload, @compile_workload
using DocStringExtensions
using LinearAlgebra
using Parameters
using Printf
using Serialization
using SHA
using CodecZlib
using Tar
using Statistics
using Suppressor
using Timers

# --- Numerical, Modeling & Scientific Computing ---
using ModelingToolkit
using ControlSystemsBase
using RecipesBase
using StaticArrays
using SymbolicIndexingInterface
using SymbolicIndexingInterface: AbstractIndexer

# --- Solvers (Nonlinear, Differential Equations) ---
using NonlinearSolve
using OrdinaryDiffEqBDF
using OrdinaryDiffEqCore
using OrdinaryDiffEqNonlinearSolve
using SteadyStateDiffEq

# --- Open Source AWE Packages ---
using AtmosphericModels
using KiteUtils
using VortexStepMethod
using WinchModels

#======================================================================#
#                  IMPORTS (for extending functions)
#======================================================================#

import KiteUtils: init!, next_step!, update_sys_state!
import ModelingToolkit: t_nounits as t, D_nounits as D
import ModelingToolkit.SciMLBase: successful_retcode

#======================================================================#
#                          EXPORTS
#                 (The Public API of this Module)
#======================================================================#

# --- KiteUtils ---
export init!, next_step!, update_sys_state!, get_data_path, set_data_path, se
export SysState, Settings, AbstractKiteModel

# --- Types ---
# Core Model
export SymbolicAWEModel
# System Structure Components
export SystemStructure, Point, Group, Segment, Pulley, Tether, Winch, Wing, Transform
# Enums
export DynamicsType, DYNAMIC, QUASI_STATIC, WING, STATIC
export SegmentType, POWER_LINE, STEERING_LINE, BRIDLE

# --- High-Level Simulation Functions (Workers) ---
export sim!, sim_oscillate!, sim_turn!

# --- Low-Level Simulation Functions ---
export find_steady_state!
export linearize!, simple_linearize!
export copy_to_simple!

# --- System Structure Creators ---
export create_ram_sys_struct
export create_tether_sys_struct
export create_simple_ram_sys_struct

# --- Getter Functions ---
export winch_force
export unstretched_length
export tether_length

# --- Helper Functions ---
export copy_examples
export copy_bin

set_zero_subnormals(true)       # required to avoid drastic slow down on Intel CPUs when numbers become very small

# Type definitions
"""
    const SimFloat = Float64

This type is used for all real variables, used in the Simulation. Possible alternatives: Float32, Double64, Dual
Other types than Float64 or Float32 do require support of Julia types by the solver. 
"""
const SimFloat = Float64

"""
   const KVec3    = MVector{3, SimFloat}

Basic 3-dimensional vector, stack allocated, mutable.
"""
const KVec3    = MVector{3, SimFloat}
const KVec4    = MVector{4, SimFloat}

"""
   const SVec3    = SVector{3, SimFloat}

Basic 3-dimensional vector, stack allocated, immutable.
"""
const SVec3    = SVector{3, SimFloat}  

# Defined in ext/SymbolicAWEModelsControlPlotsExt.jl
function plot end

function __init__()
    if isdir(joinpath(pwd(), "data")) && isfile(joinpath(pwd(), "data", "system.yaml"))
        set_data_path(joinpath(pwd(), "data"))
    end
end

include("system_structure.jl")
include("symbolic_awe_model.jl")
include("predefined_structures.jl")
include("tether_properties.jl")
include("linearize.jl")
include("generate_system.jl")
include("plot_recipe.jl")
include("simulate.jl")

function upwind_dir(v_wind_gnd)
    if v_wind_gnd[1] == 0.0 && v_wind_gnd[2] == 0.0
        return NaN
    end
    wind_dir = atan(v_wind_gnd[2], v_wind_gnd[1])
    -(wind_dir + π/2)
end

# rotate a 3d vector around the x axis in the yz plane - following the right hand rule
function rotate_around_x(vec, angle::T) where T
    result = zeros(T, 3)
    result[1] = vec[1]
    result[2] = cos(angle) * vec[2] - sin(angle) * vec[3]
    result[3] = sin(angle) * vec[2] + cos(angle) * vec[3]
    result
end

# rotate a 3d vector around the y axis in the xz plane - following the right hand rule
function rotate_around_y(vec, angle::T) where T
    result = zeros(T, 3)
    result[1] = cos(angle) * vec[1] + sin(angle) * vec[3]
    result[2] = vec[2]
    result[3] = -sin(angle) * vec[1] + cos(angle) * vec[3]
    result
end

# rotate a 3d vector around the z axis in the yx plane - following the right hand rule
function rotate_around_z(vec, angle::T) where T
    result = zeros(T, 3)
    result[1] = cos(angle) * vec[1] - sin(angle) * vec[2]
    result[2] = sin(angle) * vec[1] + cos(angle) * vec[2]
    result[3] = vec[3]
    result
end

"""
    copy_examples()

Copy all example scripts to the folder "examples"
(it will be created if it doesn't exist).
"""
function copy_examples()
    src_data_path = joinpath(dirname(pathof(@__MODULE__)), "..", "examples")
    dst_data_path = abspath(joinpath(pwd(), "examples"))
    copy_dir_no_overwrite(src_data_path, dst_data_path)
end

"""
    copy_bin()

Copy all bin scripts to the folder "bin"
(it will be created if it doesn't exist).
"""
function copy_bin()
    src_data_path = joinpath(dirname(pathof(@__MODULE__)), "..", "bin")
    dst_data_path = abspath(joinpath(pwd(), "bin"))
    copy_dir_no_overwrite(src_data_path, dst_data_path)
end

"""
    copy_data()

Copy all data scripts to the folder "data"
(it will be created if it doesn't exist).
"""
function copy_data()
    src_data_path = joinpath(dirname(pathof(@__MODULE__)), "..", "data")
    dst_data_path = abspath(joinpath(pwd(), "data"))
    copy_dir_no_overwrite(src_data_path, dst_data_path)
end

"""
    copy_dir_no_overwrite(src_dir, dst_dir)

Copies all files from `src_dir` to `dst_dir`.
Overwrites existing files if force=true.
Creates `dst_dir` if it does not exist.
"""
function copy_dir(src_dir::AbstractString, dst_dir::AbstractString; force=false)
    if !isdir(dst_dir)
        mkdir(dst_dir)
    end
    for file in readdir(src_dir)
        src_file = joinpath(src_dir, file)
        dst_file = joinpath(dst_dir, file)
        if force || (isfile(src_file) && !isfile(dst_file))
            cp(src_file, dst_file; force=true)
            chmod(dst_file, 0o774)
        elseif isdir(src_file)
            copy_dir(src_file, dst_file; force)
        end
    end
end

"""
    init_module()

Initializes the module in the current working directory.

- Copies all files from the module's `data` directory to `pwd()/data`, but does NOT overwrite files that already exist in `pwd()/data`.
- Copies all example scripts to `pwd()/examples`, creating the folder if it does not exist, and does not overwrite existing files.
- Copies the `bin` directory from the module to `pwd()/bin`, creating the folder if it does not exist, and does not overwrite existing files.
- Installs all required packages if they are not already installed.
"""
function init_module()
    copy_data()
    copy_examples()
    copy_bin()

    # Install required packages if not already present
    pkgs = ["KiteUtils", "ControlPlots", 
            "LaTeXStrings", "Timers"]
    for pkg in pkgs
        if !(pkg in keys(Pkg.project().dependencies))
            Pkg.add(pkg)
        end
    end

    println("Initialization complete! Data, examples, and binaries are prepared in the current directory.")
end

"""
    copy_bin()

Copy the scripts create_sys_image and run_julia to the folder "bin"
(it will be created if it doesn't exist).
"""
function copy_bin()
    PATH = "bin"
    if ! isdir(PATH) 
        mkdir(PATH)
    end
    src_path = joinpath(dirname(pathof(@__MODULE__)), "..", PATH)
    cp(joinpath(src_path, "create_sys_image2"), joinpath(PATH, "create_sys_image"), force=true)
    cp(joinpath(src_path, "run_julia"), joinpath(PATH, "run_julia"), force=true)
    chmod(joinpath(PATH, "create_sys_image"), 0o774)
    chmod(joinpath(PATH, "run_julia"), 0o774)
    PATH = "test"
    if ! isdir(PATH) 
        mkdir(PATH)
    end
    src_path = joinpath(dirname(pathof(@__MODULE__)), "..", PATH)
    cp(joinpath(src_path, "create_sys_image2.jl"), joinpath(PATH, "create_sys_image.jl"), force=true)
    cp(joinpath(src_path, "test_for_precompile.jl"), joinpath(PATH, "test_for_precompile.jl"), force=true)
    cp(joinpath(src_path, "update_packages.jl"), joinpath(PATH, "update_packages.jl"), force=true)
    chmod(joinpath(PATH, "create_sys_image.jl"), 0o664)
    chmod(joinpath(PATH, "test_for_precompile.jl"), 0o664)
    chmod(joinpath(PATH, "update_packages.jl"), 0o664)
end

include("precompile.jl")

end
