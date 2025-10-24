# coupled_2plate_linearized.jl
# Copyright (c) 2025
# SPDX-License-Identifier: MPL-2.0

"""
2-Plate kite: coupled aero–structural simulation using **linearized** VSM updates.

- Loads the 2-plate kite from YAML via SymbolicAWEModels
- Runs time-marching with `vsm_interval = 3` (relinearize every 3 steps)
- Optionally renders a few snapshots

Usage:
    julia --project -e 'include("coupled_2plate_linearized.jl"); Coupled2PlateLinearized.main()'
"""
module Coupled2PlateLinearized

using SymbolicAWEModels
using GLMakie
using Statistics

# ============= user settings =============
const MODEL_NAME    = "2plate_kite"
const GEOM_PATH     = joinpath("data", MODEL_NAME, "struc_geometry.yaml")
const SIM_TIME      = 10.0     # s
const N_PLOTS       = 3        # number of static snapshots
const REMAKE_CACHE  = false
const VSM_INTERVAL  = 3        # ← linearized VSM: relinearize every 3 steps
const LOG_EVERY     = 100      # steps
# =========================================

# If your repo provides a custom loader, keep this; otherwise the package export is used.
yaml_loader_path = joinpath(@__DIR__, "..", "yaml_loader.jl")
if isfile(yaml_loader_path)
    include(yaml_loader_path)  # provides load_sys_struct_from_yaml
end

"""
Return unique, sorted step indices including 0 and `n_steps`.
"""
function snapshot_steps(n_steps::Integer, n_plots::Integer)
    n     = max(n_plots, 2)
    steps = unique!(sort!(round.(Int, range(0, stop=n_steps, length=n))))
    steps[1] != 0         && pushfirst!(steps, 0)
    steps[end] != n_steps && push!(steps, n_steps)
    return steps
end

"""
Build and initialize the coupled model from YAML.
"""
function build_model(geom_path::AbstractString=GEOM_PATH; remake::Bool=REMAKE_CACHE)
    @assert isfile(geom_path) "Geometry file not found: $geom_path"
    set = SymbolicAWEModels.load_settings(MODEL_NAME)

    # Prefer project loader if included, else package function
    sys = if @isdefined load_sys_struct_from_yaml
        load_sys_struct_from_yaml(geom_path; system_name=MODEL_NAME, set=set)
    else
        SymbolicAWEModels.load_sys_struct_from_yaml(geom_path; system_name=MODEL_NAME, set=set)
    end

    sam = SymbolicAWEModel(set, sys)
    SymbolicAWEModels.init!(sam; remake=remake)
    return sam, set
end

"""
Run the simulation with linearized VSM (`vsm_interval = 3`) and plot snapshots.
"""
function main(; sim_time::Real=SIM_TIME, n_plots::Int=N_PLOTS)
    @info "Loading $MODEL_NAME from $GEOM_PATH"
    sam, set = build_model(GEOM_PATH; remake=REMAKE_CACHE)

    Δt      = 1.0 / max(1, getproperty(set, :sample_freq, 100))
    n_steps = max(1, round(Int, sim_time / Δt))
    keep    = snapshot_steps(n_steps, n_plots)

    @info "Running $(sim_time) s → $n_steps steps (Δt=$(round(Δt, digits=4)) s, vsm_interval=$(VSM_INTERVAL))…"

    # Initial snapshot
    try
        display(plot(sam.sys_struct; title="2-Plate Kite – t=0"))
    catch e
        @warn "Initial plot failed: $e"
    end

    # Time-marching loop (linearized aero between relinearizations)
    for step in 1:n_steps
        next_step!(sam; dt=Δt, vsm_interval=VSM_INTERVAL)

        if (step % LOG_EVERY == 0) || (step == n_steps)
            @info "Step $step/$n_steps (t=$(round(step*Δt, digits=2)) / $(sim_time) s)"
        end

        if step in keep
            # Render a snapshot of the current state
            try
                t = round(step*Δt, digits=2)
                display(plot(sam.sys_struct; title="2-Plate Kite – step $step (t=$t s)"))
            catch e
                @warn "Snapshot at step $step failed: $e"
            end
        end
    end

    # Quick text summary
    println("\n", "="^60)
    println("Final Simulation Results (t = $(sim_time) s)")
    println("="^60)
    if !isempty(sam.sys_struct.wings)
        w = sam.sys_struct.wings[1]
        pos = w.pos_w
        println("  Wing position: [$(round.(pos; digits=2)...)] m")
        println("  Elevation/Azimuth/Heading [deg]: ",
                "$(round(rad2deg(w.elevation), digits=2)), ",
                "$(round(rad2deg(w.azimuth),   digits=2)), ",
                "$(round(rad2deg(w.heading),   digits=2))")
    end
    dyn = filter(p -> p.type == SymbolicAWEModels.DYNAMIC, sam.sys_struct.points)
    if !isempty(dyn)
        avg = (mean(p.pos_w[1] for p in dyn),
               mean(p.pos_w[2] for p in dyn),
               mean(p.pos_w[3] for p in dyn))
        println("  Avg dynamic point pos: [$(round.(collect(avg); digits=2)...)] m")
    end
    println("="^60, "\n")

    @info "Done (linearized VSM with vsm_interval=$(VSM_INTERVAL))."
    return sam
end

end # module

# Run if executed directly
if abspath(PROGRAM_FILE) == @__FILE__
    Coupled2PlateLinearized.main()
end
