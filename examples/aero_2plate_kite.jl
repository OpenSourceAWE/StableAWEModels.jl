# Copyright (c) 2025 Bart van de Lint, Jelle Poland
# SPDX-License-Identifier: MPL-2.0

using LinearAlgebra
using VortexStepMethod
using GLMakie

# Load custom Makie polar plotting functions
include(joinpath(@__DIR__, "makie_polar_plots.jl"))

project_dir = dirname(dirname(pathof(VortexStepMethod)))  # Go up one level from src to project root

# Load VSM vsm_settings from YAML configuration file
vsm_settings = VSMSettings("2plate_kite/vsm_settings.yaml")

# Reduce panel resolution for faster sweeps
for wing_cfg in vsm_settings.wings
    wing_cfg.n_panels = min(wing_cfg.n_panels, 10)
end

# Create wing, body_aero, and solver objects using vsm_settings
wing = VortexStepMethod.Wing(vsm_settings)
body_aero = VortexStepMethod.BodyAerodynamics([wing])
solver = VortexStepMethod.Solver(body_aero, vsm_settings)

# Set flight conditions from settings
set_va!(body_aero, vsm_settings)

# Extract values for plotting (optional - for reference)
wind_speed = vsm_settings.condition.wind_speed
angle_of_attack_deg = vsm_settings.condition.alpha
sideslip_deg = vsm_settings.condition.beta
yaw_rate = vsm_settings.condition.yaw_rate

# Run the solver
results = VortexStepMethod.solve(solver, body_aero; log=true)

# Using plotting modules, to create more comprehensive plots
PLOT = true
USE_TEX = false

# Plotting polars with custom Makie function
if PLOT
    fig_polar = plot_polars_makie(
        [solver],
        [body_aero],
        ["VSM Pyramid Model"],
        angle_range=range(0, 15, length=5),
        angle_type="angle_of_attack",
        angle_of_attack=angle_of_attack_deg,
        side_slip=sideslip_deg,
        v_a=wind_speed,
        title="$(wing.n_panels) panels $(wing.spanwise_distribution) pyramid model",
        size=(1400, 1000)
    )
    display(fig_polar)
end

# Plotting geometry (TODO: Convert to Makie)
# results = VortexStepMethod.solve(solver, body_aero; log=true)
# PLOT && plot_geometry(
#     body_aero,
#     "";
#     data_type=".svg",
#     save_path="",
#     is_save=false,
#     is_show=true,
#     view_elevation=15,
#     view_azimuth=-120,
#     use_tex=USE_TEX
# )

# Use Makie-based geometry plotting instead
if PLOT
    results = VortexStepMethod.solve(solver, body_aero; log=true)
    fig_geom = Figure(size=(1200, 800))
    ax_geom = Axis3(fig_geom[1, 1], aspect=:data,
                    xlabel="X [m]", ylabel="Y [m]", zlabel="Z [m]",
                    title="Wing Geometry",
                    azimuth=9/8*π)
    plot!(ax_geom, body_aero)
    display(fig_geom)
end

# Plotting spanwise distributions (TODO: Convert to Makie)
# body_y_coordinates = [panel.aero_center[2] for panel in body_aero.panels]
# 
# PLOT && plot_distribution(
#     [body_y_coordinates],
#     [results],
#     ["VSM"];
#     title="pyramid_spanwise_distributions_alpha_$(round(angle_of_attack_deg, digits=1))_delta_$(round(sideslip_deg, digits=1))_yaw_$(round(yaw_rate, digits=1))_v_a_$(round(wind_speed, digits=1))",
#     data_type=".pdf",
#     is_save=false,
#     is_show=true,
#     use_tex=USE_TEX
# )

nothing
