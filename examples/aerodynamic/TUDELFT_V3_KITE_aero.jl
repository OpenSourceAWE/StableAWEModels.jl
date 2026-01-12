using LinearAlgebra
using VortexStepMethod
using DataFrames
using DelimitedFiles
using GLMakie
using CairoMakie
using LaTeXStrings
GLMakie.activate!()

project_dir = joinpath(@__DIR__, "..", "..")  # Go up two levels from examples to project root

# Geometry configuration
DEPOWER = 0.0
TIP_REDUCTION = 0.4
TE_FRAC = 0.95
GEOM_SUFFIX = "depower$(DEPOWER)_tip$(TIP_REDUCTION)_te$(TE_FRAC)"
AERO_YAML_PATH = joinpath(project_dir, "data", "v3",
    "aero_geometry_$(GEOM_SUFFIX).yaml")

literature_paths = [
    joinpath(project_dir, "data", "v3", "literature_results",
        "CFD_RANS_Rey_5e5_Poland2025_alpha_sweep_beta_0_NoStruts.csv"),
    joinpath(project_dir, "data", "v3", "literature_results",
        "CFD_RANS_Rey_10e5_Poland2025_alpha_sweep_beta_0.csv"),
    joinpath(project_dir, "data", "v3", "literature_results",
        "Python_VSM_Rey_5e5_Poland2025_alpha_sweep_beta_0.csv"),
    joinpath(project_dir, "data", "v3", "literature_results",
        "WindTunnel_Re_5e5_Poland2025_alpha_sweep_beta_0.csv"),
]
labels = [
    "Julia VSM",
    "CFD RANS Re=5e5",
    "CFD RANS Re=10e5 (With Struts)",
    "Python VSM Re=5e5",
    "Wind Tunnel Re=5e5 (With Struts)"
]


# Set up VSM objects for the TU Delft V3 kite
vsm_settings = VSMSettings(
    joinpath(project_dir, "data", "v3", "vsm_settings_reduced_for_coupling.yaml");
    data_prefix=false,
)
vsm_settings.wings[1].geometry_file = AERO_YAML_PATH
wing = VortexStepMethod.Wing(vsm_settings)
body_aero = VortexStepMethod.BodyAerodynamics([wing])
solver = VortexStepMethod.Solver(body_aero, vsm_settings)

# Apply flight conditions from settings
set_va!(body_aero, vsm_settings)
wind_speed = vsm_settings.condition.wind_speed
angle_of_attack_deg = vsm_settings.condition.alpha
sideslip_deg = vsm_settings.condition.beta

# Polar sweep including CMy, using solve! and calculate_results
function compute_polar_with_cm(
    solver,
    body_aero,
    angle_range;
    angle_type::String="angle_of_attack",
    angle_of_attack::Float64=0.0,
    side_slip::Float64=0.0,
    v_a::Float64=10.0
)
    n_angles = length(angle_range)
    cl = zeros(n_angles)
    cd = zeros(n_angles)
    cs = zeros(n_angles)
    cmy = fill(NaN, n_angles)
    reynolds_number = zeros(n_angles)

    gamma_prev = solver.sol.gamma_distribution
    for (i, angle_i) in enumerate(angle_range)
        if angle_type == "angle_of_attack"
            α = deg2rad(angle_i)
            β = deg2rad(side_slip)
        elseif angle_type == "side_slip"
            α = deg2rad(angle_of_attack)
            β = deg2rad(angle_i)
        else
            throw(ArgumentError("angle_type must be 'angle_of_attack' or 'side_slip'"))
        end

        set_va!(body_aero, [cos(α) * cos(β), sin(β), sin(α)] * v_a)
        solve!(solver, body_aero, gamma_prev; log=false)
        gamma_prev = solver.sol.gamma_distribution

        results_local = calculate_results(
            body_aero,
            solver.lr.gamma_new,
            zeros(MVec3),
            solver.density,
            solver.aerodynamic_model_type,
            solver.core_radius_fraction,
            solver.mu,
            solver.lr.alpha_dist,
            solver.lr.v_a_dist,
            solver.sol._chord_dist,
            solver.sol._x_airf_dist,
            solver.sol._y_airf_dist,
            solver.sol._z_airf_dist,
            solver.sol._va_dist,
            solver.br.va_norm_dist,
            solver.br.va_unit_dist,
            body_aero.panels,
            solver.is_only_f_and_gamma_output;
            correct_aoa=solver.correct_aoa
        )

        cl[i] = results_local["cl"]
        cd[i] = results_local["cd"]
        cs[i] = results_local["cs"]
        cmy[i] = get(results_local, "cmy", NaN)
        reynolds_number[i] = results_local["Rey"]
    end

    return (angle=angle_range, cl=cl, cd=cd, cs=cs, cmy=cmy, rey=reynolds_number)
end

function plot_polars_with_cmy(
    solver_list,
    body_aero_list,
    label_list;
    literature_path_list::Vector{String}=String[],
    angle_range=range(-10, 40, step=1),
    angle_type::String="angle_of_attack",
    angle_of_attack::Float64=0.0,
    side_slip::Float64=0.0,
    v_a::Float64=10.0,
    title::String="polar_with_cm",
    fig_size::Tuple{Int,Int}=(1200, 800),
    angle_xlim::Tuple{Real,Real}=(-10, 40)
)
    total_cases = length(body_aero_list) + length(literature_path_list)
    length(label_list) == total_cases || throw(ArgumentError("labels length ($(length(label_list))) must match number of cases ($total_cases)"))
    length(solver_list) == length(body_aero_list) || throw(ArgumentError("solver_list length must match body_aero_list length"))

    polar_data_list = Vector{Any}()
    labels_full = String[]

    # Computational cases
    for (solver_i, body, lbl) in zip(solver_list, body_aero_list, label_list[1:length(solver_list)])
        pd = compute_polar_with_cm(
            solver_i,
            body,
            angle_range;
            angle_type=angle_type,
            angle_of_attack=angle_of_attack,
            side_slip=side_slip,
            v_a=v_a
        )
        @info "polar sample (solver)" label=lbl first_cl=pd.cl[1] first_cd=pd.cd[1] first_cs=pd.cs[1] first_cmy=pd.cmy[1]
        push!(polar_data_list, pd)
        re_tag = round(Int, first(pd.rey) * 1e-5)
        push!(labels_full, "$(lbl) Re=$(re_tag)e5")
    end

    # Literature cases
    for (path, lbl) in zip(literature_path_list, label_list[length(solver_list)+1:end])
        data = readdlm(path, ',')
        header_raw = string.(data[1, :])
        header = lowercase.(strip.(header_raw))
        alpha_idx = findfirst(x -> occursin("alpha", x) || occursin("aoa", x), header)
        cl_idx    = findfirst(x -> occursin("cl", x), header)
        cd_idx    = findfirst(x -> occursin("cd", x), header)
        cs_idx    = findfirst(x -> occursin("cs", x), header)
        cmy_idx   = findfirst(x -> occursin("cmy", x), header)

        parse_col(col) = begin
            vals = Float64[]
            for v in col
                if v isa Real
                    push!(vals, Float64(v))
                else
                    s = strip(String(v))
                    y = tryparse(Float64, s)
                    push!(vals, isnothing(y) ? NaN : y)
                end
            end
            vals
        end

        alpha_col = parse_col(data[2:end, alpha_idx])
        cl_col    = parse_col(data[2:end, cl_idx])
        cd_col    = parse_col(data[2:end, cd_idx])
        cs_col    = cs_idx === nothing ? zeros(size(data, 1)-1) : parse_col(data[2:end, cs_idx])
        cmy_col   = cmy_idx === nothing ? nothing : parse_col(data[2:end, cmy_idx])

        push!(polar_data_list, (angle=alpha_col, cl=cl_col, cd=cd_col, cs=cs_col, cmy=cmy_col, rey=fill(NaN, length(alpha_col))))
        push!(labels_full, lbl)
    end

    fig = Figure(size=(1200, 400))

    ax_cl    = Axis(fig[1, 1], xlabel=L"\alpha \; [°]", ylabel=L"C_L \; [-]")
    ax_cd    = Axis(fig[1, 2], xlabel=L"\alpha \; [°]", ylabel=L"C_D \; [-]")
    ax_polar = Axis(fig[1, 3], xlabel=L"C_D \; [-]", ylabel=L"C_L \; [-]")

    xlims!(ax_cl, angle_xlim...)
    xlims!(ax_cd, angle_xlim...)

    colors = Makie.wong_colors()

    for (idx, (pd, lbl)) in enumerate(zip(polar_data_list, labels_full))
        color = colors[mod1(idx, length(colors))]

        lines!(ax_cl, pd.angle, pd.cl; color, label=lbl)
        scatter!(ax_cl, pd.angle, pd.cl; color, markersize=6)

        lines!(ax_cd, pd.angle, pd.cd; color, label=lbl)
        scatter!(ax_cd, pd.angle, pd.cd; color, markersize=6)

        lines!(ax_polar, pd.cd, pd.cl; color, label=lbl)
        scatter!(ax_polar, pd.cd, pd.cl; color, markersize=6)
    end

    Legend(fig[2, :], ax_cl; orientation=:horizontal)
    return fig
end


# Save PDF with CairoMakie
CairoMakie.activate!()
save_filename = joinpath(project_dir, "data", "v3", "polar_$(GEOM_SUFFIX).pdf")
fig = plot_polars_with_cmy(
    [solver],
    [body_aero],
    labels;
    literature_path_list=literature_paths,
    angle_range=range(-5, 40, step=1),
    angle_type="angle_of_attack",
    angle_of_attack=angle_of_attack_deg,
    side_slip=sideslip_deg,
    v_a=wind_speed,
    title="V3 Kite Polars ($(GEOM_SUFFIX))",
)
save(save_filename, fig)
@info "Polar plot saved to $save_filename"

# Display interactive plot with GLMakie
GLMakie.activate!()
fig = plot_polars_with_cmy(
    [solver],
    [body_aero],
    labels;
    literature_path_list=literature_paths,
    angle_range=range(-5, 40, step=1),
    angle_type="angle_of_attack",
    angle_of_attack=angle_of_attack_deg,
    side_slip=sideslip_deg,
    v_a=wind_speed,
    title="V3 Kite Polars ($(GEOM_SUFFIX))",
)
scr = display(fig)
wait(scr)



# # Plotting geometry
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


# # Plotting spanwise distributions
# body_y_coordinates = [panel.aero_center[2] for panel in body_aero.panels]

# PLOT && plot_distribution(
#     [body_y_coordinates],
#     [results],
#     ["VSM"];
#     title="CAD_spanwise_distributions_alpha_$(round(angle_of_attack_deg, digits=1))_delta_$(round(sideslip_deg, digits=1))_yaw_$(round(yaw_rate, digits=1))_v_a_$(round(wind_speed, digits=1))",
#     data_type=".pdf",
#     is_save=false,
#     is_show=true,
#     use_tex=USE_TEX
# )
