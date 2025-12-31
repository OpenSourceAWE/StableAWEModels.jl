using LinearAlgebra
using VortexStepMethod
using DataFrames
using DelimitedFiles
using GLMakie

# Compute CL/CD/CS/CM polars for a solver/body_aero sweep
function compute_polar_with_cm(solver, body_aero, angle_range;
        angle_type::String="angle_of_attack",
        angle_of_attack::Float64=0.0,
        side_slip::Float64=0.0,
        v_a::Float64=10.0)

    n_angles = length(angle_range)
    cl = zeros(n_angles)
    cd = zeros(n_angles)
    cs = zeros(n_angles)
    cmy = fill(NaN, n_angles)  # moment coefficient about body y (pitch)
    reynolds_number = zeros(n_angles)

    gamma_prev = nothing
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
        results = solve(solver, body_aero, gamma_prev)

        cl[i] = results["cl"]
        cd[i] = results["cd"]
        cs[i] = results["cs"]
        cmy[i] = get(results, "cmy", NaN)
        reynolds_number[i] = results["Rey"]
        gamma_prev = results["gamma_distribution"]
    end

    return (angle=angle_range, cl=cl, cd=cd, cs=cs, cmy=cmy, rey=reynolds_number)
end

# Plot polars with optional CM (Makie; plots CM only when available)
function plot_polars_with_cmy(solver_list, body_aero_list, label_list;
        literature_path_list::Vector{String}=String[],
        angle_range=range(-10, 40, step=1),
        angle_type::String="angle_of_attack",
        angle_of_attack::Float64=0.0,
        side_slip::Float64=0.0,
        v_a::Float64=10.0,
        title::String="polar_with_cm",
        fig_size::Tuple{Int,Int}=(1200, 800),
        angle_xlim::Tuple{Real,Real}=(-10, 40))

    total_cases = length(body_aero_list) + length(literature_path_list)
    length(label_list) == total_cases || throw(ArgumentError("labels length ($(length(label_list))) must match number of cases ($total_cases)"))
    length(solver_list) == length(body_aero_list) || throw(ArgumentError("solver_list length must match body_aero_list length"))

    polar_data_list = Vector{Any}()
    labels_full = String[]

    # Computational cases
    for (solver, body, lbl) in zip(solver_list, body_aero_list, label_list[1:length(solver_list)])
        pd = compute_polar_with_cm(solver, body, angle_range;
            angle_type=angle_type,
            angle_of_attack=angle_of_attack,
            side_slip=side_slip,
            v_a=v_a)
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

    fig = Figure(size=fig_size)
    Label(fig[0, :], title; fontsize=20, font=:bold)

    ax_cl    = Axis(fig[1, 1], title="CL vs $angle_type [deg]",  xlabel="$angle_type [deg]", ylabel="CL [-]")
    ax_cd    = Axis(fig[1, 2], title="CD vs $angle_type [deg]",  xlabel="$angle_type [deg]", ylabel="CD [-]")
    ax_cmy   = Axis(fig[1, 3], title="CMy vs $angle_type [deg]", xlabel="$angle_type [deg]", ylabel="CMy [-]")
    ax_cs    = Axis(fig[2, 1], title="CS vs $angle_type [deg]",  xlabel="$angle_type [deg]", ylabel="CS [-]")
    ax_polar = Axis(fig[2, 2], title="CL vs CD", xlabel="CD [-]", ylabel="CL [-]")

    # Fix angle x-limits on all angle-based subplots
    xlims!(ax_cl, angle_xlim...)
    xlims!(ax_cd, angle_xlim...)
    xlims!(ax_cmy, angle_xlim...)
    xlims!(ax_cs, angle_xlim...)

    colors = Makie.wong_colors()

    for (idx, (pd, lbl)) in enumerate(zip(polar_data_list, labels_full))
        color = colors[mod1(idx, length(colors))]

        lines!(ax_cl, pd.angle, pd.cl; color, label=lbl)
        scatter!(ax_cl, pd.angle, pd.cl; color, markersize=6)

        lines!(ax_cd, pd.angle, pd.cd; color, label=lbl)
        scatter!(ax_cd, pd.angle, pd.cd; color, markersize=6)

        lines!(ax_cs, pd.angle, pd.cs; color, label=lbl)
        scatter!(ax_cs, pd.angle, pd.cs; color, markersize=6)

        lines!(ax_polar, pd.cd, pd.cl; color, label=lbl)
        scatter!(ax_polar, pd.cd, pd.cl; color, markersize=6)

        if pd.cmy !== nothing
            cmy_vals = Float64.(pd.cmy)
            if !all(isnan, cmy_vals)
                lines!(ax_cmy, pd.angle, cmy_vals; color, label=lbl)
                scatter!(ax_cmy, pd.angle, cmy_vals; color, markersize=6)
            end
        end
    end

    # Legend at (2,3) (use CL axis entries)
    Legend(fig[2, 3], ax_cl)

    return fig
end

# --- Your call site (updated angle_range) ---
PLOT && display(plot_polars_with_cmy(
    [solver],
    [body_aero],
    labels;
    literature_path_list=literature_paths,
    angle_range=range(-5, 40, step=1),   # alpha-axis data sweep
    angle_type="angle_of_attack",
    angle_of_attack=angle_of_attack_deg,
    side_slip=sideslip_deg,
    v_a=wind_speed,
    title="$(wing.n_panels)_panels_$(wing.spanwise_distribution)_from_yaml_settings",
    angle_xlim=(-5, 40)                  # enforce x-limits on angle plots
))

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
