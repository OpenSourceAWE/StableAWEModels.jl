# SPDX-FileCopyrightText: 2025 Bart van de Lint <bart@vandelint.net>
# SPDX-License-Identifier: MPL-2.0

using YAML

mutable struct Settings
    dict::Vector{Dict}
    sim_settings::String
    log_file::String
    log_level::Int
    time_lapse::Float64
    sim_time::Float64
    segments::Int64
    sample_freq::Int64
    zoom::Float64
    kite_scale::Float64
    fixed_font::String
    elevations::Vector{Float64}
    elevation_rates::Vector{Float64}
    azimuths::Vector{Float64}
    azimuth_rates::Vector{Float64}
    headings::Vector{Float64}
    heading_rates::Vector{Float64}
    l_tethers::Vector{Float64}
    kite_distances::Vector{Float64}
    v_reel_outs::Vector{Float64}
    depowers::Vector{Float64}
    steerings::Vector{Float64}
    abs_tol::Float64
    rel_tol::Float64
    solver::String
    linear_solver::String
    max_order::Int64
    max_iter::Int64
    relaxation::Float64
    c0::Float64
    c_s::Float64
    c2_cor::Float64
    k_ds::Float64
    delta_st::Float64
    max_steering::Float64
    cs_4p::Float64
    alpha_d_max::Float64
    depower_offset::Float64
    model::String
    physical_model::String
    version::Int64
    mass::Float64
    area::Float64
    rel_side_area::Float64
    height_k::Float64
    alpha_cl::Vector{Float64}
    cl_list::Vector{Float64}
    alpha_cd::Vector{Float64}
    cd_list::Vector{Float64}
    cms::Float64
    width::Float64
    alpha_zero::Float64
    alpha_ztip::Float64
    m_k::Float64
    rel_nose_mass::Float64
    rel_top_mass::Float64
    smc::Float64
    cmq::Float64
    cord_length::Float64
    c_spring_kite::Float64
    damping_kite_springs::Float64
    rel_mass_p2::Float64
    rel_mass_p3::Float64
    rel_mass_p4::Float64
    foil_file::String
    top_bridle_points::Vector{Vector{Float64}}
    bridle_tether_diameter::Float64
    power_tether_diameter::Float64
    steering_tether_diameter::Float64
    crease_frac::Float64
    bridle_fracs::Vector{Float64}
    fixed_index::Int64
    quasi_static::Bool
    d_line::Float64
    h_bridle::Float64
    l_bridle::Float64
    rel_compr_stiffness::Float64
    rel_damping::Float64
    kcu_model::String
    kcu_mass::Float64
    kcu_diameter::Float64
    cd_kcu::Float64
    depower_zero::Float64
    degrees_per_percent_power::Float64
    power2steer_dist::Float64
    depower_drum_diameter::Float64
    tape_thickness::Float64
    v_depower::Float64
    v_steering::Float64
    depower_gain::Float64
    steering_gain::Float64
    d_tether::Float64
    cd_tether::Float64
    damping::Float64
    c_spring::Float64
    rho_tether::Float64
    e_tether::Float64
    winch_model::String
    max_force::Float64
    v_ro_max::Float64
    v_ro_min::Float64
    max_acc::Float64
    drum_radius::Float64
    gear_ratio::Float64
    inertia_total::Float64
    f_coulomb::Float64
    c_vf::Float64
    p_speed::Float64
    i_speed::Float64
    v_wind::Float64
    upwind_dir::Float64
    temp_ref::Float64
    height_gnd::Float64
    h_ref::Float64
    rho_0::Float64
    alpha::Float64
    z0::Float64
    profile_law::Int64
    use_turbulence::Float64
    v_wind_gnds::Vector{Float64}
    avg_height::Float64
    rel_turbs::Vector{Float64}
    i_ref::Float64
    v_ref::Float64
    grid::Vector{Int64}
    height_step::Float64
    grid_step::Float64
    g_earth::Float64
end

function local_load_settings(yaml_path::AbstractString)
    if !isfile(yaml_path)
        error("YAML file not found: $yaml_path")
    end
    dict = YAML.load_file(yaml_path)
    if isa(dict, Dict) && length(dict) == 1 && isa(first(values(dict)), Dict)
        dict = first(values(dict))
    end
    flat = Dict{Symbol,Any}()
    for v in values(dict)
        if isa(v, Dict)
            for (kk, vv) in v
                flat[Symbol(kk)] = vv
            end
        end
    end
    return Settings(
        [Dict()],
        get(flat, :sim_settings, ""),
        get(flat, :log_file, ""),
        get(flat, :log_level, 2),
        get(flat, :time_lapse, 0.0),
        get(flat, :sim_time, 0.0),
        get(flat, :segments, 0),
        get(flat, :sample_freq, 0),
        get(flat, :zoom, 0.0),
        get(flat, :kite_scale, 1.0),
        get(flat, :fixed_font, ""),
        get(flat, :elevations, [70.0]),
        get(flat, :elevation_rates, [0.0]),
        get(flat, :azimuths, [0.0]),
        get(flat, :azimuth_rates, [0.0]),
        get(flat, :headings, [0.0]),
        get(flat, :heading_rates, [0.0]),
        get(flat, :l_tethers, [0.0]),
        get(flat, :kite_distances, [0.0]),
        get(flat, :v_reel_outs, [0.0]),
        get(flat, :depowers, [0.0]),
        get(flat, :steerings, [0.0]),
        get(flat, :abs_tol, 0.0),
        get(flat, :rel_tol, 0.0),
        get(flat, :solver, "DFBDF"),
        get(flat, :linear_solver, "GMRES"),
        get(flat, :max_order, 4),
        get(flat, :max_iter, 1),
        get(flat, :relaxation, 0.0),
        get(flat, :c0, 0.0),
        get(flat, :c_s, 0.0),
        get(flat, :c2_cor, 0.0),
        get(flat, :k_ds, 0.0),
        get(flat, :delta_st, 0.0),
        get(flat, :max_steering, 0.0),
        get(flat, :cs_4p, 1.0),
        get(flat, :alpha_d_max, 0.0),
        get(flat, :depower_offset, 23.6),
        get(flat, :model, "data/kite.obj"),
        get(flat, :physical_model, ""),
        get(flat, :version, 1),
        get(flat, :mass, 0.0),
        get(flat, :area, 0.0),
        get(flat, :rel_side_area, 0.0),
        get(flat, :height_k, 0.0),
        get(flat, :alpha_cl, Float64[]),
        get(flat, :cl_list, Float64[]),
        get(flat, :alpha_cd, Float64[]),
        get(flat, :cd_list, Float64[]),
        get(flat, :cms, 0.0),
        get(flat, :width, 0.0),
        get(flat, :alpha_zero, 0.0),
        get(flat, :alpha_ztip, 0.0),
        get(flat, :m_k, 0.0),
        get(flat, :rel_nose_mass, 0.0),
        get(flat, :rel_top_mass, 0.0),
        get(flat, :smc, 0.0),
        get(flat, :cmq, 0.0),
        get(flat, :cord_length, 0.0),
        get(flat, :c_spring_kite, 0.0),
        get(flat, :damping_kite_springs, 0.0),
        get(flat, :rel_mass_p2, 0.0),
        get(flat, :rel_mass_p3, 0.0),
        get(flat, :rel_mass_p4, 0.0),
        get(flat, :foil_file, "data/ram_air_kite_foil.dat"),
        get(flat, :top_bridle_points, [[0.290199, 0.784697, -2.61305], [0.392683, 0.785271, -2.61201], [0.498202, 0.786175, -2.62148], [0.535543, 0.786175, -2.62148]]),
        get(flat, :bridle_tether_diameter, 2.0),
        get(flat, :power_tether_diameter, 2.0),
        get(flat, :steering_tether_diameter, 1.0),
        get(flat, :crease_frac, 0.82),
        get(flat, :bridle_fracs, [0.088, 0.31, 0.58, 0.93]),
        get(flat, :fixed_index, 1),
        get(flat, :quasi_static, false),
        get(flat, :d_line, 0.0),
        get(flat, :h_bridle, 0.0),
        get(flat, :l_bridle, 0.0),
        get(flat, :rel_compr_stiffness, 0.0),
        get(flat, :rel_damping, 0.0),
        get(flat, :kcu_model, "KCU1"),
        get(flat, :kcu_mass, 0.0),
        get(flat, :kcu_diameter, 0.0),
        get(flat, :cd_kcu, 0.0),
        get(flat, :depower_zero, 0.0),
        get(flat, :degrees_per_percent_power, 0.0),
        get(flat, :power2steer_dist, 0.0),
        get(flat, :depower_drum_diameter, 0.0),
        get(flat, :tape_thickness, 0.0),
        get(flat, :v_depower, 0.0),
        get(flat, :v_steering, 0.0),
        get(flat, :depower_gain, 3.0),
        get(flat, :steering_gain, 3.0),
        get(flat, :d_tether, 0.0),
        get(flat, :cd_tether, 0.0),
        get(flat, :damping, 0.0),
        get(flat, :c_spring, 0.0),
        get(flat, :rho_tether, 0.0),
        get(flat, :e_tether, 0.0),
        get(flat, :winch_model, ""),
        get(flat, :max_force, 4000.0),
        get(flat, :v_ro_max, 8.0),
        get(flat, :v_ro_min, -8.0),
        get(flat, :max_acc, 0.0),
        get(flat, :drum_radius, 0.1615),
        get(flat, :gear_ratio, 6.2),
        get(flat, :inertia_total, 0.0),
        get(flat, :f_coulomb, 122.0),
        get(flat, :c_vf, 30.6),
        get(flat, :p_speed, 0.0),
        get(flat, :i_speed, 0.0),
        get(flat, :v_wind, 0.0),
        get(flat, :upwind_dir, 0.0),
        get(flat, :temp_ref, 0.0),
        get(flat, :height_gnd, 0.0),
        get(flat, :h_ref, 0.0),
        get(flat, :rho_0, 0.0),
        get(flat, :alpha, 0.0),
        get(flat, :z0, 0.0),
        get(flat, :profile_law, 0),
        get(flat, :use_turbulence, 0.0),
        get(flat, :v_wind_gnds, Float64[]),
        get(flat, :avg_height, 0.0),
        get(flat, :rel_turbs, Float64[]),
        get(flat, :i_ref, 0.0),
        get(flat, :v_ref, 0.0),
        get(flat, :grid, Int[]),
        get(flat, :height_step, 0.0),
        get(flat, :grid_step, 0.0),

        get(flat, :g_earth, 0.0)
    )
end

