# SPDX-FileCopyrightText: 2025 Bart van de Lint, Uwe Fechner
#
# SPDX-License-Identifier: MIT

module SymbolicAWEModelsControlPlotsExt
using ControlPlots, LaTeXStrings, KiteUtils, SymbolicAWEModels

export plot

function ControlPlots.plot(sys::SystemStructure, lg::SysLog)
    sl = lg.syslog
    turn_rates_deg = rad2deg.(hcat(sl.turn_rates...))
    v_reelout_23 = [sl.v_reelout[i][2] for i in eachindex(sl.v_reelout)], [sl.v_reelout[i][3] for i in eachindex(sl.v_reelout)] # Winch 2 and 3
    aero_force_z = [sl.aero_force_b[i][3] for i in eachindex(sl.aero_force_b)]
    aero_moment_z = [sl.aero_moment_b[i][3] for i in eachindex(sl.aero_moment_b)]
    twist_angles_deg = rad2deg.(hcat(sl.twist_angles...))[eachindex(sys.groups),:]
    AoA_deg = rad2deg.(sl.AoA)
    heading_deg = rad2deg.(sl.heading)

    twist_labels = ["twist[$i]" for i in eachindex(sys.groups)]

    ControlPlots.plotx(sl.time,
        [turn_rates_deg[1,:], turn_rates_deg[2,:], turn_rates_deg[3,:]],
        v_reelout_23,
        [aero_force_z, aero_moment_z],
        [twist_angles_deg[i,:] for i in eachindex(sys.groups)],
        [AoA_deg],
        [heading_deg];
        ylabels=["turn rates [°/s]", L"v_{ro}~[m/s]", "aero F/M", "twist [°]", "AoA [°]", "heading [°]"],
        ysize=10,
        labels=[
            [L"\omega_x", L"\omega_y", L"\omega_z"],
            ["v_ro[2]", "v_ro[3]"],
            [L"F_{aero,z}", L"M_{aero,z}"],
            twist_labels,
            ["AoA"],
            ["heading"]
        ],
        fig="Oscillating Steering Input Response")
end

function ControlPlots.plot(sys::SystemStructure, reltime; l_tether=50.0, wing_pos=nothing, e_z=zeros(3), zoom=false, front=false, xy=nothing)
    pos = [sys.points[i].pos_w for i in eachindex(sys.points)]
    !isnothing(wing_pos) && (pos = [pos..., wing_pos...])
    seg = [[sys.segments[i].point_idxs[1], sys.segments[i].point_idxs[2]] for i in eachindex(sys.segments)]
    if zoom && !front
        xlim = (pos[end][1] - 6, pos[end][1]+6)
        ylim = (pos[end][3] - 10, pos[end][3]+2)
    elseif zoom && front
        xlim = (pos[end][2] - 6, pos[end][2]+6)
        ylim = (pos[end][3] - 10, pos[end][3]+2)
    elseif !zoom && !front
        xlim = (-0.1l_tether, 1.2l_tether)
        ylim = (-0.1l_tether, 1.2l_tether)
    elseif !zoom && front
        xlim = (-0.75l_tether, 0.75l_tether)
        ylim = (-0.1l_tether, 1.2l_tether)
    end
    ControlPlots.plot2d(pos, seg, reltime; zoom, front, xlim, ylim, dz_zoom=0.6, xy)
end

function ControlPlots.plot(s::SymbolicAWEModel, reltime; kwargs...)
    wings = s.sys_struct.wings
    pos = s.integrator[s.sys.pos]
    if length(wings) > 0
        wing_pos = [s.integrator[s.sys.wing_pos[i, :]] for i in eachindex(wings)]
        e_z = [s.integrator[s.sys.e_z[i, :]] for i in eachindex(wings)]
    else
        wing_pos = nothing
        e_z = zeros(3)
    end
        
    for (i, point) in enumerate(s.sys_struct.points)
        point.pos_w .= pos[:, i]
    end
    plot(s.sys_struct, reltime; s.set.l_tether, wing_pos, e_z, kwargs...)
end

end
