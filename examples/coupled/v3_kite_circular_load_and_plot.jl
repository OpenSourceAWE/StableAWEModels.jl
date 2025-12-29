using SymbolicAWEModels
using VortexStepMethod
using LinearAlgebra
using Statistics
using GLMakie
using KiteUtils
using DiscretePIDs
using Dates

function load_log_and_system(; log_name::String)
    up = parse(Float64, log_name[4:5])
    us = parse(Float64, log_name[10:11])
    v_wind = parse(Float64, log_name[16:17])
    initial_damping = 100.0

    # Load settings
    wing_type = SymbolicAWEModels.REFINE
    wing_type_str = "REFINE"
    @info "Running v3 kite simulation with REFINE wing type..."

    set_data_path("data/v3")
    set = Settings("system.yaml")
    set.v_wind = v_wind
    set.upwind_dir = -90.0

    # Load YAML structure path
    model_name = "v3_refine"
    struc_yaml_path = joinpath("data", "v3", "struc_geometry.yaml")

    # Load VSMSettings
    vsm_set_path = joinpath(get_data_path(), "vsm_settings.yaml")
    vsm_set = VortexStepMethod.VSMSettings(vsm_set_path; data_prefix=false)

    # Use 36 panels for both wing types (matches vsm_settings.yaml default)
    vsm_set.wings[1].n_panels = 36
    # Note: n_unrefined_sections is automatically inferred from YAML geometry

    # Load system structure with wing_type and vsm_set parameters
    sys = load_sys_struct_from_yaml(struc_yaml_path;
        system_name=model_name, set, wing_type, vsm_set)

    # Initialize damping
    SymbolicAWEModels.set_world_frame_damping(sys, initial_damping)

    wing_points = [p for p in sys.points if p.type == WING]
    n_unrefined = sys.wings[1].vsm_wing.n_unrefined_sections
    @info "REFINE wing setup:" n_wing_points=length(wing_points) n_groups=length(sys.groups) n_unrefined=n_unrefined n_panels=length(sys.wings[1].vsm_aero.panels) n_segments=length(sys.segments)

    # Create symbolic model
    sam = SymbolicAWEModel(set, sys)
    lg = KiteUtils.load_log(log_name; path="processed_data/v3_kite")
    return lg, sam, up, us, v_wind
end

function print_and_plot_wing(lg, sam)

   # Grab the last recorded SysState from the log
   lg_last = lg.syslog[end]

   wing = sam.sys_struct.wings[1]
   origin_idx = wing.origin_idx
   origin_w = [
      lg_last.X[origin_idx],
      lg_last.Y[origin_idx],
      lg_last.Z[origin_idx],
   ]
   R_b_w = SymbolicAWEModels.quaternion_to_rotation_matrix(lg_last.orient)
   wing_point_idxs = [p.idx for p in sam.sys_struct.points if p.type == SymbolicAWEModels.WING]

   println("Wing node positions (world frame):")
   for idx in wing_point_idxs
      println("$(lg_last.X[idx]), $(lg_last.Y[idx]), $(lg_last.Z[idx])")
   end

   println("\nWing node positions (body frame):")
   for idx in wing_point_idxs
      pos_w = [
         lg_last.X[idx],
         lg_last.Y[idx],
         lg_last.Z[idx],
      ]
      pos_b = R_b_w' * (pos_w .- origin_w)
      println("$(pos_b[1]), $(pos_b[2]), $(pos_b[3])")
   end

   # 2D scatter plots of wing nodes in body frame
   xs_b = Float64[]
   ys_b = Float64[]
   zs_b = Float64[]
   for idx in wing_point_idxs
      pos_w = [
         lg_last.X[idx],
         lg_last.Y[idx],
         lg_last.Z[idx],
      ]
      pos_b = R_b_w' * (pos_w .- origin_w)
      push!(xs_b, pos_b[1]); push!(ys_b, pos_b[2]); push!(zs_b, pos_b[3])
   end

   fig2 = Figure(resolution=(900, 300))
   ax_xy = Axis(fig2[1, 1], title="Top (x,y)", xlabel="y_b", ylabel="x_b")
   ax_xz = Axis(fig2[1, 2], title="Side (x,z)", xlabel="x_b", ylabel="z_b")
   ax_yz = Axis(fig2[1, 3], title="Front (y,z)", xlabel="y_b", ylabel="z_b")

   ax_xy.aspect = DataAspect()
   ax_xz.aspect = DataAspect()
   ax_yz.aspect = DataAspect()

   xs_b_rot = ys_b
   ys_b_rot = [-x for x in xs_b]

   scatter!(ax_xy, xs_b_rot, ys_b_rot)
   scatter!(ax_xz, xs_b, zs_b)
   scatter!(ax_yz, ys_b, zs_b)

   # enforce equal scale/step size across all plots (square span)
   spans = [
      maximum(xs_b_rot) - minimum(xs_b_rot),
      maximum(ys_b_rot) - minimum(ys_b_rot),
      maximum(xs_b) - minimum(xs_b),
      maximum(zs_b) - minimum(zs_b),
      maximum(ys_b) - minimum(ys_b),
   ]
   global_span = maximum(spans)
   global_span = global_span > 0 ? global_span : 1.0
   global_span *= 1.05  # slight padding

   function set_square_limits!(ax, xs, ys, span)
      cx = mean(xs)
      cy = mean(ys)
      half = span / 2
      xlims!(ax, cx - half, cx + half)
      ylims!(ax, cy - half*0.5, cy + half*0.5)
   end

   set_square_limits!(ax_xy, xs_b_rot, ys_b_rot, global_span)
   set_square_limits!(ax_xz, xs_b, zs_b, global_span)
   set_square_limits!(ax_yz, ys_b, zs_b, global_span)

   # Draw wing structural segments for the plotted wing (indices 1:19 in YAML)
   wing_seg_idxs = 1:19  # leading edge + struts
   lines_xy = Point2f[]
   lines_xz = Point2f[]
   lines_yz = Point2f[]
   for seg_idx in wing_seg_idxs
      seg = sam.sys_struct.segments[seg_idx]
      p1, p2 = seg.point_idxs
      pos1_w = [lg_last.X[p1], lg_last.Y[p1], lg_last.Z[p1]]
      pos2_w = [lg_last.X[p2], lg_last.Y[p2], lg_last.Z[p2]]
      pos1_b = R_b_w' * (pos1_w .- origin_w)
      pos2_b = R_b_w' * (pos2_w .- origin_w)
      x1, y1, z1 = pos1_b
      x2, y2, z2 = pos2_b
      # rotate only the xy projection by 90 deg clockwise: (x,y)->(y,-x)
      push!(lines_xy, Point2f(y1, -x1), Point2f(y2, -x2), Point2f(NaN, NaN)) # NaN breaks segments
      push!(lines_xz, Point2f(x1, z1), Point2f(x2, z2), Point2f(NaN, NaN))
      push!(lines_yz, Point2f(y1, z1), Point2f(y2, z2), Point2f(NaN, NaN))
   end

   lines!(ax_xy, lines_xy, color=:gray)
   lines!(ax_xz, lines_xz, color=:gray)
   lines!(ax_yz, lines_yz, color=:gray)

   fig2
   display(fig2)
end

function plot_time_series(lg, sam)
    fig = plot(sam.sys_struct, lg;
            plot_turn_rates=true,
            plot_reelout=false,
         #    plot_tether=true,
            plot_gk=true,
         #    plot_aero_force=true,
         #    plot_aero_moment=true,
         #    plot_tether_moment=true,
         #    plot_twist=true,
            plot_aoa=true,
            plot_heading=false,
         #    plot_old_heading=true,
         #    plot_distance=true,
         #    plot_cone_angle=true,
            plot_elevation=true,
            plot_azimuth=true,
            plot_winch_force=false,
            plot_set_values=false)
    display(fig)
end


log_name = "up_40_us_25_vw_15_date_2025_12_29_18_31"
lg, sam, up, us, v_wind = load_log_and_system(log_name=log_name)

### plot time series
# plot_time_series(lg, sam)

## show 3D animation
scene = replay(lg, sam.sys_struct; autoplay=false, loop=true)
display(scene)

### show 2D wing node plots
# print_and_plot_wing(lg, sam)

##TODO: record does not work
# record(scene, "v3_kite_circular_load_and_plot.mp4"; fps=30, duration=20)  # Adjust duration as needed
