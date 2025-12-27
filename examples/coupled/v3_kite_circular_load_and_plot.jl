lg = load_log("up_59_us_20_vw_15_2025_12_27_19_34"; path="processed_data/v3_kite")

fig = plot(sam_refine.sys_struct, syslog_refine;
           plot_turn_rates=true,
        #    plot_reelout=true,
           plot_tether=true,
           plot_aero_force=true,
        #    plot_aero_moment=true,
        #    plot_tether_moment=true,
        #    plot_twist=true,
           plot_aoa=true,
           plot_heading=true,
        #    plot_old_heading=true,
           plot_distance=true,
        #    plot_cone_angle=true,
           plot_elevation=true,
           plot_azimuth=true,
        #    plot_winch_force=true,
        #    plot_set_values=true,
           heading_setpoint=heading_setpoint_refine)
# display(fig)

scene = replay(lg, sam_refine.sys_struct; autoplay=false, loop=true)
display(scene)
# record(scene, "v3_kite_circular_load_and_plot.mp4"; fps=30, duration=20)  # Adjust duration as needed

