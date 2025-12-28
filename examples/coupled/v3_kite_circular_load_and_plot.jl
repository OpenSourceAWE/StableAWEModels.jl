using SymbolicAWEModels
using VortexStepMethod
using LinearAlgebra
using Statistics
using GLMakie
using KiteUtils
using DiscretePIDs
using Dates

log_name = "up_59_us_20_vw_15_2025_12_27_19_34"
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
lg = load_log(log_name; path="processed_data/v3_kite")

fig = plot(sam.sys_struct, lg;
           plot_turn_rates=true,
           plot_reelout=false,
        #    plot_tether=true,
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

## uncomment to show the 3D video replay
# scene = replay(lg, sam.sys_struct; autoplay=false, loop=true)
# display(scene)

##TODO: record does not work
# record(scene, "v3_kite_circular_load_and_plot.mp4"; fps=30, duration=20)  # Adjust duration as needed
