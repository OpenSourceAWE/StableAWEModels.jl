"""
CSV Replay Example

Reads flight test data from CSV and replays steering inputs through the
SymbolicAWEModel simulator. The kite is initialized to steady state, then
CSV steering commands are applied during simulation.

Usage:
    julia --project=examples examples/coupled/csv_replay.jl
"""

using SymbolicAWEModels, VortexStepMethod, KiteUtils
using SymbolicAWEModels: reposition!, rotate_around_z, rotate_around_y
using CSV, DataFrames, DiscretePIDs
using GLMakie
using Statistics
using Rotations
using UnPack
using LinearAlgebra
using OrdinaryDiffEqBDF
using KiteUtils: calc_elevation, azimuth_east

# Configuration parameters
CSV_PATH = "data/v3/2025-10-09_16-58-33_ProtoLogger_lidar.csv"
START_FRAME = 22068 + 1000        # First frame to replay
END_FRAME = START_FRAME + 200 # Last frame to replay (use nothing for all frames)
REMAKE_CACHE = false

# V3 Kite steering/depower calibration (from KCU documentation)
STEERING_L0 = 1.506  # Neutral steering tape length (m)
DEPOWER_L0 = 3.129

VIZ = true
PHYS = false

"""
    load_flight_data(csv_path::String)

Load and parse CSV flight data using space/multi-space delimiter.
Returns a named tuple with all CSV columns as time series data.
"""
function load_flight_data(csv_path::String)
    @info "Loading CSV data from: $csv_path"
    df = CSV.read(csv_path, DataFrame;
                  delim=' ',
                  silencewarnings=true,
                  normalizenames=true,
                  types=Float64,
                  strict=false)
    t0 = df.time[START_FRAME]
    df.time .= df.time .- t0
    col_names = Tuple(Symbol(name) for name in names(df))
    data = NamedTuple{col_names}(Tuple(eachcol(df)))
    return data
end

"""
    limit_frames(data; start_frame=1, end_frame=nothing)

Limit data to frame range [start_frame, end_frame].
Automatically slices all fields in the named tuple.
"""
function limit_frames(data; start_frame=1, end_frame=nothing)
    n_total = length(data.time)
    start_idx = max(1, min(start_frame, n_total))
    if isnothing(end_frame)
        end_idx = n_total
    else
        end_idx = max(start_idx, min(end_frame, n_total))
    end
    limited = NamedTuple{keys(data)}(
        Tuple(field[start_idx:end_idx] for field in data)
    )
    return limited
end

"""
    euler_to_quaternion(roll_deg, pitch_deg, yaw_deg)

Convert Euler angles (in degrees) from NED to ENU quaternion.
CSV data is in NED (North East Down) frame, but Q_b_w requires ENU (East North Up).

NED to ENU transformation:
  X_ENU = Y_NED (East)
  Y_ENU = X_NED (North)
  Z_ENU = -Z_NED (Up = -Down)
"""
function euler_to_quaternion(roll_deg, pitch_deg, yaw_deg)
    # Convert degrees to radians
    roll_rad = deg2rad(roll_deg)
    pitch_rad = deg2rad(pitch_deg)
    yaw_rad = deg2rad(yaw_deg)
    # Create rotation in NED frame using Rotations.jl (ZYX convention)
    rot_ned = RotZYX(yaw_rad, pitch_rad, roll_rad)
    R_ned_to_enu = [0.0 1.0 0.0;    # X_ENU = Y_NED (East)
                    1.0 0.0 0.0;    # Y_ENU = X_NED (North)
                    0.0 0.0 -1.0]   # Z_ENU = -Z_NED (Up)
    rot_enu = R_ned_to_enu * Matrix(rot_ned)
    q = SymbolicAWEModels.rotation_matrix_to_quaternion(rot_enu)
    return q
end

function calc_R_b_w(sys_struct::SystemStructure)
    @unpack points, wings, wind_vec_gnd = sys_struct
    wing = wings[1]
    R_b_w, origin = SymbolicAWEModels.calc_refine_wing_frame(
        points,
        wing.z_ref_points,
        wing.y_ref_points,
        wing.origin_idx
    )
    return R_b_w
end


function calc_heading(sys_struct::SystemStructure, R_b_w)
    e_x = R_b_w[:, 1]
    wind_norm = [1,0,0]
    # Project -e_x onto plane perpendicular to wind
    minus_e_x = -e_x
    proj_on_wind = dot(minus_e_x, wind_norm) * wind_norm
    e_x_perp = minus_e_x - proj_on_wind
    # Heading components in wind-perpendicular plane
    wind_cross_z = [wind_norm[2], -wind_norm[1], 0]
    heading_x = dot(e_x_perp, wind_cross_z)
    heading_z = e_x_perp[3]
    return atan(heading_x, heading_z)
end

function update_vel_from_csv!(sys, row)
    @unpack wings, points, winches, segments = sys
    wing = wings[1]

    # calc delta heading
    quat = euler_to_quaternion(row.roll,
                               row.pitch,
                               row.yaw)
    csv_heading =
        calc_heading(sys, SymbolicAWEModels.quaternion_to_rotation_matrix(quat)) + pi
    wing.R_b_w = calc_R_b_w(sys)
    curr_heading = calc_heading(sys, wing.R_b_w)
    delta_heading = csv_heading - curr_heading

    # calc needed heading and elevation
    csv_pos = [row.x, row.y, row.z]
    csv_elevation = calc_elevation(csv_pos)
    csv_azimuth = azimuth_east(csv_pos)
    curr_pos = wing.pos_w
    curr_elevation = calc_elevation(curr_pos)
    curr_azimuth = azimuth_east(curr_pos)
    delta_elevation = csv_elevation - curr_elevation
    delta_azimuth = csv_azimuth - curr_azimuth

    # apply heading and elevation
    wing.pos_w .= rotate_around_z(rotate_around_y(wing.pos_w, -delta_elevation),
                                    -delta_azimuth)
    for point in points
        point.pos_w .= rotate_around_z(rotate_around_y(point.pos_w, -delta_elevation),
                                        -delta_azimuth)
    end

    # apply vel
    csv_vel = [row.vx, row.vy, row.vz]
    points[wing.origin_idx].vel_w .= csv_vel

    # apply heading
    k = normalize(wing.pos_w)
    R_b_w = copy(wing.R_b_w)
    R_b_w[:, 1] = SymbolicAWEModels.rotate_v_around_k(R_b_w[:, 1], k, delta_heading)
    R_b_w[:, 2] = SymbolicAWEModels.rotate_v_around_k(R_b_w[:, 2], k, delta_heading)
    R_b_w[:, 3] = SymbolicAWEModels.rotate_v_around_k(R_b_w[:, 3], k, delta_heading)
    wing.R_b_w = R_b_w
    for point in points
        point.pos_w .= SymbolicAWEModels.rotate_v_around_k(point.pos_w, k, delta_heading)
    end

    # update tether length
    winches[1].tether_len = row.tether_len
    winches[1].tether_vel = row.tether_vel
    winches[1].brake = true
    segments[87].l0 = STEERING_L0 + row.steering / 100
    segments[89].l0 = STEERING_L0 - row.steering / 100
    segments[88].l0 = DEPOWER_L0 + row.depower / 100
    nothing
end

"""
    update_sys_struct_from_csv!(sys_struct, csv_row_data)

Update system structure from a single CSV row.
Updates wing orientation from Euler angles and position via transform system.
"""
function update_sys_struct_from_csv!(sys, row)
    @unpack wings, points, winches, segments = sys
    wing = wings[1]

    # calc delta heading
    quat = euler_to_quaternion(row.roll,
                               row.pitch,
                               row.yaw)
    csv_heading =
        calc_heading(sys, SymbolicAWEModels.quaternion_to_rotation_matrix(quat)) + pi
    wing.R_b_w = calc_R_b_w(sys)
    curr_heading = calc_heading(sys, wing.R_b_w)
    delta_heading = csv_heading - curr_heading

    # calc needed transform
    csv_pos = [row.x, row.y, row.z]
    curr_pos = wing.pos_w
    delta_pos = csv_pos - curr_pos
    # apply transform
    wing.pos_w .+= delta_pos
    for point in points
        transform_frac = point.pos_w ⋅ normalize(wing.pos_w) / norm(wing.pos_w)
        point.pos_w .+= transform_frac * delta_pos
    end

    # apply vel
    csv_vel = [row.vx, row.vy, row.vz]
    wing.vel_w .= csv_vel
    for point in points
        transform_frac = point.pos_w ⋅ normalize(wing.pos_w) / norm(wing.pos_w)
        point.vel_w .= transform_frac * csv_vel
    end

    # apply heading
    k = normalize(wing.pos_w)
    R_b_w = copy(wing.R_b_w)
    R_b_w[:, 1] = SymbolicAWEModels.rotate_v_around_k(R_b_w[:, 1], k, delta_heading)
    R_b_w[:, 2] = SymbolicAWEModels.rotate_v_around_k(R_b_w[:, 2], k, delta_heading)
    R_b_w[:, 3] = SymbolicAWEModels.rotate_v_around_k(R_b_w[:, 3], k, delta_heading)
    wing.R_b_w = R_b_w
    for point in points
        point.pos_w .= SymbolicAWEModels.rotate_v_around_k(point.pos_w, k, delta_heading)
    end

    # update tether length
    winches[1].tether_len = row.tether_len
    winches[1].tether_vel = row.tether_vel
    winches[1].brake = true
    segments[87].l0 = STEERING_L0 + row.steering / 100
    segments[89].l0 = STEERING_L0 - row.steering / 100
    segments[88].l0 = DEPOWER_L0 + row.depower / 100
    nothing
end

"""
    run_csv_replay(csv_path; kwargs...)

Main function to replay CSV data (visualization only, no physics).
"""
function run_csv_replay(csv_path::String;
                        start_frame=START_FRAME,
                        end_frame=END_FRAME)

    @info "="^60
    @info "CSV Data Replay (Visualization Only)"
    @info "="^60

    # Load and preprocess CSV data
    raw_data = load_flight_data(csv_path)
    csv_data = limit_frames(raw_data; start_frame, end_frame)

    # Load v3 kite system structure (without physics model)
    @info "Loading v3 kite system structure from YAML"
    set_data_path("data/v3")
    set = Settings("system.yaml")

    vsm_set_path = joinpath(get_data_path(), "vsm_settings.yaml")
    vsm_set = VortexStepMethod.VSMSettings(vsm_set_path; data_prefix=false)

    sys_struct = load_sys_struct_from_yaml("data/v3/struc_geometry.yaml";
        system_name="v3", set, wing_type=REFINE, vsm_set)
    sam = SymbolicAWEModel(set, sys_struct)
    init!(sam)

    # Create logger to store CSV replay data
    n_points = length(csv_data.time)
    @info "Creating log with $n_points timesteps"

    # Create SysState for logging
    sys_state = SysState(sam)

    # Create logger
    logger = Logger(sam, n_points)

    # Loop through CSV data and update sys_struct
    @info "Replaying CSV data..."
    replay_start = time()

    for step in 1:n_points
        # Create row data structure
        csv_row = (
            time = csv_data.time[step],
            roll = csv_data.kite_0_roll[step],
            pitch = csv_data.kite_0_pitch[step],
            yaw = csv_data.kite_0_yaw[step],
            x = csv_data.kite_pos_east[step],
            y = csv_data.kite_pos_north[step],
            z = csv_data.kite_height[step],
            vx = csv_data.kite_est_vx[step],
            vy = csv_data.kite_est_vy[step],
            vz = csv_data.kite_est_vz[step],
            tether_len = csv_data.ground_tether_length[step],
            tether_vel = csv_data.ground_tether_reelout_speed[step],
            steering = csv_data.kite_actual_steering[step],
            depower = csv_data.kite_actual_depower[step],
        )

        # Update system structure from CSV
        if step==1
            update_sys_struct_from_csv!(sam.sys_struct, csv_row)
        else
            update_vel_from_csv!(sam.sys_struct, csv_row)
        end
        SymbolicAWEModels.reinit!(sam, sam.prob, FBDF())

        # Update sys_state from sam and log it
        update_sys_state!(sys_state, sam)
        sys_state.time = csv_row.time
        log!(logger, sys_state)

        # Progress reporting
        if step % max(1, div(n_points, 10)) == 0 || step == n_points
            elapsed = time() - replay_start
            @info "  Logged $step/$n_points points " *
                  "(t = $(round(csv_row.time, digits=2)) s)"
        end
    end

    replay_elapsed = time() - replay_start
    @info "CSV replay logged in $(round(replay_elapsed, digits=2)) s"

    # Save log
    @info "Saving replay log..."
    save_log(logger, "csv_replay")
    syslog = load_log("csv_replay")

    return sam, syslog, csv_data
end

"""
    run_physics_replay(csv_path; kwargs...)

Replay CSV data using physics simulation with CSV inputs.
Updates tether length and steering/depower segments from CSV data at each step.
"""
function run_physics_replay(csv_path::String;
                             start_frame=START_FRAME,
                             end_frame=END_FRAME,
                             initial_damping=100.0,
                             decay_time=2.0)

    @info "="^60
    @info "CSV Physics Simulation Replay"
    @info "="^60

    # Load and preprocess CSV data
    raw_data = load_flight_data(csv_path)
    csv_data = limit_frames(raw_data; start_frame, end_frame)

    # Load v3 kite system structure with physics
    @info "Loading v3 kite system structure with physics from YAML"
    set_data_path("data/v3")
    set = Settings("system.yaml")

    # Update wind settings from CSV initial values
    set.v_wind = csv_data.ground_wind_velocity[1]
    set.upwind_dir = csv_data.est_upwind_direction[1]

    vsm_set_path = joinpath(get_data_path(), "vsm_settings.yaml")
    vsm_set = VortexStepMethod.VSMSettings(vsm_set_path; data_prefix=false)

    sys_struct = load_sys_struct_from_yaml("data/v3/struc_geometry.yaml";
        system_name="v3", set, wing_type=REFINE, vsm_set)

    # Create SymbolicAWEModel with physics
    @info "Creating SymbolicAWEModel with physics..."
    sam = SymbolicAWEModel(set, sys_struct)

    # Set winch to brake mode (lengths controlled directly)
    if !isempty(sam.sys_struct.winches)
        sam.sys_struct.winches[1].brake = true
    end

    # Initialize with initial CSV position
    @info "Initializing with CSV initial conditions..."
    transform = sam.sys_struct.transforms[1]
    transform.elevation = csv_data.kite_elevation[1]
    transform.azimuth = csv_data.kite_azimuth[1]

    # Calculate initial heading from CSV orientation
    quat = euler_to_quaternion(csv_data.kite_0_roll[1],
                               csv_data.kite_0_pitch[1],
                               csv_data.kite_0_yaw[1])
    sam.sys_struct.wings[1].Q_b_w .= quat
    R_b_w = calc_R_b_w(sam.sys_struct)
    csv_heading = calc_heading(sam.sys_struct, R_b_w)
    transform.heading = csv_heading + pi

    # Reinit with transform
    SymbolicAWEModels.reinit!(sam.sys_struct.transforms, sam.sys_struct)

    # Initialize the model
    init!(sam)

    # Set initial damping to help with stability at startup
    @info "Setting initial world frame damping: $initial_damping N·s/m (decay over $decay_time s)"
    SymbolicAWEModels.set_world_frame_damping(sam.sys_struct, initial_damping)

    # Create logger
    n_points = length(csv_data.time)
    @info "Creating log with $n_points timesteps"
    sys_state = SysState(sam)
    logger = Logger(sam, n_points)

    # Physics simulation loop at 100 fps
    @info "Running physics simulation at 100 fps with CSV inputs..."
    sim_start = time()

    SIM_DT = 0.01  # 100 fps = 0.01 second timesteps
    sim_time = 0.0
    csv_idx = 1
    sim_step = 0

    # Calculate total simulation steps
    total_sim_time = csv_data.time[end] - csv_data.time[1]
    total_steps = round(Int, total_sim_time / SIM_DT)

    while sim_time <= total_sim_time && csv_idx <= n_points
        sim_step += 1

        # Find the CSV frame that corresponds to current simulation time
        # Interpolate if needed, or use nearest
        target_time = csv_data.time[1] + sim_time
        while csv_idx < n_points && csv_data.time[csv_idx] < target_time
            csv_idx += 1
        end

        # Get CSV inputs for this timestep
        tether_length = csv_data.ground_tether_length[csv_idx]
        steering_cm = csv_data.kite_actual_steering[csv_idx]
        depower_cm = csv_data.kite_actual_depower[csv_idx]

        # Update tether length
        @show sam.sys_struct.winches[1].tether_len - tether_length
        # sam.sys_struct.winches[1].tether_len = tether_length

        # Update steering segments (segments 87 and 89: left and right steering tapes)
        # CSV steering is in cm, need to convert to meters and apply to differential
        steering_m = steering_cm / 100.0  # Convert cm to m

        # Steering differential: positive = right turn, negative = left turn
        # Left steering tape (segment 87, points 37->1)
        left_steering_len = L0_NEUTRAL - steering_m / 2
        @show left_steering_len
        # Right steering tape (segment 89, points 38->1)
        right_steering_len = L0_NEUTRAL + steering_m / 2

        # sam.sys_struct.segments[87].l0 = left_steering_len
        # sam.sys_struct.segments[89].l0 = right_steering_len

        # Update depower segment (segment 88: center power tape, points 35->1)
        # CSV depower is in cm
        depower_m = depower_cm / 100.0  # Convert cm to m
        # Assume nominal depower length and add depower adjustment
        # sam.sys_struct.segments[88].l0 = 3.129 + depower_m

        # Update damping (decay from initial_damping to 0 over decay_time)
        if sim_time <= decay_time
            current_damping = initial_damping * (1.0 - sim_time / decay_time)
            SymbolicAWEModels.set_world_frame_damping(sam.sys_struct, current_damping)
        elseif sim_time <= decay_time + SIM_DT
            # Ensure we set it to exactly 0 after decay period
            SymbolicAWEModels.set_world_frame_damping(sam.sys_struct, 0.0)
        end

        # Step the simulation at 100 fps
        set_values = zeros(3)  # No active torque control (brake mode)

        # Try to step, catch errors and continue
        try
            next_step!(sam; set_values, dt=SIM_DT)
            update_sys_state!(sys_state, sam)
            sys_state.time = sim_time
            log!(logger, sys_state)
        catch e
            @warn "Simulation became unstable at step $sim_step (t=$(round(sim_time, digits=2))s)" exception=e
            @info "Stopping simulation early and generating visualization with collected data..."
            break
        end

        # Progress reporting
        if sim_step % max(1, div(total_steps, 10)) == 0 || sim_step == total_steps
            elapsed = time() - sim_start
            @info "  Simulated $sim_step/$total_steps steps " *
                  "(t = $(round(sim_time, digits=2)) s, " *
                  "elapsed = $(round(elapsed, digits=2)) s)"
        end

        sim_time += SIM_DT
    end

    sim_elapsed = time() - sim_start
    @info "Physics simulation completed in $(round(sim_elapsed, digits=2)) s"

    @info "Saving physics replay log..."
    save_log(logger, "csv_physics_replay")
    syslog = load_log("csv_physics_replay")

    return sam, syslog, csv_data
end

# Main execution
if VIZ
    @info "Running CSV data replay (visualization only)..."
    sam_viz, syslog_viz, csv_data = run_csv_replay(CSV_PATH)
end
if PHYS
    @info "\nRunning physics simulation replay..."
    sam_phys, syslog_phys, csv_data_phys, fig_phys = run_physics_replay(CSV_PATH)
end

if VIZ && PHYS
    plot([sam_viz.sys_struct, sam_phys.sys_struct], [syslog_viz, syslog_phys]; plot_elevation=true)
elseif VIZ
    plot(sam_viz.sys_struct, syslog_viz)
elseif PHYS
    plot(sam_phys.sys_struct, syslog_phys)
end

