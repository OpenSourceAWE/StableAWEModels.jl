"""
CSV Replay Example

Reads flight test data from CSV and replays steering inputs through the
SymbolicAWEModel simulator. The kite is initialized to steady state, then
CSV steering commands are applied during simulation.

Usage:
    julia --project=examples examples/coupled/csv_replay.jl
"""

using SymbolicAWEModels, VortexStepMethod, KiteUtils
using SymbolicAWEModels: reposition!, rotate_around_z, rotate_around_y, calc_steady_torque
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
START_FRAME = 22068 + 1000 # First frame to replay
END_FRAME = START_FRAME + 50 # Last frame to replay (use nothing for all frames)
REMAKE_CACHE = false

# V3 Kite steering/depower calibration (from KCU documentation)
# Steering calibration
STEERING_L0 = 1.506  # Neutral steering tape length (m)
STEERING_GAIN = 1.2  # Maximum differential (m) at |u_s| = 1

# Depower calibration
DEPOWER_L0 = 0.2 # SUPPOSED TO BE 0.2
DEPOWER_GAIN = 5.0

INITIAL_DAMPING = 400.0
DECAY_TIME = 2.0
MIN_DAMPING = 300.0

# PID controller parameters for heading control
HEADING_KP = 0.1    # Low proportional gain
HEADING_TAU_I = 0.1  # Integral time constant (seconds)
HEADING_KD = 0.0     # No derivative gain
DT_CONTROL = 0.001    # Control timestep

# PI controller parameters for winch length control
WINCH_LENGTH_KP = 0.0      # Low proportional gain for length tracking
WINCH_LENGTH_TAU_I = 10.0   # Integral time constant
WINCH_LENGTH_KD = 0.0      # No derivative gain


"""
    steering_percentage_to_lengths(percentage)

Convert CSV steering percentage to left/right tape lengths (m).
Percentage convention: negative = left turn, positive = right turn.
"""
function steering_percentage_to_lengths(percentage)
    u_s = percentage / 100.0  # Convert percentage to [-1, 1]
    L_left = STEERING_L0 - (STEERING_GAIN / 2.0) * u_s
    L_right = STEERING_L0 + (STEERING_GAIN / 2.0) * u_s
    return L_left, L_right
end

"""
    depower_percentage_to_length(percentage)

Convert CSV depower percentage to tape length (m).
"""
function depower_percentage_to_length(percentage)
    u_p = percentage / 100.0  # Convert percentage to [0, 1]
    L_depower = DEPOWER_L0 + DEPOWER_GAIN * u_p
    return L_depower
end

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
    add_distance_column(data)

Add distance and cumulative_distance columns to CSV data.
Calculates 3D Euclidean distance between consecutive kite positions.
"""
function add_distance_column(data)
    n = length(data.time)
    distances = zeros(Float64, n)
    cumulative_distances = zeros(Float64, n)

    for i in 2:n
        dx = data.kite_pos_east[i] - data.kite_pos_east[i-1]
        dy = data.kite_pos_north[i] - data.kite_pos_north[i-1]
        dz = data.kite_height[i] - data.kite_height[i-1]
        distances[i] = sqrt(dx^2 + dy^2 + dz^2)
        cumulative_distances[i] = cumulative_distances[i-1] + distances[i]
    end

    # Add new fields to the named tuple
    return merge(data, (distance=distances, cumulative_distance=cumulative_distances))
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

"""
    wrap_to_pi(angle)

Wrap angle to [-π, π] range.
"""
function wrap_to_pi(angle)
    return mod(angle + π, 2π) - π
end

"""
    calc_feedforward_torque(tether_force_n, winch)

Calculate feed-forward torque from tether force.
Similar to calc_steady_torque but uses measured force.
"""
function calc_feedforward_torque(tether_force_n, winch)
    torque = -winch.drum_radius / winch.gear_ratio * tether_force_n - winch.friction
    return torque
end

"""
    calc_heading(sys_struct::SystemStructure, R_b_w)

Calculate heading angle from rotation matrix, wrapped to [-π, π].
"""
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
    heading = atan(heading_x, heading_z)
    return wrap_to_pi(heading)
end

"""
    calc_csv_heading(roll_deg, pitch_deg, yaw_deg, sys_struct)

Calculate heading from CSV Euler angles, wrapped to [-π, π].
"""
function calc_csv_heading(roll_deg, pitch_deg, yaw_deg, sys_struct)
    quat = euler_to_quaternion(roll_deg, pitch_deg, yaw_deg)
    R = SymbolicAWEModels.quaternion_to_rotation_matrix(quat)
    heading = calc_heading(sys_struct, R)
    return wrap_to_pi(heading + π)
end

function apply_force!(sys, control)
    wing = sys.wings[1]
    R_b_w = wing.R_b_w
    for point in sys.points
        distance_frac = point.pos_w ⋅ normalize(wing.pos_w) / norm(wing.pos_w)
        point.disturb .= -R_b_w[:, 1] * control * distance_frac
    end
end

function update_vel_from_csv!(sys, row, heading_pid, winch_length_pid, brake)
    @unpack wings, points, winches, segments = sys
    wing = wings[1]

    # calc delta heading
    csv_heading = calc_csv_heading(row.roll, row.pitch, row.yaw, sys)
    wing.R_b_w = calc_R_b_w(sys)
    curr_heading = calc_heading(sys, wing.R_b_w)
    delta_heading = wrap_to_pi(csv_heading - curr_heading)

    # PID control for steering based on heading error
    # Use CSV steering as feedforward
    L_left_csv, L_right_csv = steering_percentage_to_lengths(row.steering)
    u_s_csv = (L_right_csv - L_left_csv) / 2.0  # Convert back to differential (m)

    # calculate_control!(pid, r, y, uff) - r:setpoint, y:measurement, uff:feedforward
    steering_control = DiscretePIDs.calculate_control!(heading_pid, 0.0, delta_heading, u_s_csv)

    tip_push = 10
    points[2].disturb .= wing.R_b_w[:, 2] * tip_push
    points[3].disturb .= wing.R_b_w[:, 2] * tip_push
    points[20].disturb .= -wing.R_b_w[:, 2] * tip_push
    points[21].disturb .= -wing.R_b_w[:, 2] * tip_push

    # Apply steering via differential tape lengths
    segments[87].l0 = STEERING_L0 - (STEERING_GAIN/2.0)*steering_control  # Left
    segments[89].l0 = STEERING_L0 + (STEERING_GAIN/2.0)*steering_control  # Right

    # Winch length control with feed-forward torque
    winch = winches[1]

    # Calculate feed-forward torque from CSV tether force
    ff_torque = calc_feedforward_torque(row.tether_force, winch)

    # PI control for length tracking with feed-forward torque
    length_error = row.tether_len - winch.tether_len
    torque_control = DiscretePIDs.calculate_control!(
        winch_length_pid, row.tether_len, winch.tether_len, ff_torque)

    winch.brake = brake
    winch.set_value = torque_control

    # update depower (from CSV)
    L_depower = depower_percentage_to_length(row.depower)
    segments[88].l0 = L_depower

    return winch.set_value
end

"""
    update_sys_struct_from_csv!(sys_struct, csv_row_data)

Update system structure from a single CSV row.
Updates wing orientation from Euler angles and position via transform system.
"""
function update_sys_struct_from_csv!(sys, row)
    @unpack wings, points, winches, segments, transforms = sys
    wing = wings[1]
    transform = transforms[1]

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
    for (n, point_idx) in enumerate(39:44)
        points[point_idx].pos_b .= [0.0, 0.0, -n * row.tether_len / 6 * 1.01]
    end
    transform.elevation = KiteUtils.calc_elevation(csv_pos)
    transform.azimuth = KiteUtils.azimuth_east(csv_pos)
    SymbolicAWEModels.reinit!([transform], sys)

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

    # update tether length and velocity
    winches[1].tether_len = row.tether_len
    winches[1].tether_vel = row.tether_vel
    winches[1].brake = true

    # Convert CSV percentages to tape lengths
    L_left, L_right = steering_percentage_to_lengths(row.steering)
    L_depower = depower_percentage_to_length(row.depower)

    segments[87].l0 = L_left   # Left steering tape
    segments[89].l0 = L_right  # Right steering tape
    segments[88].l0 = L_depower  # Depower tape
end

"""
    run_physics_replay(csv_path; kwargs...)

Replay CSV data using physics simulation with CSV inputs.
Updates tether length and steering/depower segments from CSV data at each step.
"""
function run_physics_replay(csv_path::String;
                        start_frame=START_FRAME,
                        end_frame=END_FRAME)

    raw_data = load_flight_data(csv_path)
    csv_data = limit_frames(raw_data; start_frame, end_frame)
    csv_data = add_distance_column(csv_data)

    @info "Loading v3 kite system structure from YAML"
    set_data_path("data/v3")
    set = Settings("system.yaml")
    vsm_set_path = joinpath(get_data_path(), "vsm_settings.yaml")
    vsm_set = VortexStepMethod.VSMSettings(vsm_set_path; data_prefix=false)
    sys_struct = load_sys_struct_from_yaml("data/v3/struc_geometry.yaml";
        system_name="v3", set, wing_type=SymbolicAWEModels.REFINE, vsm_set)
    csv_sys_struct = load_sys_struct_from_yaml("data/v3/struc_geometry.yaml";
        system_name="v3", set, wing_type=SymbolicAWEModels.REFINE, vsm_set)
    sam = SymbolicAWEModel(set, sys_struct)
    init!(sam)

    max_substeps = 1000
    n_csv_steps = length(csv_data.time)
    n_steps = length(csv_data.time) * max_substeps
    @info "Creating log with $n_steps timesteps"
    sys_state = SysState(sam)
    logger = Logger(sam, n_steps)

    # Create CSV reference model for visualization
    csv_sam = SymbolicAWEModel(set, csv_sys_struct)
    init!(csv_sam)
    csv_state = SysState(csv_sam)
    csv_logger = Logger(csv_sam, n_csv_steps)

    # Loop through CSV data and update sys_struct
    @info "Replaying CSV data..."
    replay_start = time()
    sys = sam.sys_struct
    SymbolicAWEModels.set_body_frame_damping(sys, INITIAL_DAMPING)
    t = 0.0
    dt = DT_CONTROL

    # Initialize heading PID controller
    heading_pid = DiscretePID(;
        K = HEADING_KP,
        Ti = HEADING_TAU_I,
        Td = false,
        Ts = DT_CONTROL,
        umin = -0.5,
        umax = 0.5)

    # Initialize winch length PI controller
    winch_length_pid = DiscretePID(;
        K = WINCH_LENGTH_KP,
        Ti = WINCH_LENGTH_TAU_I,
        Td = false,
        Ts = DT_CONTROL,
        umin = -2000.0,
        umax = 2000.0)

    function get_row(csv_data, step)
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
            tether_force = csv_data.ground_tether_force[step] * 9.81,  # Convert kg to N
            steering = csv_data.kite_actual_steering[step],
            depower = csv_data.kite_actual_depower[step],
            distance = csv_data.distance[step],
            cumulative_distance = csv_data.cumulative_distance[step],
        )
    end

    try
        for step in 1:n_csv_steps-1
            @show step
            # Create row data structure
            csv_row = get_row(csv_data, step)
            next_row = get_row(csv_data, step+1)

            # Update CSV reference model and log
            update_sys_struct_from_csv!(csv_sam.sys_struct, csv_row)
            reinit!(csv_sam.sys_struct)
            update_sys_state!(csv_state, csv_sam)
            csv_state.time = csv_row.time
            log!(csv_logger, csv_state)

            if t <= DECAY_TIME
                current_damping = (INITIAL_DAMPING - MIN_DAMPING) * (1.0 - t / DECAY_TIME) + MIN_DAMPING
                SymbolicAWEModels.set_body_frame_damping(sam.sys_struct, current_damping)
            else
                SymbolicAWEModels.set_body_frame_damping(sam.sys_struct, MIN_DAMPING)
            end

            # Update system structure from CSV
            if step==1
                update_sys_struct_from_csv!(sam.sys_struct, csv_row)
                SymbolicAWEModels.reinit!(sam, sam.prob, FBDF())
            end

            # Distance-based stepping
            target_distance = csv_row.distance
            simulated_distance = 0.0
            last_pos_w = copy(sam.sys_struct.wings[1].pos_w)
            substep_count = 0
            min_distance_threshold = 0.01  # meters

            # if t < 0.1
                brake = true
            # else
            #     brake = false
            # end

            # Step until we match CSV distance (or hit safety limit)
            wing = sam.sys_struct.wings[1]
            last_sphere_distance = Inf
            consecutive_increases = 0
            required_increases = 3  # Need 3 consecutive increases to stop
            while (substep_count < max_substeps)
                t += dt
                set_value = update_vel_from_csv!(
                    sam.sys_struct, csv_row, heading_pid, winch_length_pid, brake)
                next_step!(sam; dt, set_values=[set_value])

                # Calculate distance moved in this substep
                next_pos = [next_row.x, next_row.y, next_row.z]
                next_elevation = KiteUtils.calc_elevation(next_pos)
                next_azimuth = KiteUtils.azimuth_east(next_pos)
                sphere_distance = norm([next_elevation, next_azimuth] -
                                       [wing.elevation, wing.azimuth])

                current_pos_w = wing.pos_w
                step_distance = norm(current_pos_w - last_pos_w)
                simulated_distance += step_distance
                last_pos_w = copy(current_pos_w)
                update_sys_state!(sys_state, sam)
                sys_state.time = t
                log!(logger, sys_state)
                substep_count += 1

                # Track consecutive increases in sphere distance
                if sphere_distance > last_sphere_distance
                    consecutive_increases += 1
                    if consecutive_increases >= required_increases
                        break
                    end
                else
                    consecutive_increases = 0
                end
                last_sphere_distance = sphere_distance
            end

            # Warn if we hit the safety limit
            if substep_count >= max_substeps
                @warn "Hit max substeps at step $step " *
                      "(simulated: $(round(simulated_distance, digits=3))m / " *
                      "target: $(round(target_distance, digits=3))m)"
            end

            # Progress reporting
            if step % max(1, div(n_steps, 10)) == 0 || step == n_steps
                elapsed = time() - replay_start
                @info "  Logged $step/$n_steps points " *
                      "(t = $(round(csv_row.time, digits=2)) s)"
            end
        end
    catch err
        if err isa AssertionError
            @warn "Still plotting"
        else
            rethrow(err)
        end
    end

    replay_elapsed = time() - replay_start
    @info "CSV replay logged in $(round(replay_elapsed, digits=2)) s"

    @info "Saving logs..."
    save_log(logger, "csv_replay")
    save_log(csv_logger, "csv_reference")
    syslog = load_log("csv_replay")
    csvlog = load_log("csv_reference")

    return sam, syslog, csv_sam, csvlog, csv_data, raw_data
end

# Main execution
sam, syslog, csv_sam, csvlog, csv_data, raw_data = run_physics_replay(CSV_PATH)
fig = plot([sam.sys_struct, csv_sam.sys_struct], [syslog, csvlog];
     plot_tether=true)

