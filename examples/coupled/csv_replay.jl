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
END_FRAME = START_FRAME + 30 # Last frame to replay (use nothing for all frames)
REMAKE_CACHE = false

# V3 Kite steering/depower calibration (from KCU documentation)
STEERING_L0 = 1.506  # Neutral steering tape length (m)
DEPOWER_L0 = 3.129

INITIAL_DAMPING = 100.0
DECAY_TIME = 2.0
MIN_DAMPING = 50.0

# PID controller parameters for heading control
HEADING_KP = 0.05    # Low proportional gain
HEADING_TAU_I = 0.1  # Integral time constant (seconds)
HEADING_KD = 0.0     # No derivative gain
DT_CONTROL = 0.01    # Control timestep

# PID controller parameters for winch control
WINCH_KP = 0.1    # Proportional gain for tether length
WINCH_TAU_I = 0.001    # Integral time constant (seconds)
WINCH_KD = 0.0       # No derivative gain
INIT_I = 300

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

function update_vel_from_csv!(sys, row, heading_pid, winch_pid)
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
    delta_heading = csv_heading - curr_heading - 2pi

    # PID control for steering based on heading error
    # calculate_control!(pid, r, y, uff) - r:setpoint, y:measurement, uff:feedforward
    steering_control = DiscretePIDs.calculate_control!(heading_pid, 0.0, delta_heading, 0.0)
    total_steering_len = 2 * STEERING_L0
    # Apply steering: positive control increases seg 87, decreases seg 89
    # segments[87].l0 = STEERING_L0 + steering_control
    # segments[89].l0 = total_steering_len - segments[87].l0
    # @show steering_control

    # PID control for winch to maintain tether length
    # Negate output: tether too long (positive error) needs negative torque (reel in)
    winch = winches[1]
    torque_control = -DiscretePIDs.calculate_control!(
        winch_pid, winch.tether_len, row.tether_len, 0.0)
    winch.brake = false
    winch.set_value = torque_control

    # update depower (from CSV)
    # segments[88].l0 = DEPOWER_L0 + row.depower / 100
    return winch.set_value
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
    return INIT_I
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

    @info "Loading v3 kite system structure from YAML"
    set_data_path("data/v3")
    set = Settings("system.yaml")
    vsm_set_path = joinpath(get_data_path(), "vsm_settings.yaml")
    vsm_set = VortexStepMethod.VSMSettings(vsm_set_path; data_prefix=false)
    sys_struct = load_sys_struct_from_yaml("data/v3/struc_geometry.yaml";
        system_name="v3", set, wing_type=REFINE, vsm_set)
    sam = SymbolicAWEModel(set, sys_struct)
    init!(sam)

    n_steps = length(csv_data.time)
    @info "Creating log with $n_steps timesteps"
    sys_state = SysState(sam)
    logger = Logger(sam, n_steps)

    # Loop through CSV data and update sys_struct
    @info "Replaying CSV data..."
    replay_start = time()
    sys = sam.sys_struct
    SymbolicAWEModels.set_world_frame_damping(sys, INITIAL_DAMPING)
    t = 0.0
    dt = 0.01

    # Initialize heading PID controller
    heading_pid = DiscretePID(;
        K = HEADING_KP,
        Ti = HEADING_TAU_I,
        Td = false,
        Ts = DT_CONTROL,
        umin = -0.5,
        umax = 0.5)

    # Initialize winch PID controller with integrator at current force
    # Negative because longer tether needs more negative torque to reel in
    initial_winch_force = norm(sys.winches[1].force)
    winch_pid = DiscretePID(;
        K = WINCH_KP,
        Ti = WINCH_TAU_I,
        Td = false,
        Ts = DT_CONTROL,
        I = INIT_I,
        umin = -50000.0,
        umax = 50000.0)
    for step in 1:n_steps
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
        if t <= DECAY_TIME
            current_damping = (INITIAL_DAMPING - MIN_DAMPING) * (1.0 - t / DECAY_TIME) + MIN_DAMPING
            SymbolicAWEModels.set_world_frame_damping(sam.sys_struct, current_damping)
        else
            SymbolicAWEModels.set_world_frame_damping(sam.sys_struct, MIN_DAMPING)
        end

        # Update system structure from CSV
        if step==1
            set_value = update_sys_struct_from_csv!(sam.sys_struct, csv_row)
            SymbolicAWEModels.reinit!(sam, sam.prob, FBDF())
        end
        if step==1
            substeps = Int(round((csv_data.time[step+1] - csv_data.time[step]) / dt))
        else
            substeps = Int(round((csv_data.time[step] - csv_data.time[step-1]) / dt))
        end
        for substep in 1:substeps
            t += dt
            set_value = update_vel_from_csv!(sam.sys_struct, csv_row, heading_pid, winch_pid)
            next_step!(sam; dt, set_values=[set_value])
        end

        # Update sys_state from sam and log it
        update_sys_state!(sys_state, sam)
        sys_state.time = t
        log!(logger, sys_state)

        # Progress reporting
        if step % max(1, div(n_steps, 10)) == 0 || step == n_steps
            elapsed = time() - replay_start
            @info "  Logged $step/$n_steps points " *
                  "(t = $(round(csv_row.time, digits=2)) s)"
        end
    end

    replay_elapsed = time() - replay_start
    @info "CSV replay logged in $(round(replay_elapsed, digits=2)) s"

    @info "Saving replay log..."
    save_log(logger, "csv_replay")
    syslog = load_log("csv_replay")

    # Calculate CSV heading for comparison
    csv_heading = zeros(n_steps)
    csv_tether_len = zeros(n_steps)
    for step in 1:n_steps
        quat = euler_to_quaternion(
            csv_data.kite_0_roll[step],
            csv_data.kite_0_pitch[step],
            csv_data.kite_0_yaw[step]
        )
        csv_pos = [csv_data.kite_pos_east[step],
                   csv_data.kite_pos_north[step],
                   csv_data.kite_height[step]]
        R = SymbolicAWEModels.quaternion_to_rotation_matrix(quat)
        csv_heading[step] = calc_heading(sys, R) - pi
        csv_tether_len[step] = csv_data.ground_tether_length[step]
    end

    return sam, syslog, csv_data, csv_heading, csv_tether_len
end

# Main execution
sam_phys, syslog_phys, csv_data_phys, csv_heading, csv_tether_len = run_physics_replay(CSV_PATH)
plot(sam_phys.sys_struct, syslog_phys;
     heading_setpoint=csv_heading,
     tether_len_setpoint=csv_tether_len,
     plot_tether=true)

