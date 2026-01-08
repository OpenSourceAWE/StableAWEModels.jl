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
using CairoMakie
GLMakie.activate!()
using Statistics
using Rotations
using UnPack
using LinearAlgebra
using OrdinaryDiffEqBDF
using KiteUtils: calc_elevation, azimuth_east
using NonlinearSolve, ADTypes
using Dates

# Configuration parameters
SECTION = "straight_right"
STRUC_YAML_PATH = "data/v3/struc_geometry_$(SECTION).yaml"
AERO_YAML_PATH = "data/v3/aero_geometry_$(SECTION).yaml"
CSV_PATH = "data/v3/2025-10-09_16-58-33_ProtoLogger_lidar.csv"
WARN_STEP = false  # Show distance warnings

# Video frame mapping (video_frame 7182 = UTC 15:36:31.0)
VIDEO_FRAME_REF = 7182
UTC_REF_SECONDS = 15*3600 + 36*60 + 31.0  # UTC 15:36:31.0 in seconds since midnight
VIDEO_FPS = 29.97

# Extra points comparison (set to nothing to disable)
EXTRA_POINTS_CSV = "data/v3/straight_flight_reelout_frame_7182.csv"
EXTRA_POINTS_FRAME = 7182
# EXTRA_POINTS_CSV = "data/v3/right_turn_reelout_frame_7362.csv"
# EXTRA_POINTS_FRAME = 7362

# Maneuver selection - specify by UTC time
if SECTION == "straight_right"
    START_UTC = "15:36:29.0"
    END_UTC = "15:36:37.1"  # Extended to include frame 7362
elseif SECTION == "straight_left"
    START_UTC = "15:36:49.0"
    END_UTC = "15:36:52.0"
elseif SECTION == "power_depower"
    START_UTC = "15:42:11.0"
    END_UTC = "15:42:22.0"
else
    error("Unknown section: $SECTION")
end

# V3 Kite steering/depower calibration (from KCU documentation)
# Steering calibration
STEERING_L0 = 1.6  # Neutral steering tape length (m)
STEERING_GAIN = 1.4  # Maximum differential (m) at |u_s| = 1 (matches v3_kite_circular)
STEERING_MULTIPLIER = 1.0

# Depower calibration
DEPOWER_L0 = 0.2 # SUPPOSED TO BE 0.2
DEPOWER_GAIN = 5.0
DEPOWER_OFFSET = 0.0

# Restabilize: update YAML with final sys_struct positions
RESTABLE = false

INITIAL_DAMPING = [0.0, 300.0, 600.0]
DECAY_TIME = 1.0
MIN_DAMPING = [0.0, 60, 120]

# PID controller parameters for heading control
HEADING_KP = 0.0
HEADING_TAU_I = 0.0
HEADING_KD = 0.0
DT_CONTROL = 0.001

# PI controller parameters for winch length control
WINCH_LENGTH_KP = 0.01      # Low proportional gain for length tracking
WINCH_LENGTH_TAU_I = 10.0   # Integral time constant
WINCH_LENGTH_KD = 0.0      # No derivative gain

"""
    parse_time_to_seconds(time_str)

Parse "HH:MM:SS" or "HH:MM:SS.sss" to seconds since midnight.
"""
function parse_time_to_seconds(time_str::String)
    parts = split(time_str, ":")
    h = parse(Int, parts[1])
    m = parse(Int, parts[2])
    s = parse(Float64, parts[3])
    return h * 3600.0 + m * 60.0 + s
end

"""
    utc_to_video_frame(utc_seconds)

Convert UTC time (seconds since midnight) to video frame number.
"""
function utc_to_video_frame(utc_seconds::Float64)
    return round(Int, VIDEO_FRAME_REF + (utc_seconds - UTC_REF_SECONDS) * VIDEO_FPS)
end

"""
    unix_to_utc_seconds(unix_timestamp)

Convert Unix timestamp to UTC seconds since midnight.
"""
function unix_to_utc_seconds(unix_timestamp::Float64)
    dt = Dates.unix2datetime(unix_timestamp)
    return Dates.hour(dt)*3600 + Dates.minute(dt)*60 +
           Dates.second(dt) + Dates.millisecond(dt)/1000
end

"""
    load_extra_points(csv_path, sys_struct)

Load extra points from CSV and transform from camera frame to simulation frame.
CSV has columns: group, idx_in_group, x, y, z.

Alignment constraints:
1. Spanwise: CSV LE[10], LE[11] align with sim points 10, 12 (center LE)
2. To-kite direction: camera→LE matches bridle(point 27)→LE
"""
function load_extra_points(csv_path::String, sys_struct; body_offset=[0.3, 0.0, 0.2])
    df = CSV.read(csv_path, DataFrame)

    # CSV reference: LE[10], LE[11] (0-indexed in CSV, so Julia indices 11, 12)
    le_pts = [[r.x, r.y, r.z] for r in eachrow(df) if r.group == "LE"]
    csv_le10, csv_le11 = le_pts[11], le_pts[12]
    csv_le_center = (csv_le10 + csv_le11) / 2

    # CSV strut centers: strut3[1]/strut4[1] are at TE, [end] are at LE
    strut3 = [[r.x, r.y, r.z] for r in eachrow(df) if r.group == "strut3"]
    strut4 = [[r.x, r.y, r.z] for r in eachrow(df) if r.group == "strut4"]
    csv_le_center = (strut3[end] + strut4[end]) / 2
    csv_te_center = (strut3[1] + strut4[1]) / 2

    # Sim reference: points 10, 12 (center LE), point 27 (bridle)
    sim_p10 = collect(sys_struct.points[10].pos_w)
    sim_p12 = collect(sys_struct.points[12].pos_w)
    sim_le_center = (sim_p10 + sim_p12) / 2
    point_27 = collect(sys_struct.points[27].pos_w)
    cam_pos = point_27 + sys_struct.wings[1].R_b_w * [0, 0.2, 0]

    # Direction vectors
    csv_span = normalize(strut4[end] - strut3[end])

    # CSV basis: y=spanwise, z from wing center geometry, x from cross
    csv_y = csv_span
    csv_wing_center = (csv_le_center + csv_te_center) / 2
    @show csv_wing_center
    csv_z = normalize(csv_wing_center - csv_y * 0.84/2)
    csv_x = cross(csv_y, csv_z)

    # Sim basis: directly from wing rotation matrix
    R_b_w = sys_struct.wings[1].R_b_w
    sim_x = R_b_w[:, 1]
    sim_y = R_b_w[:, 2]
    sim_z = R_b_w[:, 3]

    # Rotation: R * csv_basis = sim_basis
    csv_basis = hcat(csv_x, csv_y, csv_z)
    sim_basis = hcat(sim_x, sim_y, sim_z)
    R = sim_basis * csv_basis'

    # Translation: align LE centers
    T = sim_le_center - R * csv_le_center + R_b_w * body_offset

    # Transform all points (including camera origin marker)
    all_pts = [[row.x, row.y, row.z] for row in eachrow(df)]
    push!(all_pts, zeros(3))
    transformed = [Tuple(R * p + T) for p in all_pts]

    return transformed
end

"""
    steering_percentage_to_lengths(percentage)

Convert CSV steering percentage to left/right tape lengths (m).
Percentage convention: negative = left turn, positive = right turn.
"""
function steering_percentage_to_lengths(percentage)
    u_s = percentage / 100.0  # Convert percentage to [-1, 1]
    L_left = STEERING_L0 + STEERING_GAIN * u_s
    L_right = STEERING_L0 - STEERING_GAIN * u_s
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
Returns a DataFrame with all CSV columns.
"""
function load_flight_data(csv_path::String)
    @info "Loading CSV data from: $csv_path"
    df = CSV.read(csv_path, DataFrame;
                  delim=' ',
                  silencewarnings=true,
                  normalizenames=true,
                  types=Float64,
                  strict=false)
    return df
end

"""
    find_csv_indices_by_utc(df, start_utc, end_utc)

Find CSV row indices corresponding to UTC time range.
Returns (start_idx, end_idx) tuple.
"""
function find_csv_indices_by_utc(df, start_utc::String, end_utc::String)
    start_sec = parse_time_to_seconds(start_utc)
    end_sec = parse_time_to_seconds(end_utc)

    start_idx = nothing
    end_idx = nothing

    for (i, unix_t) in enumerate(df.time)
        if ismissing(unix_t)
            continue
        end
        utc_sec = unix_to_utc_seconds(unix_t)
        if isnothing(start_idx) && utc_sec >= start_sec
            start_idx = i
        end
        if utc_sec <= end_sec
            end_idx = i
        end
    end

    if isnothing(start_idx) || isnothing(end_idx)
        error("Could not find UTC range $start_utc to $end_utc in CSV")
    end

    return start_idx, end_idx
end

"""
    limit_by_utc(df, start_utc, end_utc)

Limit DataFrame to UTC time range and convert to named tuple.
Normalizes time column to start at 0. Adds video_frame column.
"""
function limit_by_utc(df, start_utc::String, end_utc::String)
    start_idx, end_idx = find_csv_indices_by_utc(df, start_utc, end_utc)
    @info "UTC range $start_utc to $end_utc -> rows $start_idx to $end_idx"

    # Slice the dataframe
    limited_df = df[start_idx:end_idx, :]

    # Calculate video_frame from unix timestamps before normalization (Float64 for interpolation)
    video_frames = [Float64(utc_to_video_frame(unix_to_utc_seconds(t)))
                    for t in limited_df.time]

    # Normalize time to start at 0
    t0 = limited_df.time[1]
    limited_df.time .= limited_df.time .- t0

    # Convert to named tuple and add video_frame
    col_names = Tuple(Symbol(name) for name in names(limited_df))
    data = NamedTuple{col_names}(Tuple(eachcol(limited_df)))
    data = merge(data, (video_frame=video_frames,))

    # Print video frame range
    @info "Video frame range: $(video_frames[1]) to $(video_frames[end])"

    return data, start_idx
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
    interpolate_csv_data(data, n_substeps)

Linearly interpolate CSV data to create n_substeps points between each pair
of original data points. Handles missing values by propagating them.
"""
function interpolate_csv_data(data, n_substeps)
    n_original = length(data.time)
    n_interp = (n_original - 1) * n_substeps + 1

    # Create interpolated arrays for each field
    interp_data = Dict{Symbol, Vector}()

    for field in keys(data)
        # Determine element type (handle Union{Float64, Missing})
        eltype_field = eltype(data[field])
        interp_values = Vector{eltype_field}(undef, n_interp)

        for i in 1:(n_original-1)
            # Starting index in interpolated array
            start_idx = (i-1) * n_substeps + 1

            val_i = data[field][i]
            val_next = data[field][i+1]

            # Interpolate between data[i] and data[i+1]
            for j in 0:(n_substeps-1)
                idx = start_idx + j
                alpha = j / n_substeps

                # Handle missing values
                if ismissing(val_i) || ismissing(val_next)
                    interp_values[idx] = missing
                else
                    interp_values[idx] = (1 - alpha) * val_i + alpha * val_next
                end
            end
        end

        # Add the last point
        interp_values[end] = data[field][end]
        interp_data[field] = interp_values
    end

    # Convert back to named tuple
    col_names = Tuple(keys(data))
    return NamedTuple{col_names}(Tuple(interp_data[k] for k in col_names))
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
    torque = -winch.drum_radius / winch.gear_ratio * tether_force_n + winch.friction
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

# Mutable state for heading spike filter
const PREV_CSV_HEADING = Ref{Float64}(NaN)
const MAX_HEADING_CHANGE = deg2rad(20.0)  # Max allowed change per step

"""
    filter_csv_heading(heading)

Filter out spikes in CSV heading data. Rejects changes > MAX_HEADING_CHANGE.
"""
function filter_csv_heading(heading)
    if isnan(PREV_CSV_HEADING[])
        PREV_CSV_HEADING[] = heading
        return heading
    end
    delta = wrap_to_pi(heading - PREV_CSV_HEADING[])
    if abs(delta) > MAX_HEADING_CHANGE
        # Spike detected, keep previous value
        return PREV_CSV_HEADING[]
    end
    PREV_CSV_HEADING[] = heading
    return heading
end

function apply_force!(sys, control)
    wing = sys.wings[1]
    R_b_w = wing.R_b_w
    for point in sys.points
        distance_frac = point.pos_w ⋅ normalize(wing.pos_w) / norm(wing.pos_w)
        point.disturb .= -R_b_w[:, 1] * control * distance_frac
    end
end

# Mutable state for simulation cumulative distance
const SIM_PREV_POS = Ref{Vector{Float64}}(zeros(3))
const SIM_CUMULATIVE_DIST = Ref{Float64}(0.0)

function reset_distance_tracker!()
    SIM_PREV_POS[] = zeros(3)
    SIM_CUMULATIVE_DIST[] = 0.0
end

function update_sim_distance!(wing_pos)
    if SIM_PREV_POS[] == zeros(3)
        SIM_PREV_POS[] = copy(wing_pos)
        return 0.0
    end
    dist = norm(wing_pos - SIM_PREV_POS[])
    SIM_CUMULATIVE_DIST[] += dist
    SIM_PREV_POS[] = copy(wing_pos)
    return SIM_CUMULATIVE_DIST[]
end

function update_vel_from_csv!(sys, row, brake, heading_pid)
    @unpack wings, points, winches, segments = sys
    wing = wings[1]

    # Calc delta heading (with spike filter on CSV heading)
    raw_csv_heading = calc_csv_heading(row.roll, row.pitch, row.yaw, sys)
    csv_heading = filter_csv_heading(raw_csv_heading)
    wing.R_b_w = calc_R_b_w(sys)
    curr_heading = calc_heading(sys, wing.R_b_w)
    delta_heading = -wrap_to_pi(csv_heading - curr_heading)

    sim_cumulative_dist = update_sim_distance!(wing.pos_w)
    sys.set.v_wind = row.wind_at_kite

    # Apply steering via differential tape lengths
    # PID control for steering based on heading error
    steering_control = DiscretePIDs.calculate_control!(
        heading_pid, 0.0, delta_heading, 0.0)
    steering = clamp(row.steering, -100.0, 100.0)
    L_left, L_right = steering_percentage_to_lengths(steering * STEERING_MULTIPLIER)
    segments[87].l0 = L_left + STEERING_GAIN * steering_control
    segments[89].l0 = L_right - STEERING_GAIN * steering_control

    # Winch length control with feed-forward torque
    winch = winches[1]

    # Calculate feed-forward torque from CSV tether force
    ff_torque = calc_feedforward_torque(row.tether_force, winch)

    # # PI control for length tracking with feed-forward torque
    # length_error = row.tether_len - winch.tether_len
    # torque_control = DiscretePIDs.calculate_control!(
    #     winch_length_pid, row.tether_len, winch.tether_len, ff_torque)

    winch.brake = brake
    winch.set_value = ff_torque

    # Update depower from CSV
    L_depower = depower_percentage_to_length(row.depower + DEPOWER_OFFSET)
    segments[88].l0 = L_depower

    # Return torque and effective tape percentages (with multiplier/offset + control)
    effective_steering = steering * STEERING_MULTIPLIER + steering_control * 100
    effective_depower = row.depower + DEPOWER_OFFSET
    return winch.set_value, effective_steering, effective_depower
end

"""
    update_sys_struct_from_csv!(sys_struct, csv_row_data)

Update system structure from a single CSV row.
Updates wing orientation from Euler angles and position via transform system.
Uses nonlinear solve to find rotation angle that achieves desired heading.
"""
function update_sys_struct_from_csv!(sys, row)
    @unpack wings, points, winches, segments, transforms = sys
    wing = wings[1]
    transform = transforms[1]

    # calc target heading from CSV
    quat = euler_to_quaternion(row.roll, row.pitch, row.yaw)
    csv_heading = calc_heading(sys,
        SymbolicAWEModels.quaternion_to_rotation_matrix(quat)) + pi
    wing.R_b_w = calc_R_b_w(sys)
    curr_heading = calc_heading(sys, wing.R_b_w)

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
    transform.heading = csv_heading
    SymbolicAWEModels.reinit!([transform], sys)

    # apply vel
    csv_vel = [row.vx, row.vy, row.vz]
    wing.vel_w .= csv_vel
    for point in points
        transform_frac = point.pos_w ⋅ normalize(wing.pos_w) / norm(wing.pos_w)
        point.vel_w .= transform_frac * csv_vel
    end

    # update tether length and velocity
    # winches[1].tether_len = row.tether_len
    # winches[1].tether_vel = row.tether_vel
    winches[1].brake = true

    # Convert CSV percentages to tape lengths
    L_left, L_right = steering_percentage_to_lengths(row.steering)
    L_depower = depower_percentage_to_length(row.depower + DEPOWER_OFFSET)

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
                        start_utc=START_UTC,
                        end_utc=END_UTC,
                        n_substeps=5)

    df = load_flight_data(csv_path)
    limited_data, start_idx = limit_by_utc(df, start_utc, end_utc)
    limited_data = add_distance_column(limited_data)

    # Interpolate CSV data for smoother control
    @info "Interpolating CSV data with $n_substeps substeps"
    csv_data = interpolate_csv_data(limited_data, n_substeps)

    @info "Loading v3 kite system structure from YAML"
    set_data_path("data/v3")
    set = Settings("system.yaml")
    set.g_earth = 9.81
    set.l_tether = 212.68
    vsm_set_path = joinpath(get_data_path(), "vsm_settings_reduced_for_coupling.yaml")
    vsm_set = VortexStepMethod.VSMSettings(vsm_set_path; data_prefix=false)
    depower_int = Int(round(40+DEPOWER_OFFSET))
    @info "Using depower int $(depower_int)"
    vsm_set.wings[1].geometry_file = AERO_YAML_PATH
    sys_struct = load_sys_struct_from_yaml(STRUC_YAML_PATH;
        system_name="v3", set, wing_type=SymbolicAWEModels.REFINE, vsm_set)
    csv_sys_struct = load_sys_struct_from_yaml(STRUC_YAML_PATH;
        system_name="v3", set, wing_type=SymbolicAWEModels.REFINE, vsm_set)
    sam = SymbolicAWEModel(set, sys_struct)
    init!(sam)

    n_steps = length(csv_data.time)
    @info "Creating log with $n_steps timesteps"
    sys_state = SysState(sam)
    logger = Logger(sam, n_steps)

    # Create CSV reference model for visualization
    csv_sam = SymbolicAWEModel(set, csv_sys_struct)
    init!(csv_sam)
    csv_state = SysState(csv_sam)
    csv_logger = Logger(csv_sam, n_steps)

    # Storage for tape percentages (for plotting)
    csv_tape_times = Float64[]
    csv_tape_steering_pct = Float64[]
    csv_tape_depower_pct = Float64[]
    phys_tape_times = Float64[]
    phys_tape_steering_pct = Float64[]
    phys_tape_depower_pct = Float64[]

    # Loop through CSV data and update sys_struct
    @info "Replaying CSV data..."
    replay_start = time()
    sys = sam.sys_struct
    SymbolicAWEModels.set_body_frame_damping(sys, INITIAL_DAMPING)

    # Calculate dt from interpolated CSV timesteps
    dt = csv_data.time[2] - csv_data.time[1]
    @info "Using timestep dt = $dt s"

    # Calculate initial delta between set.l_tether and CSV tether length
    first_csv_tether_len = csv_data.ground_tether_length[1]
    tether_len_delta = set.l_tether - first_csv_tether_len
    @info "Tether length delta" set_l_tether=set.l_tether csv_tether=first_csv_tether_len delta=tether_len_delta

    # Reset heading spike filter and distance tracker
    PREV_CSV_HEADING[] = NaN
    reset_distance_tracker!()

    # Initialize heading PID controller
    heading_pid = DiscretePID(;
        K = HEADING_KP,
        Ti = HEADING_TAU_I,
        Td = false,
        Ts = dt,
        umin = -1.0,
        umax = 1.0)

    function get_row(csv_data, step)
        csv_row = (
            time = csv_data.time[step],
            video_frame = round(Int, csv_data.video_frame[step]),
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
            wind_at_kite = coalesce(csv_data.lidar_wind_velocity_at_kite_mps[step], 10.0),
            angle_of_attack = deg2rad(coalesce(csv_data.airspeed_angle_of_attack[step], NaN))
        )
    end

    try
        for step in 1:n_steps-1
            # Create row data structure
            csv_row = get_row(csv_data, step)

            # Update CSV reference model and log
            update_sys_struct_from_csv!(csv_sam.sys_struct, csv_row)
            SymbolicAWEModels.reinit!(csv_sam, csv_sam.prob, FBDF())
            update_sys_state!(csv_state, csv_sam)
            csv_state.winch_force[1] = csv_row.tether_force
            csv_state.AoA = csv_row.angle_of_attack
            csv_state.time = csv_row.time
            csv_state.l_tether[1] = csv_row.tether_len
            csv_state.v_reelout[1] = csv_row.tether_vel
            csv_state.v_wind_gnd[1] = csv_row.wind_at_kite
            log!(csv_logger, csv_state)

            # Store CSV tape percentages for plotting
            push!(csv_tape_times, csv_row.time)
            push!(csv_tape_steering_pct, csv_row.steering)
            push!(csv_tape_depower_pct, csv_row.depower)

            # Update system structure from CSV on first step
            if step == 1
                update_sys_struct_from_csv!(sam.sys_struct, csv_row)
                SymbolicAWEModels.reinit!(sam, sam.prob, FBDF())
            end

            # Update damping
            t = csv_row.time
            if t <= DECAY_TIME
                current_damping = (INITIAL_DAMPING - MIN_DAMPING) *
                                  (1.0 - t / DECAY_TIME) + MIN_DAMPING
                SymbolicAWEModels.set_body_frame_damping(sam.sys_struct, current_damping)
            end

            # Apply control and step
            brake = true
            set_value, eff_steering, eff_depower = update_vel_from_csv!(
                sam.sys_struct, csv_row, brake, heading_pid)

            # Update winch tether length and velocity from CSV
            sam.sys_struct.winches[1].tether_len = csv_row.tether_len + tether_len_delta
            sam.sys_struct.winches[1].tether_vel = csv_row.tether_vel
            SymbolicAWEModels.reinit!(sam, sam.prob, FBDF())

            # Distance warnings
            if WARN_STEP
                csv_dist = csv_row.cumulative_distance
                sim_dist = SIM_CUMULATIVE_DIST[]
                dist_error = csv_dist - sim_dist
                if dist_error < -dt * 10
                    println("Sim ahead by $(round(-dist_error, digits=2))m")
                elseif dist_error > dt * 10
                    println("Sim behind by $(round(dist_error, digits=2))m")
                end
            end

            next_step!(sam; dt=dt, set_values=[set_value])

            # Plot comparison at specified frame
            if !isnothing(EXTRA_POINTS_CSV) &&
               Int(round(csv_row.video_frame)) == EXTRA_POINTS_FRAME
                @info "Plotting comparison at frame $(EXTRA_POINTS_FRAME)..."
                extra_pts = load_extra_points(EXTRA_POINTS_CSV, sam.sys_struct)
                comparison_scene = plot(sam.sys_struct; extra_points=extra_pts)
                scr = display(comparison_scene)
                @info "Close the plot window to continue..."
                wait(scr)
            end

            # Log state
            update_sys_state!(sys_state, sam)
            sys_state.time = t
            log!(logger, sys_state)

            # Store physics tape percentages for plotting (with multiplier/offset + control)
            push!(phys_tape_times, t)
            push!(phys_tape_steering_pct, eff_steering)
            push!(phys_tape_depower_pct, eff_depower)

            # Progress reporting
            if step % max(1, div(n_steps, 10)) == 0 || step == n_steps
                elapsed = time() - replay_start
                @info "  Step $step/$n_steps " *
                      "(t = $(round(csv_row.time, digits=2)) s, " *
                      "video_frame = $(csv_row.video_frame))"
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

    # Restabilize: update YAML with final sys_struct positions if enabled
    if RESTABLE
        @info "Restabilizing: updating YAML files with final positions..."
        SymbolicAWEModels.update_yaml_from_sys_struct!(
            sam.sys_struct, STRUC_YAML_PATH, STRUC_YAML_PATH,
            AERO_YAML_PATH, AERO_YAML_PATH)
        @info "YAML files updated: $STRUC_YAML_PATH, $AERO_YAML_PATH"
    end

    # Create tape percentages data for plotting
    phys_tape_pct = (
        time = phys_tape_times,
        steering = phys_tape_steering_pct,
        depower = phys_tape_depower_pct
    )
    csv_tape_pct = (
        time = csv_tape_times,
        steering = csv_tape_steering_pct,
        depower = csv_tape_depower_pct
    )

    return sam, syslog, csv_sam, csvlog, csv_data, phys_tape_pct, csv_tape_pct
end

# Main execution
sam, syslog, csv_sam, csvlog, csv_data, phys_tape_pct, csv_tape_pct = run_physics_replay(CSV_PATH)

# Display with GLMakie
fig = plot([sam.sys_struct, csv_sam.sys_struct], [syslog, csvlog];
     plot_tether=true, plot_aero_force=false, plot_kite_vel=true,
     plot_wind=true, plot_reelout=false, plot_v_app=false, plot_turn_rates=true,
     tape_lengths=[phys_tape_pct, csv_tape_pct],
     suffixes=["phys", "csv"])
display(fig)

# Save PDF with CairoMakie
CairoMakie.activate!()
fig_pdf = plot([sam.sys_struct, csv_sam.sys_struct], [syslog, csvlog];
     plot_tether=true, plot_aero_force=false, plot_kite_vel=true,
     plot_wind=true, plot_reelout=false, plot_v_app=false, plot_turn_rates=true,
     plot_gk=true,
     tape_lengths=[phys_tape_pct, csv_tape_pct],
     suffixes=["phys", "csv"])
CairoMakie.save("csv_replay_$(SECTION).pdf", fig_pdf)
@info "Saved plot to csv_replay_$(SECTION).pdf"
GLMakie.activate!()

sphere = plot_sphere_trajectory([syslog,csvlog])

