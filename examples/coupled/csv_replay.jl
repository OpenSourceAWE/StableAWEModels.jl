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

include("utils.jl")

# Configuration parameters
SECTION = "straight_right"
CSV_PATH = "data/v3/2025-10-09_16-58-33_ProtoLogger_lidar.csv"
WARN_STEP = false  # Show distance warnings

# Geometry configuration (must match settle_refine_wing.jl output)
REDUCE_TIP_LE = true         # Reduce tip LE segments (47,48,57,58)
REDUCE_TE = true             # Reduce TE segments
DEPOWER_PERCENTAGE = 40      # depower [0, 100]
TETHER_LENGTH = 240          # Total tether length (m)
TE_FRAC = 0.95               # Factor for TE wires (segments 20-28)

# Build geometry filename suffix
TETHER_INT = Int(round(TETHER_LENGTH))
TIP_LE_STR = REDUCE_TIP_LE ? "tipLE" : "no_tipLE"
TE_STR = REDUCE_TE ? "TE$(Int(round(TE_FRAC*100)))" : "no_TE"
GEOM_SUFFIX = "depower$(DEPOWER_PERCENTAGE)_tether$(TETHER_INT)_$(TIP_LE_STR)_$(TE_STR)"

STRUC_YAML_PATH = "data/v3/struc_geometry_$(GEOM_SUFFIX).yaml"
AERO_YAML_PATH = "data/v3/aero_geometry_$(GEOM_SUFFIX).yaml"

# Video frame mapping (video_frame 7182 = UTC 15:36:31.0)
VIDEO_FRAME_REF = 7182
UTC_REF_SECONDS = 15*3600 + 36*60 + 31.0  # UTC 15:36:31.0 in seconds since midnight
VIDEO_FPS = 29.97

# Extra points comparison (set to nothing to disable)

# Maneuver selection - specify by UTC time
if SECTION == "straight_right"
    START_UTC = "15:36:29.0"
    END_UTC = "15:36:37.1"  # Extended to include frame 7362
    EXTRA_POINTS_CSV = "data/v3/straight_flight_reelout_frame_7182.csv"
    EXTRA_POINTS_FRAME = 7182
    # EXTRA_POINTS_CSV = "data/v3/right_turn_reelout_frame_7362.csv"
    # EXTRA_POINTS_FRAME = 7362
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
STEERING_L0 = 1.6  # Neutral steering tape length (m) TODO: was 1.6
STEERING_GAIN = 1.4  # Maximum differential (m) at |u_s| = 1
DEPOWER_L0 = 0.2 # TODO: was 0.2
DEPOWER_GAIN = 5.0

STEERING_MULTIPLIER = 1.0
DEPOWER_OFFSET = 0.0

# Restabilize: update YAML with final sys_struct positions
RESTABLE = false

# Stop immediately after plotting at EXTRA_POINTS_FRAME
STOP_EARLY = false

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
    calc_feedforward_torque(tether_force_n, winch)

Calculate feed-forward torque from tether force.
Similar to calc_steady_torque but uses measured force.
"""
function calc_feedforward_torque(tether_force_n, winch)
    torque = -winch.drum_radius / winch.gear_ratio * tether_force_n + winch.friction
    return torque
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
    set.l_tether = TETHER_LENGTH
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

    # Storage for extra points (loaded at EXTRA_POINTS_FRAME)
    extra_pts = nothing
    extra_groups = nothing

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

            # Plot and save comparison at specified frame
            if !isnothing(EXTRA_POINTS_CSV) &&
               Int(round(csv_row.video_frame)) == EXTRA_POINTS_FRAME
                @info "Plotting comparison at frame $(EXTRA_POINTS_FRAME)..."
                extra_pts, extra_groups = load_extra_points(EXTRA_POINTS_CSV, sam.sys_struct)

                # Save PDFs for all three views
                CairoMakie.activate!()
                for dir in (:front, :side, :top)
                    scene = plot_body_frame_local(sam.sys_struct;
                        extra_points=extra_pts, extra_groups=extra_groups, dir)
                    pdf_filename = "data/v3/body_frame_$(dir)_frame$(EXTRA_POINTS_FRAME)_$(GEOM_SUFFIX).pdf"
                    save(pdf_filename, scene)
                    @info "Plot saved" pdf_filename
                end

                # Display interactive plot
                GLMakie.activate!()
                comparison_scene = plot_body_frame_local(sam.sys_struct;
                    extra_points=extra_pts, extra_groups=extra_groups, dir=:side)
                scr = display(comparison_scene)
                @info "Close the plot window to continue..."
                wait(scr)
                aoa_scene = plot_aoa(sam.sys_struct)
                scr = display(aoa_scene)
                @info "Close the plot window to continue..."
                wait(scr)

                if STOP_EARLY
                    @info "STOP_EARLY enabled, breaking out of loop"
                    break
                end
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

    return sam, syslog, csv_sam, csvlog, csv_data, phys_tape_pct, csv_tape_pct,
           extra_pts, extra_groups
end

# Main execution
sam, syslog, csv_sam, csvlog, csv_data, phys_tape_pct, csv_tape_pct,
    extra_pts, extra_groups = run_physics_replay(CSV_PATH)

# Display with GLMakie
fig = plot([sam.sys_struct, csv_sam.sys_struct], [syslog, csvlog];
     plot_tether=true, plot_aero_force=false, plot_kite_vel=true,
     plot_wind=false, plot_reelout=false, plot_v_app=true, plot_turn_rates=true,
     tape_lengths=[phys_tape_pct, csv_tape_pct],
     suffixes=["phys", "csv"])
# display(fig)

# # Save PDF with CairoMakie
# CairoMakie.activate!()
# fig_pdf = plot([sam.sys_struct, csv_sam.sys_struct], [syslog, csvlog];
#      plot_tether=true, plot_aero_force=false, plot_kite_vel=true,
#      plot_wind=true, plot_reelout=false, plot_v_app=false, plot_turn_rates=true,
#      plot_gk=true,
#      tape_lengths=[phys_tape_pct, csv_tape_pct],
#      suffixes=["phys", "csv"])
# CairoMakie.save("csv_replay_$(SECTION).pdf", fig_pdf)
# @info "Saved plot to csv_replay_$(SECTION).pdf"
# GLMakie.activate!()

sphere = plot_sphere_trajectory([syslog,csvlog])

