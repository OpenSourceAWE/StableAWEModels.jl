using SymbolicAWEModels
using VortexStepMethod
using LinearAlgebra
using Statistics
using Dates
using StaticArrays
using KiteUtils

const WINDOW_SEC = 200.0

"""
	parse_tags(log_name::AbstractString)

Parse tags from zenith-circle batch log names.
Expected pattern fragment:
  zenith_circle__up_<up>_us_<us>_vw_<vw>_lt_<lt>_el_<elev|yaml>_g_<g10|yaml>_date_...
Returns (up::Float64, us::Float64, vw::Int, lt::Int, set_elev::Union{Float64,Nothing}, g_earth::Union{Float64,Nothing}).
"""
function parse_tags(log_name::AbstractString)
	m = match(r"zenith_circle__up_([0-9]+)_us_([0-9._-]+)_vw_([0-9]+)_lt_([0-9]+)_el_([0-9]+|yaml)_g_([0-9]+|yaml)", log_name)
	m === nothing && return nothing
	up_raw = parse(Float64, m.captures[1])
	us_tokens = split(m.captures[2], "_")
	us_raw = parse(Float64, us_tokens[1])
	v_wind = parse(Int, m.captures[3])
	lt = parse(Int, m.captures[4])
	elev_tok = m.captures[5]
	g_tok = m.captures[6]
	set_elev = elev_tok == "yaml" ? nothing : parse(Float64, elev_tok)
	g_earth = begin
		if g_tok == "yaml"
			nothing
		else
			# stored as g*10 integer
			parse(Float64, g_tok) / 10
		end
	end
	return up_raw / 100, us_raw / 100, v_wind, lt, set_elev, g_earth
end

function adjust_tether_length!(sys::SystemStructure, set::Settings, tether_length_raw; tether_point_idxs=39:44)
	tether_length = float(tether_length_raw)

	if !isempty(set.l_tethers)
		set.l_tethers[1] = tether_length
	end

	n_points = length(tether_point_idxs)
	for (n, p_idx) in enumerate(tether_point_idxs)
		pos = (0.0, 0.0, -n * tether_length / n_points)
		sys.points[p_idx].pos_cad .= pos
		sys.points[p_idx].pos_b .= pos
	end

	if !isempty(sys.transforms)
		transform = sys.transforms[1]
		if !isempty(sys.wings) && norm(sys.wings[1].pos_w) > 0
			target_pos = normalize(sys.wings[1].pos_w) * tether_length
			transform.elevation = KiteUtils.calc_elevation(target_pos)
			transform.azimuth = KiteUtils.azimuth_east(target_pos)
		end
		SymbolicAWEModels.reinit!([transform], sys)
	end

	if !isempty(sys.winches)
		winch = sys.winches[1]
		winch.tether_len = tether_length
		winch.tether_vel = 0.0
		winch.brake = true
	end
	return nothing
end

function build_sys(; v_wind=10.0, tether_length=150.0, g_earth=nothing)
	wing_type = SymbolicAWEModels.REFINE
	set_data_path("data/v3")
	set = Settings("system.yaml")
	set.v_wind = v_wind
	set.upwind_dir = -90.0
	if g_earth !== nothing
		set.g_earth = g_earth
	end
	if !isempty(set.l_tethers)
		set.l_tethers[1] = tether_length
	end

	model_name = "v3_refine"
	struc_yaml_path = joinpath("data", "v3", "CORRECT_struc_geometry.yaml")
	vsm_set_path = joinpath(get_data_path(), "CORRECT_vsm_settings.yaml")
	vsm_set = VortexStepMethod.VSMSettings(vsm_set_path; data_prefix=false)
	vsm_set.wings[1].n_panels = 36

	sys = load_sys_struct_from_yaml(struc_yaml_path;
		system_name=model_name, set, wing_type, vsm_set)
	adjust_tether_length!(sys, set, tether_length)
	return sys
end

function mean_last_window(values::AbstractVector{<:Real}, times::AbstractVector{<:Real};
						  window_sec::Real=WINDOW_SEC)
	@assert length(values) == length(times)
	t_end = times[end]
	mask = times .>= (t_end - window_sec)
	if !any(mask)
		mask = trues(length(times))
	end
	data = values[mask]
	data = data[isfinite.(data)]
	return isempty(data) ? NaN : mean(data)
end

function analyze_log(lg, sys; window_sec::Real=WINDOW_SEC)
	sl = lg.syslog
	if length(sl.time) < 2
		return (
			aero_force=NaN, v_app=NaN, kite_vel=NaN, aoa=NaN, final_elevation=NaN, azimuth=NaN,
			proj_area=NaN, ref_area=NaN
		)
	end

	aero_force_z = [sl.aero_force_b[i][3] for i in eachindex(sl.aero_force_b)]
	aero_force = mean_last_window(aero_force_z, sl.time; window_sec)

	v_app = mean_last_window(sl.v_app, sl.time; window_sec)

	v_kite_norm = [norm(v) for v in sl.vel_kite]
	kite_vel = mean_last_window(v_kite_norm, sl.time; window_sec)

	aoa_deg = rad2deg.(sl.AoA)
	aoa = mean_last_window(aoa_deg, sl.time; window_sec)

	elevation_deg = rad2deg.(sl.elevation)
	final_elevation = mean_last_window(elevation_deg, sl.time; window_sec)

	azimuth_deg = rad2deg.(sl.azimuth)
	azimuth = mean_last_window(azimuth_deg, sl.time; window_sec)

	# Compute projected area from end-state (last sample)
	last_idx = length(sl.time)
	proj_area = compute_projected_area(sl, last_idx, sys)
	
	# Reference area from system
	ref_area = calc_ref_area(sys)

	return (
		aero_force=aero_force,
		v_app=v_app,
		kite_vel=kite_vel,
		aoa=aoa,
		final_elevation=final_elevation,
		azimuth=azimuth,
		proj_area=proj_area,
		ref_area=ref_area,
	)
end

function calc_ref_area(sys::SystemStructure)
	isempty(sys.wings) && return NaN
	wing = sys.wings[1]
	hasproperty(wing, :vsm_aero) || return NaN
	panels = wing.vsm_aero.panels
	isempty(panels) && return NaN
	return sum(p.chord * p.width for p in panels)
end

function mid_te_position(sl, k; eps=1e-12)
	"""Compute mid trailing-edge position (average of points 10 & 11)"""
	Xk = sl.X[k]; Yk = sl.Y[k]; Zk = sl.Z[k]
	if length(Xk) < 11 || length(Yk) < 11 || length(Zk) < 11
		return nothing
	end
	pte10 = SVector{3, Float64}(Xk[10], Yk[10], Zk[10])
	pte11 = SVector{3, Float64}(Xk[11], Yk[11], Zk[11])
	return (pte10 + pte11) / 2
end

function mid_le_position(sl, k; eps=1e-12)
	"""Compute mid leading-edge position (average of points 12 & 14)"""
	Xk = sl.X[k]; Yk = sl.Y[k]; Zk = sl.Z[k]
	if length(Xk) < 14 || length(Yk) < 14 || length(Zk) < 14
		return nothing
	end
	ple12 = SVector{3, Float64}(Xk[12], Yk[12], Zk[12])
	ple14 = SVector{3, Float64}(Xk[14], Yk[14], Zk[14])
	return (ple12 + ple14) / 2
end

function compute_projected_area(sl, k, sys::SystemStructure; eps=1e-12)
	"""Compute projected area of wing onto plane defined by body y-axis and chord direction.
	
	Projects wing quadrilaterals (nodes 2-7) onto plane spanned by:
	- Chord direction (mid LE to mid TE)
	- Body y-axis from wing orientation
	"""
	Xk = sl.X[k]; Yk = sl.Y[k]; Zk = sl.Z[k]
	if length(Xk) < 7 || length(Yk) < 7 || length(Zk) < 7
		return NaN
	end
	
	# Get mid LE and mid TE
	p_le = mid_le_position(sl, k; eps=eps)
	p_te = mid_te_position(sl, k; eps=eps)
	
	if p_le === nothing || p_te === nothing
		return NaN
	end
	
	# Chord direction (from LE to TE)
	chord_dir = p_te - p_le
	chord_norm = norm(chord_dir)
	if chord_norm <= eps
		return NaN
	end
	chord_unit = chord_dir / chord_norm
	
	# Body y-axis from wing orientation
	if length(sl.orient) < k
		return NaN
	end
	R_b_w = SymbolicAWEModels.quaternion_to_rotation_matrix(sl.orient[k])
	body_y = SVector{3, Float64}(R_b_w[1, 2], R_b_w[2, 2], R_b_w[3, 2])
	
	# Normalize body y-axis
	body_y_norm = norm(body_y)
	if body_y_norm <= eps
		return NaN
	end
	body_y_unit = body_y / body_y_norm
	
	# Ensure body y is orthogonal to chord (Gram-Schmidt)
	body_y_unit = body_y_unit - dot(body_y_unit, chord_unit) * chord_unit
	body_y_norm2 = norm(body_y_unit)
	if body_y_norm2 <= eps
		return NaN
	end
	body_y_unit = body_y_unit / body_y_norm2
	
	# Get wing node positions (2-7)
	wing_nodes = [SVector{3, Float64}(Xk[i], Yk[i], Zk[i]) for i in 2:7]
	
	# Project nodes onto plane: express in (chord, body_y) basis
	projected = [SVector{2, Float64}(dot(node - p_le, chord_unit), dot(node - p_le, body_y_unit)) for node in wing_nodes]
	
	# Compute area using shoelace formula for polygon
	area = 0.0
	for i in 1:length(projected)-1
		area += projected[i][1] * projected[i+1][2] - projected[i+1][1] * projected[i][2]
	end
	area = abs(area) / 2.0
	
	return area
end

function find_log_names(batch_dir::AbstractString)
	isdir(batch_dir) || error("Batch folder not found: $batch_dir")
	files = readdir(batch_dir; join=true)
	names = String[]
	for file in files
		isfile(file) || continue
		endswith(file, ".txt") && continue
		name = splitext(basename(file))[1]
		parse_tags(name) === nothing && continue
		push!(names, name)
	end
	return sort(unique(names))
end

function write_csv(path::AbstractString, rows)
	header = "vw,up,us,lt,aero_force,v_app,kite_vel,aoa,set_elevation,final_elevation,azimuth,proj_area,ref_area,g_earth"
	open(path, "w") do io
		println(io, header)
		for r in rows
			println(io, join([
				r.vw, r.up, r.us, r.lt, r.aero_force, r.v_app, r.kite_vel, r.aoa,
				r.set_elevation, r.final_elevation, r.azimuth, r.proj_area, r.ref_area, r.g_earth
			], ","))
		end
	end
end

function main()
	batch_name = isempty(ARGS) ? "" : strip(ARGS[1])
    batch_name = "zenith_2025_batch_2026_01_10_20_50_34"
	if isempty(batch_name)
		print("Enter batch folder name (e.g. batch_2026_01_07_10_04_38): ")
		batch_name = strip(readline())
	end
	isempty(batch_name) && error("Batch folder name is required.")

	batch_dir = joinpath("processed_data", "v3_kite", batch_name)
	log_names = find_log_names(batch_dir)
	isempty(log_names) && error("No logs found in: $batch_dir")

	rows = NamedTuple[]
	sys_cache = Dict{Tuple{Int, Int, Float64}, SystemStructure}()

	for log_name in log_names
		tags = parse_tags(log_name)
		tags === nothing && continue
		up, us, vw, lt, set_elev, g_earth = tags
		# Use g_earth from filename when available, else YAML default (0)
		g_eff = isnothing(g_earth) ? 0.0 : g_earth
		key = (vw, lt, g_eff)
		sys = get!(sys_cache, key) do
			build_sys(v_wind=float(vw), tether_length=float(lt), g_earth=g_earth)
		end
		lg = KiteUtils.load_log(log_name; path=batch_dir)
		metrics = analyze_log(lg, sys)
		push!(rows, (
			vw=vw,
			up=up,
			us=us,
			lt=lt,
			aero_force=metrics.aero_force,
			v_app=metrics.v_app,
			kite_vel=metrics.kite_vel,
			aoa=metrics.aoa,
			set_elevation=isnothing(set_elev) ? NaN : set_elev,
			final_elevation=metrics.final_elevation,
			azimuth=metrics.azimuth,
			proj_area=metrics.proj_area,
			ref_area=metrics.ref_area,
			g_earth=g_eff,
		))
	end

	sort!(rows, by=r -> (r.vw, r.up, r.us, r.lt, r.set_elevation, r.g_earth))

	out_path = joinpath(batch_dir, "zenith_circle_batch_analysis.csv")
	write_csv(out_path, rows)
	@info "Wrote batch analysis CSV" path=out_path rows=length(rows)
end


main()


nothing
