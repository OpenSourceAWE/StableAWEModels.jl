# run_pyramid_model.jl
#
# Purpose:
#   - Load pyramid_model structural YAML
#   - Compute panel normals & nodal loads (perpendicular to each panel)
#   - Build system, init model
#   - Run a few simulation steps and plot each iteration
#
# Requires: YAML.jl, SymbolicAWEModels, VortexStepMethod, ControlPlots
# Optional: include("../yaml_loader.jl") if your helper lives there.

using LinearAlgebra
using YAML
using SymbolicAWEModels
using VortexStepMethod
using StaticArrays

# ============= User settings =============
const MODEL_NAME = "pyramid_model"
const GEOM_PATH  = joinpath("data", MODEL_NAME, "struc_geometry_converged.yaml")
const GEOM_PATH  = joinpath("data", MODEL_NAME, "struc_geometry.yaml")
const N_PLOTS    = 3  # number of static snapshots to render after simulation

# Aerodynamic parameters
const WIND_VECTOR = SVector(10.0, 0.0, 1.0)
const WIND_DIR    = normalize(WIND_VECTOR)
const FORCE_COEFF = 2π

# Panel topology (two panels: right & left)
const PANEL_DEFINITIONS = (
    (nodes=(1, 2, 4, 3), le_pair=(1, 3), te_pair=(2, 4)),
    (nodes=(3, 4, 6, 5), le_pair=(3, 5), te_pair=(4, 6)),
)

# External aerodynamic loads (array rows match node indices 0-6)
const F_AERO_WING = [
       0.0            0.0            0.0        ;
     -49.23422944  169.91398256  225.26646516;
     -21.28999615   70.38702941   93.1314033 ;
    -136.84008626  -16.26365694  564.1619389 ;
     -53.14868866   -7.52853902  222.3234133 ;
     -46.91474061 -185.15041031  224.70657846;
     -20.7152601   -77.89786305   94.32859107
]

# simulation / plotting knobs
n_steps      = 50
use_zoom     = false
remake_cache = false
# ========================================

# ---- helpers ----
cross3(a, b) = SVector(
    a[2]*b[3] - a[3]*b[2],
    a[3]*b[1] - a[1]*b[3],
    a[1]*b[2] - a[2]*b[1],
)

function unit_normal_quad(p1::SVector{3,<:Real}, p2::SVector{3,<:Real},
                          p3::SVector{3,<:Real}, p4::SVector{3,<:Real})
    v12 = p2 .- p1
    v23 = p3 .- p2
    v34 = p4 .- p3
    n1 = cross3(v12, v23)
    n2 = cross3(v23, v34)
    n  = n1 + n2
    norm_n = norm(n)
    norm_n < 1e-9 && error("Panel is degenerate; cannot compute normal.")
    return n / norm_n
end

to_svec3(vec) = SVector{3,Float64}(Float64(vec[1]), Float64(vec[2]), Float64(vec[3]))

function point_position(point; prefer_world::Bool=true)
    if prefer_world
        pw = point.pos_w
        if sum(abs2, pw) > 1e-12
            return to_svec3(pw)
        end
    end
    return to_svec3(point.pos_cad)
end

function compute_node_forces(points; prefer_world::Bool=true)
    node_forces = Dict{Int, SVector{3,Float64}}(Int(p.idx) => SVector{3,Float64}(0.0, 0.0, 0.0) for p in points)
    positions = Dict(Int(p.idx) => point_position(p; prefer_world=prefer_world) for p in points)

    for definition in PANEL_DEFINITIONS
        panel_nodes = definition.nodes
        le_pair = definition.le_pair
        te_pair = definition.te_pair

        p1 = positions[panel_nodes[1]]
        p2 = positions[panel_nodes[2]]
        p3 = positions[panel_nodes[3]]
        p4 = positions[panel_nodes[4]]

        n̂ = unit_normal_quad(p1, p2, p3, p4)
        mid_le = (positions[le_pair[1]] + positions[le_pair[2]]) / 2
        mid_te = (positions[te_pair[1]] + positions[te_pair[2]]) / 2
        chord_vec = mid_te - mid_le
        chord_norm = norm(chord_vec)
        chord_norm < 1e-9 && error("Mid-chord vector is degenerate for panel $(panel_nodes)")
        chord_dir = chord_vec / chord_norm

        α = atan(norm(cross3(chord_dir, WIND_DIR)), dot(chord_dir, WIND_DIR))
        lift_coeff = FORCE_COEFF * α
        panel_force = lift_coeff * n̂

        force_per_node = panel_force / 4
        for nid in panel_nodes
            node_forces[nid] = node_forces[nid] + force_per_node
        end
    end

    return node_forces
end

function apply_node_forces!(points, node_forces)
    for point in points
        nid = Int(point.idx)
        f = get(node_forces, nid, SVector{3,Float64}(0.0, 0.0, 0.0))
        point.disturb[1] = f[1]
        point.disturb[2] = f[2]
        point.disturb[3] = f[3]
    end
    return nothing
end

function make_point_static!(sys, node_id::Int; fix_sphere::Bool=true)
    point = sys.points[node_id]
    sys.points[node_id] = SymbolicAWEModels.Point(
        point.idx,
        point.transform_idx,
        point.wing_idx,
        point.pos_cad,
        point.pos_b,
        point.pos_w,
        point.vel_w,
        point.disturb,
        point.force,
        SymbolicAWEModels.STATIC,
        point.mass,
        point.body_frame_damping,
        point.world_frame_damping,
        fix_sphere ? true : point.fix_sphere,
    )
    return nothing
end

# ---- load geometry YAML ----
@info "Loading geometry from: $GEOM_PATH"
isfile(GEOM_PATH) || error("Geometry file not found at: $GEOM_PATH")

geom = open(GEOM_PATH, "r") do io
    YAML.load(io)
end

# Collect wing node coordinates (ids → SVector{3,Float64})
nodes = Dict{Int, SVector{3,Float64}}()
for row in geom["wing_particles"]["data"]
    id  = Int(row[1])
    x,y,z = Float64(row[2]), Float64(row[3]), Float64(row[4])
    nodes[id] = SVector(x,y,z)
end


# ---- settings & system ----
# Prefer package helper if present (works with your repo layout).
# If your helper is in a sibling path, adjust include() accordingly.
yaml_loader_path = joinpath(@__DIR__, "..", "yaml_loader.jl")
if isfile(yaml_loader_path)
    include(yaml_loader_path)  # provides load_sys_struct_from_yaml
else
    @warn "yaml_loader.jl not found relative to this script; expecting load_sys_struct_from_yaml to be available in scope."
end

# also include the plotly utilities
plotly_helpers_path = joinpath(@__DIR__, "..", "plotly_plots.jl")
if isfile(plotly_helpers_path)
    include(plotly_helpers_path)  # provides plot helpers
else
    @warn "plotly_plots.jl not found relative to this script; expecting plot helpers to be available in scope."
end

# Load default settings for a minimal structural run.
# If you have a project-specific settings dir, swap for load_settings(MODEL_NAME).
set = try
    SymbolicAWEModels.load_settings(MODEL_NAME)  # expects data/pyramid_model/settings.yaml etc.
catch
    # Fallback: a minimal Settings file if your project uses a base/system.yaml
    Settings("base/system.yaml")
end

# Make sure recent naming is used (KiteUtils ≥ v0.11.0)
if hasfield(typeof(set), :axial_damping) == false && hasfield(typeof(set), :damping)
    @warn "Old settings field 'damping' detected; consider upgrading to 'axial_damping'."
end


# Load system structure from YAML
@info "Building system structure from YAML…"
sys = load_sys_struct_from_yaml(GEOM_PATH; system_name=MODEL_NAME, set=set)

# Make node 3 static/fixed in the structural system
# make_point_static!(sys, 3)
# @info "Marked node 3 as STATIC"

# Initial 3D plot of the structure (t=0)
# plt = plot3d_v3(sys.points, sys.segments; title="Pyramid Model (Initial)")


@info "System summary: points=$(length(sys.points)), segments=$(length(sys.segments)), pulleys=$(length(sys.pulleys))"
fixed_count   = count(p -> p.type == SymbolicAWEModels.STATIC, sys.points)
dynamic_count = count(p -> p.type == SymbolicAWEModels.DYNAMIC, sys.points)
@info "Fixed points: $fixed_count, Dynamic points: $dynamic_count"

# Apply provided static loads to system points
for (i, point) in enumerate(sys.points)
    # Provided loads are 0-indexed: row 1 is node 0 (ground), rows 2-7 are wing nodes 1-6
    if i <= size(F_AERO_WING, 1)
        point.disturb[1] = F_AERO_WING[i, 1]
        point.disturb[2] = F_AERO_WING[i, 2]
        point.disturb[3] = F_AERO_WING[i, 3]
    end
end
@info "Assigned provided static loads to $(length(sys.points)) points"


# ---- model, init, simulate, and collect states ----
sam = SymbolicAWEModel(set, sys)

# Ensure provided static loads are set on the model structure before initialization
for (i, point) in enumerate(sam.sys_struct.points)
    if i <= size(F_AERO_WING, 1)
        point.disturb[1] = F_AERO_WING[i, 1]
        point.disturb[2] = F_AERO_WING[i, 2]
        point.disturb[3] = F_AERO_WING[i, 3]
    end
end

# (Optional) if you want to set basic solver tolerances / rates:
# set.abs_tol = 1e-3
# set.rel_tol = 1e-3
# set.sample_freq = 50  # if not already set in settings

init!(sam; remake=remake_cache)

# Re-apply provided static loads after initialization
for (i, point) in enumerate(sam.sys_struct.points)
    if i <= size(F_AERO_WING, 1)
        point.disturb[1] = F_AERO_WING[i, 1]
        point.disturb[2] = F_AERO_WING[i, 2]
        point.disturb[3] = F_AERO_WING[i, 3]
    end
end

# Time loop - collect snapshots for static plots
Δt = 1.0 / max(1, hasproperty(set, :sample_freq) ? set.sample_freq : 100)
@info "Running simulation for $n_steps steps (Δt = $(round(Δt, digits=4)) s)..."

# Determine which steps to capture (include start and end)
num_samples = max(N_PLOTS, 2)
snapshot_steps = unique!(sort!(round.(Int, range(0, stop=n_steps, length=num_samples))))
snapshot_steps[1] != 0 && pushfirst!(snapshot_steps, 0)
snapshot_steps[end] != n_steps && push!(snapshot_steps, n_steps)

snapshots = Dict{Int, Vector{Point}}(0 => deepcopy(sam.sys_struct.points))

for step in 1:n_steps
    t = step * Δt

    # Re-apply provided static loads at each step
    for (i, point) in enumerate(sam.sys_struct.points)
        if i <= size(F_AERO_WING, 1)
            point.disturb[1] = F_AERO_WING[i, 1]
            point.disturb[2] = F_AERO_WING[i, 2]
            point.disturb[3] = F_AERO_WING[i, 3]
        end
    end

    # Advance simulation one step
    next_step!(sam)

    # Store current state if requested
    if step in snapshot_steps
        snapshots[step] = deepcopy(sam.sys_struct.points)
    end

    # Optional: print progress
    if step % 20 == 0
        @info "Step $step/$n_steps (t = $(round(t, digits=2)) s)"
    end
end

# Ensure final state is captured
snapshots[n_steps] = get(snapshots, n_steps, deepcopy(sam.sys_struct.points))

captured_steps = sort!(collect(keys(snapshots)))

@info "Simulation complete. Rendering $(length(captured_steps)) static plots..."

for (idx, step) in enumerate(captured_steps)
    points_snapshot = snapshots[step]
    t = step * Δt
    plot_title = "Pyramid Model – Step $(step) (t=$(round(t, digits=2)) s)"
    plot3d_v3(points_snapshot, sam.sys_struct.segments; title=plot_title)
    @info "Rendered static snapshot $(idx)/$(length(snapshot_steps)) at step $step"
end

@info "Done. Created $(length(snapshot_steps)) static Plotly plots covering the simulation window."
