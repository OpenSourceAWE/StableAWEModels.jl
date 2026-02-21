# SPDX-FileCopyrightText: 2025 Bart van de Lint
# SPDX-License-Identifier: MPL-2.0

# Benchmark ODE RHS for 2plate_kite model.
# Tests: raw RHS, registered functions directly, allocation profile.
#
# Usage: julia --project=test test/test_bench.jl

using SymbolicAWEModels
using SymbolicAWEModels: VortexStepMethod, SystemStructure, KVec3
using KiteUtils
using BenchmarkTools
using Statistics
using Printf
using Profile

function setup_sam()
    pkg_root = dirname(@__DIR__)
    src_data = joinpath(pkg_root, "data", "2plate_kite")

    tmpdir = mktempdir()
    data_path = joinpath(tmpdir, "2plate_kite")
    cp(src_data, data_path; force=true)

    set_data_path(data_path)
    set = Settings("system.yaml")

    vsm_set = VortexStepMethod.VSMSettings(
        joinpath(data_path, "vsm_settings.yaml");
        data_prefix=false)

    struc_yaml = joinpath(data_path,
                          "quat_struc_geometry.yaml")
    sys = load_sys_struct_from_yaml(struc_yaml;
        system_name="bench", set, vsm_set)
    sys.winches[:main_winch].brake = true

    sam = SymbolicAWEModel(set, sys)
    init!(sam; remake=false, prn=true)
    return sam
end

# --- 1. RHS Benchmark ---
function bench_rhs(sam)
    f = sam.integrator.f
    u = copy(sam.integrator.u)
    p = sam.integrator.p
    t = sam.integrator.t
    du = similar(u)
    f(du, u, p, t)
    return @benchmark $f($du, $u, $p, $t) samples=10
end

# --- 2. Direct registered function calls ---
function bench_registered_functions(sys::SystemStructure)
    idx = Int64(1)
    # Warmup
    SymbolicAWEModels.get_pos_w(sys, idx)
    SymbolicAWEModels.get_Q_b_to_w(sys, idx)
    SymbolicAWEModels.get_l0(sys, idx)
    SymbolicAWEModels.get_extra_mass(sys, idx)

    header = @sprintf("  %-20s %6s\n", "Function", "Allocs")
    sep = "  " * "-"^30

    # Raw field access
    println("\n  Raw field access (no @register):")
    println(sep)
    print(header)
    println(sep)

    a = @allocations sys.segments[idx].l0
    @printf("  %-20s %6d\n", "seg.l0", a)

    a = @allocations sys.points[idx].extra_mass
    @printf("  %-20s %6d\n", "pt.extra_mass", a)

    a = @allocations sys.points[idx].pos_w
    @printf("  %-20s %6d\n", "pt.pos_w", a)

    a = @allocations sys.wings[idx].Q_b_to_w
    @printf("  %-20s %6d\n", "wing.Q_b_to_w", a)

    a = @allocations sys.wings[idx].com_w
    @printf("  %-20s %6d\n", "wing.com_w", a)

    a = @allocations sys.wings[idx].R_b_to_p
    @printf("  %-20s %6d\n", "wing.R_b_to_p", a)

    a = @allocations sys.wings[idx].inertia_principal
    @printf("  %-20s %6d\n", "wing.inertia_princ", a)

    # @register_symbolic calls
    println("\n  @register_symbolic calls:")
    println(sep)
    print(header)
    println(sep)

    a = @allocations SymbolicAWEModels.get_l0(sys, idx)
    @printf("  %-20s %6d\n", "get_l0", a)

    a = @allocations SymbolicAWEModels.get_extra_mass(sys, idx)
    @printf("  %-20s %6d\n", "get_extra_mass", a)

    a = @allocations SymbolicAWEModels.get_pos_w(sys, idx)
    @printf("  %-20s %6d\n", "get_pos_w", a)

    a = @allocations SymbolicAWEModels.get_Q_b_to_w(sys, idx)
    @printf("  %-20s %6d\n", "get_Q_b_to_w", a)

    a = @allocations SymbolicAWEModels.get_com_w(sys, idx)
    @printf("  %-20s %6d\n", "get_com_w", a)

    a = @allocations SymbolicAWEModels.get_R_b_to_p(sys, idx)
    @printf("  %-20s %6d\n", "get_R_b_to_p", a)

    a = @allocations SymbolicAWEModels.get_inertia_principal(
        sys, idx)
    @printf("  %-20s %6d\n", "get_inertia_princ", a)
end

# --- 3. Allocation profile ---
function profile_rhs_allocs(sam)
    f = sam.integrator.f
    u = copy(sam.integrator.u)
    p = sam.integrator.p
    t = sam.integrator.t
    du = similar(u)
    f(du, u, p, t)

    Profile.Allocs.clear()
    Profile.Allocs.@profile sample_rate=1.0 f(du, u, p, t)
    results = Profile.Allocs.fetch()

    type_counts = Dict{Any,Int}()
    type_bytes = Dict{Any,Int}()

    # Helper: is this a meaningful Julia source frame?
    function is_julia_src(frame)
        frame.line <= 0 && return false
        file = string(frame.file)
        # Skip C runtime, GC, profiler, sysimage
        contains(file, "gc-") && return false
        contains(file, "Profile") && return false
        contains(file, "datatype.c") && return false
        endswith(file, ".c") && return false
        endswith(file, ".h") && return false
        file == ":-1" && return false
        return true
    end

    # Aggregate by source location (first Julia src frame)
    loc_counts = Dict{String,Int}()
    loc_bytes = Dict{String,Int}()

    for a in results.allocs
        T = a.type
        type_counts[T] = get(type_counts, T, 0) + 1
        type_bytes[T] = get(type_bytes, T, 0) + a.size

        loc = "unknown"
        for frame in a.stacktrace
            is_julia_src(frame) || continue
            file = string(frame.file)
            # Skip base Julia math (sys.so compiled)
            endswith(file, ".so") && continue
            basename(file) in (
                "int.jl", "float.jl", "promotion.jl",
                "number.jl", "boot.jl") && continue
            loc = "$(basename(file)):$(frame.line)"
            break
        end
        loc_counts[loc] = get(loc_counts, loc, 0) + 1
        loc_bytes[loc] = get(loc_bytes, loc, 0) + a.size
    end

    # Aggregate by call stack (top 5 Julia src frames,
    # skipping base math)
    stack_counts = Dict{String,Int}()
    stack_bytes = Dict{String,Int}()
    for a in results.allocs
        frames = String[]
        for frame in a.stacktrace
            is_julia_src(frame) || continue
            fname = basename(string(frame.file))
            func = string(frame.func)
            # Truncate long generated function names
            if length(func) > 40
                func = func[1:37] * "..."
            end
            push!(frames, "$func ($fname:$(frame.line))")
            length(frames) >= 5 && break
        end
        key = isempty(frames) ? "unknown" :
            join(frames, "\n       <- ")
        stack_counts[key] = get(stack_counts, key, 0) + 1
        stack_bytes[key] = get(stack_bytes, key, 0) + a.size
    end

    return (; type_counts, type_bytes,
              loc_counts, loc_bytes,
              stack_counts, stack_bytes)
end

# --- Reporting ---
function print_rhs_bench(rhs)
    println("\n", "="^60)
    println("  ODE RHS Benchmark")
    println("="^60)
    @printf("  median: %8.3f ms\n",
            median(rhs.times) / 1e6)
    @printf("  allocs: %8d  (%d bytes)\n",
            rhs.allocs, rhs.memory)
    println("="^60)
end

function print_alloc_profile(prof)
    (; type_counts, type_bytes,
       loc_counts, loc_bytes,
       stack_counts, stack_bytes) = prof
    sorted = sort(collect(type_counts); by=last, rev=true)
    total = sum(values(type_counts))

    println("\n", "="^60)
    println("  Allocation Profile by Type ($total total)")
    println("="^60)
    @printf("  %-35s %7s %10s\n", "Type", "Count", "Bytes")
    println("  ", "-"^56)
    for (T, count) in sorted[1:min(15, end)]
        name = string(T)
        if length(name) > 35
            name = name[1:32] * "..."
        end
        bytes = get(type_bytes, T, 0)
        @printf("  %-35s %7d %10d\n", name, count, bytes)
    end
    println("="^60)

    # --- By source location ---
    loc_sorted = sort(collect(loc_counts);
                      by=last, rev=true)
    println("\n", "="^70)
    println("  Top Allocating Source Locations")
    println("="^70)
    @printf("  %-50s %7s %10s\n",
            "Location", "Count", "Bytes")
    println("  ", "-"^67)
    for (loc, count) in loc_sorted[1:min(20, end)]
        name = loc
        if length(name) > 50
            name = "..." * name[end-46:end]
        end
        bytes = get(loc_bytes, loc, 0)
        @printf("  %-50s %7d %10d\n", name, count, bytes)
    end
    println("="^70)

    # --- By call stack ---
    stack_sorted = sort(collect(stack_counts);
                        by=last, rev=true)
    println("\n", "="^70)
    println("  Top Allocating Call Stacks (top 3 frames)")
    println("="^70)
    for (i, (stack, count)) in enumerate(
            stack_sorted[1:min(15, end)])
        bytes = get(stack_bytes, stack, 0)
        @printf("  #%-2d  %5d allocs (%d bytes)\n",
                i, count, bytes)
        for frame in split(stack, " <- ")
            println("       ", frame)
        end
        println()
    end
    println("="^70)
end

# --- Run ---
println("\n>> Setting up 2plate_kite model...")
sam = setup_sam()

println("\n>> Benchmarking ODE RHS...")
rhs = bench_rhs(sam)
print_rhs_bench(rhs)

println("\n>> Benchmarking registered functions directly...")
bench_registered_functions(sam.sys_struct)

println("\n>> Profiling RHS allocations...")
prof = profile_rhs_allocs(sam)
print_alloc_profile(prof)
