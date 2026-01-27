# SPDX-FileCopyrightText: 2025 Bart van de Lint <bart@vandelint.net>
#
# SPDX-License-Identifier: MPL-2.0

using SymbolicAWEModels, KiteUtils, Printf, LaTeXStrings

# --- Setup Models ---
F_step = 1.0
set_data_path("data")
set = Settings("system.yaml")
set.sample_freq = 600
set.abs_tol = 1e-6
set.rel_tol = 1e-6
# set.segments = 3
one_seg_sam = SymbolicAWEModel(set, "ram")
init!(one_seg_sam)
one_seg_tether_sam = SymbolicAWEModel(set, "tether")
init!(one_seg_tether_sam)

# --- Calculate Properties and Get Step Response Data ---
find_steady_state!(one_seg_sam)
unit_stiffness, unit_damping, tether_lens, dt = 
    SymbolicAWEModels.calc_spring_props(one_seg_sam, one_seg_tether_sam; F_step)

# --- Print Comparison Table ---
segments = one_seg_sam.sys_struct.segments
tethers = one_seg_sam.sys_struct.tethers
segments = [segments[tether.segment_idxs[1]] for tether in tethers]
real_unit_stiffness = [segment.unit_stiffness for segment in segments]
real_unit_damping = [segment.unit_damping for segment in segments]

println("\n--- Tether Spring Properties ---")
@printf "%-8s | %-15s %-15s %-10s | %-15s %-15s %-10s\n" "Tether" "Calc. Stiffness" "Real Stiffness" "Error (%)" "Calc. Damping" "Real Damping" "Error (%)"
println(repeat("-", 100))
for i in 1:4
    # Calculate relative errors in percent
    stiffness_err = 100 * abs(unit_stiffness[i] - real_unit_stiffness[i]) / real_unit_stiffness[i]
    damping_err   = 100 * abs(unit_damping[i] - real_unit_damping[i]) / real_unit_damping[i]
    # Print data rows
    @printf "%-8d | %-15.2f %-15.2f %-10.2f | %-15.2f %-15.2f %-10.2f\n" i unit_stiffness[i] real_unit_stiffness[i] stiffness_err unit_damping[i] real_unit_damping[i] damping_err
end

# --- Print summary ---
steps = size(tether_lens, 2) - 1
println("\nStep response computed over $(steps) steps with dt=$(dt)")
println("Final tether length: $(tether_lens[1,end]) m")
