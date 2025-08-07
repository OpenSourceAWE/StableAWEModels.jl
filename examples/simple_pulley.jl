# SPDX-FileCopyrightText: 2025 Bart van de Lint
#
# SPDX-License-Identifier: MPL-2.0

using SymbolicAWEModels, VortexStepMethod, ControlPlots

set = load_settings("base")  # Loads from data/base/settings.yaml
set.v_wind = 0
set.l_tether = 5.0 # set l_tether as it affects the plot size
dynamics_type = DYNAMIC

# pulley point was placed in the middle, as the side-to-sideo movement convergences very slowly
points = Point[]
push!(points, Point(1, [0.0, 0.0, 5.0], STATIC))
push!(points, Point(2, [5.0, 0.0, 5.0], STATIC))
push!(points, Point(3, [2.5, 0.0, 1], DYNAMIC; mass=2))

segments = Segment[]
push!(segments, Segment(1, set, (3,1), BRIDLE,; l0=3.5, compression_frac=0.01, diameter_mm=5))  # 5mm diameter, 4m rest length
push!(segments, Segment(2, set, (3,2), BRIDLE,; l0=3.5, compression_frac=0.01, diameter_mm=5))  # 5mm diameter, 4m rest length

pulleys = Pulley[]
push!(pulleys, Pulley(1, (1,2), DYNAMIC))

transforms = [Transform(1, -deg2rad(0.0), 0.0, 0.0; base_pos=[0.0, 0.0, 5.0], base_point_idx=1, rot_point_idx=2)]
sys_struct = SymbolicAWEModels.SystemStructure("pulley", set; points, segments, pulleys, transforms)
plot(sys_struct, 0.0; zoom=false, l_tether=set.l_tether)

sam = SymbolicAWEModel(set, sys_struct)

init!(sam; remake=false)
for i in 1:100
    plot(sam, i/set.sample_freq; zoom=false)
    next_step!(sam)
end
