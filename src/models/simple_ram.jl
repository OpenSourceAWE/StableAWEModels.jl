function create_simple_ram_sys_struct(set::Settings, wing::RamAirWing)
    points = Point[]
    groups = Group[]
    segments = Segment[]
    pulleys = Pulley[]
    tethers = Tether[]
    winches = Winch[]

    dynamics_type = set.quasi_static ? QUASI_STATIC : DYNAMIC
    gammas = [-3/4, -1/4, 1/4, 3/4] * wing.gamma_tip
    
    bridle_top_left = [wing.R_cad_body * (set.top_bridle_points[i] + wing.T_cad_body) for i in eachindex(set.top_bridle_points)] # cad to kite frame
    bridle_top_right = [bridle_top_left[i] .* [1, -1, 1] for i in eachindex(set.top_bridle_points)]

    # ==================== CREATE DEFORMING WING GROUPS ==================== #
    points = [
        points
        Point(1, calc_pos(wing, gammas[1], set.bridle_fracs[4]), WING)
        Point(2, calc_pos(wing, gammas[2], set.bridle_fracs[4]), WING)
        Point(3, calc_pos(wing, gammas[3], set.bridle_fracs[4]), WING)
        Point(4, calc_pos(wing, gammas[4], set.bridle_fracs[4]), WING)
    ]
    groups = [
        groups
        Group(1, [1], wing, gammas[1], DYNAMIC, set.bridle_fracs[2])
        Group(2, [2], wing, gammas[2], DYNAMIC, set.bridle_fracs[2])
        Group(3, [3], wing, gammas[3], DYNAMIC, set.bridle_fracs[2])
        Group(4, [4], wing, gammas[4], DYNAMIC, set.bridle_fracs[2])
    ]
    # ==================== CREATE PULLEY BRIDLE SYSTEM ==================== #
    points = [
        points
        Point(5, bridle_top_left[2], WING)
        Point(6, bridle_top_left[4], dynamics_type)
        Point(7, bridle_top_right[2], WING)
        Point(8, bridle_top_right[4], dynamics_type)
    ]

    segments = [
        segments
        Segment(1, (1, 6), BRIDLE)
        Segment(2, (2, 6), BRIDLE)
        Segment(3, (3, 8), BRIDLE)
        Segment(4, (4, 8), BRIDLE)
    ]

    points, segments, tethers, left_power_idx = create_tether(1, set, points, segments, tethers, points[5], POWER_LINE, dynamics_type)
    points, segments, tethers, right_power_idx = create_tether(2, set, points, segments, tethers, points[7], POWER_LINE, dynamics_type)
    points, segments, tethers, left_steering_idx = create_tether(3, set, points, segments, tethers, points[6], STEERING_LINE, dynamics_type)
    points, segments, tethers, right_steering_idx = create_tether(4, set, points, segments, tethers, points[8], STEERING_LINE, dynamics_type)

    winches = [winches; Winch(1, TorqueControlledMachine(set), [left_power_idx, right_power_idx])]
    winches = [winches; Winch(2, TorqueControlledMachine(set), [left_steering_idx])]
    winches = [winches; Winch(3, TorqueControlledMachine(set), [right_steering_idx])]

    wings = [Wing(1, [1,2,3,4], I(3), zeros(3))]
    transforms = [Transform(1, deg2rad(set.elevation), deg2rad(set.azimuth), deg2rad(set.heading);
                                    base_pos= zeros(3), base_point_idx=points[end].idx, wing_idx=1)]

    return SystemStructure(set.physical_model, set; points, groups, segments, pulleys, tethers, winches, wings, transforms)
end

