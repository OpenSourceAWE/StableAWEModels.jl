
function SystemStructure(set::Settings)
    if set.physical_model == "ram"
        return create_ram_sys_struct(set)
    elseif set.physical_model == "simple_ram"
        return create_simple_ram_sys_struct(set)
    else
        throw(ArgumentError("Undefined physical model"))
    end
end

function create_ram_sys_struct(set::Settings)
    vsm_wing = RamAirWing(set; prn=false)
    points = Point[]
    groups = Group[]
    segments = Segment[]
    pulleys = Pulley[]
    tethers = Tether[]
    winches = Winch[]
    wings = Wing[]

    attach_points = Point[]
    
    bridle_top_left = [cad_to_body_frame(vsm_wing, set.top_bridle_points[i]) for i in eachindex(set.top_bridle_points)]
    bridle_top_right = [bridle_top_left[i] .* [1, -1, 1] for i in eachindex(set.top_bridle_points)]

    dynamics_type = set.quasi_static ? QUASI_STATIC : DYNAMIC
    z = vsm_wing.R_cad_body[:,3]

    function create_bridle(bridle_top, gammas)
        i_pnt = length(points) # last point idx
        i_seg = length(segments) # last segment idx
        i_pul = length(pulleys) # last pulley idx
        i_grp = length(groups) # last group idx

        # ==================== CREATE DEFORMING WING GROUPS ==================== #
        points = [
            points
            Point(1+i_pnt, calc_pos(vsm_wing, gammas[1], set.bridle_fracs[1]), WING)
            Point(2+i_pnt, calc_pos(vsm_wing, gammas[1], set.bridle_fracs[3]), WING)
            Point(3+i_pnt, calc_pos(vsm_wing, gammas[1], set.bridle_fracs[4]), WING)

            Point(4+i_pnt, calc_pos(vsm_wing, gammas[2], set.bridle_fracs[1]), WING)
            Point(5+i_pnt, calc_pos(vsm_wing, gammas[2], set.bridle_fracs[3]), WING)
            Point(6+i_pnt, calc_pos(vsm_wing, gammas[2], set.bridle_fracs[4]), WING)
        ]
        groups = [
            groups
            Group(1+i_grp, [1+i_pnt, 2+i_pnt, 3+i_pnt], vsm_wing, gammas[1], DYNAMIC, set.bridle_fracs[2])
            Group(2+i_grp, [4+i_pnt, 5+i_pnt, 6+i_pnt], vsm_wing, gammas[2], DYNAMIC, set.bridle_fracs[2])
        ]

        # ==================== CREATE PULLEY BRIDLE SYSTEM ==================== #
        bridle_damping = 1.0
        points = [
            points
            Point(7+i_pnt, bridle_top[1], dynamics_type; bridle_damping)
            Point(8+i_pnt, bridle_top[2], WING)
            Point(9+i_pnt, bridle_top[3], dynamics_type; bridle_damping)
            Point(10+i_pnt, bridle_top[4], dynamics_type; bridle_damping)

            Point(11+i_pnt, bridle_top[2] - 1z, dynamics_type; bridle_damping)

            Point(12+i_pnt, bridle_top[1] - 2z, dynamics_type; bridle_damping)
            Point(13+i_pnt, bridle_top[3] - 2z, dynamics_type; bridle_damping)

            Point(14+i_pnt, bridle_top[1] - 4z, dynamics_type; bridle_damping)
            Point(15+i_pnt, bridle_top[3] - 4z, dynamics_type; bridle_damping)
        ]
        segments = [
            segments
            Segment(1+i_seg, set, (1+i_pnt, 7+i_pnt), BRIDLE)
            Segment(2+i_seg, set, (2+i_pnt, 9+i_pnt), BRIDLE)
            Segment(3+i_seg, set, (3+i_pnt, 10+i_pnt), BRIDLE)

            Segment(4+i_seg, set, (4+i_pnt, 7+i_pnt), BRIDLE)
            Segment(5+i_seg, set, (5+i_pnt, 9+i_pnt), BRIDLE)
            Segment(6+i_seg, set, (6+i_pnt, 10+i_pnt), BRIDLE)

            Segment(7+i_seg, set, (7+i_pnt, 12+i_pnt), BRIDLE; l0=2)
            Segment(8+i_seg, set, (8+i_pnt, 11+i_pnt), BRIDLE; l0=1)
            Segment(9+i_seg, set, (9+i_pnt, 13+i_pnt), BRIDLE; l0=2)
            Segment(10+i_seg, set, (10+i_pnt, 15+i_pnt), BRIDLE; l0=4)
            
            Segment(11+i_seg, set, (11+i_pnt, 12+i_pnt), BRIDLE; l0=1)
            Segment(12+i_seg, set, (11+i_pnt, 13+i_pnt), BRIDLE; l0=1)
            
            Segment(13+i_seg, set, (12+i_pnt, 14+i_pnt), BRIDLE; l0=2)
            Segment(14+i_seg, set, (13+i_pnt, 14+i_pnt), BRIDLE; l0=2)
            Segment(15+i_seg, set, (13+i_pnt, 15+i_pnt), BRIDLE; l0=2)
        ]
        pulleys = [
            pulleys
            Pulley(1+i_pul, (11+i_seg, 12+i_seg), dynamics_type)
            Pulley(2+i_pul, (14+i_seg, 15+i_seg), dynamics_type)
        ]
        push!(attach_points, points[end-1])
        push!(attach_points, points[end])
        return nothing
    end

    gammas = [-3/4, -1/4, 1/4, 3/4] * vsm_wing.gamma_tip
    create_bridle(bridle_top_left, gammas[[1,2]])
    create_bridle(bridle_top_right, gammas[[3,4]])

    points, segments, tethers, left_power_idx =
        create_tether(1, set, points, segments, tethers, attach_points[1], 
                      POWER_LINE, dynamics_type, z)
    points, segments, tethers, right_power_idx =
        create_tether(2, set, points, segments, tethers, attach_points[3], 
                      POWER_LINE, dynamics_type, z)
    points, segments, tethers, left_steering_idx =
        create_tether(3, set, points, segments, tethers, attach_points[2], 
                      STEERING_LINE, dynamics_type, z)
    points, segments, tethers, right_steering_idx =
        create_tether(4, set, points, segments, tethers, attach_points[4], 
                      STEERING_LINE, dynamics_type, z)

    winches = [
        Winch(1, TorqueControlledMachine(set), [left_power_idx, right_power_idx])
        Winch(2, TorqueControlledMachine(set), [left_steering_idx])
        Winch(3, TorqueControlledMachine(set), [right_steering_idx])
    ]

    vsm_aero = BodyAerodynamics([vsm_wing])
    vsm_solver = Solver(vsm_aero; solver_type=NONLIN, atol=2e-8, rtol=2e-8)
    wings = [Wing(1, vsm_aero, vsm_wing, vsm_solver, [1,2,3,4], I(3), zeros(3))]
    transforms = [Transform(1, deg2rad(set.elevation), deg2rad(set.azimuth), deg2rad(set.heading);
                                    base_pos= zeros(3), base_point_idx=points[end].idx, wing_idx=1)]
    
    return SystemStructure(set.physical_model, set; points, groups, segments, pulleys, tethers, winches, wings, transforms)
end

function create_simple_ram_sys_struct(set::Settings)
    vsm_wing = RamAirWing(set)
    gammas = [-1/2, 1/2] * wing.gamma_tip
    
    bridle_top_left = [vsm_wing.R_cad_body * (set.top_bridle_points[i] + vsm_wing.T_cad_body) for i in eachindex(set.top_bridle_points)] # cad to kite frame
    bridle_top_right = [bridle_top_left[i] .* [1, -1, 1] for i in eachindex(set.top_bridle_points)]

    points = [
        Point(1, bridle_top_left[2], WING)
        Point(2, bridle_top_right[2], WING)
        Point(3, calc_pos(vsm_wing, gammas[1], set.bridle_fracs[4]), WING)
        Point(4, calc_pos(vsm_wing, gammas[2], set.bridle_fracs[4]), WING)

        Point(5, [0, 0, -set.l_tether], STATIC)
        Point(6, [0, 0, -set.l_tether], STATIC)
        Point(7, [0, 0, -set.l_tether], STATIC)
        Point(8, [0, 0, -set.l_tether], STATIC)
    ]
    groups = [
        Group(1, [3], vsm_wing, gammas[1], DYNAMIC, set.bridle_fracs[2])
        Group(2, [4], vsm_wing, gammas[2], DYNAMIC, set.bridle_fracs[2])
    ]
    segments = [
        Segment(1, set, (1,5), POWER_LINE)
        Segment(2, set, (2,6), POWER_LINE)
        Segment(3, set, (3,7), STEERING_LINE)
        Segment(4, set, (4,8), STEERING_LINE)
    ]
    tethers = [
        Tether(1, [1], 5)
        Tether(2, [2], 6)
        Tether(3, [3], 7)
        Tether(4, [4], 8)
    ]
    winches = [
        Winch(1, TorqueControlledMachine(set), [1,2])
        Winch(2, TorqueControlledMachine(set), [3])
        Winch(3, TorqueControlledMachine(set), [4])
    ]
    vsm_aero = BodyAerodynamics([vsm_wing])
    vsm_solver = Solver(vsm_aero; solver_type=NONLIN, atol=2e-8, rtol=2e-8)
    wings = [Wing(1, vsm_aero, vsm_wing, vsm_solver, [1,2], I(3), zeros(3))]
    transforms = [
        Transform(1, deg2rad(set.elevation), deg2rad(set.azimuth), deg2rad(set.heading);
                    base_pos=zeros(3), base_point_idx=5, wing_idx=1)
    ]

    return SystemStructure(set.physical_model, set; 
        points, groups, segments, tethers, winches, wings, transforms)
end

function update_simple_sam!(ssam::SymbolicAWEModel, sam::SymbolicAWEModel)
    
end

