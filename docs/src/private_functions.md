```@meta
CurrentModule = StableAWEModels
```

# Private API


This page documents the internal functions and types of `StableAWEModels.jl`. These are not
part of the public API and may change without notice. They are listed here for
developers and for those interested in the model's internal workings.

## Core types and constructors

```@docs
StableAWEModels.SerializedModel
StableAWEModels.SimFloat
StableAWEModels.KVec3
StableAWEModels.SVec3
StableAWEModels.InplaceGetter
StableAWEModels.ScatterGroup
StableAWEModels.create_vsm_wing
StableAWEModels.build_vsm_engine
```

## State management and model simplification

```@docs
StableAWEModels.copy!
StableAWEModels.reinit!
StableAWEModels.reposition!
StableAWEModels.update_sys_struct!
StableAWEModels.get_set_hash
StableAWEModels.get_sys_struct_hash
```

## Physics and geometry helpers

```@docs
StableAWEModels.WindFactor
StableAWEModels.WindFactorReader
StableAWEModels.calc_angle_of_attack
StableAWEModels.calc_heading
StableAWEModels.calc_R_t_to_w
StableAWEModels.calc_R_v_to_w
StableAWEModels.cad_to_body_frame
StableAWEModels.calc_pos
StableAWEModels.calc_winch_force
StableAWEModels.quaternion_to_rotation_matrix
StableAWEModels.rotation_matrix_to_quaternion
StableAWEModels.rotate_v_around_k
StableAWEModels.smooth_norm
StableAWEModels.smooth_normalize
StableAWEModels.apply_heading
StableAWEModels.get_rot_pos
StableAWEModels.get_base_pos
StableAWEModels.calc_aoa
```

## Equations and system management

```@docs
StableAWEModels.create_sys!
StableAWEModels.scalar_eqs!
StableAWEModels.wing_eqs!
StableAWEModels.rigid_body_eqs!
StableAWEModels.body_eqs!
StableAWEModels.joint_eqs!
StableAWEModels.timoshenko_joint_eqs!
StableAWEModels.init_rigid_body!
StableAWEModels.n_orient_frames
StableAWEModels.aero_eqs!
StableAWEModels.point_eqs!
StableAWEModels.segment_eqs!
StableAWEModels.refresh_aero!
StableAWEModels.sync_aero_density!
StableAWEModels.jacobian
StableAWEModels.load_serialized_model!
StableAWEModels.maybe_create_lin_prob!
StableAWEModels.maybe_create_control_functions!
StableAWEModels.maybe_create_prob!
StableAWEModels.has_custom_component
StableAWEModels.generate_control_funcs
StableAWEModels.generate_lin_getters
StableAWEModels.generate_prob_getters
StableAWEModels.scatter_spec
StableAWEModels.build_inplace_getter
StableAWEModels.build_grouped_views
StableAWEModels.copy_vec!
StableAWEModels.scatter_component
StableAWEModels.scatter_groups
StableAWEModels.LinProbWithAttributes
StableAWEModels.ProbWithAttributes
StableAWEModels.ControlFuncWithAttributes
```

## Utility and internal functions

```@docs
StableAWEModels.get_model_name
StableAWEModels.calc_height
StableAWEModels.pos
StableAWEModels.spring_forces
StableAWEModels.create_model_archive
StableAWEModels.filecmp
StableAWEModels.extract_model_archive
StableAWEModels.copy_bin
StableAWEModels.copy_examples
StableAWEModels.copy_data
StableAWEModels.copy_dir
StableAWEModels.get_example_packages
StableAWEModels.make_lin_sys_state
```

## Base overloads (internal use)

```@docs
StableAWEModels.SAM_FIELDS
Base.getindex
Base.getproperty
Base.setproperty!
Serialization.serialize(::Serialization.AbstractSerializer, ::StableAWEModels.InplaceGetter)
```

## YAML loader internals

```@docs
StableAWEModels.get_field_or_nothing
StableAWEModels.convert_to_type
StableAWEModels.resolve_references
StableAWEModels.calculate_derived_properties!
StableAWEModels.extract_args
StableAWEModels.call_yaml_constructor
StableAWEModels.parse_tether_init
```

## SystemStructure internals

```@docs
StableAWEModels.segment_cad_length
StableAWEModels.segment_world_length
StableAWEModels.tether_ordered_point_idxs
StableAWEModels.tether_anchor_free
StableAWEModels.rigid_point_siblings
StableAWEModels.tether_downstream_idxs
StableAWEModels.twist_surface_tethers_by_overlap
StableAWEModels.tether_unit_stiffness
StableAWEModels.apply_cluster_init_stretched_len!
StableAWEModels.apply_tether_init_stretched_lens!
StableAWEModels.init_unstretched_len
StableAWEModels.apply_tether_init_forces!
StableAWEModels.joint_endpoint_frames
StableAWEModels.init_joint_rest!
StableAWEModels.timoshenko_element_frame
StableAWEModels.assign_indices_and_resolve!
StableAWEModels.resolve_ref
StableAWEModels.resolve_ref_spec
StableAWEModels.validate_sys_struct
StableAWEModels.build_name_dict
StableAWEModels.setup_wing_frame!
StableAWEModels.auto_create_twist_surfaces!
StableAWEModels.compute_twist_surface_geometry!
StableAWEModels.setup_particle_point_mapping!
StableAWEModels.identify_wing_segments
StableAWEModels.match_aero_sections_to_structure!
StableAWEModels.compute_spatial_twist_surface_mapping!
StableAWEModels.copy_cad_to_world!
StableAWEModels.adjust_vsm_panels_to_origin!
StableAWEModels.apply_aero_z_offset!
StableAWEModels.calc_particle_dynamics_wing_frame
StableAWEModels.principal_frame
StableAWEModels.calc_inertia_y_rotation
StableAWEModels.PrincipalFrameMethod
StableAWEModels.init_principal_state!
StableAWEModels.is_wing
StableAWEModels.wing_dynamics
StableAWEModels.WingDynamics
StableAWEModels.RigidDynamics
StableAWEModels.ParticleDynamics
StableAWEModels.rotate_vsm_sections!
StableAWEModels.AERO_SCALE_CHORD
StableAWEModels.body_vsm_engine
StableAWEModels.expand_auto_tethers!
StableAWEModels.WeightedRefPoints
StableAWEModels.resolve!
StableAWEModels.validate_weights!
StableAWEModels.SegmentType
```

## NamedCollection internals

```@docs
StableAWEModels.names
StableAWEModels.get_idx
StableAWEModels.get_name
Base.keys(::NamedCollection)
Base.values(::NamedCollection)
Base.haskey(::NamedCollection, ::Symbol)
Base.setindex!(::NamedCollection, ::Any, ::Integer)
Base.setindex!(::NamedCollection, ::Any, ::Symbol)
```

## Equation builders

```@docs
StableAWEModels.tether_eqs!
StableAWEModels.pulley_eqs!
StableAWEModels.winch_eqs!
StableAWEModels.twist_surface_eqs!
StableAWEModels.validate_twist_surface_modes
```

## Plate aerodynamics internals

```@docs
StableAWEModels.load_plate_wing
StableAWEModels.plate_corners
```

## Aero-mode interface

```@docs
StableAWEModels.vsm_engine
StableAWEModels.has_vsm_engine
StableAWEModels.require_vsm_engine
StableAWEModels.couples_to_sections
StableAWEModels.has_vsm_wing
StableAWEModels.provides_aero_override
StableAWEModels.stores_point_force
StableAWEModels.aero_mode_tag
StableAWEModels.calc_side_slip
StableAWEModels.validate_aero_component
StableAWEModels.validate_aero_structure
StableAWEModels.remake_aero!
StableAWEModels.setup_aero!
StableAWEModels.attach_engine!
StableAWEModels.resize_aero_state!
StableAWEModels.init_aero_state!
StableAWEModels.normalized_inertia
StableAWEModels.normalized_point_inertia
StableAWEModels.n_aero_log_points
StableAWEModels.write_aero_log_points!
StableAWEModels.read_aero_log_points!
StableAWEModels.restore_aero_twist!
StableAWEModels.plot_wing_aero!
StableAWEModels.update_wing_aero_plot!
StableAWEModels.load_wing
StableAWEModels.yaml_n_unrefined_sections
```

## VSM and aerodynamics internals

```@docs
StableAWEModels.refresh_rigid_aero!
StableAWEModels.refresh_particle_aero!
StableAWEModels.count_aero_log_points
StableAWEModels.build_point_to_vsm_point_mapping
StableAWEModels.update_vsm_wing_from_structure!
StableAWEModels.distribute_panel_forces_to_points!
StableAWEModels.rigid_aero_baseline!
StableAWEModels.apply_direct_forces!
StableAWEModels.vsm_aero_coeffs
StableAWEModels.vsm_solve_objects
StableAWEModels.safe_vsm_solve!
StableAWEModels.finite_full
StableAWEModels.set_particle_panel_va!
StableAWEModels.build_mesh_maps!
StableAWEModels.store_billow_offsets!
StableAWEModels.store_induced_velocity!
StableAWEModels.reconstruct_sections_b
StableAWEModels.ContinuousPolar
```

## Heading and geometry

```@docs
StableAWEModels.solve_heading_rotation
StableAWEModels.get_ref_position_from_points
StableAWEModels.sym_calc_R_t_to_w
StableAWEModels.wrap_to_pi
```

## Transform internals

```@docs
StableAWEModels.apply_azimuth_elevation!
StableAWEModels.apply_heading!
StableAWEModels.finalize_transforms!
```

## Flat parameters

```@docs
StableAWEModels.read_path
StableAWEModels.PathReader
StableAWEModels.ParamEntry
StableAWEModels.ParamRegistry
StableAWEModels.make_param
StableAWEModels.make_array_param
StableAWEModels.make_callable_param
StableAWEModels.leaf_param!
StableAWEModels.param_computed!
StableAWEModels.param_descend
StableAWEModels.ParamView
StableAWEModels.PathView
StableAWEModels.ParamGroup
StableAWEModels.ParamSync
StableAWEModels.survivor_index
StableAWEModels.build_param_sync
StableAWEModels.sync_params!
StableAWEModels.joint_stiffness_term
StableAWEModels.timoshenko_rigidity
```

## Initial conditions

```@docs
StableAWEModels.InitialEntry
StableAWEModels.InitialRegistry
StableAWEModels.InitialView
StableAWEModels.InitialPath
StableAWEModels.bind_initial!
StableAWEModels.ElementReader
StableAWEModels.InitialSync
StableAWEModels.build_initial_sync
StableAWEModels.sync_initial!
```

## Other internals

```@docs
StableAWEModels.init_principal_frame!
StableAWEModels.init_body_frame_from_ref_points!
StableAWEModels.get_rot_pos_cad
KiteUtils.Logger(::SymbolicAWEModel, ::Int64)
```