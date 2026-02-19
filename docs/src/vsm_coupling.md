# VSM coupling

This document explains how SymbolicAWEModels couples with the Vortex Step Method (VSM) for aerodynamic force computation. The coupling differs between the two wing types: [`QUATERNION` and `REFINE`](@ref WingType).

## Overview

Both wing types use **unrefined sections** as the fundamental structural element that maps to VSM geometry:

- **Unrefined section**: A structural element defined by two points (leading edge and trailing edge) along the wing span
- **Refined panel**: VSM can subdivide unrefined sections into multiple panels for higher aerodynamic fidelity
- **`refined_panel_mapping`**: Maps each refined panel back to its parent unrefined section

The VSM solver computes aerodynamic coefficients (cl, cd, cm, alpha) at both refined and unrefined levels. SymbolicAWEModels uses the **unrefined-level coefficients** to map forces back to the structural model.

## Wing types

### REFINE wing type

The `REFINE` wing type creates the most direct coupling between structure and aerodynamics.

#### Structural model

- Wing structure consists of [`WING`](@ref DynamicsType)-type points organized in leading edge (LE) and trailing edge (TE) pairs
- Each consecutive pair (point i, point i+1) forms a structural segment (strut)
- Points can move independently - the wing can deform structurally
- Number of structural segments = (number of WING points) / 2

#### VSM mapping

The structural segments map 1:1 to VSM unrefined sections:

```
Structural points:  [LE₁, TE₁]  [LE₂, TE₂]  [LE₃, TE₃]
                        ↓            ↓            ↓
Unrefined sections:   Sec₁         Sec₂         Sec₃
```

Each VSM section is defined by its LE and TE point positions, taken directly from the structural point positions. VSM can subdivide these into refined panels for higher fidelity; `refined_panel_mapping` maps each refined panel back to its parent section.

#### Force distribution

Forces are computed and distributed in several steps:

1. **VSM solve**: Computes aerodynamic coefficients and per-panel forces/moments in body frame
2. **Panel → section mapping**: Each refined panel maps to its parent unrefined section via `refined_panel_mapping`
3. **Moment-preserving LE/TE split** (`compute_aerostruc_loads`): Each panel's force and moment are split into LE and TE contributions that preserve the total moment about a reference point
4. **Accumulation at structural points**: LE/TE forces are accumulated at the corresponding structural points via the `point_to_vsm_point` mapping

#### Geometry update

Each timestep, structural point positions update VSM section geometry:

```julia
# For each structural point mapped to a VSM section point:
pos_b = R_b_w' * (point.pos_w - wing.origin)
section.LE_point = pos_b  # or section.TE_point
```

This bidirectional coupling allows structural deformation to affect aerodynamics.

### QUATERNION wing type

The `QUATERNION` wing type uses a rigid body representation with optional deformable groups.

#### Structural model

- Wing treated as a rigid body with quaternion-based orientation
- No per-point wing structure - aerodynamic forces applied to wing center of mass
- Optional **groups** represent deformable sections with twist degrees of freedom
- Groups control segment twist angles but don't affect primary aerodynamic geometry

#### VSM mapping

VSM still uses unrefined sections, but they don't correspond to individual structural points:

```
Unrefined sections: [Sec₁,  Sec₂,  Sec₃,  Sec₄,  Sec₅]
                        ╲       ╱       ╲       ╱
Groups:              [─ Group₁ ─]    [─ Group₂ ─]
                     twist DOF θ₁    twist DOF θ₂
```

Multiple unrefined sections can be combined into a single group for twist control. The mapping is configured via:

```julia
group.unrefined_section_idxs = [start_idx:end_idx]
```

#### Force computation

VSM linearization returns integrated aerodynamic coefficients per wing:

1. **VSM linearize**: Computes baseline coefficients and Jacobian around the current operating point
   - Output state `vsm_x = [C_F(3), C_M(3), section_moments(n_unrefined)]`
   - Input state `vsm_y = [va_b(3), twist_angles(n_unrefined), ω_b(3)]`
2. **Symbolic linearization**: The ODE uses `F = q∞ · A · (C₀ + J · Δstate)` where `Δstate = state - state₀`
   - `C₀[1:3]` → total force coefficient, `C₀[4:6]` → total moment coefficient
   - `C₀[7:end]` → per-section twist moment coefficients, summed per group
3. **Rigid body dynamics**: Integrated force/moment drive quaternion dynamics

Groups affect twist deformation applied to VSM sections before the solve. Each group's aerodynamic moment is the sum of its unrefined section moments, driving the twist DOF.

## Refined panel mapping

Both wing types use `refined_panel_mapping` to handle VSM mesh refinement:

### Purpose

VSM can subdivide unrefined sections into multiple refined panels for higher aerodynamic fidelity. The mapping tracks which parent unrefined section each refined panel belongs to.

### Computation

After VSM refinement, `compute_refined_panel_mapping!` finds the closest unrefined section for each refined panel by comparing center positions:

```julia
for each refined_panel in wing.refined_sections
    center = compute_center(refined_panel)
    closest_section = argmin(distance(center, unrefined_section_centers))
    refined_panel_mapping[refined_panel_idx] = closest_section
end
```

### Usage

The mapping enables:

1. **Group twist angles**: Applying the correct twist angle from groups to refined panels via their parent section
2. **Force distribution (REFINE)**: Accumulating refined panel forces at the structural points of their parent section
3. **Linearization (QUATERNION)**: Propagating state perturbations through the correct sections

## Key differences summary

| Aspect | REFINE | QUATERNION |
|--------|--------|------------|
| **Structural representation** | Individual WING points | Rigid body with quaternion |
| **Unrefined section count** | = number of structural LE/TE pairs | Independent of structure |
| **Force distribution** | Per-point via moment-preserving LE/TE split | Integrated force/moment coefficients |
| **Deformation coupling** | Direct: point motion → VSM geometry | Indirect: group twists → VSM sections |
| **Computational cost** | Higher (full VSM solve per step) | Lower (linearized) |
| **Fidelity** | Higher (aeroelastic coupling) | Lower (rigid body) |

## Implementation files

- `src/vsm_refine.jl`: REFINE wing force distribution and geometry updates
- `src/system_structure/types.jl`: Component type definitions (Point, Segment, etc.)
- `src/system_structure/wing.jl`: Wing and VSMWing type definitions, group-to-section mapping
- `src/generate_system/vsm_eqs.jl`: Symbolic VSM equation generation
- `src/generate_system/wing_eqs.jl`: Wing dynamics equation generation
- `src/linearize.jl`: VSM linearization and Jacobian updates
- VortexStepMethod.jl `src/wing_geometry.jl`: `refined_panel_mapping` computation
