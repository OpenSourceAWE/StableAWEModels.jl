# VSM Coupling

This document explains how SymbolicAWEModels couples with the Vortex Step Method (VSM) for aerodynamic force computation. The coupling differs between the two wing types: `QUATERNION` and `REFINE`.

## Overview

Both wing types use the concept of **unrefined segments** (or sections) as the fundamental structural element that maps to VSM geometry:

- **Unrefined segment/section**: A structural element defined by two points (leading edge and trailing edge) along the wing span
- **Refined panel**: VSM can subdivide unrefined segments into multiple panels for higher aerodynamic fidelity
- **`refined_panel_mapping`**: Maps each refined panel back to its parent unrefined segment

The VSM solver computes aerodynamic coefficients (cl, cd, cm, alpha) at both refined and unrefined levels. SymbolicAWEModels uses the **unrefined-level coefficients** to map forces back to the structural model.

## Wing Types

### REFINE Wing Type

The `REFINE` wing type creates the most direct coupling between structure and aerodynamics.

#### Structural Model

- Wing structure consists of `WING`-type points organized in leading edge (LE) and trailing edge (TE) pairs
- Each consecutive pair (point i, point i+1) forms a structural segment (strut)
- Points can move independently - the wing can deform structurally
- Number of structural segments = (number of WING points) / 2

#### VSM Mapping

The structural segments map 1:1 to VSM unrefined sections:

```
Structural points:  [LE₁, TE₁, LE₂, TE₂, LE₃, TE₃, ...]
                      ↓    ↓    ↓    ↓    ↓    ↓
VSM sections:       [Sec₁,   Sec₂,   Sec₃,   ...]
VSM panels:            [Panel₁,  Panel₂,  ...]
```

Each VSM section is defined by its LE and TE point positions, taken directly from the structural point positions.

#### Force Distribution

Forces are computed and distributed in several steps:

1. **VSM solve**: Computes aerodynamic coefficients for each unrefined segment
   - `cl_unrefined_array`, `cd_unrefined_array`, `cm_unrefined_array`, `alpha_unrefined_array`

2. **Panel-level force calculation**: For each unrefined panel (connects sections i and i+1):
   - Compute total lift and drag from cl, cd
   - Use pitching moment (cm) to determine center of pressure
   - Split forces between LE and TE based on moment equilibrium:
     - `x_cp = 0.25 + cm/cl` (normalized chord position)
     - `L_LE = L_total × (1 - x_cp)`, `L_TE = L_total × x_cp`
     - Drag split equally: `D_LE = D_total / 2`, `D_TE = D_total / 2`

3. **Distribution to structural points**: Each panel's forces distributed 50/50 to adjacent sections
   - Panel i connects sections i and i+1
   - Section i receives 50% of panel forces
   - Section i+1 receives 50% of panel forces
   - LE forces go to LE point, TE forces go to TE point

This creates a smooth spanwise load distribution where internal sections receive contributions from both adjacent panels.

#### Geometry Update

Each timestep, structural point positions update VSM section geometry:

```julia
# For each structural point mapped to a VSM section point:
pos_b = R_b_w' * (point.pos_w - wing.origin)
section.LE_point = pos_b  # or section.TE_point
```

This bidirectional coupling allows structural deformation to affect aerodynamics.

### QUATERNION Wing Type

The `QUATERNION` wing type uses a rigid body representation with optional deformable groups.

#### Structural Model

- Wing treated as a rigid body with quaternion-based orientation
- No per-point wing structure - aerodynamic forces applied to wing center of mass
- Optional **groups** represent deformable sections with twist degrees of freedom
- Groups control segment twist angles but don't affect primary aerodynamic geometry

#### VSM Mapping

VSM still uses unrefined sections, but they don't correspond to individual structural points:

```
VSM sections:       [Sec₁,   Sec₂,   Sec₃,   Sec₄,   Sec₅,   ...]
Groups:                [────Group₁────]  [────Group₂────]
Unrefined sections:      [───Panel₁───]  [───Panel₂───]  [...]
```

Multiple unrefined sections can be combined into a single group for twist control. The mapping is configured via:

```julia
group.unrefined_section_idxs = [start_idx:end_idx]
```

#### Force Integration

Since the wing is rigid:

1. **VSM solve**: Computes forces per panel using unrefined coefficients
2. **Force integration**: All panel forces integrated into:
   - Total force vector at wing center of mass
   - Total moment vector about wing center of mass
3. **Rigid body dynamics**: Integrated force/moment drive quaternion dynamics

Groups only affect twist deformation applied to VSM sections before the solve, not force distribution.

## Refined Panel Mapping

Both wing types use `refined_panel_mapping` to handle VSM mesh refinement:

### Purpose

VSM can subdivide panels for aerodynamic accuracy (e.g., 10 structural segments → 50 refined panels). The mapping tracks which unrefined segment each refined panel came from.

### Computation

After VSM refinement, `compute_refined_panel_mapping!` finds the closest unrefined panel for each refined panel by comparing panel center positions:

```julia
for each refined_panel in wing.refined_sections
    center = compute_center(refined_panel)
    closest_unrefined = argmin(distance(center, unrefined_centers))
    refined_panel_mapping[refined_panel_idx] = closest_unrefined
end
```

### Usage

The mapping enables:

1. **Group twist angles**: Applying the correct twist angle from unrefined groups to refined sections
2. **Deformation**: Mapping twist/deflection from structural groups to aerodynamic panels
3. **Linearization**: Propagating state perturbations to the correct panels

For REFINE wings, since structural segments map 1:1 to unrefined sections, the mapping is typically trivial (identity mapping) unless VSM refinement is used.

## Key Differences Summary

| Aspect | REFINE | QUATERNION |
|--------|--------|------------|
| **Structural representation** | Individual WING points | Rigid body with quaternion |
| **Unrefined section count** | = number of structural segments | Independent of structure |
| **Force distribution** | Per-point via LE/TE split | Integrated force/moment |
| **Deformation coupling** | Direct: point motion → VSM geometry | Indirect: group twists → VSM sections |
| **Computational cost** | 2-3× higher | Lower |
| **Fidelity** | Higher (aeroelastic coupling) | Lower (rigid body) |

## Implementation Files

- `src/vsm_refine.jl`: REFINE wing force distribution and geometry updates
- `src/system_structure.jl`: Wing type definitions and group-to-section mapping (lines 1399-1420)
- `src/generate_system.jl`: Symbolic equation generation for both wing types
- VortexStepMethod.jl `src/wing_geometry.jl`: `refined_panel_mapping` computation (lines 703-756)
