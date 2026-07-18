# Understand coordinate systems

## Explain
Can you explain commit a2aed4cb

**Summary:** Bart tried switching the no-ref-points fallback orientation of the body frame from principal-inertia frame to CAD orientation, and in parallel made the principal-frame computation itself configurable (Y-axis-constrained vs. generic eigendecomposition). Uwe's commit a2aed4cb reverts only the orientation-fallback change back to the original (principal-inertia frame), presumably because the CAD-orientation default broke something or wasn't the desired behavior, while keeping the unrelated inertia-method refactor from the intervening commits.

**Reference Frames**
1. CAD frame (c) — the fixed reference frame the wing's raw geometry/mesh comes in (points, pos_cad, the inertia tensor I_cad as returned by the aero hook). Everything else is defined relative to this.

2. Principal-inertia frame (p) — the frame in which the inertia tensor is diagonal. Computed from I_cad via either principal_frame (generic eigendecomposition) or calc_inertia_y_rotation (closed-form, Y-axis-constrained, for wings symmetric about the XZ-plane), selected by wing.principal_frame_method. The rotation R_p_to_c (principal → CAD) is stored, and wing.inertia_principal holds the diagonalized moments.

3. Body frame (b) — the frame the rigid-body dynamics is actually integrated in, origin at the COM. It's set either:

from ref points (z_ref_points/y_ref_points) via calc_particle_dynamics_wing_frame, when those are provided — letting you pin the body axes to specific geometric features, or
as a fallback equal to the principal frame (R_b_to_c .= wing.R_p_to_c) when no ref points are given — this is exactly the line commit a2aed4cb restored.
A fourth composed rotation, R_b_to_p = R_p_to_c' * R_b_to_c (body → principal), is also stored — that's what the actual equations of motion use, since the inertia tensor is diagonal in the principal frame.

So the naming convention throughout is R_x_to_y: rotates a vector expressed in frame x into frame y. The commit you asked about earlier was purely about which frame the body frame defaults to (CAD vs. principal) when ref points aren't supplied — it doesn't touch the CAD or principal frame definitions themselves.

## TODO
Use the version field of Settings to distinguish:

- Version one: Body frame is Principle frame
- Version two: Body frame is CAD frame

Version two makes sense if the direction of the x-axis of the CAD frame is 
from the leading edge point at the center to the trailing edge point of the center. 