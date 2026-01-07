# Yaw-rate vs Turn-rate in SymbolicAWEModels

## Definitions (from `src/generate_system.jl`)
- **Heading ψ**: Wind-perpendicular/view-frame heading of the body x-axis. Built from `R_v_w = [x y z]` with `z = normalize(wing_pos)`, `y = normalize(z × e_x)`, `x = y × z`; then `ψ = atan(heading_x, heading_z)` where `heading_x = e_x_perp ⋅ wind_cross_z`, `heading_z = e_x_perp[3]`, `e_x_perp` is `-e_x` projected onto the plane ⟂ wind.
- **Turn-rate vector ω_v**: Body angular velocity expressed in the view frame: `ω_v = R_v_w' * (R_b_w * ω_b)`. Logged as `sl.turn_rates`; the plotted component is `ω_v,z = sl.turn_rates[3]` (spin about tether/radial axis).

## What the plots show (`ext/SymbolicAWEModelsMakieExt.jl`)
- **Yaw-rate panel (`plot_yaw_rate`)**: Unwrap `sl.heading` to remove ±π jumps, then finite-difference in time: `ψ̇ = diff(rad2deg(ψ_unwrapped)) ./ diff(t)` (deg/s). This matches the paper’s yaw-rate in the turn-rate law `ψ̇ = g_k v_a (u_s(t − d(t)) − u_s,0)`.
- **Turn-rate panel (`plot_turn_rates`)**: Also computes `ψ̇` as above (for comparison) and plots `ω_v,z` converted to deg/s from `sl.turn_rates[3]`.

## When they match vs differ
- **Coincide** when the view/tether frame is effectively stationary about its z-axis and the body x-axis stays well in the wind-perpendicular plane (e.g., smooth circular flight with fixed tether length). Then `ψ̇ ≈ ω_v,z`.
- **Differ** when tether direction changes rapidly, view frame spins about z, or projection effects from pitch/roll/sideslip matter; then `ψ̇` (angle rate) and `ω_v,z` (gyro spin about tether axis) separate.
