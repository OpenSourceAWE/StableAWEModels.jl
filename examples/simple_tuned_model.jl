# SPDX-FileCopyrightText: 2025 Bart van de Lint <bart@vandelint.net>
#
# SPDX-License-Identifier: MPL-2.0

using Pkg
if ! ("ControlPlots" ∈ keys(Pkg.project().dependencies))
    using TestEnv; TestEnv.activate()
end

using SymbolicAWEModels, ControlPlots

set = Settings("system.yaml")
sam = SymbolicAWEModel(set, "ram")
init!(sam)

tether_set = Settings("system.yaml")
tether_sam = SymbolicAWEModel(tether_set, "tether")
init!(tether_sam)

simple_set = Settings("system.yaml")
simple_sam = SymbolicAWEModel(simple_set, "simple_ram")
init!(simple_sam)

bias = 0.3
find_steady_state!(sam; dt=3, t=10)
SymbolicAWEModels.copy_to_simple!(sam, tether_sam, simple_sam)
@show sam.sys_struct.wings[1].tether_moment
@show simple_sam.sys_struct.wings[1].tether_moment

sl, _ = sim_oscillate!(sam; total_time=5.0, prn=true, bias) # TODO: add first frac ram model
display(plot(sam.sys_struct, sl; plot_default=false, plot_elevation=true,
             plot_aoa=true, plot_heading=true, plot_aero_moment=true))

simple_sl, _ = sim_oscillate!(simple_sam; total_time=5.0, prn=true, bias)
display(plot(simple_sam.sys_struct, simple_sl; plot_default=false, plot_elevation=true,
             plot_aoa=true, plot_heading=true, plot_aero_moment=true))

# 2.5deg difference in aoa and elevation
# should be solved by changing the attach point, twist is good

