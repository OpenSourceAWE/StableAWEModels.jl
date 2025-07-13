using SymbolicAWEModels, VortexStepMethod, KiteUtils, ControlPlots, Statistics

set = Settings("system.yaml")
sam = SymbolicAWEModel(set)
SymbolicAWEModels.init!(sam)
sys = sam.sys
find_steady_state!(sam; t=10.0, dt=3.0)

simple_set = Settings("system.yaml")
simple_set.physical_model = "simple_ram"
ssam = SymbolicAWEModel(simple_set)
SymbolicAWEModels.init!(ssam)
ssys = ssam.sys
find_steady_state!(ssam)

# update_ssam!(ssam, sam)

steps = 100
winches = sam.sys_struct.winches
sam.integrator.ps[sys.fix_wing] = true
tether_lens = zeros(3, steps)
set_values = -sam.set.drum_radius .* sam.integrator[sys.winch_force] .+ 1.0
@time for i in 1:steps
    next_step!(sam; set_values, vsm_interval=0)
    [tether_lens[j,i] = winches[j].tether_len for j in 1:3]
end
plotx(1:steps, tether_lens[1,:], tether_lens[2,:], tether_lens[3,:])

