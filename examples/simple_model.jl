using SymbolicAWEModels, VortexStepMethod, KiteUtils, ControlPlots, Statistics

set = Settings("system.yaml")
sam = SymbolicAWEModel(set)
SymbolicAWEModels.init!(sam)
find_steady_state!(sam)

simple_set = Settings("system.yaml")
simple_set.physical_model = "simple_ram"
ssam = SymbolicAWEModel(simple_set)
SymbolicAWEModels.init!(ssam)
ssys = ssam.sys
find_steady_state!(ssam)

# update_ssam!(ssam, sam)

steps = 50
winches = ssam.sys_struct.winches
ssam.integrator.ps[ssys.fix_wing] = true
tether_lens = zeros(3, steps)
tether_torques = -sam.set.drum_radius .* sam.integrator[ssys.tether_spring_force]
set_values = [mean(tether_torques[1:2]), tether_torques[3], tether_torques[4]]
for i in 1:steps
    @show set_values
    next_step!(ssam; set_values)
    [tether_lens[j,i] = winches[j].tether_len for j in 1:3]
end
plot(1:steps, [tether_lens[i,:] for i in 1:3])

