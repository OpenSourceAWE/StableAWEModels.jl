using SymbolicAWEModels, VortexStepMethod, KiteUtils, Plots

set = Settings("system.yaml")
sam = SymbolicAWEModel(set)
SymbolicAWEModels.init!(sam)
find_steady_state!(sam)

simple_set = Settings("system.yaml")
simple_set.physical_model = "simple_ram"
ssam = SymbolicAWEModel(simple_set)
SymbolicAWEModels.init!(ssam)
find_steady_state!(ssam)

# update_ssam!(ssam, sam)
ssam.integrator[ssam.sys.fix_wing] = true
for t in 0:0.05:1.0
    find_steady_winch_state!(ssam)
    next_step!()

