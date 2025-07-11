using SymbolicAWEModels, VortexStepMethod, KiteUtils

set = Settings("system.yaml")
sam = SymbolicAWEModel(set)
SymbolicAWEModels.init!(sam)
find_steady_state!(sam)

simple_set = Settings("system.yaml")
simple_set.physical_model = "simple_ram"
simple_sam = SymbolicAWEModel(simple_set)
SymbolicAWEModels.init!(simple_sam)
find_steady_state!(simple_sam)

# update_simple_sam!(simple_sam, sam)

