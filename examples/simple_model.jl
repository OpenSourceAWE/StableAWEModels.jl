using SymbolicAWEModels, VortexStepMethod, KiteUtils

set = Settings("system.yaml")
sam = SymbolicAWEModel(set)
SymbolicAWEModels.init!(sam)

