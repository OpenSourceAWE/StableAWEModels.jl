using Pkg
Pkg.add(url="https://github.com/OpenSourceAWE/SymbolicAWEModels.jl.git", rev="internal_error")
using SymbolicAWEModels, KiteUtils, VortexStepMethod, ModelingToolkit

set_data_path("data")
if !ispath("data/settings.yaml")
    SymbolicAWEModels.copy_model_settings()
end

function generate_sys(segments)
    set = Settings("system.yaml")
    set.segments = segments
    set.physical_model = "ram"
    sam = SymbolicAWEModel(set)
    println("Creating sys with $segments segments")
    @time inputs = SymbolicAWEModels.create_sys!(sam, sam.sys_struct; prn=false)
    println("Simplifying sys with $segments segments")
    @time sys = mtkcompile(sam.full_sys; inputs)
end

generate_sys(3)
nothing
