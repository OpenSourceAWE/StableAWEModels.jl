```@meta
CurrentModule = SymbolicAWEModels
```

# SymbolicAWEModels
Documentation for the package [SymbolicAWEModels](https://github.com/OpenSourceAWE/SymbolicAWEModels.jl).

This package provides modular symbolic models of Airborne Wind Energy (AWE) systems, 
which consist of one or more wings, tethers, winches and a bridle system with 
or without pulleys. The kite is modeled as a deforming rigid body with orientation governed 
by quaternion dynamics. The aerodynamic forces and moments are computed using the 
Vortex Step Method. The tether is modeled as point masses connected by spring-damper 
elements, with aerodynamic drag modeled realistically. The winchs are modeled as
motors/generators that can reel in or out the tethers.

The [`SymbolicAWEModel`](@ref) has the following subcomponents implemented in separate packages:
- AtmosphericModel from [AtmosphericModels](https://github.com/aenarete/AtmosphericModels.jl)
- WinchModel from [WinchModels](https://github.com/aenarete/WinchModels.jl) 
- The aerodynamic forces and moments of some of the models are calculated using the 
  package [VortexStepMethod](https://github.com/Albatross-Kite-Transport/VortexStepMethod.jl)

This package is part of the Julia Kite Power Tools, which consist of the following packages:

![Julia Kite Power Tools](kite_power_tools.png)

## Installation
Install [Julia 1.11](https://OpenSourceAWE.github.io/2024/08/09/installing-julia-with-juliaup.html), if you haven't already. On Linux, make sure that Python3 and Matplotlib are installed:
```
sudo apt install python3-matplotlib
```
Before installing this software it is suggested to create a new project, for example like this:
```bash
mkdir test
cd test
julia --project="."
```
Then add SymbolicAWEModels from  Julia's package manager, by typing:
```julia
using Pkg
pkg"add SymbolicAWEModels"
``` 
at the Julia prompt. You can run the unit tests with the command (careful, can take 60 min):
```julia
pkg"test SymbolicAWEModels"
```
You can copy the examples to your project with:
```julia
using SymbolicAWEModels
SymbolicAWEModels.install_examples()
```
This also adds the extra packages, needed for the examples to the project. Furthermore, it creates a folder `data`
with some example input files. You can now run the examples with the command:
```julia
include("examples/menu.jl")
```
You can also run the ram-air-kite example like this:
```julia
include("examples/ram_air_kite.jl")
```
This might take two minutes. To speed up the model initialization, you can create a system image:
```bash
cd bin
./create_sys_image
```
If you now launch Julia with `./bin/run_julia` and then run the above example again, it should 
be about three times faster.

## Ram air kite model
This model represents the kite as a deforming rigid body, with orientation governed by 
quaternion dynamics. Aerodynamics are computed using the Vortex Step Method. The kite is 
controlled from the ground via four tethers.

Initialize:
```julia
using SymbolicAWEModels, ControlPlots
set = Settings("system.yaml")
sam = SymbolicAWEModel(set, "ram")
init!(sam)
```

Model and plot:
```julia
log = sim_oscillate!(sam)
plot(sam.sys_struct, log)
```

The simple ram air kite model removes the bridle system, and has 1-segment tethers. 
The tether properties and attach points can be approximated using the complex ram air kite 
model and a helper tether model.

Initialize:
```julia
init!(sam)
tether_sam = SymbolicAWEModel(set, "tether")
init!(tether_sam)
simple_sam = SymbolicAWEModel(set, "simple_ram")
init!(simple_sam)
```

Model and plot:
```julia
SymbolicAWEModels.copy_to_simple!(sam, tether_sam, simple_sam)
simple_log = sim_oscillate!(simple_sam)
plot(sam.sys_struct, simple_log)
```

## See also
- [Research Fechner](https://research.tudelft.nl/en/publications/?search=Fechner+wind&pageSize=50&ordering=rating&descending=true) for the scientic background of the winches and tethers.
- More kite models [KiteModels](https://github.com/ufechner7/KiteModels.jl)
- The meta-package [KiteSimulators](https://github.com/aenarete/KiteSimulators.jl)
- the package [KiteUtils](https://github.com/OpenSourceAWE/KiteUtils.jl)
- the packages [WinchModels](https://github.com/aenarete/WinchModels.jl) and [KitePodModels](https://github.com/aenarete/KitePodModels.jl) and [AtmosphericModels](https://github.com/aenarete/AtmosphericModels.jl)
- the packages [KiteControllers](https://github.com/aenarete/KiteControllers.jl) and [KiteViewers](https://github.com/aenarete/KiteViewers.jl)
- the [VortexStepMethod](https://github.com/Albatross-Kite-Transport/VortexStepMethod.jl)

## Questions?
If you have any questions, please ask in the Julia Discourse forum in the section 
[modelling and simulation](https://discourse.julialang.org/c/domain/models) , or in in the 
section [First steps](https://discourse.julialang.org/c/first-steps) . The Julia community 
is very friendly and responsive.
You can also send an email to Bart van de Lint (bart@vandelint.net).

Authors: Bart van de Lint (bart@vandelint.net), Uwe Fechner (uwe.fechner.msc@gmail.com)
