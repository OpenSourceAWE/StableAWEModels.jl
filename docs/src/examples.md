```@meta
CurrentModule = SymbolicAWEModels
```
# Examples for using the ram air kite model

## Visualization with GLMakie

SymbolicAWEModels provides plotting functionality through a package extension that automatically loads when you use GLMakie. This means you don't need to explicitly load any plotting package beyond GLMakie itself.

**To enable plotting:**
```julia
using SymbolicAWEModels
using GLMakie  # This automatically loads the plotting extension
```

Once GLMakie is loaded, you can use the `plot` function to visualize:
- **3D System Structure**: Interactive 3D visualization of the kite system with clickable segments
  ```julia
  plot(sam.sys_struct)  # Plot the system structure
  ```

- **Time-Series Data**: Multi-panel plots of simulation results
  ```julia
  (log, _) = sim_oscillate!(sam)
  plot(sam.sys_struct, log; plot_default=true)  # Plot simulation results
  ```

The plotting functions support many keyword arguments to customize which data is displayed. See the function documentation for details:
```julia
?plot  # In Julia REPL, shows available options
```

## Getting Started with Examples

The approach depends on how you're using SymbolicAWEModels:

### Registry Users

If you installed the package via `pkg"add SymbolicAWEModels"`:

1. **Copy examples to your project**:
   ```julia
   using SymbolicAWEModels
   SymbolicAWEModels.init_module()
   ```

   This copies all example files to an `examples/` directory and installs dependencies.

2. **Run examples**:
   ```julia
   include("examples/ram_air_kite.jl")
   include("examples/menu.jl")  # Interactive menu
   ```

### Cloned Package Users

If you cloned the repository but aren't modifying source code:

1. **Navigate to the repository**:
   ```bash
   cd path/to/SymbolicAWEModels.jl
   ```

2. **Start Julia with the examples project**:
   ```bash
   julia --project=examples
   ```

3. **Link the local package** (first time only):
   ```julia
   using Pkg
   pkg"dev ."
   ```

4. **Run examples**:
   ```julia
   include("examples/ram_air_kite.jl")
   include("examples/menu.jl")
   ```

   **Note**: `--project=examples` sets the project environment but doesn't change your working directory, so you still need the `examples/` prefix.

### Developers

If you're developing the package, see the [Developer Guide](developers.md#running-examples-during-development) for detailed instructions on using Revise.jl and running examples with live code reloading.

For more details, see the [Getting Started Guide](getting_started.md).

## Running the first example
```julia
include("examples/ram_air_kite.jl")
```
Expected output for first run:
```julia
[ Info: Loading packages
[ Info: Precompiling SymbolicAWEModelsMakieExt [ba1985c5-ebcb-5009-8fba-190d67777252] (cache misses...)
Time elapsed: 24.149010176 s
[ Info: Creating SymbolicAWEModel:
Time elapsed: 38.51030329 s
[ Info: System initialized at:
Time elapsed: 120.861073757 s
[ Info: Generating oscillating steering commands...
[ Info: Starting simulation
┌ Info: Performance Summary:
│ Component    | Speedup (×)  | Total Time
│ -------------|--------------|------------
│ Simulation   |         2.34 |       2.14
│ Step         |         5.96 |       0.84
│ Integrator   |        34.15 |       0.15
└ VSM          |         2.62 |       1.91
[ Info: Simulated at:
Time elapsed: 213.981574478 s
[ Info: Plotted at:
Time elapsed: 220.05251771 s
220.05251771
```
After the second time it runs much faster, because the simplified ODE system is cached in the 
`model*.bin` file in the `data` folder and the deserialization is precompiled:
```julia
[ Info: Loading packages 
Time elapsed: 0.017465498 s
[ Info: Creating SymbolicAWEModel:
Time elapsed: 1.85496798 s
[ Info: System initialized at:
Time elapsed: 5.364508521 s
[ Info: Generating oscillating steering commands...
[ Info: Starting simulation
┌ Info: Performance Summary:
│ Component    | Speedup (×)  | Total Time
│ -------------|--------------|------------
│ Simulation   |         2.53 |       1.97
│ Step         |         6.61 |       0.76
│ Integrator   |        41.21 |       0.12
└ VSM          |         2.79 |       1.79
[ Info: Simulated at:
Time elapsed: 7.995500858 s
[ Info: Plotted at:
Time elapsed: 8.309523118 s
8.309523118
```

In this example, the kite is first stabilized, and then a sinus-shaped steering input is 
applied such that the kite is dancing in the sky.

![Oscillating steering input response](assets/oscillating_steering.png)

## Running the second example
```julia
include("examples/simple_model.jl")
```

```julia
Tether 1: ω_n=227.039 rad/s,
                      T_s=0.1 s, 
                      ζ=0.1319, c=3.4209 Ns/m
Tether 2: ω_n=228.866 rad/s,
                      T_s=0.1 s, 
                      ζ=0.1309, c=3.4211 Ns/m
Tether 3: ω_n=235.585 rad/s,
                      T_s=0.1 s, 
                      ζ=0.1272, c=0.8549 Ns/m
Tether 4: ω_n=236.834 rad/s,
                      T_s=0.1 s, 
                      ζ=0.1265, c=0.8552 Ns/m
Summary of Results:
Tether 1: k = 2943.0749476475257 N/m, c = 3.4208593453647103 Ns/m
Tether 2: k = 2990.816045085089 N/m, c = 3.4210625600793843 Ns/m
Tether 3: k = 791.9269801315223 N/m, c = 0.8549123148071024 Ns/m
Tether 4: k = 800.5953494095929 N/m, c = 0.8551804418169453 Ns/m
[ Info: Generating oscillating steering commands...
[ Info: Starting simulation
┌ Info: Performance Summary:
│ Component    | Speedup (×)  | Total Time
│ -------------|--------------|------------
│ Simulation   |         1.99 |       1.25
│ Step         |         5.54 |       0.45
│ Integrator   |        24.27 |       0.10
└ VSM          |         2.59 |       0.97
[ Info: Generating oscillating steering commands...
[ Info: Starting simulation
┌ Info: Performance Summary:
│ Component    | Speedup (×)  | Total Time
│ -------------|--------------|------------
│ Simulation   |         0.46 |       5.44
│ Step         |        10.22 |       0.24
│ Integrator   |      1385.79 |       0.00
└ VSM          |         3.49 |       0.72
```

The simple model connects the tethers directly to the wing. In contrast, the default model 
features a more complex [speed system](https://kiteboarding.com/proddetail.asp?prod=ozone-r1v4-pro-tune-speedsystem-complete)
that uses pulleys and multiple attachment points. By 
matching key dynamic properties—such as the moment on the wing groups, stiffness, and 
damping—the simple model can be tuned to closely replicate the response of the more complex 
default model, at a much better performance.

![Oscillating steering input response, simple system](assets/oscillating_steering_simple.png)

## Linearization
The following example creates a nonlinear system model, finds a steady-state operating point, 
linearizes the model around this operating point and creates bode plots from inputs (torques)
to output (heading).
```julia
include("examples/lin_ram_model.jl")
```
See: [`lin_ram_model.jl`](https://github.com/OpenSourceAWE/SymbolicAWEModels.jl/blob/main/examples/lin_ram_model.jl)

