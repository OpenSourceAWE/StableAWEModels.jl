```@meta
CurrentModule = SymbolicAWEModels
```

# Getting Started

This guide explains how to get started with SymbolicAWEModels based on your use case:

- **Registry Users**: You want to use the package in your own project
- **Cloned Package Users**: You have cloned the repository and want to run examples
- **Developers**: You want to modify the package source code

---

## For Registry Users

This is the recommended approach for most users who want to use SymbolicAWEModels in their projects.

### Installation

1. **Install Julia** using [juliaup](https://github.com/JuliaLang/juliaup):
   ```bash
   curl -fsSL https://install.julialang.org | sh
   juliaup add release
   juliaup default release
   ```

2. **Create a new project**:
   ```bash
   mkdir my_kite_project
   cd my_kite_project
   julia --project=.
   ```

3. **Add SymbolicAWEModels**:
   ```julia
   using Pkg
   pkg"add SymbolicAWEModels"
   ```

4. **Copy examples and data files**:
   ```julia
   using SymbolicAWEModels
   SymbolicAWEModels.init_module(; force=false)
   ```

   This will:
   - Copy all example files to an `examples/` directory
   - Copy configuration files to a `data/` directory
   - Install necessary dependencies (GLMakie, KiteUtils)

### Running Examples

After running `init_module()`, you can run the examples:

```julia
include("examples/menu.jl")  # Interactive menu
```

Or run a specific example:

```julia
include("examples/ram_air_kite.jl")
```

**Important**: The first time you run an example, it will be slow due to compilation and precompilation (this can take several minutes). Run the same example a second time to see the significant speedup - subsequent runs will be much faster as the compiled code is cached.

### Testing

Run the unit tests (can take about 60 minutes):

```julia
pkg"test SymbolicAWEModels"
```

---

## For Cloned Package Users

If you've cloned the repository and want to run examples **without** modifying the source code:

### Setup

1. **Clone the repository**:
   ```bash
   git clone https://github.com/OpenSourceAWE/SymbolicAWEModels.jl
   cd SymbolicAWEModels.jl
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

   This tells Julia to use the local source code in the current directory instead of downloading the registered package.

### Running Examples

Now you can run any example file (paths are relative to your current directory, **not** the examples directory):

```julia
include("examples/ram_air_kite.jl")
include("examples/simple_tuned_model.jl")
include("examples/menu.jl")
```

**Note**: `--project=examples` tells Julia which project environment to use, but doesn't change your working directory.

---

## For Developers

If you want to contribute to the package or modify its source code:

### Initial Setup

1. **Fork and clone** (see the [Developer Guide](developers.md) for detailed instructions)
   ```bash
   git clone https://github.com/<YourUsername>/SymbolicAWEModels.jl
   cd SymbolicAWEModels.jl
   ```

2. **Install Revise.jl** globally (highly recommended):
   ```julia
   # In a Julia REPL (without --project flag)
   using Pkg
   pkg"add Revise"
   ```

   Configure it to [load automatically](https://timholy.github.io/Revise.jl/stable/config/#Using-Revise-by-default) in your `~/.julia/config/startup.jl`.

### Running Examples During Development

1. **Start Julia with the examples project**:
   ```bash
   julia --project=examples
   ```

2. **Link your local development version** (first time only):
   ```julia
   using Pkg
   pkg"dev ."
   ```

3. **Load Revise** (if not in startup.jl):
   ```julia
   using Revise
   ```

4. **Run examples**:
   ```julia
   include("examples/ram_air_kite.jl")
   ```

   With Revise.jl loaded, any changes you make to the source code in `src/` will be automatically reflected when you run the examples—no need to restart Julia!

**Important**: `--project=examples` sets the project environment but doesn't change your working directory. You still need to include the `examples/` prefix in your paths.

### Disabling Precompilation

When actively developing, you can disable precompilation to speed up Julia startup:

```bash
cp LocalPreferences.toml.default LocalPreferences.toml
```

Remember to delete `LocalPreferences.toml` if you modify the precompilation workload.

### Building Documentation Locally

To preview documentation changes:

1. **Start Julia with the docs project**:
   ```bash
   julia --project=docs
   ```

2. **Link your local development version** (first time only):
   ```julia
   using Pkg
   pkg"dev ."
   ```

3. **Serve the docs with live preview**:
   ```julia
   using LiveServer
   servedocs(launch_browser=true)
   ```

   This will build the documentation and open it in your browser. The docs will automatically rebuild when you save changes to any documentation files.

**Alternative**: Build docs once without live server:
```julia
include("docs/make.jl")
```
Then open `docs/build/index.html` in your browser.

---

## Quick Reference

| Task | Command |
|------|---------|
| Install from registry | `pkg"add SymbolicAWEModels"` |
| Copy examples (registry users) | `SymbolicAWEModels.init_module()` |
| Run examples (cloned/dev) | `julia --project=examples` then `pkg"dev ."` |
| Build docs locally | `julia --project=docs` then `pkg"dev ."` and `servedocs()` |
| Run tests | `pkg"test SymbolicAWEModels"` |

---

## Next Steps

- See the [Examples](examples.md) page for detailed example walkthroughs
- Read the [Developer Guide](developers.md) for contribution guidelines
- Check out the [API Documentation](exported_functions.md) for available functions
