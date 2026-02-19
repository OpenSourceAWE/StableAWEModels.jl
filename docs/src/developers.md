```@meta
CurrentModule = SymbolicAWEModels
```

# Developer guide

This guide provides instructions and best practices for developers contributing to `SymbolicAWEModels.jl`.

-----

## Prerequisites

Before you begin, ensure you have the following software installed:

  - **Julia**: Latest release version. Install using [juliaup](https://github.com/JuliaLang/juliaup):
    ```bash
    curl -fsSL https://install.julialang.org | sh
    juliaup add release
    juliaup default release
    ```
  - **Git**: For version control.
  - **Bash**: A Unix-like shell environment.
  - **Code editor**: Your preferred code editor with Julia support.

-----

## Getting started: development workflow

Follow these steps to set up your local development environment:

**Fork the Repository**  
[Fork](https://github.com/OpenSourceAWE/SymbolicAWEModels.jl/fork) the `SymbolicAWEModels.jl` repository on GitHub to create your own copy.

**Clone Your Fork**  
Clone your forked repository to your local machine. Replace `<UserName>` with your GitHub username.

```bash
git clone https://github.com/<UserName>/SymbolicAWEModels.jl
```

**Configure the Upstream Remote**
Add the original `OpenSourceAWE` repository as a remote named `upstream`. This allows you to pull in the latest changes from the main project.

```bash
cd SymbolicAWEModels.jl
git remote add upstream https://github.com/OpenSourceAWE/SymbolicAWEModels.jl
```

**Activate the Project**
Start a Julia session with the project environment activated:

```bash
julia --project=.
```

-----

## Contributing code: branches and pull requests

To contribute your changes, please follow this standard Git workflow:

**Sync with the Main Project**
Before starting new work, fetch the latest changes from the `upstream` repository and update your local `main` branch. This helps prevent merge conflicts.
```bash
git fetch upstream
git checkout main
git merge upstream/main
```

**Keep Your Feature Branch Up to Date**
While working on your feature branch, regularly merge the latest changes from `main` to avoid merge conflicts later:
```bash
git fetch upstream
git checkout main
git merge upstream/main
git checkout add_lei_model
git merge main
```

This is especially important for long-running feature branches. Merging frequently makes conflicts smaller and easier to resolve.

**Create a Feature Branch**
Create a new branch from your up-to-date `main` branch. Give it a short, descriptive name that summarizes your change.

```bash
# Create and switch to your new branch
git checkout -b add_lei_model
```

Good branch names include `add_lei_model`, `improve_plot_recipe`, or `fix_winch_dynamics`.

**Make and Commit Your Changes**
Work on your feature and commit your changes as you go. Write clear and concise commit messages.

```bash
git add .
git commit -m "Add initial structure for LEI kite model"
```

**Push to Your Fork**
Push your new branch to your forked repository on GitHub.

```bash
git push -u origin add_lei_model
```

**Create a Pull Request**
Go to the GitHub page for your fork. You should see a prompt to create a pull request from your new branch. Create a pull request that targets the `main` branch of the original `OpenSourceAWE/SymbolicAWEModels.jl` repository. Provide a clear title and a detailed description of your changes.

-----

## Improving the development experience

### Use Revise.jl for faster workflow

We strongly recommend adding **[Revise.jl](https://timholy.github.io/Revise.jl/stable/)** to your global Julia environment. It allows you to modify source code without restarting your Julia session, which is essential for efficient development.

**Install Revise.jl globally:**

```julia
# Start Julia without a project
julia

# In the REPL
using Pkg
pkg"add Revise"
```

**Configure Revise to auto-load on startup:**

Create or edit `~/.julia/config/startup.jl` (on Linux/Mac) or `%USERPROFILE%\.julia\config\startup.jl` (on Windows):

```julia
try
    @eval using Revise
catch e
    @warn "Error initializing Revise" exception=(e, catch_backtrace())
end
```

This will automatically load Revise every time you start Julia. The `try`/`catch` block ensures Julia will still start even if Revise encounters an issue.

**Verify it works:**

Start a new Julia session and you should see Revise load automatically. You can verify by checking:

```julia
julia> @which Revise
```

Now any changes you make to package source code will be automatically reflected in your Julia session!

### Running examples during development

When developing the package, you'll want to test your changes with the examples. Here's how to set up the examples to use your local development version:

#### Setup

1. **From the package root directory**, start Julia with the examples project:
   ```bash
   julia --project=examples
   ```

2. **Link your local development version**:
   ```julia
   ]  # Press ] to enter Pkg mode - prompt shows (examples) pkg>
   dev .
   ```

   This command tells Julia to use the local source code in the current directory (`.`) instead of the registered package version. Use `]st` to verify the package is linked to your local path.

#### Running examples

Now any changes you make to the source code will be immediately reflected when you run the examples (thanks to Revise.jl):

```julia
include("examples/ram_air_kite.jl")
include("examples/simple_tuned_model.jl")
include("examples/menu.jl")
```

**Important**: `--project=examples` sets which project environment to use, but doesn't change your current working directory. You still need to use `examples/` in the include paths.

The `examples/Project.toml` file already contains the necessary dependencies:
- `GLMakie` - for visualization
- `KiteUtils` - for utility functions
- `SymbolicAWEModels` - the package itself

#### Managing package dependencies

**Understanding the Package Manager:**

Press `]` in the Julia REPL to enter package manager (Pkg) mode. The prompt changes to show your current project:
```julia
julia> ]  # Press ] to enter Pkg mode
(examples) pkg>  # Prompt shows you're in the examples project
```

Press backspace to exit Pkg mode and return to the Julia REPL.

**Common Pkg commands:**
- `add PackageName` - Add a package to the current project
- `rm PackageName` - Remove a package
- `dev .` or `dev ..` - Use local source code instead of registered version
- `st` - Show status (list all packages and their versions)
- `up` - Update all packages
- `instantiate` - Install all packages from Project.toml

**Adding packages to the examples:**
```julia
# Start Julia with examples project
julia --project=examples

]  # Enter Pkg mode - prompt shows (examples) pkg>
add YourPackage
st  # Verify the package was added
```

**Adding packages to SymbolicAWEModels itself:**
```julia
# Start Julia with the main project
julia --project=.

]  # Enter Pkg mode - prompt shows (SymbolicAWEModels) pkg>
add YourPackage
st  # Verify the package was added
```

The prompt `(ProjectName) pkg>` always tells you which project you're modifying.

**Tip**: Create a shell alias to quickly start the development environment:
```bash
alias jl-ex='julia --project=examples'
```

### Building documentation locally

To preview documentation changes as you work:

#### Using LiveServer (recommended)

1. **Start Julia with the docs project**:
   ```bash
   julia --project=docs
   ```

2. **Link your local development version** (first time only):
   ```julia
   ]  # Press ] to enter Pkg mode - prompt shows (docs) pkg>
   dev .
   ```

3. **Serve the docs with live reload**:
   ```julia
   using LiveServer
   servedocs(launch_browser=true)
   ```

   This will:
   - Build the documentation
   - Open it in your default browser
   - Watch for changes to documentation files
   - Automatically rebuild and refresh when you save changes

#### Manual build

Alternatively, you can build the documentation once without the live server:

```bash
julia --project=docs
```

```julia
include("docs/make.jl")
```

Then open `docs/build/index.html` in your browser.

**Note**: If you make changes to the package source code (not just documentation), you'll need to reload Julia or use Revise.jl for the changes to be reflected in the built documentation.

-----

## Coding style guidelines

Please adhere to the following style guidelines to maintain code quality and readability:

  - **Environment:** Add packages like `Revise` to your global Julia environment, not to the project's `Project.toml`.
  - **No Magic Numbers:** Avoid hard-coded values (e.g., `9.81`). Define them as constants (e.g., `G_EARTH`) or read them from a configuration file.
  - **Line Length:** Keep lines under 100 characters, including in documentation.
  - **Operators:**
      - Use the tilde `~` for scalar equations in `ModelingToolkit` instead of the broadcasted `.~`.
      - Use the `\cdot` operator for the dot product (`⋅`) for improved readability.
      - Enclose binary operators (`+`, `*`, `=`) with single spaces (e.g., `y = a * x + b`).
  - **Spacing:** Use a space after a comma (e.g., `my_function(x, y)`).
  - **Alignment:** Align assignment operators (`=`) in blocks of related assignments to improve readability:
    ```julia
    tether_rhs = [force_eqs[j, i].rhs for j in 1:3]
    kite_rhs   = [force_eqs[j, i+3].rhs for j in 1:3]
    f_xy       = dot(tether_rhs, e_z) * e_z
    ```
  - **Settings:** Use the `Settings()` constructor to load the settings for the active project. You can specify a file with `set = Settings("my_settings.yaml")`. Use `set = Settings("")` to load the default settings file.

-----

## Source code organization

The source code is organized into modular directories:

- **`src/system_structure/`** — component types and assembly
  - `types.jl`: [`Point`](@ref), [`Segment`](@ref), [`Pulley`](@ref), [`Tether`](@ref), [`Winch`](@ref), [`Group`](@ref), [`Transform`](@ref)
  - `wing.jl`: [`AbstractWing`](@ref) and [`VSMWing`](@ref) types, aerodynamic setup
  - `system_structure_core.jl`: [`SystemStructure`](@ref) constructor, reference resolution
  - `named_collection.jl`: Symbol-based indexing ([`NamedCollection`](@ref))
  - `transforms.jl`: Spherical coordinate transforms
  - `utilities.jl`: Validation, tether creation, state management
- **`src/generate_system/`** — symbolic equation generation
  - `create_sys.jl`: Top-level orchestrator
  - `point_eqs.jl`, `segment_eqs.jl`, `wing_eqs.jl`, etc.: per-subsystem equations
- **`src/yaml_loader.jl`** — YAML configuration file parser ([`load_sys_struct_from_yaml`](@ref))
- **`src/linearize.jl`** — VSM linearization ([`linearize!`](@ref))
- **`src/simulate.jl`** — high-level simulation functions ([`sim!`](@ref), [`sim_oscillate!`](@ref))

## Outlook

Current development goals include:

  - Adding a Leading Edge Inflatable (LEI) kite model
  - YAML-based model validation packages for ram air and V3 kites

