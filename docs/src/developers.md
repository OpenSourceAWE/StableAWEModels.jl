```@meta
CurrentModule = SymbolicAWEModels
```

# Developer Guide

This guide provides instructions and best practices for developers contributing to `SymbolicAWEModels.jl`.

-----

## Prerequisites

Before you begin, ensure you have the following software installed:

  - **Julia**: Version 1.10 or 1.11.
  - **Git**: For version control.
  - **Bash**: A Unix-like shell environment.
  - **Visual Studio Code**: Recommended for its excellent Julia support via the [julia-vscode.org](https://www.julia-vscode.org/) extension.

For detailed setup instructions, see the [Julia and VSCode installation guide](https://www.google.com/search?q=https://OpenSourceAWE.github.io/2024/08/09/installing-julia-with-juliaup.html).

-----

## Getting Started: Development Workflow

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

## Contributing Code: Branches and Pull Requests

To contribute your changes, please follow this standard Git workflow:

**Sync with the Main Project**
Before starting new work, fetch the latest changes from the `upstream` repository and update your local `main` branch. This helps prevent merge conflicts.
```bash
git fetch upstream
git checkout main
git merge upstream/main
```

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

## Improving the Development Experience

### Use Revise.jl for Faster Workflow

We strongly recommend adding **[Revise.jl](https://timholy.github.io/Revise.jl/stable/)** to your global Julia environment. It allows you to modify source code without restarting your Julia session, which is essential for efficient development.

Consider [configuring Revise to run by default](https://timholy.github.io/Revise.jl/stable/config/#Using-Revise-by-default) in your `startup.jl` file.

### Disable Precompilation for Core Development

When actively developing, you can temporarily disable the heavy precompilation workload to speed up restarts. To do this, copy the provided default preferences file:

```bash
cp LocalPreferences.toml.default LocalPreferences.toml
```

**Note:** Remember to delete `LocalPreferences.toml` if you make changes to the precompilation workload itself.

-----

## Coding Style Guidelines

Please adhere to the following style guidelines to maintain code quality and readability:

  - **Environment:** Add packages like `TestEnv` and `Revise` to your global Julia environment, not to the project's `Project.toml`.
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
  - **Settings:** Use the `se()` function from `KiteUtils` to load the settings for the active project. You can specify a file with `set = se("my_settings.yaml")`.

-----

## Outlook

Current development goals include:

  - Adding a Leading Edge Inflatable (LEI) kite model.
  - Implementing a swinging arm system for ground-based testing.
  - Adding a rigid wing model as an alternative to the ram-air kite.

