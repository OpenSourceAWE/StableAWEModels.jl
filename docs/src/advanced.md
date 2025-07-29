```@meta
CurrentModule = SymbolicAWEModels
```
# Advanced usage
For advanced users it is suggested to install git, bash and vscode or vscodium in addition to 
Julia. vscode and vscodium both have a very good plugin for Julia support, see [https://www.julia-vscode.org](https://www.julia-vscode.org/).
Installation instructions: [Julia and VSCode](https://OpenSourceAWE.github.io/2024/08/09/installing-julia-with-juliaup.html) .

## Forking the repository and creating a custom system image
1. Add [Revise](https://timholy.github.io/Revise.jl/stable/config/#Using-Revise-by-default)
    to avoid having to restart your Julia session. This makes development a lot easier and 
    faster.
2. Go to the website https://github.com/OpenSourceAWE/SymbolicAWEModels.jl and click on the 
    **Fork** button at the top right.
3. clone the new repository which is owned by you with a command similar to this one: 
    ```git clone https://github.com/OpenSourceAWE/SymbolicAWEModels.jl``` 
    Your own git user name must appear in the URL, otherwise you will not be able to push 
    your changes.

After cloning the repo you can launch julia with the command:
```bash
julia --project=.
```
And run an example:
```julia
julia> @time include("examples/ram_air_kite.jl")
```

Running the example after starting your Julia session and assuming everything is precompiled
takes around 30 seconds, while it takes less than 10 seconds the second time. This is why 
Revise is so useful: you can edit the source code without having to restart Julia, and save 
20 seconds (or a lot more if precompilation is triggered).

## Hints for Developers
### Coding style

- add the packages `TestEnv` and `Revise` to your global environment, not to any project

- avoid hard-coded numeric values like `9.81` in the code, instead define a global constant 
    `G_EARTH` or read this value from a configuration file

- stick to a line length limit of 100 characters (also in docs)

- try to avoid dot operators unless you have to. 
Bad: `norm1        .~ norm(segment)`
Good: `norm1        ~ norm(segment)`

- if you need to refer to the settings you can use `se()` which will load the settings of the 
    active project. To define the active project use a line like `set = se("system_3l.yaml")` 
    at the beginning of your program.
- use the `\cdot` operator for the dot product for improved readability
- use a space after a comma, e.g. `force_eqs[j, i]`
- enclose operators like `+` and `*` in single spaces, like `0.5 * (s.pos[s.i_C] + s.pos[s.i_D])`;  
  exception: `mass_tether_particle[i-1]`
- try to align the equation signs for improved readability like this:
```julia
    tether_rhs        = [force_eqs[j, i].rhs for j in 1:3]
    kite_rhs          = [force_eqs[j, i+3].rhs for j in 1:3]
    f_xy              = dot(tether_rhs, e_z) * e_z
```

## Outlook

The next steps:
- add LEI kite
- add swinging arm system
- add a rigid wing model

