# Copyright (c) 2022, 2024 Uwe Fechner
# SPDX-License-Identifier: MIT
using REPL.TerminalMenus

options = [
        "catenary_line = include(\"catenary_line.jl\")",
        "hanging_mass = include(\"hanging_mass.jl\")",
        "simple_pulley = include(\"simple_pulley.jl\")",
        "lin_ram_model = include(\"lin_ram_model.jl\")",
        "lin_simple_tuned_model = include(\"lin_simple_tuned_model.jl\")",
        "pulley = include(\"pulley.jl\")",
        "ram_air_kite = include(\"ram_air_kite.jl\")",
        "sam_tutorial = include(\"sam_tutorial.jl\")",
        "simple_lin_model = include(\"simple_lin_model.jl\")",
        "simple_tuned_model = include(\"simple_tuned_model.jl\")",
        "tether_props = include(\"tether_props.jl\")",
        "quit"
]

function example_menu()
    active = true
    while active
        menu = RadioMenu(options, pagesize=8)
        choice = request("\nChoose function to execute or `q` to quit: ", menu)

        if choice != -1 && choice != length(options)
            eval(Meta.parse(options[choice]))
        else
            println("Left menu. Press <ctrl><d> to quit Julia!")
            active = false
        end
    end
end

example_menu()
