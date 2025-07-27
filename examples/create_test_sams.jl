# SPDX-FileCopyrightText: 2025 Bart van de Lint <bart@vandelint.net>
#
# SPDX-License-Identifier: MPL-2.0

using SymbolicAWEModels

set = Settings("system.yaml")
prn = true
sam = SymbolicAWEModel(set, "ram")
init!(sam; prn)
tsam = SymbolicAWEModel(set, "tether")
init!(tsam; prn)
ssam = SymbolicAWEModel(set, "simple_ram")
init!(ssam; prn)
set.segments = 1
sam = SymbolicAWEModel(set, "ram")
init!(sam; prn)
tsam = SymbolicAWEModel(set, "tether")
init!(tsam; prn)
ssam = SymbolicAWEModel(set, "simple_ram")
init!(ssam; prn)
