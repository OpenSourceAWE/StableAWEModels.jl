# SPDX-FileCopyrightText: 2025 Bart van de Lint <bart@vandelint.net>
#
# SPDX-License-Identifier: MPL-2.0

using SymbolicAWEModels

set = Settings("system.yaml")
prn = true
sam = SymbolicAWEModel(set, "ram")
init!(sam; prn)
tether_sam = SymbolicAWEModel(set, "tether")
init!(tether_sam; prn)
simple_sam = SymbolicAWEModel(set, "simple_ram")
init!(simple_sam; prn)
set.segments = 1
sam = SymbolicAWEModel(set, "ram")
init!(sam; prn)
tether_sam = SymbolicAWEModel(set, "tether")
init!(tether_sam; prn)
