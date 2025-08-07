## Step 1. Include necessary packages
using SymbolicAWEModels
using KiteUtils
using WinchModels
using ControlPlots

## Step 2. Define settings, in .yaml format
# This file should be stotrd in data/<simulation_name>/settings.yaml
# We will use the example, shown in data/tutorial_1
# The settings are loaded using KiteUtils/src/settings.jl
# Which looks in the data directory for the settings file
set = Settings("tutorial_1/system.yaml")

### Step 3. Initialize the model

# set.l_tether will set the firt index of set.l_tethers