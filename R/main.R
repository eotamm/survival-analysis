# Utility functions
source("funcs.R")

# Load and preprocess the data
source("data.R")

# TabPfn function and fit the survival models
source("tabpfnfunc.R")
source("model.R")

# Render the Quarto report
quarto::quarto_render("survival.qmd")

