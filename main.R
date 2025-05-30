# Load and preprocess the data
source("data.R")

# Fit the survival model
source("model.R")

# Render the Quarto report
quarto::quarto_render("survival.qmd")
