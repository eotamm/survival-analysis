# Clean column names
clean_column_names <- function(df) {
  clean_names <- names(df)
  clean_names <- gsub("__+", "_", clean_names)            
  clean_names <- gsub("_$", "", clean_names)              
  clean_names <- gsub("[^A-Za-z0-9_]", "", clean_names)     
  clean_names <- ifelse(grepl("^[A-Za-z]", clean_names), 
                        clean_names, 
                        paste0("X", clean_names))         
  clean_names <- make.unique(clean_names)               
  names(df) <- clean_names
  return(df)
}

# Model fitting function for ratios
fit_model <- function(df, model_name) {
  # Define predictor variables
  predictors <- setdiff(names(df), c("Event", "Event_time"))
  
  # Build Cox model formula with censoring specification
  formula_str <- paste("Event_time | cens(1 - Event) ~", paste(predictors, collapse = " + "))
  message("Fitting model: ", model_name)
  
  # Fit the model using Bayesian Cox regression via brms
  fit <- brm(
    formula = as.formula(formula_str),
    data = df,
    family = brmsfamily("cox"),
    prior = set_prior("normal(0, 1)", class = "b"),
    chains = 4,
    iter = 4000,
    control = list(max_treedepth = 15),
    cores = parallel::detectCores()
  )
  
  # Save the fitted model to disk
  saveRDS(fit, paste0(model_name, ".rds"))
  return(fit)
}

