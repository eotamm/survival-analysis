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



# Function to prepare data for analysis
prepare_logratios_from_clr <- function(tse_subset, method = c("clr", "rclr")) {
  method <- match.arg(method)
  
  # Perform CLR or rCLR transformation using mia
  tse_subset <- transformAssay(
    tse_subset,
    method = method,
    assay.type = "counts",
    pseudocount = 1e-6,
    name = method
  )
  
  # Extract the transformed matrix (samples x microbes)
  mat <- t(assay(tse_subset, method))
  
  # Calculate log-ratios: log(x/y) = clr(x) - clr(y) for each microbe pair
  microbe_names <- colnames(mat)
  combs <- combn(microbe_names, 2, simplify = FALSE)
  
  df <- purrr::map_dfc(combs, function(pair) {
    col1 <- pair[1]
    col2 <- pair[2]
    ratio_name <- paste0("logratio_", col1, "_vs_", col2)
    tibble(!!ratio_name := mat[, col1] - mat[, col2])
  })
  
  # Standardize features (mean = 0, sd = 1)
  df <- as.data.frame(scale(df))
  
  # Add response variables: event indicator and survival time
  df$Event <- colData(tse_subset)$Event
  df$Event_time <- colData(tse_subset)$Event_time
  
  # Clean up column names
  df <- clean_column_names(df)

  # Remove rows with missing or invalid survival times
  df <- df %>%
    filter(Event_time >= 0) %>%
    drop_na()
  
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
