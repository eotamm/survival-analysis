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

# Model fitting function
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


# Prepare data for analysis
extract_top_features_by_transformation <- function(tse, method = c("clr", "rclr", "log_abund", "lra"),
                                                   top_n = 20, pseudocount = 1e-6) {
  method <- match.arg(method)
  
  # Apply selected transformation
  tse_tx <- switch(method,
                   clr = transformAssay(tse, method = "clr", assay.type = "counts", pseudocount = pseudocount, name = "clr"),
                   rclr = transformAssay(tse, method = "rclr", assay.type = "counts", pseudocount = pseudocount, name = "rclr"),
                   log_abund = tse |>
                     transformAssay(method = "relabundance", assay.type = "counts", name = "relabundance") |>
                     transformAssay(method = "log", assay.type = "relabundance", pseudocount = pseudocount, name = "log_abund"),
                   lra = tse |>
                     transformAssay(method = "relabundance", assay.type = "counts", name = "relabundance") |>
                     transformAssay(method = "log", assay.type = "relabundance", pseudocount = pseudocount, name = "log_abund") |>
                     transformAssay(method = "difference", assay.type = "log_abund", name = "logratios")
  )
  
  # Extract transformed assay
  df <- if (method == "lra") {
    as.data.frame(t(as.matrix(assay(altExp(tse_tx, "logratios")))))
  } else {
    as.data.frame(t(assay(tse_tx, method)))
  }
  
  # Add survival data
  df$Event <- colData(tse_tx)$Event
  df$Event_time <- colData(tse_tx)$Event_time
  df <- clean_column_names(df)
  
  # Univariate Cox per feature
  predictors <- setdiff(names(df), c("Event", "Event_time"))
  stats <- lapply(predictors, function(p) {
    tryCatch({
      s <- summary(coxph(as.formula(paste0("Surv(Event_time, Event) ~ ", p)), data = df))
      data.frame(feature = p, p = s$logtest["pvalue"])
    }, error = function(e) data.frame(feature = p, p = 1))
  }) %>%
    bind_rows() %>%
    arrange(p)
  
  # Select top N features
  top_feats <- stats$feature[1:min(top_n, nrow(stats))]
  df_top <- df[, c(top_feats, "Event", "Event_time")]
  
  return(df_top)
}



# Feature selection and survival modeling across data transformations
compare_transforms_with_brm <- function(tse) {
  # Compute prevalence of each taxon
  prevalence <- getPrevalence(tse, detection = 1)
  results <- list()
  
  # Determine number of taxa and max feature count
  ntaxa <- nrow(tse)
  max_features <- min(ntaxa, 2000)
  
  # Define feature counts (powers of 2) to test
  feature_counts <- 2^(0:floor(log2(max_features)))
  if (tail(feature_counts, 1) < max_features) {
    feature_counts <- c(feature_counts, max_features)
  }
  
  # Split data into 80% training and 20% testing
  set.seed(2025)
  idx <- sample(seq_len(ncol(tse)), size = floor(ncol(tse) * 0.8))
  tse_train <- tse[, idx]
  tse_test  <- tse[, -idx]
  
  # Define transformations to test
  transforms <- list(
    log_abund = function(x) {
      x |>
        transformAssay(method = "relabundance", assay.type = "counts", name = "relabundance") |>
        transformAssay(method = "log", assay.type = "relabundance", pseudocount = 1e-6, name = "log_abund")
    },
    clr = function(x) {
      transformAssay(x, method = "clr", assay.type = "counts", pseudocount = 1e-6, name = "clr")
    },
    rclr = function(x) {
      transformAssay(x, method = "rclr", assay.type = "counts", pseudocount = 1e-6, name = "rclr")
    },
    lra = function(x) {
      x |>
        transformAssay(method = "relabundance", assay.type = "counts", name = "relabundance") |>
        transformAssay(method = "log", assay.type = "relabundance", pseudocount = 1e-6, name = "log_abund") |>
        transformAssay(method = "difference", assay.type = "log_abund", name = "logratios")
    }
  )
  
  result_ntaxa <- list()
  
  for (tname in names(transforms)) {
    cat("→", tname, "\n")
    
    # Apply transformation
    tse_train_tx <- transforms[[tname]](tse_train)
    tse_test_tx  <- transforms[[tname]](tse_test)
    
    # Extract transformed matrices
    df_train <- if (tname == "lra") as.data.frame(t(as.matrix(assay(altExp(tse_train_tx, "logratios"))))) else
      as.data.frame(t(assay(tse_train_tx, tname)))
    df_test  <- if (tname == "lra") as.data.frame(t(as.matrix(assay(altExp(tse_test_tx, "logratios"))))) else
      as.data.frame(t(assay(tse_test_tx, tname)))
    
    # Add survival data
    df_train$Event <- colData(tse_train_tx)$Event
    df_train$Event_time <- colData(tse_train_tx)$Event_time
    df_test$Event <- colData(tse_test_tx)$Event
    df_test$Event_time <- colData(tse_test_tx)$Event_time
    
    # Clean column names
    df_train <- clean_column_names(df_train)
    df_test <- clean_column_names(df_test)
    
    # Get predictors
    predictors <- setdiff(names(df_train), c("Event", "Event_time"))
    
    # Univariate Cox regression for feature ranking
    stats <- lapply(predictors, function(p) {
      tryCatch({
        s <- summary(coxph(as.formula(paste0("Surv(Event_time, Event) ~ ", p)), data = df_train))
        data.frame(feature = p, p = s$logtest["pvalue"])
      }, error = function(e) data.frame(feature = p, p = 1))
    }) %>%
      bind_rows() %>%
      arrange(p)
    
    # Fit BRMS models with increasing number of features
    res_feature_counts <- data.frame(N = feature_counts, C_index = NA, features = I(vector("list", length(feature_counts))))
    
    for (i in seq_along(feature_counts)) {
      N <- feature_counts[i]
      feats <- stats$feature[1:min(N, nrow(stats))]
      
      formula_str <- paste("Event_time | cens(1 - Event) ~", paste(feats, collapse = " + "))
      message("   Fitting brms model with ", N, " features (", tname, ")...")
      
      fit <- brm(
        formula = as.formula(formula_str),
        data = df_train[, c(feats, "Event", "Event_time")],
        family = brmsfamily("cox"),
        prior = set_prior("normal(0, 1)", class = "b"),
        chains = 4,
        iter = 4000,
        cores = 2,
        silent = 2
      )
      
      # Predict on test set and compute C-index
      pred <- posterior_linpred(fit, newdata = df_test[, c(feats, "Event", "Event_time")], transform = FALSE)
      lp_mean <- colMeans(pred)
      cidx <- survConcordance(Surv(df_test$Event_time, df_test$Event) ~ lp_mean)$concordance
      res_feature_counts$C_index[i] <- cidx
      res_feature_counts$features[[i]] <- feats
    }
    
    # Save results for this transformation
    result_ntaxa[[tname]] <- res_feature_counts
  }
  
  results[[paste0("ntaxa_", ntaxa)]] <- result_ntaxa
  return(results)
}

