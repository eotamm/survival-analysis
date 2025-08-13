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

# Feature selection and survival modeling across data transformations
compare_transforms_with_brm <- function(tse, pseudocount = 1e-6, alr_ref = "g__Sutterella") {
  prevalence <- getPrevalence(tse, detection = 1)
  results <- list()
  
  ntaxa <- nrow(tse)
  max_features <- min(ntaxa, 1000)
  feature_counts <- 2^(0:floor(log2(max_features)))
  if (tail(feature_counts, 1) < max_features) {
    feature_counts <- c(feature_counts, max_features)
  }
  
  set.seed(2025)
  idx <- sample(seq_len(ncol(tse)), size = floor(ncol(tse) * 0.8))
  tse_train <- tse[, idx]
  tse_test  <- tse[, -idx]
  
  #  Define transformations
  transforms <- list(
    clr = function(x) transformAssay(x, method = "clr", assay.type = "counts", pseudocount = pseudocount, name = "clr"),
    rclr = function(x) transformAssay(x, method = "rclr", assay.type = "counts", pseudocount = pseudocount, name = "rclr"),
    lra = function(x) {
      x |>
        transformAssay(method = "relabundance", assay.type = "counts", name = "relabundance") |>
        transformAssay(method = "log", assay.type = "relabundance", pseudocount = pseudocount, name = "log_abund") |>
        transformAssay(method = "difference", assay.type = "log_abund", name = "logratios", MARGIN = 1L)
    },
    pa = function(x) {
      x <- transformAssay(x, method = "relabundance", assay.type = "counts", name = "tss")
      ra <- assay(x, "tss")
      assay(x, "pa") <- (ra >= 0.001) * 1  # 0.1%
      x
    },
    tss = function(x) transformAssay(x, method = "relabundance", assay.type = "counts", name = "tss"),
    logtss = function(x) {
      x |>
        transformAssay(method = "relabundance", assay.type = "counts", name = "tss") |>
        transformAssay(method = "log", assay.type = "tss", pseudocount = pseudocount, name = "logtss")
    },
    asin = function(x) {
      x_tss <- transformAssay(x, method = "relabundance", assay.type = "counts", name = "tss")
      asin_mat <- asin(sqrt(assay(x_tss, "tss")))
      altExp(x_tss, "asin") <- SummarizedExperiment(assays = list(asin = asin_mat))
      return(x_tss)
    },
    alr = function(x) transformAssay(x, method = "alr", assay.type = "counts", pseudocount = pseudocount, name = "alr", ref_taxa = alr_ref)
  )
  
  result_ntaxa <- list()
  
  for (tname in names(transforms)) {
    cat("→", tname, "\n")
    
    tse_train_tx <- transforms[[tname]](tse_train)
    tse_test_tx  <- transforms[[tname]](tse_test)
    
    # Convert assay data to data frame
    if (tname == "lra") {
      df_train <- as.data.frame(t(as.matrix(assay(altExp(tse_train_tx, "logratios")))))
      df_test  <- as.data.frame(t(as.matrix(assay(altExp(tse_test_tx, "logratios")))))
    } else if (tname == "asin") {
      df_train <- as.data.frame(t(as.matrix(assay(altExp(tse_train_tx, "asin")))))
      df_test  <- as.data.frame(t(as.matrix(assay(altExp(tse_test_tx, "asin")))))
    } else {
      df_train <- as.data.frame(t(assay(tse_train_tx, tname)))
      df_test  <- as.data.frame(t(assay(tse_test_tx, tname)))
    }
    
    # Add surviavl outcome
    df_train$Event <- colData(tse_train_tx)$Event
    df_train$Event_time <- colData(tse_train_tx)$Event_time
    df_test$Event <- colData(tse_test_tx)$Event
    df_test$Event_time <- colData(tse_test_tx)$Event_time
    
    df_train <- clean_column_names(df_train)
    df_test  <- clean_column_names(df_test)
    
    predictors <- setdiff(names(df_train), c("Event", "Event_time"))
    
    stats <- lapply(predictors, function(p) {
      tryCatch({
        s <- summary(coxph(as.formula(paste0("Surv(Event_time, Event) ~ ", p)), data = df_train))
        data.frame(feature = p, p = s$logtest["pvalue"])
      }, error = function(e) data.frame(feature = p, p = 1))
    }) %>%
      bind_rows() %>%
      arrange(p)
    
    res_feature_counts <- data.frame(N = feature_counts, C_index = NA, features = I(vector("list", length(feature_counts))))
    
    for (i in seq_along(feature_counts)) {
      N <- feature_counts[i]
      feats <- stats$feature[1:min(N, nrow(stats))]
      df_subset <- df_train[, c(feats, "Event", "Event_time")]
      
      if (nrow(na.omit(df_subset)) == 0) {
        warning("→ Ei kelvollisia havaintoja transformaatiossa ", tname, " ja N = ", N, ". Hypätään yli.")
        next
      }
      
      formula_str <- paste("Event_time | cens(1 - Event) ~", paste(feats, collapse = " + "))
      message("   Fitting brms model with ", N, " features (", tname, ")...")
      
      fit <- brm(
        formula = as.formula(formula_str),
        data = df_subset,
        family = brmsfamily("cox"),
        prior = c(
          set_prior("normal(0, 1)", class = "Intercept"),
          set_prior("normal(0, 1)", class = "b")
        ),
        chains = 4,
        iter = 4000,
        control = list(adapt_delta = 0.99,
                       max_treedepth = 14),
        cores = 4,
        silent = 2
      )
      
      pred <- posterior_linpred(fit, newdata = df_test[, c(feats, "Event", "Event_time")], transform = FALSE)
      lp_mean <- colMeans(pred)
      cidx <- survConcordance(Surv(df_test$Event_time, df_test$Event) ~ lp_mean)$concordance
      res_feature_counts$C_index[i] <- cidx
      res_feature_counts$features[[i]] <- feats
    }
    
    result_ntaxa[[tname]] <- res_feature_counts
  }
  
  results[[paste0("ntaxa_", ntaxa)]] <- result_ntaxa
  return(results)
}



# Model fitting function with selectable prior
fit_model <- function(df, model_name, prior_type = c("normal", "horseshoe")) {
  prior_type <- match.arg(prior_type)
  
  # Define predictor variables
  predictors <- setdiff(names(df), c("Event", "Event_time"))
  
  # Build Cox model formula with censoring specification
  formula_str <- paste("Event_time | cens(1 - Event) ~", paste(predictors, collapse = " + "))
  message("Fitting model with ", prior_type, " prior: ", model_name)
  
  # Define prior based on type
  prior <- switch(prior_type,
                  normal = c(
                    set_prior("normal(0, 1)", class = "Intercept"),
                    set_prior("normal(0, 1)", class = "b")
                  ),
                  horseshoe = c(
                    set_prior("normal(0, 1)", class = "Intercept"),
                    set_prior(horseshoe(df = 1, par_ratio = 0.1), class = "b")
                  )
  )
  
  # Fit the model
  fit <- brm(
    formula = as.formula(formula_str),
    data = df,
    family = brmsfamily("cox"),
    prior = prior,
    chains = 4,
    iter = 4000,
    control = list(adapt_delta = 0.99, max_treedepth = 16),
    cores = parallel::detectCores(),
    save_pars = save_pars(all = TRUE)
  )
  
  # Save the model
  saveRDS(fit, file.path("model_result", paste0(model_name, ".rds")))
  return(fit)
}


# Prepare data for analysis
extract_top_features_by_transformation <- function(tse,
                                                   method = c("clr", "rclr", "log_abund", "lra", "pa", "tss", "logtss", "asin", "alr"),
                                                   top_n = 20,
                                                   pseudocount = 1e-6,
                                                   alr_ref = "g__Sutterella") {
  method <- match.arg(method)
  
  # Apply transformation
  tse_tx <- switch(method,
                   clr = transformAssay(tse, method = "clr", assay.type = "counts", pseudocount = pseudocount, name = "clr"),
                   rclr = transformAssay(tse, method = "rclr", assay.type = "counts", pseudocount = pseudocount, name = "rclr"),
                   log_abund = tse |>
                     transformAssay(method = "relabundance", assay.type = "counts", name = "relabundance") |>
                     transformAssay(method = "log", assay.type = "relabundance", pseudocount = pseudocount, name = "log_abund"),
                   lra = tse |>
                     transformAssay(method = "relabundance", assay.type = "counts", name = "relabundance") |>
                     transformAssay(method = "log", assay.type = "relabundance", pseudocount = pseudocount, name = "log_abund") |>
                     transformAssay(method = "difference", assay.type = "log_abund", name = "logratios", MARGIN = 1L),
                   pa = {
                     out <- transformAssay(tse, method = "relabundance", assay.type = "counts", name = "tss")
                     ra  <- assay(out, "tss")
                     assay(out, "pa") <- (ra >= 0.001) * 1 #0.1%
                     out
                   },
                   tss = transformAssay(tse, method = "relabundance", assay.type = "counts", name = "tss"),
                   logtss = tse |>
                     transformAssay(method = "relabundance", assay.type = "counts", name = "tss") |>
                     transformAssay(method = "log", assay.type = "tss", pseudocount = pseudocount, name = "logtss"),
                   asin = {
                     tse_tmp <- transformAssay(tse, method = "relabundance", assay.type = "counts", name = "tss")
                     asin_mat <- asin(sqrt(assay(tse_tmp, "tss")))
                     attr(asin_mat, "Event") <- colData(tse_tmp)$Event
                     attr(asin_mat, "Event_time") <- colData(tse_tmp)$Event_time
                     asin_mat
                   },
                   alr = transformAssay(tse, method = "alr", assay.type = "counts", pseudocount = pseudocount, name = "alr", ref_taxa = alr_ref)
  )
  
  # Extract assay data into dataframe
  df <- if (method == "lra") {
    as.data.frame(t(as.matrix(assay(altExp(tse_tx, "logratios")))))
  } else if (method == "asin") {
    as.data.frame(t(tse_tx))
  } else {
    as.data.frame(t(assay(tse_tx, method)))
  }
  
  # Add survival metadata
  if (method == "asin") {
    df$Event <- attr(tse_tx, "Event")
    df$Event_time <- attr(tse_tx, "Event_time")
  } else {
    df$Event <- colData(tse_tx)$Event
    df$Event_time <- colData(tse_tx)$Event_time
  }
  
  df <- clean_column_names(df)
  
  # Univariate Cox regression
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

# Plot top posterior hazard ratios with 95% credible intervals
plot_top_hr <- function(fit, title = "Posterior estimates (HR)", top_k = 10) {
  
  ps <- posterior::as_draws_df(fit) %>%
    select(starts_with("b_")) %>%
    select(-any_of("b_Intercept")) %>%
    pivot_longer(everything(), names_to = "variable", values_to = "value") %>%
    group_by(variable) %>%
    summarise(
      median = median(value),
      lower  = quantile(value, 0.025),
      upper  = quantile(value, 0.975),
      .groups = "drop"
    ) %>%
    mutate(
      variable_clean = str_remove(variable, "^b_"),
      significant    = ifelse(lower > 0 | upper < 0, "yes", "no"),
      hr_median      = exp(median),
      hr_lower       = exp(lower),
      hr_upper       = exp(upper)
    ) %>%
    arrange(desc(median)) %>% 
    slice_head(n = top_k) %>%
    mutate(variable_clean = factor(variable_clean, levels = rev(variable_clean)))
  
  ggplot(ps, aes(
    y = variable_clean,
    x = hr_median,
    xmin = hr_lower,
    xmax = hr_upper,
    color = significant
  )) +
    geom_point() +
    geom_errorbarh(height = 0.2) +
    geom_vline(xintercept = 1, linetype = "dashed", color = "gray") +
    scale_color_manual(values = c("yes" = "red", "no" = "black")) +
    labs(
      x = "Hazard Ratio (exp(median))",
      y = "",
      title = title,
      color = "Significant (95% CI excludes 0)"
    ) +
    theme_minimal()
}

# Combine LOO results from both priors into one data frame
get_elpd_df <- function(loo_list, prior_label) {
  data.frame(
    model = names(loo_list),
    elpd = sapply(loo_list, function(x) x$estimates["elpd_loo", "Estimate"]),
    se   = sapply(loo_list, function(x) x$estimates["elpd_loo", "SE"]),
    prior = prior_label
  )
}

# Function to extract median risk score per model
get_pred_df <- function(model, method_name) {
  pred <- posterior_linpred(model, transform = FALSE)
  risk_score <- apply(pred, 2, median)
  data.frame(
    risk_score = risk_score,
    Event = model$data$Event,
    Event_time = model$data$Event_time,
    Method = method_name
  )
}

# Compute AUROC from all posterior samples
get_auc_posterior <- function(model, method_name, prior_label) {
  pred_mat <- posterior_linpred(model, transform = FALSE)
  labels <- model$data$Event
  
  auc_vals <- apply(pred_mat, 1, function(scores) {
    roc_obj <- pROC::roc(response = labels, predictor = scores, quiet = TRUE)
    as.numeric(roc_obj$auc)
  })
  
  tibble(
    model = method_name,
    prior = prior_label,
    AUROC = auc_vals
  )
}


# Compute C-index
get_c_posterior <- function(model, method_name, prior_label) {
  pred_mat <- posterior_linpred(model, transform = FALSE)
  time  <- model$data$Event_time
  event <- model$data$Event
  c_vals <- apply(pred_mat, 1, function(lp) {
    ok <- is.finite(lp) & is.finite(time) & is.finite(event)
    if (!any(ok)) return(NA_real_)
    cf <- survival::concordance(Surv(time[ok], event[ok]) ~ lp[ok], reverse = TRUE)
    as.numeric(cf$concordance)
  })
  tibble(
    model  = method_name,
    prior  = prior_label,
    Cindex = c_vals
  )
}


# Count Más-o-menos riskscore
calculate_masomenos <- function(df) {
  predictors <- setdiff(names(df), c("Event", "Event_time"))
  signs <- sapply(predictors, function(var) {
    mod <- coxph(Surv(Event_time, Event) ~ df[[var]], data = df)
    sign(coef(mod))
  })
  score <- as.matrix(df[, predictors]) %*% signs / length(signs)
  return(as.numeric(score))
}
