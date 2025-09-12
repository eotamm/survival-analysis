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
extract_all_features_by_transformation <- function(tse,
                                                   method = c("clr","rclr","log_abund","lra","pa","tss","logtss","asin","alr"),
                                                   pseudocount = 1e-6,
                                                   alr_ref = "g_Turicibacter") {
  method <- match.arg(method)
  
  # Apply transformation
  tse_tx <- switch(
    method,
    clr  = transformAssay(tse, method = "clr",  assay.type = "counts", pseudocount = pseudocount, name = "clr"),
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
      assay(out, "pa") <- (ra >= 0.001) * 1  # 0.1%
      out
    },
    tss    = transformAssay(tse, method = "relabundance", assay.type = "counts", name = "tss"),
    logtss = tse |>
      transformAssay(method = "relabundance", assay.type = "counts", name = "tss") |>
      transformAssay(method = "log", assay.type = "tss", pseudocount = pseudocount, name = "logtss"),
    asin = {
      tmp <- transformAssay(tse, method = "relabundance", assay.type = "counts", name = "tss")
      asin_mat <- asin(sqrt(assay(tmp, "tss")))
      attr(asin_mat, "Event")      <- colData(tmp)$Event
      attr(asin_mat, "Event_time") <- colData(tmp)$Event_time
      asin_mat
    },
    alr = transformAssay(tse, method = "alr", assay.type = "counts",
                         pseudocount = pseudocount, name = "alr", ref_taxa = alr_ref)
  )
  
  # Extract assay matrix
  M <- if (method == "lra") {
    assay(altExp(tse_tx, "logratios"))
  } else if (method == "asin") {
    tse_tx
  } else {
    assay(tse_tx, method)
  }
  
  df <- as.data.frame(t(as.matrix(M)))
  
  # Clean column names if helper exists
  if (exists("clean_column_names")) df <- clean_column_names(df)
  
  # ALR: drop the reference taxon column so it is not used as a feature
  if (method == "alr") {
    canon <- function(x) gsub("[^a-z0-9]+", "", tolower(x))
    ref_variants <- unique(c(
      alr_ref,
      sub("__", "_", alr_ref, fixed = TRUE),
      sub("^g__", "g_", alr_ref),
      make.names(alr_ref)
    ))
    ref_canon <- canon(ref_variants)
    drop_idx <- which(canon(colnames(df)) %in% ref_canon)
    if (length(drop_idx)) df <- df[, -drop_idx, drop = FALSE]
  }
  
  # Remove all-NA columns
  keep <- colSums(!is.na(df)) > 0
  df <- df[, keep, drop = FALSE]
  
  # Add survival metadata
  if (method == "asin") {
    df$Event      <- attr(M, "Event")
    df$Event_time <- attr(M, "Event_time")
  } else {
    df$Event      <- colData(tse_tx)$Event
    df$Event_time <- colData(tse_tx)$Event_time
  }
  
  df
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



# Compute Harrell's C-index for given times, events, and risk scores
harrell_c <- function(time, event, score, reverse = TRUE) {
  as.numeric(
    survival::concordance(
      survival::Surv(time, event) ~ score,
      reverse = reverse
    )$concordance
  )
}

# Create stratified K-fold indices based on binary Event (0/1)
make_stratified_folds <- function(event, K = 5, seed = 1) {
  set.seed(seed)
  n <- length(event)
  f <- integer(n)
  idx0 <- which(event == 0L)
  idx1 <- which(event == 1L)
  f[idx0] <- sample(rep(seq_len(K), length.out = length(idx0)))
  f[idx1] <- sample(rep(seq_len(K), length.out = length(idx1)))
  f
}

# Run generic K-fold CV for a model that returns out-of-fold risk scores
cv5_cindex <- function(df, risk_fun, seed = 1, reverse = TRUE, K = 5, fold_id = NULL) {
  n <- nrow(df)
  if (!is.null(fold_id) && length(fold_id) != n)
    stop("'fold_id' length must equal nrow(df)")
  if (is.null(fold_id))
    fold_id <- make_stratified_folds(df$Event, K = K, seed = seed)
  
  preds <- numeric(n)
  for (k in seq_len(K)) {
    train <- df[fold_id != k, , drop = FALSE]
    test  <- df[fold_id == k, , drop = FALSE]
    preds[fold_id == k] <- as.numeric(risk_fun(train, test))
  }
  
  c_overall <- harrell_c(df$Event_time, df$Event, preds, reverse = reverse)
  c_folds <- vapply(
    seq_len(K),
    function(k) harrell_c(
      df$Event_time[fold_id == k],
      df$Event[fold_id == k],
      preds[fold_id == k],
      reverse = reverse
    ),
    numeric(1)
  )
  
  list(
    c_overall = c_overall,
    c_folds   = c_folds,
    preds     = preds,
    fold_id   = fold_id,
    data      = df,
    reverse   = reverse,
    K         = K
  )
}

# Summarize a CV run to a tibble with estimate and normal-approx CI
summarize_cv <- function(model_name, method_name, cv_res, level = 0.95) {
  cc <- survival::concordance(
    survival::Surv(cv_res$data$Event_time, cv_res$data$Event) ~ cv_res$preds,
    reverse = isTRUE(cv_res$reverse)
  )
  est <- as.numeric(cc$concordance)
  se  <- if (!is.null(cc$std.err)) as.numeric(cc$std.err) else sqrt(as.numeric(cc$var))
  z   <- stats::qnorm(1 - (1 - level) / 2)
  lo  <- est - z * se
  hi  <- est + z * se
  
  tibble::tibble(
    model    = model_name,
    method   = method_name,
    metric   = "C",
    estimate = est,
    lower    = lo,
    upper    = hi,
    se       = se
  )
}

# Random Survival Forest 5-fold CV
rsf_cindex_cv5 <- function(
    df, method_name,
    num.trees = 1500, min.node.size = 10, mtry = NULL,
    seed = 1, K = 5, fold_id = NULL
) {
  message(sprintf("[RSF] %s: starting %d-fold CV", method_name, K))
  
  risk_fun <- function(train, test) {
    p_train  <- max(1L, ncol(train) - 2L)
    mtry_use <- if (is.null(mtry)) max(1L, min(p_train, floor(sqrt(p_train)))) else as.integer(mtry)
    fit <- ranger::ranger(
      survival::Surv(Event_time, Event) ~ .,
      data = train,
      num.trees = num.trees,
      mtry = mtry_use,
      min.node.size = min.node.size,
      splitrule = "logrank",
      write.forest = TRUE
    )
    chf <- predict(fit, data = test)$chf
    if (is.matrix(chf)) chf[, ncol(chf), drop = TRUE] else as.numeric(chf)
  }
  
  res <- cv5_cindex(df, risk_fun, seed = seed, reverse = TRUE, K = K, fold_id = fold_id)
  message(sprintf("[RSF] %s: done. CV%d C = %.3f", method_name, K, res$c_overall))
  summarize_cv(method_name, sprintf("RSF_%dCV", K), res)
}

# Logistic regression (binomial) 5-fold CV
logit_cindex_cv5 <- function(
    df, method_name, maxit = 200,
    seed = 1, K = 5, fold_id = NULL
) {
  message(sprintf("[Logit] %s: starting %d-fold CV", method_name, K))
  
  risk_fun <- function(train, test) {
    xtr_df <- train[, setdiff(names(train), c("Event_time","Event")), drop = FALSE]
    xte_df <- test[,  setdiff(names(test),  c("Event_time","Event")), drop = FALSE]
    MM     <- stats::model.matrix(~ . - 1, data = rbind(xtr_df, xte_df))
    ntr    <- nrow(xtr_df)
    Xtr    <- MM[seq_len(ntr), , drop = FALSE]
    Xte    <- MM[(ntr + 1L):nrow(MM), , drop = FALSE]
    ytr    <- as.integer(train$Event)
    
    fit <- stats::glm(
      ytr ~ .,
      data   = data.frame(Xtr, check.names = FALSE),
      family = stats::binomial(),
      control = list(maxit = maxit)
    )
    as.numeric(stats::predict(fit, newdata = data.frame(Xte, check.names = FALSE), type = "response"))
  }
  
  res <- cv5_cindex(df, risk_fun, seed = seed, reverse = TRUE, K = K, fold_id = fold_id)
  message(sprintf("[Logit] %s: done. CV%d C = %.3f", method_name, K, res$c_overall))
  summarize_cv(method_name, sprintf("Logit_%dCV", K), res)
}

# Cox proportional hazards 5-fold CV
coxph_cindex_cv5 <- function(
    df, method_name,
    ties = c("efron","breslow","exact"),
    iter_max = 50,
    seed = 1, K = 5, fold_id = NULL
) {
  ties <- match.arg(ties)
  message(sprintf("[CoxPH] %s: starting %d-fold CV (ties=%s)", method_name, K, ties))
  
  risk_fun <- function(train, test) {
    fit <- survival::coxph(
      survival::Surv(Event_time, Event) ~ .,
      data = train,
      ties = ties,
      control = survival::coxph.control(iter.max = iter_max),
      x = FALSE, y = FALSE
    )
    as.numeric(stats::predict(fit, newdata = test, type = "lp"))
  }
  
  res <- cv5_cindex(df, risk_fun, seed = seed, reverse = TRUE, K = K, fold_id = fold_id)
  message(sprintf("[CoxPH] %s: done. CV%d C = %.3f", method_name, K, res$c_overall))
  summarize_cv(method_name, sprintf("CoxPH_%dCV", K), res)
}

# DeepSurv, 5-fold CV
deepsurv_cindex_cv5 <- function(
    df, method_name,
    hidden = c(64, 32), dropout = 0.0, l2 = 1e-4,
    lr = 1e-3, epochs = 300, patience = 30,
    verbose = 0, run_eagerly = TRUE,
    seed = 1, K = 5, fold_id = NULL
) {
  message(sprintf("[DeepSurv] %s: starting %d-fold CV", method_name, K))
  
  risk_fun <- function(train, test) {
    xtr_df <- train[, setdiff(names(train), c("Event_time","Event")), drop = FALSE]
    xte_df <- test[,  setdiff(names(test),  c("Event_time","Event")), drop = FALSE]
    mm_all <- stats::model.matrix(~ . - 1, data = rbind(xtr_df, xte_df))
    ntr    <- nrow(xtr_df)
    Xtr    <- mm_all[seq_len(ntr), , drop = FALSE]
    Xte    <- mm_all[(ntr + 1L):nrow(mm_all), , drop = FALSE]
    p      <- ncol(Xtr)
    
    mu <- matrixStats::colMeans2(Xtr)
    sd <- matrixStats::colSds(Xtr); sd[!is.finite(sd) | sd < 1e-8] <- 1
    Xtr_s <- sweep(sweep(Xtr, 2, mu, "-"), 2, sd, "/")
    Xte_s <- sweep(sweep(Xte, 2, mu, "-"), 2, sd, "/")
    
    ord   <- order(train$Event_time)
    Xtr_s <- Xtr_s[ord, , drop = FALSE]
    ev_tr <- as.numeric(train$Event[ord])
    
    cox_ph_loss_safe <- function(y_true, y_pred) {
      tf <- tensorflow::tf
      y_pred <- tf$reshape(y_pred, shape = c(-1L))
      events <- tf$reshape(y_true, shape = c(-1L))
      haz <- tf$math$exp(y_pred)
      rev_csum  <- tf$math$cumsum(tf$reverse(haz, list(0L)))
      risk_csum <- tf$reverse(rev_csum, list(0L))
      log_risk  <- tf$math$log(risk_csum + 1e-8)
      -tf$math$reduce_sum((y_pred - log_risk) * events) /
        (tf$math$reduce_sum(events) + 1e-8)
    }
    
    inp <- keras::layer_input(shape = p, dtype = "float32")
    x <- inp
    for (u in hidden) {
      x <- keras::layer_dense(
        x, units = u, activation = "relu",
        kernel_regularizer = keras::regularizer_l2(l = l2)
      )
      if (dropout > 0) x <- keras::layer_dropout(x, rate = dropout)
    }
    out <- keras::layer_dense(
      x, units = 1, activation = "linear",
      kernel_regularizer = keras::regularizer_l2(l = l2)
    )
    model <- keras::keras_model(inp, out)
    model %>% keras::compile(
      optimizer = keras::optimizer_adam(learning_rate = lr, clipnorm = 1.0, clipvalue = 0.5),
      loss = cox_ph_loss_safe,
      run_eagerly = run_eagerly
    )
    model %>% keras::fit(
      x = Xtr_s, y = matrix(ev_tr, ncol = 1),
      batch_size = nrow(Xtr_s), epochs = epochs,
      shuffle = FALSE, verbose = verbose,
      callbacks = list(keras::callback_early_stopping(
        monitor = "loss", patience = patience, restore_best_weights = TRUE))
    )
    as.numeric(model$predict(Xte_s, verbose = as.integer(verbose)))
  }
  
  res <- cv5_cindex(df, risk_fun, seed = seed, reverse = TRUE, K = K, fold_id = fold_id)
  message(sprintf("[DeepSurv] %s: done. CV%d C = %.3f", method_name, K, res$c_overall))
  summarize_cv(method_name, sprintf("DeepSurv_%dCV", K), res)
}

# CatBoost (binary)
catboost_bin_cindex_cv5 <- function(
    df, method_name,
    seed = 1, K = 5, fold_id = NULL,
    iterations = 1000, depth = 6, learning_rate = 0.05,
    l2_leaf_reg = 3.0, border_count = 254, verbose = 0
) {
  message(sprintf("[CatBoost Bin→C] %s: starting %d-fold CV", method_name, K))
  
  risk_fun <- function(train, test) {
    xtr <- train[, setdiff(names(train), c("Event_time","Event")), drop = FALSE]
    xte <- test[,  setdiff(names(test),  c("Event_time","Event")), drop = FALSE]
    
    xtr[] <- lapply(xtr, function(v) if (is.character(v) || is.logical(v)) factor(v) else v)
    xte[] <- lapply(xte, function(v) if (is.character(v) || is.logical(v)) factor(v) else v)
    
    for (nm in intersect(names(xtr), names(xte))) {
      if (is.factor(xtr[[nm]])) xte[[nm]] <- factor(xte[[nm]], levels = levels(xtr[[nm]]))
    }
    
    ytr <- as.integer(train$Event)
    pool_tr <- catboost::catboost.load_pool(xtr, label = ytr)
    pool_te <- catboost::catboost.load_pool(xte)
    
    fit <- catboost::catboost.train(
      learn_pool = pool_tr,
      params = list(
        loss_function = "Logloss",
        random_seed   = seed,
        iterations    = iterations,
        depth         = depth,
        learning_rate = learning_rate,
        l2_leaf_reg   = l2_leaf_reg,
        border_count  = border_count,
        verbose       = as.integer(verbose)
      )
    )
    
    as.numeric(catboost::catboost.predict(fit, pool_te, prediction_type = "Probability"))
  }
  
  res <- cv5_cindex(df, risk_fun, seed = seed, reverse = TRUE, K = K, fold_id = fold_id)
  message(sprintf("[CatBoost Bin→C] %s: done. CV%d C = %.3f", method_name, K, res$c_overall))
  summarize_cv(method_name, sprintf("CatBoost_Binomial_%dCV_toC", K), res)
}

# XGBoost Cox
xgb_cox_cindex_cv5_grid_tune_once <- function(
    df, method_name,
    seed = 1, K = 5, k_inner = 3,
    depth_grid     = c(3L, 4L),
    min_child_grid = c(4L, 8L),
    subsample_grid = c(0.70, 0.85, 1.00),
    colsample_grid = c(0.50, 0.70, 0.90),
    eta_grid       = c(0.05, 0.08, 0.12),
    lambda_grid    = c(0, 2, 5, 10),
    alpha_grid     = c(0, 0.10, 0.50, 1.00),
    nrounds_max = 2000, early_stopping_rounds = 50,
    fold_id = NULL, use_gpu = FALSE, nthread = NULL
) {
  if (!requireNamespace("xgboost", quietly = TRUE))
    stop("Please install.packages('xgboost')")
  if (!is.null(fold_id) && length(fold_id) != nrow(df))
    stop("'fold_id' length must equal nrow(df)")
  
  message(sprintf(
    "[XGB GRID-ONCE] %s: tuning once (k_inner=%d, grid=%d), then %d-fold CV",
    method_name, k_inner,
    length(depth_grid) * length(min_child_grid) * length(subsample_grid) *
      length(colsample_grid) * length(eta_grid) * length(lambda_grid) * length(alpha_grid),
    K
  ))
  
  x_df <- df[, setdiff(names(df), c("Event_time","Event")), drop = FALSE]
  MM   <- stats::model.matrix(~ . - 1, data = x_df)
  eps  <- .Machine$double.eps
  time <- pmax(df$Event_time, eps)
  y    <- ifelse(df$Event == 1, time, -time)
  
  set.seed(seed)
  inner_id <- make_stratified_folds(df$Event, K = k_inner, seed = seed)
  
  feval_cindex <- function(preds, dmat) {
    ylab <- xgboost::getinfo(dmat, "label")
    t    <- abs(ylab)
    e    <- as.integer(ylab > 0)
    cval <- harrell_c(t, e, preds, reverse = TRUE)
    list(metric = "C", value = cval)
  }
  
  tree_method      <- if (use_gpu) "gpu_hist" else "hist"
  sampling_method  <- "uniform"
  grow_policy      <- "depthwise"
  max_bin          <- 256
  colsample_bynode <- 1.0
  
  grid <- expand.grid(
    max_depth        = depth_grid,
    min_child_weight = min_child_grid,
    subsample        = subsample_grid,
    colsample_bytree = colsample_grid,
    eta              = eta_grid,
    lambda           = lambda_grid,
    alpha            = alpha_grid,
    KEEP.OUT.ATTRS = FALSE, stringsAsFactors = FALSE
  )
  
  best_mean_c <- -Inf
  best_pars   <- NULL
  best_iters  <- NULL
  
  for (gi in seq_len(nrow(grid))) {
    g <- grid[gi, ]
    c_fold <- numeric(k_inner)
    iters  <- integer(k_inner)
    
    for (ki in seq_len(k_inner)) {
      tr <- inner_id != ki
      va <- inner_id == ki
      
      dtr <- xgboost::xgb.DMatrix(MM[tr, , drop = FALSE], label = y[tr])
      dva <- xgboost::xgb.DMatrix(MM[va, , drop = FALSE], label = y[va])
      
      params <- list(
        booster          = "gbtree",
        objective        = "survival:cox",
        tree_method      = tree_method,
        sampling_method  = sampling_method,
        grow_policy      = grow_policy,
        max_bin          = max_bin,
        max_depth        = as.integer(g$max_depth),
        min_child_weight = as.integer(g$min_child_weight),
        subsample        = g$subsample,
        colsample_bytree = g$colsample_bytree,
        colsample_bynode = colsample_bynode,
        eta              = g$eta,
        lambda           = g$lambda,
        alpha            = g$alpha,
        verbosity        = 0
      )
      if (!is.null(nthread)) params$nthread <- as.integer(nthread)
      
      bst <- xgboost::xgb.train(
        params = params, data = dtr, nrounds = nrounds_max,
        watchlist = list(valid = dva),
        early_stopping_rounds = early_stopping_rounds,
        feval = feval_cindex, maximize = TRUE,
        verbose = 0
      )
      
      c_fold[ki] <- as.numeric(bst$best_score)
      bi <- bst$best_iteration; if (is.null(bi) || is.na(bi) || bi < 1) bi <- nrounds_max
      iters[ki] <- as.integer(bi)
    }
    
    mc <- mean(c_fold)
    if (mc > best_mean_c) {
      best_mean_c <- mc
      best_pars   <- g
      best_iters  <- iters
    }
  }
  
  chosen_nrounds <- as.integer(stats::median(best_iters))
  message(sprintf(
    "[XGB GRID-ONCE] %s: best -> depth=%d, min_child=%d, subs=%.2f, col_tree=%.2f, eta=%.3f, lambda=%.2f, alpha=%.2f, nrounds=%d (mean C=%.4f)",
    method_name,
    as.integer(best_pars$max_depth),
    as.integer(best_pars$min_child_weight),
    best_pars$subsample, best_pars$colsample_bytree,
    best_pars$eta, best_pars$lambda, best_pars$alpha,
    chosen_nrounds, best_mean_c
  ))
  
  risk_fun <- function(train, test) {
    xtr_df <- train[, setdiff(names(train), c("Event_time","Event")), drop = FALSE]
    xte_df <- test[,  setdiff(names(test),  c("Event_time","Event")), drop = FALSE]
    MM_all <- stats::model.matrix(~ . - 1, data = rbind(xtr_df, xte_df))
    ntr    <- nrow(xtr_df)
    Xtr    <- MM_all[seq_len(ntr), , drop = FALSE]
    Xte    <- MM_all[(ntr + 1L):nrow(MM_all), , drop = FALSE]
    
    eps <- .Machine$double.eps
    ttr <- pmax(train$Event_time, eps)
    ytr <- ifelse(train$Event == 1, ttr, -ttr)
    
    dtr <- xgboost::xgb.DMatrix(Xtr, label = ytr)
    dte <- xgboost::xgb.DMatrix(Xte)
    
    params_best <- list(
      booster          = "gbtree",
      objective        = "survival:cox",
      tree_method      = tree_method,
      sampling_method  = sampling_method,
      grow_policy      = grow_policy,
      max_bin          = max_bin,
      max_depth        = as.integer(best_pars$max_depth),
      min_child_weight = as.integer(best_pars$min_child_weight),
      subsample        = best_pars$subsample,
      colsample_bytree = best_pars$colsample_bytree,
      colsample_bynode = colsample_bynode,
      eta              = best_pars$eta,
      lambda           = best_pars$lambda,
      alpha            = best_pars$alpha,
      verbosity        = 0
    )
    if (!is.null(nthread)) params_best$nthread <- as.integer(nthread)
    
    bst <- xgboost::xgb.train(
      params = params_best, data = dtr, nrounds = chosen_nrounds, verbose = 0
    )
    as.numeric(predict(bst, dte))
  }
  
  res <- cv5_cindex(df, risk_fun, seed = seed, reverse = TRUE, K = K, fold_id = fold_id)
  summarize_cv(method_name, sprintf("XGB_Cox_%dCV_Grid_TunedOnce", K), res)
}


