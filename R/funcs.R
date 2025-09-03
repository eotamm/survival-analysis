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



# C-INDEX
harrell_c <- function(time, event, score, reverse = TRUE) {
  ok <- is.finite(time) & is.finite(event) & is.finite(score)
  ok[is.na(ok)] <- FALSE
  if (!any(ok)) return(NA_real_)
  as.numeric(
    survival::concordance(
      survival::Surv(time[ok], event[ok]) ~ score[ok],
      reverse = reverse
    )$concordance
  )
}
# Time-dependent AUC 
time_auc_at <- function(time, event, marker, tau) {
  keep <- is.finite(time) & is.finite(event) & is.finite(marker)
  if (!any(keep) || !is.finite(tau)) return(NA_real_)
  t <- time[keep]; d <- event[keep]; m <- marker[keep]
  
  if (!(any(d == 1 & t <= tau) && any(t > tau))) return(NA_real_)
  if (stats::sd(m, na.rm = TRUE) <= 0) return(NA_real_)
  
  tr <- tryCatch(
    timeROC::timeROC(T = t, delta = d, marker = m, cause = 1, times = tau, iid = FALSE),
    error = function(e) NULL
  )
  if (is.null(tr)) return(NA_real_)
  as.numeric(tr$AUC[which.min(abs(tr$times - tau))])
}

# Ggeneric OOB bootstrap
oob_bootstrap <- function(df, B = 500,
                          risk_fun,
                          require_train_events = TRUE,
                          require_oob_events   = FALSE,
                          verbose = TRUE) {
  df_cc <- df[stats::complete.cases(df), , drop = FALSE]
  n <- nrow(df_cc); if (n < 5) stop("Too few observations for OOB bootstrap.")
  tau <- stats::median(df_cc$Event_time[df_cc$Event == 1], na.rm = TRUE)
  
  oob_c   <- rep(NA_real_, B)
  oob_auc <- rep(NA_real_, B)
  
  if (isTRUE(verbose)) message(sprintf("[OOB] B=%d; tau=%.4f", B, tau))
  
  for (b in seq_len(B)) {
    # Resample until OOB non-empty
    idx <- sample.int(n, replace = TRUE)
    oob <- setdiff(seq_len(n), unique(idx))
    while (length(oob) == 0L) { idx <- sample.int(n, replace = TRUE); oob <- setdiff(seq_len(n), unique(idx)) }
    
    train <- df_cc[idx, , drop = FALSE]
    test  <- df_cc[oob, , drop = FALSE]
    
    # Event guards
    if (require_train_events && sum(train$Event == 1, na.rm = TRUE) == 0L) next
    if (require_oob_events   && sum(test$Event  == 1, na.rm = TRUE) == 0L) next
    
    # Fit model + get test risk
    risk <- tryCatch(risk_fun(train, test), error = function(e) rep(NA_real_, nrow(test)))
    risk <- as.numeric(risk)
    
    # C-index
    oob_c[b] <- harrell_c(test$Event_time, test$Event, risk, reverse = TRUE)
    
    # AUC
    oob_auc[b] <- time_auc_at(test$Event_time, test$Event, risk, tau)
  }
  
  list(tau = tau, c = oob_c, auc = oob_auc)
}

# Summarize bootstrap vectors
summarize_boot <- function(model_name, method_name, tau, c_vec, auc_vec) {
  c_has   <- any(is.finite(c_vec))
  auc_has <- any(is.finite(auc_vec))
  
  c_est <- if (c_has) stats::median(c_vec, na.rm = TRUE) else NA_real_
  c_lo  <- if (c_has) as.numeric(stats::quantile(c_vec,  0.025, na.rm = TRUE)) else NA_real_
  c_hi  <- if (c_has) as.numeric(stats::quantile(c_vec,  0.975, na.rm = TRUE)) else NA_real_
  
  auc_est <- if (auc_has) stats::median(auc_vec, na.rm = TRUE) else NA_real_
  auc_lo  <- if (auc_has) as.numeric(stats::quantile(auc_vec, 0.025, na.rm = TRUE)) else NA_real_
  auc_hi  <- if (auc_has) as.numeric(stats::quantile(auc_vec, 0.975, na.rm = TRUE)) else NA_real_
  
  tibble::tibble(
    model    = model_name,
    method   = method_name,
    metric   = c("C", "AUC"),
    time     = c(NA_real_, tau),
    estimate = c(c_est, auc_est),
    lower    = c(c_lo,  auc_lo),
    upper    = c(c_hi,  auc_hi)
  )
}

# RSF: OOB bootstrap
rsf_cindex_boot_oob <- function(
    df, method_name,
    B = 500,
    num.trees = 1500,
    min.node.size = 10,
    mtry = NULL
) {
  risk_fun <- function(train, test) {
    p_train <- max(1, ncol(train) - 2L)
    mtry_use <- if (is.null(mtry)) max(1, min(p_train, floor(sqrt(p_train)))) else mtry
    fit <- ranger::ranger(
      survival::Surv(Event_time, Event) ~ .,
      data = train,
      num.trees = num.trees,
      mtry = mtry_use,
      min.node.size = min.node.size,
      splitrule = "logrank",
      write.forest = TRUE
    )
    pr <- predict(fit, data = test)
    chf <- pr$chf
    if (is.matrix(chf)) chf[, ncol(chf), drop = TRUE] else as.numeric(chf)
  }
  message(sprintf("[RSF %s] B=%d", method_name, B))
  res <- oob_bootstrap(df, B = B, risk_fun = risk_fun,
                       require_train_events = TRUE, require_oob_events = FALSE, verbose = FALSE)
  tbl <- summarize_boot(model_name = method_name, method_name = "RSF_OOB",
                        tau = res$tau, c_vec = res$c, auc_vec = res$auc)
}


# XGBoost: OOB bootstrap
xgb_cox_cindex_boot_oob <- function(
    df, method_name,
    B = 500,
    nrounds = 500,
    max_depth = 3, eta = 0.04,
    subsample = 0.65, colsample_bytree = NULL,
    min_child_weight = 3, reg_lambda = 1, reg_alpha = 0
) {
  if (!requireNamespace("xgboost", quietly = TRUE))
    stop("Please install.packages('xgboost')")
  # risk_fun: fit XGB-Cox on 'train', return risk scores for 'test'
  risk_fun <- function(train, test) {
    xtr_df <- train[, setdiff(names(train), c("Event_time", "Event")), drop = FALSE]
    xte_df <- test[,  setdiff(names(test),  c("Event_time", "Event")), drop = FALSE]
    mm_all <- stats::model.matrix(~ . - 1, data = rbind(xtr_df, xte_df))
    ntr    <- nrow(xtr_df)
    Xtr    <- mm_all[seq_len(ntr), , drop = FALSE]
    Xte    <- mm_all[(ntr + 1L):nrow(mm_all), , drop = FALSE]
    p <- ncol(Xtr); if (p == 0L) return(rep(NA_real_, nrow(test)))
    colsample_use <- if (is.null(colsample_bytree)) sqrt(p) / p else colsample_bytree
    
    # XGB Cox labels: +time for events, -time for censored
    eps <- .Machine$double.eps
    ttr <- pmax(train$Event_time, eps)
    ytr <- ifelse(train$Event == 1, ttr, -ttr)
    dtr <- xgboost::xgb.DMatrix(data = Xtr, label = ytr)
    dte <- xgboost::xgb.DMatrix(data = Xte)
    params <- list(objective="survival:cox", eval_metric="cox-nloglik",
                   max_depth=max_depth, eta=eta, subsample=subsample,
                   colsample_bytree=colsample_use, min_child_weight=min_child_weight,
                   lambda=reg_lambda, alpha=reg_alpha)
    fit <- xgboost::xgb.train(params=params, data=dtr, nrounds=nrounds, verbose=0)
    as.numeric(predict(fit, dte))
  }
  message(sprintf("[XGB%s] B=%d, nrounds=%d", method_name, B, nrounds))
  res <- oob_bootstrap(df, B = B, risk_fun = risk_fun,
                       require_train_events = TRUE, require_oob_events = FALSE, verbose = FALSE)
  tbl <- summarize_boot(model_name = method_name, method_name = "XGB_Cox_OOB",
                        tau = res$tau, c_vec = res$c, auc_vec = res$auc)
}


# DeepSurv: OOB bootstrap
deepsurv_cindex_boot_oob <- function(
    df, method_name,
    B = 500,
    hidden = c(64, 32),
    dropout = 0.0,
    l2 = 1e-4,
    lr = 1e-3,
    epochs = 300, patience = 30,
    verbose = 0,
    run_eagerly = TRUE,
    require_oob_events = TRUE
) {
  
  # risk_fun: fit DeepSurv on 'train', return risk scores for 'test'
  risk_fun <- function(train, test) {
    xtr_df <- train[, setdiff(names(train), c("Event_time","Event")), drop = FALSE]
    xte_df <- test[,  setdiff(names(test),  c("Event_time","Event")), drop = FALSE]
    mm_all <- stats::model.matrix(~ . - 1, data = rbind(xtr_df, xte_df))
    ntr    <- nrow(xtr_df)
    Xtr    <- mm_all[seq_len(ntr), , drop = FALSE]
    Xte    <- mm_all[(ntr + 1L):nrow(mm_all), , drop = FALSE]
    p      <- ncol(Xtr); if (p == 0L) return(rep(NA_real_, nrow(test)))
    
    # Standardize using training statistics
    mu <- matrixStats::colMeans2(Xtr)
    sd <- matrixStats::colSds(Xtr); sd[!is.finite(sd) | sd < 1e-8] <- 1
    Xtr_s <- sweep(sweep(Xtr, 2, mu, "-"), 2, sd, "/")
    Xte_s <- sweep(sweep(Xte, 2, mu, "-"), 2, sd, "/")
    
    # Sort training rows by time for stable partial-likelihood risk sets
    ord   <- order(train$Event_time)
    Xtr_s <- Xtr_s[ord, , drop = FALSE]
    ev_tr <- as.numeric(train$Event[ord])
    
    # Cox partial likelihood loss
    cox_ph_loss_safe <- function(y_true, y_pred) {
      tf <- tensorflow::tf
      y_pred <- tf$reshape(y_pred, shape = c(-1L))
      events <- tf$reshape(y_true, shape = c(-1L))
      haz <- tf$math$exp(y_pred)
      rev_csum  <- tf$math$cumsum(tf$reverse(haz, list(0L)))
      risk_csum <- tf$reverse(rev_csum, list(0L))
      log_risk  <- tf$math$log(risk_csum + 1e-8)
      -tf$math$reduce_sum((y_pred - log_risk) * events) / (tf$math$reduce_sum(events) + 1e-8)
    }
    # Network
    inp <- keras::layer_input(shape = p, dtype = "float32")
    x <- inp
    for (u in hidden) {
      x <- keras::layer_dense(x, units = u, activation = "relu",
                              kernel_regularizer = keras::regularizer_l2(l = l2))
      if (dropout > 0) x <- keras::layer_dropout(x, rate = dropout)
    }
    out <- keras::layer_dense(x, units = 1, activation = "linear",
                              kernel_regularizer = keras::regularizer_l2(l = l2))
    model <- keras::keras_model(inp, out)
    model %>% keras::compile(
      optimizer = keras::optimizer_adam(learning_rate = lr, clipnorm = 1.0, clipvalue = 0.5),
      loss = cox_ph_loss_safe,
      run_eagerly = run_eagerly
    )
    # Fit
    ok <- TRUE
    tryCatch({
      model %>% keras::fit(
        x = Xtr_s, y = matrix(ev_tr, ncol = 1),
        batch_size = nrow(Xtr_s), epochs = epochs,
        shuffle = FALSE, verbose = verbose,
        callbacks = list(keras::callback_early_stopping(monitor = "loss",
                                                        patience = patience, restore_best_weights = TRUE))
      )
    }, error = function(e) ok <<- FALSE)
    if (!ok) return(rep(NA_real_, nrow(test)))
    # OOB risk scores
    risk <- as.numeric(model$predict(Xte_s, verbose = as.integer(verbose)))
    if (!all(is.finite(risk))) {
      if (any(is.finite(risk))) {
        med <- stats::median(risk[is.finite(risk)], na.rm = TRUE)
        risk[!is.finite(risk)] <- med
      } else {
        risk <- rep(NA_real_, length(risk))
      }
    }
    risk
  }
  message(sprintf("[DeepSurv %s] B=%d", method_name, B))
  res <- oob_bootstrap(df, B = B, risk_fun = risk_fun,
                       require_train_events = TRUE, require_oob_events = require_oob_events, verbose = FALSE)
  tbl <- summarize_boot(model_name = method_name, method_name = "DeepSurv_OOB",
                        tau = res$tau, c_vec = res$c, auc_vec = res$auc)
}



# Logistic (censoring ignored): OOB bootstrap
logit_cindex_boot_oob <- function(
    df, method_name,
    B = 500,
    class_weights = TRUE,
    maxit = 200
) {
  
  # risk_fun: fit logistic regression on 'train' (Event ~ X), predict P(Event=1) for 'test'
  risk_fun <- function(train, test) {
    xtr_df <- train[, setdiff(names(train), c("Event_time","Event")), drop = FALSE]
    xte_df <- test[,  setdiff(names(test),  c("Event_time","Event")), drop = FALSE]
    MM     <- stats::model.matrix(~ . - 1, data = rbind(xtr_df, xte_df))
    ntr    <- nrow(xtr_df)
    Xtr    <- MM[seq_len(ntr), , drop = FALSE]
    Xte    <- MM[(ntr + 1L):nrow(MM), , drop = FALSE]
    ytr    <- as.integer(train$Event)
    if (length(unique(ytr)) < 2L) return(rep(mean(ytr), nrow(test)))
    keep <- apply(Xtr, 2, function(z) stats::sd(z, na.rm = TRUE) > 0)
    if (!any(keep)) return(rep(mean(ytr), nrow(test)))
    Xtr <- Xtr[, keep, drop = FALSE]
    Xte <- Xte[, keep, drop = FALSE]
    w <- NULL
    if (isTRUE(class_weights)) {
      n_pos <- sum(ytr == 1L); n_neg <- sum(ytr == 0L)
      if (n_pos > 0 && n_neg > 0) w <- ifelse(ytr == 1L, n_neg / n_pos, 1)
    }
    
    # Fit logistic regression
    dtr <- data.frame(y = ytr, Xtr, check.names = FALSE)
    fit <- tryCatch(stats::glm(y ~ ., data = dtr, family = stats::binomial(), weights = w,
                               control = list(maxit = maxit)),
                    error = function(e) NULL)
    if (is.null(fit)) return(rep(mean(ytr), nrow(test)))
    pr <- tryCatch(as.numeric(stats::predict(fit, newdata = data.frame(Xte, check.names = FALSE), type = "response")),
                   error = function(e) rep(mean(ytr), nrow(test)))
    if (!all(is.finite(pr))) {
      if (any(is.finite(pr))) {
        med <- stats::median(pr[is.finite(pr)], na.rm = TRUE)
        pr[!is.finite(pr)] <- med
      } else {
        pr <- rep(mean(ytr), length(pr))
      }
    }
    pr
  }
  message(sprintf("[Logit %s] B=%d", method_name, B))
  res <- oob_bootstrap(df, B = B, risk_fun = risk_fun,
                       require_train_events = FALSE, require_oob_events = FALSE, verbose = FALSE)
  tbl <- summarize_boot(model_name = method_name, method_name = "Logit_OOB",
                        tau = res$tau, c_vec = res$c, auc_vec = res$auc)
}


# CoxPH: OOB bootstrap
coxph_cindex_boot_oob <- function(
    df, method_name,
    B = 500,
    ties = c("efron","breslow","exact"),
    iter_max = 50,
    quiet = TRUE
) {
  ties <- match.arg(ties)
  
  # risk_fun: fit Cox on 'train', return linear predictor for 'test'
  risk_fun <- function(train, test) {
    pred_names <- setdiff(names(train), c("Event_time","Event"))
    keep <- vapply(pred_names, function(nm) stats::sd(train[[nm]], na.rm = TRUE) > 0, logical(1))
    preds_kept <- pred_names[keep]
    if (length(preds_kept) == 0L) return(rep(NA_real_, nrow(test)))
    fml <- stats::as.formula(
      paste0("survival::Surv(Event_time, Event) ~ ", paste(preds_kept, collapse = " + "))
    )
    fit_call <- function() survival::coxph(
      fml, data = train, ties = ties,
      control = survival::coxph.control(iter.max = iter_max),
      x = FALSE, y = FALSE
    )
    fit <- tryCatch(if (quiet) suppressWarnings(fit_call()) else fit_call(),
                    error = function(e) NULL)
    if (is.null(fit)) return(rep(NA_real_, nrow(test)))
    
    # LP on OOB (replace non-finite with finite median if needed)
    lp <- tryCatch(as.numeric(stats::predict(fit, newdata = test, type = "lp")),
                   error = function(e) rep(NA_real_, nrow(test)))
    if (!all(is.finite(lp))) {
      if (any(is.finite(lp))) {
        med <- stats::median(lp[is.finite(lp)], na.rm = TRUE)
        lp[!is.finite(lp)] <- med
      }
    }
    lp
  }
  message(sprintf("[CoxPH %s] B=%d, ties=%s", method_name, B, ties))
  res <- oob_bootstrap(df, B = B, risk_fun = risk_fun,
                       require_train_events = TRUE, require_oob_events = FALSE, verbose = FALSE)
  tbl <- summarize_boot(model_name = method_name, method_name = "CoxPH_OOB",
                        tau = res$tau, c_vec = res$c, auc_vec = res$auc)
}


