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


# Prepare data for analysis
extract_all_features_by_transformation <- function(
    tse,
    method = c("clr","rclr","log_abund","lra","pa","tss","logtss","asin","alr"),
    assay_in = "counts",         # set to your assay
    pseudocount = 1e-6,
    alr_ref = NULL,              # required when method == "alr"
    event_col = "Event",
    time_col  = "Event_time"
) {
  method <- match.arg(method)
  
  # Validate inputs
  if (!assay_in %in% SummarizedExperiment::assayNames(tse)) {
    stop(sprintf("Assay '%s' not found. Set 'assay_in' correctly.", assay_in))
  }
  if (method == "alr" && (is.null(alr_ref) || !nzchar(alr_ref))) {
    stop("ALR requires 'alr_ref', e.g. 'g_Turicibacter'.")
  }
  
  # Transform
  tse_tx <- switch(
    method,
    clr  = mia::transformAssay(tse, method = "clr",  assay.type = assay_in,
                               pseudocount = pseudocount, name = "clr"),
    rclr = mia::transformAssay(tse, method = "rclr", assay.type = assay_in,
                               pseudocount = pseudocount, name = "rclr"),
    log_abund = tse |>
      mia::transformAssay(method = "relabundance", assay.type = assay_in, name = "relabundance") |>
      mia::transformAssay(method = "log", assay.type = "relabundance",
                          pseudocount = pseudocount, name = "log_abund"),
    lra = tse |>
      mia::transformAssay(method = "relabundance", assay.type = assay_in, name = "relabundance") |>
      mia::transformAssay(method = "log", assay.type = "relabundance",
                          pseudocount = pseudocount, name = "log_abund") |>
      mia::transformAssay(method = "difference", assay.type = "log_abund",
                          name = "logratios", MARGIN = 1L),
    pa     = mia::transformAssay(tse, method = "pa",  assay.type = assay_in, name = "pa"),
    tss    = mia::transformAssay(tse, method = "relabundance", assay.type = assay_in, name = "tss"),
    logtss = tse |>
      mia::transformAssay(method = "relabundance", assay.type = assay_in, name = "tss") |>
      mia::transformAssay(method = "log", assay.type = "tss",
                          pseudocount = pseudocount, name = "logtss"),
    asin = mia::transformAssay(tse, method = "relabundance", assay.type = assay_in, name = "tss"),
    alr = mia::transformAssay(tse, method = "alr", assay.type = assay_in,
                              pseudocount = pseudocount, name = "alr", ref_taxa = alr_ref)
  )
  
  # Extract matrix
  M <- if (method == "lra") {
    assay(altExp(tse_tx, "logratios"))
  } else if (method == "asin") {
    asin(sqrt(assay(tse_tx, "tss")))
  } else {
    assay(tse_tx, method)
  }
  
  # Build feature frame
  df <- as.data.frame(t(as.matrix(M)))
  rownames(df) <- colnames(M)
  
  # Clean names
  df <- clean_column_names(df)
  
  # Drop ALR reference (clean with your cleaner)
  if (method == "alr") {
    refs <- unique(c(alr_ref,
                     sub("__","_", alr_ref, fixed = TRUE),
                     sub("^g__","g_", alr_ref),
                     make.names(alr_ref)))
    ref_clean <- names(
      clean_column_names(
        setNames(as.data.frame(matrix(NA, nrow = 1, ncol = length(refs))), refs)
      )
    )
    df <- df[, !(colnames(df) %in% ref_clean), drop = FALSE]
  }
  
  # Add metadata
  cd <- SummarizedExperiment::colData(tse_tx)
  if (!event_col %in% colnames(cd)) stop(sprintf("Missing '%s' in colData.", event_col))
  if (!time_col  %in% colnames(cd)) stop(sprintf("Missing '%s' in colData.", time_col))
  df$Event      <- as.vector(cd[[event_col]])
  df$Event_time <- as.vector(cd[[time_col]])
  
  df
}


# Compute Harrell's C-index
harrell_c <- function(time, event, score, reverse = TRUE) {
  as.numeric(
    survival::concordance(
      survival::Surv(time, event) ~ score,
      reverse = reverse
    )$concordance
  )
}

# Stratified K-fold indices on binary Event (0/1)
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

# Generic K-fold CV (model returns OOF risk scores)
cv5_cindex <- function(df, risk_fun, seed = 1, reverse = TRUE, K = 5, fold_id = NULL,
                       event_col = "Event", time_col = "Event_time") {
  n <- nrow(df)
  if (!is.null(fold_id) && length(fold_id) != n)
    stop("'fold_id' length must equal nrow(df)")
  if (is.null(fold_id))
    fold_id <- make_stratified_folds(df[[event_col]], K = K, seed = seed)
  
  preds <- numeric(n)
  for (k in seq_len(K)) {
    train <- df[fold_id != k, , drop = FALSE]
    test  <- df[fold_id == k, , drop = FALSE]
    preds[fold_id == k] <- as.numeric(risk_fun(train, test))
  }
  
  c_overall <- harrell_c(df[[time_col]], df[[event_col]], preds, reverse = reverse)
  c_folds <- vapply(
    seq_len(K),
    function(k) harrell_c(
      df[[time_col]][fold_id == k],
      df[[event_col]][fold_id == k],
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

# Summarize CV to tibble (normal-approx CI)
summarize_cv <- function(model_name, method_name, cv_res, level = 0.95,
                         event_col = "Event", time_col = "Event_time") {
  cc <- survival::concordance(
    survival::Surv(cv_res$data[[time_col]], cv_res$data[[event_col]]) ~ cv_res$preds,
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


# Random Survival Forest K-fold CV
rsf_cindex_cv5 <- function(
    df, method_name,
    num.trees = 1500, min.node.size = 10, mtry = NULL,
    seed = 1, K = 5, fold_id = NULL,
    event_col = "Event", time_col = "Event_time"
) {
  message(sprintf("[RSF] %s: starting %d-fold CV", method_name, K))
  
  risk_fun <- function(train, test) {
    # build Surv() ~ . with given column names
    fml <- stats::as.formula(sprintf("survival::Surv(%s, %s) ~ .", time_col, event_col))
    # ranger needs factors, not characters
    train[] <- lapply(train, function(x) if (is.character(x)) factor(x) else x)
    test[]  <- lapply(test,  function(x) if (is.character(x)) factor(x) else x)
    
    p_train  <- max(1L, ncol(train) - 2L)
    mtry_use <- if (is.null(mtry)) max(1L, min(p_train, floor(sqrt(p_train)))) else as.integer(mtry)
    
    fit <- ranger::ranger(
      formula = fml,
      data = train,
      num.trees = num.trees,
      mtry = mtry_use,
      min.node.size = min.node.size,
      splitrule = "logrank",
      write.forest = TRUE,
      seed = seed
    )
    chf <- predict(fit, data = test)$chf
    if (is.matrix(chf)) chf[, ncol(chf), drop = TRUE] else as.numeric(chf)
  }
  
  res <- cv5_cindex(df, risk_fun, seed = seed, reverse = TRUE, K = K,
                    fold_id = fold_id, event_col = event_col, time_col = time_col)
  
  message(sprintf("[RSF] %s: done. CV%d C = %.3f", method_name, K, res$c_overall))
  summarize_cv(method_name, sprintf("RSF"), res,
               event_col = event_col, time_col = time_col)
}

# Logistic regression (binomial) K-fold CV
logit_cindex_cv5 <- function(
    df, method_name, maxit = 200,
    seed = 1, K = 5, fold_id = NULL,
    event_col = "Event", time_col = "Event_time"
) {
  message(sprintf("[Logit] %s: starting %d-fold CV", method_name, K))
  
  risk_fun <- function(train, test) {
    # drop survival cols
    xtr <- train[, setdiff(names(train), c(time_col, event_col)), drop = FALSE]
    xte <- test[,  setdiff(names(test),  c(time_col, event_col)), drop = FALSE]
    # handle characters as factors
    xtr[] <- lapply(xtr, function(z) if (is.character(z)) factor(z) else z)
    xte[] <- lapply(xte, function(z) if (is.character(z)) factor(z) else z)
    # shared design
    MM  <- stats::model.matrix(~ . - 1, data = rbind(xtr, xte))
    ntr <- nrow(xtr)
    Xtr <- MM[seq_len(ntr), , drop = FALSE]
    Xte <- MM[(ntr + 1L):nrow(MM), , drop = FALSE]
    ytr <- as.integer(train[[event_col]])
    
    fit <- stats::glm(
      ytr ~ .,
      data    = data.frame(Xtr, check.names = FALSE),
      family  = stats::binomial(),
      control = list(maxit = maxit)
    )
    as.numeric(stats::predict(fit, newdata = data.frame(Xte, check.names = FALSE), type = "response"))
  }
  
  res <- cv5_cindex(df, risk_fun, seed = seed, reverse = TRUE, K = K,
                    fold_id = fold_id, event_col = event_col, time_col = time_col)
  
  message(sprintf("[Logit] %s: done. CV%d C = %.3f", method_name, K, res$c_overall))
  summarize_cv(method_name, sprintf("Logit"), res,
               event_col = event_col, time_col = time_col)
}


# Cox proportional hazards K-fold CV
coxph_cindex_cv5 <- function(
    df, method_name,
    ties = c("efron","breslow","exact"),
    iter_max = 50,
    seed = 1, K = 5, fold_id = NULL,
    event_col = "Event", time_col = "Event_time"
) {
  ties <- match.arg(ties)
  message(sprintf("[CoxPH] %s: starting %d-fold CV (ties=%s)", method_name, K, ties))
  
  risk_fun <- function(train, test) {
    # build Surv() ~ . with given column names
    fml <- stats::as.formula(sprintf("survival::Surv(%s, %s) ~ .", time_col, event_col))
    # ensure characters are factors
    train[] <- lapply(train, function(x) if (is.character(x)) factor(x) else x)
    test[]  <- lapply(test,  function(x) if (is.character(x)) factor(x) else x)
    
    fit <- survival::coxph(
      formula = fml,
      data    = train,
      ties    = ties,
      control = survival::coxph.control(iter.max = iter_max),
      x = FALSE, y = FALSE
    )
    as.numeric(stats::predict(fit, newdata = test, type = "lp"))
  }
  
  res <- cv5_cindex(df, risk_fun, seed = seed, reverse = TRUE, K = K,
                    fold_id = fold_id, event_col = event_col, time_col = time_col)
  
  message(sprintf("[CoxPH] %s: done. CV%d C = %.3f", method_name, K, res$c_overall))
  summarize_cv(method_name, sprintf("CoxPH"), res,
               event_col = event_col, time_col = time_col)
}


# DeepSurv, K-fold CV
deepsurv_cindex_cv5 <- function(
    df, method_name,
    hidden = c(64, 32), dropout = 0.0, l2 = 1e-4,
    lr = 1e-3, epochs = 300, patience = 30,
    verbose = 0, run_eagerly = TRUE,
    seed = 1, K = 5, fold_id = NULL,
    event_col = "Event", time_col = "Event_time"
) {
  message(sprintf("[DeepSurv] %s: starting %d-fold CV", method_name, K))
  
  risk_fun <- function(train, test) {
    # features only
    xtr <- train[, setdiff(names(train), c(time_col, event_col)), drop = FALSE]
    xte <- test[,  setdiff(names(test),  c(time_col, event_col)), drop = FALSE]
    # design matrix (shared)
    MM  <- stats::model.matrix(~ . - 1, data = rbind(xtr, xte))
    ntr <- nrow(xtr)
    Xtr <- MM[seq_len(ntr), , drop = FALSE]
    Xte <- MM[(ntr + 1L):nrow(MM), , drop = FALSE]
    p   <- ncol(Xtr)
    
    # standardize by train stats
    mu <- matrixStats::colMeans2(Xtr)
    sd <- matrixStats::colSds(Xtr); sd[!is.finite(sd) | sd < 1e-8] <- 1
    Xtr_s <- sweep(sweep(Xtr, 2, mu, "-"), 2, sd, "/")
    Xte_s <- sweep(sweep(Xte, 2, mu, "-"), 2, sd, "/")
    
    # sort by time (ascending)
    ord   <- order(train[[time_col]])
    Xtr_s <- Xtr_s[ord, , drop = FALSE]
    ev_tr <- as.numeric(train[[event_col]][ord])
    
    # Cox partial likelihood loss
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
    
    # model
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
    model %>% keras::fit(
      x = Xtr_s, y = matrix(ev_tr, ncol = 1),
      batch_size = nrow(Xtr_s), epochs = epochs,
      shuffle = FALSE, verbose = verbose,
      callbacks = list(keras::callback_early_stopping(
        monitor = "loss", patience = patience, restore_best_weights = TRUE))
    )
    as.numeric(model$predict(Xte_s, verbose = as.integer(verbose)))
  }
  
  res <- cv5_cindex(df, risk_fun, seed = seed, reverse = TRUE, K = K,
                    fold_id = fold_id, event_col = event_col, time_col = time_col)
  
  message(sprintf("[DeepSurv] %s: done. CV%d C = %.3f", method_name, K, res$c_overall))
  summarize_cv(method_name, sprintf("DeepSurv"), res,
               event_col = event_col, time_col = time_col)
}


# CatBoost (binomial)
catboost_bin_cindex_cv5 <- function(
    df, method_name,
    seed = 1, K = 5, fold_id = NULL,
    iterations = 1000, depth = 6, learning_rate = 0.05,
    l2_leaf_reg = 3.0, border_count = 254, verbose = 0,
    event_col = "Event", time_col = "Event_time"
) {
  message(sprintf("[CatBoost Bin] %s: starting %d-fold CV", method_name, K))
  
  risk_fun <- function(train, test) {
    # features only
    xtr <- train[, setdiff(names(train), c(time_col, event_col)), drop = FALSE]
    xte <- test[,  setdiff(names(test),  c(time_col, event_col)), drop = FALSE]
    
    # catboost: characters/logicals as factors; align levels
    xtr[] <- lapply(xtr, function(v) if (is.character(v) || is.logical(v)) factor(v) else v)
    xte[] <- lapply(xte, function(v) if (is.character(v) || is.logical(v)) factor(v) else v)
    for (nm in intersect(names(xtr), names(xte))) {
      if (is.factor(xtr[[nm]])) xte[[nm]] <- factor(xte[[nm]], levels = levels(xtr[[nm]]))
    }
    
    ytr <- as.integer(train[[event_col]])
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
    
    # probability of event (higher = higher risk)
    as.numeric(catboost::catboost.predict(fit, pool_te, prediction_type = "Probability"))
  }
  
  res <- cv5_cindex(df, risk_fun, seed = seed, reverse = TRUE, K = K,
                    fold_id = fold_id, event_col = event_col, time_col = time_col)
  
  message(sprintf("[CatBoost Bin] %s: done. CV%d C = %.3f", method_name, K, res$c_overall))
  summarize_cv(method_name, sprintf("CatBoost"), res,
               event_col = event_col, time_col = time_col)
}

# XGBoost Cox: tune once with inner CV, then K-fold CV (print only outer C)
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
    fold_id = NULL, use_gpu = FALSE, nthread = NULL,
    event_col = "Event", time_col = "Event_time"
) {
  message(sprintf("[XGB] %s: starting %d-fold CV", method_name, K))
  
  # design for inner tuning
  x_df <- df[, setdiff(names(df), c(time_col, event_col)), drop = FALSE]
  MM   <- stats::model.matrix(~ . - 1, data = x_df)
  eps  <- .Machine$double.eps
  time <- pmax(df[[time_col]], eps)
  y    <- ifelse(df[[event_col]] == 1, time, -time)
  
  set.seed(seed)
  inner_id <- make_stratified_folds(df[[event_col]], K = k_inner, seed = seed)
  
  feval_cindex <- function(preds, dmat) {
    ylab <- xgboost::getinfo(dmat, "label")
    t    <- abs(ylab); e <- as.integer(ylab > 0)
    list(metric = "C", value = harrell_c(t, e, preds, reverse = TRUE))
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
      tr <- inner_id != ki; va <- inner_id == ki
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
  
  risk_fun <- function(train, test) {
    xtr_df <- train[, setdiff(names(train), c(time_col, event_col)), drop = FALSE]
    xte_df <- test[,  setdiff(names(test),  c(time_col, event_col)), drop = FALSE]
    MM_all <- stats::model.matrix(~ . - 1, data = rbind(xtr_df, xte_df))
    ntr    <- nrow(xtr_df)
    Xtr    <- MM_all[seq_len(ntr), , drop = FALSE]
    Xte    <- MM_all[(ntr + 1L):nrow(MM_all), , drop = FALSE]
    
    eps <- .Machine$double.eps
    ttr <- pmax(train[[time_col]], eps)
    ytr <- ifelse(train[[event_col]] == 1, ttr, -ttr)
    
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
    
    bst <- xgboost::xgb.train(params = params_best, data = dtr,
                              nrounds = chosen_nrounds, verbose = 0)
    as.numeric(predict(bst, dte))
  }
  
  res <- cv5_cindex(df, risk_fun, seed = seed, reverse = TRUE, K = K,
                    fold_id = fold_id, event_col = event_col, time_col = time_col)
  
  message(sprintf("[XGB] %s: done. CV%d C = %.3f", method_name, K, res$c_overall))
  summarize_cv(method_name, sprintf("XGB_Cox"), res,
               event_col = event_col, time_col = time_col)
}


# TabPFN: init + wrappers
tabpfn_init <- function(conda_env = "r-tabpfn") {
  reticulate::use_condaenv(conda_env, required = TRUE)
  if (!reticulate::py_module_available("tabpfn"))
    stop("Python module 'tabpfn' not found in conda env '", conda_env, "'.")
  reticulate::py_run_string("
import importlib
def _make_tabpfn_clf(device='cpu', ensemble=16, seed=1):
    tabpfn = importlib.import_module('tabpfn')
    clf = tabpfn.TabPFNClassifier(device=device, random_state=int(seed))
    if hasattr(clf, 'set_ensemble_size'):
        clf.set_ensemble_size(int(ensemble))
    elif hasattr(clf, 'N_ensemble_configurations'):
        clf.N_ensemble_configurations = int(ensemble)
    elif hasattr(clf, 'n_ensemble_configurations'):
        clf.n_ensemble_configurations = int(ensemble)
    if hasattr(clf, 'ignore_pretraining_limits'):
        clf.ignore_pretraining_limits = True
    return clf
")
  invisible(TRUE)
}

one_hot_mm <- function(df) {
  X <- stats::model.matrix(~ . - 1, data = df)
  storage.mode(X) <- "double"
  X
}

tabpfn_classify <- function(X_train_df, y_train, X_test_df,
                            device = "cpu", ensemble = 32, seed = 1) {
  Xtr <- one_hot_mm(X_train_df)
  Xte <- one_hot_mm(X_test_df)
  miss <- setdiff(colnames(Xtr), colnames(Xte))
  if (length(miss))
    Xte <- cbind(Xte, matrix(0, nrow(Xte), length(miss), dimnames = list(NULL, miss)))
  Xte <- Xte[, colnames(Xtr), drop = FALSE]
  clf <- reticulate::py$`_make_tabpfn_clf`(
    device = device, ensemble = as.integer(ensemble), seed = as.integer(seed)
  )
  clf$fit(Xtr, as.integer(y_train))
  p <- clf$predict_proba(Xte)
  reticulate::py_to_r(p)
}

tabpfn_bin_cindex_cv5 <- function(
    df, method_name,
    seed = 1, K = 5, fold_id = NULL,
    device = "cpu", ensemble = 32,
    event_col = "Event", time_col = "Event_time"
) {
  message(sprintf("[TabPFN Bin] %s: starting %d-fold CV", method_name, K))
  risk_fun <- function(train, test) {
    xtr <- train[, setdiff(names(train), c(time_col, event_col)), drop = FALSE]
    xte <- test[,  setdiff(names(test),  c(time_col, event_col)), drop = FALSE]
    ytr <- as.integer(train[[event_col]])
    pmat <- tabpfn_classify(xtr, ytr, xte, device = device, ensemble = ensemble, seed = seed)
    if (is.matrix(pmat) && ncol(pmat) >= 2) as.numeric(pmat[, 2]) else as.numeric(pmat)
  }
  res <- cv5_cindex(df, risk_fun, seed = seed, reverse = TRUE, K = K,
                    fold_id = fold_id, event_col = event_col, time_col = time_col)
  message(sprintf("[TabPFN Bin] %s: done. CV%d C = %.3f", method_name, K, res$c_overall))
  summarize_cv(method_name, sprintf("TabPFN"), res,
               event_col = event_col, time_col = time_col)
}

