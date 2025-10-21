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
    assay_in = "counts",
    pseudocount = 1e-6,
    event_col = "Event",
    time_col  = "Event_time"
) {
  method <- match.arg(method)
  
  # Validate inputs
  if (!assay_in %in% SummarizedExperiment::assayNames(tse)) {
    stop(sprintf("Assay '%s' not found. Set 'assay_in' correctly.", assay_in))
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
                              pseudocount = pseudocount, name = "alr")
  )
  
  # Extract matrix
  M <- if (method == "lra") {
    assay(altExp(tse_tx, "logratios"))
  } else if (method == "asin") {
    asin(sqrt(assay(tse_tx, "tss")))
  } else {
    assay(tse_tx, method)
  }
  
  # Drop ALR reference
  if (method == "alr") {
    idx <- tryCatch(attributes(M)$parameters$index, error = function(e) NULL)
    if (is.null(idx) || !is.numeric(idx) || length(idx) != 1L) idx <- 1L
    ref_used <- rownames(M)[idx]
    M <- M[setdiff(rownames(M), ref_used), , drop = FALSE]
  }
  
  # Build feature frame
  df <- as.data.frame(t(as.matrix(M)))
  rownames(df) <- colnames(M)
  
  # Clean names
  df <- clean_column_names(df)
  
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

# Generic K-fold CV
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

# Summarize CV
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
  
  out <- tibble::tibble(
    model    = model_name,
    method   = method_name,
    metric   = "C",
    estimate = est,
    lower    = lo,
    upper    = hi,
    se       = se
  )
  
  if (!is.null(cv_res$c_folds)) {
    foldC_tbl <- tibble::tibble(
      model  = model_name,
      method = method_name,
      fold   = seq_along(cv_res$c_folds),
      C      = as.numeric(cv_res$c_folds)
    )
    attr(out, "foldC") <- foldC_tbl
  }
  out
}


# Logistic regression (binomial) K-fold CV + SHAP
logit_cindex_cv5 <- function(
    df, method_name, maxit = 200,
    seed = 1, K = 5, fold_id = NULL,
    event_col = "Event", time_col = "Event_time",
    shap = FALSE, shap_nsim = 64, shap_bg_max = 128, shap_verbose = TRUE
) {
  set.seed(seed)
  
  # folds
  n <- nrow(df)
  if (!is.null(fold_id) && length(fold_id) != n) stop("'fold_id' length must equal nrow(df)")
  if (is.null(fold_id)) fold_id <- make_stratified_folds(df[[event_col]], K = K, seed = seed)
  
  message(sprintf("[Logit] %s: starting %d-fold CV", method_name, K))
  
  preds <- rep(NA_real_, n)
  feat_cols <- setdiff(names(df), c(event_col, time_col))
  shap_all <- if (isTRUE(shap)) list() else NULL
  
  fit_imp <- function(d) {
    is_num <- vapply(d, is.numeric, TRUE)
    list(
      num_cols = names(d)[is_num],
      cat_cols = names(d)[!is_num],
      num_med  = vapply(d[is_num], function(z) median(z, na.rm = TRUE), numeric(1)),
      cat_lvls = lapply(d[!is_num], function(z) levels(addNA(as.factor(z), ifany = TRUE)))
    )
  }
  apply_imp <- function(d, imp) {
    out <- d
    for (nm in imp$num_cols) {
      v <- out[[nm]]
      v[!is.finite(v)] <- NA
      v[is.na(v)] <- imp$num_med[[nm]]
      out[[nm]] <- v
    }
    for (nm in imp$cat_cols) {
      out[[nm]] <- addNA(as.factor(out[[nm]]), ifany = TRUE)
      out[[nm]] <- factor(out[[nm]], levels = imp$cat_lvls[[nm]])
    }
    out
  }
  mm_train_test <- function(xtr, xte) {
    Xtr <- stats::model.matrix(~ . - 1, data = xtr)
    fn  <- colnames(Xtr)
    Xte <- stats::model.matrix(~ . - 1, data = xte)
    miss <- setdiff(fn, colnames(Xte))
    if (length(miss)) Xte <- cbind(Xte, matrix(0, nrow(Xte), length(miss), dimnames = list(NULL, miss)))
    Xte <- Xte[, fn, drop = FALSE]
    list(Xtr = Xtr, Xte = Xte)
  }
  
  for (k in seq_len(K)) {
    if (shap_verbose) message(sprintf("[Logit][%s] fold %d/%d", method_name, k, K))
    
    tr <- fold_id != k; te <- fold_id == k
    train0 <- df[tr, , drop = FALSE]
    test0  <- df[te, , drop = FALSE]
    cast_cols <- function(d) {
      d[] <- lapply(d, function(z) if (is.character(z) || is.logical(z)) factor(z) else z)
      d
    }
    train0 <- cast_cols(train0); test0 <- cast_cols(test0)
    
    # keep only complete cases on features
    cc_tr <- rep(TRUE, nrow(train0))
    cc_te <- rep(TRUE, nrow(test0))
    train <- train0
    test  <- test0
    
    if (nrow(train) < 2L) next
    
    # align dummies using only complete rows
    xtr_raw <- train[, feat_cols, drop = FALSE]
    xte_raw <- test[,  feat_cols, drop = FALSE]
    imp <- fit_imp(xtr_raw)
    xtr_imp <- apply_imp(xtr_raw, imp)
    xte_imp <- apply_imp(xte_raw, imp)
    MM  <- mm_train_test(xtr_imp, xte_imp)
    Xtr <- MM$Xtr
    Xte <- MM$Xte
    
    # drop zero-variance cols by train
    keep0 <- if (requireNamespace("matrixStats", quietly = TRUE)) {
      matrixStats::colSds(Xtr) > 0
    } else {
      apply(Xtr, 2, sd) > 0
    }
    if (!all(keep0)) {
      Xtr <- Xtr[, keep0, drop = FALSE]
      Xte <- Xte[, keep0, drop = FALSE]
    }
    
    # drop collinear columns via QR on TRAIN
    if (ncol(Xtr) > 0) {
      qrX <- qr(Xtr)
      rk  <- qrX$rank
      if (rk < ncol(Xtr)) {
        keep_qr <- qrX$pivot[seq_len(rk)]
        Xtr <- Xtr[, keep_qr, drop = FALSE]
        Xte <- Xte[, keep_qr, drop = FALSE]
      }
    }
    
    # fit GLM (binomial)
    ytr <- as.integer(train[[event_col]])
    df_fit <- data.frame(y = ytr, Xtr, check.names = FALSE)
    fit <- suppressWarnings(stats::glm(
      y ~ .,
      data    = df_fit,
      family  = stats::binomial(),
      control = stats::glm.control(maxit = maxit)
    ))
    
    # OOF risk (probability)
    pr_te <- suppressWarnings(stats::predict(fit, newdata = as.data.frame(Xte), type = "response"))
    pr_te[!is.finite(pr_te)] <- NA_real_
    preds[te][cc_te] <- as.numeric(pr_te)
    
    # SHAP
    if (isTRUE(shap) && nrow(Xtr) > 0) {
      ttX  <- stats::delete.response(stats::terms(fit))
      beta <- stats::coef(fit)
      cols <- setdiff(names(beta), "(Intercept)")
      
      if (length(cols) > 0) {
        bg <- as.data.frame(Xtr)
        if (nrow(bg) > shap_bg_max) bg <- bg[sample.int(nrow(bg), shap_bg_max), , drop = FALSE]
        
        pred_wrap <- function(object, newdata) {
          X <- stats::model.matrix(ttX, data = newdata)
          miss <- setdiff(cols, colnames(X))
          if (length(miss)) {
            X <- cbind(X, matrix(0, nrow(X), length(miss), dimnames = list(NULL, miss)))
          }
          X <- X[, cols, drop = FALSE]
          b <- stats::coef(object)[cols]
          b[!is.finite(b)] <- 0
          eta <- drop(X %*% b); eta[!is.finite(eta)] <- 0
          as.numeric(eta)
        }
        
        shap_raw <- fastshap::explain(
          object = fit,
          X = bg,
          pred_wrapper = pred_wrap,
          newdata = as.data.frame(as.matrix(Xte)),
          nsim = shap_nsim,
          adjust = TRUE
        )
        
        shap_mat <- if (is.matrix(shap_raw) || inherits(shap_raw, "explain")) {
          dm <- dim(shap_raw); dmn <- dimnames(shap_raw)
          out <- unclass(shap_raw); dim(out) <- dm; dimnames(out) <- dmn; out
        } else if (is.data.frame(shap_raw)) {
          as.matrix(shap_raw)
        } else stop("Unexpected SHAP object type: ", paste(class(shap_raw), collapse = "/"))
        
        shap_mat[!is.finite(shap_mat)] <- 0
        
        shap_tbl <- tibble::as_tibble(shap_mat, .name_repair = "minimal")
        idx_rows <- which(te)[cc_te]
        shap_tbl$.row_id <- idx_rows
        shap_tbl <- tidyr::pivot_longer(
          shap_tbl, cols = - .row_id,
          names_to = "feature", values_to = "shap"
        )
        shap_all[[k]] <- shap_tbl
      } else {
        shap_all[[k]] <- tibble::tibble(.row_id = integer(0), feature = character(0), shap = numeric(0))
      }
    }
  }
  
  # fold-wise C
  c_folds <- vapply(
    seq_len(K),
    function(k) {
      idx <- (fold_id == k) & is.finite(preds)
      if (!any(idx)) return(NA_real_)
      harrell_c(df[[time_col]][idx], df[[event_col]][idx], preds[idx], reverse = TRUE)
    },
    numeric(1)
  )
  
  # overall C
  idx_all <- is.finite(preds)
  c_overall <- if (any(idx_all)) harrell_c(df[[time_col]][idx_all], df[[event_col]][idx_all], preds[idx_all], reverse = TRUE) else NA_real_
  
  res <- summarize_cv(
    model_name  = method_name,
    method_name = "Logit",
    cv_res = list(
      data = df, preds = preds, reverse = TRUE,
      K = K, fold_id = fold_id, c_folds = c_folds
    ),
    event_col = event_col, time_col = time_col
  )
  message(sprintf("[Logit] %s: done. CV%d C = %.3f", method_name, K, c_overall))
  
  # attach SHAP
  if (isTRUE(shap)) {
    shap_long <- dplyr::bind_rows(shap_all)
    shap_long$model_key <- "logit"
    shap_long$transform <- method_name
    
    shap_agg <- dplyr::summarise(
      dplyr::group_by(shap_long, model_key, transform, feature),
      mean_abs    = mean(abs(shap), na.rm = TRUE),
      mean_signed = mean(shap,      na.rm = TRUE),
      .groups = "drop"
    ) |>
      dplyr::group_by(model_key, transform) |>
      dplyr::mutate(rel_abs = mean_abs / pmax(sum(mean_abs), .Machine$double.eps)) |>
      dplyr::ungroup()
    
    attr(res, "shap_long") <- shap_long
    attr(res, "shap_agg")  <- shap_agg
  }
  
  return(res)
}


# Random Survival Forest K-fold CV + SHAP
rsf_cindex_cv5 <- function(
    df, method_name,
    num.trees = 1500, min.node.size = 10, mtry = NULL,
    seed = 1, K = 5, fold_id = NULL,
    event_col = "Event", time_col = "Event_time",
    shap = FALSE, shap_nsim = 64, shap_bg_max = 128, shap_verbose = TRUE
) {
  set.seed(seed)
  n <- nrow(df)
  if (!is.null(fold_id) && length(fold_id) != n)
    stop("'fold_id' length must equal nrow(df)")
  if (is.null(fold_id))
    fold_id <- make_stratified_folds(df[[event_col]], K = K, seed = seed)
  
  message(sprintf("[RSF] %s: starting %d-fold CV", method_name, K))
  
  preds <- rep(NA_real_, n)
  feat_cols <- setdiff(names(df), c(event_col, time_col))
  shap_all <- if (isTRUE(shap)) list() else NULL
  
  cast_cols <- function(d) { d[] <- lapply(d, function(x) if (is.character(x)) factor(x) else x); d }
  
  fit_imp <- function(d) {
    is_num <- vapply(d, is.numeric, TRUE)
    list(
      num_cols = names(d)[is_num],
      cat_cols = names(d)[!is_num],
      num_med  = vapply(d[is_num], function(z) median(z, na.rm = TRUE), numeric(1)),
      cat_mode = vapply(d[!is_num], function(z) {
        z <- as.character(z); if (all(is.na(z))) NA_character_ else names(sort(table(z), TRUE))[1]
      }, character(1))
    )
  }
  apply_imp <- function(d, imp) {
    out <- d
    for (nm in imp$num_cols) {
      v <- out[[nm]]; v[!is.finite(v)] <- NA; v[is.na(v)] <- imp$num_med[[nm]]; out[[nm]] <- v
    }
    for (nm in imp$cat_cols) {
      z <- as.character(out[[nm]]); z[is.na(z)] <- imp$cat_mode[[nm]]; out[[nm]] <- factor(z)
    }
    out
  }
  
  for (k in seq_len(K)) {
    if (shap_verbose) message(sprintf("[RSF][%s] fold %d/%d: fit + predict%s",
                                      method_name, k, K, if (shap) " + SHAP" else ""))
    
    tr <- fold_id != k; te <- fold_id == k
    train <- df[tr, , drop = FALSE]
    test  <- df[te, , drop = FALSE]
    
    # Factors for ranger
    train <- cast_cols(train)
    test  <- cast_cols(test)
    
    # Choose mtry
    p_train  <- max(1L, ncol(train) - 2L)
    mtry_use <- if (is.null(mtry)) max(1L, min(p_train, floor(sqrt(p_train)))) else as.integer(mtry)
    
    # fit RSF
    fml <- stats::as.formula(sprintf("survival::Surv(%s, %s) ~ .", time_col, event_col))
    
    cc_tr <- rep(TRUE, nrow(train))
    cc_te <- rep(TRUE, nrow(test))
    train_cc <- train
    test_cc  <- test
    xtr_imp <- apply_imp(train_cc[, feat_cols, drop = FALSE], fit_imp(train_cc[, feat_cols, drop = FALSE]))
    xte_imp <- apply_imp(test_cc[,  feat_cols, drop = FALSE], fit_imp(train_cc[, feat_cols, drop = FALSE]))
    train_cc[, feat_cols] <- xtr_imp
    test_cc[,  feat_cols] <- xte_imp
    if (nrow(train_cc) < 2L || nrow(test_cc) == 0L) next
    
    fit <- ranger::ranger(
      formula = fml,
      data = train_cc,
      num.trees = num.trees,
      mtry = mtry_use,
      min.node.size = min.node.size,
      splitrule = "logrank",
      write.forest = TRUE,
      seed = seed
    )
    
    # OOF risk on test
    chf <- predict(fit, data = test_cc)$chf
    preds[te][cc_te] <- if (is.matrix(chf)) chf[, ncol(chf), drop = TRUE] else as.numeric(chf)
    
    # SHAP
    if (isTRUE(shap)) {
      # Background from TRAIN
      bg <- train_cc[, feat_cols, drop = FALSE]
      if (nrow(bg) > shap_bg_max) bg <- bg[sample.int(nrow(bg), shap_bg_max), , drop = FALSE]
      
      # Prediction wrapper
      pred_wrap <- function(object, newdata) {
        newdata[] <- lapply(newdata, function(z) if (is.character(z)) factor(z) else z)
        chf <- predict(object, data = newdata)$chf
        r   <- if (is.matrix(chf)) chf[, ncol(chf), drop = TRUE] else as.numeric(chf)
        as.numeric(log(pmax(r, 1e-8)))
      }
      
      set.seed(seed)
      shap_raw <- fastshap::explain(
        object = fit,
        X = bg,
        pred_wrapper = pred_wrap,
        newdata = test_cc[, feat_cols, drop = FALSE],
        nsim = shap_nsim,
        adjust = TRUE
      )
      
      shap_mat <- if (is.matrix(shap_raw) || inherits(shap_raw, "explain")) {
        dm  <- dim(shap_raw); dmn <- dimnames(shap_raw)
        shap_mat <- unclass(shap_raw)
        dim(shap_mat) <- dm; dimnames(shap_mat) <- dmn
        shap_mat
      } else if (is.data.frame(shap_raw)) {
        as.matrix(shap_raw)
      } else {
        stop("Unexpected SHAP object type: ", paste(class(shap_raw), collapse = "/"))
      }
      
      shap_tbl <- tibble::as_tibble(shap_mat, .name_repair = "minimal") %>%
        dplyr::mutate(.row_id = which(te)[cc_te]) %>%
        tidyr::pivot_longer(-.row_id, names_to = "feature", values_to = "shap")
      
      shap_all[[k]] <- shap_tbl
    }
  }
  
  # fold-wise C
  c_folds <- vapply(
    seq_len(K),
    function(k) {
      idx <- (fold_id == k) & is.finite(preds)
      if (!any(idx)) return(NA_real_)
      harrell_c(
        df[[time_col]][idx],
        df[[event_col]][idx],
        preds[idx],
        reverse = TRUE
      )
    },
    numeric(1)
  )
  
  # CV metric
  idx_all <- is.finite(preds)
  c_overall <- if (any(idx_all)) harrell_c(df[[time_col]][idx_all], df[[event_col]][idx_all], preds[idx_all], reverse = TRUE) else NA_real_
  res <- summarize_cv(
    model_name = method_name,
    method_name = "RSF",
    cv_res = list(
      data = df, preds = preds, reverse = TRUE,
      K = K, fold_id = fold_id, c_folds = c_folds
    ),
    event_col = event_col, time_col = time_col
  )
  message(sprintf("[RSF] %s: done. CV%d C = %.3f", method_name, K, c_overall))
  
  # Attach SHAP
  if (isTRUE(shap)) {
    shap_all <- lapply(shap_all, function(x) {
      tibble::as_tibble(x)
    })
    
    shap_long <- dplyr::bind_rows(shap_all) %>%
      dplyr::mutate(model_key = "rsf", transform = method_name)
    
    shap_agg <- shap_long %>%
      dplyr::group_by(model_key, transform, feature) %>%
      dplyr::summarise(
        mean_abs    = mean(abs(shap), na.rm = TRUE),
        mean_signed = mean(shap,      na.rm = TRUE),
        .groups = "drop"
      ) %>%
      dplyr::group_by(model_key, transform) %>%
      dplyr::mutate(rel_abs = mean_abs / pmax(sum(mean_abs), .Machine$double.eps)) %>%
      dplyr::ungroup()
    
    attr(res, "shap_long") <- shap_long
    attr(res, "shap_agg")  <- shap_agg
  }
  
  return(res)
}

# DeepSurv K-fold CV + SHAP
deepsurv_cindex_cv5 <- function(
    df, method_name,
    hidden = c(64, 32), dropout = 0.0, l2 = 1e-4,
    lr = 1e-3, epochs = 300, patience = 30,
    verbose = 0, run_eagerly = TRUE,
    seed = 1, K = 5, fold_id = NULL,
    event_col = "Event", time_col = "Event_time",
    shap = FALSE, shap_nsim = 64, shap_bg_max = 128, shap_verbose = TRUE
) {
  
  set.seed(seed)
  n <- nrow(df)
  if (!is.null(fold_id) && length(fold_id) != n) stop("'fold_id' length must equal nrow(df)")
  if (is.null(fold_id)) fold_id <- make_stratified_folds(df[[event_col]], K = K, seed = seed)
  
  message(sprintf("[DeepSurv] %s: starting %d-fold CV", method_name, K))
  
  preds <- rep(NA_real_, n)
  feat_cols <- setdiff(names(df), c(event_col, time_col))
  shap_all <- if (isTRUE(shap)) list() else NULL
  
  cast_cols <- function(d) { d[] <- lapply(d, function(z) if (is.character(z) || is.logical(z)) factor(z) else z); d }
  
  fit_imp <- function(d) {
    is_num <- vapply(d, is.numeric, TRUE)
    list(
      num_cols = names(d)[is_num],
      cat_cols = names(d)[!is_num],
      num_med  = vapply(d[is_num], function(z) median(z, na.rm = TRUE), numeric(1)),
      cat_lvls = lapply(d[!is_num], function(z) levels(addNA(as.factor(z), ifany = TRUE)))
    )
  }
  apply_imp <- function(d, imp) {
    out <- d
    for (nm in imp$num_cols) {
      v <- out[[nm]]; v[!is.finite(v)] <- NA; v[is.na(v)] <- imp$num_med[[nm]]; out[[nm]] <- v
    }
    for (nm in imp$cat_cols) {
      out[[nm]] <- addNA(as.factor(out[[nm]]), ifany = TRUE)
      out[[nm]] <- factor(out[[nm]], levels = imp$cat_lvls[[nm]])
    }
    out
  }
  mm_train_test <- function(xtr, xte) {
    Xtr <- stats::model.matrix(~ . - 1, data = xtr)
    fn  <- colnames(Xtr)
    Xte <- stats::model.matrix(~ . - 1, data = xte)
    miss <- setdiff(fn, colnames(Xte))
    if (length(miss)) Xte <- cbind(Xte, matrix(0, nrow(Xte), length(miss), dimnames = list(NULL, miss)))
    Xte <- Xte[, fn, drop = FALSE]
    list(Xtr = Xtr, Xte = Xte)
  }
  
  for (k in seq_len(K)) {
    if (shap_verbose) message(sprintf("[DeepSurv][%s] fold %d/%d: fit + predict%s",
                                      method_name, k, K, if (shap) " + SHAP" else ""))
    
    tr <- fold_id != k; te <- fold_id == k
    train0 <- df[tr, , drop = FALSE]
    test0  <- df[te, , drop = FALSE]
    
    train0 <- cast_cols(train0)
    test0  <- cast_cols(test0)
    
    cc_tr <- rep(TRUE, nrow(train0))
    cc_te <- rep(TRUE, nrow(test0))
    train <- train0
    test  <- test0
    if (nrow(train) < 2L || nrow(test) == 0L) next
    
    xtr_raw <- train[, feat_cols, drop = FALSE]
    xte_raw <- test[,  feat_cols, drop = FALSE]
    imp <- fit_imp(xtr_raw)
    xtr_imp <- apply_imp(xtr_raw, imp)
    xte_imp <- apply_imp(xte_raw, imp)
    MM      <- mm_train_test(xtr_imp, xte_imp)
    Xtr     <- MM$Xtr
    Xte     <- MM$Xte
    p       <- ncol(Xtr)
    
    mu <- matrixStats::colMeans2(Xtr)
    sd <- matrixStats::colSds(Xtr); sd[!is.finite(sd) | sd < 1e-8] <- 1
    Xtr_s <- sweep(sweep(Xtr, 2, mu, "-"), 2, sd, "/")
    Xte_s <- sweep(sweep(Xte, 2, mu, "-"), 2, sd, "/")
    
    ord       <- order(train[[time_col]])
    Xtr_s_ord <- Xtr_s[ord, , drop = FALSE]
    ev_tr     <- as.numeric(train[[event_col]][ord])
    
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
      x = Xtr_s_ord, y = matrix(ev_tr, ncol = 1),
      batch_size = nrow(Xtr_s_ord), epochs = epochs,
      shuffle = FALSE, verbose = verbose,
      callbacks = list(keras::callback_early_stopping(
        monitor = "loss", patience = patience, restore_best_weights = TRUE))
    )
    
    preds[te][cc_te] <- as.numeric(model$predict(Xte_s, verbose = as.integer(verbose)))
    
    # SHAP
    if (isTRUE(shap)) {
      bg_raw <- train[, feat_cols, drop = FALSE]
      if (nrow(bg_raw) > shap_bg_max) {
        bg_raw <- bg_raw[sample.int(nrow(bg_raw), shap_bg_max), , drop = FALSE]
      }
      
      feat_names <- colnames(Xtr)
      pred_wrap <- function(object, newdata) {
        X <- stats::model.matrix(~ . - 1, data = newdata)
        miss <- setdiff(feat_names, colnames(X))
        if (length(miss)) {
          X <- cbind(X, matrix(0, nrow(X), length(miss), dimnames = list(NULL, miss)))
        }
        X <- X[, feat_names, drop = FALSE]
        Xs <- sweep(sweep(X, 2, mu, "-"), 2, sd, "/")
        as.numeric(object$predict(Xs, verbose = 0))
      }
      
      set.seed(seed)
      shap_raw <- fastshap::explain(
        object = model,
        X = {
          Xbg <- stats::model.matrix(~ . - 1, data = bg_raw)
          miss <- setdiff(feat_names, colnames(Xbg))
          if (length(miss)) {
            Xbg <- cbind(Xbg, matrix(0, nrow(Xbg), length(miss),
                                     dimnames = list(NULL, miss)))
          }
          Xbg <- Xbg[, feat_names, drop = FALSE]
          bg_raw
        },
        pred_wrapper = function(object, newdata) pred_wrap(object, newdata),
        newdata = test[, feat_cols, drop = FALSE],
        nsim = shap_nsim,
        adjust = TRUE
      )
      
      shap_mat <- if (is.matrix(shap_raw) || inherits(shap_raw, "explain")) {
        dm  <- dim(shap_raw); dmn <- dimnames(shap_raw)
        out <- unclass(shap_raw); dim(out) <- dm; dimnames(out) <- dmn; out
      } else if (is.data.frame(shap_raw)) {
        as.matrix(shap_raw)
      } else stop("Unexpected SHAP object type: ", paste(class(shap_raw), collapse = "/"))
      
      shap_tbl <- tibble::as_tibble(shap_mat, .name_repair = "minimal") %>%
        dplyr::mutate(.row_id = which(te)[cc_te]) %>%
        tidyr::pivot_longer(-.row_id, names_to = "feature", values_to = "shap")
      
      shap_all[[k]] <- shap_tbl
    }
  }
  
  # fold-wise C
  c_folds <- vapply(
    seq_len(K),
    function(k) {
      idx <- (fold_id == k) & is.finite(preds)
      if (!any(idx)) return(NA_real_)
      harrell_c(
        df[[time_col]][idx],
        df[[event_col]][idx],
        preds[idx],
        reverse = TRUE
      )
    },
    numeric(1)
  )
  
  # CV metric
  idx_all <- is.finite(preds)
  c_overall <- if (any(idx_all)) harrell_c(df[[time_col]][idx_all], df[[event_col]][idx_all], preds[idx_all], reverse = TRUE) else NA_real_
  res <- summarize_cv(
    model_name = method_name,
    method_name = "DeepSurv",
    cv_res = list(
      data = df, preds = preds, reverse = TRUE,
      K = K, fold_id = fold_id, c_folds = c_folds
    ),
    event_col = event_col, time_col = time_col
  )
  message(sprintf("[DeepSurv] %s: done. CV%d C = %.3f", method_name, K, c_overall))
  
  # Attach SHAP
  if (isTRUE(shap)) {
    shap_all <- lapply(shap_all, tibble::as_tibble)
    shap_long <- dplyr::bind_rows(shap_all) %>%
      dplyr::mutate(model_key = "deepsurv", transform = method_name)
    
    shap_agg <- shap_long %>%
      dplyr::group_by(model_key, transform, feature) %>%
      dplyr::summarise(
        mean_abs    = mean(abs(shap), na.rm = TRUE),
        mean_signed = mean(shap,      na.rm = TRUE),
        .groups = "drop"
      ) %>%
      dplyr::group_by(model_key, transform) %>%
      dplyr::mutate(rel_abs = mean_abs / pmax(sum(mean_abs), .Machine$double.eps)) %>%
      dplyr::ungroup()
    
    attr(res, "shap_long") <- shap_long
    attr(res, "shap_agg")  <- shap_agg
  }
  
  return(res)
}


# Penalized Cox (glmnet) K-fold CV + SHAP
coxnet_cindex_cv5 <- function(
    df, method_name,
    alpha = 0.5,                 # 1=LASSO, 0=Ridge, (0,1)=Elastic Net
    seed = 1, K = 5, fold_id = NULL,
    event_col = "Event", time_col = "Event_time",
    shap = FALSE, shap_nsim = 64, shap_bg_max = 128, shap_verbose = TRUE
) {
  stopifnot(event_col %in% names(df), time_col %in% names(df))
  set.seed(seed)
  
  df <- as.data.frame(df, check.names = FALSE)
  df[] <- lapply(df, function(z) if (is.character(z) || is.logical(z)) factor(z) else z)
  n <- nrow(df)
  if (!is.null(fold_id) && length(fold_id) != n) stop("'fold_id' length must equal nrow(df)")
  if (is.null(fold_id)) fold_id <- make_stratified_folds(df[[event_col]], K = K, seed = seed)
  
  feat_cols <- setdiff(names(df), c(event_col, time_col))
  preds <- rep(NA_real_, n)
  shap_all <- if (isTRUE(shap)) list() else NULL
  
  fit_imp <- function(d) {
    is_num <- vapply(d, is.numeric, TRUE)
    list(
      num_cols = names(d)[is_num],
      cat_cols = names(d)[!is_num],
      num_med  = vapply(d[is_num], function(z) median(z, na.rm = TRUE), numeric(1)),
      cat_lvls = lapply(d[!is_num], function(z) levels(addNA(as.factor(z), ifany = TRUE)))
    )
  }
  apply_imp <- function(d, imp) {
    out <- d
    for (nm in imp$num_cols) {
      v <- out[[nm]]; v[!is.finite(v)] <- NA; v[is.na(v)] <- imp$num_med[[nm]]; out[[nm]] <- v
    }
    for (nm in imp$cat_cols) {
      out[[nm]] <- addNA(as.factor(out[[nm]]), ifany = TRUE)
      out[[nm]] <- factor(out[[nm]], levels = imp$cat_lvls[[nm]])
    }
    out
  }
  mm_train_test <- function(xtr, xte) {
    Xtr <- stats::model.matrix(~ . - 1, data = xtr)
    fn  <- colnames(Xtr)
    Xte <- stats::model.matrix(~ . - 1, data = xte)
    miss <- setdiff(fn, colnames(Xte))
    if (length(miss)) Xte <- cbind(Xte, matrix(0, nrow(Xte), length(miss), dimnames = list(NULL, miss)))
    Xte <- Xte[, fn, drop = FALSE]
    list(Xtr = Xtr, Xte = Xte)
  }
  
  for (k in seq_len(K)) {
    if (shap_verbose) message(sprintf("[CoxNet][%s] fold %d/%d", method_name, k, K))
    
    tr <- fold_id != k; te <- fold_id == k
    train0 <- df[tr, , drop = FALSE]
    test0  <- df[te, , drop = FALSE]
    
    # complete cases on features
    cc_tr <- rep(TRUE, nrow(train0))
    cc_te <- rep(TRUE, nrow(test0))
    train <- train0
    test  <- test0
    if (nrow(train) < 2L) next
    
    Xtr_raw <- train[, feat_cols, drop = FALSE]
    Xte_raw <- test[,  feat_cols, drop = FALSE]
    imp <- fit_imp(Xtr_raw)
    Xtr_i <- apply_imp(Xtr_raw, imp)
    Xte_i <- apply_imp(Xte_raw, imp)
    MM_all  <- mm_train_test(Xtr_i, Xte_i)
    Xtr_mm  <- MM_all$Xtr
    Xte_mm  <- MM_all$Xte
    
    # drop zero-variance cols by train
    keep <- if (requireNamespace("matrixStats", quietly = TRUE)) {
      matrixStats::colSds(Xtr_mm) > 0
    } else {
      apply(Xtr_mm, 2, sd) > 0
    }
    if (!all(keep)) {
      Xtr_mm <- Xtr_mm[, keep, drop = FALSE]
      Xte_mm <- Xte_mm[, keep, drop = FALSE]
    }
    feat_names <- colnames(Xtr_mm)
    
    y_surv <- survival::Surv(train[[time_col]], train[[event_col]])
    
    # inner CV to choose lambda
    cvfit <- glmnet::cv.glmnet(
      x = Xtr_mm, y = y_surv,
      family = "cox",
      alpha  = alpha,
      nfolds = 5,
      standardize = TRUE,
      type.measure = "deviance",
      grouped = TRUE,
      intercept = FALSE
    )
    lam <- cvfit$lambda.min
    
    # refit on full train at chosen lambda
    fit <- glmnet::glmnet(
      x = Xtr_mm, y = y_surv,
      family = "cox",
      alpha  = alpha,
      lambda = lam,
      standardize = TRUE,
      intercept = FALSE
    )
    
    # OOF risk (linear predictor)
    lp_te <- as.numeric(stats::predict(fit, newx = Xte_mm, s = lam, type = "link"))
    lp_te[!is.finite(lp_te)] <- NA_real_
    preds[te][cc_te] <- lp_te
    
    # SHAP via fastshap
    if (isTRUE(shap)) {
      bg_raw <- Xtr_raw
      if (nrow(bg_raw) > shap_bg_max)
        bg_raw <- bg_raw[sample.int(nrow(bg_raw), shap_bg_max), , drop = FALSE]
      
      pred_wrap <- function(object, newdata) {
        X <- stats::model.matrix(~ . - 1, data = newdata)
        miss <- setdiff(feat_names, colnames(X))
        if (length(miss)) {
          X <- cbind(X, matrix(0, nrow(X), length(miss),
                               dimnames = list(NULL, miss)))
        }
        X <- X[, feat_names, drop = FALSE]
        as.numeric(stats::predict(object, newx = X, s = lam, type = "link"))
      }
      
      shap_raw <- fastshap::explain(
        object = fit,
        X = bg_raw,
        pred_wrapper = pred_wrap,
        newdata = Xte_raw,
        nsim = shap_nsim,
        adjust = TRUE
      )
      
      shap_mat <- if (is.matrix(shap_raw) || inherits(shap_raw, "explain")) {
        dm <- dim(shap_raw); dmn <- dimnames(shap_raw)
        out <- unclass(shap_raw); dim(out) <- dm; dimnames(out) <- dmn; out
      } else if (is.data.frame(shap_raw)) {
        as.matrix(shap_raw)
      } else stop("Unexpected SHAP object type: ", paste(class(shap_raw), collapse = "/"))
      
      shap_tbl <- tibble::as_tibble(shap_mat, .name_repair = "minimal") |>
        dplyr::mutate(.row_id = which(te)[cc_te]) |>
        tidyr::pivot_longer(-.row_id, names_to = "feature", values_to = "shap")
      
      shap_all[[k]] <- shap_tbl
    }
  }
  
  # fold-wise C
  c_folds <- vapply(
    seq_len(K),
    function(k) {
      idx <- (fold_id == k) & is.finite(preds)
      if (!any(idx)) return(NA_real_)
      harrell_c(df[[time_col]][idx], df[[event_col]][idx], preds[idx], reverse = TRUE)
    },
    numeric(1)
  )
  
  # overall C
  idx_all <- is.finite(preds)
  c_overall <- if (any(idx_all)) harrell_c(df[[time_col]][idx_all], df[[event_col]][idx_all], preds[idx_all], reverse = TRUE) else NA_real_
  
  res <- summarize_cv(
    model_name  = method_name,
    method_name = "CoxNet",
    cv_res = list(
      data = df, preds = preds, reverse = TRUE,
      K = K, fold_id = fold_id, c_folds = c_folds
    ),
    event_col = event_col, time_col = time_col
  )
  message(sprintf("[CoxNet] %s: done. CV%d C = %.3f", method_name, K, c_overall))
  
  # attach SHAP
  if (isTRUE(shap)) {
    shap_all <- lapply(shap_all, tibble::as_tibble)
    shap_long <- dplyr::bind_rows(shap_all) |>
      dplyr::mutate(model_key = "coxnet", transform = method_name)
    
    shap_agg <- shap_long |>
      dplyr::group_by(model_key, transform, feature) |>
      dplyr::summarise(
        mean_abs    = mean(abs(shap), na.rm = TRUE),
        mean_signed = mean(shap,      na.rm = TRUE),
        .groups = "drop"
      ) |>
      dplyr::group_by(model_key, transform) |>
      dplyr::mutate(rel_abs = mean_abs / pmax(sum(mean_abs), .Machine$double.eps)) |>
      dplyr::ungroup()
    
    attr(res, "shap_long") <- shap_long
    attr(res, "shap_agg")  <- shap_agg
  }
  
  return(res)
}


# CatBoost (binomial) K-fold CV + SHAP
catboost_bin_cindex_cv5 <- function(
    df, method_name,
    seed = 1, K = 5, fold_id = NULL,
    iterations = 1000, depth = 6, learning_rate = 0.05,
    l2_leaf_reg = 3.0, border_count = 254, verbose = 0,
    event_col = "Event", time_col = "Event_time",
    shap = FALSE, shap_nsim = 64, shap_bg_max = 128, shap_verbose = TRUE
) {
  set.seed(seed)
  n <- nrow(df)
  if (!is.null(fold_id) && length(fold_id) != n) stop("'fold_id' length must equal nrow(df)")
  if (is.null(fold_id)) fold_id <- make_stratified_folds(df[[event_col]], K = K, seed = seed)
  
  message(sprintf("[CatBoost Bin] %s: starting %d-fold CV", method_name, K))
  
  preds <- rep(NA_real_, n)
  feat_cols <- setdiff(names(df), c(event_col, time_col))
  shap_all <- if (isTRUE(shap)) list() else NULL
  
  for (k in seq_len(K)) {
    if (shap_verbose) message(sprintf("[CatBoost][%s] fold %d/%d: fit + predict%s",
                                      method_name, k, K, if (shap) " + SHAP" else ""))
    
    tr <- fold_id != k; te <- fold_id == k
    train0 <- df[tr, , drop = FALSE]
    test0  <- df[te, , drop = FALSE]
    
    xtr0 <- train0[, feat_cols, drop = FALSE]
    xte0 <- test0[,  feat_cols, drop = FALSE]
    xtr0[] <- lapply(xtr0, function(v) if (is.character(v) || is.logical(v)) factor(v) else v)
    xte0[] <- lapply(xte0, function(v) if (is.character(v) || is.logical(v)) factor(v) else v)
    
    cc_tr <- rep(TRUE, nrow(xtr0))
    cc_te <- rep(TRUE, nrow(xte0))
    xtr <- xtr0
    xte <- xte0
    train <- train0
    test  <- test0
    
    for (nm in intersect(names(xtr), names(xte))) {
      if (is.factor(xtr[[nm]])) xte[[nm]] <- factor(xte[[nm]], levels = levels(xtr[[nm]]))
    }
    xte <- xte[, names(xtr), drop = FALSE]
    
    if (nrow(xtr) < 2L || nrow(xte) == 0L) next
    
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
    
    preds[te][cc_te] <- as.numeric(catboost::catboost.predict(fit, pool_te, prediction_type = "Probability"))
    
    # SHAP via fastshap
    if (isTRUE(shap)) {
      bg <- xtr
      if (nrow(bg) > shap_bg_max) bg <- bg[sample.int(nrow(bg), shap_bg_max), , drop = FALSE]
      bg <- bg[, names(xtr), drop = FALSE]
      
      pred_wrap <- function(object, newdata) {
        newdata[] <- lapply(newdata, function(v) if (is.character(v) || is.logical(v)) factor(v) else v)
        for (nm in intersect(names(xtr), names(newdata))) {
          if (is.factor(xtr[[nm]])) newdata[[nm]] <- factor(newdata[[nm]], levels = levels(xtr[[nm]]))
        }
        newdata <- newdata[, names(xtr), drop = FALSE]
        pool <- catboost::catboost.load_pool(newdata)
        as.numeric(catboost::catboost.predict(object, pool, prediction_type = "RawFormulaVal"))
      }
      
      set.seed(seed)
      shap_raw <- fastshap::explain(
        object = fit,
        X = bg,
        pred_wrapper = pred_wrap,
        newdata = xte,
        nsim = shap_nsim,
        adjust = TRUE
      )
      
      shap_mat <- if (is.matrix(shap_raw) || inherits(shap_raw, "explain")) {
        dm <- dim(shap_raw); dmn <- dimnames(shap_raw)
        out <- unclass(shap_raw); dim(out) <- dm; dimnames(out) <- dmn; out
      } else if (is.data.frame(shap_raw)) {
        as.matrix(shap_raw)
      } else stop("Unexpected SHAP object type: ", paste(class(shap_raw), collapse = "/"))
      
      shap_tbl <- tibble::as_tibble(shap_mat, .name_repair = "minimal") %>%
        dplyr::mutate(.row_id = which(te)[cc_te]) %>%
        tidyr::pivot_longer(-.row_id, names_to = "feature", values_to = "shap")
      
      shap_all[[k]] <- shap_tbl
    }
  }
  
  # fold-wise C
  c_folds <- vapply(
    seq_len(K),
    function(k) {
      idx <- (fold_id == k) & is.finite(preds)
      if (!any(idx)) return(NA_real_)
      harrell_c(
        df[[time_col]][idx],
        df[[event_col]][idx],
        preds[idx],
        reverse = TRUE
      )
    },
    numeric(1)
  )
  
  # overall C
  idx_all <- is.finite(preds)
  c_overall <- if (any(idx_all)) harrell_c(df[[time_col]][idx_all], df[[event_col]][idx_all], preds[idx_all], reverse = TRUE) else NA_real_
  res <- summarize_cv(
    model_name = method_name,
    method_name = "CatBoost",
    cv_res = list(
      data = df, preds = preds, reverse = TRUE,
      K = K, fold_id = fold_id, c_folds = c_folds
    ),
    event_col = event_col, time_col = time_col
  )
  message(sprintf("[CatBoost Bin] %s: done. CV%d C = %.3f", method_name, K, c_overall))
  
  # attach SHAP
  if (isTRUE(shap)) {
    shap_all <- lapply(shap_all, tibble::as_tibble)
    shap_long <- dplyr::bind_rows(shap_all) %>%
      dplyr::mutate(model_key = "catboost", transform = method_name)
    
    shap_agg <- shap_long %>%
      dplyr::group_by(model_key, transform, feature) %>%
      dplyr::summarise(
        mean_abs    = mean(abs(shap), na.rm = TRUE),
        mean_signed = mean(shap,      na.rm = TRUE),
        .groups = "drop"
      ) %>%
      dplyr::group_by(model_key, transform) %>%
      dplyr::mutate(rel_abs = mean_abs / pmax(sum(mean_abs), .Machine$double.eps)) %>%
      dplyr::ungroup()
    
    attr(res, "shap_long") <- shap_long
    attr(res, "shap_agg")  <- shap_agg
  }
  
  return(res)
}


# XGBoost Cox: inner-grid once, then K-fold CV + SHAP
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
    event_col = "Event", time_col = "Event_time",
    shap = FALSE, shap_nsim = 64, shap_bg_max = 128, shap_verbose = TRUE
) {
  message(sprintf("[XGB] %s: starting %d-fold CV", method_name, K))
  set.seed(seed)
  
  # inner data for tuning
  x_df <- df[, setdiff(names(df), c(time_col, event_col)), drop = FALSE]
  MM   <- stats::model.matrix(~ . - 1, data = x_df)
  eps  <- .Machine$double.eps
  time <- pmax(df[[time_col]], eps)
  y    <- ifelse(df[[event_col]] == 1, time, -time)
  
  inner_id <- make_stratified_folds(df[[event_col]], K = k_inner, seed = seed)
  
  feval_cindex <- function(preds, dmat) {
    ylab <- xgboost::getinfo(dmat, "label")
    t <- abs(ylab); e <- as.integer(ylab > 0)
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
  
  # inner CV
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
  
  # outer CV
  n <- nrow(df)
  if (!is.null(fold_id) && length(fold_id) != n) stop("'fold_id' length must equal nrow(df)")
  if (is.null(fold_id)) fold_id <- make_stratified_folds(df[[event_col]], K = K, seed = seed)
  
  preds <- rep(NA_real_, n)
  feat_cols <- setdiff(names(df), c(event_col, time_col))
  shap_all <- if (isTRUE(shap)) list() else NULL
  
  mm_train_test <- function(xtr, xte) {
    Xtr <- stats::model.matrix(~ . - 1, data = xtr)
    fn  <- colnames(Xtr)
    Xte <- stats::model.matrix(~ . - 1, data = xte)
    miss <- setdiff(fn, colnames(Xte))
    if (length(miss)) Xte <- cbind(Xte, matrix(0, nrow(Xte), length(miss), dimnames = list(NULL, miss)))
    Xte <- Xte[, fn, drop = FALSE]
    list(Xtr = Xtr, Xte = Xte)
  }
  
  for (k in seq_len(K)) {
    if (shap_verbose) message(sprintf("[XGB][%s] fold %d/%d: fit + predict%s",
                                      method_name, k, K, if (shap) " + SHAP" else ""))
    
    tr <- fold_id != k; te <- fold_id == k
    train0 <- df[tr, , drop = FALSE]
    test0  <- df[te, , drop = FALSE]
    
    xtr_df0 <- train0[, feat_cols, drop = FALSE]
    xte_df0 <- test0[,  feat_cols, drop = FALSE]
    cc_tr <- rep(TRUE, nrow(xtr_df0))
    cc_te <- rep(TRUE, nrow(xte_df0))
    xtr_df <- xtr_df0
    xte_df <- xte_df0
    train  <- train0
    test   <- test0
    
    if (nrow(xtr_df) < 2L || nrow(xte_df) == 0L) next
    
    MM_all <- mm_train_test(xtr_df, xte_df)
    Xtr    <- MM_all$Xtr
    Xte    <- MM_all$Xte
    
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
    
    preds[te][cc_te] <- as.numeric(predict(bst, dte))
    
    # SHAP via fastshap
    if (isTRUE(shap)) {
      bg_raw <- xtr_df
      if (nrow(bg_raw) > shap_bg_max)
        bg_raw <- bg_raw[sample.int(nrow(bg_raw), shap_bg_max), , drop = FALSE]
      feat_names <- colnames(Xtr)
      
      pred_wrap <- function(object, newdata) {
        X <- stats::model.matrix(~ . - 1, data = newdata)
        miss <- setdiff(feat_names, colnames(X))
        if (length(miss)) {
          X <- cbind(X, matrix(0, nrow(X), length(miss), dimnames = list(NULL, miss)))
        }
        X <- X[, feat_names, drop = FALSE]
        as.numeric(predict(object, xgboost::xgb.DMatrix(X)))
      }
      
      shap_raw <- fastshap::explain(
        object = bst,
        X = bg_raw,
        pred_wrapper = function(object, newdata) pred_wrap(object, newdata),
        newdata = xte_df,
        nsim = shap_nsim,
        adjust = TRUE
      )
      
      shap_mat <- if (is.matrix(shap_raw) || inherits(shap_raw, "explain")) {
        dm <- dim(shap_raw); dmn <- dimnames(shap_raw)
        out <- unclass(shap_raw); dim(out) <- dm; dimnames(out) <- dmn; out
      } else if (is.data.frame(shap_raw)) {
        as.matrix(shap_raw)
      } else stop("Unexpected SHAP object type: ", paste(class(shap_raw), collapse = "/"))
      
      shap_tbl <- tibble::as_tibble(shap_mat, .name_repair = "minimal") %>%
        dplyr::mutate(.row_id = which(te)[cc_te]) %>%
        tidyr::pivot_longer(-.row_id, names_to = "feature", values_to = "shap")
      
      shap_all[[k]] <- shap_tbl
    }
  }
  
  # fold-wise C
  c_folds <- vapply(
    seq_len(K),
    function(k) harrell_c(
      df[[time_col]][fold_id == k],
      df[[event_col]][fold_id == k],
      preds[fold_id == k],
      reverse = TRUE
    ),
    numeric(1)
  )
  
  # Overall C
  c_overall <- {
    idx <- is.finite(preds)
    if (any(idx)) harrell_c(df[[time_col]][idx], df[[event_col]][idx], preds[idx], reverse = TRUE) else NA_real_
  }
  
  res <- summarize_cv(
    model_name = method_name,
    method_name = "XGB_Cox",
    cv_res = list(
      data = df, preds = preds, reverse = TRUE,
      K = K, fold_id = fold_id, c_folds = c_folds
    ),
    event_col = event_col, time_col = time_col
  )
  message(sprintf("[XGB] %s: done. CV%d C = %.3f", method_name, K, c_overall))
  
  # attach SHAP
  if (isTRUE(shap)) {
    shap_all <- lapply(shap_all, tibble::as_tibble)
    shap_long <- dplyr::bind_rows(shap_all) %>%
      dplyr::mutate(model_key = "xgb_cox", transform = method_name)
    
    shap_agg <- shap_long %>%
      dplyr::group_by(model_key, transform, feature) %>%
      dplyr::summarise(
        mean_abs    = mean(abs(shap), na.rm = TRUE),
        mean_signed = mean(shap,      na.rm = TRUE),
        .groups = "drop"
      ) %>%
      dplyr::group_by(model_key, transform) %>%
      dplyr::mutate(rel_abs = mean_abs / pmax(sum(mean_abs), .Machine$double.eps)) %>%
      dplyr::ungroup()
    
    attr(res, "shap_long") <- shap_long
    attr(res, "shap_agg")  <- shap_agg
  }
  
  return(res)
}


# TabPFN: init helper
tabpfn_init <- function(conda_env = "r-tabpfn") {
  reticulate::use_condaenv(conda_env, required = TRUE)
  if (!reticulate::py_module_available("tabpfn"))
    stop("Python module 'tabpfn' not found in env '", conda_env, "'.")
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

# One-hot encoder
one_hot_mm <- function(df) {
  X <- stats::model.matrix(~ . - 1, data = df)
  storage.mode(X) <- "double"
  X
}

# TabPFN (binomial) K-fold CV + SHAP
tabpfn_bin_cindex_cv5 <- function(
    df, method_name,
    seed = 1, K = 5, fold_id = NULL,
    device = "cpu", ensemble = 32,
    event_col = "Event", time_col = "Event_time",
    shap = FALSE, shap_nsim = 64, shap_bg_max = 128, shap_verbose = TRUE
) {
  message(sprintf("[TabPFN Bin] %s: starting %d-fold CV", method_name, K))
  set.seed(seed)
  
  n <- nrow(df)
  if (!is.null(fold_id) && length(fold_id) != n) stop("'fold_id' length must equal nrow(df)")
  if (is.null(fold_id)) fold_id <- make_stratified_folds(df[[event_col]], K = K, seed = seed)
  
  preds <- rep(NA_real_, n)
  feat_cols <- setdiff(names(df), c(event_col, time_col))
  shap_all <- if (isTRUE(shap)) list() else NULL
  
  fit_imp <- function(d) {
    is_num <- vapply(d, is.numeric, TRUE)
    list(
      num_cols = names(d)[is_num],
      cat_cols = names(d)[!is_num],
      num_med  = vapply(d[is_num], function(z) median(z, na.rm = TRUE), numeric(1)),
      cat_lvls = lapply(d[!is_num], function(z) levels(addNA(as.factor(z), ifany = TRUE)))
    )
  }
  apply_imp <- function(d, imp) {
    out <- d
    for (nm in imp$num_cols) {
      v <- out[[nm]]; v[!is.finite(v)] <- NA; v[is.na(v)] <- imp$num_med[[nm]]; out[[nm]] <- v
    }
    for (nm in imp$cat_cols) {
      out[[nm]] <- addNA(as.factor(out[[nm]]), ifany = TRUE)
      out[[nm]] <- factor(out[[nm]], levels = imp$cat_lvls[[nm]])
    }
    out
  }
  
  for (k in seq_len(K)) {
    if (shap_verbose) message(sprintf("[TabPFN][%s] fold %d/%d: fit + predict%s",
                                      method_name, k, K, if (shap) " + SHAP" else ""))
    
    tr <- fold_id != k; te <- fold_id == k
    train0 <- df[tr, , drop = FALSE]
    test0  <- df[te, , drop = FALSE]
    
    # one-hot + align on complete cases only (features)
    xtr_df0 <- train0[, feat_cols, drop = FALSE]
    xte_df0 <- test0[,  feat_cols, drop = FALSE]
    cc_tr <- rep(TRUE, nrow(xtr_df0))
    cc_te <- rep(TRUE, nrow(xte_df0))
    xtr_df <- xtr_df0
    xte_df <- xte_df0
    train  <- train0
    test   <- test0
    
    if (nrow(xtr_df) < 2L || nrow(xte_df) == 0L) next
    
    imp <- fit_imp(xtr_df)
    xtr_imp <- apply_imp(xtr_df, imp)
    xte_imp <- apply_imp(xte_df, imp)
    
    Xtr <- one_hot_mm(xtr_imp)
    Xte <- one_hot_mm(xte_imp)
    miss <- setdiff(colnames(Xtr), colnames(Xte))
    if (length(miss)) Xte <- cbind(Xte, matrix(0, nrow(Xte), length(miss), dimnames = list(NULL, miss)))
    Xte <- Xte[, colnames(Xtr), drop = FALSE]
    
    ytr <- as.integer(train[[event_col]])
    
    # fit TabPFN
    clf <- reticulate::py$`_make_tabpfn_clf`(
      device = device, ensemble = as.integer(ensemble), seed = as.integer(seed)
    )
    clf$fit(Xtr, as.integer(ytr))
    
    # OOF prob
    p <- reticulate::py_to_r(clf$predict_proba(Xte))
    p1 <- if (is.matrix(p) && ncol(p) >= 2) p[, 2] else as.numeric(p)
    preds[te][cc_te] <- as.numeric(p1)
    
    # SHAP
    if (isTRUE(shap)) {
      bg_raw <- xtr_df
      if (nrow(bg_raw) > shap_bg_max)
        bg_raw <- bg_raw[sample.int(nrow(bg_raw), shap_bg_max), , drop = FALSE]
      feat_names <- colnames(Xtr)
      
      pred_wrap <- function(object, newdata) {
        X <- one_hot_mm(apply_imp(newdata, imp))
        miss <- setdiff(feat_names, colnames(X))
        if (length(miss)) X <- cbind(X, matrix(0, nrow(X), length(miss), dimnames = list(NULL, miss)))
        X <- X[, feat_names, drop = FALSE]
        pr <- reticulate::py_to_r(object$predict_proba(X))
        p1 <- if (is.matrix(pr) && ncol(pr) >= 2) pr[, 2] else as.numeric(pr)
        p1 <- pmin(pmax(p1, 1e-6), 1 - 1e-6)
        qlogis(p1)
      }
      
      shap_raw <- fastshap::explain(
        object = clf,
        X = bg_raw,
        pred_wrapper = pred_wrap,
        newdata = xte_df,
        nsim = shap_nsim,
        adjust = TRUE
      )
      
      shap_mat <- if (is.matrix(shap_raw) || inherits(shap_raw, "explain")) {
        dm  <- dim(shap_raw); dmn <- dimnames(shap_raw)
        out <- unclass(shap_raw); dim(out) <- dm; dimnames(out) <- dmn; out
      } else if (is.data.frame(shap_raw)) {
        as.matrix(shap_raw)
      } else stop("Unexpected SHAP object type: ", paste(class(shap_raw), collapse = "/"))
      
      shap_tbl <- tibble::as_tibble(shap_mat, .name_repair = "minimal") %>%
        dplyr::mutate(.row_id = which(te)[cc_te]) %>%
        tidyr::pivot_longer(-.row_id, names_to = "feature", values_to = "shap")
      
      shap_all[[k]] <- shap_tbl
    }
  }
  
  # fold-wise C
  c_folds <- vapply(
    seq_len(K),
    function(k) harrell_c(
      df[[time_col]][fold_id == k],
      df[[event_col]][fold_id == k],
      preds[fold_id == k],
      reverse = TRUE
    ),
    numeric(1)
  )
  
  # CV summary
  c_overall <- {
    idx <- is.finite(preds)
    if (any(idx)) harrell_c(df[[time_col]][idx], df[[event_col]][idx], preds[idx], reverse = TRUE) else NA_real_
  }
  res <- summarize_cv(
    model_name = method_name,
    method_name = "TabPFN",
    cv_res = list(
      data = df, preds = preds, reverse = TRUE,
      K = K, fold_id = fold_id, c_folds = c_folds
    ),
    event_col = event_col, time_col = time_col
  )
  message(sprintf("[TabPFN Bin] %s: done. CV%d C = %.3f", method_name, K, c_overall))
  
  # attach SHAP
  if (isTRUE(shap)) {
    shap_all <- lapply(shap_all, tibble::as_tibble)
    shap_long <- dplyr::bind_rows(shap_all) %>%
      dplyr::mutate(model_key = "tabpfn", transform = method_name)
    
    shap_agg <- shap_long %>%
      dplyr::group_by(model_key, transform, feature) %>%
      dplyr::summarise(
        mean_abs    = mean(abs(shap), na.rm = TRUE),
        mean_signed = mean(shap,      na.rm = TRUE),
        .groups = "drop"
      ) %>%
      dplyr::group_by(model_key, transform) %>%
      dplyr::mutate(rel_abs = mean_abs / pmax(sum(mean_abs), .Machine$double.eps)) %>%
      dplyr::ungroup()
    
    attr(res, "shap_long") <- shap_long
    attr(res, "shap_agg")  <- shap_agg
  }
  
  return(res)
}



# PERMANOVA
permanova_r2_cv5 <- function(
    df, method_name,
    group_col   = "Event",
    event_col   = "Event",
    time_col    = "Event_time",
    harmonize   = c("zscore","ranknorm","none"),
    metric      = c("euclidean","manhattan","cosine"),
    permutations = 999,
    K           = 5,
    fold_id     = NULL,
    run_on      = c("train","test"),
    shapley     = FALSE,
    shapley_M   = 64,
    shapley_pmax= 64,
    seed        = 1,
    messages    = FALSE
){
  harmonize <- match.arg(harmonize)
  metric    <- match.arg(metric)
  run_on    <- match.arg(run_on)
  set.seed(seed)
  
  msg <- function(...) if (isTRUE(messages)) message(sprintf(...))
  
  # Harmonize matrix columns
  harmonize_matrix <- function(M) {
    M <- as.matrix(M)
    keep <- apply(M, 2, function(v) sd(v, na.rm = TRUE) > 0)
    if (!all(keep)) M <- M[, keep, drop = FALSE]
    if (harmonize == "none") return(M)
    if (harmonize == "zscore") {
      mu  <- if (requireNamespace("matrixStats", quietly = TRUE)) matrixStats::colMeans2(M) else colMeans(M)
      sdv <- if (requireNamespace("matrixStats", quietly = TRUE)) matrixStats::colSds(M)   else apply(M,2,sd)
      sdv[!is.finite(sdv) | sdv < 1e-12] <- 1
      return( sweep(sweep(M,2,mu,"-"),2,sdv,"/") )
    }
    if (harmonize == "ranknorm") {
      rn <- function(x){
        r <- rank(x, ties.method = "average", na.last = "keep"); n <- sum(is.finite(r))
        if (n <= 2) return(rep(0,length(x)))
        p <- (r - 3/8)/(n + 1/4)
        qnorm(p)
      }
      return(apply(M, 2, rn))
    }
  }
  
  # Distance builder
  make_dist <- function(M){
    if (metric %in% c("euclidean","manhattan")) return(stats::dist(M, method = metric))
    if (metric == "cosine") {
      if (!requireNamespace("proxy", quietly = TRUE)) stop("Package 'proxy' is required for cosine.")
      return(proxy::dist(M, method = "cosine"))
    }
  }
  
  # Feature set (exclude time/event and .fold)
  feat_cols_all <- setdiff(names(df), c(event_col, time_col, ".fold"))
  
  # Compute R^2 for a subset
  compute_r2 <- function(df_sub){
    grp <- factor(df_sub[[group_col]])
    X   <- df_sub[, intersect(feat_cols_all, names(df_sub)), drop = FALSE]
    if (!is.matrix(X) || !all(vapply(X, is.numeric, logical(1)))) X <- stats::model.matrix(~ . - 1, data = X)
    X[is.na(X)] <- 0
    Xh <- harmonize_matrix(X)
    D  <- make_dist(Xh)
    fit <- vegan::adonis2(D ~ grp, permutations = permutations, by = "terms")
    rr <- which(rownames(fit) %in% c("grp","grp ")); if (length(rr) != 1L) rr <- 1L
    as.numeric(fit$R2[rr])
  }
  
  # Build folds (stratified 0/1)
  n <- nrow(df)
  if (!is.null(fold_id) && length(fold_id) != n) stop("'fold_id' length must equal nrow(df).")
  if (is.null(fold_id)) {
    ev <- df[[group_col]]
    idx0 <- which(ev == 0L); idx1 <- which(ev == 1L)
    fold_id <- integer(n)
    fold_id[idx0] <- sample(rep(seq_len(K), length.out = length(idx0)))
    fold_id[idx1] <- sample(rep(seq_len(K), length.out = length(idx1)))
  }
  df$.fold <- fold_id
  
  msg("[PERMANOVA][%s] Start CV (K=%d, run_on=%s, metric=%s, harmonize=%s)", method_name, K, run_on, metric, harmonize)
  
  # Fold-wise R^2 on train/test
  fold_tbl <- lapply(seq_len(K), function(k){
    idx <- if (run_on == "train") df$.fold != k else df$.fold == k
    dfk <- df[idx, , drop = FALSE]
    if (nrow(dfk) < 10 || length(unique(dfk[[group_col]])) < 2)
      return(data.frame(model = method_name, method = "PERMANOVA", fold = k, R2 = NA_real_))
    R2 <- tryCatch(compute_r2(dfk), error = function(e) NA_real_)
    data.frame(model = method_name, method = "PERMANOVA", fold = k, R2 = R2)
  }) %>% dplyr::bind_rows() %>% dplyr::filter(is.finite(R2))
  
  # Shapley for overall R^2 (no .row_id, exclude .fold)
  shap_long_tbl <- shap_agg_tbl <- NULL
  if (isTRUE(shapley)) {
    msg("[PERMANOVA][%s] Shapley: M=%d, pmax=%d", method_name, shapley_M, shapley_pmax)
    X <- df[, feat_cols_all, drop = FALSE]
    if (!is.matrix(X) || !all(vapply(X, is.numeric, logical(1)))) X <- stats::model.matrix(~ . - 1, data = X)
    X[is.na(X)] <- 0
    Xh <- harmonize_matrix(X)
    grp <- factor(df[[group_col]])
    feat_names <- colnames(Xh); if (is.null(feat_names)) feat_names <- paste0("V", seq_len(ncol(Xh)))
    
    r2_subset <- function(cols){
      if (length(cols) == 0) return(0)
      Dsub <- make_dist(Xh[, cols, drop = FALSE])
      fit  <- vegan::adonis2(Dsub ~ grp, permutations = permutations, by = "terms")
      rr <- which(rownames(fit) %in% c("grp","grp ")); if (length(rr) != 1L) rr <- 1L
      as.numeric(fit$R2[rr])
    }
    
    contr_list <- vector("list", shapley_M)
    for (m in seq_len(shapley_M)) {
      p <- ncol(Xh)
      idx_sub <- if (p > shapley_pmax) sort(sample.int(p, shapley_pmax)) else seq_len(p)
      ord <- sample(idx_sub, length(idx_sub), replace = FALSE)
      prev <- integer(0); prev_r2 <- 0
      rows <- vector("list", length(ord))
      for (i in seq_along(ord)) {
        j <- ord[i]
        cur <- c(prev, j)
        cur_r2 <- r2_subset(cur)
        rows[[i]] <- data.frame(feature = feat_names[j], shap = cur_r2 - prev_r2)
        prev <- cur; prev_r2 <- cur_r2
      }
      contr_list[[m]] <- dplyr::bind_rows(rows)
    }
    shap_long_tbl <- dplyr::bind_rows(contr_list) %>%
      dplyr::mutate(model_key = "permanova", transform = method_name, .before = 1) %>%
      dplyr::relocate(model_key, transform, feature, shap)
    
    shap_agg_tbl <- shap_long_tbl %>%
      dplyr::group_by(model_key, transform, feature) %>%
      dplyr::summarise(
        mean_abs    = mean(abs(shap), na.rm = TRUE),
        mean_signed = mean(shap,      na.rm = TRUE),
        .groups = "drop"
      ) %>%
      dplyr::group_by(model_key, transform) %>%
      dplyr::mutate(rel_abs = mean_abs / pmax(sum(mean_abs), .Machine$double.eps)) %>%
      dplyr::ungroup()
  }
  
  # CV summary (no n_folds in output)
  m <- nrow(fold_tbl)
  if (isTRUE(is.finite(m)) && m >= 2) {
    est <- mean(fold_tbl$R2, na.rm = TRUE)
    sdv <- stats::sd(fold_tbl$R2, na.rm = TRUE)
    se  <- sdv / sqrt(m)
    tcr <- stats::qt(0.975, df = m - 1)
    lo  <- est - tcr * se
    hi  <- est + tcr * se
  } else if (m == 1) {
    est <- mean(fold_tbl$R2, na.rm = TRUE)
    se  <- NA_real_; lo <- NA_real_; hi <- NA_real_
  } else {
    est <- NA_real_; se <- NA_real_; lo <- NA_real_; hi <- NA_real_
  }
  
  res <- tibble::tibble(
    model    = method_name,
    method   = "PERMANOVA",
    metric   = "R2",
    estimate = est,
    lower    = lo,
    upper    = hi,
    se       = se
  )
  
  # Attach attributes
  attr(res, "foldR2")  <- tibble::as_tibble(fold_tbl)
  if (!is.null(shap_long_tbl)) attr(res, "shap_long") <- shap_long_tbl
  if (!is.null(shap_agg_tbl))  attr(res, "shap_agg")  <- shap_agg_tbl
  
  msg("[PERMANOVA][%s] Done. CV mean R2=%.4f (95%% CI %.4f–%.4f)", method_name, est, lo, hi)
  res
}





##### FUNCTIONS FOR SUBSAMPLING AND TRUNCATION

truncate_followup <- function(df, time_col, event_col, tmax){
  df2 <- df; over <- df2[[time_col]] > tmax
  df2[[time_col]] <- pmin(df2[[time_col]], tmax)
  df2[[event_col]][over] <- 0L
  df2
}
subset_folds <- function(folds, idx_keep) folds[idx_keep]


# CoxNet
run_coxnet_once <- function(p_grid, t_grid, event_col, time_col){
  set.seed(CV_SEED)
  sub_runs <- purrr::map(p_grid, function(pct){
    message(sprintf("[CoxNet] Subsample %.2f", pct))
    n <- length(global_folds); k <- max(2L, floor(pct*n)); idx <- sample.int(n, k)
    trans_s <- purrr::imap(transforms, ~ .x[idx, , drop = FALSE])
    folds_s <- subset_folds(global_folds, idx)
    lst <- purrr::imap(trans_s, function(X, tr){
      res <- coxnet_cindex_cv5(X, tr, alpha=0.5,
                               seed=CV_SEED, K=K_FOLDS, fold_id=folds_s, event_col=event_col, time_col=time_col, shap=FALSE)
      list(metrics = dplyr::mutate(res, transform = tr),
           foldC   = dplyr::mutate(attr(res,"foldC"), transform = tr))
    })
    list(
      metrics = purrr::map(lst,"metrics") %>% dplyr::bind_rows() %>% dplyr::mutate(model="CoxNet", sample_pct=pct),
      foldC   = purrr::map(lst,"foldC")   %>% purrr::compact() %>% dplyr::bind_rows() %>% dplyr::mutate(model="CoxNet", sample_pct=pct)
    )
  })
  metrics_sub <- purrr::map(sub_runs,"metrics") %>% dplyr::bind_rows()
  foldC_sub   <- purrr::map(sub_runs,"foldC")   %>% dplyr::bind_rows()
  saveRDS(metrics_sub, file.path(out_dir, "metrics_subsample_coxnet.rds"))
  saveRDS(foldC_sub,   file.path(out_dir, "foldC_subsample_coxnet.rds"))
  
  trunc_runs <- purrr::map(t_grid, function(tmax){
    message(sprintf("[CoxNet] Trunc %.4g", tmax))
    trans_t <- purrr::imap(transforms, ~ truncate_followup(.x, time_col, event_col, tmax))
    lst <- purrr::imap(trans_t, function(X, tr){
      res <- coxnet_cindex_cv5(X, tr, alpha=0.5,
                               seed=CV_SEED, K=K_FOLDS, fold_id=global_folds, event_col=event_col, time_col=time_col, shap=FALSE)
      list(metrics = dplyr::mutate(res, transform = tr),
           foldC   = dplyr::mutate(attr(res,"foldC"), transform = tr))
    })
    list(
      metrics = purrr::map(lst,"metrics") %>% dplyr::bind_rows() %>% dplyr::mutate(model="CoxNet", followup_time=tmax),
      foldC   = purrr::map(lst,"foldC")   %>% purrr::compact() %>% dplyr::bind_rows() %>% dplyr::mutate(model="CoxNet", followup_time=tmax)
    )
  })
  metrics_tr <- purrr::map(trunc_runs,"metrics") %>% dplyr::bind_rows()
  foldC_tr   <- purrr::map(trunc_runs,"foldC")   %>% dplyr::bind_rows()
  saveRDS(metrics_tr, file.path(out_dir, "metrics_truncation_coxnet.rds"))
  saveRDS(foldC_tr,   file.path(out_dir, "foldC_truncation_coxnet.rds"))
  invisible(NULL)
}

# RSF
run_rsf_once <- function(p_grid, t_grid, event_col, time_col){
  message("[RSF] Start")
  set.seed(CV_SEED)
  sub_runs <- purrr::map(p_grid, function(pct){
    message(sprintf("[RSF] Subsample %.2f", pct))
    n <- length(global_folds); k <- max(2L, floor(pct*n)); idx <- sample.int(n, k)
    trans_s <- purrr::imap(transforms, ~ .x[idx, , drop = FALSE])
    folds_s <- subset_folds(global_folds, idx)
    lst <- purrr::imap(trans_s, function(X, tr){
      res <- rsf_cindex_cv5(X, tr, seed=CV_SEED, K=K_FOLDS, fold_id=folds_s, event_col=event_col, time_col=time_col, shap=FALSE)
      list(metrics = dplyr::mutate(res, transform = tr),
           foldC   = dplyr::mutate(attr(res,"foldC"), transform = tr))
    })
    list(
      metrics = purrr::map(lst,"metrics") |> dplyr::bind_rows() |> dplyr::mutate(model="RSF", sample_pct=pct),
      foldC   = purrr::map(lst,"foldC")   |> purrr::compact() |> dplyr::bind_rows() |> dplyr::mutate(model="RSF", sample_pct=pct)
    )
  })
  metrics_sub <- purrr::map(sub_runs,"metrics") |> dplyr::bind_rows()
  foldC_sub   <- purrr::map(sub_runs,"foldC")   |> dplyr::bind_rows()
  saveRDS(metrics_sub, file.path(out_dir, "metrics_subsample_rsf.rds"))
  saveRDS(foldC_sub,   file.path(out_dir, "foldC_subsample_rsf.rds"))
  
  trunc_runs <- purrr::map(t_grid, function(tmax){
    message(sprintf("[RSF] Trunc %.4g", tmax))
    trans_t <- purrr::imap(transforms, ~ truncate_followup(.x, time_col, event_col, tmax))
    lst <- purrr::imap(trans_t, function(X, tr){
      res <- rsf_cindex_cv5(X, tr, seed=CV_SEED, K=K_FOLDS, fold_id=global_folds, event_col=event_col, time_col=time_col, shap=FALSE)
      list(metrics = dplyr::mutate(res, transform = tr),
           foldC   = dplyr::mutate(attr(res,"foldC"), transform = tr))
    })
    list(
      metrics = purrr::map(lst,"metrics") |> dplyr::bind_rows() |> dplyr::mutate(model="RSF", followup_time=tmax),
      foldC   = purrr::map(lst,"foldC")   |> purrr::compact() |> dplyr::bind_rows() |> dplyr::mutate(model="RSF", followup_time=tmax)
    )
  })
  metrics_tr <- purrr::map(trunc_runs,"metrics") |> dplyr::bind_rows()
  foldC_tr   <- purrr::map(trunc_runs,"foldC")   |> dplyr::bind_rows()
  saveRDS(metrics_tr, file.path(out_dir, "metrics_truncation_rsf.rds"))
  saveRDS(foldC_tr,   file.path(out_dir, "foldC_truncation_rsf.rds"))
  invisible(NULL)
}

# Logit
run_logit_once <- function(p_grid, t_grid, event_col, time_col){
  message("[Logit] Start")
  set.seed(CV_SEED)
  sub_runs <- purrr::map(p_grid, function(pct){
    message(sprintf("[Logit] Subsample %.2f", pct))
    n <- length(global_folds); k <- max(2L, floor(pct*n)); idx <- sample.int(n, k)
    trans_s <- purrr::imap(transforms, ~ .x[idx, , drop = FALSE])
    folds_s <- subset_folds(global_folds, idx)
    lst <- purrr::imap(trans_s, function(X, tr){
      res <- logit_cindex_cv5(X, tr, maxit=200, seed=CV_SEED, K=K_FOLDS, fold_id=folds_s, event_col=event_col, time_col=time_col, shap=FALSE)
      list(metrics = dplyr::mutate(res, transform = tr),
           foldC   = dplyr::mutate(attr(res,"foldC"), transform = tr))
    })
    list(
      metrics = purrr::map(lst,"metrics") |> dplyr::bind_rows() |> dplyr::mutate(model="Logit", sample_pct=pct),
      foldC   = purrr::map(lst,"foldC")   |> purrr::compact() |> dplyr::bind_rows() |> dplyr::mutate(model="Logit", sample_pct=pct)
    )
  })
  metrics_sub <- purrr::map(sub_runs,"metrics") |> dplyr::bind_rows()
  foldC_sub   <- purrr::map(sub_runs,"foldC")   |> dplyr::bind_rows()
  saveRDS(metrics_sub, file.path(out_dir, "metrics_subsample_logit.rds"))
  saveRDS(foldC_sub,   file.path(out_dir, "foldC_subsample_logit.rds"))
  
  trunc_runs <- purrr::map(t_grid, function(tmax){
    message(sprintf("[Logit] Trunc %.4g", tmax))
    trans_t <- purrr::imap(transforms, ~ truncate_followup(.x, time_col, event_col, tmax))
    lst <- purrr::imap(trans_t, function(X, tr){
      res <- logit_cindex_cv5(X, tr, maxit=200, seed=CV_SEED, K=K_FOLDS, fold_id=global_folds, event_col=event_col, time_col=time_col, shap=FALSE)
      list(metrics = dplyr::mutate(res, transform = tr),
           foldC   = dplyr::mutate(attr(res,"foldC"), transform = tr))
    })
    list(
      metrics = purrr::map(lst,"metrics") |> dplyr::bind_rows() |> dplyr::mutate(model="Logit", followup_time=tmax),
      foldC   = purrr::map(lst,"foldC")   |> purrr::compact() |> dplyr::bind_rows() |> dplyr::mutate(model="Logit", followup_time=tmax)
    )
  })
  metrics_tr <- purrr::map(trunc_runs,"metrics") |> dplyr::bind_rows()
  foldC_tr   <- purrr::map(trunc_runs,"foldC")   |> dplyr::bind_rows()
  saveRDS(metrics_tr, file.path(out_dir, "metrics_truncation_logit.rds"))
  saveRDS(foldC_tr,   file.path(out_dir, "foldC_truncation_logit.rds"))
  invisible(NULL)
}

# CatBoost
run_catboost_once <- function(p_grid, t_grid, event_col, time_col){
  message("[CatBoost] Start")
  set.seed(CV_SEED)
  sub_runs <- purrr::map(p_grid, function(pct){
    message(sprintf("[CatBoost] Subsample %.2f", pct))
    n <- length(global_folds); k <- max(2L, floor(pct*n)); idx <- sample.int(n, k)
    trans_s <- purrr::imap(transforms, ~ .x[idx, , drop = FALSE])
    folds_s <- subset_folds(global_folds, idx)
    lst <- purrr::imap(trans_s, function(X, tr){
      res <- catboost_bin_cindex_cv5(X, tr, seed=CV_SEED, K=K_FOLDS, fold_id=folds_s, event_col=event_col, time_col=time_col, shap=FALSE)
      list(metrics = dplyr::mutate(res, transform = tr),
           foldC   = dplyr::mutate(attr(res,"foldC"), transform = tr))
    })
    list(
      metrics = purrr::map(lst,"metrics") |> dplyr::bind_rows() |> dplyr::mutate(model="CatBoost", sample_pct=pct),
      foldC   = purrr::map(lst,"foldC")   |> purrr::compact() |> dplyr::bind_rows() |> dplyr::mutate(model="CatBoost", sample_pct=pct)
    )
  })
  metrics_sub <- purrr::map(sub_runs,"metrics") |> dplyr::bind_rows()
  foldC_sub   <- purrr::map(sub_runs,"foldC")   |> dplyr::bind_rows()
  saveRDS(metrics_sub, file.path(out_dir, "metrics_subsample_catboost.rds"))
  saveRDS(foldC_sub,   file.path(out_dir, "foldC_subsample_catboost.rds"))
  
  trunc_runs <- purrr::map(t_grid, function(tmax){
    message(sprintf("[CatBoost] Trunc %.4g", tmax))
    trans_t <- purrr::imap(transforms, ~ truncate_followup(.x, time_col, event_col, tmax))
    lst <- purrr::imap(trans_t, function(X, tr){
      res <- catboost_bin_cindex_cv5(X, tr, seed=CV_SEED, K=K_FOLDS, fold_id=global_folds, event_col=event_col, time_col=time_col, shap=FALSE)
      list(metrics = dplyr::mutate(res, transform = tr),
           foldC   = dplyr::mutate(attr(res,"foldC"), transform = tr))
    })
    list(
      metrics = purrr::map(lst,"metrics") |> dplyr::bind_rows() |> dplyr::mutate(model="CatBoost", followup_time=tmax),
      foldC   = purrr::map(lst,"foldC")   |> purrr::compact() |> dplyr::bind_rows() |> dplyr::mutate(model="CatBoost", followup_time=tmax)
    )
  })
  metrics_tr <- purrr::map(trunc_runs,"metrics") |> dplyr::bind_rows()
  foldC_tr   <- purrr::map(trunc_runs,"foldC")   |> dplyr::bind_rows()
  saveRDS(metrics_tr, file.path(out_dir, "metrics_truncation_catboost.rds"))
  saveRDS(foldC_tr,   file.path(out_dir, "foldC_truncation_catboost.rds"))
  invisible(NULL)
}

# XGB Cox
run_xgb_once <- function(p_grid, t_grid, event_col, time_col){
  message("[XGB] Start")
  set.seed(CV_SEED)
  sub_runs <- purrr::map(p_grid, function(pct){
    message(sprintf("[XGB] Subsample %.2f", pct))
    n <- length(global_folds); k <- max(2L, floor(pct*n)); idx <- sample.int(n, k)
    trans_s <- purrr::imap(transforms, ~ .x[idx, , drop = FALSE])
    folds_s <- subset_folds(global_folds, idx)
    lst <- purrr::imap(trans_s, function(X, tr){
      res <- xgb_cox_cindex_cv5_grid_tune_once(df=X, method_name=tr, seed=CV_SEED, K=K_FOLDS, k_inner=3,
                                               fold_id=folds_s, use_gpu=FALSE, event_col=event_col, time_col=time_col, shap=FALSE)
      list(metrics = dplyr::mutate(res, transform = tr),
           foldC   = dplyr::mutate(attr(res,"foldC"), transform = tr))
    })
    list(
      metrics = purrr::map(lst,"metrics") |> dplyr::bind_rows() |> dplyr::mutate(model="XGB_Cox", sample_pct=pct),
      foldC   = purrr::map(lst,"foldC")   |> purrr::compact() |> dplyr::bind_rows() |> dplyr::mutate(model="XGB_Cox", sample_pct=pct)
    )
  })
  metrics_sub <- purrr::map(sub_runs,"metrics") |> dplyr::bind_rows()
  foldC_sub   <- purrr::map(sub_runs,"foldC")   |> dplyr::bind_rows()
  saveRDS(metrics_sub, file.path(out_dir, "metrics_subsample_xgb.rds"))
  saveRDS(foldC_sub,   file.path(out_dir, "foldC_subsample_xgb.rds"))
  
  trunc_runs <- purrr::map(t_grid, function(tmax){
    message(sprintf("[XGB] Trunc %.4g", tmax))
    trans_t <- purrr::imap(transforms, ~ truncate_followup(.x, time_col, event_col, tmax))
    lst <- purrr::imap(trans_t, function(X, tr){
      res <- xgb_cox_cindex_cv5_grid_tune_once(df=X, method_name=tr, seed=CV_SEED, K=K_FOLDS, k_inner=3,
                                               fold_id=global_folds, use_gpu=FALSE, event_col=event_col, time_col=time_col, shap=FALSE)
      list(metrics = dplyr::mutate(res, transform = tr),
           foldC   = dplyr::mutate(attr(res,"foldC"), transform = tr))
    })
    list(
      metrics = purrr::map(lst,"metrics") |> dplyr::bind_rows() |> dplyr::mutate(model="XGB_Cox", followup_time=tmax),
      foldC   = purrr::map(lst,"foldC")   |> purrr::compact() |> dplyr::bind_rows() |> dplyr::mutate(model="XGB_Cox", followup_time=tmax)
    )
  })
  metrics_tr <- purrr::map(trunc_runs,"metrics") |> dplyr::bind_rows()
  foldC_tr   <- purrr::map(trunc_runs,"foldC")   |> dplyr::bind_rows()
  saveRDS(metrics_tr, file.path(out_dir, "metrics_truncation_xgb.rds"))
  saveRDS(foldC_tr,   file.path(out_dir, "foldC_truncation_xgb.rds"))
  invisible(NULL)
}

# TabPFN (before DeepSurv)
run_tabpfn_once <- function(p_grid, t_grid, event_col, time_col){
  message("[TabPFN] Start")
  set.seed(CV_SEED)
  sub_runs <- purrr::map(p_grid, function(pct){
    message(sprintf("[TabPFN] Subsample %.2f", pct))
    n <- length(global_folds); k <- max(2L, floor(pct*n)); idx <- sample.int(n, k)
    trans_s <- purrr::imap(transforms, ~ .x[idx, , drop = FALSE])
    folds_s <- subset_folds(global_folds, idx)
    lst <- purrr::imap(trans_s, function(X, tr){
      res <- tabpfn_bin_cindex_cv5(X, tr, seed=CV_SEED, K=K_FOLDS, fold_id=folds_s, device="cpu", ensemble=1,
                                   event_col=event_col, time_col=time_col, shap=FALSE)
      list(metrics = dplyr::mutate(res, transform = tr),
           foldC   = dplyr::mutate(attr(res,"foldC"), transform = tr))
    })
    list(
      metrics = purrr::map(lst,"metrics") |> dplyr::bind_rows() |> dplyr::mutate(model="TabPFN", sample_pct=pct),
      foldC   = purrr::map(lst,"foldC")   |> purrr::compact() |> dplyr::bind_rows() |> dplyr::mutate(model="TabPFN", sample_pct=pct)
    )
  })
  metrics_sub <- purrr::map(sub_runs,"metrics") |> dplyr::bind_rows()
  foldC_sub   <- purrr::map(sub_runs,"foldC")   |> dplyr::bind_rows()
  saveRDS(metrics_sub, file.path(out_dir, "metrics_subsample_tabpfn.rds"))
  saveRDS(foldC_sub,   file.path(out_dir, "foldC_subsample_tabpfn.rds"))
  
  trunc_runs <- purrr::map(t_grid, function(tmax){
    message(sprintf("[TabPFN] Trunc %.4g", tmax))
    trans_t <- purrr::imap(transforms, ~ truncate_followup(.x, time_col, event_col, tmax))
    lst <- purrr::imap(trans_t, function(X, tr){
      res <- tabpfn_bin_cindex_cv5(X, tr, seed=CV_SEED, K=K_FOLDS, fold_id=global_folds, device="cpu", ensemble=1,
                                   event_col=event_col, time_col=time_col, shap=FALSE)
      list(metrics = dplyr::mutate(res, transform = tr),
           foldC   = dplyr::mutate(attr(res,"foldC"), transform = tr))
    })
    list(
      metrics = purrr::map(lst,"metrics") |> dplyr::bind_rows() |> dplyr::mutate(model="TabPFN", followup_time=tmax),
      foldC   = purrr::map(lst,"foldC")   |> purrr::compact() |> dplyr::bind_rows() |> dplyr::mutate(model="TabPFN", followup_time=tmax)
    )
  })
  metrics_tr <- purrr::map(trunc_runs,"metrics") |> dplyr::bind_rows()
  foldC_tr   <- purrr::map(trunc_runs,"foldC")   |> dplyr::bind_rows()
  saveRDS(metrics_tr, file.path(out_dir, "metrics_truncation_tabpfn.rds"))
  saveRDS(foldC_tr,   file.path(out_dir, "foldC_truncation_tabpfn.rds"))
  invisible(NULL)
}

# DeepSurv (after TabPFN)
run_deepsurv_once <- function(p_grid, t_grid, event_col, time_col){
  message("[DeepSurv] Start")
  set.seed(CV_SEED)
  sub_runs <- purrr::map(p_grid, function(pct){
    message(sprintf("[DeepSurv] Subsample %.2f", pct))
    n <- length(global_folds); k <- max(2L, floor(pct*n)); idx <- sample.int(n, k)
    trans_s <- purrr::imap(transforms, ~ .x[idx, , drop = FALSE])
    folds_s <- subset_folds(global_folds, idx)
    lst <- purrr::imap(trans_s, function(X, tr){
      res <- deepsurv_cindex_cv5(X, tr, hidden=c(32,16), dropout=0.1, l2=1e-4, lr=1e-3, epochs=120, patience=12,
                                 verbose=0, run_eagerly=TRUE, seed=CV_SEED, K=K_FOLDS, fold_id=folds_s,
                                 event_col=event_col, time_col=time_col, shap=FALSE)
      list(metrics = dplyr::mutate(res, transform = tr),
           foldC   = dplyr::mutate(attr(res,"foldC"), transform = tr))
    })
    list(
      metrics = purrr::map(lst,"metrics") |> dplyr::bind_rows() |> dplyr::mutate(model="DeepSurv", sample_pct=pct),
      foldC   = purrr::map(lst,"foldC")   |> purrr::compact() |> dplyr::bind_rows() |> dplyr::mutate(model="DeepSurv", sample_pct=pct)
    )
  })
  metrics_sub <- purrr::map(sub_runs,"metrics") |> dplyr::bind_rows()
  foldC_sub   <- purrr::map(sub_runs,"foldC")   |> dplyr::bind_rows()
  saveRDS(metrics_sub, file.path(out_dir, "metrics_subsample_deepsurv.rds"))
  saveRDS(foldC_sub,   file.path(out_dir, "foldC_subsample_deepsurv.rds"))
  
  trunc_runs <- purrr::map(t_grid, function(tmax){
    message(sprintf("[DeepSurv] Trunc %.4g", tmax))
    trans_t <- purrr::imap(transforms, ~ truncate_followup(.x, time_col, event_col, tmax))
    lst <- purrr::imap(trans_t, function(X, tr){
      res <- deepsurv_cindex_cv5(X, tr, hidden=c(32,16), dropout=0.1, l2=1e-4, lr=1e-3, epochs=120, patience=12,
                                 verbose=0, run_eagerly=TRUE, seed=CV_SEED, K=K_FOLDS, fold_id=global_folds,
                                 event_col=event_col, time_col=time_col, shap=FALSE)
      list(metrics = dplyr::mutate(res, transform = tr),
           foldC   = dplyr::mutate(attr(res,"foldC"), transform = tr))
    })
    list(
      metrics = purrr::map(lst,"metrics") |> dplyr::bind_rows() |> dplyr::mutate(model="DeepSurv", followup_time=tmax),
      foldC   = purrr::map(lst,"foldC")   |> purrr::compact() |> dplyr::bind_rows() |> dplyr::mutate(model="DeepSurv", followup_time=tmax)
    )
  })
  metrics_tr <- purrr::map(trunc_runs,"metrics") |> dplyr::bind_rows()
  foldC_tr   <- purrr::map(trunc_runs,"foldC")   |> dplyr::bind_rows()
  saveRDS(metrics_tr, file.path(out_dir, "metrics_truncation_deepsurv.rds"))
  saveRDS(foldC_tr,   file.path(out_dir, "foldC_truncation_deepsurv.rds"))
  invisible(NULL)
}

# PERMANOVA (R2)
run_permanova_once <- function(p_grid, t_grid, event_col, time_col){
  message("[PERMANOVA] Start")
  set.seed(CV_SEED)
  sub_runs <- purrr::map(p_grid, function(pct){
    message(sprintf("[PERMANOVA] Subsample %.2f", pct))
    n <- length(global_folds); k <- max(2L, floor(pct*n)); idx <- sample.int(n, k)
    trans_s <- purrr::imap(transforms, ~ .x[idx, , drop = FALSE])
    folds_s <- subset_folds(global_folds, idx)
    lst <- purrr::imap(trans_s, function(X, tr){
      res <- permanova_r2_cv5(X, tr, group_col=event_col, event_col=event_col, time_col=time_col,
                              harmonize="zscore", metric="euclidean", permutations=999,
                              K=K_FOLDS, fold_id=folds_s, run_on="train", shapley=FALSE, seed=CV_SEED, messages=FALSE)
      list(metrics = dplyr::mutate(res, transform = tr, method="PERMANOVA"),
           foldR2  = dplyr::mutate(attr(res,"foldR2"), transform = tr))
    })
    list(
      metrics = purrr::map(lst,"metrics") |> dplyr::bind_rows() |> dplyr::mutate(model="PERMANOVA", sample_pct=pct),
      foldR2  = purrr::map(lst,"foldR2")  |> purrr::compact() |> dplyr::bind_rows() |> dplyr::mutate(model="PERMANOVA", sample_pct=pct)
    )
  })
  metrics_sub <- purrr::map(sub_runs,"metrics") |> dplyr::bind_rows()
  foldR2_sub  <- purrr::map(sub_runs,"foldR2")  |> dplyr::bind_rows()
  saveRDS(metrics_sub, file.path(out_dir, "metrics_subsample_permanova.rds"))
  saveRDS(foldR2_sub,  file.path(out_dir, "foldR2_subsample_permanova.rds"))
  
  trunc_runs <- purrr::map(t_grid, function(tmax){
    message(sprintf("[PERMANOVA] Trunc %.4g", tmax))
    trans_t <- purrr::imap(transforms, ~ truncate_followup(.x, time_col, event_col, tmax))
    lst <- purrr::imap(trans_t, function(X, tr){
      res <- permanova_r2_cv5(X, tr, group_col=event_col, event_col=event_col, time_col=time_col,
                              harmonize="zscore", metric="euclidean", permutations=999,
                              K=K_FOLDS, fold_id=global_folds, run_on="train", shapley=FALSE, seed=CV_SEED, messages=FALSE)
      list(metrics = dplyr::mutate(res, transform = tr, method="PERMANOVA"),
           foldR2  = dplyr::mutate(attr(res,"foldR2"), transform = tr))
    })
    list(
      metrics = purrr::map(lst,"metrics") |> dplyr::bind_rows() |> dplyr::mutate(model="PERMANOVA", followup_time=tmax),
      foldR2  = purrr::map(lst,"foldR2")  |> purrr::compact() |> dplyr::bind_rows() |> dplyr::mutate(model="PERMANOVA", followup_time=tmax)
    )
  })
  metrics_tr <- purrr::map(trunc_runs,"metrics") |> dplyr::bind_rows()
  foldR2_tr  <- purrr::map(trunc_runs,"foldR2")  |> dplyr::bind_rows()
  saveRDS(metrics_tr, file.path(out_dir, "metrics_truncation_permanova.rds"))
  saveRDS(foldR2_tr,  file.path(out_dir, "foldR2_truncation_permanova.rds"))
  invisible(NULL)
}
