# Run Bayesian survival analysis across different data transformations and feature set sizes,
# with feature counts selected automatically based on the number of taxa in the input data (tse)
set.seed(1)
if (!dir.exists("model_result")) {
  dir.create("model_result")
}

res_brm_result <- compare_transforms_with_brm(tse)

# Save the results to file
saveRDS(res_brm_result, file.path("model_result", "feature_survival_result2.rds"))


# Fit probabilistic Cox models for each transformation output
set.seed(1)
# Normal prior
fit_clr_normal    <- fit_model(df_clr,    model_name = "fit_clr_normal",    prior_type = "normal")
fit_rclr_normal   <- fit_model(df_rclr,   model_name = "fit_rclr_normal",   prior_type = "normal")
fit_logtss_normal <- fit_model(df_logtss, model_name = "fit_logtss_normal", prior_type = "normal")
fit_lra_normal    <- fit_model(df_lra,    model_name = "fit_lra_normal",    prior_type = "normal")
fit_pa_normal     <- fit_model(df_pa,     model_name = "fit_pa_normal",     prior_type = "normal")
fit_tss_normal    <- fit_model(df_tss,    model_name = "fit_tss_normal",    prior_type = "normal")
fit_asin_normal   <- fit_model(df_asin,   model_name = "fit_asin_normal",   prior_type = "normal")
fit_alr_normal    <- fit_model(df_alr,    model_name = "fit_alr_normal",    prior_type = "normal")

# Horseshoe prior
fit_clr_hs    <- fit_model(df_clr,    model_name = "fit_clr_hs",    prior_type = "horseshoe")
fit_rclr_hs   <- fit_model(df_rclr,   model_name = "fit_rclr_hs",   prior_type = "horseshoe")
fit_logtss_hs <- fit_model(df_logtss, model_name = "fit_logtss_hs", prior_type = "horseshoe")
fit_lra_hs    <- fit_model(df_lra,    model_name = "fit_lra_hs",    prior_type = "horseshoe")
fit_pa_hs     <- fit_model(df_pa,     model_name = "fit_pa_hs",     prior_type = "horseshoe")
fit_tss_hs    <- fit_model(df_tss,    model_name = "fit_tss_hs",    prior_type = "horseshoe")
fit_asin_hs   <- fit_model(df_asin,   model_name = "fit_asin_hs",   prior_type = "horseshoe")
fit_alr_hs    <- fit_model(df_alr,    model_name = "fit_alr_hs",    prior_type = "horseshoe")


# LOO for normal-prior models
models_normal <- list(
  clr    = fit_clr_normal,
  rclr   = fit_rclr_normal,
  logtss = fit_logtss_normal,
  lra    = fit_lra_normal,
  pa     = fit_pa_normal,
  tss    = fit_tss_normal,
  asin   = fit_asin_normal,
  alr    = fit_alr_normal
)

loo_normal <- lapply(models_normal, function(f) {
  loo(f, moment_match = TRUE)
})

# Comapre
loo_comp_normal <- loo_compare(list(
  clr_normal    = loo_normal$clr,
  rclr_normal   = loo_normal$rclr,
  logtss_normal = loo_normal$logtss,
  lra_normal    = loo_normal$lra,
  pa_normal     = loo_normal$pa,
  tss_normal    = loo_normal$tss,
  asin_normal   = loo_normal$asin,
  alr_normal    = loo_normal$alr
))
saveRDS(loo_normal, file = file.path("model_result", "loo_normal.rds"))
saveRDS(loo_comp_normal, file = file.path("model_result", "loo_compare_normal.rds"))


# LOO for horseshoe-prior models
models_hs <- list(
  clr    = fit_clr_hs,
  rclr   = fit_rclr_hs,
  logtss = fit_logtss_hs,
  lra    = fit_lra_hs,
  pa     = fit_pa_hs,
  tss    = fit_tss_hs,
  asin   = fit_asin_hs,
  alr    = fit_alr_hs
)

loo_hs <- lapply(models_hs, function(f) {
  loo(f, moment_match = TRUE)
})

# Compare
loo_comp_hs <- loo_compare(list(
  clr_hs    = loo_hs$clr,
  rclr_hs   = loo_hs$rclr,
  logtss_hs = loo_hs$logtss,
  lra_hs    = loo_hs$lra,
  pa_hs     = loo_hs$pa,
  tss_hs    = loo_hs$tss,
  asin_hs   = loo_hs$asin,
  alr_hs    = loo_hs$alr
))
saveRDS(loo_hs, file = file.path("model_result", "loo_hs.rds"))
saveRDS(loo_comp_hs, file = file.path("model_result", "loo_compare_hs.rds"))


# Collect AUC results for all models
auc_all <- bind_rows(
  get_auc_posterior(fit_rclr_normal,   "rCLR",              "Normal"),
  get_auc_posterior(fit_clr_normal,    "CLR",               "Normal"),
  get_auc_posterior(fit_logtss_normal, "logTSS",            "Normal"),
  get_auc_posterior(fit_lra_normal,    "LRA",               "Normal"),
  get_auc_posterior(fit_pa_normal,     "PA",                "Normal"),
  get_auc_posterior(fit_tss_normal,    "TSS",               "Normal"),
  get_auc_posterior(fit_asin_normal,   "Asin",              "Normal"),
  get_auc_posterior(fit_alr_normal,    "ALR",               "Normal"),
  
  get_auc_posterior(fit_rclr_hs,   "rCLR",              "Horseshoe"),
  get_auc_posterior(fit_clr_hs,    "CLR",               "Horseshoe"),
  get_auc_posterior(fit_logtss_hs, "logTSS",            "Horseshoe"),
  get_auc_posterior(fit_lra_hs,    "LRA",               "Horseshoe"),
  get_auc_posterior(fit_pa_hs,     "PA",                "Horseshoe"),
  get_auc_posterior(fit_tss_hs,    "TSS",               "Horseshoe"),
  get_auc_posterior(fit_asin_hs,   "Asin",              "Horseshoe"),
  get_auc_posterior(fit_alr_hs,    "ALR",               "Horseshoe")
)

# Save
saveRDS(auc_all, file.path("model_result", "posterior_auc_all.rds"))


# Collect C-index results
c_all <- dplyr::bind_rows(
  get_c_posterior(fit_rclr_normal,   "rCLR",              "Normal"),
  get_c_posterior(fit_clr_normal,    "CLR",               "Normal"),
  get_c_posterior(fit_logtss_normal, "logTSS",            "Normal"),
  get_c_posterior(fit_lra_normal,    "LRA",               "Normal"),
  get_c_posterior(fit_pa_normal,     "PA",                "Normal"),
  get_c_posterior(fit_tss_normal,    "TSS",               "Normal"),
  get_c_posterior(fit_asin_normal,   "Asin",              "Normal"),
  get_c_posterior(fit_alr_normal,    "ALR",               "Normal"),
  
  get_c_posterior(fit_rclr_hs,   "rCLR",              "Horseshoe"),
  get_c_posterior(fit_clr_hs,    "CLR",               "Horseshoe"),
  get_c_posterior(fit_logtss_hs, "logTSS",            "Horseshoe"),
  get_c_posterior(fit_lra_hs,    "LRA",               "Horseshoe"),
  get_c_posterior(fit_pa_hs,     "PA",                "Horseshoe"),
  get_c_posterior(fit_tss_hs,    "TSS",               "Horseshoe"),
  get_c_posterior(fit_asin_hs,   "Asin",              "Horseshoe"),
  get_c_posterior(fit_alr_hs,    "ALR",               "Horseshoe")
)

# Save
saveRDS(c_all, file.path("model_result", "cindex_posterior_all.rds"))


# Input datasets for each transformation
transforms <- list(
  CLR    = df_clr_all,
  rCLR   = df_rclr_all,
  logTSS = df_logtss_all,
  LRA    = df_lra_all,
  PA     = df_pa_all,
  TSS    = df_tss_all,
  Asin   = df_asin_all,
  ALR    = df_alr_all
)

# CV setup
CV_SEED <- 1
K_FOLDS <- 5
set.seed(CV_SEED)

# Shared folds for all transforms (same order/length across data frames)
global_folds <- make_stratified_folds(df_clr$Event, K = K_FOLDS, seed = CV_SEED)

# CoxPH
cox_res <- purrr::imap_dfr(
  transforms,
  ~ coxph_cindex_cv5(.x, .y, ties = "efron",
                     seed = CV_SEED, K = K_FOLDS, fold_id = global_folds)
)

# Random Survival Forest
rsf_res <- purrr::imap_dfr(
  transforms,
  ~ rsf_cindex_cv5(.x, .y, seed = CV_SEED, K = K_FOLDS, fold_id = global_folds)
)

# Logistic (binary)
logit_res <- purrr::imap_dfr(
  transforms,
  ~ logit_cindex_cv5(.x, .y, maxit = 200,
                     seed = CV_SEED, K = K_FOLDS, fold_id = global_folds)
)

# DeepSurv
deepsurv_res <- purrr::imap_dfr(
  transforms,
  ~ deepsurv_cindex_cv5(
    .x, .y,
    hidden = c(32, 16), dropout = 0.1, l2 = 1e-4,
    lr = 1e-3, epochs = 120, patience = 12,
    verbose = 0, run_eagerly = TRUE,
    seed = CV_SEED, K = K_FOLDS, fold_id = global_folds
  )
)

# CatBoost (binomial)
cb_res <- purrr::imap_dfr(
  transforms,
  ~ catboost_bin_cindex_cv5(.x, .y, seed = CV_SEED, K = K_FOLDS, fold_id = global_folds)
)

# XGBoost Cox with one-time grid tuning (inner 3-fold), then outer K-fold
xgb_res <- purrr::imap_dfr(
  transforms,
  ~ xgb_cox_cindex_cv5_grid_tune_once(
    df = .x, method_name = .y,
    seed = CV_SEED, K = K_FOLDS, k_inner = 3,
    depth_grid     = c(3L, 4L),
    min_child_grid = c(4L, 8L),
    subsample_grid = c(0.70, 0.85, 1.00),
    colsample_grid = c(0.50, 0.70, 0.90),
    eta_grid       = c(0.05, 0.08, 0.12),
    lambda_grid    = c(0, 2, 5, 10),
    alpha_grid     = c(0, 0.10, 0.50, 1.00),
    nrounds_max = 2000, early_stopping_rounds = 50,
    fold_id = global_folds,
    use_gpu = FALSE
  )
)

# TabPFN (binary)
tabpfn_res <- purrr::imap_dfr(
  transforms,
  ~ tabpfn_bin_cindex_cv5(
    df = .x, method_name = .y,
    seed = CV_SEED, K = K_FOLDS, fold_id = global_folds,
    device = "cpu", ensemble = 32
  )
)

# Combine all model summaries into one
all_models <- dplyr::bind_rows(
  cox_res, rsf_res, logit_res, deepsurv_res, xgb_res, cb_res, tabpfn_res
)

# Save 
out_dir <- file.path("model_result", "results")
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)
tag <- paste0("cv", K_FOLDS)
saveRDS(cox_res,       file.path(out_dir, paste0("coxph_cindex_",       tag, "_summary.rds")))
saveRDS(rsf_res,       file.path(out_dir, paste0("rsf_cindex_",         tag, "_summary.rds")))
saveRDS(logit_res,     file.path(out_dir, paste0("logit_cindex_",       tag, "_summary.rds")))
saveRDS(deepsurv_res,  file.path(out_dir, paste0("deepsurv_cindex_",    tag, "_summary.rds")))
saveRDS(xgb_res,       file.path(out_dir, paste0("xgb_cindex_",         tag, "_summary.rds")))
saveRDS(cb_res,        file.path(out_dir, paste0("catboost_bin_cindex_",tag, "_summary.rds")))
saveRDS(tabpfn_res,    file.path(out_dir, paste0("tabpfn_bin_cindex_",  tag, "_summary.rds")))

# Combined 
saveRDS(all_models, file.path(out_dir, paste0("all_models_", tag, "_with_ci.rds")))
