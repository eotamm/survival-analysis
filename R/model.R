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

