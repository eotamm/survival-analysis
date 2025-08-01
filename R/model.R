# Run Bayesian survival analysis across different data transformations and feature set sizes,
# with feature counts selected automatically based on the number of taxa in the input data (tse)
set.seed(2025)
res_brm_result <- compare_transforms_with_brm(tse)

# Save the results to file
saveRDS(res_brm_result, "feature_survival_result.rds")

# Fit probabilistic Cox models for each transformation output
fit_clr   <- fit_model(df_clr,   model_name = "fit_clr")
fit_rclr  <- fit_model(df_rclr,  model_name = "fit_rclr")
fit_abund <- fit_model(df_abund, model_name = "fit_abund")
fit_lra   <- fit_model(df_lra,   model_name = "fit_lra")

# Calculate LOO for each model individually with moment matching
loo_lra   <- loo(fit_lra,   moment_match = TRUE)
loo_clr   <- loo(fit_clr,   moment_match = TRUE)
loo_rclr  <- loo(fit_rclr,  moment_match = TRUE)
loo_abund <- loo(fit_abund, moment_match = TRUE)

# List of all LOO objects for convenient access
loos <- list(
  lra   = loo_lra,
  clr   = loo_clr,
  rclr  = loo_rclr,
  abund = loo_abund
)

# Compare models using LOO
loo_comp <- loo_compare(
  loo_lra,
  loo_clr,
  loo_rclr,
  loo_abund
)

# Save results
saveRDS(loos, file = "loos.rds")
saveRDS(loo_comp, file = "loo_comparison.rds")
