# Fit probabilistic Cox model with brms using helper function
fit_brms <- fit_model(df, model_name = "fit_brms")

# Fit the model using log-transformed relative abundances (no ratios)
fit_abund_top20 <- fit_model(df_abund, "fit_abund_top20")

# Fit the model using log-ratios from log-relative-abundancesgi
fit_lra_ratios <- fit_model(df_lra_ratios, "fit_lra_ratios")

# Fit the model using CLR-transformed features
fit_clr <- fit_model(df_clr, "fit_clr")

# Fit the model using rCLR-transformed features
fit_rclr <- fit_model(df_rclr, "fit_rclr")

# Calculate LOO for each model individually with moment matching
loo_lra   <- loo(fit_lra_ratios, moment_match = TRUE)
loo_clr   <- loo(fit_clr, moment_match = TRUE)
loo_rclr  <- loo(fit_rclr, moment_match = TRUE)
loo_abund <- loo(fit_abund_top20, moment_match = TRUE)

# List of all LOO objects
loos <- list(
  lra = loo_lra,
  clr = loo_clr,
  rclr = loo_rclr,
  abund = loo_abund
)

# Compare models using LOO
loo_comp <- loo_compare(
  loo_lra,
  loo_clr,
  loo_rclr,
  loo_abund
)

# Save
saveRDS(loos, file = "loos.rds")
saveRDS(loo_comp, file = "loo_comparison.rds")
