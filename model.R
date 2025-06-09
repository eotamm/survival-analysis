# Fit probabilistic Cox model with brms using helper function
fit_brms <- fit_model(dff, model_name = "fit_brms")

# Fit the model of top15 microbes
fit_abund_top15 <- fit_model(df_abund, "fit_abund_top15")

# Fit log-relative-abundance model
fit_lra_ratios  <- fit_model(df_lra_ratios, "fit_lra_ratios")

# Clr vs rclr models
fit_clr_ratios <- fit_model(df_clr_ratios, "fit_clr_ratios")
fit_rclr_ratios <- fit_model(df_rclr_ratios, "fit_rclr_ratios")