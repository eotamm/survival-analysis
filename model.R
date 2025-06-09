# Fit probabilistic Cox model with brms using helper function
fit_brms <- fit_model(dff, model_name = "fit_brms")

# Fit log-relative-abundance model
fit_lra_ratios  <- fit_model(df_lra_ratios, "fit_lra_ratios")

# Clr vs rclr models
fit_clr_ratios <- fit_model(df_clr_ratios, "fit_clr_ratios")
fit_rclr_ratios <- fit_model(df_rclr_ratios, "fit_rclr_ratios")