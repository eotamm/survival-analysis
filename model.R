# Fit probabilistic Cox model with brms
fit_brms <- brm(
  formula = as.formula(
    paste("Event_time | cens(1 - Event) ~", paste(predictors, collapse = " + "))),
  data = dff,
  family = brmsfamily("cox"),
  prior = set_prior("normal(0, 1)", class = "b"),
  chains = 4,
  iter = 4000,
  cores = parallel::detectCores()
)

# Save model
saveRDS(fit_brms, file = "fit_brms.rds")
