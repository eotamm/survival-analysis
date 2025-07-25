# Survivalanalysis using probabilistic models


# Introduction

Let’s download all the necessary libraries:

``` r
library(tidyverse)
library(brms)
library(tidybayes)
library(ggplot2)
library(TreeSummarizedExperiment)
library(posterior)
library(dplyr)
library(tidyr)
library(stringr)
library(survival)
library(survminer)
library(vegan)
library(mia)
library(bayesboot)
library(IDPSurvival)
library(cmdstanr)
```

# Survival model

## Data manipulation

All data preprocessing was carried out in a separate script (`data.R`) to ensure clarity and modularity in the analysis workflow. The script loads the raw data, applies appropriate transformations, removes invalid samples, and derives variables needed for modeling.

To reduce dimensionality and avoid overfitting in probabilistic survival modeling, only the 5 most prevalent microbial taxa were selected based on their presence across samples. This targeted feature selection improves the robustness and interpretability of the model, especially given the relatively limited sample size.

## Modeling

We construct a Cox model formula using the log-transformed predictors already available in dff. A probabilistic Cox proportional hazards model is then fitted using brm(), with a normal(0, 1) prior on each bacterial coefficient to reflect weak prior knowledge. All of this is implemented in a separate script (`model.R`), and we now load the results.

``` r
# Load model
fit_brms <- readRDS("fit_brms.rds")
```

We take summary of the model.

``` r
# Print model summary
summary(fit_brms)
```

     Family: cox 
      Links: mu = log 
    Formula: Event_time | cens(1 - Event) ~ g_Bacteroides + f_Lachnospiraceae_g + f_Enterobacteriaceae_g + g_Ruminococcus + o_Clostridiales_g 
       Data: df (Number of observations: 150) 
      Draws: 4 chains, each with iter = 4000; warmup = 2000; thin = 1;
             total post-warmup draws = 8000

    Regression Coefficients:
                           Estimate Est.Error l-95% CI u-95% CI Rhat Bulk_ESS
    Intercept                 -0.91      1.16    -3.31     1.24 1.00     8528
    g_Bacteroides             -0.09      0.04    -0.17    -0.00 1.00     8012
    f_Lachnospiraceae_g        0.14      0.10    -0.06     0.33 1.00     6400
    f_Enterobacteriaceae_g     0.08      0.04    -0.00     0.17 1.00     8524
    g_Ruminococcus             0.03      0.06    -0.09     0.14 1.00     7648
    o_Clostridiales_g         -0.07      0.07    -0.20     0.07 1.00     7200
                           Tail_ESS
    Intercept                  6689
    g_Bacteroides              5963
    f_Lachnospiraceae_g        5656
    f_Enterobacteriaceae_g     5755
    g_Ruminococcus             5951
    o_Clostridiales_g          5313

    Draws were sampled using sampling(NUTS). For each parameter, Bulk_ESS
    and Tail_ESS are effective sample size measures, and Rhat is the potential
    scale reduction factor on split chains (at convergence, Rhat = 1).

The model shows good convergence diagnostics. All Rhat values are equal to 1.00, indicating that the chains have converged well. Additionally, both Bulk and Tail Effective Sample Sizes (ESS) are high across all parameters, suggesting that the posterior distributions are well estimated and the sampling was efficient.

Posterior samples of the coefficients are extracted and summarized by computing the median, 2.5%, and 97.5% quantiles, forming the 95% credible interval for each predictor.

``` r
# Extract raw draws as dataframe
draws_df <- as_draws_df(fit_brms)

# Select only b_ variables, excluding intercept
posterior_summary <- draws_df %>%
  select(starts_with("b_")) %>%
  select(-b_Intercept) %>%
  pivot_longer(everything(), names_to = "variable", values_to = "value") %>%
  group_by(variable) %>%
  summarise(
    median = median(value),
    lower = quantile(value, 0.025),
    upper = quantile(value, 0.975)
  ) %>%
  mutate(variable_clean = str_remove(variable, "b_"))


# Add significance column (if 95% CI does not include 0, (log(HR)=0))
posterior_summary <- posterior_summary %>%
  mutate(significant = ifelse(lower > 0 | upper < 0, "yes", "no"))
```

## Visualization

The posterior medians and credible intervals are plotted as hazard ratios. Predictors whose intervals do not cross 1 (indicating statistical relevance) are highlighted in red, and others in black. The plot is sorted by effect size.

``` r
# Plot hazard ratios
ggplot(posterior_summary, aes(y = reorder(variable_clean, median), 
                              x = exp(median), 
                              xmin = exp(lower), 
                              xmax = exp(upper), 
                              color = significant)) +
  geom_point() +
  geom_errorbarh(height = 0.2) +
  geom_vline(xintercept = 1, linetype = "dashed", color = "gray") +
  scale_color_manual(values = c("yes" = "red", "no" = "black")) +
  labs(x = "Hazard Ratio", y = "", title = "Posterior estimates for all predictors (log-transformed)")
```

![](figures/unnamed-chunk-6-1.png)

This figure displays the posterior estimates of hazard ratios for all log-transformed bacterial predictors. Each point represents the median of the posterior distribution, and the horizontal lines show the 95% credible intervals. Predictors highlighted in red are statistically significant, meaning their credible intervals do not include 1.

Hazard ratios below 1 (to the left of the dashed line) are associated with a reduced risk of mortality, while those above 1 suggest an increased risk. In this model, most predictors show no statistically significant association with mortality; only g_Bacteroides demonstrates a significant effect, indicating a potential reduce in risk.

## CLR-transformed modeling

Next, we built a separate survival model using the five most common microbial taxa, as before. Because microbiome data is compositional, we applied a centered log-ratio (CLR) transformation to the relative abundances. This method expresses each taxon’s abundance in relation to the others, helping to reduce misleading correlations and making the results easier to compare across samples.

``` r
# Load model
fit_clr5 <- readRDS("fit_clr5.rds")
```

We take summary of the model.

``` r
# Print model summary
summary(fit_clr5)
```

     Family: cox 
      Links: mu = log 
    Formula: Event_time | cens(1 - Event) ~ g_Bacteroides + f_Lachnospiraceae_g + f_Enterobacteriaceae_g + g_Ruminococcus + o_Clostridiales_g 
       Data: df (Number of observations: 150) 
      Draws: 4 chains, each with iter = 4000; warmup = 2000; thin = 1;
             total post-warmup draws = 8000

    Regression Coefficients:
                           Estimate Est.Error l-95% CI u-95% CI Rhat Bulk_ESS
    Intercept                  0.06      0.20    -0.36     0.45 1.00     6728
    g_Bacteroides             -0.15      0.45    -1.04     0.75 1.00     2275
    f_Lachnospiraceae_g        0.19      0.47    -0.75     1.12 1.00     2442
    f_Enterobacteriaceae_g     0.09      0.45    -0.80     1.00 1.00     2268
    g_Ruminococcus             0.03      0.46    -0.86     0.96 1.00     2311
    o_Clostridiales_g         -0.11      0.46    -1.00     0.82 1.00     2352
                           Tail_ESS
    Intercept                  6684
    g_Bacteroides              3015
    f_Lachnospiraceae_g        3490
    f_Enterobacteriaceae_g     2953
    g_Ruminococcus             3121
    o_Clostridiales_g          3139

    Draws were sampled using sampling(NUTS). For each parameter, Bulk_ESS
    and Tail_ESS are effective sample size measures, and Rhat is the potential
    scale reduction factor on split chains (at convergence, Rhat = 1).

The model shows good convergence diagnostics. All Rhat values are equal to 1.00, indicating that the chains have converged well. Additionally, both Bulk and Tail Effective Sample Sizes (ESS) are high across all parameters, suggesting that the posterior distributions are well estimated and the sampling was efficient.

Posterior samples of the coefficients are extracted and summarized by computing the median, 2.5%, and 97.5% quantiles, forming the 95% credible interval for each predictor.

``` r
# Extract raw draws as dataframe
draws_df_clr <- as_draws_df(fit_clr5)

# Select only b_ variables, excluding intercept
posterior_summary_clr <- draws_df_clr %>%
  select(starts_with("b_")) %>%
  select(-b_Intercept) %>%
  pivot_longer(everything(), names_to = "variable", values_to = "value") %>%
  group_by(variable) %>%
  summarise(
    median = median(value),
    lower = quantile(value, 0.025),
    upper = quantile(value, 0.975)
  ) %>%
  mutate(variable_clean = str_remove(variable, "b_")) %>%
  mutate(significant = ifelse(lower > 0 | upper < 0, "yes", "no"))
```

The posterior medians and credible intervals are plotted as hazard ratios as earlier.

``` r
# Plot
ggplot(posterior_summary_clr, aes(y = reorder(variable_clean, median), 
                              x = exp(median), 
                              xmin = exp(lower), 
                              xmax = exp(upper), 
                              color = significant)) +
  geom_point() +
  geom_errorbarh(height = 0.2) +
  geom_vline(xintercept = 1, linetype = "dashed", color = "gray") +
  scale_color_manual(values = c("yes" = "red", "no" = "black")) +
  labs(x = "Hazard Ratio", y = "", title = "Posterior estimates for all predictors (CLR-transformed)")
```

![](figures/unnamed-chunk-10-1.png)

Both modeling approaches identified g_Bacteroides as potentially associated with survival. In the log-transformed model, its hazard ratio was significantly below 1, suggesting a protective effect. However, this association was not statistically significant in the CLR-transformed model. Overall, results were directionally consistent across models, but only the log-transformed version yielded strong statistical evidence.

This discrepancy is most likely explained by the difference in interpretation: the log-transformed model estimates each taxon’s effect more independently, while the CLR-transformed model evaluates each taxon relative to the others. This relative structure introduces greater variability and typically results in wider credible intervals.

## LRA-based modeling

Lastly, we built a survival model using pairwise log-ratios (LRA) between the five most common taxa. This approach is well-suited for microbiome data, as it looks at how the balance between two taxa relates to the outcome, rather than their individual abundance. It also handles the relative nature of the data and avoids scale issues.

``` r
# Load LRA model
fit_lra5 <- readRDS("fit_lra5.rds")
```

We take summary of the model.

``` r
# Print model summary
summary(fit_lra5)
```

     Family: cox 
      Links: mu = log 
    Formula: Event_time | cens(1 - Event) ~ diff_g_Bacteroidesf_Lachnospiraceae_g + diff_g_Bacteroidesf_Enterobacteriaceae_g + diff_g_Bacteroidesg_Ruminococcus + diff_g_Bacteroideso_Clostridiales_g + diff_f_Lachnospiraceae_g_f_Enterobacteriaceae_g + diff_f_Lachnospiraceae_g_g_Ruminococcus + diff_f_Lachnospiraceae_g_o_Clostridiales_g + diff_f_Enterobacteriaceae_g_g_Ruminococcus + diff_f_Enterobacteriaceae_g_o_Clostridiales_g + diff_g_Ruminococcuso_Clostridiales_g 
       Data: df (Number of observations: 150) 
      Draws: 4 chains, each with iter = 4000; warmup = 2000; thin = 1;
             total post-warmup draws = 8000

    Regression Coefficients:
                                                    Estimate Est.Error l-95% CI
    Intercept                                          -0.12      0.14    -0.39
    diff_g_Bacteroidesf_Lachnospiraceae_g              -0.15      0.76    -1.62
    diff_g_Bacteroidesf_Enterobacteriaceae_g           -0.11      0.82    -1.70
    diff_g_Bacteroidesg_Ruminococcus                   -0.07      0.83    -1.71
    diff_g_Bacteroideso_Clostridiales_g                -0.05      0.80    -1.60
    diff_f_Lachnospiraceae_g_f_Enterobacteriaceae_g    -0.02      0.85    -1.70
    diff_f_Lachnospiraceae_g_g_Ruminococcus             0.03      0.66    -1.28
    diff_f_Lachnospiraceae_g_o_Clostridiales_g          0.10      0.56    -1.01
    diff_f_Enterobacteriaceae_g_g_Ruminococcus          0.06      0.80    -1.49
    diff_f_Enterobacteriaceae_g_o_Clostridiales_g       0.07      0.87    -1.66
    diff_g_Ruminococcuso_Clostridiales_g                0.06      0.78    -1.48
                                                    u-95% CI Rhat Bulk_ESS Tail_ESS
    Intercept                                           0.15 1.00     9496     6383
    diff_g_Bacteroidesf_Lachnospiraceae_g               1.35 1.00     6429     5505
    diff_g_Bacteroidesf_Enterobacteriaceae_g            1.52 1.00     7156     5297
    diff_g_Bacteroidesg_Ruminococcus                    1.56 1.00     7478     4988
    diff_g_Bacteroideso_Clostridiales_g                 1.48 1.00     6158     5417
    diff_f_Lachnospiraceae_g_f_Enterobacteriaceae_g     1.64 1.00     7611     5552
    diff_f_Lachnospiraceae_g_g_Ruminococcus             1.34 1.00     5626     5627
    diff_f_Lachnospiraceae_g_o_Clostridiales_g          1.20 1.00     4921     5116
    diff_f_Enterobacteriaceae_g_g_Ruminococcus          1.65 1.00     7026     5434
    diff_f_Enterobacteriaceae_g_o_Clostridiales_g       1.78 1.00     7475     4720
    diff_g_Ruminococcuso_Clostridiales_g                1.57 1.00     5177     5017

    Draws were sampled using sampling(NUTS). For each parameter, Bulk_ESS
    and Tail_ESS are effective sample size measures, and Rhat is the potential
    scale reduction factor on split chains (at convergence, Rhat = 1).

The model shows good convergence diagnostics. All Rhat values are equal to 1.00, indicating that the chains have converged well. Bulk and Tail Effective Sample Sizes are also high, suggesting that the posterior distributions are well estimated and sampling was efficient.

Posterior samples of the log-ratio coefficients are extracted and summarized using the median and 95% credible intervals.

``` r
# Extract posterior draws
draws_df_lra <- as_draws_df(fit_lra5)

# Summarize LRA coefficients and clean names
posterior_summary_lra <- draws_df_lra %>%
  select(starts_with("b_")) %>%
  select(-b_Intercept) %>%
  pivot_longer(cols = everything(), names_to = "variable", values_to = "value") %>%
  group_by(variable) %>%
  summarise(
    median = median(value),
    lower = quantile(value, 0.025),
    upper = quantile(value, 0.975)
  ) %>%
  mutate(
    variable_clean = str_remove(variable, "b_") %>% 
                     str_replace("^diff_", "logratio_"),
    significant = ifelse(lower > 0 | upper < 0, "yes", "no")
  )
```

The log-ratio coefficients are plotted below as hazard ratios.

``` r
# Plot
ggplot(posterior_summary_lra, aes(
  y = reorder(variable_clean, median),
  x = exp(median),
  xmin = exp(lower),
  xmax = exp(upper),
  color = significant)) +
  geom_point() +
  geom_errorbarh(height = 0.2) +
  geom_vline(xintercept = 1, linetype = "dashed", color = "gray") +
  scale_color_manual(values = c("yes" = "red", "no" = "black")) +
  labs(
    x = "Hazard Ratio",
    y = "",
    title = "Posterior estimates for all predictors (LRA-based)"
  )
```

![](figures/unnamed-chunk-14-1.png)

None of the log-ratio predictors reached statistical significance, as all 95% credible intervals included 1. This likely reflects limited statistical power due to the relatively small sample size and the increased uncertainty that comes from modeling many pairwise ratios. Additionally, the structure of the LRA model plays a role: it estimates the effect of one taxon always in relation to another, which spreads the information across multiple comparisons and tends to produce wider credible intervals than models focusing on individual taxa.

# Survival curves

## Kaplan-Meier curve

First classical Kaplan–Meier estimator is applied to provide a robust and interpretable nonparametric summary of cumulative mortality over time.

``` r
# Fit Kaplan-Meier survival model
surv_fit <- survfit(Surv(Event_time, Event) ~ 1, data = df)

# Plot cumulative mortality curve
ggsurvplot(
  surv_fit,
  data = df,
  conf.int = TRUE,
  fun = "event",  
  palette = "blue",
  xlab = "Time (years)",
  ylab = "Cumulative mortality (%)",
  title = "Overall cumulative mortality"
)
```

![](figures/unnamed-chunk-15-1.png)

This plot shows the overall cumulative mortality over a 10-year follow-up period. The curve represents the estimated probability of death over time in the entire cohort. The shaded area around the curve indicates the 95% confidence interval, and the small vertical ticks represent censored observations. By the end of the follow-up, cumulative mortality reaches approximately 60%, with a fairly steady increase over time.

## Probabilistic Survival Curve with IDPSurvival

We now estimate a probabilistic survival curve using the `isurvfit()` function from the `IDPSurvival` package. This method is based on the Imprecise Dirichlet Process, which provides robust survival estimates with uncertainty bounds, especially useful for small or uncertain datasets.

``` r
# Create the Surv object
surv_obj <- Surv(time = df$Event_time, event = df$Event)

# Estimate the IDP survival curve
fit <- isurvfit(surv_obj ~ 1, data = df, s = 1, 
                conf.type = "exact", nsamples = 2000, display = FALSE)

# Plot the survival curve
plot(fit)
title("Probabilistic Survival Curve (IDP)")
mtext("Time (years)", side = 1, line = 2)
mtext("Survival probability", side = 2, line = 2)
legend('bottomleft', c("Lower expectation",
          "Upper expectation","Confidence intervals"), lty=c(1,1,2),lwd=c(1,2,1))
```

![](figures/unnamed-chunk-16-1.png)

The parameter `s = 1` controls the strength of the prior. A higher `s` reflects more confidence in the prior and results in narrower uncertainty bands. Conversely, lower values of `s` allow for more imprecision, widening the interval between lower and upper expectations. The parameter `nsamples = 2000` specifies the number of posterior samples used to construct the credible intervals. A larger value gives a smoother and more stable estimate of the uncertainty region. The `conf.type = "exact"` option determines how the uncertainty bounds are calculated. When set to `"exact"`, the intervals are computed using the full posterior distribution.

## Comparison of Classical and Probabilistic Survival Estimation Methods

To better understand the differences between classical and probabilistic survival estimates, we compare the Kaplan–Meier curve with the survival curve derived from the Imprecise Dirichlet Process approach.

``` r
# Fit the classical Kaplan–Meier survival model
km_fit <- survfit(Surv(Event_time, Event) ~ 1, data = df)

# Fit the probabilistic IDP survival model (s = 1)
idp_fit <- isurvfit(Surv(Event_time, Event) ~ 1, data = df, s = 1,
                    conf.type = "exact", nsamples = 2000, display = FALSE)

# Prepare Kaplan–Meier survival curve data
km_df <- data.frame(
  time = km_fit$time,
  surv = km_fit$surv,
  lower = km_fit$lower,
  upper = km_fit$upper,
  method = "Kaplan–Meier"
)

# Prepare IDP survival curve data using midpoint between lower and upper expectations
idp_df <- data.frame(
  time = idp_fit$time,
  surv = (idp_fit$survUP + idp_fit$survLOW) / 2,
  lower = idp_fit$lower,
  upper = idp_fit$upper,
  method = "IDP (s = 1)"
)

# Combine the datasets for plotting
combined_df <- bind_rows(km_df, idp_df)

# Plot the survival curves
ggplot(combined_df, aes(x = time, y = surv, color = method, fill = method)) +
  geom_step(linewidth = 1.2) +
  geom_ribbon(aes(ymin = lower, ymax = upper), alpha = 0.2, color = NA) +
  labs(
    title = "Survival Curve Comparison: Kaplan–Meier vs. IDP (s = 1)",
    x = "Time",
    y = "Survival probability",
    color = "Method",
    fill = "Method"
  ) +
  theme_minimal(base_size = 14)
```

![](figures/unnamed-chunk-17-1.png)

The point estimates from the IDP and Kaplan–Meier methods are very similar, but the uncertainty band from the IDP model lies slightly lower than uncertainty band from the Kaplan-Meier.

# Comparison of CLR, rCLR, and Log-Relative-Abundance Transformations Using Microbial Log-Ratio Features in Survival Analysis

To evaluate the impact of compositional preprocessing choices, we compare three common transformations: centered log-ratio (CLR), robust CLR (rCLR), and log-relative-abundance (LRA). Using the 20 most prevalent microbial taxa, we compute all pairwise log(x/y) ratios as input features for survival modeling.

The log(x/y) transformation expresses the relative abundance of one taxon compared to another on a log scale, making the data more statistically appropriate for modeling and removing compositional constraints. These features are used in Bayesian Cox models, and the predictive performance of all three models is compared using leave-one-out cross-validation (LOO).

As a baseline reference, we also include a model using log-transformed absolute abundances (without any ratio transformation) of the same 20 taxa. This provides a point of comparison to assess whether ratio-based features offer an advantage over simpler log-abundance inputs in survival prediction.

``` r
# Count how many observations have pareto_k > 0.7 per model
sapply(loos, function(x) sum(x$diagnostics$pareto_k > 0.7))
```

      lra   clr  rclr abund 
        1     0     0     0 

``` r
# Print the results
print(loo_comp)
```

                    elpd_diff se_diff
    fit_lra_ratios   0.0       0.0   
    fit_clr         -1.2       0.8   
    fit_rclr        -1.5       0.8   
    fit_abund_top20 -3.3       1.1   

After applying moment matching to improve the reliability of the leave-one-out estimates, all Pareto $k$ values were below the diagnostic threshold of $0.7$, except for a single observation in the LRA-based model. This indicates that the LOO estimates are reliable and can be used for model comparison.

The LRA model showed the best predictive performance. Both the CLR and rCLR models performed slightly worse, but the differences were small and well within one standard error, suggesting no meaningful drop in accuracy. This indicates that the specific choice of log-ratio transformation has little impact on predictive performance in this setting.

In contrast, the model based on raw log-transformed abundances (without log-ratios) performed notably worse. Its predictive performance was clearly lower than the LRA model, with a difference that exceeded one standard error. Indicating a meaningful decline when ratio-based features are not used.

The CLR- and rCLR-based models included 20 transformed taxa as predictors, whereas the LRA-based model incorporated 190 pairwise log-ratio features derived from the same 20 taxa. This led to a substantial increase in dimensionality for the LRA model, although model fitting remained relatively fast and stable.
