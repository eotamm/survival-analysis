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
library(patchwork)
library(pROC)
library(RColorBrewer) 
library(ranger)
library(timeROC)
library(xgboost)
library(matrixStats)
library(keras)
library(tensorflow)
library(catboost)
library(reticulate)
```

## Data manipulation

All data preprocessing was carried out in a separate script (`data.R`) to ensure clarity and modularity in the analysis workflow. The script loads the raw data, applies appropriate transformations, removes invalid samples, and derives variables needed for modeling.

## Data transformations

To evaluate the effect of different data preprocessing strategies on survival model performance, several commonly used microbial feature transformations were applied:

-   **CLR (Centered Log-Ratio)**  
    Transforms counts to log-ratios relative to the geometric mean of all features. Suitable for compositional data.

-   **rCLR (Robust CLR)**  
    A robust version of CLR that downweights outliers. Helps reduce the influence of extreme values.

-   **log-TSS (Log-transformed Total Sum Scaling)**  
    Counts are converted to relative abundances and then log-transformed to stabilize variance and reduce skewness.

-   **LRA (Log-Ratio Analysis)**  
    Computes pairwise log-ratios between features after log-transformation. Emphasizes relative differences.

-   **PA (Presence/Absence)**  
    Binary encoding of whether each taxon is detected. Ignores abundance but captures occurrence. The detection threshold for PA was relative (0.1% TSS)

-   **TSS (Total Sum Scaling)**  
    Normalizes counts to relative abundances by dividing each sample by its total count.

-   **Arcsin (Arcsin-sqrt transformation of TSS)**  
    Applies arcsin square-root transformation to relative abundances, useful for proportions near 0 or 1.

-   **ALR (Additive Log-Ratio)**  
    Computes log-ratios relative to a selected reference taxon. Sensitive to choice of reference. In this analysis reference taxon is g\_\_Sutterella.

Each transformation highlights different aspects of the microbial composition and is evaluated independently in downstream survival models.

# Optimizing Feature Count for Survival Prediction

In this analysis, we systematically evaluate how the number of input features affects the predictive performance of Bayesian survival models. Starting from a single feature and increasing the count in powers of two up to a predefined maximum, we compare model performance using the concordance index (C-index) across multiple data transformations. The models are fit using a weakly informative normal prior for the regression coefficients, which encourages regularization while allowing flexibility. This approach helps identify a potential “sweet spot” where the model achieves high performance with a minimal and interpretable feature set, before diminishing returns or overfitting may occur.

To analyze how predictive performance varies across transformations and different numbers of selected features, we reorganize the model results into a long-format data frame. This structure facilitates visualization and comparison across transformation methods and feature subset sizes. Each row represents a single model fit, annotated with the number of features used (`N`), the transformation method applied, and the number of taxa (`ntaxa`) considered in the subset.

``` r
# Reshape results into long format
res_brm_df <- purrr::imap_dfr(res_brm_result, function(method_list, ntaxa_name) {
  ntaxa_val <- as.numeric(gsub("ntaxa_", "", ntaxa_name))
  purrr::imap_dfr(method_list, function(df, method_name) {
    df %>%
      mutate(
        ntaxa = ntaxa_val,
        method = method_name
      )
  })
})

# Visualize C-index values across transformation methods and feature counts
ggplot(res_brm_df, aes(x = N, y = C_index, color = method)) +
  geom_line(linewidth = 1) +
  geom_point(size = 2) +
  facet_wrap(~ ntaxa, labeller = label_both, scales = "free_y") +
  labs(
    title = "C-index of Bayesian Cox models across transformation methods",
    x = "Number of selected features (N)",
    y = "C-index (posterior mean)",
    color = "Transformation"
  ) +
  scale_x_log10() +
  theme_minimal(base_size = 14)
```

![](figures/unnamed-chunk-4-1.png)

Model performance, measured by the C-index, improves as more features are included but levels off around 20 features. After that, additional features offer no benefit for most of the models. This suggests that using about 20 top features is a good balance between simplicity and performance.

# Survival model

In survival analysis, Cox models are fitted using microbial features as predictors. To reduce dimensionality and avoid overfitting, the 20 most informative features are selected for each data transformation based on univariate Cox regression. This improves model performance and interpretability.

Bayesian Cox models are then fitted using brm() with two alternative priors on the coefficients: a weakly informative normal(0, 1) prior and a sparsity-inducing horseshoe prior. The modeling is implemented in model.R, and results are loaded for analysis.

Posterior hazard ratio estimates are computed and visualized only for selected transformations: log-transformed relative abundance, centered log-ratio (CLR), and log-ratio analysis (LRA). This allows clearer comparison of interpretability and the direction of feature effects.

## Relative log-abundance modeling

The log-transformed total sum scaling (log-TSS) approach is a commonly used normalization method where relative abundances are scaled and log-transformed to stabilize variance and reduce skewness. Both normal and horseshoe priors were used when fitting the Bayesian Cox models to compare how the choice of prior affects the estimated hazard ratios. This transformation makes microbial data more suitable for downstream regression modeling.

We take summaries of the models.

``` r
# Print model summary
summary(fit_logtss_normal)
```

     Family: cox 
      Links: mu = log 
    Formula: Event_time | cens(1 - Event) ~ g_Streptococcus + g_Bacteroides + g_Roseburia + g_Dialister + g_Lachnospira + f_Gemellaceae_g + g_Faecalibacterium + g_Actinomyces + g_Granulicatella + g_Prevotella + g_Conchiformibius + f_Rikenellaceae_g + o_Clostridiales_g + o_Lactobacillales_g + g_Holdemania + g_Coprococcus + g_Parabacteroides + g_Aggregatibacter + g_Odoribacter + g_Veillonella 
       Data: df (Number of observations: 150) 
      Draws: 4 chains, each with iter = 4000; warmup = 2000; thin = 1;
             total post-warmup draws = 8000

    Regression Coefficients:
                        Estimate Est.Error l-95% CI u-95% CI Rhat Bulk_ESS Tail_ESS
    Intercept               0.74      1.28    -1.77     3.21 1.00     8489     6899
    g_Streptococcus         0.11      0.08    -0.05     0.26 1.00     6753     5591
    g_Bacteroides          -0.09      0.06    -0.20     0.02 1.00     7507     6034
    g_Roseburia            -0.10      0.05    -0.21    -0.00 1.00     9367     6370
    g_Dialister             0.09      0.04     0.01     0.17 1.00     8825     5580
    g_Lachnospira          -0.13      0.07    -0.26     0.00 1.00     9276     6404
    f_Gemellaceae_g         0.09      0.08    -0.06     0.25 1.00     6456     6174
    g_Faecalibacterium      0.03      0.05    -0.07     0.13 1.00     8041     6139
    g_Actinomyces           0.09      0.09    -0.09     0.26 1.00     7621     6180
    g_Granulicatella       -0.08      0.09    -0.26     0.09 1.00     6902     5891
    g_Prevotella           -0.15      0.06    -0.27    -0.05 1.00     8184     6391
    g_Conchiformibius       0.09      0.05    -0.02     0.19 1.00     8717     5998
    f_Rikenellaceae_g      -0.01      0.05    -0.10     0.09 1.00     8058     5832
    o_Clostridiales_g      -0.00      0.07    -0.15     0.14 1.00     7690     5923
    o_Lactobacillales_g    -0.02      0.05    -0.12     0.08 1.00     8299     6420
    g_Holdemania            0.02      0.10    -0.17     0.20 1.00     8351     6462
    g_Coprococcus           0.03      0.06    -0.08     0.14 1.00     7843     5684
    g_Parabacteroides      -0.02      0.05    -0.11     0.07 1.00     8084     6405
    g_Aggregatibacter       0.02      0.05    -0.09     0.12 1.00     7840     5800
    g_Odoribacter           0.01      0.06    -0.11     0.12 1.00     7338     6205
    g_Veillonella          -0.02      0.06    -0.14     0.10 1.00     7990     5609

    Draws were sampled using sampling(NUTS). For each parameter, Bulk_ESS
    and Tail_ESS are effective sample size measures, and Rhat is the potential
    scale reduction factor on split chains (at convergence, Rhat = 1).

``` r
summary(fit_logtss_hs)
```

     Family: cox 
      Links: mu = log 
    Formula: Event_time | cens(1 - Event) ~ g_Streptococcus + g_Bacteroides + g_Roseburia + g_Dialister + g_Lachnospira + f_Gemellaceae_g + g_Faecalibacterium + g_Actinomyces + g_Granulicatella + g_Prevotella + g_Conchiformibius + f_Rikenellaceae_g + o_Clostridiales_g + o_Lactobacillales_g + g_Holdemania + g_Coprococcus + g_Parabacteroides + g_Aggregatibacter + g_Odoribacter + g_Veillonella 
       Data: df (Number of observations: 150) 
      Draws: 4 chains, each with iter = 4000; warmup = 2000; thin = 1;
             total post-warmup draws = 8000

    Regression Coefficients:
                        Estimate Est.Error l-95% CI u-95% CI Rhat Bulk_ESS Tail_ESS
    Intercept               0.43      0.98    -1.46     2.35 1.00     5770     6720
    g_Streptococcus         0.07      0.07    -0.01     0.21 1.00     3842     5483
    g_Bacteroides          -0.03      0.04    -0.13     0.01 1.00     4627     7080
    g_Roseburia            -0.05      0.05    -0.16     0.01 1.00     3934     5608
    g_Dialister             0.07      0.04    -0.00     0.15 1.00     3455     3008
    g_Lachnospira          -0.05      0.06    -0.18     0.01 1.00     4368     6403
    f_Gemellaceae_g         0.03      0.04    -0.03     0.14 1.00     6073     7179
    g_Faecalibacterium     -0.01      0.02    -0.06     0.04 1.00     7794     7321
    g_Actinomyces           0.03      0.05    -0.03     0.16 1.00     6247     7244
    g_Granulicatella        0.01      0.03    -0.06     0.09 1.00     8290     6037
    g_Prevotella           -0.09      0.06    -0.21     0.00 1.00     3734     2301
    g_Conchiformibius       0.03      0.04    -0.02     0.14 1.00     5283     5437
    f_Rikenellaceae_g      -0.01      0.02    -0.07     0.03 1.00     8095     7241
    o_Clostridiales_g      -0.01      0.03    -0.09     0.05 1.00     7721     6858
    o_Lactobacillales_g     0.00      0.02    -0.04     0.06 1.00     9282     7146
    g_Holdemania           -0.00      0.03    -0.08     0.06 1.00     9020     6834
    g_Coprococcus          -0.00      0.02    -0.05     0.05 1.00     8959     7552
    g_Parabacteroides      -0.01      0.02    -0.06     0.03 1.00     8687     7671
    g_Aggregatibacter       0.01      0.03    -0.03     0.08 1.00     7240     7523
    g_Odoribacter          -0.00      0.02    -0.06     0.04 1.00     8668     7394
    g_Veillonella           0.00      0.02    -0.05     0.05 1.00     8756     6894

    Draws were sampled using sampling(NUTS). For each parameter, Bulk_ESS
    and Tail_ESS are effective sample size measures, and Rhat is the potential
    scale reduction factor on split chains (at convergence, Rhat = 1).

The models show good convergence diagnostics. All Rhat values are equal to 1.00, indicating that the chains have converged well. Additionally, both Bulk and Tail Effective Sample Sizes (ESS) are high across all parameters, suggesting that the posterior distributions are well estimated and the sampling was efficient.

Posterior samples of the coefficients are extracted and summarized by computing the median, 2.5%, and 97.5% quantiles, forming the 95% credible interval for each predictor. The posterior medians and credible intervals are plotted as hazard ratios. Predictors whose intervals do not cross 1 (indicating statistical relevance) are highlighted in red, and others in black. The plot is sorted by effect size.

``` r
# Plot
p1 <- plot_top_hr(fit_logtss_normal,
                  title = "logTSS – Normal prior",
                  top_k = 20)
p2 <- plot_top_hr(fit_logtss_hs,
                  title = "logTSS – Horseshoe prior",
                  top_k = 20)
p1 + p2
```

![](figures/unnamed-chunk-6-1.png)

In the logTSS-transformed models, the normal prior identified several taxa with credible intervals excluding 1, including g_Dialister (higher risk) and g_Roseburia and g_Prevotella (lower risk). Under the horseshoe prior, all coefficients were shrunk toward zero and no taxa remained significant, reflecting the stronger regularization of this prior.

## CLR-transformed modeling

In another model, we used the top microbial features based on their individual associations with the outcome. Before fitting the model, we applied a centered log-ratio (CLR) transformation to express each microbe’s level in relation to the others. Both models were fitted with different priors (normal and horseshoe) to compare the effect of prior choice on coefficient estimates. This helps make the results easier to interpret and reduces the chance of misleading correlations.

We take summaries of the models.

``` r
# Print model summary
summary(fit_clr_normal)
```

     Family: cox 
      Links: mu = log 
    Formula: Event_time | cens(1 - Event) ~ g_Streptococcus + g_Bacteroides + g_Roseburia + g_Dialister + g_Lachnospira + f_Gemellaceae_g + g_Actinomyces + g_Faecalibacterium + g_Granulicatella + o_Clostridiales_g + g_Prevotella + g_Conchiformibius + f_Rikenellaceae_g + g_Coprococcus + o_Lactobacillales_g + g_Holdemania + g_Parabacteroides + g_Aggregatibacter + g_Veillonella + g_Odoribacter 
       Data: df (Number of observations: 150) 
      Draws: 4 chains, each with iter = 4000; warmup = 2000; thin = 1;
             total post-warmup draws = 8000

    Regression Coefficients:
                        Estimate Est.Error l-95% CI u-95% CI Rhat Bulk_ESS Tail_ESS
    Intercept               0.40      0.57    -0.74     1.49 1.00    10460     6247
    g_Streptococcus         0.15      0.12    -0.09     0.37 1.00     7474     5864
    g_Bacteroides          -0.15      0.09    -0.32     0.02 1.00     7936     6173
    g_Roseburia            -0.16      0.08    -0.31    -0.02 1.00    10460     6175
    g_Dialister             0.11      0.06    -0.01     0.23 1.00     9379     6385
    g_Lachnospira          -0.20      0.10    -0.40    -0.01 1.00    10385     5938
    f_Gemellaceae_g         0.11      0.11    -0.10     0.33 1.00     7888     6744
    g_Actinomyces           0.08      0.14    -0.19     0.34 1.00     7854     6757
    g_Faecalibacterium      0.03      0.07    -0.11     0.17 1.00     8756     6207
    g_Granulicatella       -0.13      0.13    -0.39     0.11 1.00     7825     6249
    o_Clostridiales_g      -0.04      0.12    -0.27     0.18 1.00     8256     6489
    g_Prevotella           -0.23      0.08    -0.40    -0.07 1.00     8420     5926
    g_Conchiformibius       0.11      0.08    -0.05     0.25 1.00     8138     6397
    f_Rikenellaceae_g      -0.03      0.07    -0.17     0.12 1.00     8800     6540
    g_Coprococcus           0.02      0.09    -0.16     0.19 1.00     7678     5921
    o_Lactobacillales_g    -0.05      0.08    -0.21     0.10 1.00     7041     5984
    g_Holdemania            0.02      0.13    -0.24     0.27 1.00     9082     6577
    g_Parabacteroides      -0.04      0.07    -0.17     0.10 1.00     8964     6251
    g_Aggregatibacter       0.01      0.08    -0.14     0.16 1.00     8248     6744
    g_Veillonella          -0.04      0.08    -0.20     0.13 1.00     9204     6066
    g_Odoribacter           0.00      0.09    -0.17     0.17 1.00     7715     6206

    Draws were sampled using sampling(NUTS). For each parameter, Bulk_ESS
    and Tail_ESS are effective sample size measures, and Rhat is the potential
    scale reduction factor on split chains (at convergence, Rhat = 1).

``` r
summary(fit_clr_hs)
```

     Family: cox 
      Links: mu = log 
    Formula: Event_time | cens(1 - Event) ~ g_Streptococcus + g_Bacteroides + g_Roseburia + g_Dialister + g_Lachnospira + f_Gemellaceae_g + g_Actinomyces + g_Faecalibacterium + g_Granulicatella + o_Clostridiales_g + g_Prevotella + g_Conchiformibius + f_Rikenellaceae_g + g_Coprococcus + o_Lactobacillales_g + g_Holdemania + g_Parabacteroides + g_Aggregatibacter + g_Veillonella + g_Odoribacter 
       Data: df (Number of observations: 150) 
      Draws: 4 chains, each with iter = 4000; warmup = 2000; thin = 1;
             total post-warmup draws = 8000

    Regression Coefficients:
                        Estimate Est.Error l-95% CI u-95% CI Rhat Bulk_ESS Tail_ESS
    Intercept               0.29      0.48    -0.64     1.24 1.00     5968     6458
    g_Streptococcus         0.08      0.09    -0.02     0.29 1.00     5222     7235
    g_Bacteroides          -0.07      0.08    -0.24     0.02 1.00     4455     7300
    g_Roseburia            -0.09      0.08    -0.26     0.01 1.00     4147     5499
    g_Dialister             0.09      0.06    -0.01     0.21 1.00     4577     5379
    g_Lachnospira          -0.10      0.10    -0.32     0.02 1.00     4266     6138
    f_Gemellaceae_g         0.03      0.06    -0.04     0.19 1.00     7232     7258
    g_Actinomyces           0.03      0.06    -0.06     0.21 1.00     8359     7207
    g_Faecalibacterium     -0.01      0.03    -0.10     0.05 1.00     9085     7372
    g_Granulicatella        0.01      0.05    -0.08     0.12 1.00     9827     7133
    o_Clostridiales_g      -0.02      0.05    -0.15     0.06 1.00     8164     6334
    g_Prevotella           -0.15      0.09    -0.33     0.00 1.00     3922     2905
    g_Conchiformibius       0.03      0.05    -0.03     0.16 1.00     6316     7486
    f_Rikenellaceae_g      -0.02      0.04    -0.12     0.04 1.00     8895     7248
    g_Coprococcus          -0.00      0.03    -0.08     0.07 1.00     9778     7581
    o_Lactobacillales_g     0.00      0.03    -0.06     0.08 1.00     9601     7385
    g_Holdemania           -0.01      0.05    -0.13     0.07 1.00    10109     7450
    g_Parabacteroides      -0.01      0.03    -0.10     0.04 1.00     8841     7577
    g_Aggregatibacter       0.01      0.03    -0.04     0.10 1.00     8601     7160
    g_Veillonella           0.00      0.03    -0.07     0.07 1.00     9756     7675
    g_Odoribacter          -0.01      0.03    -0.10     0.06 1.00    10190     7430

    Draws were sampled using sampling(NUTS). For each parameter, Bulk_ESS
    and Tail_ESS are effective sample size measures, and Rhat is the potential
    scale reduction factor on split chains (at convergence, Rhat = 1).

The models show good convergence diagnostics. All Rhat values are equal to 1.00, indicating that the chains have converged well. Additionally, both Bulk and Tail Effective Sample Sizes (ESS) are high across all parameters, suggesting that the posterior distributions are well estimated and the sampling was efficient.

The posterior medians and credible intervals are plotted as hazard ratios as earlier.

``` r
# Plot
p1 <- plot_top_hr(fit_clr_normal,
                  title = "clr – Normal prior",
                  top_k = 20)
p2 <- plot_top_hr(fit_clr_hs,
                  title = "clr – Horseshoe prior",
                  top_k = 20)
p1 + p2
```

![](figures/unnamed-chunk-8-1.png)

In the CLR-transformed models, the normal prior identified g_Roseburia, g_Lachnospira, and g_Prevotella as significantly associated with lower risk. With the horseshoe prior, all effects were shrunk toward zero and no taxa remained significant, indicating stronger regularization.

## LRA-based modeling

We also built survival models using pairwise log-ratios between selected microbial taxa. Both models were fitted with different priors (normal and horseshoe) to compare the effect of prior choice on coefficient estimates. This method focuses on the balance between microbes instead of looking at their individual levels. It helps reduce misleading associations and gives a clearer picture of how microbes relate to the outcome relative to one another.

We take summaries of the models.

``` r
# Print model summary
summary(fit_lra_normal)
```

     Family: cox 
      Links: mu = log 
    Formula: Event_time | cens(1 - Event) ~ g_Lachnospirag_Streptococcus + g_Dialisterg_Roseburia + g_Prevotellag_Streptococcus + g_Roseburiag_Streptococcus + g_Dialistero_Clostridiales_g + g_Dialisterg_Lachnospira + g_Bacteroidesg_Streptococcus + o_Clostridiales_g_g_Streptococcus + f_Gemellaceae_g_g_Lachnospira + g_Bacteroidesf_Gemellaceae_g + g_Faecalibacteriumg_Streptococcus + g_Streptococcusg_Holdemania + g_Dialisterg_Prevotella + g_Bacteroidesg_Dialister + f_Gemellaceae_g_g_Prevotella + g_Faecalibacteriumg_Dialister + g_Bacteroidesg_Conchiformibius + g_Actinomycesg_Roseburia + g_Bacteroidesg_Granulicatella + f_Peptostreptococcaceae_g_g_Streptococcus 
       Data: df (Number of observations: 150) 
      Draws: 4 chains, each with iter = 4000; warmup = 2000; thin = 1;
             total post-warmup draws = 8000

    Regression Coefficients:
                                              Estimate Est.Error l-95% CI u-95% CI
    Intercept                                     0.47      0.54    -0.61     1.50
    g_Lachnospirag_Streptococcus                 -0.02      0.72    -1.44     1.40
    g_Dialisterg_Roseburia                       -0.02      0.64    -1.28     1.20
    g_Prevotellag_Streptococcus                  -0.05      0.73    -1.47     1.38
    g_Roseburiag_Streptococcus                    0.05      0.64    -1.21     1.29
    g_Dialistero_Clostridiales_g                 -0.02      0.64    -1.28     1.23
    g_Dialisterg_Lachnospira                      0.03      0.72    -1.41     1.48
    g_Bacteroidesg_Streptococcus                 -0.06      0.74    -1.52     1.37
    o_Clostridiales_g_g_Streptococcus             0.02      0.64    -1.24     1.29
    f_Gemellaceae_g_g_Lachnospira                 0.03      0.68    -1.29     1.36
    g_Bacteroidesf_Gemellaceae_g                 -0.02      0.67    -1.31     1.32
    g_Faecalibacteriumg_Streptococcus             0.01      0.65    -1.27     1.28
    g_Streptococcusg_Holdemania                  -0.03      0.07    -0.17     0.11
    g_Dialisterg_Prevotella                       0.07      0.73    -1.35     1.49
    g_Bacteroidesg_Dialister                     -0.03      0.73    -1.48     1.41
    f_Gemellaceae_g_g_Prevotella                  0.05      0.67    -1.28     1.37
    g_Faecalibacteriumg_Dialister                 0.00      0.65    -1.26     1.29
    g_Bacteroidesg_Conchiformibius               -0.08      0.04    -0.16     0.01
    g_Actinomycesg_Roseburia                      0.14      0.08    -0.03     0.30
    g_Bacteroidesg_Granulicatella                 0.09      0.08    -0.08     0.26
    f_Peptostreptococcaceae_g_g_Streptococcus    -0.11      0.05    -0.21    -0.01
                                              Rhat Bulk_ESS Tail_ESS
    Intercept                                 1.00     9189     7300
    g_Lachnospirag_Streptococcus              1.00     8028     5442
    g_Dialisterg_Roseburia                    1.00     6835     4961
    g_Prevotellag_Streptococcus               1.00     8230     5462
    g_Roseburiag_Streptococcus                1.00     6905     5227
    g_Dialistero_Clostridiales_g              1.00     6773     5306
    g_Dialisterg_Lachnospira                  1.00     8013     6195
    g_Bacteroidesg_Streptococcus              1.00     7368     5397
    o_Clostridiales_g_g_Streptococcus         1.00     6782     5398
    f_Gemellaceae_g_g_Lachnospira             1.00     7252     6046
    g_Bacteroidesf_Gemellaceae_g              1.00     7579     5938
    g_Faecalibacteriumg_Streptococcus         1.00     6270     5139
    g_Streptococcusg_Holdemania               1.00    11795     6053
    g_Dialisterg_Prevotella                   1.00     8347     5485
    g_Bacteroidesg_Dialister                  1.00     7929     5481
    f_Gemellaceae_g_g_Prevotella              1.00     7504     5464
    g_Faecalibacteriumg_Dialister             1.00     6301     5272
    g_Bacteroidesg_Conchiformibius            1.00    11819     5901
    g_Actinomycesg_Roseburia                  1.00    10992     6212
    g_Bacteroidesg_Granulicatella             1.00    11343     5956
    f_Peptostreptococcaceae_g_g_Streptococcus 1.00    12437     6200

    Draws were sampled using sampling(NUTS). For each parameter, Bulk_ESS
    and Tail_ESS are effective sample size measures, and Rhat is the potential
    scale reduction factor on split chains (at convergence, Rhat = 1).

``` r
summary(fit_lra_hs)
```

     Family: cox 
      Links: mu = log 
    Formula: Event_time | cens(1 - Event) ~ g_Lachnospirag_Streptococcus + g_Dialisterg_Roseburia + g_Prevotellag_Streptococcus + g_Roseburiag_Streptococcus + g_Dialistero_Clostridiales_g + g_Dialisterg_Lachnospira + g_Bacteroidesg_Streptococcus + o_Clostridiales_g_g_Streptococcus + f_Gemellaceae_g_g_Lachnospira + g_Bacteroidesf_Gemellaceae_g + g_Faecalibacteriumg_Streptococcus + g_Streptococcusg_Holdemania + g_Dialisterg_Prevotella + g_Bacteroidesg_Dialister + f_Gemellaceae_g_g_Prevotella + g_Faecalibacteriumg_Dialister + g_Bacteroidesg_Conchiformibius + g_Actinomycesg_Roseburia + g_Bacteroidesg_Granulicatella + f_Peptostreptococcaceae_g_g_Streptococcus 
       Data: df (Number of observations: 150) 
      Draws: 4 chains, each with iter = 4000; warmup = 2000; thin = 1;
             total post-warmup draws = 8000

    Regression Coefficients:
                                              Estimate Est.Error l-95% CI u-95% CI
    Intercept                                     0.39      0.43    -0.46     1.22
    g_Lachnospirag_Streptococcus                 -0.02      0.04    -0.12     0.03
    g_Dialisterg_Roseburia                        0.02      0.04    -0.02     0.11
    g_Prevotellag_Streptococcus                  -0.03      0.04    -0.14     0.02
    g_Roseburiag_Streptococcus                   -0.02      0.03    -0.11     0.03
    g_Dialistero_Clostridiales_g                  0.00      0.02    -0.05     0.07
    g_Dialisterg_Lachnospira                      0.02      0.03    -0.03     0.11
    g_Bacteroidesg_Streptococcus                 -0.01      0.02    -0.07     0.04
    o_Clostridiales_g_g_Streptococcus            -0.00      0.02    -0.05     0.05
    f_Gemellaceae_g_g_Lachnospira                 0.02      0.03    -0.02     0.11
    g_Bacteroidesf_Gemellaceae_g                 -0.01      0.03    -0.08     0.03
    g_Faecalibacteriumg_Streptococcus            -0.00      0.02    -0.04     0.04
    g_Streptococcusg_Holdemania                   0.00      0.02    -0.05     0.05
    g_Dialisterg_Prevotella                       0.04      0.05    -0.02     0.14
    g_Bacteroidesg_Dialister                     -0.01      0.03    -0.08     0.03
    f_Gemellaceae_g_g_Prevotella                  0.03      0.04    -0.02     0.14
    g_Faecalibacteriumg_Dialister                -0.00      0.02    -0.05     0.04
    g_Bacteroidesg_Conchiformibius               -0.03      0.04    -0.12     0.01
    g_Actinomycesg_Roseburia                      0.03      0.04    -0.02     0.12
    g_Bacteroidesg_Granulicatella                -0.00      0.02    -0.05     0.05
    f_Peptostreptococcaceae_g_g_Streptococcus    -0.03      0.04    -0.13     0.01
                                              Rhat Bulk_ESS Tail_ESS
    Intercept                                 1.00     7409     6727
    g_Lachnospirag_Streptococcus              1.00     5072     6876
    g_Dialisterg_Roseburia                    1.00     4094     6445
    g_Prevotellag_Streptococcus               1.00     3480     6263
    g_Roseburiag_Streptococcus                1.00     4655     5668
    g_Dialistero_Clostridiales_g              1.00     6964     6332
    g_Dialisterg_Lachnospira                  1.00     4793     5716
    g_Bacteroidesg_Streptococcus              1.00     6322     6911
    o_Clostridiales_g_g_Streptococcus         1.00     7281     6524
    f_Gemellaceae_g_g_Lachnospira             1.00     4767     6748
    g_Bacteroidesf_Gemellaceae_g              1.00     5732     6341
    g_Faecalibacteriumg_Streptococcus         1.00     8181     6996
    g_Streptococcusg_Holdemania               1.00     8274     6991
    g_Dialisterg_Prevotella                   1.00     2965     6268
    g_Bacteroidesg_Dialister                  1.00     5648     5587
    f_Gemellaceae_g_g_Prevotella              1.00     3561     6540
    g_Faecalibacteriumg_Dialister             1.00     7443     7012
    g_Bacteroidesg_Conchiformibius            1.00     4324     6517
    g_Actinomycesg_Roseburia                  1.00     4047     6729
    g_Bacteroidesg_Granulicatella             1.00     7617     6334
    f_Peptostreptococcaceae_g_g_Streptococcus 1.00     3984     6720

    Draws were sampled using sampling(NUTS). For each parameter, Bulk_ESS
    and Tail_ESS are effective sample size measures, and Rhat is the potential
    scale reduction factor on split chains (at convergence, Rhat = 1).

The models show good convergence diagnostics. All Rhat values are equal to 1.00, indicating that the chains have converged well. Bulk and Tail Effective Sample Sizes are also high, suggesting that the posterior distributions are well estimated and sampling was efficient.

The log-ratio coefficients are plotted below as hazard ratios.

``` r
# Plot
p1 <- plot_top_hr(fit_lra_normal,
                  title = "LRA – Normal prior",
                  top_k = 20)
p2 <- plot_top_hr(fit_lra_hs,
                  title = "LRA – Horseshoe prior",
                  top_k = 20)
p1 + p2
```

![](figures/unnamed-chunk-10-1.png)

In the LRA-transformed models, the normal prior produced some wide credible intervals, likely reflecting collinearity among the log-ratio features. Only one feature (f_Peptostreptococcaceae_g_g_Streptococcus) was significant. Under the horseshoe prior, all effects were more strongly shrunk toward zero, and no features remained significant, consistent with its stronger regularization.

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

![](figures/unnamed-chunk-11-1.png)

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

![](figures/unnamed-chunk-12-1.png)

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

![](figures/unnamed-chunk-13-1.png)

The point estimates from the IDP and Kaplan–Meier methods are very similar, but the uncertainty band from the IDP model lies slightly lower than uncertainty band from the Kaplan-Meier.

# Comparison of transformations Using Microbial Features in Survival Analysis

To evaluate how different data transformations affect survival modeling, we compared eight preprocessing approaches:
centered log-ratio (CLR), robust CLR (rCLR), pairwise log-ratio analysis (LRA), additive log-ratio (ALR), presence/absence (PA), total sum scaling (TSS), log-transformed TSS (logTSS), and arcsin square root (ASIN).

For each transformation, we selected the 20 microbial features most strongly associated with survival time based on univariate Cox regression. This consistent feature selection ensures a fair comparison across transformations.

Bayesian Cox proportional hazards models were fitted for each transformation using two different priors on the coefficients:

-   Normal(0,1) prior, representing weakly informative assumptions.

-   Horseshoe prior, promoting sparsity and reducing overfitting in high-dimensional settings.

Model performance for each prior–transformation combination was evaluated using leave-one-out cross-validation (LOO) to estimate out-of-sample predictive accuracy.

``` r
# Check Pareto k diagnostics
cat("Normal prior models – Pareto k > 0.7:\n")
```

    Normal prior models – Pareto k > 0.7:

``` r
print(sapply(loo_normal, function(x) sum(x$diagnostics$pareto_k > 0.7)))
```

       clr   rclr logtss    lra     pa    tss   asin    alr 
         0      0      0      0      0      0      0      0 

``` r
cat("\nHorseshoe prior models – Pareto k > 0.7:\n")
```


    Horseshoe prior models – Pareto k > 0.7:

``` r
print(sapply(loo_hs, function(x) sum(x$diagnostics$pareto_k > 0.7)))
```

       clr   rclr logtss    lra     pa    tss   asin    alr 
         0      0      0      0      0      0      0      0 

``` r
# Print LOO comparisons
cat("\nLOO comparison – Normal prior:\n")
```


    LOO comparison – Normal prior:

``` r
print(loo_comp_normal)
```

                  elpd_diff se_diff
    lra_normal      0.0       0.0  
    asin_normal    -1.8       5.6  
    alr_normal     -7.1       2.2  
    tss_normal     -8.9       6.5  
    pa_normal      -9.4       5.6  
    rclr_normal   -10.6       2.9  
    clr_normal    -10.9       3.0  
    logtss_normal -11.0       2.5  

``` r
cat("\nLOO comparison – Horseshoe prior:\n")
```


    LOO comparison – Horseshoe prior:

``` r
print(loo_comp_hs)
```

              elpd_diff se_diff
    lra_hs      0.0       0.0  
    asin_hs    -3.6       2.8  
    alr_hs     -5.1       1.2  
    rclr_hs    -5.8       1.6  
    clr_hs     -5.9       1.5  
    logtss_hs  -6.5       1.3  
    tss_hs     -8.1       3.9  
    pa_hs     -14.5       3.8  

``` r
# Create combined data frame
dff_all <- rbind(
  get_elpd_df(loo_normal, "Normal"),
  get_elpd_df(loo_hs, "Horseshoe")
)

# Keep model order consistent by average performance
model_order <- dff_all %>%
  group_by(model) %>%
  summarise(mean_elpd = mean(elpd, na.rm = TRUE)) %>%
  arrange(desc(mean_elpd)) %>%
  pull(model)

dff_all$model <- factor(dff_all$model, levels = model_order)

# Plot
ggplot(dff_all, aes(x = model, y = elpd, color = model, shape = prior)) +
  geom_point(size = 4, position = position_dodge(width = 0.5)) +
  geom_errorbar(
    aes(ymin = elpd - se, ymax = elpd + se),
    width = 0.15,
    size = 1,
    position = position_dodge(width = 0.5)
  ) +
  geom_hline(yintercept = max(dff_all$elpd), linetype = "dashed", color = "gray40") +
  scale_color_brewer(palette = "Set2") +
  labs(
    title = "Model Comparison via LOO-ELPD (Normal vs Horseshoe)",
    y = "Expected Log Predictive Density (ELPD)",
    x = "Transformation",
    shape = "Prior"
  ) +
  theme_bw(base_size = 14) +
  theme(
    plot.title = element_text(hjust = 0.5, face = "bold")
  )
```

![](figures/unnamed-chunk-16-1.png)

After applying moment matching, all Pareto $k$ values were below the recommended threshold of $0.7$, indicating that the LOO estimates are reliable for model comparison.

Across both prior specifications, the LRA (log-ratio analysis) model achieved the highest expected log predictive density (ELPD), showing the best out-of-sample predictive performance. Under the normal prior, performance differences between transformations were more pronounced, with CLR, rCLR, logTSS. With the horseshoe prior, differences between models narrowed, but LRA remained the top performer.The detection threshold for PA was relative (0.1% TSS), and it was the only method that performed better with the normal prior than with the horseshoe prior.

This pattern suggests that the stronger regularization of the horseshoe prior reduces performance gaps between transformations. Overall, these results reinforce that transformations, particularly LRA, yield more predictive survival models in this context.

# Model Predictions vs. Known Events

We plot the predicted risk scores from all transformations under the normal prior to examine how well each model distinguishes between events and non-events, and whether higher predicted risk corresponds to shorter observed survival times.

``` r
# Normal prior
df_rclr_n   <- get_pred_df(fit_rclr_normal,   method_name = "rCLR")
df_clr_n    <- get_pred_df(fit_clr_normal,    method_name = "CLR")
df_logtss_n <- get_pred_df(fit_logtss_normal, method_name = "logTSS")
df_lra_n    <- get_pred_df(fit_lra_normal,    method_name = "LRA")
df_pa_n     <- get_pred_df(fit_pa_normal,     method_name = "PA")
df_tss_n    <- get_pred_df(fit_tss_normal,    method_name = "TSS")
df_asin_n   <- get_pred_df(fit_asin_normal,   method_name = "Arcsin sqrt(TSS)")
df_alr_n    <- get_pred_df(fit_alr_normal,    method_name = "ALR")

df_all_n <- dplyr::bind_rows(
  df_clr_n, df_rclr_n, df_logtss_n, df_lra_n, df_pa_n, df_tss_n, df_asin_n, df_alr_n
)

ggplot(df_all_n, aes(x = risk_score, y = Event_time, color = factor(Event))) +
  geom_point(alpha = 0.7, size = 2) +
  geom_smooth(method = "loess", se = FALSE, color = "black", span = 1) +
  scale_color_manual(values = c("0" = "#00BFC4", "1" = "#F8766D")) +
  facet_wrap(~Method, scales = "free_x") +
  labs(
    title = "Predicted Risk Score vs. Survival Time (All Transformations) – Normal prior",
    x = "Predicted Risk Score (Posterior Median)",
    y = "Observed Survival Time",
    color = "Event"
  ) +
  theme_minimal(base_size = 14)
```

![](figures/unnamed-chunk-17-1.png)

The plot shows a consistent negative association between predicted risk score and survival time across most transformations. Higher risk scores are linked to shorter survival, with event cases (red) clustering toward high-risk, high-survival regions and censored cases (blue) concentrated at lower risk scores and longer survival times. The fitted trend lines reinforce this relationship, suggesting that the models capture meaningful patterns in outcome timing.

Next, we repeat the visualization for models fitted with the horseshoe prior.

``` r
# Horseshoe prior
df_rclr_h   <- get_pred_df(fit_rclr_hs,   method_name = "rCLR")
df_clr_h    <- get_pred_df(fit_clr_hs,    method_name = "CLR")
df_logtss_h <- get_pred_df(fit_logtss_hs, method_name = "logTSS")
df_lra_h    <- get_pred_df(fit_lra_hs,    method_name = "LRA")
df_pa_h     <- get_pred_df(fit_pa_hs,     method_name = "PA")
df_tss_h    <- get_pred_df(fit_tss_hs,    method_name = "TSS")
df_asin_h   <- get_pred_df(fit_asin_hs,   method_name = "Arcsin sqrt(TSS)")
df_alr_h    <- get_pred_df(fit_alr_hs,    method_name = "ALR")

df_all_h <- dplyr::bind_rows(
  df_clr_h, df_rclr_h, df_logtss_h, df_lra_h, df_pa_h, df_tss_h, df_asin_h, df_alr_h
)

ggplot(df_all_h, aes(x = risk_score, y = Event_time, color = factor(Event))) +
  geom_point(alpha = 0.7, size = 2) +
  geom_smooth(method = "loess", se = FALSE, color = "black", span = 1) +
  scale_color_manual(values = c("0" = "#00BFC4", "1" = "#F8766D")) +
  facet_wrap(~Method, scales = "free_x") +
  labs(
    title = "Predicted Risk Score vs. Survival Time (All Transformations) – Horseshoe prior",
    x = "Predicted Risk Score (Posterior Median)",
    y = "Observed Survival Time",
    color = "Event"
  ) +
  theme_minimal(base_size = 14)
```

![](figures/unnamed-chunk-18-1.png)

The overall patterns are visually very similar to those observed with the normal prior, indicating that the choice of prior has little effect on the qualitative relationship between predicted risk and survival time in this dataset.

# Posterior-based AUROC and C-index Comparison of Survival Models

We assess model discrimination with two complementary metrics computed from the posterior. First, we calculate AUROC for every posterior draw of the linear predictor against the observed events. This gives a distribution of AUROC values that naturally reflects parameter uncertainty. Second, we compute Harrell’s C-index for each draw using the observed survival times and censoring, interpreting higher scores as higher risk (earlier events). We report both metrics as posterior medians with 95% credible intervals.

``` r
# Read RDS
auc_all <- readRDS(file.path("model_result", "posterior_auc_all.rds"))

# Summarise: median + 95% credible interval for AUROC
auc_summary <- auc_all %>%
  group_by(model, prior) %>%
  summarise(
    median_auc = median(AUROC),
    lower = quantile(AUROC, 0.025),
    upper = quantile(AUROC, 0.975),
    .groups = "drop"
  )

# Order models by mean AUROC across priors
model_order <- auc_summary %>%
  group_by(model) %>%
  summarise(mean_auc = mean(median_auc), .groups = "drop") %>%
  arrange(desc(mean_auc)) %>%
  pull(model)

auc_summary$model <- factor(auc_summary$model, levels = model_order)

# Plot
ggplot(auc_summary, aes(x = model, y = median_auc, color = model, shape = prior)) +
  geom_point(size = 4, position = position_dodge(width = 0.5)) +
  geom_errorbar(
    aes(ymin = lower, ymax = upper),
    width = 0.15,
    size = 1,
    position = position_dodge(width = 0.5)
  ) +
  geom_hline(yintercept = max(auc_summary$median_auc), linetype = "dashed", color = "gray40") +
  scale_color_brewer(palette = "Set2") +
  labs(
    title = "Model Comparison via Posterior AUROC (Normal vs Horseshoe)",
    y = "AUROC",
    x = "Transformation",
    shape = "Prior"
  ) +
  theme_bw(base_size = 14) +
  theme(
    plot.title = element_text(hjust = 0.5, face = "bold")
  )
```

![](figures/unnamed-chunk-19-1.png)

``` r
# Read RDS
c_all <- readRDS(file.path("model_result", "cindex_posterior_all.rds"))

# Summarise: median + 95% credible interval
c_summary <- c_all %>%
  filter(is.finite(Cindex)) %>%
  group_by(model, prior) %>%
  summarise(
    median_c = median(Cindex),
    lower    = quantile(Cindex, 0.025),
    upper    = quantile(Cindex, 0.975),
    .groups  = "drop"
  )

# Order models by average median C across priors
model_order <- c_summary %>%
  group_by(model) %>%
  summarise(mean_med = mean(median_c), .groups = "drop") %>%
  arrange(desc(mean_med)) %>%
  pull(model)

c_summary$model <- factor(c_summary$model, levels = model_order)

# Plot
ggplot(c_summary, aes(x = model, y = median_c, color = model, shape = prior)) +
  geom_point(size = 4, position = position_dodge(width = 0.5)) +
  geom_errorbar(
    aes(ymin = lower, ymax = upper),
    width = 0.15, size = 1,
    position = position_dodge(width = 0.5)
  ) +
  geom_hline(yintercept = max(c_summary$median_c), linetype = "dashed", color = "gray40") +
  scale_color_brewer(palette = "Set2") +
  labs(
    title = "Model Comparison via Posterior Harrell's C (Normal vs Horseshoe)",
    y = "Harrell's C",
    x = "Transformation",
    shape = "Prior"
  ) +
  theme_bw(base_size = 14) +
  theme(plot.title = element_text(hjust = 0.5, face = "bold"))
```

![](figures/unnamed-chunk-20-1.png)

In both the AUROC and Harrell’s C panels the story is consistent. Log-ratio-based models perform best, LRA leads, with ALR next, while CLR/rCLR are competitive but a step down. Abundance transforms (logTSS, asin) sit in the middle. PA and especially TSS are clearly weakest and show the widest uncertainty. Across priors, Normal is generally on par with or slightly ahead of Horseshoe, suggesting strong shrinkage trims some discriminative signal, the clearest Normal advantage appears for PA. The time-aware C-index corroborates the AUROC ranking.

# Kaplan–Meier curves by median risk (Normal-prior models)

For each model, we compute a risk score as the posterior median of its linear predictor and split individuals at the model-specific median into Low vs High risk. The panels show Kaplan–Meier curves with 95% CIs for each transformation.

``` r
# Collect risk scores from Normal-prior models
df_all_n <- dplyr::bind_rows(
  get_pred_df(fit_rclr_normal,   "rCLR"),
  get_pred_df(fit_clr_normal,    "CLR"),
  get_pred_df(fit_logtss_normal, "logTSS"),
  get_pred_df(fit_lra_normal,    "LRA"),
  get_pred_df(fit_pa_normal,     "PA"),
  get_pred_df(fit_tss_normal,    "TSS"),
  get_pred_df(fit_asin_normal,   "Arcsin sqrt(TSS)"),
  get_pred_df(fit_alr_normal,    "ALR")
)

# Median split within each model: Low vs High
df_all_n <- df_all_n %>%
  dplyr::group_by(Method) %>%
  dplyr::mutate(
    RiskGroup = if_else(risk_score >= median(risk_score, na.rm = TRUE), "High", "Low"),
    RiskGroup = factor(RiskGroup, levels = c("Low","High"))
  ) %>%
  dplyr::ungroup()

# Kaplan–Meier data per model and risk group
km_df <- df_all_n %>%
  dplyr::group_split(Method) %>%
  purrr::map_dfr(function(d){
    meth <- as.character(d$Method[1])
    sf <- survival::survfit(Surv(Event_time, Event) ~ RiskGroup, data = d,
                            conf.int = 0.95, conf.type = "log-log")
    survminer::surv_summary(sf, data = d) %>%
      dplyr::mutate(Method = meth, RiskGroup = sub("^RiskGroup=", "", strata)) %>%
      dplyr::select(Method, time, surv, lower, upper, RiskGroup)
  })

# Plot
ggplot(km_df, aes(x = time, y = surv, color = RiskGroup, fill = RiskGroup)) +
  geom_step(linewidth = 1) +
  geom_ribbon(aes(ymin = lower, ymax = upper), alpha = 0.15, color = NA) +
  facet_wrap(~ Method, ncol = 4) +
  labs(
    title = "Kaplan–Meier by median risk (High vs Low) — Normal prior models",
    x = "Time",
    y = "Survival probability",
    color = "Risk group",
    fill  = "Risk group"
  ) +
  theme_bw(base_size = 14) +
  theme(plot.title = element_text(hjust = 0.5, face = "bold"),
        panel.spacing.x = unit(0.8, "lines"))
```

![](figures/unnamed-chunk-21-1.png)

The panels show clear risk separation across all Normal-prior models. The High group (red) consistently has lower survival than the Low group (blue), with gaps opening early and often widening over time. Log-ratio transforms (ALR, LRA, CLR, rCLR) perform best, logTSS and asin are close behind, and TSS is weakest. The split is model specific.

We will also repeat the analysis with a four-group split using the 25th, 50th, and 75th percentiles (Q1–Q4)

``` r
# Quartile split (25/50/75%) within each model: Q1 (Low), Q4 (High)
df_q <- df_all_n %>%
  dplyr::group_by(Method) %>%
  dplyr::mutate(
    RiskQ = dplyr::ntile(risk_score, 4),
    RiskQ = factor(RiskQ, levels = 1:4,
                   labels = c("Q1 (Low)", "Q2", "Q3", "Q4 (High)"))
  ) %>%
  dplyr::ungroup()

# KM curves by quartile for each model
km_df_q <- df_q %>%
  dplyr::group_split(Method) %>%
  purrr::map_dfr(function(d){
    meth <- as.character(d$Method[1])
    sf <- survival::survfit(Surv(Event_time, Event) ~ RiskQ, data = d,
                            conf.int = 0.95, conf.type = "log-log")
    survminer::surv_summary(sf, data = d) %>%
      dplyr::mutate(Method = meth, RiskQ = sub("^RiskQ=", "", strata)) %>%
      dplyr::select(Method, time, surv, lower, upper, RiskQ)
  })

# Plot
ggplot(km_df_q, aes(x = time, y = surv, color = RiskQ, fill = RiskQ)) +
  geom_step(linewidth = 1) +
  geom_ribbon(aes(ymin = lower, ymax = upper), alpha = 0.15, color = NA) +
  facet_wrap(~ Method, ncol = 4) +
  scale_color_brewer(palette = "Set1") +
  scale_fill_brewer(palette  = "Set1") +
  labs(
    title = "Kaplan–Meier by risk quartiles (Q1–Q4) — Normal prior models",
    x = "Time",
    y = "Survival probability",
    color = "Risk quartile",
    fill  = "Risk quartile"
  ) +
  theme_bw(base_size = 14) +
  theme(plot.title = element_text(hjust = 0.5, face = "bold"),
        panel.spacing.x = unit(0.8, "lines"))
```

![](figures/unnamed-chunk-22-1.png)

The quartile panels show a clear risk gradient. Q4 (High) has the lowest survival and Q1 (Low) the highest, with gaps widening over time. Log ratio transforms separate best, and TSS is weakest.

# Comparing Más-o-menos scoring with the LRA model

To explore alternative strategies for risk prediction, we implemented the [Más-o-menos method](https://academic.oup.com/bioinformatics/article/30/21/3062/2422201) (Armañanzas et al., 2014), a simple sign-based approach that averages the direction of effect from univariate Cox models. Each feature contributes either positively or negatively based on the sign of its estimated coefficient, and the overall risk score is calculated as the mean of these directional contributions.

We compared this score against the log-ratio-based (LRA) Cox models (normal and horseshoe), which had previously shown the best predictive performance among the tested transformations. Both methods were evaluated on the same set of top-ranked microbial features, and their predicted risk scores were compared in relation to observed event outcomes.

``` r
# LRA, normal prior
pred_normal <- posterior_linpred(fit_lra_normal, transform = FALSE)
df_lra_normal <- df_lra
df_lra_normal$risk_score <- apply(pred_normal, 2, median)
df_lra_normal$method <- "LRA – Normal"

# LRA, horseshoe prior
pred_hs <- posterior_linpred(fit_lra_hs, transform = FALSE)
df_lra_hs <- df_lra
df_lra_hs$risk_score <- apply(pred_hs, 2, median)
df_lra_hs$method <- "LRA – Horseshoe"

# Más-o-menos
df_masomenos <- df_lra
df_masomenos$risk_score <- calculate_masomenos(df_masomenos)
df_masomenos$method <- "Más-o-menos"

df_all <- rbind(
  df_lra_normal[, c("risk_score", "Event", "method")],
  df_lra_hs[, c("risk_score", "Event", "method")],
  df_masomenos[, c("risk_score", "Event", "method")]
)

ggplot(df_all, aes(x = factor(Event), y = risk_score, fill = method)) +
  geom_boxplot(position = position_dodge(width = 0.75)) +
  scale_fill_manual(values = c(
    "LRA – Normal"    = "#1f78b4",
    "LRA – Horseshoe" = "#33a02c",
    "Más-o-menos"     = "#e31a1c"
  )) +
  labs(
    title = "Predicted Risk Scores by Event Status",
    x = "Event (0 = censored, 1 = occurred)",
    y = "Risk Score",
    fill = "Method"
  ) +
  theme_minimal(base_size = 14)
```

![](figures/unnamed-chunk-23-1.png)

The figure indicates that all three methods produce higher predicted risk scores for individuals who experienced the event compared to those who were censored, showing that they capture the expected risk–outcome relationship reasonably well. The separation between groups is visible for each method, but given the degree of overlap in score distributions, no firm conclusions can be drawn about differences in performance between the methods based on this plot alone.

# Comparing different survival models across transformations (Updated)

We compare seven survival models: Cox proportional hazards model, logistic regression, Random Survival Forests, XGBoost, DeepSurv, CatBoost, and TabPFN. To ensure a fair comparison, we evaluate all available taxonomy/feature-coding schemes and apply the same choice across methods. Transformations and models are evaluated with 5-fold cross-validation using Harrell’s C-index as the primary metric. For computational reasons, the Cox model is implemented in a frequentist rather than Bayesian form.

``` r
# Combine all models
combined_metrics <- bind_rows(cox, rsf, xgb, deep, logit, cb, tab)
```

``` r
# Method label order
method_levels <- c(
  "CoxPH_5CV",
  "RSF_5CV",
  "XGB_Cox_5CV",
  "DeepSurv_5CV",
  "Logit_Binomial_5CV",
  "CatBoost_Binomial_5CV",
  "TabPFN_Binomial_5CV"
)

# Use C-index rows for ordering
cindex_metrics <- combined_metrics %>%
  dplyr::filter(metric == "C") %>%
  dplyr::mutate(method = as.character(method))
present <- intersect(method_levels, unique(cindex_metrics$method))
others  <- setdiff(unique(cindex_metrics$method), present)
method_levels_use <- c(present, others)

cindex_metrics <- cindex_metrics %>%
  dplyr::mutate(method = factor(method, levels = method_levels_use))

# Order transformations by mean C across methods
model_order <- cindex_metrics %>%
  dplyr::group_by(model) %>%
  dplyr::summarise(mean_cindex = mean(estimate, na.rm = TRUE), .groups = "drop") %>%
  dplyr::arrange(dplyr::desc(mean_cindex)) %>%
  dplyr::pull(model)

cindex_metrics$model <- factor(cindex_metrics$model, levels = model_order)

# Plot 
ggplot(cindex_metrics, aes(x = model, y = estimate, color = method)) +
  geom_point(size = 4, stroke = 1.1,
             position = position_dodge(width = 0.7), na.rm = TRUE) +
  geom_errorbar(aes(ymin = lower, ymax = upper),
                width = 0.15, size = 1,
                position = position_dodge(width = 0.7), na.rm = TRUE) +
  geom_hline(yintercept = 0.5, linetype = "dashed") +
  labs(
    title = "C-index: CoxPH, RSF, XGB-Cox, DeepSurv, Logit, CatBoost, TabPFN",
    x = "Transformation",
    y = "Harrell's C",
    color = "Method"
  ) +
  coord_cartesian(ylim = c(0.45, 1.00)) +
  theme_bw(base_size = 14) +
  theme(
    plot.title = element_text(hjust = 0.5, face = "bold"),
    axis.text.x = element_text(angle = 20, hjust = 1)
  )
```

![](figures/unnamed-chunk-26-1.png)
