# Packages
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
library(Matrix)
library(patchwork)
library(pROC)
library(RColorBrewer) 
library(cmdstanr)
library(ranger)
library(timeROC)
library(xgboost)
library(matrixStats)
library(catboost)
library(reticulate)

# Load funcs
source("Funcs_cleaned.R")

# Load the TreeSummarizedExperiment object
tse <- readRDS("survival_tse.rds")

# All features data
# Set once
# - Change ASSAY_IN if your counts are not under 'counts'.
# - Change EVENT_COL/TIME_COL if your metadata uses different names.
ASSAY_IN    <- "counts"         # input assay name
EVENT_COL   <- "Event"          # event (0/1)
TIME_COL    <- "Event_time"     # follow-up time
PSEUDOCOUNT <- 1e-6             # for log-based transforms
ALR_REF     <- "g_Turicibacter" # required for ALR


df_clr      <- extract_all_features_by_transformation(
  tse, method = "clr",
  assay_in = ASSAY_IN, pseudocount = PSEUDOCOUNT,
  event_col = EVENT_COL, time_col = TIME_COL
)
df_rclr     <- extract_all_features_by_transformation(
  tse, method = "rclr",
  assay_in = ASSAY_IN, pseudocount = PSEUDOCOUNT,
  event_col = EVENT_COL, time_col = TIME_COL
)
df_logabund <- extract_all_features_by_transformation(
  tse, method = "log_abund",
  assay_in = ASSAY_IN, pseudocount = PSEUDOCOUNT,
  event_col = EVENT_COL, time_col = TIME_COL
)
df_lra      <- extract_all_features_by_transformation(
  tse, method = "lra",
  assay_in = ASSAY_IN, pseudocount = PSEUDOCOUNT,
  event_col = EVENT_COL, time_col = TIME_COL
)
df_pa       <- extract_all_features_by_transformation(
  tse, method = "pa",
  assay_in = ASSAY_IN,
  event_col = EVENT_COL, time_col = TIME_COL
)
df_tss      <- extract_all_features_by_transformation(
  tse, method = "tss",
  assay_in = ASSAY_IN,
  event_col = EVENT_COL, time_col = TIME_COL
)
df_logtss   <- extract_all_features_by_transformation(
  tse, method = "logtss",
  assay_in = ASSAY_IN, pseudocount = PSEUDOCOUNT,
  event_col = EVENT_COL, time_col = TIME_COL
)
df_asin     <- extract_all_features_by_transformation(
  tse, method = "asin",
  assay_in = ASSAY_IN,
  event_col = EVENT_COL, time_col = TIME_COL
)
# For ALR set `alr_ref` to the reference taxon ID used in your features.
df_alr      <- extract_all_features_by_transformation(
  tse, method = "alr",
  assay_in = ASSAY_IN, pseudocount = PSEUDOCOUNT, alr_ref = ALR_REF,
  event_col = EVENT_COL, time_col = TIME_COL
)


# Transforms
transforms <- list(
  CLR    = df_clr,
  rCLR   = df_rclr,
  logTSS = df_logtss,
  LRA    = df_lra,
  PA     = df_pa,
  TSS    = df_tss,
  Asin   = df_asin,
  ALR    = df_alr
)

# CV setup
CV_SEED  <- 1
K_FOLDS  <- 5
set.seed(CV_SEED)
global_folds <- make_stratified_folds(df_clr[[EVENT_COL]], K = K_FOLDS, seed = CV_SEED)

# CoxPH
cox_res <- purrr::imap_dfr(
  transforms,
  ~ coxph_cindex_cv5(.x, .y,
                     ties = "efron",
                     seed = CV_SEED, K = K_FOLDS, fold_id = global_folds,
                     event_col = EVENT_COL, time_col = TIME_COL)
)

# Random Survival Forest
rsf_res <- purrr::imap_dfr(
  transforms,
  ~ rsf_cindex_cv5(.x, .y,
                   seed = CV_SEED, K = K_FOLDS, fold_id = global_folds,
                   event_col = EVENT_COL, time_col = TIME_COL)
)

# Logistic (binomial)
logit_res <- purrr::imap_dfr(
  transforms,
  ~ logit_cindex_cv5(.x, .y, maxit = 200,
                     seed = CV_SEED, K = K_FOLDS, fold_id = global_folds,
                     event_col = EVENT_COL, time_col = TIME_COL)
)

# TabPFN (binomial)
# Lock Python to 'r-tabpfn' (do this before any keras/tensorflow usage)
if (!"reticulate" %in% .packages(TRUE)) install.packages("reticulate")
library(reticulate)
Sys.setenv(RETICULATE_PYTHON = reticulate::conda_python("r-tabpfn"))
reticulate::py_config()

# Ensure required Python pkgs exist in r-tabpfn
if (!reticulate::py_module_available("tabpfn")) {
  reticulate::conda_install("r-tabpfn", c("pip"))
  reticulate::conda_install("r-tabpfn", pip = TRUE,
                            packages = c("torch", "tabpfn @ git+https://github.com/PriorLabs/TabPFN.git"))
}

# Initialize TabPFN helpers (your function defined earlier)
tabpfn_init("r-tabpfn")

tabpfn_res <- purrr::imap_dfr(
  transforms,
  ~ tabpfn_bin_cindex_cv5(
    df = .x, method_name = .y,
    seed = CV_SEED, K = K_FOLDS, fold_id = global_folds,
    device = "cpu", ensemble = 32,
    event_col = EVENT_COL, time_col = TIME_COL
  )
)

# DeepSurv
library(keras)
library(tensorflow)
deepsurv_res <- purrr::imap_dfr(
  transforms,
  ~ deepsurv_cindex_cv5(
    .x, .y,
    hidden = c(32, 16), dropout = 0.1, l2 = 1e-4,
    lr = 1e-3, epochs = 120, patience = 12,
    verbose = 0, run_eagerly = TRUE,
    seed = CV_SEED, K = K_FOLDS, fold_id = global_folds,
    event_col = EVENT_COL, time_col = TIME_COL
  )
)

# CatBoost (binomial)
cb_res <- purrr::imap_dfr(
  transforms,
  ~ catboost_bin_cindex_cv5(.x, .y,
                            seed = CV_SEED, K = K_FOLDS, fold_id = global_folds,
                            event_col = EVENT_COL, time_col = TIME_COL)
)

# XGBoost Cox
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
    fold_id = global_folds, use_gpu = FALSE,
    event_col = EVENT_COL, time_col = TIME_COL
  )
)


# Combine
combined_metrics <- bind_rows(
  cox_res, rsf_res, xgb_res, deepsurv_res, logit_res, cb_res, tabpfn_res
)

# Fix order by mean
cindex <- combined_metrics %>%
  filter(metric == "C") %>%
  mutate(method = factor(method,
                         levels = c("CoxPH","RSF","XGB_Cox","DeepSurv","Logit","CatBoost","TabPFN"))
  )
model_order <- cindex %>%
  group_by(model) %>%
  summarise(meanC = mean(estimate, na.rm = TRUE), .groups = "drop") %>%
  arrange(desc(meanC)) %>%
  pull(model)

cindex <- cindex %>% mutate(model = factor(model, levels = model_order))

# Plot
ggplot(cindex, aes(x = model, y = estimate, color = method)) +
  geom_point(position = position_dodge(0.6), size = 3) +
  geom_errorbar(aes(ymin = lower, ymax = upper),
                position = position_dodge(0.6), width = 0.15) +
  geom_hline(yintercept = 0.5, linetype = "dashed") +
  labs(title = "C-index across transformations (ordered by mean C)",
       x = "Transformation", y = "Harrell's C", color = "Method") +
  theme_bw() +
  theme(axis.text.x = element_text(angle = 20, hjust = 1))
