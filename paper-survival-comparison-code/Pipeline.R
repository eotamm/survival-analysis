# Packages
library(tidyverse)
library(brms)
library(tidybayes)
library(ggplot2)
library(TreeSummarizedExperiment)
library(SingleCellExperiment)
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
library(fastshap)
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
ASSAY_IN    <- "counts"                   # input assay name
EVENT_COL   <- "Event"                    # event (0/1)
TIME_COL    <- "Event_time"               # follow-up time
PSEUDOCOUNT <- 1e-6                       # for log-based transforms

#  Core settings
PREV_MIN      <- 0.30   # min fraction of samples where taxon is "present"
# presence if relative abundance >= 0.1% for this data as there is no zeros, but if it has this can be 0
REL_PRES_MIN  <- 1e-3

# Relative abundance
tse <- transformAssay(tse, method = "relabundance",
                           assay.type = ASSAY_IN, name = "tss")

rel <- SummarizedExperiment::assay(tse, "tss")

# Prevalence on relative scale
# If REL_PRES_MIN == 0, uses strict '>' so zeros are not counted present
prev <- if (REL_PRES_MIN == 0) {
  rowMeans(rel > 0)
} else {
  rowMeans(rel >= REL_PRES_MIN)
}

# Core taxa by prevalence
core_taxa <- names(prev)[prev >= PREV_MIN]
if (length(core_taxa) < 2L) {
  stop("Core is too small. Lower PREV_MIN / REL_PRES_MIN.")
}

# Subset TSE to core only
tse_core <- tse[core_taxa, , drop = FALSE]

# Core summary
message(sprintf("Core taxa: %d | PREV_MIN=%.2f | REL_PRES_MIN=%g",
                length(core_taxa), PREV_MIN, REL_PRES_MIN))

# Build data frames for transformations
df_clr <- extract_all_features_by_transformation(
  tse_core, method = "clr",
  assay_in = ASSAY_IN, pseudocount = PSEUDOCOUNT,
  event_col = EVENT_COL, time_col = TIME_COL
)

df_rclr <- extract_all_features_by_transformation(
  tse_core, method = "rclr",
  assay_in = ASSAY_IN, pseudocount = PSEUDOCOUNT,
  event_col = EVENT_COL, time_col = TIME_COL
)

df_logabund <- extract_all_features_by_transformation(
  tse_core, method = "log_abund",
  assay_in = ASSAY_IN, pseudocount = PSEUDOCOUNT,
  event_col = EVENT_COL, time_col = TIME_COL
)

df_lra <- extract_all_features_by_transformation(
  tse_core, method = "lra",
  assay_in = ASSAY_IN, pseudocount = PSEUDOCOUNT,
  event_col = EVENT_COL, time_col = TIME_COL
)

df_pa <- extract_all_features_by_transformation(
  tse_core, method = "pa",
  assay_in = ASSAY_IN,
  event_col = EVENT_COL, time_col = TIME_COL
)

df_tss <- extract_all_features_by_transformation(
  tse_core, method = "tss",
  assay_in = ASSAY_IN,
  event_col = EVENT_COL, time_col = TIME_COL
)

df_logtss <- extract_all_features_by_transformation(
  tse_core, method = "logtss",
  assay_in = ASSAY_IN, pseudocount = PSEUDOCOUNT,
  event_col = EVENT_COL, time_col = TIME_COL
)

df_asin <- extract_all_features_by_transformation(
  tse_core, method = "asin",
  assay_in = ASSAY_IN,
  event_col = EVENT_COL, time_col = TIME_COL
)

# No reference, so the transform uses the first feature as the reference
df_alr <- extract_all_features_by_transformation(
  tse_core, method = "alr",
  assay_in = ASSAY_IN, pseudocount = PSEUDOCOUNT,
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


out_dir <- file.path("model_result", "shap_outputs")
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

# Run models
# Coxnet
coxnet_list <- purrr::imap(
  transforms,
  ~ coxnet_cindex_cv5(
    df = .x, method_name = .y,
    alpha = 0,
    seed = CV_SEED, K = K_FOLDS, fold_id = global_folds,
    event_col = EVENT_COL, time_col = TIME_COL,
    shap = TRUE, shap_nsim = 16, shap_bg_max = 32, shap_verbose = TRUE
  )
)

coxnet_res <- dplyr::bind_rows(coxnet_list)
coxnet_foldC <- coxnet_list %>% purrr::map(~ attr(.x, "foldC")) %>% purrr::compact() %>% dplyr::bind_rows()
coxnet_shap_long <- coxnet_list %>% purrr::map(~ attr(.x, "shap_long")) %>% purrr::compact() %>% dplyr::bind_rows()
coxnet_shap_agg  <- coxnet_list %>% purrr::map(~ attr(.x, "shap_agg"))  %>% purrr::compact() %>% dplyr::bind_rows()

saveRDS(coxnet_res,       file.path(out_dir, "coxnet_metrics.rds"))
saveRDS(coxnet_foldC,     file.path(out_dir, "coxnet_foldC.rds"))
saveRDS(coxnet_shap_long, file.path(out_dir, "coxnet_shap_long.rds"))
saveRDS(coxnet_shap_agg,  file.path(out_dir, "coxnet_shap_agg.rds"))


# RSF
rsf_list <- purrr::imap(
  transforms,
  ~ rsf_cindex_cv5(
    df = .x, method_name = .y,
    seed = CV_SEED, K = K_FOLDS, fold_id = global_folds,
    event_col = EVENT_COL, time_col = TIME_COL,
    shap = TRUE, shap_nsim = 16, shap_bg_max = 32
  )
)
rsf_res <- dplyr::bind_rows(rsf_list)
rsf_foldC <- rsf_list %>% purrr::map(~ attr(.x, "foldC")) %>% purrr::compact() %>% dplyr::bind_rows()
rsf_shap_long <- rsf_list %>% purrr::map(~ attr(.x, "shap_long")) %>% purrr::compact() %>% dplyr::bind_rows()
rsf_shap_agg  <- rsf_list %>% purrr::map(~ attr(.x, "shap_agg"))  %>% purrr::compact() %>% dplyr::bind_rows()

saveRDS(rsf_res,       file.path(out_dir, "rsf_metrics.rds"))
saveRDS(rsf_foldC,     file.path(out_dir, "rsf_foldC.rds"))
saveRDS(rsf_shap_long, file.path(out_dir, "rsf_shap_long.rds"))
saveRDS(rsf_shap_agg,  file.path(out_dir, "rsf_shap_agg.rds"))


# Logistic (binomial)
logit_list <- purrr::imap(
  .x = transforms,
  .f = ~ logit_cindex_cv5(
    df = .x, method_name = .y, maxit = 200,
    seed = CV_SEED, K = K_FOLDS, fold_id = global_folds,
    event_col = EVENT_COL, time_col = TIME_COL,
    shap = TRUE, shap_nsim = 16, shap_bg_max = 32, shap_verbose = TRUE
  )
)
logit_foldC <- logit_list %>% purrr::map(~ attr(.x, "foldC")) %>% purrr::compact() %>% dplyr::bind_rows()
logit_res <- dplyr::bind_rows(logit_list)
logit_shap_long <- logit_list %>% purrr::map(~ attr(.x, "shap_long")) %>% purrr::compact() %>% dplyr::bind_rows()
logit_shap_agg  <- logit_list %>% purrr::map(~ attr(.x, "shap_agg"))  %>% purrr::compact() %>% dplyr::bind_rows()

saveRDS(logit_res,       file.path(out_dir, "logit_metrics.rds"))
saveRDS(logit_foldC,     file.path(out_dir, "logit_foldC.rds"))
saveRDS(logit_shap_long, file.path(out_dir, "logit_shap_long.rds"))
saveRDS(logit_shap_agg,  file.path(out_dir, "logit_shap_agg.rds"))

# CatBoost (binomial)
cb_list <- purrr::imap(
  transforms,
  ~ catboost_bin_cindex_cv5(
    df = .x, method_name = .y,
    seed = CV_SEED, K = K_FOLDS, fold_id = global_folds,
    event_col = EVENT_COL, time_col = TIME_COL,
    shap = TRUE, shap_nsim = 16, shap_bg_max = 32
  )
)
cb_res <- dplyr::bind_rows(cb_list)
cb_foldC <- cb_list %>% purrr::map(~ attr(.x, "foldC")) %>% purrr::compact() %>% dplyr::bind_rows()
cb_shap_long <- cb_list %>% purrr::map(~ attr(.x, "shap_long")) %>% purrr::compact() %>% dplyr::bind_rows()
cb_shap_agg  <- cb_list %>% purrr::map(~ attr(.x, "shap_agg"))  %>% purrr::compact() %>% dplyr::bind_rows()

saveRDS(cb_res,       file.path(out_dir, "cb_metrics.rds"))
saveRDS(cb_foldC,     file.path(out_dir, "cb_foldC.rds"))
saveRDS(cb_shap_long, file.path(out_dir, "cb_shap_long.rds"))
saveRDS(cb_shap_agg,  file.path(out_dir, "cb_shap_agg.rds"))

# TabPFN (binomial)
# reticulate env + deps
library(reticulate)
Sys.setenv(RETICULATE_PYTHON = reticulate::conda_python("r-tabpfn"))
reticulate::py_config()

if (!reticulate::py_module_available("tabpfn")) {
  reticulate::conda_install("r-tabpfn", c("pip"))
  reticulate::conda_install("r-tabpfn", pip = TRUE,
                            packages = c("torch", "tabpfn @ git+https://github.com/PriorLabs/TabPFN.git"))
}
tabpfn_init("r-tabpfn")
tabpfn_list <- purrr::imap(
  transforms,
  ~ tabpfn_bin_cindex_cv5(
    df = .x, method_name = .y,
    seed = CV_SEED, K = K_FOLDS, fold_id = global_folds,
    device = "cpu", ensemble = 1,
    event_col = EVENT_COL, time_col = TIME_COL,
    shap = TRUE, shap_nsim = 2, shap_bg_max = 16
  )
)
tabpfn_res <- dplyr::bind_rows(tabpfn_list)
tabpfn_foldC <- tabpfn_list %>% purrr::map(~ attr(.x, "foldC")) %>% purrr::compact() %>% dplyr::bind_rows()
tabpfn_shap_long <- tabpfn_list %>% purrr::map(~ attr(.x, "shap_long")) %>% purrr::compact() %>% dplyr::bind_rows()
tabpfn_shap_agg  <- tabpfn_list %>% purrr::map(~ attr(.x, "shap_agg"))  %>% purrr::compact() %>% dplyr::bind_rows()

saveRDS(tabpfn_res,       file.path(out_dir, "tabpfn_metrics.rds"))
saveRDS(tabpfn_foldC,     file.path(out_dir, "tabpfn_foldC.rds"))
saveRDS(tabpfn_shap_long, file.path(out_dir, "tabpfn_shap_long.rds"))
saveRDS(tabpfn_shap_agg,  file.path(out_dir, "tabpfn_shap_agg.rds"))

# DeepSurv
library(keras); library(tensorflow)
deepsurv_list <- purrr::imap(
  transforms,
  ~ deepsurv_cindex_cv5(
    df = .x, method_name = .y,
    hidden = c(32, 16), dropout = 0.1, l2 = 1e-4,
    lr = 1e-3, epochs = 120, patience = 12,
    verbose = 0, run_eagerly = TRUE,
    seed = CV_SEED, K = K_FOLDS, fold_id = global_folds,
    event_col = EVENT_COL, time_col = TIME_COL,
    shap = TRUE, shap_nsim = 16, shap_bg_max = 32
  )
)
deepsurv_res <- dplyr::bind_rows(deepsurv_list)
deepsurv_foldC <- deepsurv_list %>% purrr::map(~ attr(.x, "foldC")) %>% purrr::compact() %>% dplyr::bind_rows()
deepsurv_shap_long <- deepsurv_list %>% purrr::map(~ attr(.x, "shap_long")) %>% purrr::compact() %>% dplyr::bind_rows()
deepsurv_shap_agg  <- deepsurv_list %>% purrr::map(~ attr(.x, "shap_agg"))  %>% purrr::compact() %>% dplyr::bind_rows()

saveRDS(deepsurv_res,       file.path(out_dir, "deepsurv_metrics.rds"))
saveRDS(deepsurv_foldC,     file.path(out_dir, "deepsurv_foldC.rds"))
saveRDS(deepsurv_shap_long, file.path(out_dir, "deepsurv_shap_long.rds"))
saveRDS(deepsurv_shap_agg,  file.path(out_dir, "deepsurv_shap_agg.rds"))

# XGBoost Cox
xgb_list <- purrr::imap(
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
    event_col = EVENT_COL, time_col = TIME_COL,
  shap = TRUE, shap_nsim = 16, shap_bg_max = 32
  )
)
xgb_res <- dplyr::bind_rows(xgb_list)
xgb_foldC <- xgb_list %>% purrr::map(~ attr(.x, "foldC")) %>% purrr::compact() %>% dplyr::bind_rows()
xgb_shap_long <- xgb_list %>% purrr::map(~ attr(.x, "shap_long")) %>% purrr::compact() %>% dplyr::bind_rows()
xgb_shap_agg  <- xgb_list %>% purrr::map(~ attr(.x, "shap_agg"))  %>% purrr::compact() %>% dplyr::bind_rows()

saveRDS(xgb_res,       file.path(out_dir, "xgb_metrics.rds"))
saveRDS(xgb_foldC,     file.path(out_dir, "xgb_foldC.rds"))
saveRDS(xgb_shap_long, file.path(out_dir, "xgb_shap_long.rds"))
saveRDS(xgb_shap_agg,  file.path(out_dir, "xgb_shap_agg.rds"))


# PERMANOVA + SHAP
permanova_list <- purrr::imap(
  transforms,
  ~ permanova_r2_cv5(
    df = .x, method_name = .y,
    group_col = EVENT_COL, event_col = EVENT_COL, time_col = TIME_COL,
    harmonize = "zscore", metric = "euclidean",
    permutations = 999,
    K = K_FOLDS, fold_id = global_folds, run_on = "train",
    shapley = TRUE, shapley_M = 64, shapley_pmax = 64,
    seed = CV_SEED, messages = TRUE
  )
)

permanova_res <- dplyr::bind_rows(permanova_list)
permanova_foldR2 <- permanova_list %>% purrr::map(~ attr(.x, "foldR2")) %>% purrr::compact() %>% dplyr::bind_rows()
permanova_shap_long <- permanova_list %>% purrr::map(~ attr(.x, "shap_long")) %>% purrr::compact() %>% dplyr::bind_rows()
permanova_shap_agg  <- permanova_list %>% purrr::map(~ attr(.x, "shap_agg"))  %>% purrr::compact() %>% dplyr::bind_rows()

saveRDS(permanova_res,       file.path(out_dir, "permanova_metrics.rds"))
saveRDS(permanova_foldR2,    file.path(out_dir, "permanova_foldR2.rds"))
saveRDS(permanova_shap_long, file.path(out_dir, "permanova_shap_long.rds"))
saveRDS(permanova_shap_agg,  file.path(out_dir, "permanova_shap_agg.rds"))


out_dir <- file.path("model_result", "shap_outputs")

# CoxNet
coxnet_metrics    <- readRDS(file.path(out_dir, "coxnet_metrics.rds"))
coxnet_shap_agg   <- readRDS(file.path(out_dir, "coxnet_shap_agg.rds"))
coxnet_shap_long  <- readRDS(file.path(out_dir, "coxnet_shap_long.rds"))
coxnet_foldC      <- readRDS(file.path(out_dir, "coxnet_foldC.rds"))

# RSF
rsf_metrics    <- readRDS(file.path(out_dir, "rsf_metrics.rds"))
rsf_shap_agg   <- readRDS(file.path(out_dir, "rsf_shap_agg.rds"))
rsf_shap_long  <- readRDS(file.path(out_dir, "rsf_shap_long.rds"))
rsf_foldC      <- readRDS(file.path(out_dir, "rsf_foldC.rds"))

# Logistic
logit_metrics   <- readRDS(file.path(out_dir, "logit_metrics.rds"))
logit_shap_agg  <- readRDS(file.path(out_dir, "logit_shap_agg.rds"))
logit_shap_long <- readRDS(file.path(out_dir, "logit_shap_long.rds"))
logit_foldC     <- readRDS(file.path(out_dir, "logit_foldC.rds"))

# CatBoost
cb_metrics    <- readRDS(file.path(out_dir, "cb_metrics.rds"))
cb_shap_agg   <- readRDS(file.path(out_dir, "cb_shap_agg.rds"))
cb_shap_long  <- readRDS(file.path(out_dir, "cb_shap_long.rds"))
cb_foldC      <- readRDS(file.path(out_dir, "cb_foldC.rds"))

# DeepSurv
deepsurv_metrics   <- readRDS(file.path(out_dir, "deepsurv_metrics.rds"))
deepsurv_shap_agg  <- readRDS(file.path(out_dir, "deepsurv_shap_agg.rds"))
deepsurv_shap_long <- readRDS(file.path(out_dir, "deepsurv_shap_long.rds"))
deepsurv_foldC     <- readRDS(file.path(out_dir, "deepsurv_foldC.rds"))

# XGBoost Cox
xgb_metrics    <- readRDS(file.path(out_dir, "xgb_metrics.rds"))
xgb_shap_agg   <- readRDS(file.path(out_dir, "xgb_shap_agg.rds"))
xgb_shap_long  <- readRDS(file.path(out_dir, "xgb_shap_long.rds"))
xgb_foldC      <- readRDS(file.path(out_dir, "xgb_foldC.rds"))

# TabPFN
tabpfn_metrics   <- readRDS(file.path(out_dir, "tabpfn_metrics.rds"))
tabpfn_shap_agg  <- readRDS(file.path(out_dir, "tabpfn_shap_agg.rds"))
tabpfn_shap_long <- readRDS(file.path(out_dir, "tabpfn_shap_long.rds"))
tabpfn_foldC     <- readRDS(file.path(out_dir, "tabpfn_foldC.rds"))

# PERMANOVA
permanova_metrics    <- readRDS(file.path(out_dir, "permanova_metrics.rds"))
permanova_foldR2     <- readRDS(file.path(out_dir, "permanova_foldR2.rds"))
permanova_shap_agg   <- readRDS(file.path(out_dir, "permanova_shap_agg.rds"))
permanova_shap_long  <- readRDS(file.path(out_dir, "permanova_shap_long.rds"))

# Combine overall metrics
metrics_all <- dplyr::bind_rows(
  coxnet_metrics, rsf_metrics, logit_metrics,
  cb_metrics, deepsurv_metrics, xgb_metrics, tabpfn_metrics,
  permanova_metrics
)

# Combine fold-wise C only
# PERMANOVA folds kept separate (used for R^2)
# permanova_foldR2
foldC_all <- dplyr::bind_rows(
  coxnet_foldC, rsf_foldC, logit_foldC,
  cb_foldC, deepsurv_foldC, xgb_foldC, tabpfn_foldC
)

# Combine SHAP (long/agg) across all models
shap_long_all <- dplyr::bind_rows(
  coxnet_shap_long, rsf_shap_long, logit_shap_long,
  cb_shap_long, deepsurv_shap_long, xgb_shap_long, tabpfn_shap_long,
  permanova_shap_long
) %>% dplyr::select(-dplyr::any_of(".row_id"))


shap_agg_all <- dplyr::bind_rows(
  coxnet_shap_agg, rsf_shap_agg, logit_shap_agg,
  cb_shap_agg, deepsurv_shap_agg, xgb_shap_agg, tabpfn_shap_agg,
  permanova_shap_agg
)




## Subsample & truncation (restart R first; TabPFN may conflict with tensorflow from DeepSurv)
# Subsample percentages 
p_grid <- c(0.20, 0.40, 0.60, 0.80, 1.00)

# Match truncation times to the same percent levels of follow-up distribution
t_vec  <- df_clr[[TIME_COL]]
t_grid <- as.numeric(stats::quantile(t_vec, probs = p_grid, na.rm = TRUE, type = 7))

# Run all models
run_coxnet_once(p_grid, t_grid, EVENT_COL, TIME_COL)
run_rsf_once(p_grid, t_grid, EVENT_COL, TIME_COL)
run_logit_once(p_grid, t_grid, EVENT_COL, TIME_COL)
run_catboost_once(p_grid, t_grid, EVENT_COL, TIME_COL)
run_tabpfn_once(p_grid, t_grid, EVENT_COL, TIME_COL)
run_deepsurv_once(p_grid, t_grid, EVENT_COL, TIME_COL)
run_xgb_once(p_grid, t_grid, EVENT_COL, TIME_COL)
run_permanova_once(p_grid, t_grid, EVENT_COL, TIME_COL)


# Folder
out_dir <- file.path("model_result", "shap_outputs")

# Read all 4 files per model
load_one_model <- function(key) {
  fold_tag <- if (tolower(key) == "permanova") "foldR2" else "foldC"
  list(
    sub_metrics = readRDS(file.path(out_dir, sprintf("metrics_subsample_%s.rds",  key))),
    sub_folds   = readRDS(file.path(out_dir, sprintf("%s_subsample_%s.rds",       fold_tag, key))),
    tr_metrics  = readRDS(file.path(out_dir, sprintf("metrics_truncation_%s.rds", key))),
    tr_folds    = readRDS(file.path(out_dir, sprintf("%s_truncation_%s.rds",      fold_tag, key)))
  )
}

model_keys <- c(
  coxnet    = "coxnet",
  rsf       = "rsf",
  logit     = "logit",
  catboost  = "catboost",
  tabpfn    = "tabpfn",
  deepsurv  = "deepsurv",
  xgb       = "xgb",
  permanova = "permanova"
)

# Download to the environment
for (nm in names(model_keys)) {
  assign(paste0(nm, "_grids"), load_one_model(model_keys[[nm]]), envir = .GlobalEnv)
}

# Save each model as a single grid
for (nm in names(model_keys)) {
  key   <- model_keys[[nm]]
  grid  <- load_one_model(key)
  assign(paste0(nm, "_grids"), grid, envir = .GlobalEnv)
  saveRDS(grid, file.path(out_dir, sprintf("grid_%s.rds", key)))
}


