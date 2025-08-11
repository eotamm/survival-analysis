# Load required packages
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
library(cmdstanr)
# Install CmdStan
if (is.null(cmdstanr::cmdstan_path())) {
  message("CmdStan is not installed, installing now...")
  cmdstanr::install_cmdstan()
}

set.seed(1)

# Load the TreeSummarizedExperiment object
tse <- readRDS("survival_tse.rds")

# Apply relative abundance and log transformation to raw count data
tse <- tse |>
  transformAssay(
    method = "relabundance",
    assay.type = "counts",
    name = "relabundance"
  ) |>
  transformAssay(
    method = "log",
    assay.type = "relabundance",
    pseudocount = 1e-6,
    name = "log_abund"
  )

# Extract log-transformed relative abundance matrix and convert to data frame
df <- as.data.frame(t(assay(tse, "log_abund")))

# Add survival outcome metadata (event indicator and survival time)
df$Event <- colData(tse)$Event
df$Event_time <- colData(tse)$Event_time

# Clean column names
df <- clean_column_names(df)

# Top 20 features for all transformations
df_clr     <- extract_top_features_by_transformation(tse, method = "clr",     top_n = 20)
df_rclr    <- extract_top_features_by_transformation(tse, method = "rclr",    top_n = 20)
df_logabund <- extract_top_features_by_transformation(tse, method = "log_abund", top_n = 20)
df_lra     <- extract_top_features_by_transformation(tse, method = "lra",     top_n = 20)
df_pa      <- extract_top_features_by_transformation(tse, method = "pa",      top_n = 20)
df_tss     <- extract_top_features_by_transformation(tse, method = "tss",     top_n = 20)
df_logtss  <- extract_top_features_by_transformation(tse, method = "logtss",  top_n = 20)
df_asin    <- extract_top_features_by_transformation(tse, method = "asin",    top_n = 20)
df_alr     <- extract_top_features_by_transformation(tse, method = "alr",     top_n = 20)

