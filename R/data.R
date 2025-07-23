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

# Load the preprocessed TreeSummarizedExperiment object
tse <- readRDS("survival_tse.rds")

# Calculate prevalence (how many samples each microbe appears in, detection threshold = 1)
prevalence <- getPrevalence(tse, detection = 1)

# Select the 5 most prevalent microbes
top5_microbes <- names(sort(prevalence, decreasing = TRUE))[1:5]

# Subset the TSE object to keep only the top 5 microbes
tse_top5 <- tse[top5_microbes, ]

# Transform raw count data to log-transformed relative abundances
tse_top5 <- tse_top5 |>
  transformAssay(
    method = "relabundance",
    assay.type = "counts",
    name = "relabundance"
  ) |>
  transformAssay(
    method = "log",
    assay.type = "relabundance",
    pseudocount = 1e-6,
    name = "logcounts"
  )

# Convert the log-transformed assay data to a data frame
df <- as.data.frame(t(assay(tse_top5, "logcounts")))

# Add survival outcome variables to the data frame
df$Event <- colData(tse_top5)$Event
df$Event_time <- colData(tse_top5)$Event_time

# Clean up column names for model fitting
df <- clean_column_names(df)

# Identify predictor variable names (exclude outcome variables)
remove_vars <- c("Event", "Event_time")
predictors <- setdiff(names(df), remove_vars)

# Select the 20 most prevalent microbes for log-ratio analysis
top20_microbes <- names(sort(prevalence, decreasing = TRUE))[1:20]
tse_top20 <- tse[top20_microbes, ]


# Transform counts to log-ratio features
tse_top20 <- tse_top20 |>
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
  ) |>
  transformAssay(
    method = "difference",
    assay.type = "log_abund",
    name = "logratios"
  )


# Convert sparse logratio matrix to dense matrix, then to data.frame
logratio_mat <- as.matrix(assay(altExp(tse_top20, "logratios")))
df_lra_ratios <- as.data.frame(t(logratio_mat))

# Standardize features and add survival metadata
df_lra_ratios <- as.data.frame(scale(df_lra_ratios))
df_lra_ratios$Event <- colData(tse_top20)$Event
df_lra_ratios$Event_time <- colData(tse_top20)$Event_time
df_lra_ratios <- clean_column_names(df_lra_ratios)


# CLR
tse_top20 <- transformAssay(
  tse_top20,
  method = "clr",
  assay.type = "counts",
  pseudocount = 1e-6,
  name = "clr"
)
df_clr <- as.data.frame(t(assay(tse_top20, "clr")))
df_clr$Event <- colData(tse_top20)$Event
df_clr$Event_time <- colData(tse_top20)$Event_time
df_clr <- clean_column_names(df_clr)

# rCLR
tse_top20 <- transformAssay(
  tse_top20,
  method = "rclr",
  assay.type = "counts",
  pseudocount = 1e-6,
  name = "rclr"
)
df_rclr <- as.data.frame(t(assay(tse_top20, "rclr")))
df_rclr$Event <- colData(tse_top20)$Event
df_rclr$Event_time <- colData(tse_top20)$Event_time
df_rclr <- clean_column_names(df_rclr)

# Direct log-transformed abundances (not ratios)
df_abund <- as.data.frame(t(assay(tse_top20, "log_abund")))
df_abund$Event <- colData(tse_top20)$Event
df_abund$Event_time <- colData(tse_top20)$Event_time
df_abund <- clean_column_names(df_abund)
