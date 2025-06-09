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

# Load the preprocessed TreeSummarizedExperiment object
tse <- readRDS("survival_tse.rds")

# Calculate prevalence (how many samples each microbe appears in, detection threshold = 1)
prevalence <- getPrevalence(tse, detection = 1)

# Select the 5 most prevalent microbes
top5_microbes <- names(sort(prevalence, decreasing = TRUE))[1:5]

# Subset the TSE object to keep only the top 5 microbes
tse_top5 <- tse[top5_microbes, ]

# Log-transform the counts to reduce skewness
tse_top5 <- transformAssay(
  tse_top5,
  assay.type = "counts",
  method = "log",
  name = "logcounts"
)

# Convert the log-transformed assay data to a data frame
df <- as.data.frame(t(assay(tse_top5, "logcounts")))

# Add survival outcome variables to the data frame
df$Event <- colData(tse_top5)$Event
df$Event_time <- colData(tse_top5)$Event_time

# Clean up column names for model fitting
df <- clean_column_names(df)

# Remove rows with missing or invalid survival times
dff <- df %>%
  filter(Event_time >= 0) %>%
  drop_na()

# Identify predictor variable names (exclude outcome variables)
remove_vars <- c("Event", "Event_time")
predictors <- setdiff(names(dff), remove_vars)

# Select the 15 most prevalent microbes for log-ratio analysis
prevalence <- getPrevalence(tse, detection = 1)
top15_microbes <- names(sort(prevalence, decreasing = TRUE))[1:15]
tse_top15 <- tse[top15_microbes, ]

# Log-transform counts of top 15 most prevalent microbes
tse_log_abund <- transformAssay(
  tse_top15,
  assay.type = "counts",
  method = "log",
  name = "logcounts"
)

# Convert assay matrix to data frame
df_abund <- as.data.frame(t(assay(tse_log_abund, "logcounts")))

# Add survival outcome variables
df_abund$Event <- colData(tse_log_abund)$Event
df_abund$Event_time <- colData(tse_log_abund)$Event_time

# Clean column names
df_abund <- clean_column_names(df_abund)

# Compute the ratios
df_clr_ratios <- prepare_logratios_from_clr(tse_top15, method = "clr")
df_rclr_ratios <- prepare_logratios_from_clr(tse_top15, method = "rclr")
df_lra_ratios   <- prepare_logratios_from_clr(tse_top15, method = "lra")