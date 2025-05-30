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
clean_names <- names(df)
clean_names <- gsub("__+", "_", clean_names)
clean_names <- gsub("_$", "", clean_names)
clean_names <- gsub("[^A-Za-z0-9_]", "", clean_names)
clean_names <- ifelse(grepl("^[A-Za-z]", clean_names), clean_names, paste0("X", clean_names))
clean_names <- make.unique(clean_names)
names(df) <- clean_names

# Remove rows with missing or invalid survival times
dff <- df %>%
  filter(Event_time >= 0) %>%
  drop_na()

# Identify predictor variable names (exclude outcome variables)
remove_vars <- c("Event", "Event_time")
predictors <- setdiff(names(dff), remove_vars)
