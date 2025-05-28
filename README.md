# Survival-analysis

This repository contains code for a Bayesian survival analysis focused on microbiome-derived predictors.

## Overview

The project implements a probabilistic Cox proportional hazards model using the `brms` package in R. Microbial abundance data are used as predictors of time-to-event outcomes, with a focus on the most prevalent taxa.

Key steps include:
- Data preprocessing and transformation (log-transformed microbial counts)
- Feature selection based on taxon prevalence
- Model fitting with censored survival data
- Posterior summarization and interpretation of hazard ratios
- Visualization of survival using Kaplan–Meier curves

## Files

- `survival_tse.rds`: Contains the preprocessed TreeSummarizedExperiment (TSE) object used in the survival analysis.
- `data.R`: Loads and preprocesses microbiome and survival data.
- `model.R`: Fits a Bayesian Cox model using selected microbial predictors.
- `survival.qmd`: Generates the final report (Quarto).
- `main.R`: Wrapper script to execute all of the above in correct order.

## Running the analysis

Please ensure that all required files are located in the working directory.
 
```r
source("main.R")
