# Survival-analysis

This repository contains code for a Bayesian survival analysis focused on microbiome-derived predictors.

## Overview

The project implements a probabilistic Cox proportional hazards model using the `brms` package in R. Microbial abundance data are used as predictors of time-to-event outcomes, with a focus on the most prevalent taxa.

Key steps include:
- Data preprocessing and transformation
- Feature selection based on taxon prevalence
- Model fitting with censored survival data: Bayesian survival models are fit to the data using different priors (Normal and Horseshoe).
- Posterior summarization and interpretation of hazard ratios
- Visualization of survival using Kaplan–Meier curves
- Probabilistic survival estimation using the Imprecise Dirichlet Process (IDP)
- Construction of centered log-ratio (CLR), robust CLR (rCLR), pairwise log-ratio analysis (LRA), additive log-ratio (ALR), presence/absence (PA), total sum scaling (TSS), log-transformed TSS (logTSS), and arcsin square root (ASIN) and comparing them.
- Multiple survival modeling approaches are compared (CoxPH, Random Survival Forests, XGB-Cox, DeepSurv, and a logistic-model that ignores censoring).

## Files

- `survival_tse.rds`: Contains the preprocessed TreeSummarizedExperiment (TSE) object used in the survival analysis.
- `data.R`: Loads and preprocesses microbiome and survival data.
- `model.R`: Fits a Bayesian Cox model using selected microbial predictors.
- `funcs.R`: Defines reusable functions for data transformation, log-ratio construction and model fitting.
- `survival.qmd`: Generates the final report (Quarto).
- `survival.md`: Generated report.
- `main.R`: Wrapper script to execute all of the above in correct order.

## How to use

Please ensure that all required files are located in the working directory. And then run the following code in R:
 
```r
source("main.R")
