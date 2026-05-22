# ------------------------------------------------------------------------------
#
# apply_within_MI.R
#
# This file applies multiple imputation to the BMI and Smoking covariates
# MI is conducted in "within" fashion, meaning that the entire analysis
# (all subsequent YAML actions) are conducted in this file in parallel
# across each of the 10 imputed datasets
#
# Arguments:
#  - cohort - string, defines which of three opensafely cohorts to describe
#             (prevax, vax, unvax)
#  - preex - boolean/string, defines preexisting conditions
#            for the replication preex = FALSE always
#            ("All", TRUE, or FALSE)
#
# Returns:
#  - MANY
#
# Authors: Emma Tarmey
#
# ------------------------------------------------------------------------------


# Load libraries ---------------------------------------------------------------
print("Load libraries")

library(magrittr)
library(mice)
library(here)
library(dplyr)
library(fs)
library(survival)


# Define apply_across_MI output folder ---------------------------------------
print("Creating output/apply_across_MI output folder")

apply_across_MI_dir <- "output/apply_across_MI/"
dir_create(here::here(apply_across_MI_dir))


# Source common functions ------------------------------------------------------
print("Source common functions")

source("analysis/utility.R")


# Specify arguments ------------------------------------------------------------
print("Specify arguments")

args <- commandArgs(trailingOnly = TRUE)
print(length(args))

if (length(args) == 0) {
  # default argument values
  cohort  <- "prevax"
  preex   <- "All"
} else {
  # YAML arguments
  cohort  <- args[[1]]

  # optional argument
  if (length(args) < 2) {
    preex <- "All"
  } else {
    preex <- args[[2]]
  } # allow an empty input for the preex variable
}


# Load data --------------------------------------------------------------------
print("Load data")

df <- readr::read_rds(paste0(
  "output/dataset_clean/input_",
  cohort,
  "_clean.rds"
))

df <- as.data.frame(df)


