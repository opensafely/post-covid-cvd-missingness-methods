# ------------------------------------------------------------------------------
#
# table_missingness.R
#
# This file generates the missingness table
#
# Arguments:
#  - cohort - string, defines which of three opensafely cohorts to describe
#             (prevax, vax, unvax)
#
# Returns:
#  - table!
#
# Authors: Emma Tarmey, Venexia Walker, UoB ehrQL Team
#
# ------------------------------------------------------------------------------


# Refresh local R session ------------------------------------------------------
print("Refresh local R session")

rm(list=ls())


# Load libraries ---------------------------------------------------------------
print("Load libraries")

library(magrittr)
library(here)
library(dplyr)
library(fs)


# Define table_missingness output folder ---------------------------------------------------------
print("Creating output/table_missingness output folder")

table_missingness_dir <- "output/table_missingness/"
dir_create(here::here(table_missingness_dir))


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
} else {
  # YAML arguments
  cohort  <- args[[1]]
}


# Load data --------------------------------------------------------------------
print("Load data")

df <- readr::read_rds(paste0(
  "output/dataset_clean/input_",
  cohort,
  "_clean_prehoc.rds"
))

# df <- readr::read_rds(paste0(
#   "output/dataset_definition/input_",
#   cohort,
#   ".rds"
# ))

df <- as.data.frame(df)


# Convert missingness categories to NA -----------------------------------------

df$cov_cat_ethnicity[df$cov_cat_ethnicity == "Missing"] <- NA
df$cov_cat_smoking[df$cov_cat_smoking == "Missing"]     <- NA


# Make table -------------------------------------------------------------------

all_var_names <- c(
  "cov_num_age", "cov_cat_sex", "cov_num_bmi", "cov_cat_ethnicity", "cov_cat_imd",
  "cov_cat_smoking", "cov_bin_carehome", "cov_bin_hcworker", "cov_bin_dementia",
  "cov_bin_liver_disease", "cov_bin_ckd", "cov_bin_cancer", "cov_bin_hypertension",
  "cov_bin_diabetes", "cov_bin_obesity", "cov_bin_copd", "cov_bin_depression", "cov_bin_stroke_all",
  "cov_bin_other_ae", "cov_bin_vte", "cov_bin_hf", "cov_bin_angina", "cov_bin_lipidmed",
  "cov_bin_antiplatelet", "cov_bin_anticoagulant", "cov_bin_cocp", "cov_bin_hrt", "strat_cat_region"
)

all_var_missingness        <- rep(0, length.out = length(all_var_names))
names(all_var_missingness) <- all_var_names

all_var_missingness["cov_num_bmi"]       <- signif(100 * (sum(is.na(df$cov_num_bmi)) / length(df$cov_num_bmi)), digits = 4)
all_var_missingness["cov_cat_ethnicity"] <- signif(100 * (sum(is.na(df$cov_cat_ethnicity)) / length(df$cov_cat_ethnicity)), digits = 4)
all_var_missingness["cov_cat_smoking"]   <- signif(100 * (sum(is.na(df$cov_cat_smoking)) / length(df$cov_cat_smoking)), digits = 4)

missingness_table <- data.frame(
  Covariates  = all_var_names,
  Missingness = all_var_missingness
)


# Save results -----------------------------------------------------------------

write.csv(
  missingness_table,
  paste0(table_missingness_dir, "table_missingness_cohort_", cohort, ".csv"),
  row.names = TRUE
)
