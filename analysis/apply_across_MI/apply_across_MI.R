# ------------------------------------------------------------------------------
#
# apply_across_MI.R
#
# This file applies multiple imputation to the BMI and Smoking covariates
# MI is conducted in "across" fashion, meaning that the result is a dataframe
# 10x larger than the origina containing all datasets
# 
# Arguments:
#  - cohort - string, defines which of three opensafely cohorts to describe
#             (prevax, vax, unvax)
#  - preex - boolean/string, defines preexisting conditions
#            for the replication preex = FALSE always
#            ("All", TRUE, or FALSE)
#
# Returns:
#  - output/apply_across_MI/input_<cohort>_clean_across_MI_ami.rds
#  - output/apply_across_MI/input_<cohort>_clean_across_MI_sahhs.rds
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


# Applying multiple imputation to BMI and smoking covariates -------------------
print("Applying multiple imputation to BMI and smoking covariates")

# set random seed
set.seed(2026)

# convert all dates to numeric
df <- (df %>% mutate_if(lubridate::is.Date, as.numeric))

# censorship
df$cens_status <- (!is.na(df$cens_date_dereg)) | (!is.na(df$cens_date_death))

# define binary covid19 exposure status
df$cov_bin_covid19 <- !is.na(df$exp_date_covid)

# define binary sahhs outcome
df$cov_bin_sahhs <- !is.na(df$out_date_stroke_sahhs)

# remove dates
date_vars   <- grep("vax_date", colnames(df))
df_dates    <- df[, date_vars]
df_no_vax_dates <- df[, -date_vars]

# re-cast missing smoking to NA
smoking_missing <- df_no_vax_dates$cov_cat_smoking == "Missing"
df_no_vax_dates$cov_cat_smoking[smoking_missing] <- NA

# check missingness of smoking and bmi variables
percent_smoking_missing <- signif(100 * (sum(is.na(df_no_vax_dates$cov_cat_smoking)) / length(df_no_vax_dates$cov_cat_smoking)), digits = 4)
percent_bmi_missing     <- signif(100 * (sum(is.na(df_no_vax_dates$cov_num_bmi))     / length(df_no_vax_dates$cov_num_bmi)),     digits = 4)

print(paste0("The variable smoking is ", percent_smoking_missing, "% missing"))
print(paste0("The variable bmi is ",     percent_bmi_missing,     "% missing"))


# Specify imputation methods for each outcome (ami and sahhs) ------------------
print("Specify imputation methods for each outcome (ami and sahhs)")

# The below ensures that bmi and smoking are handled with specific
# imputation methods and that all other covariates are left alone
# See: https://www.rdocumentation.org/packages/mice/versions/3.17.0/topics/mice
imp_method                    <- rep("", length.out = length(colnames(df_no_vax_dates)))
names(imp_method)             <- colnames(df_no_vax_dates)
imp_method["cov_cat_smoking"] <- "polyreg" # smoking is categorical with 3 levels, Polytomous logistic regression
imp_method["cov_num_bmi"]     <- "norm"    # bmi is numerical, Bayesian linear regression


# Specify imputation formulas for each outcome (ami and sahhs) -----------------
print("Specify imputation formulas for each outcome (ami and sahhs)")

# Specify imputation formulas, exclude variable such as index date
all_var_names <- c(
  "cov_num_age", "cov_cat_sex", "cov_cat_ethnicity",
  "cov_cat_imd", "cov_bin_carehome", "cov_bin_hcworker",
  "cov_bin_dementia", "cov_bin_liver_disease", "cov_bin_ckd",
  "cov_bin_cancer", "cov_bin_hypertension", "cov_bin_diabetes",
  "cov_bin_obesity", "cov_bin_copd",
  "cov_bin_depression", "cov_bin_stroke_all", "cov_bin_other_ae",
  "cov_bin_vte", "cov_bin_hf", "cov_bin_angina",
  "cov_bin_lipidmed", "cov_bin_antiplatelet", "cov_bin_anticoagulant",
  "cov_bin_cocp", "cov_bin_hrt", "cens_date_dereg",
  "strat_cat_region", "vax_cat_jcvi_group",
  "cens_status", "cov_bin_covid19",
  "cov_bin_sahhs", "cov_bin_ami"
)

my_formulas_ami <- list(
  cov_cat_smoking = as.formula(paste0("cov_cat_smoking ~ ", paste(all_var_names, collapse = " + "), " + H0_ami")),
  cov_num_bmi     = as.formula(paste0("cov_num_bmi ~ ",     paste(all_var_names, collapse = " + "), " + H0_ami"))
)

my_formulas_sahhs <- list(
  cov_cat_smoking = as.formula(paste0("cov_cat_smoking ~ ", paste(all_var_names, collapse = " + "), " + H0_sahhs")),
  cov_num_bmi     = as.formula(paste0("cov_num_bmi ~ ",     paste(all_var_names, collapse = " + "), " + H0_sahhs"))
)


# Calculate Nelson-Aalen Estimator for both outcomes (ami and sahhs) -----------
print("Calculate Nelson-Aalen Estimator for both outcomes (ami and sahhs)")

# ami
H0_ami          <- (survfit(Surv(out_date_ami, cens_status) ~ 1, data = df_no_vax_dates) %>% summary(times = unique(df_no_vax_dates$out_date_ami)))
H0_ami          <- H0_ami[c("time", "surv")]
names(H0_ami)   <- c("out_date_ami", "surv")
H0_ami          <- as.data.frame(H0_ami)
df_no_vax_dates <- merge(df_no_vax_dates, H0_ami, all.x = TRUE, by = "out_date_ami")
df_no_vax_dates <- rename(df_no_vax_dates, H0_ami = surv)

# sahhs
H0_stroke_sahhs        <- (survfit(Surv(out_date_stroke_sahhs, cens_status) ~ 1, data = df_no_vax_dates) %>% summary(times = unique(df_no_vax_dates$out_date_stroke_sahhs)))
H0_stroke_sahhs        <- H0_stroke_sahhs[c("time", "surv")]
names(H0_stroke_sahhs) <- c("out_date_stroke_sahhs", "surv")
H0_stroke_sahhs        <- as.data.frame(H0_stroke_sahhs)
df_no_vax_dates        <- merge(df_no_vax_dates, H0_stroke_sahhs, all.x = TRUE, by = "out_date_stroke_sahhs")
df_no_vax_dates        <- rename(df_no_vax_dates, H0_sahhs = surv)


# Applying multiple imputation to BMI and smoking covariates for ami outcome ---
print("Applying multiple imputation to BMI and smoking covariates for ami outcome")

# Apply multiple imputation
num_datasets <- 10
imp_ami <- mice::mice(
  data       = df_no_vax_dates,
  m          = num_datasets,
  maxit      = 20,
  formulas   = my_formulas_ami,
  imp_method = unname(imp_method)
)

df_post_imputation_ami <- mice::complete(
  imp_ami,
  action  = "long",
  include = FALSE
)

df_post_imputation_ami <- subset(
  df_post_imputation_ami,
  select = -c(.imp, .id)
)

df_dates_stacked_ami <- rbind(
  df_dates, df_dates, df_dates, df_dates, df_dates,
  df_dates, df_dates, df_dates, df_dates, df_dates
)

df_post_imputation_ami <- cbind(
  df_post_imputation_ami, df_dates_stacked_ami
)


# Applying multiple imputation to BMI and smoking covariates for sahhs outcome ---
print("Applying multiple imputation to BMI and smoking covariates for sahhs outcome")

# Apply multiple imputation
num_datasets <- 10
imp_sahhs <- mice::mice(
  data       = df_no_vax_dates,
  m          = num_datasets,
  maxit      = 20,
  formulas   = my_formulas_sahhs,
  imp_method = unname(imp_method)
)

df_post_imputation_sahhs <- mice::complete(
  imp_sahhs,
  action  = "long",
  include = FALSE
)

df_post_imputation_sahhs <- subset(
  df_post_imputation_sahhs,
  select = -c(.imp, .id)
)

df_dates_stacked_sahhs <- rbind(
  df_dates, df_dates, df_dates, df_dates, df_dates,
  df_dates, df_dates, df_dates, df_dates, df_dates
)

df_post_imputation_sahhs <- cbind(
  df_post_imputation_sahhs, df_dates_stacked_sahhs
)


# Remove now unused level 'missing' from smoking covariate for both outcomes (ami and sahhs) ---
print("Remove now unused level 'missing' from smoking covariate")

df_post_imputation_ami$cov_cat_smoking <- factor(
  df_post_imputation_ami$cov_cat_smoking,
  levels = levels(droplevels(df_post_imputation_ami$cov_cat_smoking))
)

df_post_imputation_sahhs$cov_cat_smoking <- factor(
  df_post_imputation_sahhs$cov_cat_smoking,
  levels = levels(droplevels(df_post_imputation_sahhs$cov_cat_smoking))
)


# Re-assign unique identifiers to imputed dataset for both outcomes (ami and sahhs)  ---
print("Re-assign unique identifiers to imputed dataset for both outcomes (ami and sahhs) ")

# NB: IDs are identical for both ami and sahhs, hence this is only calculated once

patient_id_1 <- as.numeric(unique(df_post_imputation_ami$patient_id))
shift        <- max(patient_id_1) + 1

patient_id_2  <- patient_id_1 + shift
patient_id_3  <- patient_id_2 + shift
patient_id_4  <- patient_id_3 + shift
patient_id_5  <- patient_id_4 + shift
patient_id_6  <- patient_id_5 + shift
patient_id_7  <- patient_id_6 + shift
patient_id_8  <- patient_id_7 + shift
patient_id_9  <- patient_id_8 + shift
patient_id_10 <- patient_id_9 + shift

new_patient_id <- c(
  patient_id_1, patient_id_2, patient_id_3, patient_id_4, patient_id_5,
  patient_id_6, patient_id_7, patient_id_8, patient_id_9, patient_id_10
)

df_post_imputation_ami$patient_id   <- new_patient_id
df_post_imputation_sahhs$patient_id <- new_patient_id


# Save data after 'across' multiple imputation for both outcomes (ami and sahhs)  ---
print("Save data after 'across' multiple imputation for both outcomes (ami and sahhs) ")

saveRDS(
  df_post_imputation_ami,
  file = paste0(apply_across_MI_dir, "input_", cohort, "_clean_across_MI_ami.rds"),
  compress = TRUE
)

saveRDS(
  df_post_imputation_sahhs,
  file = paste0(apply_across_MI_dir, "input_", cohort, "_clean_across_MI_sahhs.rds"),
  compress = TRUE
)
