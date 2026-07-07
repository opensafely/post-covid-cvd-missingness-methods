# ------------------------------------------------------------------------------
#
# nelson_aalen_plots.R
#
# This file generates plots of the nelson_aalen estimator for both outcomes
#
# Arguments:
#  - cohort - string, defines which of three opensafely cohorts to describe
#             (prevax, vax, unvax)
#
# Returns:
#  - plots!
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

df_ami_nelsonaalen   <- readRDS(paste0("output/dataset_clean/nelson_aalen_", cohort, "_ami.rds"))
df_sahhs_nelsonaalen <- readRDS(paste0("output/dataset_clean/nelson_aalen_", cohort, "_sahhs.rds"))


# Take random 100-point sample -------------------------------------------------
print("Take random 100-point sample")

sample_size   <- nrow(df_ami_nelsonaalen)
random_sample <- sample(x = c(1:sample_size), size = 100, replace = FALSE)

df_ami_nelsonaalen   <- df_ami_nelsonaalen[random_sample, ]
df_sahhs_nelsonaalen <- df_sahhs_nelsonaalen[random_sample, ]


# Generate Nelson-Aalen plots -----
print("Generate Nelson-Aalen plots")

png(
  paste0("output/dataset_clean/nelson_aalen_", cohort, "_ami.png"),
  width     = 800,
  height    = 600
)
plot(
  x = df_ami_nelsonaalen$time,
  y = df_ami_nelsonaalen$H0,
  main = "Nelson-Aalen Estimator (H0) Plot \nfor Acute Myocardial Infarction (ami)",
  xlab = "Time",
  ylab = "Cumulative hazard (H0)"
)
dev.off()

png(
  paste0("output/dataset_clean/nelson_aalen_", cohort, "_sahhs.png"),
  width     = 800,
  height    = 600
)
plot(
  x = df_sahhs_nelsonaalen$time,
  y = df_sahhs_nelsonaalen$H0,
  main = "Nelson-Aalen Estimator (H0) Plot for \nSubarachnoid haemhorrhage and haemhorrhage stroke (sahhs)",
  xlab = "Time",
  ylab = "Cumulative hazard (H0)"
)
dev.off()
