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
library(ggplot2)


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


# Sort
# df_ami_nelsonaalen   <- df_ami_nelsonaalen[order(df_ami_nelsonaalen$time, df_ami_nelsonaalen$H0), ]
# df_sahhs_nelsonaalen <- df_sahhs_nelsonaalen[order(df_sahhs_nelsonaalen$time, df_sahhs_nelsonaalen$H0), ]

# print(head(df_ami_nelsonaalen, n = 30))

# stop("?")


# Take random 100-point sample -------------------------------------------------
print("Take random 100-point sample")

sample_size   <- nrow(df_ami_nelsonaalen)
random_sample <- sample(x = c(1:sample_size), size = 100, replace = FALSE)

df_ami_nelsonaalen   <- df_ami_nelsonaalen[random_sample, ]
df_sahhs_nelsonaalen <- df_sahhs_nelsonaalen[random_sample, ]


# Sort by x-axis value (time) --------------------------------------------------
print("Sort by x-axis value (time)")

df_ami_nelsonaalen   <- df_ami_nelsonaalen[order(df_ami_nelsonaalen$time), ]
df_sahhs_nelsonaalen <- df_sahhs_nelsonaalen[order(df_sahhs_nelsonaalen$time), ]


# Convert time to date format --------------------------------------------------
print("Convert time to date format")

df_ami_nelsonaalen$time   <- as.Date(df_ami_nelsonaalen$time,   origin = lubridate::origin)
df_sahhs_nelsonaalen$time <- as.Date(df_sahhs_nelsonaalen$time, origin = lubridate::origin)


# Calculate 95% Confidence Interval --------------------------------------------
print("Calculate 95% Confidence Interval")

# z value is 1.96 for 95% confidence interval
z_value <- qnorm(p=(0.05/2), lower.tail = FALSE)

df_ami_nelsonaalen$upper_CI <- df_ami_nelsonaalen$H0 + (z_value * df_ami_nelsonaalen$se)
df_ami_nelsonaalen$lower_CI <- df_ami_nelsonaalen$H0 - (z_value * df_ami_nelsonaalen$se)

df_sahhs_nelsonaalen$upper_CI <- df_sahhs_nelsonaalen$H0 + (z_value * df_sahhs_nelsonaalen$se)
df_sahhs_nelsonaalen$lower_CI <- df_sahhs_nelsonaalen$H0 - (z_value * df_sahhs_nelsonaalen$se)


# Generate Nelson-Aalen plots -----
print("Generate Nelson-Aalen plots")

png(
  paste0("output/dataset_clean/nelson_aalen_", cohort, "_ami_ggplot.png"),
  width     = 1000,
  height    = 800
)
ggplot(df_ami_nelsonaalen, aes(x=time)) +
  geom_step(mapping = aes(y=H0,       colour = "Nelson Aalen Estimate"),         linetype = 1, linewidth = 0.5) +
  geom_step(mapping = aes(y=lower_CI, colour = "Lower 95% Confidence Interval"), linetype = 2, linewidth = 0.5) +
  geom_step(mapping = aes(y=upper_CI, colour = "Upper 95% Confidence Interval"), linetype = 2, linewidth = 0.5) +
  scale_colour_manual(
    name   = 'Legend',
    breaks = c('Nelson Aalen Estimate', 'Upper 95% Confidence Interval', 'Lower 95% Confidence Interval'),
    values = c('Nelson Aalen Estimate'='black', 'Upper 95% Confidence Interval'='blue', 'Lower 95% Confidence Interval'='red')
  ) +
  labs(
    title    = "Nelson-Aalen Estimator (H0) Plot",
    subtitle = "Acute Myocardial Infarction (ami) (sample 100 points)",
    x        = "Time",
    y        = "Cumulative hazard (H0)"
  ) +
  theme(legend.position = "right")
dev.off()

png(
  paste0("output/dataset_clean/nelson_aalen_", cohort, "_ami.png"),
  width     = 1000,
  height    = 800
)
plot(
  x = df_ami_nelsonaalen$time,
  y = df_ami_nelsonaalen$H0,
  type = "s",
  main = "Nelson-Aalen Estimator (H0) Plot\nAcute Myocardial Infarction (ami) (sample 100 points)",
  xlab = "Time",
  ylab = "Cumulative hazard (H0)"
)
lines(
  x = df_ami_nelsonaalen$time,
  y = df_ami_nelsonaalen$upper_CI,
  type = "s",
  lty = 2,
  col = "blue"
)
lines(
  x = df_ami_nelsonaalen$time,
  y = df_ami_nelsonaalen$lower_CI,
  type = "s",
  lty = 2,
  col = "red"
)
legend(
  "bottomright",
  legend = c("Nelson Aalen Estimate", "Upper 95% Confidence Interval", "Lower 95% Confidence Interval"),
  col    = c("black", "blue", "red"),
  lty    = c(1, 2, 2)
)
dev.off()

png(
  paste0("output/dataset_clean/nelson_aalen_", cohort, "_sahhs_ggplot.png"),
  width     = 1000,
  height    = 800
)
ggplot(df_sahhs_nelsonaalen, aes(x=time)) +
  geom_step(mapping = aes(y=H0,       colour = "Nelson Aalen Estimate"),         linetype = 1, linewidth = 0.5) +
  geom_step(mapping = aes(y=lower_CI, colour = "Lower 95% Confidence Interval"), linetype = 2, linewidth = 0.5) +
  geom_step(mapping = aes(y=upper_CI, colour = "Upper 95% Confidence Interval"), linetype = 2, linewidth = 0.5) +
  scale_colour_manual(
    name   = 'Legend',
    breaks = c('Nelson Aalen Estimate', 'Upper 95% Confidence Interval', 'Lower 95% Confidence Interval'),
    values = c('Nelson Aalen Estimate'='black', 'Upper 95% Confidence Interval'='blue', 'Lower 95% Confidence Interval'='red')
  ) +
  labs(
    title    = "Nelson-Aalen Estimator (H0) Plot",
    subtitle = "Subarachnoid haemorrhage and haemorrhage stroke (sahhs) (sample 100 points)",
    x        = "Time",
    y        = "Cumulative hazard (H0)"
  ) +
  theme(legend.position = "right")
dev.off()

png(
  paste0("output/dataset_clean/nelson_aalen_", cohort, "_sahhs.png"),
  width     = 1000,
  height    = 800
)
plot(
  x = df_sahhs_nelsonaalen$time,
  y = df_sahhs_nelsonaalen$H0,
  type = "s",
  main = "Nelson-Aalen Estimator (H0) Plot\nSubarachnoid haemhorrhage & haemhorrhage stroke (sahhs) (sample 100 points)",
  xlab = "Time",
  ylab = "Cumulative hazard (H0)"
)
lines(
  x = df_sahhs_nelsonaalen$time,
  y = df_sahhs_nelsonaalen$upper_CI,
  type = "s",
  lty = 2,
  col = "blue"
)
lines(
  x = df_sahhs_nelsonaalen$time,
  y = df_sahhs_nelsonaalen$lower_CI,
  type = "s",
  lty = 2,
  col = "red"
)
legend(
  "bottomright",
  legend = c("Nelson Aalen Estimate", "Upper 95% Confidence Interval", "Lower 95% Confidence Interval"),
  col    = c("black", "blue", "red"),
  lty    = c(1, 2, 2)
)
dev.off()
