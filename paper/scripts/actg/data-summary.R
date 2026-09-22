#!/usr/bin/env Rscript

# Reviewer-facing entry point for the public ACTG descriptive analysis.

script_argument <- grep(
  "^--file=", commandArgs(trailingOnly = FALSE), value = TRUE
)
script_path <- normalizePath(
  sub("^--file=", "", script_argument[[1L]]), mustWork = TRUE
)
root <- normalizePath(
  file.path(dirname(script_path), "..", "..", ".."),
  winslash = "/"
)
setwd(root)

if (!requireNamespace("gamlss.longitudinal", quietly = TRUE)) {
  stop(
    "Install the checked-out gamlss.longitudinal package before running this script."
  )
}
source(file.path("paper", "R", "05-actg-workflow.R"))
data <- jss_actg_load_data()
summary <- jss_actg_data_summary(data)
print(summary$time_summary)
print(summary$cohort_summary)
print(round(summary$pearson, 3))
print(round(summary$kendall, 3))
print(round(summary$upper_tail, 3))
print(round(summary$lower_tail, 3))
print(summary$missingness)
