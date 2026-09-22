#!/usr/bin/env Rscript

# Reviewer-facing entry point. The targets workflow is preferred because it
# checkpoints screening, fitting, post-fit simulation, and rendering separately.

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
source(file.path("paper", "R", "replication-helpers.R"))
source(file.path("paper", "R", "05-actg-workflow.R"))
Sys.setenv(GAMLSS_LONGITUDINAL_JSS_ROOT = root)
settings <- jss_settings(create = TRUE)
data <- jss_actg_load_data()
summaries <- jss_actg_data_summary(data)
screening <- jss_actg_screen(data, settings)
primary <- jss_actg_fit_primary(data)
families <- jss_actg_family_comparison(data, primary)
sensitivity <- jss_actg_fit_sensitivity(data)
postfit <- jss_actg_postfit(data, primary)
outputs <- jss_run_05_actg(
  settings, data, summaries, screening, primary,
  families, sensitivity, postfit
)
print(summary(primary$value))
message("ACTG outputs written under: ", settings$out_dir)
invisible(outputs)
