source(file.path("paper", "R", "09-simulation-multivariate-longitudinal", "00-multivariate-setup.R"))

mvt_load_package()
run_dir <- mvt_read_run_dir()

summary_paths <- file.path(
  run_dir,
  mvt_expected_output_files()
)
if (any(!file.exists(summary_paths))) {
  mvt_summarise_results(run_dir)
}

audit <- mvt_audit_run_dir(run_dir)
message(
  "Review audit written to: ",
  file.path(run_dir, "review_audit.csv"),
  " and ",
  file.path(run_dir, "review_audit.md")
)
message(
  "Audit checks: ",
  sum(audit$status == "pass"),
  " pass, ",
  sum(audit$status == "warn"),
  " warn, ",
  sum(audit$status == "fail"),
  " fail."
)
