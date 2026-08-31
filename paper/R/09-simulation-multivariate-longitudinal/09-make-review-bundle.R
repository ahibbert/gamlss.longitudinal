source(file.path("paper", "R", "09-simulation-multivariate-longitudinal", "00-multivariate-setup.R"))

mvt_load_package()
run_dir <- mvt_read_run_dir()

if (!file.exists(file.path(run_dir, "review_audit.md"))) {
  mvt_audit_run_dir(run_dir)
}

bundle_index <- mvt_write_review_bundle(run_dir)

message("Review bundle index written to: ", bundle_index)
