library(targets)

source(file.path("paper", "R", "replication-helpers.R"))

tar_option_set(
  packages = c("gamlss.longitudinal", "gamlss.dist", "ggplot2"),
  format = "rds"
)

list(
  tar_target(settings, jss_settings()),
  tar_target(coverage_results, jss_run_coverage_suite(settings)),
  tar_target(coverage_summary_csv, jss_write_coverage_summary(coverage_results, settings), format = "file"),
  tar_target(convergence_summary_csv, jss_write_convergence_summary(coverage_results, settings), format = "file"),
  tar_target(runtime_figure, jss_write_runtime_figure(coverage_results, settings), format = "file"),
  tar_target(convergence_figure, jss_write_convergence_figure(coverage_results, settings), format = "file"),
  tar_target(session_info, jss_write_session_info(settings), format = "file"),
  tar_target(manifest, jss_write_manifest(settings), format = "file"),
  tar_target(
    output_hashes,
    jss_write_output_hashes(
      settings,
      manifest,
      coverage_summary_csv,
      convergence_summary_csv,
      runtime_figure,
      convergence_figure,
      session_info
    ),
    format = "file"
  ),
  tar_target(validation, jss_validate_manifest(manifest, output_hashes))
)
