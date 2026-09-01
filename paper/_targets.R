library(targets)

source(file.path("paper", "R", "replication-helpers.R"))
source(file.path("paper", "R", "01-simulation-bcpe-t.R"))
source(file.path("paper", "R", "02-simulation-delaporte-clayton.R"))
source(file.path("paper", "R", "03-joint-vs-separate-optimization.R"))
source(file.path("paper", "R", "04-missingness-dropout-sensitivity.R"))
source(file.path("paper", "R", "07-gamma-copula-misspecification.R"))
source(file.path("paper", "R", "08-simulation-sensitivity-correlation-misspecification.R"))
source(file.path("paper", "R", "public-paper-producers.R"))

tar_option_set(packages = c("gamlss.longitudinal", "gamlss.dist", "ggplot2"), format = "rds")
profile <- Sys.getenv("GAMLSS_LONGITUDINAL_JSS_PROFILE", "smoke")
if (identical(profile, "expanded")) profile <- "paper"

common <- list(
  tar_target(settings, jss_settings(), cue = tar_cue(mode = "always")),
  tar_target(public_workflow_figures, jss_run_public_workflow_figures(settings)),
  tar_target(
    session_info,
    jss_write_session_info(settings),
    format = "file",
    cue = tar_cue(mode = "always")
  )
)
if (!identical(profile, "smoke")) {
  common <- c(common, list(tar_target(
    public_input_files,
    list.files(settings$public_data_dir, full.names = TRUE, recursive = TRUE, pattern = "[.](csv|tex|md)$"),
    format = "file"
  )))
}

modules <- if (identical(profile, "smoke")) {
  list(
    tar_target(coverage_results, jss_run_coverage_suite(settings)),
    tar_target(coverage_summary_csv, jss_write_coverage_summary(coverage_results, settings), format = "file"),
    tar_target(convergence_summary_csv, jss_write_convergence_summary(coverage_results, settings), format = "file"),
    tar_target(fit_events, jss_write_fit_event_audit(coverage_results, settings), format = "file")
  )
} else {
  list(
    tar_target(module_01_bcpe_t, { public_input_files; jss_run_01_bcpe_t(settings) }),
    tar_target(module_02_nbi_clayton, { public_input_files; jss_run_02_delaporte_clayton(settings) }),
    tar_target(module_03_joint_vs_separate, { public_input_files; jss_run_03_joint_vs_separate(settings) }),
    tar_target(module_04_missingness, { public_input_files; jss_run_04_missingness_dropout(settings) }),
    tar_target(module_07_copula_misspecification, {
      public_input_files
      if (identical(settings$profile, "full")) jss_run_07_gamma_copula_misspecification(settings, stage = "full") else jss_run_07_from_public_results(settings)
    }),
    tar_target(module_08_correlation_misspecification, { public_input_files; jss_run_08_public(settings) })
  )
}

output_target <- if (identical(profile, "smoke")) {
  tar_target(public_outputs, list(public_workflow_figures, coverage_summary_csv, convergence_summary_csv, fit_events))
} else {
  tar_target(public_outputs, list(public_workflow_figures, module_01_bcpe_t, module_02_nbi_clayton,
    module_03_joint_vs_separate, module_04_missingness,
    module_07_copula_misspecification, module_08_correlation_misspecification))
}

tolerance_target <- if (identical(profile, "full")) {
  tar_target(full_tolerance_audit, { public_outputs; jss_validate_full_tolerances(settings) }, format = "file")
} else {
  tar_target(full_tolerance_audit, TRUE)
}

final <- list(
  output_target,
  tar_target(manifest, jss_write_manifest(settings), format = "file"),
  tar_target(artifact_files, {
    public_outputs
    x <- utils::read.csv(manifest, stringsAsFactors = FALSE)
    paths <- jss_manifest_output_files(x, manifest)
    paths <- paths[x$access == "public" & x$publication_status == "active"]
    paths[file.exists(paths)]
  }, format = "file"),
  tolerance_target,
  tar_target(output_hashes, { session_info; artifact_files; full_tolerance_audit; jss_write_output_hashes(settings, manifest) }, format = "file"),
  tar_target(validation, { full_tolerance_audit; jss_validate_manifest(manifest, output_hashes) })
)
c(common, modules, final)
