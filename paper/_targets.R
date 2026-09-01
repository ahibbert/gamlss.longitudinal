library(targets)

source(file.path("paper", "R", "replication-helpers.R"))
source(file.path("paper", "R", "phase2-evidence-contracts.R"))
source(file.path("paper", "R", "main-recovery-evidence.R"))
source(file.path("paper", "R", "01-simulation-bcpe-t.R"))
source(file.path("paper", "R", "02-simulation-delaporte-clayton.R"))
source(file.path("paper", "R", "03-joint-vs-separate-optimization.R"))
source(file.path("paper", "R", "missingness-study-helpers.R"))
source(file.path("paper", "R", "04-missingness-dropout-sensitivity.R"))
source(file.path("paper", "R", "07-gamma-copula-misspecification.R"))
source(file.path("paper", "R", "09-simulation-multivariate-longitudinal", "00-multivariate-setup.R"))
source(file.path("paper", "R", "10-fit-scaling.R"))
source(file.path("paper", "R", "phase2-central-integration.R"))
source(file.path("paper", "R", "phase2-paper-evidence.R"))
source(file.path("paper", "R", "public-paper-producers.R"))

tar_option_set(packages = c("gamlss.longitudinal", "gamlss.dist", "ggplot2"), format = "rds")
profile <- Sys.getenv("GAMLSS_LONGITUDINAL_JSS_PROFILE", "smoke")
if (identical(profile, "expanded")) profile <- "paper"

common <- list(
  tar_target(settings, {
    value <- jss_settings()
    value$target_integrated <- TRUE
    value
  }, cue = tar_cue(mode = "always")),
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
    tar_target(module_phase2_main_recovery, { public_input_files; jss_run_phase2_main_recovery(settings) }),
    tar_target(module_03_joint_vs_separate, { public_input_files; jss_run_phase2_optimizer(settings) }),
    tar_target(module_04_missingness, { public_input_files; jss_run_04_missingness_dropout(settings) }),
    tar_target(module_07_copula_misspecification, {
      public_input_files
      jss_run_07_from_public_results(settings)
    }),
    tar_target(module_09_multivariate_benchmark, { public_input_files; jss_run_phase2_multivariate_benchmark(settings) }),
    tar_target(module_10_fit_scaling, { public_input_files; jss_run_10_fit_scaling(settings) }),
    tar_target(phase2_production_modules, {
      registered <- list(
        main_recovery = module_phase2_main_recovery,
        optimizer = module_03_joint_vs_separate,
        missingness = module_04_missingness,
        copula_selection = module_07_copula_misspecification,
        multivariate_benchmark = module_09_multivariate_benchmark,
        fit_scaling = module_10_fit_scaling
      )
      jss_phase2_validate_production_modules(settings, registered)
      registered
    }),
    tar_target(phase2_claim_evidence, {
      phase2_production_modules
      jss_write_phase2_claim_evidence(settings)
    }, format = "file"),
    tar_target(phase2_gate_audit, jss_write_phase2_gate_audit(
      settings, phase2_production_modules, phase2_claim_evidence
    ), format = "file")
  )
}

output_target <- if (identical(profile, "smoke")) {
  tar_target(public_outputs, list(public_workflow_figures, coverage_summary_csv, convergence_summary_csv, fit_events))
} else {
  tar_target(public_outputs, list(public_workflow_figures, module_phase2_main_recovery, module_03_joint_vs_separate, module_04_missingness,
    module_07_copula_misspecification, module_09_multivariate_benchmark,
    module_10_fit_scaling, phase2_claim_evidence, phase2_gate_audit))
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
