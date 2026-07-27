library(targets)

source(file.path("paper", "R", "replication-helpers.R"))
source(file.path("paper", "R", "01-simulation-bcpe-t.R"))
source(file.path("paper", "R", "02-simulation-delaporte-clayton.R"))
source(file.path("paper", "R", "03-joint-vs-separate-optimization.R"))
source(file.path("paper", "R", "04-missingness-dropout-sensitivity.R"))
source(file.path("paper", "R", "05-application-lipid.R"))
source(file.path("paper", "R", "06-application-rand-doctor-visits.R"))
source(file.path("paper", "R", "07-gamma-copula-misspecification.R"))
source(file.path("paper", "R", "08-simulation-sensitivity-correlation-misspecification.R"))

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
  tar_target(module_01_bcpe_t, jss_run_01_bcpe_t(settings)),
  tar_target(module_02_delaporte_clayton, jss_run_02_delaporte_clayton(settings)),
  tar_target(module_03_joint_vs_separate, jss_run_03_joint_vs_separate(settings)),
  tar_target(module_04_missingness_dropout, jss_run_04_missingness_dropout(settings)),
  tar_target(module_05_lipid_application, jss_run_05_lipid_application(settings)),
  tar_target(module_06_rand_doctor_visits, jss_run_06_rand_doctor_visits(settings)),
  tar_target(module_07_gamma_copula_misspecification, jss_run_07_gamma_copula_misspecification(settings)),
  tar_target(
    module_08_simulation_sensitivity_correlation_misspecification,
    jss_run_08_simulation_sensitivity_correlation_misspecification(settings)
  ),
  tar_target(
    module_files,
    jss_collect_module_files(
      module_01_bcpe_t,
      module_02_delaporte_clayton,
      module_03_joint_vs_separate,
      module_04_missingness_dropout,
      module_05_lipid_application,
      module_06_rand_doctor_visits,
      module_07_gamma_copula_misspecification,
      module_08_simulation_sensitivity_correlation_misspecification
    ),
    format = "file"
  ),
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
      module_files,
      session_info
    ),
    format = "file"
  ),
  tar_target(validation, jss_validate_manifest(manifest, output_hashes))
)
