jss_run_01_bcpe_t <- function(settings) {
  module_id <- "01-simulation-bcpe-t"
  source_dir <- file.path(
    settings$root,
    "results",
    "jss-exploratory",
    "01-continuous-bcpe-t",
    "bcpe_t_current_defaults_rep100_comparison_p05_p2_fits"
  )

  artifacts <- data.frame(
    source_file = c(
      "paper_simulation_bcpe_t_fit_characteristics.csv",
      "paper_simulation_bcpe_t_fit_characteristics.tex",
      "paper_simulation_bcpe_t_fixed_parameter_bias_rmse.csv",
      "paper_simulation_bcpe_t_fixed_parameter_bias_rmse.tex",
      "paper_simulation_bcpe_t_gamlss2_variogram_p2_trim_exclusions.csv",
      "paper_simulation_bcpe_t_fixed_effect_recovery.png",
      "paper_simulation_bcpe_t_smooth_effect_recovery.png"
    ),
    output_file = c(
      "01-simulation-bcpe-t-fit-characteristics.csv",
      "01-simulation-bcpe-t-fit-characteristics.tex",
      "01-simulation-bcpe-t-fixed-parameter-bias-rmse.csv",
      "01-simulation-bcpe-t-fixed-parameter-bias-rmse.tex",
      "01-simulation-bcpe-t-gamlss2-variogram-p2-trim-exclusions.csv",
      "01-simulation-bcpe-t-fixed-effect-recovery.png",
      "01-simulation-bcpe-t-smooth-effect-recovery.png"
    ),
    artifact_type = c(
      "table",
      "table",
      "table",
      "table",
      "data",
      "figure",
      "figure"
    ),
    role = c(
      "fit_characteristics_csv",
      "fit_characteristics_latex",
      "parameter_recovery_csv",
      "parameter_recovery_latex",
      "variogram_trim_audit",
      "fixed_effect_recovery_figure",
      "smooth_effect_recovery_figure"
    ),
    stringsAsFactors = FALSE
  )

  jss_copy_final_artifacts(
    settings = settings,
    module_id = module_id,
    title = "BCPE margin with t-copula simulation",
    source_dir = source_dir,
    artifacts = artifacts,
    notes = paste(
      "Final 100-replicate BCPE/t paper artifacts with p = 0.5 and p = 2",
      "variogram diagnostics, SE calibration, held-out predictive scores,",
      "and the audited gamlss2 p = 2 variogram trim."
    )
  )
}
