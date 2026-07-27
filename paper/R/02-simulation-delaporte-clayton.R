jss_run_02_delaporte_clayton <- function(settings) {
  module_id <- "02-simulation-delaporte-clayton"
  source_dir <- file.path(
    settings$root,
    "results",
    "jss-exploratory",
    "02-discrete-delaporte-clayton",
    "nbi_highsignal_highcorr_n500_rep100_p05_p2_fits"
  )

  artifacts <- data.frame(
    source_file = c(
      "paper_simulation_nbi_clayton_highsignal_fit_characteristics.csv",
      "paper_simulation_nbi_clayton_highsignal_fit_characteristics.tex",
      "paper_simulation_nbi_clayton_highsignal_fixed_parameter_bias_rmse.csv",
      "paper_simulation_nbi_clayton_highsignal_fixed_parameter_bias_rmse.tex",
      "paper_simulation_nbi_clayton_highsignal_convergence_counts.csv",
      "predictive_scores_summary.csv",
      "joint_distribution_metrics_summary.csv",
      "paper_simulation_nbi_clayton_highsignal_fixed_effect_recovery.png",
      "paper_simulation_nbi_clayton_highsignal_smooth_effect_recovery.png",
      "paper_simulation_nbi_clayton_highsignal_runtime_comparison.png",
      "paper_simulation_nbi_clayton_highsignal_iteration_comparison.png",
      "nbi_sigma_compare_recovery.png"
    ),
    output_file = c(
      "02-simulation-delaporte-clayton-fit-characteristics.csv",
      "02-simulation-delaporte-clayton-fit-characteristics.tex",
      "02-simulation-delaporte-clayton-fixed-parameter-bias-rmse.csv",
      "02-simulation-delaporte-clayton-fixed-parameter-bias-rmse.tex",
      "02-simulation-delaporte-clayton-convergence-counts.csv",
      "02-simulation-delaporte-clayton-predictive-scores-summary.csv",
      "02-simulation-delaporte-clayton-joint-distribution-metrics-summary.csv",
      "02-simulation-delaporte-clayton-fixed-effect-recovery.png",
      "02-simulation-delaporte-clayton-smooth-effect-recovery.png",
      "02-simulation-delaporte-clayton-runtime-comparison.png",
      "02-simulation-delaporte-clayton-iteration-comparison.png",
      "02-simulation-delaporte-clayton-sigma-compare-recovery.png"
    ),
    artifact_type = c(
      "table",
      "table",
      "table",
      "table",
      "table",
      "data",
      "data",
      "figure",
      "figure",
      "figure",
      "figure",
      "figure"
    ),
    role = c(
      "fit_characteristics_csv",
      "fit_characteristics_latex",
      "parameter_recovery_csv",
      "parameter_recovery_latex",
      "convergence_counts",
      "predictive_scores_summary",
      "joint_distribution_metrics_summary",
      "fixed_effect_recovery_figure",
      "smooth_effect_recovery_figure",
      "runtime_comparison_figure",
      "iteration_comparison_figure",
      "sigma_compare_recovery_figure"
    ),
    stringsAsFactors = FALSE
  )

  jss_copy_final_artifacts(
    settings = settings,
    module_id = module_id,
    title = "Negative binomial type I margin with Clayton copula simulation",
    source_dir = source_dir,
    artifacts = artifacts,
    notes = paste(
      "Final 100-replicate high-signal, high-correlation NBI/Clayton paper",
      "artifacts with p = 0.5 and p = 2 variogram diagnostics, SE calibration,",
      "held-out predictive scores, and Rosenblatt residual diagnostics."
    )
  )
}
