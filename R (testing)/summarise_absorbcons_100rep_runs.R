dirs <- c(
  rs_separate = "results/bcpe_t_gamlss_comparison_rs_separate_absorbcons_100rep",
  rs_dlcopdpar = "results/bcpe_t_gamlss_comparison_rs_dlcopdpar_absorbcons_100rep",
  cg_dlcopdpar = "results/bcpe_t_gamlss_comparison_cg_dlcopdpar_absorbcons_delta025_100rep",
  cg_dlcopdpar_lambda1 =
    "results/bcpe_t_gamlss_comparison_cg_dlcopdpar_absorbcons_delta025_lambda1_100rep"
)

read_result <- function(dir, file) {
  read.csv(file.path(dir, file), stringsAsFactors = FALSE)
}

get_smooth <- function(smooth, model, parameter, column) {
  value <- smooth[smooth$model == model & smooth$parameter == parameter, column]
  if (length(value) == 0L) {
    return(NA_real_)
  }
  value[1]
}

summarise_dir <- function(label, dir) {
  run <- read_result(dir, "fit_run_log.csv")
  joint <- read_result(dir, "joint_distribution_metrics_summary.csv")
  fixed <- read_result(dir, "fixed_effects_bias_rmse_table.csv")
  smooth <- read_result(dir, "smooth_integrated_metrics.csv")
  se <- read_result(dir, "fixed_effects_se_calibration.csv")

  long_joint <- joint[joint$model == "gamlss.longitudinal", ]
  gamlss2_joint <- joint[joint$model == "gamlss2", ]

  marginal_parameters <- c("mu", "sigma", "nu", "tau")
  long_fixed <- fixed[
    fixed$model == "gamlss.longitudinal" &
      fixed$parameter %in% marginal_parameters,
  ]
  gamlss2_fixed <- fixed[
    fixed$model == "gamlss2" &
      fixed$parameter %in% marginal_parameters,
  ]

  long_nonintercept <- long_fixed[long_fixed$term != "intercept", ]
  gamlss2_nonintercept <- gamlss2_fixed[gamlss2_fixed$term != "intercept", ]
  long_intercept <- long_fixed[long_fixed$term == "intercept", ]

  se_nonintercept <- se[
    se$term != "intercept" &
      se$parameter %in% marginal_parameters &
      is.finite(se$se_to_empirical_sd),
  ]
  long_se <- se_nonintercept[se_nonintercept$model == "gamlss.longitudinal", ]
  gamlss2_se <- se_nonintercept[se_nonintercept$model == "gamlss2", ]

  data.frame(
    run = label,
    longitudinal_success = sum(
      run$model == "gamlss.longitudinal" & run$success
    ),
    gamlss2_success = sum(run$model == "gamlss2" & run$success),
    longitudinal_logLik = long_joint$mean_logLik,
    gamlss2_logLik = gamlss2_joint$mean_logLik,
    delta_logLik = long_joint$mean_logLik - gamlss2_joint$mean_logLik,
    longitudinal_rosenblatt_z_lag1 =
      long_joint$mean_abs_rosenblatt_normal_lag1_cor,
    gamlss2_rosenblatt_z_lag1 =
      gamlss2_joint$mean_abs_rosenblatt_normal_lag1_cor,
    longitudinal_rosenblatt_cvm = long_joint$mean_rosenblatt_cvm,
    gamlss2_rosenblatt_cvm = gamlss2_joint$mean_rosenblatt_cvm,
    longitudinal_mean_seconds = mean(
      run$elapsed_sec[run$model == "gamlss.longitudinal"],
      na.rm = TRUE
    ),
    gamlss2_mean_seconds = mean(
      run$elapsed_sec[run$model == "gamlss2"],
      na.rm = TRUE
    ),
    longitudinal_margin_intercept_rmse = mean(long_intercept$rmse, na.rm = TRUE),
    longitudinal_margin_nonintercept_rmse =
      mean(long_nonintercept$rmse, na.rm = TRUE),
    gamlss2_margin_nonintercept_rmse =
      mean(gamlss2_nonintercept$rmse, na.rm = TRUE),
    longitudinal_mu_smooth_irmse =
      get_smooth(smooth, "gamlss.longitudinal", "mu", "irmse"),
    gamlss2_mu_smooth_irmse =
      get_smooth(smooth, "gamlss2", "mu", "irmse"),
    longitudinal_sigma_smooth_irmse =
      get_smooth(smooth, "gamlss.longitudinal", "sigma", "irmse"),
    gamlss2_sigma_smooth_irmse =
      get_smooth(smooth, "gamlss2", "sigma", "irmse"),
    longitudinal_theta_smooth_irmse =
      get_smooth(smooth, "gamlss.longitudinal", "theta", "irmse"),
    longitudinal_se_to_sd = mean(long_se$se_to_empirical_sd, na.rm = TRUE),
    gamlss2_se_to_sd = mean(gamlss2_se$se_to_empirical_sd, na.rm = TRUE),
    longitudinal_coverage_95 = mean(long_se$coverage_95, na.rm = TRUE),
    gamlss2_coverage_95 = mean(gamlss2_se$coverage_95, na.rm = TRUE)
  )
}

summary_table <- do.call(
  rbind,
  Map(summarise_dir, names(dirs), dirs)
)

output_csv <- "results/bcpe_t_gamlss_comparison_absorbcons_100rep_summary.csv"
write.csv(summary_table, output_csv, row.names = FALSE)

print(summary_table)
