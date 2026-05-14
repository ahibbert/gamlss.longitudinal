dirs <- c(
  cg_prplus = "results/bcpe_t_gamlss_comparison_cg_prplus_dlcopdpar_absorbcons_10rep",
  cg_prplus_wolfe = "results/bcpe_t_gamlss_comparison_cg_prplus_wolfe_dlcopdpar_absorbcons_10rep",
  cg_prplus_adaptcap = "results/bcpe_t_gamlss_comparison_cg_prplus_adaptcap_dlcopdpar_absorbcons_10rep"
)

summarise_probe <- function(label, dir) {
  runs <- read.csv(file.path(dir, "fit_run_log.csv"), stringsAsFactors = FALSE)
  joint <- read.csv(file.path(dir, "joint_distribution_metrics_summary.csv"), stringsAsFactors = FALSE)
  fixed <- read.csv(file.path(dir, "fixed_effects_bias_rmse_table.csv"), stringsAsFactors = FALSE)
  smooth <- read.csv(file.path(dir, "smooth_integrated_metrics.csv"), stringsAsFactors = FALSE)
  se <- read.csv(file.path(dir, "fixed_effects_se_calibration.csv"), stringsAsFactors = FALSE)

  long_joint <- joint[joint$model == "gamlss.longitudinal", ]
  fixed_long <- fixed[
    fixed$model == "gamlss.longitudinal" &
      fixed$parameter %in% c("mu", "sigma", "nu", "tau"),
  ]
  nonint <- fixed_long[fixed_long$term != "intercept", ]
  inter <- fixed_long[fixed_long$term == "intercept", ]
  se_long <- se[
    se$model == "gamlss.longitudinal" &
      se$parameter %in% c("mu", "sigma", "nu", "tau") &
      is.finite(se$se_to_empirical_sd),
  ]

  get_irmse <- function(parameter) {
    z <- smooth[
      smooth$model == "gamlss.longitudinal" &
        smooth$parameter == parameter,
      "irmse"
    ]
    mean(z, na.rm = TRUE)
  }

  data.frame(
    run = label,
    longitudinal_success = sum(runs$model == "gamlss.longitudinal" & runs$success),
    mean_logLik = long_joint$mean_logLik,
    rosenblatt_z_lag1 = long_joint$mean_abs_rosenblatt_normal_lag1_cor,
    rosenblatt_cvm = long_joint$mean_rosenblatt_cvm,
    mean_seconds = mean(runs$elapsed_sec[runs$model == "gamlss.longitudinal"], na.rm = TRUE),
    intercept_rmse = mean(inter$rmse, na.rm = TRUE),
    nonintercept_rmse = mean(nonint$rmse, na.rm = TRUE),
    mu_irmse = get_irmse("mu"),
    sigma_irmse = get_irmse("sigma"),
    theta_irmse = get_irmse("theta"),
    se_to_sd = mean(se_long$se_to_empirical_sd, na.rm = TRUE),
    coverage_95 = mean(se_long$coverage_95, na.rm = TRUE)
  )
}

out <- do.call(rbind, Map(summarise_probe, names(dirs), dirs))
write.csv(out, "results/bcpe_t_gamlss_comparison_cg_probe_summary.csv", row.names = FALSE)
print(out)
