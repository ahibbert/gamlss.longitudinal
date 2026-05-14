dirs <- c(
  rs_separate = "results/bcpe_t_gamlss_comparison_rs_separate_absorbcons_100rep",
  rs_dlcopdpar = "results/bcpe_t_gamlss_comparison_rs_dlcopdpar_absorbcons_100rep",
  cg_none = "results/bcpe_t_gamlss_comparison_cg_dlcopdpar_absorbcons_delta025_100rep",
  cg_prplus = "results/bcpe_t_gamlss_comparison_cg_prplus_dlcopdpar_absorbcons_10rep",
  cg_prplus_wolfe = "results/bcpe_t_gamlss_comparison_cg_prplus_wolfe_dlcopdpar_absorbcons_10rep",
  cg_prplus_adaptcap = "results/bcpe_t_gamlss_comparison_cg_prplus_adaptcap_dlcopdpar_absorbcons_10rep"
)

summarise_dir <- function(label, dir) {
  runs <- read.csv(file.path(dir, "fit_run_log.csv"), stringsAsFactors = FALSE)
  fixed <- read.csv(file.path(dir, "fixed_effects_bias_rmse_table.csv"), stringsAsFactors = FALSE)
  smooth_by_rep <- read.csv(file.path(dir, "smooth_integrated_metrics.csv"), stringsAsFactors = FALSE)
  joint_by_rep_path <- file.path(dir, "joint_distribution_metrics_by_rep.csv")
  joint_by_rep <- read.csv(joint_by_rep_path, stringsAsFactors = FALSE)

  runs <- runs[runs$rep <= 10, ]
  joint_by_rep <- joint_by_rep[joint_by_rep$rep <= 10, ]
  smooth_by_rep <- smooth_by_rep[smooth_by_rep$rep <= 10, ]

  long_runs <- runs[runs$model == "gamlss.longitudinal" & runs$success, ]
  long_joint <- joint_by_rep[joint_by_rep$model == "gamlss.longitudinal", ]
  fixed_long <- fixed[
    fixed$model == "gamlss.longitudinal" &
      fixed$parameter %in% c("mu", "sigma", "nu", "tau"),
  ]
  inter <- fixed_long[fixed_long$term == "intercept", ]
  nonint <- fixed_long[fixed_long$term != "intercept", ]

  get_irmse <- function(parameter) {
    z <- smooth_by_rep[
      smooth_by_rep$model == "gamlss.longitudinal" &
        smooth_by_rep$parameter == parameter,
      "irmse"
    ]
    mean(z, na.rm = TRUE)
  }

  data.frame(
    run = label,
    success = nrow(long_runs),
    mean_logLik = mean(long_runs$logLik, na.rm = TRUE),
    rosenblatt_z_lag1 = mean(long_joint$abs_rosenblatt_normal_lag1_cor, na.rm = TRUE),
    rosenblatt_cvm = mean(long_joint$rosenblatt_cvm, na.rm = TRUE),
    mean_seconds = mean(long_runs$elapsed_sec, na.rm = TRUE),
    intercept_rmse = mean(inter$rmse, na.rm = TRUE),
    nonintercept_rmse = mean(nonint$rmse, na.rm = TRUE),
    mu_irmse = get_irmse("mu"),
    sigma_irmse = get_irmse("sigma"),
    theta_irmse = get_irmse("theta")
  )
}

out <- do.call(rbind, Map(summarise_dir, names(dirs), dirs))
write.csv(out, "results/bcpe_t_gamlss_comparison_cg_probe_first10_summary.csv", row.names = FALSE)
print(out)
