#!/usr/bin/env Rscript

fmt_num <- function(x, digits = 3) {
  x <- suppressWarnings(as.numeric(x))
  ifelse(is.na(x), "--", formatC(x, format = "f", digits = digits))
}

safe_read_csv <- function(path) {
  if (!file.exists(path)) return(NULL)
  read.csv(path, check.names = FALSE)
}

tex_escape <- function(x) {
  x <- as.character(x)
  x <- gsub("\\\\", "\\\\textbackslash{}", x)
  x <- gsub("_", "\\\\_", x)
  x
}

tex_parameter <- function(x) {
  out <- switch(
    as.character(x),
    mu = "$\\mu$",
    sigma = "$\\sigma$",
    nu = "$\\nu$",
    tau = "$\\tau$",
    theta = "$\\theta$",
    zeta = "$\\zeta$",
    tex_escape(x)
  )
  out
}

out_dir <- Sys.getenv("OUT_DIR", unset = file.path("results", "bcpe_t_gamlss_comparison"))
fixed <- read.csv(file.path(out_dir, "fixed_effects_bias_rmse_compact.csv"), check.names = FALSE)
smooth <- read.csv(file.path(out_dir, "smooth_integrated_metrics.csv"), check.names = FALSE)
runs <- read.csv(file.path(out_dir, "fit_run_log.csv"), check.names = FALSE)
se_calibration <- safe_read_csv(file.path(out_dir, "fixed_effects_se_calibration.csv"))
joint_metrics <- safe_read_csv(file.path(out_dir, "joint_distribution_metrics_summary.csv"))

fixed$abs_bias_long <- abs(fixed[["bias.gamlss.longitudinal"]])
fixed$abs_bias_g2 <- abs(fixed[["bias.gamlss2"]])
fixed$delta_abs_bias <- fixed$abs_bias_long - fixed$abs_bias_g2
fixed$delta_rmse <- fixed[["rmse.gamlss.longitudinal"]] - fixed[["rmse.gamlss2"]]
fixed <- fixed[order(fixed$parameter, match(fixed$term, c("intercept", "x1", "x2", "t"))), ]

smooth_summary <- aggregate(
  cbind(bias_abs_integrated, irmse) ~ scenario + model + n + d + parameter,
  data = smooth,
  FUN = mean,
  na.rm = TRUE
)
smooth_wide <- reshape(
  smooth_summary,
  idvar = c("scenario", "n", "d", "parameter"),
  timevar = "model",
  direction = "wide"
)
for (nm in c(
  "bias_abs_integrated.gamlss.longitudinal", "irmse.gamlss.longitudinal",
  "bias_abs_integrated.gamlss2", "irmse.gamlss2"
)) {
  if (!nm %in% names(smooth_wide)) smooth_wide[[nm]] <- NA_real_
}
smooth_wide$delta_irmse <- smooth_wide[["irmse.gamlss.longitudinal"]] - smooth_wide[["irmse.gamlss2"]]
smooth_wide <- smooth_wide[order(smooth_wide$parameter), ]

run_summary <- aggregate(
  cbind(logLik, elapsed_sec) ~ scenario + model + n + d,
  data = runs[runs$success == TRUE, ],
  FUN = function(x) c(mean = mean(x, na.rm = TRUE), sd = stats::sd(x, na.rm = TRUE))
)
run_df <- data.frame(
  scenario = run_summary$scenario,
  model = run_summary$model,
  n = run_summary$n,
  d = run_summary$d,
  mean_logLik = run_summary$logLik[, "mean"],
  sd_logLik = run_summary$logLik[, "sd"],
  mean_elapsed_sec = run_summary$elapsed_sec[, "mean"],
  stringsAsFactors = FALSE
)
run_wide <- reshape(run_df, idvar = c("scenario", "n", "d"), timevar = "model", direction = "wide")
run_wide$delta_mean_logLik <- run_wide[["mean_logLik.gamlss.longitudinal"]] - run_wide[["mean_logLik.gamlss2"]]

if (!is.null(se_calibration)) {
  se_wide <- reshape(
    se_calibration,
    idvar = c("scenario", "n", "d", "parameter", "term"),
    timevar = "model",
    direction = "wide"
  )
  for (nm in c(
    "mean_std_error.gamlss.longitudinal", "sd_estimate.gamlss.longitudinal",
    "se_to_empirical_sd.gamlss.longitudinal", "coverage_95.gamlss.longitudinal",
    "mean_std_error.gamlss2", "sd_estimate.gamlss2",
    "se_to_empirical_sd.gamlss2", "coverage_95.gamlss2"
  )) {
    if (!nm %in% names(se_wide)) se_wide[[nm]] <- NA_real_
  }
  se_wide <- se_wide[order(se_wide$parameter, match(se_wide$term, c("x1", "x2", "t"))), ]
}

if (!is.null(joint_metrics)) {
  for (nm in c(
    "mean_logLik", "sd_logLik", "mean_marginal_pit_ks", "mean_marginal_pit_cvm",
    "mean_conditional_pit_ks", "mean_conditional_pit_cvm",
    "mean_abs_adjacent_pit_cor", "mean_abs_conditional_pit_cor",
    "mean_rosenblatt_ks", "mean_rosenblatt_cvm",
    "mean_abs_rosenblatt_lag1_cor", "mean_abs_rosenblatt_normal_lag1_cor",
    "mean_rosenblatt_mean_abs_time_mean", "mean_rosenblatt_normal_mean_abs_time_mean"
  )) {
    if (!nm %in% names(joint_metrics)) joint_metrics[[nm]] <- NA_real_
  }
  joint_wide <- reshape(
    joint_metrics,
    idvar = c("scenario", "n", "d"),
    timevar = "model",
    direction = "wide"
  )
  joint_wide$delta_mean_logLik <- joint_wide[["mean_logLik.gamlss.longitudinal"]] - joint_wide[["mean_logLik.gamlss2"]]
}

fixed_rows <- apply(fixed, 1, function(r) {
  paste(
    tex_parameter(r[["parameter"]]),
    tex_escape(r[["term"]]),
    fmt_num(r[["bias.gamlss.longitudinal"]]),
    fmt_num(r[["rmse.gamlss.longitudinal"]]),
    fmt_num(r[["bias.gamlss2"]]),
    fmt_num(r[["rmse.gamlss2"]]),
    fmt_num(r[["delta_abs_bias"]]),
    fmt_num(r[["delta_rmse"]]),
    sep = " & "
  )
})

smooth_rows <- apply(smooth_wide, 1, function(r) {
  paste(
    tex_parameter(r[["parameter"]]),
    fmt_num(r[["bias_abs_integrated.gamlss.longitudinal"]]),
    fmt_num(r[["irmse.gamlss.longitudinal"]]),
    fmt_num(r[["bias_abs_integrated.gamlss2"]]),
    fmt_num(r[["irmse.gamlss2"]]),
    fmt_num(r[["delta_irmse"]]),
    sep = " & "
  )
})

loglik_rows <- apply(run_wide, 1, function(r) {
  paste(
    tex_escape(r[["scenario"]]),
    r[["n"]],
    r[["d"]],
    fmt_num(r[["mean_logLik.gamlss.longitudinal"]], 1),
    fmt_num(r[["sd_logLik.gamlss.longitudinal"]], 1),
    fmt_num(r[["mean_logLik.gamlss2"]], 1),
    fmt_num(r[["sd_logLik.gamlss2"]], 1),
    fmt_num(r[["delta_mean_logLik"]], 1),
    sep = " & "
  )
})

if (exists("se_wide")) {
  se_rows <- apply(se_wide, 1, function(r) {
    paste(
      tex_parameter(r[["parameter"]]),
      tex_escape(r[["term"]]),
      fmt_num(r[["mean_std_error.gamlss.longitudinal"]]),
      fmt_num(r[["sd_estimate.gamlss.longitudinal"]]),
      fmt_num(r[["se_to_empirical_sd.gamlss.longitudinal"]]),
      fmt_num(r[["coverage_95.gamlss.longitudinal"]], 2),
      fmt_num(r[["mean_std_error.gamlss2"]]),
      fmt_num(r[["sd_estimate.gamlss2"]]),
      fmt_num(r[["se_to_empirical_sd.gamlss2"]]),
      fmt_num(r[["coverage_95.gamlss2"]], 2),
      sep = " & "
    )
  })
}

if (exists("joint_wide")) {
  joint_rows <- apply(joint_wide, 1, function(r) {
    paste(
      tex_escape(r[["scenario"]]),
      r[["n"]],
      r[["d"]],
      fmt_num(r[["mean_logLik.gamlss.longitudinal"]], 1),
      fmt_num(r[["mean_logLik.gamlss2"]], 1),
      fmt_num(r[["delta_mean_logLik"]], 1),
      fmt_num(r[["mean_rosenblatt_cvm.gamlss.longitudinal"]]),
      fmt_num(r[["mean_rosenblatt_cvm.gamlss2"]]),
      fmt_num(r[["mean_abs_rosenblatt_normal_lag1_cor.gamlss.longitudinal"]]),
      fmt_num(r[["mean_abs_rosenblatt_normal_lag1_cor.gamlss2"]]),
      fmt_num(r[["mean_rosenblatt_normal_mean_abs_time_mean.gamlss.longitudinal"]]),
      fmt_num(r[["mean_rosenblatt_normal_mean_abs_time_mean.gamlss2"]]),
      sep = " & "
    )
  })
}

lines <- c(
  paste0("% Auto-generated from ", tex_escape(out_dir), " CSV files."),
  "% Requires \\usepackage{booktabs}.",
  "",
  "\\begin{table}[ht]",
  "\\centering",
  "\\caption{Fixed-effect recovery for the BCPE-$t$ simulation ($n = 500$, $d = 4$, 10 replications). Bias and RMSE are on the linear predictor scale. Positive $\\Delta$ values indicate smaller error for gamlss2 than for gamlss.longitudinal.}",
  "\\label{tab:bcpe-t-fixed-recovery}",
  "\\begin{tabular}{llrrrrrr}",
  "\\toprule",
  "Parameter & Term & \\multicolumn{2}{c}{gamlss.longitudinal} & \\multicolumn{2}{c}{gamlss2} & \\multicolumn{2}{c}{Difference} \\\\",
  "\\cmidrule(lr){3-4} \\cmidrule(lr){5-6} \\cmidrule(lr){7-8}",
  " & & Bias & RMSE & Bias & RMSE & $\\Delta |\\mathrm{Bias}|$ & $\\Delta$ RMSE \\\\",
  "\\midrule",
  paste0(fixed_rows, " \\\\"),
  "\\bottomrule",
  "\\end{tabular}",
  "\\end{table}",
  "",
  "\\begin{table}[ht]",
  "\\centering",
  "\\caption{Integrated smooth-term recovery. Integrated absolute bias and integrated RMSE are averaged over the 10 replications. The standard gamlss2 marginal model has no $\\theta$ smooth, so those entries are not applicable.}",
  "\\label{tab:bcpe-t-smooth-recovery}",
  "\\begin{tabular}{lrrrrr}",
  "\\toprule",
  "Parameter & \\multicolumn{2}{c}{gamlss.longitudinal} & \\multicolumn{2}{c}{gamlss2} & Difference \\\\",
  "\\cmidrule(lr){2-3} \\cmidrule(lr){4-5}",
  " & Int. abs. bias & IRMSE & Int. abs. bias & IRMSE & $\\Delta$ IRMSE \\\\",
  "\\midrule",
  paste0(smooth_rows, " \\\\"),
  "\\bottomrule",
  "\\end{tabular}",
  "\\end{table}",
  "",
  "\\begin{table}[ht]",
  "\\centering",
  "\\caption{Overall log-likelihood comparison across the 10 successful replications. $\\Delta$ log-likelihood is gamlss.longitudinal minus gamlss2.}",
  "\\label{tab:bcpe-t-loglik}",
  "\\begin{tabular}{lrrrrrrr}",
  "\\toprule",
  "Scenario & $n$ & $d$ & \\multicolumn{2}{c}{gamlss.longitudinal} & \\multicolumn{2}{c}{gamlss2} & Difference \\\\",
  "\\cmidrule(lr){4-5} \\cmidrule(lr){6-7}",
  " & & & Mean logLik & SD & Mean logLik & SD & $\\Delta$ mean logLik \\\\",
  "\\midrule",
  paste0(loglik_rows, " \\\\"),
  "\\bottomrule",
  "\\end{tabular}",
  "\\end{table}",
  ""
)

if (exists("se_rows")) {
  lines <- c(
    lines,
    "\\begin{table}[ht]",
    "\\centering",
    "\\caption{Fixed-effect standard-error calibration for non-intercept terms. The SE/SD ratio compares the mean model-based standard error with the empirical standard deviation of the estimates over replications; values near one indicate calibrated standard errors. Coverage is the empirical coverage of nominal 95\\% intervals.}",
    "\\label{tab:bcpe-t-se-calibration}",
    "\\begin{tabular}{llrrrrrrrr}",
    "\\toprule",
    "Parameter & Term & \\multicolumn{4}{c}{gamlss.longitudinal} & \\multicolumn{4}{c}{gamlss2} \\\\",
    "\\cmidrule(lr){3-6} \\cmidrule(lr){7-10}",
    " & & Mean SE & Emp. SD & SE/SD & Coverage & Mean SE & Emp. SD & SE/SD & Coverage \\\\",
    "\\midrule",
    paste0(se_rows, " \\\\"),
    "\\bottomrule",
    "\\end{tabular}",
    "\\end{table}",
    ""
  )
}

if (exists("joint_rows")) {
  lines <- c(
    lines,
    "\\begin{table}[ht]",
    "\\centering",
    "\\caption{Joint-distribution diagnostics based on the Rosenblatt transform. Larger log-likelihood is better. Smaller Rosenblatt Cramer--von Mises statistics and residual normal-score correlations indicate better calibration of the fitted joint distribution.}",
    "\\label{tab:bcpe-t-joint-diagnostics}",
    "\\begin{tabular}{lrrrrrrrrrrr}",
    "\\toprule",
    "Scenario & $n$ & $d$ & \\multicolumn{3}{c}{Log-likelihood} & \\multicolumn{2}{c}{Rosenblatt CvM} & \\multicolumn{2}{c}{$|\\mathrm{corr}(Z_t,Z_{t-1})|$} & \\multicolumn{2}{c}{Mean $|\\bar{Z}_t|$} \\\\",
    "\\cmidrule(lr){4-6} \\cmidrule(lr){7-8} \\cmidrule(lr){9-10} \\cmidrule(lr){11-12}",
    " & & & Long. & gamlss2 & $\\Delta$ & Long. & gamlss2 & Long. & gamlss2 & Long. & gamlss2 \\\\",
    "\\midrule",
    paste0(joint_rows, " \\\\"),
    "\\bottomrule",
    "\\end{tabular}",
    "\\end{table}",
    ""
  )
}

tex_path <- file.path(out_dir, "model_comparison_tables.tex")
writeLines(lines, tex_path, useBytes = TRUE)
cat("Wrote LaTeX tables to ", tex_path, "\n", sep = "")
