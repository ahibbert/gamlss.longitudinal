# Consolidated support module for reviewer navigation.
# Source blocks are grouped mechanically from the files named below.

# ---- benchmark-adoption-scenarios.R ----

#' Return named adoption benchmark scenarios

#'

#' `adoption_benchmark_scenarios()` defines a small, opinionated benchmark plan

#' for comparing `gamlss.longitudinal` with GEE/GLMM/GAM defaults. The scenarios

#' are designed to be executable through [run_adoption_benchmarks()] while also

#' documenting the applied claim each scenario is meant to test.

#'

#' @param scenarios Optional character vector of scenario names to keep.

#'

#' @return A data frame with one row per benchmark scenario.

#' @export

adoption_benchmark_scenarios <- function(scenarios = NULL) {
  out <- data.frame(
    scenario = c(
      "gaussian_heteroskedastic",
      "gamma_positive",
      "poisson_count",
      "time_varying_dependence",
      "missing_visits"
    ),
    label = c(
      "Gaussian outcome with heteroskedasticity",
      "Positive skewed outcome",
      "Count outcome with longitudinal dependence",
      "Gaussian outcome with time-varying dependence",
      "Gaussian outcome with missing visits"
    ),
    family = c("NO", "GA", "PO", "NO", "NO"),
    copula = c("N", "N", "N", "N", "N"),
    design = c("scale", "covariate", "covariate", "time_dependence", "covariate"),
    n_subject = c(80L, 80L, 80L, 80L, 80L),
    n_time = c(3L, 3L, 3L, 4L, 4L),
    dependence = c("moderate", "moderate", "moderate", "moderate", "moderate"),
    missingness = c("none", "none", "none", "none", "drop_rows"),
    claim = c(
      "Tests whether a true scale-varying GAMLSS margin improves calibration over mean-only longitudinal baselines.",
      "Tests whether a positive GAMLSS margin improves mean, quantile, and tail behaviour versus standard mean baselines.",
      "Tests whether count margins and copula dependence improve dispersion and upper-tail behaviour.",
      "Tests whether theta formulas recover changing adjacent-time dependence that standard exchangeable baselines cannot represent.",
      "Tests whether the workflow remains stable when common follow-up visits are absent."
    ),
    primary_metrics = I(list(
      c("benchmark_mean_rmse", "benchmark_neg_log_score", "benchmark_q90_mae", "benchmark_upper_tail_error_90", "benchmark_interval_coverage_95", "benchmark_interval_width_95", "elapsed_sec"),
      c("benchmark_mean_rmse", "benchmark_neg_log_score", "benchmark_q90_mae", "benchmark_upper_tail_error_90", "elapsed_sec"),
      c("benchmark_mean_rmse", "benchmark_neg_log_score", "benchmark_q90_mae", "benchmark_upper_tail_error_90", "elapsed_sec"),
      c("benchmark_theta_time_abs_error", "elapsed_sec"),
      c("benchmark_mean_rmse", "benchmark_neg_log_score", "benchmark_interval_coverage_95", "benchmark_interval_width_95", "elapsed_sec")
    )),
    methods = I(list(
      c("rs_separate", "gee", "glmm", "gam"),
      c("rs_separate", "gee", "glmm", "gam"),
      c("rs_separate", "gee", "glmm", "gam"),
      c("rs_separate", "gee", "glmm", "gam"),
      c("rs_separate", "gee", "glmm", "gam")
    )),
    stringsAsFactors = FALSE
  )


  if (!is.null(scenarios)) {
    scenarios <- as.character(scenarios)

    missing <- setdiff(scenarios, out$scenario)

    if (length(missing) > 0L) {
      stop("Unknown adoption benchmark scenario(s): ", paste(missing, collapse = ", "), call. = FALSE)
    }

    out <- out[match(scenarios, out$scenario), , drop = FALSE]
  }

  rownames(out) <- NULL

  out
}
