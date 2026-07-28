.gl_check_table <- function(summary_obj, scores, pit_stats, tail_stats, lag1_cor,
                            dependence_cor_cutoff = 0.25, vcov_method = NA_character_) {
  converged <- isTRUE(summary_obj$convergence$converged)

  marginal_fail <- is.finite(pit_stats$ks_p_value) && pit_stats$ks_p_value < 0.05

  tail_fail <- is.finite(tail_stats$tail_ratio_max) && tail_stats$tail_ratio_max > 2

  copula_fail <- is.finite(lag1_cor) && abs(lag1_cor) > dependence_cor_cutoff

  variance_review <- identical(vcov_method, "numderiv")

  do.call(rbind, list(
    .gl_check_row(
      area = "Convergence",
      quantity_checked = "object$convergence$converged",
      value = if (converged) "TRUE" else "not TRUE",
      threshold_condition = "Not TRUE",
      default = "n/a",
      status = if (converged) "PASS" else "FAIL",
      message = if (converged) "Convergence was confirmed." else "Convergence was not confirmed.",
      action = if (converged) {
        "Continue with broader diagnostics."
      } else {
        "Refit with more iterations, different starts, or a simpler specification."
      }
    ),
    .gl_check_row(
      area = "Marginal fit",
      quantity_checked = "PIT Kolmogorov-Smirnov p-value vs Uniform(0, 1)",
      value = if (is.finite(pit_stats$ks_p_value)) formatC(pit_stats$ks_p_value, digits = 4, format = "fg") else NA_character_,
      threshold_condition = "ks_p_value < 0.05",
      default = "0.05",
      status = if (marginal_fail) "FAIL" else "PASS",
      message = if (marginal_fail) {
        "The marginal distribution is off by the PIT uniformity screen."
      } else {
        "The PIT uniformity screen did not flag marginal misfit."
      },
      action = if (marginal_fail) {
        "Inspect PIT, QQ, worm, and rootogram diagnostics; try a richer margin or covariate specification."
      } else {
        "Continue with visual marginal diagnostics."
      }
    ),
    .gl_check_row(
      area = "Tail fit",
      quantity_checked = "Maximum lower/upper PIT tail ratio over thresholds 0.05 and 0.10",
      value = if (is.finite(tail_stats$tail_ratio_max)) formatC(tail_stats$tail_ratio_max, digits = 4, format = "fg") else NA_character_,
      threshold_condition = "max(lower_ratio, upper_ratio) > 2",
      default = "2",
      status = if (tail_fail) "FAIL" else "PASS",
      message = if (tail_fail) {
        "Tail observations occur more often than the fitted margin expects."
      } else {
        "The basic PIT tail-ratio screen did not flag tail misfit."
      },
      action = if (tail_fail) {
        "Inspect lower/upper PIT tails and consider heavier-tailed or asymmetric margins."
      } else {
        "Continue with tail-sensitive diagnostics when tails are substantively important."
      }
    ),
    .gl_check_row(
      area = "Copula fit",
      quantity_checked = "Absolute lag-1 Rosenblatt normal-score residual correlation after fitted copula",
      value = if (is.finite(lag1_cor)) formatC(abs(lag1_cor), digits = 4, format = "fg") else NA_character_,
      threshold_condition = "abs(lag1_cor) > dependence_cor_cutoff",
      default = formatC(dependence_cor_cutoff, digits = 4, format = "fg"),
      status = if (copula_fail) "FAIL" else "PASS",
      message = if (copula_fail) {
        paste0(
          "Dependence remains after the copula in Rosenblatt normal-score residuals (|lag-1 cor| > ",
          dependence_cor_cutoff,
          ")."
        )
      } else {
        "The lag-1 Rosenblatt residual correlation screen did not flag residual dependence."
      },
      action = if (copula_fail) {
        "Consider a different copula family, time-varying dependence, richer serial structure, or a sensitivity refit before treating this as a failure."
      } else {
        "Continue with broader copula diagnostics."
      }
    ),
    .gl_check_row(
      area = "Variance calculation",
      quantity_checked = "Variance-covariance method from summary",
      value = if (is.na(vcov_method) || !nzchar(vcov_method)) NA_character_ else vcov_method,
      threshold_condition = 'vcov_method == "numderiv"',
      default = "n/a",
      status = if (variance_review) "REVIEW" else "PASS",
      message = if (variance_review) {
        "Variance-covariance inference used the numerical Hessian path."
      } else {
        "The variance-covariance method did not trigger the numerical-Hessian review screen."
      },
      action = if (variance_review) {
        "Cite intervals and tests as approximate numerical-Hessian inference."
      } else {
        "Continue with inference checks when reporting intervals or tests."
      }
    )
  ))
}
