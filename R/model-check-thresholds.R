.gl_check_table <- function(summary_obj, scores, pit_stats, tail_stats, lag1_cor,
                            dependence_cor_cutoff = NULL, vcov_method = NA_character_) {
  converged <- isTRUE(summary_obj$convergence$converged)
  has_user_cutoff <- !is.null(dependence_cor_cutoff)
  dependence_flagged <- has_user_cutoff && is.finite(lag1_cor) &&
    abs(lag1_cor) > dependence_cor_cutoff
  variance_review <- identical(vcov_method, "numderiv")

  do.call(rbind, list(
    .gl_check_row(
      area = "Convergence",
      quantity_checked = "object$convergence$converged",
      value = if (converged) "TRUE" else "not TRUE",
      threshold_condition = "Fitted-object convergence contract",
      default = "n/a",
      status = if (converged) "converged" else "not_converged",
      message = if (converged) "Convergence was confirmed." else "Convergence was not confirmed.",
      action = if (converged) {
        "Continue with graphical and subject-matter diagnostics."
      } else {
        "Refit with more iterations, different starts, or a simpler specification before inference."
      }
    ),
    .gl_check_row(
      area = "Marginal fit",
      quantity_checked = "PIT Kolmogorov-Smirnov p-value vs Uniform(0, 1)",
      value = if (is.finite(pit_stats$ks_p_value)) formatC(pit_stats$ks_p_value, digits = 4, format = "fg") else NA_character_,
      threshold_condition = "Descriptive only; no package decision threshold",
      default = "none",
      status = "descriptive",
      message = "The PIT uniformity statistic is reported as a descriptive diagnostic, not a model verdict.",
      action = "Interpret with PIT, QQ, worm, and rootogram plots and the analysis context."
    ),
    .gl_check_row(
      area = "Tail fit",
      quantity_checked = "Maximum lower/upper PIT tail ratio over thresholds 0.05 and 0.10",
      value = if (is.finite(tail_stats$tail_ratio_max)) formatC(tail_stats$tail_ratio_max, digits = 4, format = "fg") else NA_character_,
      threshold_condition = "Descriptive only; no package decision threshold",
      default = "none",
      status = "descriptive",
      message = "Tail ratios are reported descriptively and are not converted into pass/fail decisions.",
      action = "Inspect lower and upper tails and compare substantively plausible alternative margins."
    ),
    .gl_check_row(
      area = "Copula fit",
      quantity_checked = "Absolute lag-1 Rosenblatt normal-score residual correlation after fitted copula",
      value = if (is.finite(lag1_cor)) formatC(abs(lag1_cor), digits = 4, format = "fg") else NA_character_,
      threshold_condition = if (has_user_cutoff) {
        paste0("User-supplied flag: abs(lag1_cor) > ", dependence_cor_cutoff)
      } else {
        "Descriptive only; no package decision threshold"
      },
      default = "none",
      status = if (!has_user_cutoff) "descriptive" else if (dependence_flagged) "flagged" else "not_flagged",
      message = if (!has_user_cutoff) {
        "Residual dependence is reported descriptively; higher-lag availability is recorded separately."
      } else if (dependence_flagged) {
        "The lag-1 residual correlation exceeded the user-supplied review threshold."
      } else {
        "The lag-1 residual correlation did not exceed the user-supplied review threshold."
      },
      action = "Review lag summaries and graphical copula diagnostics; threshold flags are not formal tests."
    ),
    .gl_check_row(
      area = "Variance calculation",
      quantity_checked = "Variance-covariance method from summary",
      value = if (is.na(vcov_method) || !nzchar(vcov_method)) NA_character_ else vcov_method,
      threshold_condition = "Method provenance only",
      default = "n/a",
      status = if (variance_review) "review" else if (is.na(vcov_method) || !nzchar(vcov_method)) "not_computed" else "recorded",
      message = if (variance_review) {
        "Variance-covariance inference used the numerical Hessian path."
      } else if (is.na(vcov_method) || !nzchar(vcov_method)) {
        "Variance-covariance inference was not requested for this diagnostic summary."
      } else {
        "The variance-covariance method is recorded for provenance."
      },
      action = if (variance_review) {
        "Report intervals and tests as approximate numerical-Hessian inference."
      } else {
        "Use the recorded inference contract when reporting intervals or tests."
      }
    )
  ))
}
