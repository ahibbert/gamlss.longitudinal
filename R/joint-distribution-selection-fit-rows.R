.joint_selection_success_row <- function(
    margin_family,
    copula_family,
    fit,
    fit_metrics,
    elapsed,
    warnings) {
  data.frame(
    margin_family = margin_family,
    copula_family = copula_family,
    logLik = as.numeric(fit_metrics$logLik),
    AIC = as.numeric(fit_metrics$AIC),
    BIC = as.numeric(fit_metrics$BIC),
    EDF = as.numeric(fit_metrics$model_selection["EDF", "joint"]),
    converged = isTRUE(fit$convergence$converged),
    hit_outer_limit = isTRUE(fit$convergence$hit_outer_limit),
    elapsed_sec = elapsed,
    warnings = paste(unique(warnings), collapse = "\n"),
    error = NA_character_,
    stringsAsFactors = FALSE
  )
}

.joint_selection_failed_row <- function(margin_family, copula_family, elapsed, warnings, error) {
  data.frame(
    margin_family = margin_family,
    copula_family = copula_family,
    logLik = NA_real_,
    AIC = NA_real_,
    BIC = NA_real_,
    EDF = NA_real_,
    converged = FALSE,
    hit_outer_limit = NA,
    elapsed_sec = elapsed,
    warnings = paste(unique(warnings), collapse = "\n"),
    error = error,
    stringsAsFactors = FALSE
  )
}
