.copula_v2_gaussian_limit_warning <- function(object, threshold = 7.5) {
  if (is.null(object) || is.null(object$model$copula_dist)) {
    return(NULL)
  }

  copula_dist <- object$model$copula_dist
  if (!identical(copula_dist, "t")) {
    return(NULL)
  }

  coef_tbl <- object$coefficients
  if (is.null(coef_tbl) || nrow(coef_tbl) == 0 || !("parameter" %in% names(coef_tbl))) {
    return(NULL)
  }

  zeta_rows <- coef_tbl[coef_tbl$parameter == "zeta" & is.finite(coef_tbl$estimate), , drop = FALSE]
  if (nrow(zeta_rows) == 0) {
    return(NULL)
  }

  zeta_link <- zeta_rows$estimate
  near_limit <- is.finite(zeta_link) & zeta_link >= threshold
  if (!any(near_limit)) {
    return(NULL)
  }

  zeta_nat <- exp(zeta_link[near_limit]) + 2
  zeta_label <- paste0(formatC(zeta_link[near_limit], format = "f", digits = 2), collapse = ", ")
  df_label <- paste0(formatC(zeta_nat, format = "f", digits = 1), collapse = ", ")

  paste0(
    "WARNING: t-copula zeta is near the Gaussian-limit regime (link-scale zeta = ",
    zeta_label,
    "; implied degrees of freedom ",
    df_label,
    "). The t-copula may be collapsing toward Gaussian dependence."
  )
}

#' @export
print.summary.gamlss.longitudinal <- function(x, digits = max(3, getOption("digits") - 3), ...) {
  cat("\nGAMLSS Longitudinal Model Summary\n")
  cat("--------------------------------\n")
  cat("Margin distribution:", x$model$margin_dist, "\n")
  cat("Copula distribution:", x$model$copula_dist, "\n")
  if (identical(x$convergence$converged, FALSE)) {
    cat(
      "Convergence: NOT CONVERGED (stop reason:",
      x$convergence$stop_reason %||% "unknown", ")\n"
    )
    cat("Inference is disabled; likelihood criteria below are provisional and must not be used for model selection.\n")
  } else {
    cat("Convergence: CONVERGED\n")
  }

  copula_warning <- .copula_v2_gaussian_limit_warning(x)
  if (!is.null(copula_warning)) {
    cat(copula_warning, "\n")
  }

  cat(
    "Observations:", x$model$n_obs,
    " | Subjects:", x$model$n_subjects,
    " | Time points:", x$model$n_timepoints, "\n"
  )
  if (identical(x$fit$likelihood_contract$objective, "segmented")) {
    cat(
      "Missingness: SEGMENTED LIKELIHOOD (",
      x$fit$missingness$n_subjects_with_gaps %||% NA_integer_,
      " subject(s) with intermittent gaps; ",
      x$fit$missingness$n_segments %||% NA_integer_,
      " contiguous segments)\n",
      sep = ""
    )
    cat("Between-gap assumption: different observed segments are treated as independent.\n")
    cat("AIC and BIC are available under this explicit segmented-independence assumption.\n")
    cat("Numerical integration of dependence across gaps is not currently implemented.\n")
  }
  cat(
    "Fixed coefficients:", x$model$n_fixed,
    " | Smooth terms:", x$model$n_smooth_terms,
    " | Smooth EDF:", format(round(x$model$edf_smooth, digits), nsmall = 2), "\n"
  )
  if (isTRUE(x$fit$vcov_included)) {
    cat(
      "VCOV:", x$fit$vcov_method %||% "unknown",
      " | Requested:", x$fit$vcov_method_requested %||% "unknown",
      "\n"
    )
    hd <- x$fit$hessian_diagnostics
    if (!is.null(hd) && is.finite(hd$scaled_condition_number %||% NA_real_)) {
      cat(
        "Scaled information condition number:",
        formatC(hd$scaled_condition_number, digits = digits, format = "fg"),
        " | Validation:", hd$validation_profile %||% "unknown", "\n"
      )
    }
    contract <- x$fit$inference_contract %||% attr(x, "inference_contract")
    if (!is.null(contract)) {
      cat("Inference target: fixed coefficients conditional on fitted smooth structure.\n")
      cat("Omitted uncertainty: fixed-smooth covariance and smoothing-parameter uncertainty.\n")
    }
  }

  cat("\nFixed coefficients:\n")
  cat("--------------------\n")
  coef_disp <- .gl_summary_coefficient_display(x$coefficients, digits = digits)
  .gl_print_summary_coefficient_blocks(coef_disp)

  cat("  Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1\n")

  .gl_print_summary_smooth_terms(x$smooth_terms, digits = digits)
  .gl_print_summary_model_selection(x$fit, digits = digits)
  cat("--------------------------------\n")

  invisible(x)
}
