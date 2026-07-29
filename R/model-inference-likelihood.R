.gl_model_edf <- function(object) {
  fixed <- length(object$par)

  smooth <- 0

  if (!is.null(object$df_s) && length(object$df_s) > 0L) {
    vals <- suppressWarnings(as.numeric(unlist(object$df_s, use.names = FALSE)))

    vals <- vals[is.finite(vals)]

    smooth <- sum(vals)
  }

  fixed + smooth
}

.gl_joint_loglik <- function(object) {
  ll <- NULL

  if (!is.null(object$calc_lik_out_end) && !is.null(object$calc_lik_out_end$log_lik)) {
    ll <- object$calc_lik_out_end$log_lik
  } else if (!is.null(object$calc_lik_out) && !is.null(object$calc_lik_out$log_lik)) {
    ll <- object$calc_lik_out$log_lik
  }

  if (!is.null(ll) && "joint" %in% names(ll)) {
    return(as.numeric(ll[["joint"]]))
  }

  NA_real_
}

#' Compare fitted models with likelihood-ratio summaries

#'

#' `likelihood_compare()` gives a compact sequential likelihood comparison for

#' nested or approximately nested `gamlss.longitudinal` models. It reports joint

#' log-likelihood, effective degrees of freedom, AIC, BIC, likelihood-ratio

#' increments, and chi-square reference p-values.

#'

#' @param ... Fitted `gamlss.longitudinal` objects, or a single list of fitted

#'   objects.

#' @param sort Logical; order models by AIC, lowest first, before

#'   computing sequential comparisons.

#'

#' @return An object of class `gamlss_longitudinal_likelihood_compare`.

#' @export

likelihood_compare <- function(..., sort = TRUE) {
  labels <- .gl_likelihood_compare_call_labels(substitute(list(...)))

  models <- .gl_likelihood_compare_models(list(...), labels = labels)

  .gl_likelihood_compare_table(models, sort = sort)
}

#' @export

print.gamlss_longitudinal_likelihood_compare <- function(x, digits = max(3, getOption("digits") - 3), ...) {
  cat("\nLikelihood Comparison for gamlss.longitudinal\n")

  cat("--------------------------------------------\n")

  cat("Sequential LR rows compare each model with the previous row.\n\n")

  print.data.frame(x, digits = digits, row.names = FALSE, ...)

  invisible(x)
}
