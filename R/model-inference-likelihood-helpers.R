#' Normalize and validate likelihood-comparison model inputs
#'
#' @param models List collected from `...`.
#' @return Named list of fitted models.
#' @noRd
.gl_likelihood_compare_models <- function(models) {
  if (length(models) == 1L && is.list(models[[1L]]) && !inherits(models[[1L]], "gamlss.longitudinal")) {
    models <- models[[1L]]
  }

  if (length(models) < 2L) {
    stop("At least two fitted models are required.", call. = FALSE)
  }

  ok <- vapply(models, inherits, logical(1), what = "gamlss.longitudinal")

  if (!all(ok)) {
    stop("All inputs must be fitted 'gamlss.longitudinal' objects.", call. = FALSE)
  }

  labels <- names(models)

  if (is.null(labels) || any(labels == "")) {
    labels <- paste0("model_", seq_along(models))
  }

  names(models) <- labels
  models
}

#' Build the likelihood-comparison result table
#'
#' @param models Named list of fitted `gamlss.longitudinal` objects.
#' @param sort Logical; order rows by effective degrees of freedom and log-likelihood.
#' @return Classed likelihood-comparison data frame.
#' @noRd
.gl_likelihood_compare_table <- function(models, sort = TRUE) {
  labels <- names(models)

  n_obs <- vapply(models, function(x) length(x$response), integer(1))

  if (length(unique(n_obs)) > 1L) {
    warning(
      "Models have different observation counts; likelihood-ratio comparisons may not be valid.",
      call. = FALSE
    )
  }

  df <- vapply(models, .gl_model_edf, numeric(1))

  loglik <- vapply(models, .gl_joint_loglik, numeric(1))

  if (isTRUE(sort)) {
    ord <- order(df, loglik)

    labels <- labels[ord]

    n_obs <- n_obs[ord]

    df <- df[ord]

    loglik <- loglik[ord]
  }

  aic <- -2 * loglik + 2 * df

  bic <- -2 * loglik + log(pmax(1, n_obs)) * df

  delta_df <- c(NA_real_, diff(df))

  lr <- c(NA_real_, 2 * diff(loglik))

  p_value <- rep(NA_real_, length(models))

  valid <- is.finite(lr) & is.finite(delta_df) & delta_df > 0

  p_value[valid] <- stats::pchisq(lr[valid], df = delta_df[valid], lower.tail = FALSE)

  out <- data.frame(
    model = labels,
    n_obs = as.integer(n_obs),
    df = as.numeric(df),
    logLik = as.numeric(loglik),
    AIC = as.numeric(aic),
    BIC = as.numeric(bic),
    delta_df = as.numeric(delta_df),
    LR_statistic = as.numeric(lr),
    p_value = as.numeric(p_value),
    stringsAsFactors = FALSE
  )

  class(out) <- c("gamlss_longitudinal_likelihood_compare", "data.frame")

  out
}
