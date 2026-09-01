#' Normalize and validate likelihood-comparison model inputs
#'
#' @param models List collected from `...`.
#' @param labels Optional labels captured from the unevaluated call.
#' @return Named list of fitted models.
#' @noRd
.gl_likelihood_compare_models <- function(models, labels = NULL) {
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

  model_names <- names(models)

  if (is.null(model_names)) {
    model_names <- rep("", length(models))
  }

  if (!is.null(labels) && length(labels) == length(models)) {
    missing_names <- model_names == ""

    model_names[missing_names] <- labels[missing_names]
  }

  if (any(model_names == "")) {
    model_names[model_names == ""] <- paste0("model_", which(model_names == ""))
  }

  names(models) <- model_names
  models
}

#' Capture likelihood-comparison labels from the user call
#'
#' @param expr Unevaluated `list(...)` call from `likelihood_compare()`.
#' @return Character vector of labels, with empty strings for non-symbol inputs.
#' @noRd
.gl_likelihood_compare_call_labels <- function(expr) {
  args <- as.list(expr)[-1L]

  if (length(args) == 1L && is.call(args[[1L]]) && identical(args[[1L]][[1L]], as.name("list"))) {
    args <- as.list(args[[1L]])[-1L]
  }

  labels <- rep("", length(args))

  arg_names <- names(args)

  if (is.null(arg_names)) {
    arg_names <- rep("", length(args))
  }

  for (i in seq_along(args)) {
    if (nzchar(arg_names[[i]])) {
      labels[[i]] <- arg_names[[i]]
    } else if (is.symbol(args[[i]])) {
      labels[[i]] <- as.character(args[[i]])
    }
  }

  labels
}

#' Build the likelihood-comparison result table
#'
#' @param models Named list of fitted `gamlss.longitudinal` objects.
#' @param sort Logical; order rows by AIC, lowest first.
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

  aic <- -2 * loglik + 2 * df

  bic <- -2 * loglik + log(pmax(1, n_obs)) * df

  if (isTRUE(sort)) {
    ord <- order(aic, df, -loglik)

    labels <- labels[ord]

    n_obs <- n_obs[ord]

    df <- df[ord]

    loglik <- loglik[ord]

    aic <- aic[ord]

    bic <- bic[ord]
  }

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

  attr(out, "inference_contract") <- .gl_inference_contract(
    "likelihood_ratio_nested",
    validity_status = if (all(valid[-1L])) "reference_available" else "partly_or_wholly_unavailable"
  )

  out
}
