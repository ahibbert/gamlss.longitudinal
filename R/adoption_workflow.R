#' Extract the best-fitting candidate from a selection result
#'
#' @param x A selection result, such as from [select_margin()] or
#'   [select_copula()].
#' @param ... Reserved for methods.
#'
#' @return `best_fit()` returns a list of selected-fit metadata.
#'   `best_fit_family()` returns the family value that can be supplied directly
#'   to fitting helpers.
#' @export
best_fit <- function(x, ...) {
  UseMethod("best_fit")
}

#' @rdname best_fit
#' @export
best_fit_family <- function(x, ...) {
  UseMethod("best_fit_family")
}

.margin_family_object <- function(family_name) {
  if (length(family_name) != 1L || is.na(family_name) || !nzchar(family_name)) {
    return(NULL)
  }
  family_fun <- tryCatch(
    get(family_name, envir = asNamespace("gamlss.dist"), mode = "function", inherits = FALSE),
    error = function(e) NULL
  )
  if (is.null(family_fun)) {
    return(NULL)
  }
  tryCatch(do.call(family_fun, list()), error = function(e) NULL)
}

#' @export
best_fit.margin_selection <- function(x, ...) {
  if (nrow(x) == 0L) {
    return(list(
      family_name = NA_character_,
      family = NULL,
      rank = NA_integer_,
      AIC = NA_real_,
      type = attr(x, "response_type")
    ))
  }
  row <- as.data.frame(x)[1L, , drop = FALSE]
  family_name <- as.character(row$family[[1L]])
  row_meta <- as.list(row)
  row_meta$family <- NULL
  c(
    list(
      family_name = family_name,
      family = .margin_family_object(family_name)
    ),
    row_meta
  )
}

#' @export
best_fit.margin_screen <- best_fit.margin_selection

#' @export
best_fit_family.margin_selection <- function(x, ...) {
  best_fit(x)$family
}

#' @export
best_fit_family.margin_screen <- best_fit_family.margin_selection

#' @export
`$.margin_selection` <- function(x, name) {
  if (identical(name, "best_fit")) {
    return(best_fit(x))
  }
  .subset2(as.data.frame(x), name, exact = FALSE)
}

#' @export
`$.margin_screen` <- `$.margin_selection`

#' Select candidate marginal distributions
#'
#' `select_margin()` is a lightweight wrapper around [gamlss::fitDist()] for
#' the recommended longitudinal workflow: choose a plausible marginal family,
#' then fit dependence with [gamlss_longitudinal()].
#'
#' @param response Numeric response vector. For the common
#'   `select_margin(dat, response_var = "y")` call, a data frame supplied here
#'   is treated as `data`.
#' @param data Optional data frame containing the response.
#' @param response_var Optional response column name in `data`.
#' @param type Optional `gamlss::fitDist()` type. If `NULL`, a simple heuristic
#'   uses `"counts"` for non-negative integer responses, `"realplus"` for
#'   positive continuous responses, and `"realAll"` otherwise.
#' @param families Optional character vector used to filter the returned table.
#' @param time_intercepts Logical; if `TRUE`, rank retained families by a
#'   GAMLSS fit with factor time intercepts for each distribution parameter.
#'   This is useful for longitudinal screening where the final model is allowed
#'   to vary marginal location, scale, or shape over time.
#' @param time_var Optional time column name in `data`, required when
#'   `time_intercepts = TRUE`.
#' @param try.gamlss,trace,... Passed to [gamlss::fitDist()].
#'
#' @return A data frame ordered by AIC and class `margin_selection`. The
#'   selected family is also stored in the `"selected"` attribute for backward
#'   compatibility. Use [best_fit()] or [best_fit_family()] to extract the
#'   selected family in a fitting-friendly form.
#' @export
select_margin <- function(
  response = NULL,
  data = NULL,
  response_var = NULL,
  type = NULL,
  families = NULL,
  time_intercepts = FALSE,
  time_var = NULL,
  try.gamlss = FALSE,
  trace = FALSE,
  ...
) {
  if (is.null(data) && is.data.frame(response) && !is.null(response_var)) {
    data <- response
    response <- NULL
  }

  if (!is.null(data)) {
    data <- as.data.frame(data, stringsAsFactors = FALSE)
    if (is.null(response_var) || !is.character(response_var) || length(response_var) != 1L) {
      stop("'response_var' must be a single column name when 'data' is supplied.", call. = FALSE)
    }
    if (!response_var %in% names(data)) {
      stop("response_var='", response_var, "' not found in 'data'.", call. = FALSE)
    }
    response <- data[[response_var]]
  }

  if (is.null(response)) {
    stop("Supply either 'response' or both 'data' and 'response_var'.", call. = FALSE)
  }
  response_all <- as.numeric(response)
  if (isTRUE(time_intercepts)) {
    if (is.null(data)) {
      stop("'time_intercepts = TRUE' requires 'data' and 'time_var'.", call. = FALSE)
    }
    if (is.null(time_var) || !is.character(time_var) || length(time_var) != 1L) {
      stop("'time_var' must be a single column name when 'time_intercepts = TRUE'.", call. = FALSE)
    }
    if (!time_var %in% names(data)) {
      stop("time_var='", time_var, "' not found in 'data'.", call. = FALSE)
    }
  }
  response <- response_all[is.finite(response_all)]
  if (length(response) < 3L) {
    stop("Need at least three finite response values to screen margins.", call. = FALSE)
  }

  if (is.null(type)) {
    is_count <- all(response >= 0) && all(abs(response - round(response)) < .Machine$double.eps^0.5)
    type <- if (is_count) {
      "counts"
    } else if (all(response > 0)) {
      "realplus"
    } else {
      "realAll"
    }
  }

  fit <- NULL
  if (isTRUE(trace)) {
    fit <- gamlss::fitDist(
      response,
      type = type,
      try.gamlss = try.gamlss,
      trace = trace,
      ...
    )
  } else {
    invisible(utils::capture.output({
      fit <- suppressWarnings(suppressMessages(gamlss::fitDist(
        response,
        type = type,
        try.gamlss = try.gamlss,
        trace = trace,
        ...
      )))
    }))
  }

  if (is.null(fit$fits)) {
    stop("gamlss::fitDist() did not return a fits table.", call. = FALSE)
  }

  out <- data.frame(
    family = names(fit$fits),
    AIC = as.numeric(fit$fits),
    type = type,
    stringsAsFactors = FALSE,
    row.names = NULL
  )
  out <- out[is.finite(out$AIC), , drop = FALSE]
  if (!is.null(families)) {
    out <- out[out$family %in% families, , drop = FALSE]
  }

  out$screen_model <- "pooled_intercept"

  if (isTRUE(time_intercepts) && nrow(out) > 0L) {
    time_aic <- .select_margin_time_intercept_aic(
      response = response_all,
      time = data[[time_var]],
      families = out$family,
      trace = trace
    )
    out$pooled_AIC <- out$AIC
    out$AIC <- time_aic[out$family]
    out$screen_model <- "time_intercepts"
    out <- out[is.finite(out$AIC), , drop = FALSE]
  }

  out <- out[order(out$AIC), , drop = FALSE]
  rownames(out) <- NULL
  out$rank <- seq_len(nrow(out))
  out$delta_AIC <- if (nrow(out) > 0L) out$AIC - min(out$AIC, na.rm = TRUE) else numeric(0)
  out$supported_by_longitudinal <- vapply(out$family, function(family) {
    all(vapply(
      paste0(c("d", "p", "q"), family),
      exists,
      logical(1),
      envir = asNamespace("gamlss.dist"),
      inherits = FALSE
    ))
  }, logical(1), USE.NAMES = FALSE)

  attr(out, "selected") <- if (nrow(out) > 0L) out$family[[1L]] else NA_character_
  attr(out, "response_type") <- type
  attr(out, "time_intercepts") <- isTRUE(time_intercepts)
  attr(out, "time_var") <- if (isTRUE(time_intercepts)) time_var else NULL
  class(out) <- c("margin_selection", "margin_screen", class(out))
  out
}

.select_margin_time_intercept_aic <- function(response, time, families, trace = FALSE) {
  response <- as.numeric(response)
  time <- as.character(time)
  keep <- is.finite(response) & !is.na(time)
  if (sum(keep) < 3L) {
    stop("Need at least three finite response values with non-missing time values to screen time-intercept margins.", call. = FALSE)
  }

  fit_data <- data.frame(
    response = response[keep],
    time_intercept = factor(time[keep], levels = unique(time[keep])),
    stringsAsFactors = FALSE
  )
  has_time_contrast <- length(unique(fit_data$time_intercept)) > 1L
  mu_formula <- if (has_time_contrast) {
    stats::as.formula("response ~ time_intercept")
  } else {
    stats::as.formula("response ~ 1")
  }
  par_formula <- if (has_time_contrast) {
    stats::as.formula("~ time_intercept")
  } else {
    stats::as.formula("~ 1")
  }

  out <- vapply(families, function(family_name) {
    family <- .margin_family_object(family_name)
    if (is.null(family) || is.null(family$parameters)) {
      return(NA_real_)
    }

    fit_args <- list(
      formula = mu_formula,
      family = family,
      data = fit_data,
      trace = trace
    )
    for (par_name in setdiff(names(family$parameters), "mu")) {
      fit_args[[paste0(par_name, ".formula")]] <- par_formula
    }

    fit <- tryCatch(
      {
        if (isTRUE(trace)) {
          do.call(gamlss::gamlss, fit_args)
        } else {
          fit <- NULL
          invisible(utils::capture.output({
            fit <- suppressWarnings(suppressMessages(do.call(gamlss::gamlss, fit_args)))
          }))
          fit
        }
      },
      error = function(e) NULL
    )
    if (is.null(fit)) {
      return(NA_real_)
    }
    if (identical(fit$converged, FALSE)) {
      return(NA_real_)
    }
    aic <- tryCatch(stats::AIC(fit), error = function(e) NA_real_)
    as.numeric(aic)[1L]
  }, numeric(1), USE.NAMES = TRUE)

  out
}

#' @rdname select_margin
#' @export
screen_margin <- function(...) {
  select_margin(...)
}

#' @export
print.margin_screen <- function(x, ..., n = 10L) {
  cat("\nMarginal Distribution Screen\n")
  cat("----------------------------\n")
  if (nrow(x) == 0L) {
    cat("No candidate families were retained.\n")
    return(invisible(x))
  }
  cat("Selected:", attr(x, "selected"), "\n\n")
  display <- utils::head(as.data.frame(x), n = n)
  unsupported <- "supported_by_longitudinal" %in% names(display) &
    any(!display$supported_by_longitudinal, na.rm = TRUE)
  display$supported_by_longitudinal <- NULL
  print(display, row.names = FALSE)
  if (isTRUE(unsupported)) {
    cat("\nNote: one or more printed families are not currently supported by gamlss_longitudinal().\n")
  }
  invisible(x)
}

.gl_prediction_frame <- function(object, newdata = NULL, require_response = FALSE) {
  diag_data <- .gl_fitted_distribution(object, newdata = newdata, require_response = require_response)
  out <- data.frame(
    subject = diag_data$subject,
    time = diag_data$time,
    response = diag_data$response,
    stringsAsFactors = FALSE
  )
  for (par_name in names(diag_data$params)) {
    out[[par_name]] <- as.numeric(diag_data$params[[par_name]])
  }
  out
}

.gl_quantile_columns <- function(params, family, probs) {
  out <- list()
  for (prob in probs) {
    label <- paste0("q", gsub("[^0-9]+", "", format(prob, trim = TRUE, scientific = FALSE)))
    if (identical(label, "q")) {
      label <- paste0("q", seq_along(out) + 1L)
    }
    out[[label]] <- .gl_call_family_fun("q", family, prob, params)
  }
  as.data.frame(out, stringsAsFactors = FALSE)
}

.gl_mu_values <- function(params) {
  if ("mu" %in% names(params)) {
    return(as.numeric(params$mu))
  }
  as.numeric(params[[1L]])
}

.gl_distribution_mean <- function(params, family, probs = seq(0.001, 0.999, length.out = 199L)) {
  family <- as.character(family)[1L]
  if (identical(family, "NO") && "mu" %in% names(params)) {
    return(as.numeric(params$mu))
  }
  if (identical(family, "GA") && "mu" %in% names(params)) {
    return(as.numeric(params$mu))
  }
  if (identical(family, "LOGNO") && all(c("mu", "sigma") %in% names(params))) {
    return(exp(as.numeric(params$mu) + 0.5 * as.numeric(params$sigma)^2))
  }

  probs <- as.numeric(probs)
  probs <- probs[is.finite(probs) & probs > 0 & probs < 1]
  if (length(probs) < 3L) {
    stop("'probs' must contain at least three finite probabilities inside (0, 1).", call. = FALSE)
  }

  q_mat <- vapply(probs, function(prob) {
    as.numeric(.gl_call_family_fun("q", family, prob, params))
  }, numeric(length(params[[1L]])))
  as.numeric(rowMeans(q_mat, na.rm = TRUE))
}

.gl_prediction_values <- function(type, params, family) {
  switch(
    type,
    response = .gl_mu_values(params),
    mu = .gl_mu_values(params),
    mean = .gl_distribution_mean(params, family),
    median = as.numeric(.gl_call_family_fun("q", family, 0.5, params)),
    stop("Unsupported prediction value type.", call. = FALSE)
  )
}

#' Predict from a longitudinal GAMLSS-copula model
#'
#' @param object A fitted `gamlss.longitudinal` object.
#' @param newdata Optional new data. If omitted, fitted rows are used.
#' @param type Prediction type: `"mean"` returns the fitted marginal response
#'   mean, `"median"` returns the fitted marginal median, `"mu"` returns the
#'   fitted GAMLSS `mu` parameter, `"response"` is retained as a compatibility
#'   alias for `"mu"`, `"parameters"` returns all fitted marginal distribution
#'   parameters, `"quantile"` returns fitted marginal quantiles,
#'   `"cdf"`/`"density"` evaluate the fitted marginal CDF or density, and
#'   `"probability"` returns probabilities below or above a threshold.
#' @param probs Quantile probabilities when `type = "quantile"`.
#' @param q Threshold value for `type = "cdf"` or `type = "probability"`.
#'   Defaults to observed responses for `type = "cdf"` when available.
#' @param y Evaluation value for `type = "density"`. Defaults to observed
#'   responses when available.
#' @param direction Probability direction for `type = "probability"`.
#' @param se.fit Logical; for `type = "response"`, `"mu"`, or `"mean"`, return
#'   approximate delta-method standard errors for the fitted `mu` linear
#'   predictor contribution. For non-`mu` response means this is a first-order
#'   approximation and should be treated as exploratory.
#' @param interval Interval type. `"confidence"` adds response-scale confidence
#'   limits when `type = "response"`, `"mu"`, or `"mean"`.
#' @param level Confidence level for `interval = "confidence"`.
#' @param vcov_method Variance-covariance method passed to [vcov.gamlss.longitudinal()].
#' @param ... Additional arguments reserved for future methods.
#'
#' @details
#' `predict()` returns marginal summaries from the fitted distribution at each
#' requested row. The fitted copula/dependence structure is not used to
#' condition a row's prediction on that subject's other observed responses.
#' Instead, dependence affects prediction indirectly through the coefficients
#' estimated by the joint copula likelihood, and through `se.fit`/confidence
#' intervals when the covariance matrix is computed from the joint model. Use
#' [simulate.gamlss.longitudinal()] for fitted-data trajectory simulation that
#' preserves the fitted copula dependence structure. Copula-preserving
#' simulation for `newdata` is not currently implemented.
#'
#' `type = "response"` is a soft-deprecated compatibility alias for `type =
#' "mu"` because GAMLSS `mu` is not the response mean for every family. New code
#' should use `type = "mean"` for response-mean estimands or `type = "mu"` for
#' the distribution parameter.
#'
#' @return A numeric vector for `type = "response"`, `"mu"`, `"mean"`, or
#'   `"median"` unless standard errors or intervals are requested; a data frame
#'   otherwise.
#' @export
predict.gamlss.longitudinal <- function(
  object,
  newdata = NULL,
  type = c("response", "mean", "mu", "median", "parameters", "quantile", "cdf", "density", "probability"),
  probs = c(0.025, 0.5, 0.975),
  q = NULL,
  y = NULL,
  direction = c("below", "above"),
  se.fit = FALSE,
  interval = c("none", "confidence"),
  level = 0.95,
  vcov_method = "analytical",
  ...
) {
  type <- match.arg(type)
  interval <- match.arg(interval)
  direction <- match.arg(direction)
  if (!inherits(object, "gamlss.longitudinal")) {
    stop("'object' must be a fitted 'gamlss.longitudinal' object.", call. = FALSE)
  }

  require_response <- (type == "cdf" && is.null(q)) || (type == "density" && is.null(y))
  diag_data <- .gl_fitted_distribution(object, newdata = newdata, require_response = require_response)
  params <- diag_data$params

  if (type %in% c("response", "mu", "mean", "median")) {
    fit_values <- .gl_prediction_values(type, params, diag_data$family)
    if (isTRUE(se.fit) || !identical(interval, "none")) {
      pred <- .gl_prediction_frame(object, newdata = newdata, require_response = FALSE)
      pred$fit <- fit_values
      se_values <- .gl_predict_response_se(object, newdata = newdata, method = vcov_method, ...)
      pred$se.fit <- se_values
      if (!identical(interval, "none")) {
        z <- stats::qnorm((1 + level) / 2)
        pred$conf.low <- pred$fit - z * pred$se.fit
        pred$conf.high <- pred$fit + z * pred$se.fit
      }
      return(pred[c("subject", "time", "response", "fit", "se.fit", intersect(c("conf.low", "conf.high"), names(pred)))])
    }
    return(fit_values)
  }

  pred <- .gl_prediction_frame(object, newdata = newdata, require_response = require_response)
  if (type == "parameters") {
    return(pred)
  }

  if (type == "quantile") {
    q_df <- .gl_quantile_columns(params, diag_data$family, probs)
    return(cbind(pred[c("subject", "time", "response")], q_df))
  }

  expand_eval <- function(value, label) {
    if (is.null(value)) {
      if (!require_response) {
        stop("'", label, "' is required for type = '", type, "' when no response is available.", call. = FALSE)
      }
      return(as.numeric(diag_data$response))
    }
    value <- as.numeric(value)
    n <- length(diag_data$subject)
    if (length(value) == 1L) {
      return(rep(value, n))
    }
    if (length(value) != n) {
      stop("'", label, "' must be length 1 or match the number of prediction rows.", call. = FALSE)
    }
    value
  }

  if (type == "cdf") {
    q_use <- expand_eval(q, "q")
    pred$q <- q_use
    pred$cdf <- .gl_call_family_fun("p", diag_data$family, q_use, params)
    return(pred[c("subject", "time", "response", "q", "cdf")])
  }

  if (type == "probability") {
    if (is.null(q)) {
      stop("'q' is required for type = 'probability'.", call. = FALSE)
    }
    q_use <- expand_eval(q, "q")
    cdf <- .gl_call_family_fun("p", diag_data$family, q_use, params)
    pred$q <- q_use
    pred$direction <- direction
    pred$probability <- if (direction == "below") cdf else 1 - cdf
    return(pred[c("subject", "time", "response", "q", "direction", "probability")])
  }

  y_use <- expand_eval(y, "y")
  pred$y <- y_use
  pred$density <- .gl_call_family_fun("d", diag_data$family, y_use, params)
  pred[c("subject", "time", "response", "y", "density")]
}

.gl_prediction_model_matrix <- function(object, newdata = NULL) {
  if (is.null(newdata)) {
    return(object$model_matrix)
  }
  copula_link <- get_copula_dist(object$copula_dist)$copula_link
  nd <- .gl_prepare_newdata_internal(object, newdata, require_response = FALSE)
  mm_use <- do.call(
    create_model_matrices,
    list(
      mu.formula = object$formulas_int$mu,
      sigma.formula = object$formulas_int$sigma,
      nu.formula = object$formulas_int$nu,
      tau.formula = object$formulas_int$tau,
      theta.formula = object$formulas_int$theta,
      zeta.formula = object$formulas_int$zeta,
      margin.family = object$margin_dist,
      copula.family = object$copula_dist,
      copula.link = copula_link,
      dataset = nd,
      quiet_gamlss2 = TRUE,
      preserve_factor_levels = TRUE
    )
  )
  .gl_align_model_matrix_columns(mm_use, object$model_matrix)
}

.gl_predict_response_se <- function(object, newdata = NULL, method = "analytical", ...) {
  mm_use <- .gl_prediction_model_matrix(object, newdata = newdata)
  if (is.null(mm_use$x$mu) || ncol(mm_use$x$mu) == 0L) {
    return(rep(NA_real_, length(predict(object, newdata = newdata, type = "response"))))
  }
  X <- as.matrix(mm_use$x$mu)
  beta_names <- colnames(X)
  beta_names <- ifelse(startsWith(beta_names, "mu."), beta_names, paste0("mu.", beta_names))

  vc <- .resolve_vcov(object, extra_args = list(method = method, ...))
  V <- vc$vcov$overall
  V_names <- colnames(V) %||% rownames(V)
  idx <- match(beta_names, V_names)
  if (any(is.na(idx))) {
    return(rep(NA_real_, nrow(X)))
  }
  V_mu <- V[idx, idx, drop = FALSE]
  se_eta <- sqrt(pmax(0, rowSums((X %*% V_mu) * X)))

  copula_link <- get_copula_dist(object$copula_dist)$copula_link
  eta_out <- calc_eta(
    par_cov = object$par,
    mm = mm_use,
    margin_dist = object$margin_dist,
    copula_link = copula_link,
    par_s = object$par_s
  )
  mu_dr <- eta_out$eta_dr$mu %||% rep(1, length(se_eta))
  as.numeric(abs(mu_dr[seq_along(se_eta)]) * se_eta)
}

#' Confidence intervals for fixed coefficients
#'
#' @param object A fitted `gamlss.longitudinal` object.
#' @param parm Optional coefficient names or numeric indices.
#' @param level Confidence level.
#' @param method Variance-covariance method passed to [vcov.gamlss.longitudinal()].
#' @param ... Additional arguments passed to [vcov.gamlss.longitudinal()].
#'
#' @return A matrix with lower and upper confidence limits.
#' @importFrom stats confint
#' @export
confint.gamlss.longitudinal <- function(
  object,
  parm = NULL,
  level = 0.95,
  method = "analytical",
  ...
) {
  if (!inherits(object, "gamlss.longitudinal")) {
    stop("'object' must be a fitted 'gamlss.longitudinal' object.", call. = FALSE)
  }
  estimates <- stats::coef(object)
  vc <- .resolve_vcov(object, extra_args = list(method = method, ...))
  se <- vc$se$overall
  if (is.null(names(se))) {
    names(se) <- names(estimates)
  }
  idx <- if (is.null(parm)) {
    seq_along(estimates)
  } else if (is.character(parm)) {
    match(parm, names(estimates))
  } else {
    parm
  }
  if (any(is.na(idx)) || any(idx < 1L) || any(idx > length(estimates))) {
    stop("'parm' contains unknown coefficient names or indices.", call. = FALSE)
  }
  z <- stats::qnorm((1 + level) / 2)
  est <- estimates[idx]
  se_use <- se[names(est)]
  out <- cbind(
    lower = as.numeric(est) - z * as.numeric(se_use),
    upper = as.numeric(est) + z * as.numeric(se_use)
  )
  colnames(out) <- c(
    paste0(round((1 - level) / 2 * 100, 1), " %"),
    paste0(round((1 + level) / 2 * 100, 1), " %")
  )
  rownames(out) <- names(est)
  out
}

.fixed_term_coefficient_names <- function(object) {
  if (!inherits(object, "gamlss.longitudinal") || is.null(object$model_matrix$x)) {
    return(NULL)
  }

  out <- list()
  for (parameter in names(object$model_matrix$x)) {
    X <- object$model_matrix$x[[parameter]]
    if (is.null(X)) next
    col_names <- colnames(X)
    if (length(col_names) == 0L) next

    assign <- attr(X, "assign")
    term_labels <- attr(X, "term.labels")
    if (length(assign) != length(col_names) || length(term_labels) == 0L) next

    for (term_idx in seq_along(term_labels)) {
      cols <- col_names[assign == term_idx]
      if (length(cols) == 0L) next
      out[[paste(parameter, term_labels[[term_idx]], sep = ".")]] <- paste(parameter, cols, sep = ".")
    }
  }
  out
}

.resolve_coefficient_terms <- function(terms, coefficient_names, arg = "terms", term_map = NULL) {
  if (is.null(terms)) {
    return(seq_along(coefficient_names))
  }

  if (!is.character(terms)) {
    idx <- terms
    if (length(idx) == 0L || any(is.na(idx)) || any(idx < 1L) || any(idx > length(coefficient_names))) {
      stop(sprintf("'%s' contains unknown coefficient names or indices.", arg), call. = FALSE)
    }
    return(as.integer(idx))
  }

  if (length(terms) == 0L || any(is.na(terms)) || any(!nzchar(terms))) {
    stop(sprintf("'%s' contains unknown coefficient names or indices.", arg), call. = FALSE)
  }

  idx_list <- lapply(terms, function(term) {
    exact <- match(term, coefficient_names)
    if (!is.na(exact)) {
      return(exact)
    }
    if (!is.null(term_map) && term %in% names(term_map)) {
      mapped_idx <- match(term_map[[term]], coefficient_names)
      mapped_idx <- mapped_idx[!is.na(mapped_idx)]
      if (length(mapped_idx) > 0L) {
        return(mapped_idx)
      }
    }
    prefix_idx <- which(startsWith(coefficient_names, term))
    if (length(prefix_idx) == 0L) {
      return(prefix_idx)
    }
    suffix <- substring(coefficient_names[prefix_idx], nchar(term) + 1L)
    main_idx <- prefix_idx[!grepl(":", suffix, fixed = TRUE)]
    if (length(main_idx) > 0L) {
      return(main_idx)
    }
    prefix_idx
  })
  missing <- vapply(idx_list, length, integer(1)) == 0L
  if (any(missing)) {
    stop(sprintf("'%s' contains unknown coefficient names or indices.", arg), call. = FALSE)
  }
  unique(as.integer(unlist(idx_list, use.names = FALSE)))
}

#' Wald tests for fixed coefficients
#'
#' `wald_test()` provides a small reporting-friendly hypothesis-test surface for
#' fitted `gamlss.longitudinal` models. It uses the same variance-covariance
#' route as [summary.gamlss.longitudinal()] and [confint.gamlss.longitudinal()],
#' so numerical-Hessian tests should be reported as approximate.
#'
#' @param object A fitted `gamlss.longitudinal` object.
#' @param terms Optional coefficient names, formula-term names such as
#'   `"mu.treatment"`, coefficient-name prefixes, or numeric indices. When `L`
#'   is `NULL`, these select coefficients for individual tests or a joint test.
#' @param L Optional contrast matrix. Columns must either be named with
#'   coefficient names or have one column per fixed coefficient in model order.
#' @param rhs Null-hypothesis value. Either a scalar or one value per tested row.
#' @param joint Logical; when `TRUE`, test selected `terms` jointly. Contrast
#'   matrices supplied through `L` are always tested jointly.
#' @param method Variance-covariance method passed to [vcov.gamlss.longitudinal()].
#' @param ... Additional arguments passed to [vcov.gamlss.longitudinal()].
#'
#' @return An object of class `gamlss_longitudinal_wald_test`.
#' @export
wald_test <- function(
  object,
  terms = NULL,
  L = NULL,
  rhs = 0,
  joint = FALSE,
  method = "analytical",
  ...
) {
  if (!inherits(object, "gamlss.longitudinal")) {
    stop("'object' must be a fitted 'gamlss.longitudinal' object.", call. = FALSE)
  }
  estimates <- stats::coef(object)
  if (length(estimates) == 0L) {
    stop("No fixed coefficients are available to test.", call. = FALSE)
  }
  vc <- .resolve_vcov(object, extra_args = list(method = method, ...))
  V <- vc$vcov$overall
  if (is.null(V) || !is.matrix(V)) {
    stop("A fixed-effect variance-covariance matrix is required for Wald tests.", call. = FALSE)
  }
  V_names <- colnames(V) %||% rownames(V)
  if (is.null(V_names)) {
    V_names <- names(estimates)
  }
  idx_v <- match(names(estimates), V_names)
  if (any(is.na(idx_v))) {
    stop("Variance-covariance matrix names do not match fitted coefficients.", call. = FALSE)
  }
  V <- V[idx_v, idx_v, drop = FALSE]
  rownames(V) <- colnames(V) <- names(estimates)

  if (!is.null(L)) {
    L <- as.matrix(L)
    if (nrow(L) == 0L || ncol(L) == 0L) {
      stop("'L' must contain at least one contrast row and one coefficient column.", call. = FALSE)
    }
    if (!is.null(colnames(L))) {
      idx <- match(colnames(L), names(estimates))
      if (any(is.na(idx))) {
        stop("Column names in 'L' must match coefficient names.", call. = FALSE)
      }
      L_full <- matrix(0, nrow = nrow(L), ncol = length(estimates))
      colnames(L_full) <- names(estimates)
      rownames(L_full) <- rownames(L)
      L_full[, idx] <- L
      L <- L_full
    } else if (ncol(L) != length(estimates)) {
      stop("Unnamed 'L' must have one column per fixed coefficient.", call. = FALSE)
    } else {
      colnames(L) <- names(estimates)
    }
    if (is.null(rownames(L))) {
      rownames(L) <- paste0("H", seq_len(nrow(L)))
    }
    joint <- TRUE
  } else {
    idx <- .resolve_coefficient_terms(
      terms,
      names(estimates),
      arg = "terms",
      term_map = .fixed_term_coefficient_names(object)
    )
    L <- diag(length(estimates))[idx, , drop = FALSE]
    colnames(L) <- names(estimates)
    rownames(L) <- names(estimates)[idx]
  }

  rhs <- rep(rhs, length.out = nrow(L))
  estimate <- as.numeric(L %*% estimates)
  diff <- estimate - rhs
  LVL <- L %*% V %*% t(L)

  if (isTRUE(joint)) {
    stat <- tryCatch(
      as.numeric(t(diff) %*% solve(LVL, diff)),
      error = function(e) NA_real_
    )
    out <- data.frame(
      hypothesis = paste(rownames(L), collapse = ", "),
      df = nrow(L),
      statistic = stat,
      p_value = stats::pchisq(stat, df = nrow(L), lower.tail = FALSE),
      method = vc$method %||% method,
      stringsAsFactors = FALSE
    )
  } else {
    se <- sqrt(pmax(0, diag(LVL)))
    z <- diff / se
    out <- data.frame(
      term = rownames(L),
      estimate = estimate,
      rhs = rhs,
      std_error = se,
      statistic = z,
      p_value = 2 * stats::pnorm(abs(z), lower.tail = FALSE),
      method = vc$method %||% method,
      stringsAsFactors = FALSE
    )
  }
  attr(out, "joint") <- isTRUE(joint)
  attr(out, "method_requested") <- method
  attr(out, "vcov_method") <- vc$method %||% method
  class(out) <- c("gamlss_longitudinal_wald_test", "data.frame")
  out
}

#' @export
print.gamlss_longitudinal_wald_test <- function(x, digits = max(3, getOption("digits") - 3), ...) {
  cat("\nWald Test for gamlss.longitudinal\n")
  cat("---------------------------------\n")
  cat("Test type:", if (isTRUE(attr(x, "joint"))) "joint" else "individual", "\n")
  cat("VCOV method:", attr(x, "vcov_method") %||% "unknown", "\n\n")
  print.data.frame(x, digits = digits, row.names = FALSE, ...)
  invisible(x)
}

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
#' @param sort Logical; order models by effective degrees of freedom before
#'   computing sequential comparisons.
#'
#' @return An object of class `gamlss_longitudinal_likelihood_compare`.
#' @export
likelihood_compare <- function(..., sort = TRUE) {
  models <- list(...)
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

#' @export
print.gamlss_longitudinal_likelihood_compare <- function(x, digits = max(3, getOption("digits") - 3), ...) {
  cat("\nLikelihood Comparison for gamlss.longitudinal\n")
  cat("--------------------------------------------\n")
  cat("Sequential LR rows compare each model with the previous row.\n\n")
  print.data.frame(x, digits = digits, row.names = FALSE, ...)
  invisible(x)
}

.gl_bootstrap_refit <- function(object, response, fit_args) {
  if (is.null(object$dataset) || is.null(object$formulas_int)) {
    stop(
      "This fitted object does not store the internal dataset/formulas needed for refitting. Refit the model with the current package version.",
      call. = FALSE
    )
  }
  dat <- object$dataset
  if (!"response" %in% names(dat)) {
    stop("Stored model dataset does not contain a 'response' column.", call. = FALSE)
  }
  if (length(response) != nrow(dat)) {
    stop("Bootstrap response length does not match the stored model dataset.", call. = FALSE)
  }
  dat$response <- response
  args <- c(
    list(
      dataset = dat,
      margin_dist = object$margin_dist,
      copula_dist = object$copula_dist,
      time_var = "time",
      subject_var = "subject",
      mu.formula = object$formulas_int$mu,
      sigma.formula = object$formulas_int$sigma,
      nu.formula = object$formulas_int$nu,
      tau.formula = object$formulas_int$tau,
      theta.formula = object$formulas_int$theta,
      zeta.formula = object$formulas_int$zeta,
      start_from = object$par,
      include_dlcopdpar = object$include_dlcopdpar,
      compute_vcov = FALSE,
      verbose = 0
    ),
    fit_args
  )
  do.call(gamlss_longitudinal, args)
}

#' Parametric bootstrap inference for fitted models
#'
#' `bootstrap_inference()` simulates responses from a fitted
#' `gamlss.longitudinal` model, refits the same model to each simulated response,
#' and summarizes the bootstrap distribution of selected fixed coefficients.
#' It is intended for opt-in applied uncertainty checks and should be run with
#' enough replicates outside CRAN-time tests for final reporting.
#'
#' @param object A fitted `gamlss.longitudinal` object.
#' @param R Number of bootstrap replicates.
#' @param terms Optional coefficient names, formula-term names such as
#'   `"mu.treatment"`, coefficient-name prefixes, or numeric indices to
#'   summarize.
#' @param level Confidence level for percentile intervals.
#' @param seed Optional random seed.
#' @param fit_args Optional named list of arguments passed to each refit, such
#'   as `max_outer_iter`, `max_inner_iter`, or convergence tolerances.
#' @param keep_fits Logical; keep successful refitted model objects.
#' @param ... Additional arguments passed to [simulate.gamlss.longitudinal()].
#'
#' @return An object of class `gamlss_longitudinal_bootstrap`.
#' @export
bootstrap_inference <- function(
  object,
  R = 100,
  terms = NULL,
  level = 0.95,
  seed = NULL,
  fit_args = list(),
  keep_fits = FALSE,
  ...
) {
  if (!inherits(object, "gamlss.longitudinal")) {
    stop("'object' must be a fitted 'gamlss.longitudinal' object.", call. = FALSE)
  }
  R <- as.integer(R)
  if (length(R) != 1L || is.na(R) || R < 1L) {
    stop("'R' must be a positive integer.", call. = FALSE)
  }
  if (!is.list(fit_args) || (length(fit_args) > 0L && (is.null(names(fit_args)) || any(names(fit_args) == "")))) {
    stop("'fit_args' must be a named list.", call. = FALSE)
  }
  dots <- list(...)
  if ("simulation_type" %in% names(dots)) {
    simulation_type <- dots$simulation_type
    if (!identical(simulation_type, "copula")) {
      stop("'simulation_type' is no longer supported; bootstrap_inference() simulates from the fitted copula model.", call. = FALSE)
    }
    dots$simulation_type <- NULL
  }
  estimates <- stats::coef(object)
  idx <- .resolve_coefficient_terms(
    terms,
    names(estimates),
    arg = "terms",
    term_map = .fixed_term_coefficient_names(object)
  )
  terms_use <- names(estimates)[idx]
  if (!is.null(seed)) {
    old_seed <- if (exists(".Random.seed", envir = .GlobalEnv, inherits = FALSE)) {
      get(".Random.seed", envir = .GlobalEnv, inherits = FALSE)
    } else {
      NULL
    }
    on.exit({
      if (is.null(old_seed)) {
        if (exists(".Random.seed", envir = .GlobalEnv, inherits = FALSE)) {
          rm(".Random.seed", envir = .GlobalEnv)
        }
      } else {
        assign(".Random.seed", old_seed, envir = .GlobalEnv)
      }
    }, add = TRUE)
    set.seed(seed)
  }

  sim <- do.call(simulate, c(list(object = object, nsim = R), dots))
  boot_coef <- matrix(NA_real_, nrow = R, ncol = length(terms_use))
  colnames(boot_coef) <- terms_use
  errors <- rep(NA_character_, R)
  fits <- if (isTRUE(keep_fits)) vector("list", R) else NULL
  missing_mask <- is.na(object$response)

  for (b in seq_len(R)) {
    response_b <- sim[[b]]
    response_b[missing_mask] <- NA_real_
    fit_b <- tryCatch(
      suppressWarnings(.gl_bootstrap_refit(object, response_b, fit_args = fit_args)),
      error = function(e) e
    )
    if (inherits(fit_b, "error")) {
      errors[[b]] <- conditionMessage(fit_b)
      next
    }
    coef_b <- stats::coef(fit_b)
    boot_coef[b, ] <- coef_b[terms_use]
    if (isTRUE(keep_fits)) {
      fits[[b]] <- fit_b
    }
  }

  alpha <- 1 - level
  q_na <- function(x, prob) {
    x <- x[is.finite(x)]
    if (length(x) == 0L) return(NA_real_)
    as.numeric(stats::quantile(x, probs = prob, names = FALSE))
  }
  summary <- data.frame(
    term = terms_use,
    estimate = as.numeric(estimates[terms_use]),
    bootstrap_mean = colMeans(boot_coef, na.rm = TRUE),
    bootstrap_se = apply(boot_coef, 2, stats::sd, na.rm = TRUE),
    conf.low = apply(boot_coef, 2, q_na, prob = alpha / 2),
    conf.high = apply(boot_coef, 2, q_na, prob = 1 - alpha / 2),
    reps = colSums(is.finite(boot_coef)),
    stringsAsFactors = FALSE
  )

  out <- list(
    summary = summary,
    replicates = as.data.frame(boot_coef, stringsAsFactors = FALSE),
    errors = errors,
    R = R,
    successful_replicates = sum(stats::complete.cases(boot_coef)),
    failed_replicates = sum(!is.na(errors)),
    level = level,
    simulation_type = "copula",
    fits = fits
  )
  class(out) <- "gamlss_longitudinal_bootstrap"
  out
}

#' @export
print.gamlss_longitudinal_bootstrap <- function(x, digits = max(3, getOption("digits") - 3), ...) {
  cat("\nParametric Bootstrap for gamlss.longitudinal\n")
  cat("-------------------------------------------\n")
  cat("Replicates:", x$R, "\n")
  cat("Simulation:", "fitted copula model", "\n")
  cat("Failed refits:", x$failed_replicates, "\n\n")
  print(x$summary, digits = digits, row.names = FALSE)
  invisible(x)
}

.gl_simulate_copula_matrix <- function(object, diag_data, nsim) {
  fit_data <- .copula_v2_fit_data(object)
  n <- nrow(fit_data)
  qfun <- get(paste0("q", diag_data$family), envir = asNamespace("gamlss.dist"), inherits = FALSE)
  out <- matrix(NA_real_, nrow = n, ncol = nsim)

  time_levels <- if (is.factor(fit_data$time)) {
    levels(fit_data$time)[levels(fit_data$time) %in% as.character(unique(fit_data$time))]
  } else {
    u <- unique(fit_data$time)
    if (is.numeric(u) || is.integer(u)) sort(u) else sort(as.character(u))
  }
  time_lookup <- stats::setNames(seq_along(time_levels), as.character(time_levels))
  fit_data$.time_idx <- unname(time_lookup[as.character(fit_data$time)])
  fit_data$.row_id <- seq_len(n)

  split_rows <- split(fit_data[order(fit_data$subject, fit_data$.time_idx), , drop = FALSE], fit_data$subject)

  for (j in seq_len(nsim)) {
    u <- rep(NA_real_, n)
    for (subject_rows in split_rows) {
      row_ids <- subject_rows$.row_id
      if (length(row_ids) == 0L) next
      u[row_ids[[1L]]] <- stats::runif(1L)
      if (length(row_ids) > 1L) {
        for (k in 2:length(row_ids)) {
          left_row <- row_ids[[k - 1L]]
          theta <- fit_data$theta_fit[[left_row]]
          zeta <- fit_data$zeta_fit[[left_row]]
          target <- stats::runif(1L)
          if (!is.finite(theta)) {
            u[row_ids[[k]]] <- target
          } else {
            u[row_ids[[k]]] <- .sim_invert_hfunc1(
              u1 = u[left_row],
              target = target,
              family = object$copula_dist,
              par = theta,
              par2 = if (is.finite(zeta)) zeta else 0
            )
          }
        }
      }
    }
    args <- c(list(p = u), diag_data$params)
    args <- args[names(args) %in% formalArgs(qfun)]
    out[, j] <- do.call(qfun, args)
  }
  out
}

#' Simulate responses from a longitudinal GAMLSS-copula model
#'
#' @param object A fitted `gamlss.longitudinal` object.
#' @param nsim Number of simulated response columns.
#' @param seed Optional random seed.
#' @param newdata Optional new data. Simulation from the fitted copula model is
#'   currently available for fitted data only.
#' @param ... Additional arguments reserved for future methods.
#'
#' @return A data frame with one column per simulation.
#' @importFrom stats simulate
#' @export
simulate.gamlss.longitudinal <- function(
  object,
  nsim = 1,
  seed = NULL,
  newdata = NULL,
  ...
) {
  if (!inherits(object, "gamlss.longitudinal")) {
    stop("'object' must be a fitted 'gamlss.longitudinal' object.", call. = FALSE)
  }
  nsim <- as.integer(nsim)
  if (!is.finite(nsim) || nsim < 1L) {
    stop("'nsim' must be a positive integer.", call. = FALSE)
  }
  dots <- list(...)
  if ("type" %in% names(dots)) {
    type <- dots$type
    if (!identical(type, "copula")) {
      stop("'type' is no longer supported; simulate() always uses the fitted copula model.", call. = FALSE)
    }
  }
  if (!is.null(newdata)) {
    stop("Copula-preserving simulation with 'newdata' is not yet implemented.", call. = FALSE)
  }

  if (!is.null(seed)) {
    old_seed <- if (exists(".Random.seed", envir = .GlobalEnv, inherits = FALSE)) {
      get(".Random.seed", envir = .GlobalEnv, inherits = FALSE)
    } else {
      NULL
    }
    on.exit({
      if (is.null(old_seed)) {
        if (exists(".Random.seed", envir = .GlobalEnv, inherits = FALSE)) {
          rm(".Random.seed", envir = .GlobalEnv)
        }
      } else {
        assign(".Random.seed", old_seed, envir = .GlobalEnv)
      }
    }, add = TRUE)
    set.seed(seed)
  }

  diag_data <- .gl_fitted_distribution(object, newdata = NULL, require_response = TRUE)
  sim_mat <- .gl_simulate_copula_matrix(object, diag_data, nsim)
  if (nrow(sim_mat) != length(object$response)) {
    full_sim_mat <- matrix(NA_real_, nrow = length(object$response), ncol = nsim)
    full_sim_mat[diag_data$keep_index, ] <- sim_mat
    sim_mat <- full_sim_mat
  }

  out <- as.data.frame(sim_mat, stringsAsFactors = FALSE)
  names(out) <- paste0("sim_", seq_len(nsim))
  out
}

.gl_residual_dependence_summary <- function(object, residual_lags = 1L) {
  out <- data.frame(lag = as.integer(residual_lags), normal_score_cor = NA_real_, n_pairs = 0L)
  if (length(residual_lags) == 0L) {
    return(data.frame(lag = integer(0), normal_score_cor = numeric(0), n_pairs = integer(0)))
  }

  copula_spec <- get_copula_dist(object$copula_dist)
  family_name <- .copula_family_code(copula_spec$copula_dist)
  family_num <- tryCatch(.copula_family_code(family_name), error = function(e) NA_character_)

  lag_summary <- tryCatch({
    fit_data <- .copula_v2_fit_data(object)
    rosenblatt_df <- .copula_v2_rosenblatt_series(fit_data, family_num)
    .copula_v2_rosenblatt_lag_summary(rosenblatt_df, lag_values = residual_lags)
  }, error = function(e) data.frame())

  if (nrow(lag_summary) == 0L) {
    return(out)
  }

  data.frame(
    lag = as.integer(lag_summary$lag),
    normal_score_cor = as.numeric(lag_summary$cor_z),
    n_pairs = as.integer(lag_summary$n_pairs),
    stringsAsFactors = FALSE
  )
}

.gl_pit_tail_summary <- function(pit, thresholds = c(0.05, 0.10)) {
  n <- sum(is.finite(pit))
  rows <- lapply(thresholds, function(threshold) {
    lower <- mean(pit <= threshold, na.rm = TRUE)
    upper <- mean(pit >= 1 - threshold, na.rm = TRUE)
    data.frame(
      threshold = threshold,
      lower = lower,
      upper = upper,
      expected = threshold,
      lower_ratio = lower / threshold,
      upper_ratio = upper / threshold,
      stringsAsFactors = FALSE
    )
  })
  out <- do.call(rbind, rows)
  if (n == 0L) {
    out[, c("lower", "upper", "lower_ratio", "upper_ratio")] <- NA_real_
  }
  tail_ratios <- c(out$lower_ratio, out$upper_ratio)
  attr(out, "tail_ratio_max") <- if (any(is.finite(tail_ratios))) {
    max(tail_ratios, na.rm = TRUE)
  } else {
    NA_real_
  }
  out
}

.gl_check_table <- function(summary_obj, scores, pit_stats, tail_stats, lag1_cor,
                            dependence_cor_cutoff = 0.25, vcov_method = NA_character_) {
  row <- function(area, quantity_checked, value, threshold_condition, default,
                  status, message, action) {
    severity <- switch(
      status,
      FAIL = "concern",
      REVIEW = "review",
      PASS = "ok",
      "ok"
    )
    data.frame(
      area = area,
      quantity_checked = quantity_checked,
      value = value,
      threshold_condition = threshold_condition,
      default = default,
      status = status,
      severity = severity,
      message = message,
      action = action,
      stringsAsFactors = FALSE
    )
  }

  converged <- isTRUE(summary_obj$convergence$converged)
  marginal_fail <- is.finite(pit_stats$ks_p_value) && pit_stats$ks_p_value < 0.05
  tail_fail <- is.finite(tail_stats$tail_ratio_max) && tail_stats$tail_ratio_max > 2
  copula_fail <- is.finite(lag1_cor) && abs(lag1_cor) > dependence_cor_cutoff
  variance_review <- identical(vcov_method, "numderiv")

  do.call(rbind, list(
    row(
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
    row(
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
    row(
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
    row(
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
    row(
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

.gl_basic_checks <- function(checks) {
  checks[, c("area", "status", "value", "threshold_condition", "message"), drop = FALSE]
}

.gl_basic_checks_result <- function(checks) {
  if (any(checks$status == "FAIL", na.rm = TRUE)) {
    return("failed")
  }
  if (any(checks$status == "REVIEW", na.rm = TRUE)) {
    return("review")
  }
  "passed"
}

#' Check a fitted longitudinal GAMLSS-copula model
#'
#' `check_model()` turns diagnostics into a compact set of basic automated
#' checks for broad applied use. It does not replace visual inspection, but it
#' provides a stable first pass over convergence, marginal calibration, residual
#' dependence, and scoring summaries.
#'
#' @param object A fitted `gamlss.longitudinal` object.
#' @param include_vcov Logical; include variance-covariance inference metadata
#'   via [summary.gamlss.longitudinal()].
#' @param include_plots Logical; include standard diagnostic plot objects in
#'   `check$plots`. Visual review should usually use `plot(object)` or the
#'   explicit diagnostic helpers instead.
#' @param dependence_cor_cutoff Absolute lag-1 Rosenblatt normal-score residual
#'   correlation above which the dependence check is flagged. The default is a
#'   review threshold rather than a formal hypothesis test.
#' @param ... Passed to [summary.gamlss.longitudinal()] when `include_vcov` is
#'   `TRUE`.
#'
#' @return An object of class `gamlss_longitudinal_check`, including a compact
#'   `basic_checks` table, a full `checks` table, a `warnings` table containing
#'   failed checks, and overall `basic_checks_passed` and `basic_checks_result`
#'   fields.
#' @export
check_model <- function(
  object,
  include_vcov = FALSE,
  include_plots = FALSE,
  dependence_cor_cutoff = 0.25,
  ...
) {
  if (!inherits(object, "gamlss.longitudinal")) {
    stop("'object' must be a fitted 'gamlss.longitudinal' object.", call. = FALSE)
  }
  dependence_cor_cutoff <- as.numeric(dependence_cor_cutoff)
  if (length(dependence_cor_cutoff) != 1L || !is.finite(dependence_cor_cutoff) ||
      dependence_cor_cutoff <= 0 || dependence_cor_cutoff >= 1) {
    stop("'dependence_cor_cutoff' must be a single number between 0 and 1.", call. = FALSE)
  }

  s <- summary(object, include_vcov = include_vcov, ...)
  pit_out <- .gl_pit(object, randomize = FALSE)
  pit <- pmin(pmax(pit_out$pit, 0), 1)
  ks_p <- tryCatch(stats::ks.test(pit, "punif")$p.value, error = function(e) NA_real_)
  pit_stats <- data.frame(
    n = length(pit),
    mean = mean(pit, na.rm = TRUE),
    sd = stats::sd(pit, na.rm = TRUE),
    expected_sd = sqrt(1 / 12),
    ks_p_value = as.numeric(ks_p),
    stringsAsFactors = FALSE
  )
  tail_summary <- .gl_pit_tail_summary(pit)
  tail_stats <- data.frame(
    tail_ratio_max = attr(tail_summary, "tail_ratio_max"),
    stringsAsFactors = FALSE
  )

  scores <- as.data.frame(as.list(proscore(object, type = c("logs", "mae", "mse", "dss"))), stringsAsFactors = FALSE)
  residual_dependence <- .gl_residual_dependence_summary(object, residual_lags = 1L)
  lag1_cor <- residual_dependence$normal_score_cor[match(1L, residual_dependence$lag)]
  if (length(lag1_cor) == 0L) lag1_cor <- NA_real_
  copula_summary <- tryCatch(copula_time_summary(object), error = function(e) NULL)
  vcov_method <- s$vcov$method %||% s$fit$vcov_method %||% NA_character_
  checks <- .gl_check_table(
    summary_obj = list(fit = s$fit, convergence = object$convergence),
    scores = scores,
    pit_stats = pit_stats,
    tail_stats = tail_stats,
    lag1_cor = lag1_cor,
    dependence_cor_cutoff = dependence_cor_cutoff,
    vcov_method = vcov_method
  )
  warnings <- checks[checks$status == "FAIL", , drop = FALSE]
  basic_checks_result <- .gl_basic_checks_result(checks)
  if (nrow(warnings) > 0L) {
    warning(
      paste0(
        "Basic model checks failed: ",
        paste(warnings$area, collapse = ", "),
        ". Review check$warnings and broader diagnostics."
      ),
      call. = FALSE
    )
  }

  out <- list(
    model = s$model,
    fit = s$fit,
    convergence = object$convergence,
    scores = scores,
    pit = pit_stats,
    tail = tail_summary,
    residual_dependence = transform(residual_dependence, cutoff = dependence_cor_cutoff),
    copula = copula_summary,
    basic_checks = .gl_basic_checks(checks),
    basic_checks_passed = !any(checks$status == "FAIL", na.rm = TRUE),
    basic_checks_result = basic_checks_result,
    checks = checks,
    warnings = warnings,
    plots = if (isTRUE(include_plots)) {
      list(
        pithist = pithist(object, plot = TRUE),
        qqrplot = qqrplot(object, plot = TRUE),
        wormplot = wormplot(object, plot = TRUE),
        rootogram = rootogram(object, plot = TRUE)
      )
    } else {
      NULL
    }
  )
  class(out) <- "gamlss_longitudinal_check"
  out
}

#' @export
print.gamlss_longitudinal_check <- function(x, digits = max(3, getOption("digits") - 3), ...) {
  fmt_num <- function(value) {
    if (length(value) == 0L || is.null(value) || !is.finite(value)) {
      return("n/a")
    }
    formatC(value, digits = digits, format = "fg")
  }
  section <- function(title) {
    cat("\n", title, "\n", paste(rep("-", nchar(title)), collapse = ""), "\n", sep = "")
  }

  cat("\n")
  cat("GAMLSS Longitudinal Model Check\n")
  cat("===============================\n")
  cat("Margin: ", x$model$margin_dist, "    Copula: ", x$model$copula_dist, "\n", sep = "")
  cat("LogLik: ", fmt_num(x$fit$logLik), "    Converged: ", if (isTRUE(x$convergence$converged)) "yes" else "no", "\n", sep = "")

  if (!is.null(x$basic_checks)) {
    section("Basic Checks")
    basic_display <- x$basic_checks[, c("area", "status"), drop = FALSE]
    names(basic_display) <- c("Area", "Status")
    print(basic_display, row.names = FALSE, right = FALSE)
    cat("\nResult: ", toupper(x$basic_checks_result), "\n", sep = "")
    cat("Note: these are basic automated checks; broader model diagnostics should also be reviewed.\n")
  }

  section("Scores")
  print(x$scores, digits = digits, row.names = FALSE)

  section("PIT")
  print(x$pit, digits = digits, row.names = FALSE)

  if (!is.null(x$tail)) {
    section("Tail Calibration")
    print(x$tail, digits = digits, row.names = FALSE)
  }

  section("Residual Dependence")
  print(x$residual_dependence, digits = digits, row.names = FALSE)
  cat("\nUse check$checks for thresholds and check$warnings for failed-check details.\n")
  invisible(x)
}

#' Counterfactual marginal effects from fitted distributional parameters
#'
#' @param object A fitted `gamlss.longitudinal` object.
#' @param newdata Data used as the counterfactual baseline.
#' @param variable Single variable to vary.
#' @param values Values to assign to `variable`. Defaults to observed factor
#'   levels for factors/characters or the 25th, 50th, and 75th percentiles for
#'   numeric variables.
#' @param parameter Distributional parameter to summarize, usually `"mu"`.
#' @param reference Optional reference value. Defaults to the first value.
#' @param se.fit Logical; when `TRUE` and `parameter = "mu"`, attach
#'   approximate delta-method standard errors for response-scale averages.
#' @param level Confidence level used when `se.fit = TRUE`.
#' @param vcov_method Variance-covariance method passed to [vcov.gamlss.longitudinal()].
#' @param ... Additional arguments passed to [predict.gamlss.longitudinal()].
#'
#' @return A data frame with average fitted parameter values and contrasts.
#' @importFrom stats predict
#' @export
marginal_effects <- function(
  object,
  newdata,
  variable,
  values = NULL,
  parameter = "mu",
  reference = NULL,
  se.fit = FALSE,
  level = 0.95,
  vcov_method = "analytical",
  ...
) {
  if (!inherits(object, "gamlss.longitudinal")) {
    stop("'object' must be a fitted 'gamlss.longitudinal' object.", call. = FALSE)
  }
  if (missing(newdata) || is.null(newdata)) {
    stop("'newdata' is required for counterfactual marginal effects.", call. = FALSE)
  }
  newdata <- as.data.frame(newdata, stringsAsFactors = FALSE)
  if (!is.character(variable) || length(variable) != 1L || !variable %in% names(newdata)) {
    stop("'variable' must be a single column name in 'newdata'.", call. = FALSE)
  }

  x <- newdata[[variable]]
  if (is.null(values)) {
    values <- if (is.factor(x)) {
      levels(x)
    } else if (is.character(x)) {
      sort(unique(x))
    } else {
      as.numeric(stats::quantile(x, probs = c(0.25, 0.5, 0.75), na.rm = TRUE, names = FALSE))
    }
  }
  if (length(values) == 0L) {
    stop("'values' must contain at least one counterfactual value.", call. = FALSE)
  }
  if (is.null(reference)) {
    reference <- values[[1L]]
  }

  baseline_n <- nrow(newdata)
  add_factor_calibration_rows <- function(nd, template) {
    factor_cols <- names(nd)[vapply(nd, is.factor, logical(1))]
    extras <- list()
    for (fc in factor_cols) {
      levs <- levels(nd[[fc]])
      if (length(levs) < 2L) next
      present <- unique(as.character(nd[[fc]]))
      missing <- setdiff(levs, present)
      if (length(missing) == 0L) next
      for (lev in missing) {
        row <- template[1L, , drop = FALSE]
        row[[fc]] <- factor(lev, levels = levs, ordered = is.ordered(nd[[fc]]))
        for (other_fc in factor_cols) {
          if (other_fc == fc) next
          row[[other_fc]] <- factor(
            as.character(row[[other_fc]]),
            levels = levels(nd[[other_fc]]),
            ordered = is.ordered(nd[[other_fc]])
          )
        }
        extras[[length(extras) + 1L]] <- row
      }
    }
    if (length(extras) == 0L) {
      return(nd)
    }
    rbind(nd, do.call(rbind, extras))
  }

  rows <- lapply(values, function(value) {
    nd <- newdata
    if (is.factor(nd[[variable]])) {
      nd[[variable]] <- factor(
        rep(value, nrow(nd)),
        levels = levels(nd[[variable]]),
        ordered = is.ordered(nd[[variable]])
      )
    } else {
      nd[[variable]] <- value
    }
    nd_pred <- add_factor_calibration_rows(nd, newdata)
    if (isTRUE(se.fit) && identical(parameter, "mu")) {
      pred <- predict(
        object,
        newdata = nd_pred,
        type = "response",
        se.fit = TRUE,
        interval = "none",
        vcov_method = vcov_method,
        ...
      )
      pred <- pred[seq_len(baseline_n), , drop = FALSE]
      estimate <- mean(pred$fit, na.rm = TRUE)
      n_finite <- sum(is.finite(pred$se.fit))
      std_error <- if (n_finite > 0L) {
        sqrt(sum(pred$se.fit^2, na.rm = TRUE)) / n_finite
      } else {
        NA_real_
      }
    } else {
      pred <- predict(object, newdata = nd_pred, type = "parameters", ...)
      pred <- pred[seq_len(baseline_n), , drop = FALSE]
      if (!parameter %in% names(pred)) {
        stop("Parameter '", parameter, "' is not available in model predictions.", call. = FALSE)
      }
      estimate <- mean(pred[[parameter]], na.rm = TRUE)
      std_error <- NA_real_
    }
    data.frame(
      variable = variable,
      value = as.character(value),
      parameter = parameter,
      estimate = estimate,
      std_error = std_error,
      stringsAsFactors = FALSE
    )
  })
  out <- do.call(rbind, rows)
  ref_idx <- match(as.character(reference), out$value)
  if (is.na(ref_idx)) {
    stop("'reference' must be one of the counterfactual values.", call. = FALSE)
  }
  out$reference <- out$value[[ref_idx]]
  out$contrast <- out$estimate - out$estimate[[ref_idx]]
  if (isTRUE(se.fit) && any(is.finite(out$std_error))) {
    alpha <- 1 - level
    z <- stats::qnorm(1 - alpha / 2)
    out$conf.low <- out$estimate - z * out$std_error
    out$conf.high <- out$estimate + z * out$std_error
  }
  rownames(out) <- NULL
  out
}
