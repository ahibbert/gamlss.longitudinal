#' Screen candidate marginal distributions
#'
#' `screen_margin()` is a lightweight wrapper around [gamlss::fitDist()] for
#' the recommended longitudinal workflow: choose a plausible marginal family,
#' then fit dependence with [gamlss_longitudinal()].
#'
#' @param response Numeric response vector. For the common
#'   `screen_margin(dat, response_var = "y")` call, a data frame supplied here
#'   is treated as `data`.
#' @param data Optional data frame containing the response.
#' @param response_var Optional response column name in `data`.
#' @param type Optional `gamlss::fitDist()` type. If `NULL`, a simple heuristic
#'   uses `"counts"` for non-negative integer responses, `"realplus"` for
#'   positive continuous responses, and `"realAll"` otherwise.
#' @param families Optional character vector used to filter the returned table.
#' @param try.gamlss,trace,... Passed to [gamlss::fitDist()].
#'
#' @return A data frame ordered by AIC. The selected family is also stored in
#'   the `"selected"` attribute.
#' @export
screen_margin <- function(
  response = NULL,
  data = NULL,
  response_var = NULL,
  type = NULL,
  families = NULL,
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
  response <- as.numeric(response)
  response <- response[is.finite(response)]
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

  fit <- gamlss::fitDist(
    response,
    type = type,
    try.gamlss = try.gamlss,
    trace = trace,
    ...
  )

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
  out <- out[order(out$AIC), , drop = FALSE]
  rownames(out) <- NULL
  out$rank <- seq_len(nrow(out))
  out$delta_AIC <- out$AIC - min(out$AIC, na.rm = TRUE)
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
  class(out) <- c("margin_screen", class(out))
  out
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
  print(utils::head(as.data.frame(x), n = n), row.names = FALSE)
  invisible(x)
}

#' Fit the recommended longitudinal GAMLSS-copula workflow
#'
#' `fit_longitudinal()` is a golden-path wrapper around [gamlss_longitudinal()]
#' with adoption-friendly defaults. It optionally screens the copula family
#' from a pseudo-observation column and stores workflow metadata on the fitted
#' object.
#'
#' @param dataset Long-format data frame.
#' @param margin_dist Marginal `gamlss.dist` family object.
#' @param time_var,subject_var Column names identifying time and subject.
#' @param mu.formula Formula for the marginal location/mean parameter.
#' @param sigma.formula,nu.formula,tau.formula Marginal parameter formulas.
#' @param theta.formula,zeta.formula Copula parameter formulas.
#' @param copula_dist Optional copula family code. If `NULL`,
#'   `fit_longitudinal()` selects from `u_var` when available, otherwise uses
#'   `"N"` as a conservative default.
#' @param auto_select_copula Logical; when `TRUE` and pseudo-observations are
#'   available, call [select_copula()].
#' @param u_var Optional pseudo-observation column used for copula screening.
#'   If omitted and `dataset` has a column named `"u"`, that column is used.
#' @param copula_families Candidate copula family codes for screening.
#' @param run_checks Logical; attach `check_model(..., include_plots = FALSE)`
#'   to the returned object.
#' @param warm_start_joint Logical; default `FALSE` for the high-level workflow
#'   because it favors reliability over a faster but more fragile warm-start
#'   phase. Advanced users can set this to `TRUE`.
#' @param method,include_dlcopdpar,compute_vcov,verbose Defaults passed to
#'   [gamlss_longitudinal()].
#' @param ... Additional arguments passed to [gamlss_longitudinal()].
#'
#' @return A fitted `gamlss.longitudinal` object with a `workflow` component.
#' @export
fit_longitudinal <- function(
  dataset,
  margin_dist,
  time_var,
  subject_var,
  mu.formula,
  sigma.formula = ~ 1,
  nu.formula = ~ 1,
  tau.formula = ~ 1,
  theta.formula = ~ 1,
  zeta.formula = ~ 1,
  copula_dist = NULL,
  auto_select_copula = TRUE,
  u_var = NULL,
  copula_families = c("N", "C", "F", "G", "J", "t"),
  run_checks = TRUE,
  method = "RS",
  include_dlcopdpar = TRUE,
  compute_vcov = FALSE,
  warm_start_joint = FALSE,
  verbose = 0,
  ...
) {
  dataset <- as.data.frame(dataset, stringsAsFactors = FALSE)
  if (is.null(u_var) && "u" %in% names(dataset)) {
    u_var <- "u"
  }

  copula_screen <- NULL
  copula_source <- "user"
  if (is.null(copula_dist)) {
    if (isTRUE(auto_select_copula) && !is.null(u_var) && u_var %in% names(dataset)) {
      copula_screen <- select_copula(
        data = dataset,
        u_var = u_var,
        subject_var = subject_var,
        time_var = time_var,
        families = copula_families
      )
      copula_dist <- attr(copula_screen, "selected")
      copula_source <- paste0("select_copula:", u_var)
    } else {
      copula_dist <- "N"
      copula_source <- "default:N"
    }
  }

  fit <- gamlss_longitudinal(
    dataset = dataset,
    margin_dist = margin_dist,
    copula_dist = copula_dist,
    time_var = time_var,
    subject_var = subject_var,
    mu.formula = mu.formula,
    sigma.formula = sigma.formula,
    nu.formula = nu.formula,
    tau.formula = tau.formula,
    theta.formula = theta.formula,
    zeta.formula = zeta.formula,
    include_dlcopdpar = include_dlcopdpar,
    method = method,
    compute_vcov = compute_vcov,
    warm_start_joint = warm_start_joint,
    verbose = verbose,
    ...
  )

  fit$workflow <- list(
    interface = "fit_longitudinal",
    copula_source = copula_source,
    copula_screen = copula_screen,
    u_var = u_var,
    run_checks = isTRUE(run_checks)
  )
  if (isTRUE(run_checks)) {
    fit$workflow$check <- check_model(fit, include_vcov = FALSE, include_plots = FALSE)
  }
  fit
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

#' Predict from a longitudinal GAMLSS-copula model
#'
#' @param object A fitted `gamlss.longitudinal` object.
#' @param newdata Optional new data. If omitted, fitted rows are used.
#' @param type Prediction type: `"response"` returns the fitted `mu` parameter,
#'   `"parameters"` returns all fitted marginal distribution parameters,
#'   `"quantile"` returns fitted marginal quantiles, `"cdf"`/`"density"`
#'   evaluate the fitted marginal CDF or density, and `"probability"` returns
#'   probabilities below or above a threshold.
#' @param probs Quantile probabilities when `type = "quantile"`.
#' @param q Threshold value for `type = "cdf"` or `type = "probability"`.
#'   Defaults to observed responses for `type = "cdf"` when available.
#' @param y Evaluation value for `type = "density"`. Defaults to observed
#'   responses when available.
#' @param direction Probability direction for `type = "probability"`.
#' @param se.fit Logical; for `type = "response"`, return approximate
#'   delta-method standard errors for the fitted response mean.
#' @param interval Interval type. `"confidence"` adds response-scale confidence
#'   limits when `type = "response"`.
#' @param level Confidence level for `interval = "confidence"`.
#' @param vcov_method Variance-covariance method passed to [vcov.gamlss.longitudinal()].
#' @param ... Additional arguments reserved for future methods.
#'
#' @return A numeric vector for `type = "response"` unless standard errors or
#'   intervals are requested; a data frame otherwise.
#' @export
predict.gamlss.longitudinal <- function(
  object,
  newdata = NULL,
  type = c("response", "parameters", "quantile", "cdf", "density", "probability"),
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

  if (type == "response") {
    fit_values <- if ("mu" %in% names(params)) as.numeric(params$mu) else as.numeric(params[[1L]])
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
    if ("mu" %in% names(params)) {
      return(as.numeric(params$mu))
    }
    return(as.numeric(params[[1L]]))
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

#' Wald tests for fixed coefficients
#'
#' `wald_test()` provides a small reporting-friendly hypothesis-test surface for
#' fitted `gamlss.longitudinal` models. It uses the same variance-covariance
#' route as [summary.gamlss.longitudinal()] and [confint.gamlss.longitudinal()],
#' so numerical-Hessian tests should be reported as approximate.
#'
#' @param object A fitted `gamlss.longitudinal` object.
#' @param terms Optional coefficient names or numeric indices. When `L` is
#'   `NULL`, these select coefficients for individual tests or a joint test.
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
    idx <- if (is.null(terms)) {
      seq_along(estimates)
    } else if (is.character(terms)) {
      match(terms, names(estimates))
    } else {
      terms
    }
    if (any(is.na(idx)) || any(idx < 1L) || any(idx > length(estimates))) {
      stop("'terms' contains unknown coefficient names or indices.", call. = FALSE)
    }
    idx <- as.integer(idx)
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
#' @param terms Optional coefficient names or numeric indices to summarize.
#' @param level Confidence level for percentile intervals.
#' @param simulation_type Passed to [simulate.gamlss.longitudinal()]. Use
#'   `"copula"` to preserve fitted dependence or `"marginal"` for independent
#'   fitted-margin simulation.
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
  simulation_type = c("copula", "marginal"),
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
  simulation_type <- match.arg(simulation_type)
  estimates <- stats::coef(object)
  idx <- if (is.null(terms)) {
    seq_along(estimates)
  } else if (is.character(terms)) {
    match(terms, names(estimates))
  } else {
    terms
  }
  if (any(is.na(idx)) || any(idx < 1L) || any(idx > length(estimates))) {
    stop("'terms' contains unknown coefficient names or indices.", call. = FALSE)
  }
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

  sim <- simulate(object, nsim = R, type = simulation_type, ...)
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
    successful_replicates = colSums(is.finite(boot_coef)),
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
    simulation_type = simulation_type,
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
  cat("Simulation type:", x$simulation_type, "\n")
  cat("Failed refits:", x$failed_replicates, "\n\n")
  print(x$summary, digits = digits, row.names = FALSE)
  invisible(x)
}

.gl_simulate_marginal_matrix <- function(object, diag_data, nsim) {
  qfun <- get(paste0("q", diag_data$family), envir = asNamespace("gamlss.dist"), inherits = FALSE)
  out <- matrix(NA_real_, nrow = length(diag_data$subject), ncol = nsim)
  for (j in seq_len(nsim)) {
    args <- c(list(p = stats::runif(nrow(out))), diag_data$params)
    args <- args[names(args) %in% formalArgs(qfun)]
    out[, j] <- do.call(qfun, args)
  }
  out
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
#' @param newdata Optional new data. Copula-preserving simulation is currently
#'   available for fitted data only; `newdata` uses marginal simulation.
#' @param type `"copula"` preserves the fitted first-order dependence on fitted
#'   data. `"marginal"` simulates independently from fitted marginal
#'   distributions.
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
  type = c("copula", "marginal"),
  ...
) {
  type <- match.arg(type)
  if (!inherits(object, "gamlss.longitudinal")) {
    stop("'object' must be a fitted 'gamlss.longitudinal' object.", call. = FALSE)
  }
  nsim <- as.integer(nsim)
  if (!is.finite(nsim) || nsim < 1L) {
    stop("'nsim' must be a positive integer.", call. = FALSE)
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

  diag_data <- .gl_fitted_distribution(object, newdata = newdata, require_response = FALSE)
  sim_mat <- if (type == "copula" && is.null(newdata)) {
    .gl_simulate_copula_matrix(object, diag_data, nsim)
  } else {
    if (type == "copula" && !is.null(newdata)) {
      warning("Copula-preserving simulation with 'newdata' is not yet implemented; using marginal simulation.", call. = FALSE)
    }
    .gl_simulate_marginal_matrix(object, diag_data, nsim)
  }

  out <- as.data.frame(sim_mat, stringsAsFactors = FALSE)
  names(out) <- paste0("sim_", seq_len(nsim))
  out
}

.gl_lag1_residual_cor <- function(object) {
  pit_out <- .gl_pit(object, randomize = FALSE)
  z <- stats::qnorm(pmin(pmax(pit_out$pit, .Machine$double.eps), 1 - .Machine$double.eps))
  dat <- data.frame(
    subject = pit_out$diag$subject,
    time = pit_out$diag$time,
    z = z,
    stringsAsFactors = FALSE
  )
  dat <- dat[order(dat$subject, dat$time), , drop = FALSE]
  pairs <- lapply(split(dat, dat$subject), function(x) {
    if (nrow(x) < 2L) return(NULL)
    data.frame(z_prev = x$z[-nrow(x)], z_curr = x$z[-1L])
  })
  pairs <- do.call(rbind, pairs[!vapply(pairs, is.null, logical(1))])
  if (is.null(pairs) || nrow(pairs) < 3L) {
    return(NA_real_)
  }
  stats::cor(pairs$z_prev, pairs$z_curr, use = "complete.obs")
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

.gl_check_table <- function(summary_obj, scores, pit_stats, tail_stats, lag1_cor) {
  checks <- list()
  if (!isTRUE(summary_obj$convergence$converged)) {
    checks[[length(checks) + 1L]] <- data.frame(
      area = "convergence",
      severity = "concern",
      message = "Convergence was not confirmed.",
      action = "Refit with more iterations, different starts, or a simpler specification.",
      stringsAsFactors = FALSE
    )
  }
  if (is.finite(pit_stats$ks_p_value) && pit_stats$ks_p_value < 0.05) {
    checks[[length(checks) + 1L]] <- data.frame(
      area = "marginal_calibration",
      severity = "concern",
      message = "The marginal distribution is off by the PIT uniformity screen.",
      action = "Inspect PIT, QQ, worm, and rootogram diagnostics; try a richer margin or covariate specification.",
      stringsAsFactors = FALSE
    )
  }
  if (is.finite(tail_stats$tail_ratio_max) && tail_stats$tail_ratio_max > 2) {
    checks[[length(checks) + 1L]] <- data.frame(
      area = "tail_calibration",
      severity = "review",
      message = "Tail observations occur more often than the fitted margin expects.",
      action = "Inspect lower/upper PIT tails and consider heavier-tailed or asymmetric margins.",
      stringsAsFactors = FALSE
    )
  }
  if (is.finite(lag1_cor) && abs(lag1_cor) > 0.15) {
    checks[[length(checks) + 1L]] <- data.frame(
      area = "dependence",
      severity = "concern",
      message = "Dependence remains after the copula in normal-score residuals.",
      action = "Consider a different copula family, time-varying dependence, or richer serial structure.",
      stringsAsFactors = FALSE
    )
  }
  if (!is.null(summary_obj$fit$vcov_method) && identical(summary_obj$fit$vcov_method, "numderiv")) {
    checks[[length(checks) + 1L]] <- data.frame(
      area = "inference",
      severity = "note",
      message = "Variance-covariance inference used the numerical Hessian path.",
      action = "Cite intervals and tests as approximate numerical-Hessian inference.",
      stringsAsFactors = FALSE
    )
  }
  if (length(checks) == 0L) {
    checks[[1L]] <- data.frame(
      area = "overall",
      severity = "ok",
      message = "No major automatic diagnostic warnings were triggered.",
      action = "Continue with visual diagnostics and subject-matter checks.",
      stringsAsFactors = FALSE
    )
  }
  do.call(rbind, checks)
}

.gl_decision_messages <- function(checks) {
  flagged <- checks[checks$severity %in% c("concern", "review", "note"), , drop = FALSE]
  if (nrow(flagged) == 0L) {
    return("No major automatic diagnostic warnings were triggered.")
  }
  paste(flagged$message, flagged$action)
}

.gl_model_recommendation <- function(checks) {
  if (any(checks$severity == "concern", na.rm = TRUE)) {
    return(data.frame(
      role = "revise_before_primary",
      message = "At least one diagnostic concern should be addressed before using this model as the primary analysis.",
      stringsAsFactors = FALSE
    ))
  }
  if (any(checks$severity %in% c("review", "note"), na.rm = TRUE)) {
    return(data.frame(
      role = "primary_candidate_with_caveats",
      message = "No automatic diagnostic concern was triggered, but review or inference notes should be reported.",
      stringsAsFactors = FALSE
    ))
  }
  data.frame(
    role = "primary_candidate",
    message = "No major automatic diagnostic warnings were triggered; continue with visual and subject-matter checks.",
    stringsAsFactors = FALSE
  )
}

#' Check a fitted longitudinal GAMLSS-copula model
#'
#' `check_model()` turns diagnostics into a compact set of decisions for broad
#' applied use. It does not replace visual inspection, but it provides a stable
#' first pass over convergence, marginal calibration, residual dependence, and
#' scoring summaries.
#'
#' @param object A fitted `gamlss.longitudinal` object.
#' @param include_vcov Logical; include variance-covariance inference metadata
#'   via [summary.gamlss.longitudinal()].
#' @param include_plots Logical; include standard diagnostic plot objects.
#' @param ... Passed to [summary.gamlss.longitudinal()] when `include_vcov` is
#'   `TRUE`.
#'
#' @return An object of class `gamlss_longitudinal_check`, including a full
#'   `checks` table, a filtered `warnings` table for diagnostics that should be
#'   reviewed or reported, and a compact `recommendation` row for model role
#'   decisions.
#' @export
check_model <- function(object, include_vcov = FALSE, include_plots = TRUE, ...) {
  if (!inherits(object, "gamlss.longitudinal")) {
    stop("'object' must be a fitted 'gamlss.longitudinal' object.", call. = FALSE)
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
  lag1_cor <- .gl_lag1_residual_cor(object)
  copula_summary <- tryCatch(copula_time_summary(object), error = function(e) NULL)
  checks <- .gl_check_table(
    summary_obj = list(fit = s$fit, convergence = object$convergence),
    scores = scores,
    pit_stats = pit_stats,
    tail_stats = tail_stats,
    lag1_cor = lag1_cor
  )
  warnings <- checks[checks$severity %in% c("concern", "review", "note"), , drop = FALSE]

  out <- list(
    model = s$model,
    fit = s$fit,
    convergence = object$convergence,
    scores = scores,
    pit = pit_stats,
    tail = tail_summary,
    residual_dependence = data.frame(lag = 1L, normal_score_cor = lag1_cor),
    copula = copula_summary,
    checks = checks,
    warnings = warnings,
    recommendation = .gl_model_recommendation(checks),
    decisions = .gl_decision_messages(checks),
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
  cat("\nGAMLSS Longitudinal Model Check\n")
  cat("------------------------------\n")
  cat("Margin distribution:", x$model$margin_dist, "\n")
  cat("Copula distribution:", x$model$copula_dist, "\n")
  cat("Converged:", if (isTRUE(x$convergence$converged)) "yes" else "no", "\n")
  if (!is.null(x$fit$logLik) && is.finite(x$fit$logLik)) {
    cat("Joint logLik:", formatC(x$fit$logLik, digits = digits, format = "fg"), "\n")
  }
  cat("\nScoring summary:\n")
  print(x$scores, digits = digits, row.names = FALSE)
  cat("\nPIT summary:\n")
  print(x$pit, digits = digits, row.names = FALSE)
  if (!is.null(x$tail)) {
    cat("\nTail calibration:\n")
    print(x$tail, digits = digits, row.names = FALSE)
  }
  cat("\nResidual dependence:\n")
  print(x$residual_dependence, digits = digits, row.names = FALSE)
  if (!is.null(x$checks)) {
    cat("\nDiagnostic checks:\n")
    print(x$checks, digits = digits, row.names = FALSE)
  }
  if (!is.null(x$warnings) && nrow(x$warnings) > 0L) {
    cat("\nWarnings to report:\n")
    print(x$warnings, digits = digits, row.names = FALSE)
  }
  if (!is.null(x$recommendation)) {
    cat("\nRecommendation:\n")
    print(x$recommendation, digits = digits, row.names = FALSE)
  }
  cat("\nDecisions:\n")
  for (msg in x$decisions) {
    cat("* ", msg, "\n", sep = "")
  }
  invisible(x)
}

#' @export
plot.gamlss_longitudinal_check <- function(
  x,
  which = c("pithist", "qqrplot", "wormplot", "rootogram"),
  ...
) {
  if (!inherits(x, "gamlss_longitudinal_check")) {
    stop("'x' must be a 'gamlss_longitudinal_check' object.", call. = FALSE)
  }
  which <- match.arg(which, several.ok = TRUE)
  if (is.null(x$plots) || length(x$plots) == 0L) {
    message("No plot objects stored. Re-run check_model(..., include_plots = TRUE).")
    return(invisible(x))
  }
  available <- intersect(which, names(x$plots))
  missing <- setdiff(which, names(x$plots))
  if (length(missing) > 0L) {
    warning(
      "Requested diagnostic plot(s) not stored: ",
      paste(missing, collapse = ", "),
      call. = FALSE
    )
  }
  if (length(available) == 0L) {
    message("No requested diagnostic plot objects are available.")
    return(invisible(x))
  }
  for (nm in available) {
    if (!is.null(x$plots[[nm]])) {
      print(x$plots[[nm]])
    }
  }
  invisible(x)
}

#' Effects for fitted longitudinal GAMLSS-copula models
#'
#' `effects.gamlss.longitudinal()` is a convenience alias for
#' [marginal_effects()] so users can reach for the familiar `effects()` generic
#' when interpreting fitted distributional parameters.
#'
#' @param object A fitted `gamlss.longitudinal` object.
#' @param ... Arguments passed to [marginal_effects()], including `newdata`,
#'   `variable`, `values`, `parameter`, and `se.fit`.
#'
#' @return A data frame from [marginal_effects()].
#' @export
effects.gamlss.longitudinal <- function(object, ...) {
  marginal_effects(object, ...)
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
