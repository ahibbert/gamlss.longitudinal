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

#' preserves the fitted copula dependence structure. With `newdata`,

#' simulation is unconditional and uses the model-implied dependence evaluated

#' on the supplied panel.

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


