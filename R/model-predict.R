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
    ...) {
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
      return(.gl_prediction_interval_frame(
        object = object,
        newdata = newdata,
        fit_values = fit_values,
        interval = interval,
        level = level,
        vcov_method = vcov_method,
        require_response = FALSE,
        ...
      ))
    }

    return(fit_values)
  }

  pred <- .gl_prediction_frame(object, newdata = newdata, require_response = require_response)

  if (type == "parameters") {
    return(pred)
  }

  if (type == "quantile") {
    return(.gl_prediction_quantile_frame(pred, params, diag_data$family, probs))
  }

  if (type == "cdf") {
    return(.gl_prediction_cdf_frame(pred, diag_data, q, require_response))
  }

  if (type == "probability") {
    return(.gl_prediction_probability_frame(pred, diag_data, q, direction, require_response))
  }

  .gl_prediction_density_frame(pred, diag_data, y, require_response)
}
