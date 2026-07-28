#' Build fitted-value prediction frame with standard errors and intervals
#'
#' @param object Fitted model object.
#' @param newdata Optional prediction data.
#' @param fit_values Fitted values already computed for the requested estimand.
#' @param interval Interval type.
#' @param level Confidence level.
#' @param vcov_method Variance-covariance method for prediction standard errors.
#' @param require_response Logical passed to `.gl_prediction_frame()`.
#' @param se_fn Prediction standard-error helper.
#' @param ... Additional arguments passed to `se_fn`.
#' @return Prediction data frame with fit, standard error, and optional limits.
#' @noRd
.gl_prediction_interval_frame <- function(
    object,
    newdata,
    fit_values,
    interval,
    level,
    vcov_method,
    require_response = FALSE,
    se_fn = .gl_predict_response_se,
    ...) {
  pred <- .gl_prediction_frame(object, newdata = newdata, require_response = require_response)

  pred$fit <- fit_values

  se_values <- se_fn(object, newdata = newdata, method = vcov_method, ...)

  pred$se.fit <- se_values

  if (!identical(interval, "none")) {
    z <- stats::qnorm((1 + level) / 2)

    pred$conf.low <- pred$fit - z * pred$se.fit

    pred$conf.high <- pred$fit + z * pred$se.fit
  }

  pred[c("subject", "time", "response", "fit", "se.fit", intersect(c("conf.low", "conf.high"), names(pred)))]
}
