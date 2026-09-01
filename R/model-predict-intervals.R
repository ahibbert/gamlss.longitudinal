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

  covariance_contract <- attr(se_values, "inference_contract")
  if (!any(is.finite(se_values))) {
    stop(.gl_inference_unavailable(list(
      status = "unavailable",
      failure_codes = "prediction_standard_errors_nonfinite",
      message = "Prediction inference is unavailable because no finite standard errors were produced."
    )))
  }

  pred$se.fit <- se_values

  if (!identical(interval, "none")) {
    z <- stats::qnorm((1 + level) / 2)

    pred$conf.low <- pred$fit - z * pred$se.fit

    pred$conf.high <- pred$fit + z * pred$se.fit
  }

  out <- pred[c("subject", "time", "response", "fit", "se.fit", intersect(c("conf.low", "conf.high"), names(pred)))]
  mu_names <- names(object$par)[startsWith(names(object$par), "mu.")]
  se_complete <- all(is.finite(se_values))
  covariance_status <- covariance_contract$validity_status %||% "unverified_covariance_contract"
  contract <- .gl_inference_contract(
    "prediction_mu_delta",
    coefficient_names = mu_names,
    method = covariance_contract$method_used %||% vcov_method,
    validity_status = if (se_complete) covariance_status else "partial_standard_error_coverage",
    failure_states = unique(c(
      covariance_contract$observed_failures %||% character(),
      if (!se_complete) "some_prediction_standard_errors_nonfinite" else character()
    ))
  )
  contract$covariance_contract <- covariance_contract
  contract$method_requested <- covariance_contract$method_requested %||% vcov_method
  contract$method_used <- covariance_contract$method_used %||% vcov_method
  contract$fallback_used <- covariance_contract$fallback_used %||% FALSE
  contract$diagnostics <- covariance_contract$diagnostics %||% NULL
  contract$standard_error_rows <- list(
    total = length(se_values), finite = sum(is.finite(se_values))
  )
  .gl_attach_inference_contract(out, contract)
}
