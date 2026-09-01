#' Evaluate one counterfactual marginal-effect row
#'
#' @param object Fitted model object.
#' @param newdata Baseline counterfactual data.
#' @param variable Single variable to set.
#' @param value Counterfactual value for `variable`.
#' @param parameter Distributional parameter to summarize.
#' @param se.fit Logical; compute standard-error summary for `parameter = "mu"`.
#' @param vcov_method Variance-covariance method passed to prediction.
#' @param predict_fn Prediction function, dependency-injected for tests.
#' @param ... Additional prediction arguments.
#' @return One-row data frame with estimate and optional standard error.
#' @noRd
.gl_effect_counterfactual_row <- function(
    object,
    newdata,
    variable,
    value,
    parameter,
    se.fit,
    vcov_method,
    predict_fn = predict,
    ...) {
  baseline_n <- nrow(newdata)
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

  nd_pred <- .gl_effect_add_factor_calibration_rows(nd, newdata)

  if (isTRUE(se.fit) && identical(parameter, "mu")) {
    pred <- predict_fn(
      object,
      newdata = nd_pred,
      type = "response",
      se.fit = TRUE,
      interval = "none",
      vcov_method = vcov_method,
      ...
    )

    prediction_contract <- attr(pred, "inference_contract")

    pred <- pred[seq_len(baseline_n), , drop = FALSE]

    estimate <- mean(pred$fit, na.rm = TRUE)

    n_finite <- sum(is.finite(pred$se.fit))

    if (n_finite == 0L) {
      stop(.gl_inference_unavailable(list(
        status = "unavailable",
        failure_codes = "marginal_effect_standard_errors_nonfinite",
        message = "Marginal-effect inference is unavailable because no finite rowwise prediction standard errors were produced."
      )))
    }
    std_error <- sqrt(sum(pred$se.fit^2, na.rm = TRUE)) / n_finite
  } else {
    pred <- predict_fn(object, newdata = nd_pred, type = "parameters", ...)

    pred <- pred[seq_len(baseline_n), , drop = FALSE]

    if (!parameter %in% names(pred)) {
      stop("Parameter '", parameter, "' is not available in model predictions.", call. = FALSE)
    }

    estimate <- mean(pred[[parameter]], na.rm = TRUE)

    std_error <- NA_real_
  }

  out <- data.frame(
    variable = variable,
    value = as.character(value),
    parameter = parameter,
    estimate = estimate,
    std_error = std_error,
    stringsAsFactors = FALSE
  )
  if (isTRUE(se.fit) && identical(parameter, "mu")) {
    attr(out, "prediction_inference_contract") <- prediction_contract
  }
  out
}
