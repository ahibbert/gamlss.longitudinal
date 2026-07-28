#' Build fitted marginal quantile prediction frame
#'
#' @param pred Base prediction frame from `.gl_prediction_frame()`.
#' @param params Fitted marginal distribution parameters.
#' @param family Marginal family name.
#' @param probs Quantile probabilities.
#' @return Data frame with subject, time, response, and quantile columns.
#' @noRd
.gl_prediction_quantile_frame <- function(pred, params, family, probs) {
  q_df <- .gl_quantile_columns(params, family, probs)

  cbind(pred[c("subject", "time", "response")], q_df)
}

#' Build fitted marginal CDF prediction frame
#'
#' @param pred Base prediction frame from `.gl_prediction_frame()`.
#' @param diag_data Fitted distribution data from `.gl_fitted_distribution()`.
#' @param q Evaluation values.
#' @param require_response Whether observed responses are required.
#' @return Data frame with subject, time, response, q, and CDF columns.
#' @noRd
.gl_prediction_cdf_frame <- function(pred, diag_data, q, require_response) {
  q_use <- .gl_prediction_eval_values(q, "q", "cdf", diag_data, require_response)

  pred$q <- q_use

  pred$cdf <- .gl_call_family_fun("p", diag_data$family, q_use, diag_data$params)

  pred[c("subject", "time", "response", "q", "cdf")]
}

#' Build fitted marginal threshold-probability prediction frame
#'
#' @param pred Base prediction frame from `.gl_prediction_frame()`.
#' @param diag_data Fitted distribution data from `.gl_fitted_distribution()`.
#' @param q Threshold values.
#' @param direction Probability direction, `"below"` or `"above"`.
#' @param require_response Whether observed responses are required.
#' @return Data frame with subject, time, response, threshold, direction, and probability.
#' @noRd
.gl_prediction_probability_frame <- function(pred, diag_data, q, direction, require_response) {
  if (is.null(q)) {
    stop("'q' is required for type = 'probability'.", call. = FALSE)
  }

  q_use <- .gl_prediction_eval_values(q, "q", "probability", diag_data, require_response)

  cdf <- .gl_call_family_fun("p", diag_data$family, q_use, diag_data$params)

  pred$q <- q_use

  pred$direction <- direction

  pred$probability <- if (direction == "below") cdf else 1 - cdf

  pred[c("subject", "time", "response", "q", "direction", "probability")]
}

#' Build fitted marginal density prediction frame
#'
#' @param pred Base prediction frame from `.gl_prediction_frame()`.
#' @param diag_data Fitted distribution data from `.gl_fitted_distribution()`.
#' @param y Evaluation values.
#' @param require_response Whether observed responses are required.
#' @return Data frame with subject, time, response, y, and density columns.
#' @noRd
.gl_prediction_density_frame <- function(pred, diag_data, y, require_response) {
  y_use <- .gl_prediction_eval_values(y, "y", "density", diag_data, require_response)

  pred$y <- y_use

  pred$density <- .gl_call_family_fun("d", diag_data$family, y_use, diag_data$params)

  pred[c("subject", "time", "response", "y", "density")]
}
