# Prediction evaluates coefficients in their fitted smooth coordinates.
# No new smoothCon call, legacy guessing, refit or covariance is used here.
.gl_prediction_basis_stop <- function(code, message) {
  stop(structure(list(message = message, call = NULL, failure_code = code),
    class = c("gamlss_longitudinal_prediction_unavailable", "error", "condition")))
}

.gl_build_prediction_smooth_model_matrices <- function(parameter, formula, data_for_par, reference) {
  if (!is.list(reference) || !is.list(reference$x) || !is.list(reference$s) ||
      !(is.matrix(reference$x[[parameter]]) || is.data.frame(reference$x[[parameter]])))
    .gl_prediction_basis_stop("prediction_fitted_reference_unavailable",
      paste("A retained fitted parameter design is required:", parameter))
  labels <- attr(stats::terms(formula), "term.labels")
  expected <- labels[grepl("^\\s*s\\(", labels)]
  retained <- reference$s[[parameter]]
  if (!length(expected)) {
    if (length(retained)) .gl_prediction_basis_stop("prediction_smooth_term_mismatch",
      paste("Fitted smooth terms do not match the prediction formula:", parameter))
    return(NULL)
  }
  if (!is.list(retained) || !identical(names(retained), expected))
    .gl_prediction_basis_stop("prediction_smooth_term_mismatch",
      paste("An exact fitted smooth-term map is required:", parameter))
  result <- setNames(vector("list", length(expected)), expected)
  for (term in expected) {
    B <- retained[[term]]
    spec <- attr(B, "gl_smooth_spec", exact = TRUE)
    if (!is.matrix(B) || !is.numeric(B) || any(!is.finite(B)) ||
        !identical(attr(B, "gl_smooth_spec_version", exact = TRUE), 1L) ||
        !identical(attr(B, "gl_smooth_construction_count", exact = TRUE), 1L) ||
        !inherits(spec, "mgcv.smooth") || !is.matrix(spec$X) ||
        !identical(dim(spec$X), dim(B)) ||
        !identical(as.vector(spec$X), as.vector(B))) {
      .gl_prediction_basis_stop("prediction_fitted_smooth_spec_unavailable",
        paste("Prediction requires the retained fitted smooth construction; legacy or changed bases remain unavailable:",
          parameter, term))
    }
    value <- mgcv::PredictMat(spec, data_for_par)
    if (!is.matrix(value) || !is.numeric(value) ||
        !identical(dim(value), c(nrow(data_for_par), ncol(B))) || any(!is.finite(value)))
      .gl_prediction_basis_stop("prediction_fitted_smooth_matrix_invalid",
        paste("The fitted smooth construction returned an invalid prediction matrix:", parameter, term))
    colnames(value) <- colnames(B)
    result[[term]] <- value
  }
  result
}
