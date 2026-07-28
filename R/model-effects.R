#' Marginal effects from fitted distributional parameters

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
    ...) {
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
    values <- .gl_effect_counterfactual_values(x)
  }

  if (length(values) == 0L) {
    stop("'values' must contain at least one counterfactual value.", call. = FALSE)
  }

  if (is.null(reference)) {
    reference <- values[[1L]]
  }

  rows <- lapply(values, function(value) {
    .gl_effect_counterfactual_row(
      object = object,
      newdata = newdata,
      variable = variable,
      value = value,
      parameter = parameter,
      se.fit = se.fit,
      vcov_method = vcov_method,
      ...
    )
  })

  .gl_finalize_marginal_effects(
    rows = rows,
    reference = reference,
    se.fit = se.fit,
    level = level
  )
}
