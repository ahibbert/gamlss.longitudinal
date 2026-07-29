#' Check response missingness against observed predictors
#'
#' `check_missingness()` screens whether response missingness is associated with
#' observed variables by fitting a logistic model for the missing-response
#' indicator. This is a practical MCAR/MAR review: associations with observed
#' predictors are compatible with missing at random after conditioning on those
#' predictors, while missing not at random cannot be confirmed or ruled out from
#' observed data alone.
#'
#' @param data Data frame containing the response and candidate predictors.
#' @param response_var Response column name.
#' @param predictors Optional predictor column names. Defaults to observed
#'   columns other than the response, subject identifier, and stored simulation
#'   truth columns.
#' @param time_var,subject_var Optional time and subject identifier column
#'   names. `time_var` is included as a predictor when present; `subject_var` is
#'   excluded by default.
#' @param alpha Significance threshold used to flag predictor associations.
#' @param max_factor_levels Maximum number of levels allowed for default factor
#'   predictors.
#'
#' @return An object of class `gamlss_longitudinal_missingness_check`. The
#'   `terms` element contains likelihood-ratio screening rows with `term`,
#'   `statistic`, `df`, `p_value`, and `method`. The default method is a
#'   multivariable likelihood-ratio test; terms that are aliased in the
#'   multivariable model fall back to a single-predictor likelihood-ratio screen.
#' @export
check_missingness <- function(
    data,
    response_var,
    predictors = NULL,
    time_var = NULL,
    subject_var = NULL,
    alpha = 0.05,
    max_factor_levels = 20) {
  args <- .missingness_validate_args(
    data = data,
    response_var = response_var,
    time_var = time_var,
    subject_var = subject_var,
    alpha = alpha,
    max_factor_levels = max_factor_levels
  )
  data <- args$data
  alpha <- args$alpha
  max_factor_levels <- args$max_factor_levels

  response <- .missingness_response_summary(data, response_var)
  missing_response <- response$missing_response
  response_summary <- response$response_summary

  predictors <- .missingness_resolve_predictors(
    data = data,
    response_var = response_var,
    predictors = predictors,
    time_var = time_var,
    subject_var = subject_var
  )

  predictor_summary <- .missingness_predictor_summary(data, predictors, max_factor_levels)
  predictors_used <- predictor_summary$predictor[predictor_summary$used]

  early <- .missingness_early_result(
    response_summary = response_summary,
    predictor_summary = predictor_summary,
    alpha = alpha
  )
  if (!is.null(early)) {
    return(early)
  }

  model_data <- .missingness_model_data(
    data = data,
    response_var = response_var,
    predictors_used = predictors_used,
    missing_response = missing_response
  )
  fit <- .missingness_fit_model(model_data, predictors_used)
  if (is.null(fit)) {
    return(.missingness_result(
      response_summary = response_summary,
      predictor_summary = predictor_summary,
      term_tests = data.frame(),
      model = NULL,
      alpha = alpha,
      assessment = "not_estimable",
      message = "The missingness model was not estimable after removing rows with missing predictors."
    ))
  }
  if (inherits(fit, "error")) {
    return(.missingness_result(
      response_summary = response_summary,
      predictor_summary = predictor_summary,
      term_tests = data.frame(),
      model = NULL,
      alpha = alpha,
      assessment = "not_estimable",
      message = paste("The missingness model failed:", conditionMessage(fit))
    ))
  }

  term_tests <- .missingness_term_tests(fit)
  assessment <- .missingness_assessment(term_tests, alpha)

  .missingness_result(
    response_summary = response_summary,
    predictor_summary = predictor_summary,
    term_tests = term_tests,
    model = fit,
    alpha = alpha,
    assessment = assessment$assessment,
    message = assessment$message
  )
}
