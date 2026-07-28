.missingness_validate_args <- function(data, response_var, time_var, subject_var, alpha, max_factor_levels) {
  if (!is.data.frame(data)) {
    data <- as.data.frame(data, stringsAsFactors = FALSE)
  } else {
    data <- as.data.frame(data, stringsAsFactors = FALSE)
  }
  .missingness_check_column(data, response_var, "response_var")
  if (!is.null(time_var)) .missingness_check_column(data, time_var, "time_var")
  if (!is.null(subject_var)) .missingness_check_column(data, subject_var, "subject_var")

  alpha <- as.numeric(alpha)
  if (length(alpha) != 1L || !is.finite(alpha) || alpha <= 0 || alpha >= 1) {
    stop("'alpha' must be a single number between 0 and 1.", call. = FALSE)
  }
  max_factor_levels <- as.integer(max_factor_levels)
  if (length(max_factor_levels) != 1L || !is.finite(max_factor_levels) || max_factor_levels < 2L) {
    stop("'max_factor_levels' must be an integer of at least 2.", call. = FALSE)
  }

  list(data = data, alpha = alpha, max_factor_levels = max_factor_levels)
}

.missingness_resolve_predictors <- function(data, response_var, predictors, time_var, subject_var) {
  if (is.null(predictors)) {
    excluded <- unique(c(response_var, subject_var, "u", names(data)[grepl("^true_", names(data))]))
    predictors <- setdiff(names(data), excluded)
  } else {
    if (!is.character(predictors) || any(is.na(predictors)) || any(!nzchar(predictors))) {
      stop("'predictors' must be NULL or a character vector of column names.", call. = FALSE)
    }
    missing_predictors <- setdiff(predictors, names(data))
    if (length(missing_predictors) > 0L) {
      stop("Predictor column(s) not found: ", paste(missing_predictors, collapse = ", "), call. = FALSE)
    }
  }
  if (!is.null(time_var) && time_var %in% names(data) && !time_var %in% predictors) {
    predictors <- c(time_var, predictors)
  }
  unique(setdiff(predictors, response_var))
}

.missingness_response_summary <- function(data, response_var) {
  missing_response <- is.na(data[[response_var]])
  list(
    missing_response = missing_response,
    response_summary = data.frame(
      n = nrow(data),
      n_missing = sum(missing_response),
      prop_missing = mean(missing_response),
      stringsAsFactors = FALSE
    )
  )
}

.missingness_early_result <- function(response_summary, predictor_summary, alpha) {
  if (response_summary$n_missing == 0L) {
    return(.missingness_result(
      response_summary = response_summary,
      predictor_summary = predictor_summary,
      term_tests = data.frame(),
      model = NULL,
      alpha = alpha,
      assessment = "no_missing_responses",
      message = "No response missingness was detected."
    ))
  }
  if (response_summary$n_missing == response_summary$n) {
    return(.missingness_result(
      response_summary = response_summary,
      predictor_summary = predictor_summary,
      term_tests = data.frame(),
      model = NULL,
      alpha = alpha,
      assessment = "all_responses_missing",
      message = "All responses are missing, so missingness cannot be modelled."
    ))
  }
  if (!any(predictor_summary$used)) {
    return(.missingness_result(
      response_summary = response_summary,
      predictor_summary = predictor_summary,
      term_tests = data.frame(),
      model = NULL,
      alpha = alpha,
      assessment = "no_usable_predictors",
      message = "No usable observed predictors were available for the missingness model."
    ))
  }
  NULL
}

.missingness_model_data <- function(data, response_var, predictors_used, missing_response) {
  model_data <- data[c(response_var, predictors_used)]
  model_data$.response_missing <- missing_response
  keep <- stats::complete.cases(model_data[predictors_used])
  model_data[keep, , drop = FALSE]
}

.missingness_fit_model <- function(model_data, predictors_used) {
  if (nrow(model_data) < 3L || length(unique(model_data$.response_missing)) < 2L) {
    return(NULL)
  }
  form <- stats::reformulate(predictors_used, response = ".response_missing")
  tryCatch(
    stats::glm(form, data = model_data, family = stats::binomial()),
    error = function(e) e
  )
}

.missingness_assessment <- function(term_tests, alpha) {
  flagged <- term_tests[is.finite(term_tests$p_value) & term_tests$p_value < alpha, , drop = FALSE]
  if (nrow(flagged) > 0L) {
    return(list(
      assessment = "covariate_related_missingness",
      message = paste(
        "Response missingness is associated with observed predictor(s):",
        paste(flagged$term, collapse = ", "),
        ". This is compatible with MAR conditional on the observed predictors,",
        "but it does not rule out MNAR."
      )
    ))
  }
  list(
    assessment = "no_detected_covariate_association",
    message = paste(
      "No observed predictor association was detected at alpha =",
      format(alpha),
      ". This is compatible with MCAR, but MNAR cannot be ruled out from observed data alone."
    )
  )
}
