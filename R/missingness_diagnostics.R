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
#' @return An object of class `gamlss_longitudinal_missingness_check`.
#' @export
check_missingness <- function(
  data,
  response_var,
  predictors = NULL,
  time_var = NULL,
  subject_var = NULL,
  alpha = 0.05,
  max_factor_levels = 20
) {
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

  missing_response <- is.na(data[[response_var]])
  response_summary <- data.frame(
    n = nrow(data),
    n_missing = sum(missing_response),
    prop_missing = mean(missing_response),
    stringsAsFactors = FALSE
  )

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
  predictors <- unique(setdiff(predictors, response_var))

  predictor_summary <- .missingness_predictor_summary(data, predictors, max_factor_levels)
  predictors_used <- predictor_summary$predictor[predictor_summary$used]

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
  if (length(predictors_used) == 0L) {
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

  model_data <- data[c(response_var, predictors_used)]
  model_data$.response_missing <- missing_response
  keep <- stats::complete.cases(model_data[predictors_used])
  model_data <- model_data[keep, , drop = FALSE]
  if (nrow(model_data) < 3L || length(unique(model_data$.response_missing)) < 2L) {
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

  form <- stats::reformulate(predictors_used, response = ".response_missing")
  fit <- tryCatch(
    stats::glm(form, data = model_data, family = stats::binomial()),
    error = function(e) e
  )
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
  flagged <- term_tests[is.finite(term_tests$p_value) & term_tests$p_value < alpha, , drop = FALSE]
  if (nrow(flagged) > 0L) {
    assessment <- "covariate_related_missingness"
    message <- paste(
      "Response missingness is associated with observed predictor(s):",
      paste(flagged$term, collapse = ", "),
      ". This is compatible with MAR conditional on the observed predictors,",
      "but it does not rule out MNAR."
    )
  } else {
    assessment <- "no_detected_covariate_association"
    message <- paste(
      "No observed predictor association was detected at alpha =",
      format(alpha),
      ". This is compatible with MCAR, but MNAR cannot be ruled out from observed data alone."
    )
  }

  .missingness_result(
    response_summary = response_summary,
    predictor_summary = predictor_summary,
    term_tests = term_tests,
    model = fit,
    alpha = alpha,
    assessment = assessment,
    message = message
  )
}

.missingness_check_column <- function(data, column, arg) {
  if (!is.character(column) || length(column) != 1L || is.na(column) || !nzchar(column)) {
    stop("'", arg, "' must be a single column name.", call. = FALSE)
  }
  if (!column %in% names(data)) {
    stop(arg, "='", column, "' not found in 'data'.", call. = FALSE)
  }
  invisible(TRUE)
}

.missingness_predictor_summary <- function(data, predictors, max_factor_levels) {
  if (length(predictors) == 0L) {
    return(data.frame(
      predictor = character(0),
      class = character(0),
      n_observed = integer(0),
      n_unique = integer(0),
      used = logical(0),
      reason = character(0),
      stringsAsFactors = FALSE
    ))
  }
  rows <- lapply(predictors, function(predictor) {
    x <- data[[predictor]]
    n_observed <- sum(!is.na(x))
    n_unique <- length(unique(x[!is.na(x)]))
    cls <- class(x)[1L]
    used <- TRUE
    reason <- "used"
    if (n_observed == 0L) {
      used <- FALSE
      reason <- "all values missing"
    } else if (n_unique < 2L) {
      used <- FALSE
      reason <- "fewer than two observed values"
    } else if ((is.factor(x) || is.character(x)) && n_unique > max_factor_levels) {
      used <- FALSE
      reason <- paste0("more than ", max_factor_levels, " levels")
    }
    data.frame(
      predictor = predictor,
      class = cls,
      n_observed = n_observed,
      n_unique = n_unique,
      used = used,
      reason = reason,
      stringsAsFactors = FALSE
    )
  })
  do.call(rbind, rows)
}

.missingness_term_tests <- function(fit) {
  out <- tryCatch(
    stats::drop1(fit, test = "Chisq"),
    error = function(e) NULL
  )
  if (!is.null(out)) {
    tab <- as.data.frame(out)
    term <- rownames(tab)
    tab <- tab[term != "<none>", , drop = FALSE]
    term <- term[term != "<none>"]
    p_col <- intersect(c("Pr(>Chi)", "Pr(>Chisq)", "Pr(>LRT)"), names(tab))
    if (length(p_col) > 0L) {
      return(data.frame(
        term = term,
        statistic = as.numeric(tab$LRT %||% tab$Chisq %||% NA_real_),
        df = as.numeric(tab$Df %||% NA_real_),
        p_value = as.numeric(tab[[p_col[[1L]]]]),
        stringsAsFactors = FALSE
      ))
    }
  }

  coef_tab <- as.data.frame(summary(fit)$coefficients)
  coef_tab$term <- rownames(coef_tab)
  coef_tab <- coef_tab[coef_tab$term != "(Intercept)", , drop = FALSE]
  data.frame(
    term = coef_tab$term,
    statistic = as.numeric(coef_tab$`z value`),
    df = NA_real_,
    p_value = as.numeric(coef_tab$`Pr(>|z|)`),
    stringsAsFactors = FALSE
  )
}

.missingness_result <- function(response_summary, predictor_summary, term_tests, model, alpha, assessment, message) {
  out <- list(
    response = response_summary,
    predictors = predictor_summary,
    terms = term_tests,
    model = model,
    alpha = alpha,
    assessment = assessment,
    message = message
  )
  class(out) <- "gamlss_longitudinal_missingness_check"
  out
}

#' @export
print.gamlss_longitudinal_missingness_check <- function(x, digits = max(3, getOption("digits") - 3), ...) {
  cat("\nResponse Missingness Check\n")
  cat("--------------------------\n")
  print(x$response, digits = digits, row.names = FALSE)
  cat("\nAssessment:", x$assessment, "\n")
  cat(x$message, "\n", sep = "")
  if (!is.null(x$terms) && nrow(x$terms) > 0L) {
    cat("\nPredictor tests:\n")
    print(x$terms, digits = digits, row.names = FALSE)
  }
  invisible(x)
}
