.missingness_check_column <- function(data, column, arg) {
  if (!is.character(column) || length(column) != 1L || is.na(column) || !nzchar(column)) {
    stop("'", arg, "' must be a single column name.", call. = FALSE)
  }
  if (!column %in% names(data)) {
    stop(arg, "='", column, "' not found in 'data'.", call. = FALSE)
  }
  invisible(TRUE)
}

.missingness_formula <- function(response, predictors) {
  terms <- vapply(c(response, predictors), .missingness_formula_name, character(1))
  stats::as.formula(
    paste(terms[[1L]], "~", paste(terms[-1L], collapse = " + ")),
    env = parent.frame()
  )
}

.missingness_formula_name <- function(name) {
  reserved <- c(
    "if", "else", "repeat", "while", "function", "for", "in", "next", "break",
    "TRUE", "FALSE", "NULL", "Inf", "NaN", "NA", "NA_integer_", "NA_real_",
    "NA_complex_", "NA_character_"
  )
  if (identical(make.names(name), name) && !name %in% reserved) {
    return(name)
  }
  paste0("`", gsub("`", "\\\\`", name), "`")
}

.missingness_term_label <- function(term) {
  if (grepl("^`.*`$", term)) {
    term <- substring(term, 2L, nchar(term) - 1L)
    term <- gsub("\\\\`", "`", term)
  }
  term
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
    term <- unname(vapply(rownames(tab), .missingness_term_label, character(1)))
    keep_term <- term != "<none>"
    tab <- tab[keep_term, , drop = FALSE]
    term <- term[keep_term]
    p_col <- intersect(c("Pr(>Chi)", "Pr(>Chisq)", "Pr(>LRT)"), names(tab))
    if (length(p_col) > 0L) {
      tests <- data.frame(
        term = term,
        statistic = as.numeric(tab$LRT %||% tab$Chisq %||% NA_real_),
        df = as.numeric(tab$Df %||% NA_real_),
        p_value = as.numeric(tab[[p_col[[1L]]]]),
        method = "multivariable_lrt",
        stringsAsFactors = FALSE
      )
      not_estimable <- !is.finite(tests$p_value) |
        !is.finite(tests$statistic) |
        !is.finite(tests$df) |
        tests$df <= 0
      if (any(not_estimable)) {
        fallback <- .missingness_single_term_tests(fit, tests$term[not_estimable])
        idx <- match(fallback$term, tests$term)
        tests[idx, c("statistic", "df", "p_value", "method")] <-
          fallback[c("statistic", "df", "p_value", "method")]
      }
      rownames(tests) <- NULL
      return(tests)
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
    method = "wald_coefficient",
    stringsAsFactors = FALSE
  )
}

.missingness_single_term_tests <- function(fit, terms) {
  terms <- unname(terms)
  model_data <- stats::model.frame(fit)
  rows <- lapply(terms, function(term) {
    .missingness_single_term_test(model_data, term)
  })
  out <- do.call(rbind, rows)
  rownames(out) <- NULL
  out
}

.missingness_single_term_test <- function(model_data, term) {
  not_estimable <- function() {
    data.frame(
      term = term,
      statistic = NA_real_,
      df = NA_real_,
      p_value = NA_real_,
      method = "not_estimable",
      stringsAsFactors = FALSE
    )
  }

  if (!term %in% names(model_data)) {
    return(not_estimable())
  }

  term_data <- model_data[c(".response_missing", term)]
  term_data <- term_data[stats::complete.cases(term_data), , drop = FALSE]
  term_data <- droplevels(term_data)
  if (
    nrow(term_data) < 3L ||
      length(unique(term_data$.response_missing)) < 2L ||
      length(unique(term_data[[term]])) < 2L
  ) {
    return(not_estimable())
  }

  null_fit <- tryCatch(
    suppressWarnings(stats::glm(.response_missing ~ 1, data = term_data, family = stats::binomial())),
    error = function(e) NULL
  )
  term_fit <- tryCatch(
    suppressWarnings(stats::glm(
      .missingness_formula(".response_missing", term),
      data = term_data,
      family = stats::binomial()
    )),
    error = function(e) NULL
  )
  if (is.null(null_fit) || is.null(term_fit)) {
    return(not_estimable())
  }

  statistic <- stats::deviance(null_fit) - stats::deviance(term_fit)
  df <- stats::df.residual(null_fit) - stats::df.residual(term_fit)
  if (!is.finite(statistic) || !is.finite(df) || df <= 0) {
    return(not_estimable())
  }

  statistic <- max(0, statistic)
  data.frame(
    term = term,
    statistic = statistic,
    df = df,
    p_value = stats::pchisq(statistic, df = df, lower.tail = FALSE),
    method = "univariable_lrt",
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
