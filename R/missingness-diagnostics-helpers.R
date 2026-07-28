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
