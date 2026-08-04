#' Broom-style table methods for fitted longitudinal GAMLSS models
#'
#' These methods expose fitted `gamlss.longitudinal` objects in the table shapes
#' expected by reporting packages. `tidy()` returns fixed-coefficient summaries,
#' `glance()` returns one-row model-fit summaries, and `augment()` returns
#' row-level fitted values and residuals.
#'
#' @param x A fitted `gamlss.longitudinal` object.
#' @param conf.int Logical; include Wald confidence intervals for coefficients.
#' @param conf.level Confidence level for intervals.
#' @param parameter Optional distributional parameter(s) to retain.
#' @param include_vcov Logical; compute or reuse variance-covariance output.
#' @param vcov_method Variance-covariance method passed to
#'   [summary.gamlss.longitudinal()].
#' @param data Optional data frame to augment. Defaults to the fitted expanded
#'   model frame.
#' @param newdata Optional new data to augment. When supplied, residuals are only
#'   added if a response column is present.
#' @param type Prediction type passed to [predict.gamlss.longitudinal()].
#' @param residual_type Residual type used when augmenting fitted rows.
#' @param se.fit Logical; include prediction standard errors when supported.
#' @param interval Interval type passed to [predict.gamlss.longitudinal()].
#' @param level Confidence level for prediction intervals.
#' @param ... Additional arguments passed to downstream methods.
#'
#' @return A data frame.
#' @name gamlss_longitudinal_table_methods
NULL

#' @rdname gamlss_longitudinal_table_methods
#' @export
tidy.gamlss.longitudinal <- function(
    x,
    conf.int = TRUE,
    conf.level = 0.95,
    parameter = NULL,
    include_vcov = TRUE,
    vcov_method = "analytical",
    ...) {
  if (!inherits(x, "gamlss.longitudinal")) {
    stop("'x' must be a fitted 'gamlss.longitudinal' object.", call. = FALSE)
  }

  s <- summary(
    x,
    include_vcov = include_vcov,
    ci_level = conf.level,
    vcov_method = vcov_method,
    ...
  )

  coef_tbl <- as.data.frame(s$coefficients, stringsAsFactors = FALSE)
  if (nrow(coef_tbl) == 0L) {
    return(.gl_empty_tidy_table(conf.int = conf.int))
  }

  split_terms <- .gl_split_coefficient_terms(coef_tbl$term)
  out <- data.frame(
    component = ifelse(coef_tbl$parameter %in% c("theta", "zeta"), "copula", "margin"),
    parameter = coef_tbl$parameter,
    term = split_terms$term,
    coefficient = coef_tbl$term,
    estimate = coef_tbl$estimate,
    std.error = coef_tbl$std_error,
    statistic = coef_tbl$estimate / coef_tbl$std_error,
    p.value = coef_tbl$p_value,
    signif = coef_tbl$signif,
    stringsAsFactors = FALSE
  )

  if (!is.null(parameter)) {
    parameter <- as.character(parameter)
    out <- out[out$parameter %in% parameter, , drop = FALSE]
  }

  if (isTRUE(conf.int)) {
    z <- stats::qnorm((1 + conf.level) / 2)
    out$conf.low <- out$estimate - z * out$std.error
    out$conf.high <- out$estimate + z * out$std.error
  }

  rownames(out) <- NULL
  out
}

#' @rdname gamlss_longitudinal_table_methods
#' @export
glance.gamlss.longitudinal <- function(
    x,
    include_vcov = FALSE,
    vcov_method = "analytical",
    ...) {
  if (!inherits(x, "gamlss.longitudinal")) {
    stop("'x' must be a fitted 'gamlss.longitudinal' object.", call. = FALSE)
  }

  s <- summary(
    x,
    include_vcov = include_vcov,
    vcov_method = vcov_method,
    ...
  )

  data.frame(
    margin_dist = s$model$margin_dist,
    copula_dist = s$model$copula_dist,
    nobs = .gl_safe_nobs(x, "observed"),
    nobs_expanded = .gl_safe_nobs(x, "expanded"),
    n_subjects = s$model$n_subjects,
    n_timepoints = s$model$n_timepoints,
    n_fixed = s$model$n_fixed,
    n_smooth_terms = s$model$n_smooth_terms,
    edf_smooth = s$model$edf_smooth,
    logLik = s$fit$logLik,
    AIC = s$fit$AIC,
    BIC = s$fit$BIC,
    vcov_method = s$fit$vcov_method,
    converged = x$convergence$converged %||% NA,
    outer_iterations = x$convergence$outer_iterations %||% NA_integer_,
    stringsAsFactors = FALSE
  )
}

#' @rdname gamlss_longitudinal_table_methods
#' @export
augment.gamlss.longitudinal <- function(
    x,
    data = NULL,
    newdata = NULL,
    type = c("mean", "mu", "median"),
    residual_type = c("response", "pearson", "quantile"),
    se.fit = FALSE,
    interval = c("none", "confidence"),
    level = 0.95,
    vcov_method = "analytical",
    ...) {
  if (!inherits(x, "gamlss.longitudinal")) {
    stop("'x' must be a fitted 'gamlss.longitudinal' object.", call. = FALSE)
  }

  type <- match.arg(type)
  residual_type <- match.arg(residual_type)
  interval <- match.arg(interval)

  supplied_newdata <- !is.null(newdata)
  out <- if (supplied_newdata) {
    as.data.frame(newdata, stringsAsFactors = FALSE)
  } else if (!is.null(data)) {
    as.data.frame(data, stringsAsFactors = FALSE)
  } else {
    stats::model.frame(x, type = "expanded")
  }

  predict_newdata <- if (supplied_newdata || !is.null(data)) out else NULL
  pred <- predict(
    x,
    newdata = predict_newdata,
    type = type,
    se.fit = se.fit,
    interval = interval,
    level = level,
    vcov_method = vcov_method,
    ...
  )

  if (is.data.frame(pred)) {
    out$.fitted <- .gl_first_nonnull_column(pred, c("fit", type, "mu", "mean", "median"))
    if ("se.fit" %in% names(pred)) out$.se.fit <- pred$se.fit
    if ("conf.low" %in% names(pred)) out$.lower <- pred$conf.low
    if ("conf.high" %in% names(pred)) out$.upper <- pred$conf.high
  } else {
    out$.fitted <- as.numeric(pred)
  }

  response_col <- x$response_var %||% "response"
  if (response_col %in% names(out)) {
    out$.resid <- as.numeric(out[[response_col]]) - out$.fitted
  } else if (!supplied_newdata && is.null(data)) {
    out$.resid <- stats::residuals(x, type = residual_type, finite = FALSE)
  }

  rownames(out) <- NULL
  out
}

#' Build a publication-ready table from a fitted model
#'
#' `publication_table()` formats the model summaries most often needed in
#' papers and reports. Use `tidy()`, `glance()`, or `reporting_table()` when you
#' need numeric columns for further calculation; use `publication_table()` when
#' the next step is rendering.
#'
#' @param object A fitted `gamlss.longitudinal` object.
#' @param table Table to build: fixed coefficients, model fit statistics, or
#'   grouped predictions.
#' @param newdata Data used for prediction tables. Defaults to the fitted
#'   observed model frame when omitted.
#' @param by Optional grouping columns for prediction tables.
#' @param probs Quantiles to include in prediction tables.
#' @param threshold Optional threshold for prediction probabilities.
#' @param direction Probability direction when `threshold` is supplied.
#' @param conf.level Confidence level for coefficient intervals.
#' @param include_vcov Logical; include coefficient uncertainty when possible.
#' @param vcov_method Variance-covariance method passed to summary methods.
#' @param digits Number of decimal places for estimates and intervals.
#' @param p_digits Number of decimal places for p-values.
#' @param output Output format. `"data.frame"` returns a formatted data frame;
#'   `"latex"` returns LaTeX source, and `"kable"`, `"gt"`, and `"flextable"`
#'   use optional reporting packages.
#' @param caption Optional table caption for rendered outputs.
#' @param ... Additional arguments passed to the underlying table builder.
#'
#' @return A formatted data frame, LaTeX `knitr_kable`, `knitr_kable`,
#'   `gt_tbl`, or `flextable` object depending on `output`.
#' @export
publication_table <- function(
    object,
    table = c("coefficients", "model", "predictions"),
    newdata = NULL,
    by = NULL,
    probs = c(0.1, 0.5, 0.9),
    threshold = NULL,
    direction = c("above", "below"),
    conf.level = 0.95,
    include_vcov = TRUE,
    vcov_method = "analytical",
    digits = 3,
    p_digits = 3,
    output = c("data.frame", "latex", "kable", "gt", "flextable"),
    caption = NULL,
    ...) {
  if (!inherits(object, "gamlss.longitudinal")) {
    stop("'object' must be a fitted 'gamlss.longitudinal' object.", call. = FALSE)
  }

  table <- match.arg(table)
  output <- match.arg(output)
  direction <- match.arg(direction)

  out <- switch(table,
    coefficients = {
      raw <- tidy.gamlss.longitudinal(
        object,
        conf.int = TRUE,
        conf.level = conf.level,
        include_vcov = include_vcov,
        vcov_method = vcov_method,
        ...
      )
      .gl_publication_coefficients(raw, conf.level, digits, p_digits)
    },
    model = {
      raw <- glance.gamlss.longitudinal(
        object,
        include_vcov = include_vcov,
        vcov_method = vcov_method,
        ...
      )
      .gl_publication_model(raw, digits)
    },
    predictions = {
      if (is.null(newdata)) {
        newdata <- stats::model.frame(object, type = "observed")
      }
      raw <- reporting_table(
        object,
        newdata = newdata,
        by = by,
        probs = probs,
        threshold = threshold,
        direction = direction
      )
      .gl_publication_predictions(raw, digits)
    }
  )

  .gl_render_publication_table(out, output = output, caption = caption)
}

#' @export
print.gamlss_longitudinal_publication_table <- function(x, ...) {
  print.data.frame(x, row.names = FALSE, ...)
  invisible(x)
}

.gl_empty_tidy_table <- function(conf.int = TRUE) {
  out <- data.frame(
    component = character(0),
    parameter = character(0),
    term = character(0),
    coefficient = character(0),
    estimate = numeric(0),
    std.error = numeric(0),
    statistic = numeric(0),
    p.value = numeric(0),
    signif = character(0),
    stringsAsFactors = FALSE
  )
  if (isTRUE(conf.int)) {
    out$conf.low <- numeric(0)
    out$conf.high <- numeric(0)
  }
  out
}

.gl_split_coefficient_terms <- function(terms) {
  parameter <- sub("\\..*$", "", terms)
  term <- sub("^[^.]+\\.", "", terms)
  term[term %in% c("intercept", "(Intercept)")] <- "(Intercept)"
  data.frame(parameter = parameter, term = term, stringsAsFactors = FALSE)
}

.gl_safe_nobs <- function(object, type) {
  tryCatch(stats::nobs(object, type = type), error = function(e) NA_integer_)
}

.gl_first_nonnull_column <- function(data, candidates) {
  for (candidate in candidates) {
    if (candidate %in% names(data)) {
      return(as.numeric(data[[candidate]]))
    }
  }
  as.numeric(data[[ncol(data)]])
}

.gl_publication_coefficients <- function(raw, conf.level, digits, p_digits) {
  ci_col <- paste0(round(conf.level * 100), "% CI")
  out <- data.frame(
    Component = .gl_title_case(raw$component),
    Parameter = raw$parameter,
    Term = raw$term,
    Estimate = .gl_format_number(raw$estimate, digits),
    `Std. Error` = .gl_format_number(raw$std.error, digits),
    check.names = FALSE,
    stringsAsFactors = FALSE
  )
  out[[ci_col]] <- .gl_format_interval(raw$conf.low, raw$conf.high, digits)
  out[["p-value"]] <- .gl_format_p_value(raw$p.value, p_digits)
  out[["Significance"]] <- ifelse(is.na(raw$signif), "", raw$signif)
  out
}

.gl_publication_model <- function(raw, digits) {
  data.frame(
    Statistic = c(
      "Margin distribution",
      "Copula distribution",
      "Observed rows",
      "Subjects",
      "Time points",
      "Fixed coefficients",
      "Smooth terms",
      "Smooth EDF",
      "Log likelihood",
      "AIC",
      "BIC",
      "VCOV method",
      "Converged",
      "Outer iterations"
    ),
    Value = c(
      raw$margin_dist,
      raw$copula_dist,
      .gl_format_count(raw$nobs),
      .gl_format_count(raw$n_subjects),
      .gl_format_count(raw$n_timepoints),
      .gl_format_count(raw$n_fixed),
      .gl_format_count(raw$n_smooth_terms),
      .gl_format_number(raw$edf_smooth, digits),
      .gl_format_number(raw$logLik, digits),
      .gl_format_number(raw$AIC, digits),
      .gl_format_number(raw$BIC, digits),
      raw$vcov_method,
      as.character(raw$converged),
      .gl_format_count(raw$outer_iterations)
    ),
    stringsAsFactors = FALSE
  )
}

.gl_publication_predictions <- function(raw, digits) {
  out <- as.data.frame(raw, stringsAsFactors = FALSE)
  for (nm in names(out)) {
    if (is.numeric(out[[nm]]) && !identical(nm, "n")) {
      out[[nm]] <- .gl_format_number(out[[nm]], digits)
    }
  }
  names(out) <- .gl_publication_prediction_names(names(out))
  out
}

.gl_publication_prediction_names <- function(nms) {
  vapply(nms, function(nm) {
    if (identical(nm, "n")) return("N")
    if (grepl("^q[0-9]+$", nm)) return(toupper(nm))
    if (grepl("^prob_(above|below)_", nm)) {
      parts <- strsplit(nm, "_", fixed = TRUE)[[1L]]
      return(paste0("Pr(", parts[2L], " ", paste(parts[-c(1L, 2L)], collapse = "_"), ")"))
    }
    .gl_title_case(gsub("_", " ", nm, fixed = TRUE))
  }, character(1))
}

.gl_render_publication_table <- function(out, output, caption = NULL) {
  if (identical(output, "data.frame")) {
    attr(out, "caption") <- caption
    class(out) <- c("gamlss_longitudinal_publication_table", class(out))
    return(out)
  }

  if (identical(output, "kable")) {
    .gl_require_reporting_package("knitr", output)
    return(knitr::kable(out, caption = caption, booktabs = TRUE))
  }

  if (identical(output, "latex")) {
    .gl_require_reporting_package("knitr", output)
    return(knitr::kable(out, format = "latex", caption = caption, booktabs = TRUE))
  }

  if (identical(output, "gt")) {
    .gl_require_reporting_package("gt", output)
    gt_out <- gt::gt(out)
    if (!is.null(caption)) {
      gt_out <- gt::tab_header(gt_out, title = caption)
    }
    return(gt_out)
  }

  .gl_require_reporting_package("flextable", output)
  ft <- flextable::flextable(out)
  if (!is.null(caption)) {
    ft <- flextable::set_caption(ft, caption = caption)
  }
  flextable::autofit(ft)
}

.gl_require_reporting_package <- function(package, output) {
  if (!requireNamespace(package, quietly = TRUE)) {
    stop(
      "Package '", package, "' is required for output = '", output, "'. ",
      "Use output = 'data.frame' or install '", package, "'.",
      call. = FALSE
    )
  }
}

.gl_format_number <- function(x, digits) {
  ifelse(
    is.na(x),
    "",
    formatC(as.numeric(x), format = "f", digits = digits)
  )
}

.gl_format_count <- function(x) {
  ifelse(is.na(x), "", formatC(as.integer(x), format = "d"))
}

.gl_format_p_value <- function(x, digits) {
  threshold <- 10^(-digits)
  ifelse(
    is.na(x),
    "",
    ifelse(
      x > 0 & x < threshold,
      paste0("<", formatC(threshold, format = "f", digits = digits)),
      formatC(x, format = "f", digits = digits)
    )
  )
}

.gl_format_interval <- function(lower, upper, digits) {
  ifelse(
    is.na(lower) | is.na(upper),
    "",
    paste0(
      "[",
      .gl_format_number(lower, digits),
      ", ",
      .gl_format_number(upper, digits),
      "]"
    )
  )
}

.gl_title_case <- function(x) {
  x <- as.character(x)
  paste0(toupper(substr(x, 1L, 1L)), substr(x, 2L, nchar(x)))
}
