# Consolidated support module for reviewer navigation.
# Source blocks are grouped mechanically from the files named below.

# ---- benchmark-comparator.R ----

#' Fit standard longitudinal comparator models for adoption benchmarks

#'

#' `benchmark_standard_models()` is an opt-in scaffold for comparing

#' `gamlss.longitudinal` with models users already know. It is intentionally

#' narrow: optionally score an already-fitted primary model, fit common

#' mean-model baselines, record whether they ran, and return timing,

#' response-scale prediction metrics, and distributional diagnostics when the

#' fitted model family supports them.

#'

#' @param data Long-format data frame.

#' @param formula Mean-model formula used by the comparator models.

#' @param subject_var Subject identifier column.

#' @param family A base R family object or one of `"gaussian"`, `"poisson"`,

#'   `"binomial"`, or `"gamma"`.

#' @param comparators Character vector containing any of `"glm"`, `"gee"`,

#'   `"glmm"`, `"gam"`, `"gamm"`, and `"glmmTMB"`.

#' @param correlation Working correlation passed to `geepack::geeglm()`.

#' @param add_subject_re_to_gamm Logical; add `s(subject, bs = "re")` to the

#'   GAMM comparator.

#' @param add_subject_re_to_gam Deprecated alias for `add_subject_re_to_gamm`.

#' @param fit Optional already-fitted primary model to score beside the standard

#'   comparators. For a `gamlss.longitudinal` fit, response-mean predictions

#'   are obtained with `predict(fit, newdata = data, type = "mean")`; quantile,

#'   probability, PIT, density, and interval metrics still use the fitted

#'   distribution.

#' @param fit_name Label used for the supplied `fit` row.

#' @param distributional_metrics Logical; compute interval coverage, PIT

#'   calibration, tail-frequency diagnostics, and truth-aware quantile/tail

#'   metrics when enough information is available.

#' @param truth_family Optional `gamlss.dist` family name used to derive true

#'   quantiles from `true_*` columns. Defaults are inferred from `family`.

#' @param quantile_prob Probability used for the benchmark quantile and upper

#'   tail metrics.

#' @param interval_level Prediction interval level used for empirical coverage.

#' @param ... Additional arguments passed to each comparator fit.

#'

#' @return An object of class `gamlss_longitudinal_benchmark` with `results`,

#'   `model_status`, `fits`, `coefficients`, and `interpretation` components. The

#'   `coefficients` component contains `mu` coefficient estimates and

#'   uncertainty comparisons for fixed/parametric terms available across the

#'   fitted models.

#' @export

benchmark_standard_models <- function(
    data,
    formula,
    subject_var,
    family = "gaussian",
    comparators = c("glm", "gee", "glmm", "gam", "gamm"),
    correlation = "exchangeable",
    add_subject_re_to_gamm = TRUE,
    add_subject_re_to_gam = NULL,
    fit = NULL,
    fit_name = "gamlss.longitudinal",
    distributional_metrics = TRUE,
    truth_family = NULL,
    quantile_prob = 0.9,
    interval_level = 0.95,
    ...) {
  inputs <- .benchmark_standard_models_inputs(
    data = data,
    formula = formula,
    subject_var = subject_var,
    family = family,
    comparators = comparators,
    add_subject_re_to_gamm = add_subject_re_to_gamm,
    add_subject_re_to_gam = add_subject_re_to_gam,
    fit = fit,
    fit_name = fit_name,
    quantile_prob = quantile_prob,
    interval_level = interval_level
  )

  data <- inputs$data
  formula <- inputs$formula
  subject_var <- inputs$subject_var
  family <- inputs$family
  comparators <- inputs$comparators
  add_subject_re_to_gamm <- inputs$add_subject_re_to_gamm
  fit_name <- inputs$fit_name
  quantile_prob <- inputs$quantile_prob
  interval_level <- inputs$interval_level


  comparator_runs <- lapply(comparators, function(comparator) {
    .benchmark_fit_one(
      data = data,
      formula = formula,
      subject_var = subject_var,
      family = family,
      comparator = comparator,
      correlation = correlation,
      add_subject_re_to_gamm = add_subject_re_to_gamm,
      distributional_metrics = distributional_metrics,
      truth_family = truth_family,
      quantile_prob = quantile_prob,
      interval_level = interval_level,
      ...
    )
  })

  supplied_run <- .benchmark_supplied_fit_one(
    fit = fit,
    data = data,
    formula = formula,
    family = family,
    fit_name = fit_name,
    distributional_metrics = distributional_metrics,
    truth_family = truth_family,
    quantile_prob = quantile_prob,
    interval_level = interval_level
  )

  runs <- c(if (!is.null(supplied_run)) list(supplied_run), comparator_runs)

  .benchmark_standard_models_result(
    runs = runs,
    supplied_run = supplied_run,
    comparators = comparators,
    fit_name = fit_name,
    formula = formula,
    subject_var = subject_var,
    family = family,
    interval_level = interval_level
  )
}

# ---- benchmark-comparator-workflow.R ----

#' Normalize and validate benchmark_standard_models() inputs
#'
#' @noRd
.benchmark_standard_models_inputs <- function(
    data,
    formula,
    subject_var,
    family,
    comparators,
    add_subject_re_to_gamm,
    add_subject_re_to_gam,
    fit,
    fit_name,
    quantile_prob,
    interval_level) {
  data <- as.data.frame(data, stringsAsFactors = FALSE)

  if (!is.character(subject_var) || length(subject_var) != 1L || !subject_var %in% names(data)) {
    stop("'subject_var' must be a single column name in 'data'.", call. = FALSE)
  }

  if (!is.null(add_subject_re_to_gam)) {
    warning(
      "'add_subject_re_to_gam' is deprecated; use 'add_subject_re_to_gamm' instead.",
      call. = FALSE
    )
    add_subject_re_to_gamm <- add_subject_re_to_gam
  }

  formula <- stats::as.formula(formula)
  family <- .benchmark_family(family)
  quantile_prob <- as.numeric(quantile_prob)[1L]
  interval_level <- as.numeric(interval_level)[1L]

  if (!is.finite(quantile_prob) || quantile_prob <= 0 || quantile_prob >= 1) {
    stop("'quantile_prob' must be a probability strictly between 0 and 1.", call. = FALSE)
  }

  if (!is.finite(interval_level) || interval_level <= 0 || interval_level >= 1) {
    stop("'interval_level' must be a probability strictly between 0 and 1.", call. = FALSE)
  }

  comparators <- unique(as.character(comparators))
  valid <- benchmark_comparator_status()$comparator
  bad <- setdiff(comparators, valid)

  if (length(bad) > 0L) {
    stop("Unknown comparator(s): ", paste(bad, collapse = ", "), call. = FALSE)
  }

  if (!is.null(fit)) {
    if (!is.character(fit_name) || length(fit_name) != 1L || is.na(fit_name) || !nzchar(fit_name)) {
      stop("'fit_name' must be a single non-empty character value.", call. = FALSE)
    }

    if (fit_name %in% comparators) {
      stop("'fit_name' must not duplicate a standard comparator name.", call. = FALSE)
    }
  }

  list(
    data = data,
    formula = formula,
    subject_var = subject_var,
    family = family,
    comparators = comparators,
    add_subject_re_to_gamm = add_subject_re_to_gamm,
    fit_name = fit_name,
    quantile_prob = quantile_prob,
    interval_level = interval_level
  )
}


#' Assemble a benchmark_standard_models() return object
#'
#' @noRd
.benchmark_standard_models_result <- function(
    runs,
    supplied_run,
    comparators,
    fit_name,
    formula,
    subject_var,
    family,
    interval_level) {
  fit_names <- c(if (!is.null(supplied_run)) fit_name, comparators)
  results <- do.call(rbind, lapply(runs, `[[`, "row"))
  primary_method <- if (!is.null(supplied_run)) fit_name else NULL
  fits <- stats::setNames(lapply(runs, `[[`, "fit"), fit_names)
  coefficients <- .benchmark_coefficient_comparison(
    fits = fits,
    results = results,
    primary_method = primary_method,
    parameter = "mu",
    level = interval_level
  )

  out <- list(
    results = results,
    model_status = .benchmark_model_status_table(results),
    fits = fits,
    formula = formula,
    subject_var = subject_var,
    family = family$family,
    primary_method = primary_method,
    coefficients = coefficients,
    interpretation = .benchmark_interpretation(results, primary_method = primary_method)
  )

  class(out) <- "gamlss_longitudinal_benchmark"
  out
}

# ---- benchmark-comparator-fit-dispatch.R ----

.benchmark_run_comparator_fit <- function(
    data,
    formula,
    subject_var,
    family,
    comparator,
    correlation,
    add_subject_re_to_gamm,
    extra_args = list()) {
  if (identical(comparator, "glm")) {
    do.call(
      stats::glm,
      c(
        list(
          formula = formula,
          data = data,
          family = family
        ),
        extra_args
      )
    )
  } else if (identical(comparator, "gee")) {
    data$.benchmark_subject_id <- data[[subject_var]]
    do.call(
      geepack::geeglm,
      c(
        list(
          formula = formula,
          id = as.name(".benchmark_subject_id"),
          data = data,
          family = family,
          corstr = correlation
        ),
        extra_args
      )
    )
  } else if (identical(comparator, "glmm")) {
    f_re <- .benchmark_formula_with_random_intercept(formula, subject_var)
    if (identical(family$family, "gaussian")) {
      do.call(
        lme4::lmer,
        c(
          list(
            formula = f_re,
            data = data
          ),
          extra_args
        )
      )
    } else {
      do.call(
        lme4::glmer,
        c(
          list(
            formula = f_re,
            data = data,
            family = family
          ),
          extra_args
        )
      )
    }
  } else if (identical(comparator, "gam")) {
    do.call(
      mgcv::gam,
      c(
        list(
          formula = formula,
          data = data,
          family = family,
          method = "REML"
        ),
        extra_args
      )
    )
  } else if (identical(comparator, "gamm")) {
    f_gamm <- if (isTRUE(add_subject_re_to_gamm)) {
      .benchmark_formula_with_subject_re_smooth(formula, subject_var)
    } else {
      formula
    }

    do.call(
      mgcv::gam,
      c(
        list(
          formula = f_gamm,
          data = data,
          family = family,
          method = "REML"
        ),
        extra_args
      )
    )
  } else {
    f_tmb <- .benchmark_formula_with_random_intercept(formula, subject_var)
    glmmTMB_fit <- getExportedValue("glmmTMB", "glmmTMB")
    do.call(
      glmmTMB_fit,
      c(
        list(
          formula = f_tmb,
          data = data,
          family = family
        ),
        extra_args
      )
    )
  }
}

# ---- benchmark-comparator-fitters.R ----

.benchmark_fit_one <- function(
    data,
    formula,
    subject_var,
    family,
    comparator,
    correlation,
    add_subject_re_to_gamm,
    distributional_metrics = TRUE,
    truth_family = NULL,
    quantile_prob = 0.9,
    interval_level = 0.95,
    ...) {
  status <- benchmark_comparator_status()

  status_row <- status[match(comparator, status$comparator), , drop = FALSE]

  if (nrow(status_row) != 1L || is.na(status_row$comparator)) {
    stop("Unknown comparator: ", comparator, call. = FALSE)
  }


  start <- Sys.time()

  warnings <- character(0)

  fit <- NULL

  error <- NULL


  if (!isTRUE(status_row$available)) {
    elapsed <- as.numeric(difftime(Sys.time(), start, units = "secs"))

    return(list(
      fit = NULL,
      row = .benchmark_empty_result_row(
        status_row = status_row,
        comparator = comparator,
        available = FALSE,
        success = FALSE,
        elapsed_sec = elapsed,
        nobs = nrow(data),
        error = paste0("Package '", status_row$package, "' is not installed.")
      )
    ))
  }


  smooth_comparators <- c("gam", "gamm")

  if (!comparator %in% smooth_comparators && .benchmark_formula_has_smooth(formula)) {
    elapsed <- as.numeric(difftime(Sys.time(), start, units = "secs"))

    return(list(
      fit = NULL,
      row = .benchmark_empty_result_row(
        status_row = status_row,
        comparator = comparator,
        available = TRUE,
        success = FALSE,
        elapsed_sec = elapsed,
        nobs = nrow(data),
        error = paste0(
          "Smooth terms via s(...) are only supported for 'gam' and 'gamm' comparators; ",
          "comparator '", comparator, "' cannot fit them."
        )
      )
    ))
  }


  captured <- withCallingHandlers(
    tryCatch(
      {
        .benchmark_run_comparator_fit(
          data = data,
          formula = formula,
          subject_var = subject_var,
          family = family,
          comparator = comparator,
          correlation = correlation,
          add_subject_re_to_gamm = add_subject_re_to_gamm,
          extra_args = list(...)
        )
      },
      error = function(e) {
        error <<- conditionMessage(e)

        NULL
      }
    ),
    warning = function(w) {
      warnings <<- c(warnings, conditionMessage(w))

      invokeRestart("muffleWarning")
    }
  )

  fit <- captured

  elapsed <- as.numeric(difftime(Sys.time(), start, units = "secs"))

  success <- !is.null(fit) && is.null(error)

  y <- .benchmark_response(formula, data)

  pred <- if (success) .benchmark_predict_response(fit, data) else rep(NA_real_, length(y))

  scores <- .benchmark_score_predictions(as.numeric(y), pred)

  distribution_scores <- if (success && isTRUE(distributional_metrics)) {
    .benchmark_distribution_summaries(
      y = as.numeric(y),
      fitted = pred,
      family = family,
      data = data,
      formula = formula,
      truth_family = truth_family,
      p = quantile_prob,
      interval_level = interval_level
    )
  } else {
    .benchmark_distribution_metric_empty()
  }

  ll <- if (success && !identical(comparator, "gee")) {
    as.numeric(tryCatch(stats::logLik(fit), error = function(e) NA_real_))
  } else {
    NA_real_
  }

  aic <- if (success && !identical(comparator, "gee")) {
    as.numeric(tryCatch(stats::AIC(fit), error = function(e) NA_real_))
  } else {
    NA_real_
  }


  list(
    fit = fit,
    row = data.frame(
      method = comparator,
      comparator = comparator,
      comparator_class = status_row$comparator_class,
      estimator = status_row$estimator,
      package = status_row$package,
      available = TRUE,
      success = success,
      elapsed_sec = elapsed,
      nobs = nrow(data),
      logLik = ll,
      AIC = aic,
      mae = unname(scores[["mae"]]),
      rmse = unname(scores[["rmse"]]),
      benchmark_mae = unname(scores[["mae"]]),
      benchmark_rmse = unname(scores[["rmse"]]),
      as.list(distribution_scores),
      warning = if (length(warnings)) paste(unique(warnings), collapse = " | ") else NA_character_,
      error = if (is.null(error)) NA_character_ else error,
      stringsAsFactors = FALSE
    )
  )
}

# ---- benchmark-comparator-supplied.R ----

.benchmark_scalar_fit_stat <- function(expr, preferred = "joint") {
  value <- tryCatch(expr, error = function(e) NA_real_)

  if (is.data.frame(value)) {
    if (preferred %in% names(value)) {
      value <- value[[preferred]]
    } else {
      value <- value[[1L]]
    }
  } else if (is.matrix(value)) {
    if (preferred %in% colnames(value)) {
      value <- value[, preferred]
    } else if (preferred %in% rownames(value)) {
      value <- value[preferred, ]
    } else {
      value <- value[[1L]]
    }
  }

  if (length(value) == 0L) {
    return(NA_real_)
  }

  if (!is.null(names(value)) && preferred %in% names(value)) {
    value <- value[[preferred]]
  } else {
    value <- value[[1L]]
  }

  value <- suppressWarnings(as.numeric(value))

  if (length(value) == 0L || !is.finite(value[[1L]])) {
    return(NA_real_)
  }

  value[[1L]]
}

.benchmark_supplied_fit_aic <- function(fit) {
  aic <- .benchmark_scalar_fit_stat(stats::AIC(fit))

  if (is.finite(aic)) {
    return(aic)
  }

  if (inherits(fit, "gamlss.longitudinal")) {
    summary_obj <- tryCatch(summary(fit, include_vcov = FALSE), error = function(e) NULL)

    value <- summary_obj$fit$AIC %||% NA_real_

    value <- suppressWarnings(as.numeric(value))

    if (length(value) > 0L && is.finite(value[[1L]])) {
      return(value[[1L]])
    }
  }

  NA_real_
}

.benchmark_supplied_fit_one <- function(
    fit,
    data,
    formula,
    family,
    fit_name,
    distributional_metrics = TRUE,
    truth_family = NULL,
    quantile_prob = 0.9,
    interval_level = 0.95) {
  if (is.null(fit)) {
    return(NULL)
  }

  fit_name <- as.character(fit_name)[1L]

  class_name <- class(fit)[1L]

  comparator_class <- if (inherits(fit, "gamlss.longitudinal")) {
    "gamlss_longitudinal"
  } else {
    class_name
  }

  estimator <- if (inherits(fit, "gamlss.longitudinal")) {
    "gamlss.longitudinal::gamlss_longitudinal"
  } else {
    class_name
  }

  package <- if (inherits(fit, "gamlss.longitudinal")) {
    "gamlss.longitudinal"
  } else {
    NA_character_
  }


  y <- .benchmark_response(formula, data)

  pred <- .benchmark_predict_response(fit, data)

  scores <- .benchmark_score_predictions(as.numeric(y), pred)

  success <- any(is.finite(pred))

  distribution_scores <- if (success && isTRUE(distributional_metrics) && inherits(fit, "gamlss.longitudinal")) {
    .benchmark_gamlss_distribution_summaries(
      fit = fit,
      data = data,
      formula = formula,
      truth_family = truth_family,
      p = quantile_prob,
      interval_level = interval_level
    )
  } else if (success && isTRUE(distributional_metrics)) {
    .benchmark_distribution_summaries(
      y = as.numeric(y),
      fitted = pred,
      family = family,
      data = data,
      formula = formula,
      truth_family = truth_family,
      p = quantile_prob,
      interval_level = interval_level
    )
  } else {
    .benchmark_distribution_metric_empty()
  }

  ll <- if (success) .benchmark_scalar_fit_stat(stats::logLik(fit)) else NA_real_

  aic <- if (success) .benchmark_supplied_fit_aic(fit) else NA_real_


  list(
    fit = fit,
    row = data.frame(
      method = fit_name,
      comparator = fit_name,
      comparator_class = comparator_class,
      estimator = estimator,
      package = package,
      available = TRUE,
      success = success,
      elapsed_sec = NA_real_,
      nobs = nrow(data),
      logLik = ll,
      AIC = aic,
      mae = unname(scores[["mae"]]),
      rmse = unname(scores[["rmse"]]),
      benchmark_mae = unname(scores[["mae"]]),
      benchmark_rmse = unname(scores[["rmse"]]),
      as.list(distribution_scores),
      warning = NA_character_,
      error = if (success) NA_character_ else "Prediction failed or returned no finite response-scale values.",
      stringsAsFactors = FALSE
    )
  )
}

# ---- benchmark-comparator-coefficients.R ----

.benchmark_coefficient_comparison <- function(fits, results, primary_method = NULL, parameter = "mu", level = 0.95) {
  if (length(fits) == 0L) {
    return(list(
      parameter = parameter,
      level = level,
      reference_method = NA_character_,
      long = data.frame(),
      estimates = data.frame(),
      uncertainty = data.frame()
    ))
  }

  results <- as.data.frame(results, stringsAsFactors = FALSE)

  success <- if ("success" %in% names(results)) {
    stats::setNames(results$success == TRUE & !is.na(results$success), results$method)
  } else {
    stats::setNames(rep(TRUE, nrow(results)), results$method)
  }

  rows <- lapply(names(fits), function(method) {
    if (!isTRUE(success[[method]] %||% TRUE)) {
      return(data.frame())
    }

    .benchmark_coef_table_one(fits[[method]], method = method, parameter = parameter, level = level)
  })

  long <- do.call(rbind, rows)

  long <- as.data.frame(long, stringsAsFactors = FALSE)

  if (nrow(long) == 0L) {
    return(list(
      parameter = parameter,
      level = level,
      reference_method = NA_character_,
      long = long,
      estimates = data.frame(),
      uncertainty = data.frame()
    ))
  }

  reference_method <- primary_method

  if (is.null(reference_method) || length(reference_method) == 0L || !reference_method %in% long$method) {
    reference_method <- long$method[[1L]]
  }

  estimates <- .benchmark_wide_values(long, "estimate")


  reference <- long[long$method == reference_method, , drop = FALSE]

  uncertainty <- long

  uncertainty$reference_method <- reference_method

  uncertainty$estimate_diff_vs_reference <- NA_real_

  uncertainty$se_ratio_vs_reference <- NA_real_

  uncertainty$ci_width <- uncertainty$conf.high - uncertainty$conf.low

  uncertainty$ci_width_ratio_vs_reference <- NA_real_

  uncertainty$ci_overlap_with_reference <- NA

  for (i in seq_len(nrow(uncertainty))) {
    ref <- reference[reference$term == uncertainty$term[[i]], , drop = FALSE]

    if (nrow(ref) != 1L) next

    uncertainty$estimate_diff_vs_reference[[i]] <- uncertainty$estimate[[i]] - ref$estimate[[1L]]

    uncertainty$se_ratio_vs_reference[[i]] <- uncertainty$std_error[[i]] / ref$std_error[[1L]]

    ref_width <- ref$conf.high[[1L]] - ref$conf.low[[1L]]

    uncertainty$ci_width_ratio_vs_reference[[i]] <- uncertainty$ci_width[[i]] / ref_width

    uncertainty$ci_overlap_with_reference[[i]] <- (

      uncertainty$conf.low[[i]] <= ref$conf.high[[1L]] &&

        uncertainty$conf.high[[i]] >= ref$conf.low[[1L]]

    )
  }


  list(
    parameter = parameter,
    level = level,
    reference_method = reference_method,
    long = long,
    estimates = estimates,
    uncertainty = uncertainty
  )
}

# ---- benchmark-comparator-coefficient-tables.R ----

.benchmark_coef_matrix_table <- function(mat, method, parameter = "mu", level = 0.95) {
  mat <- as.matrix(mat)

  if (nrow(mat) == 0L || ncol(mat) == 0L) {
    return(data.frame())
  }

  estimate_col <- grep("^Estimate$", colnames(mat), ignore.case = TRUE)

  if (length(estimate_col) == 0L) estimate_col <- 1L

  se_col <- grep("Std\\.? Error|Std\\.err|Std Error", colnames(mat), ignore.case = TRUE)

  if (length(se_col) == 0L && ncol(mat) >= 2L) se_col <- 2L

  estimate <- suppressWarnings(as.numeric(mat[, estimate_col[[1L]]]))

  std_error <- if (length(se_col) > 0L) {
    suppressWarnings(as.numeric(mat[, se_col[[1L]]]))
  } else {
    rep(NA_real_, length(estimate))
  }

  z <- stats::qnorm((1 + level) / 2)

  out <- data.frame(
    method = method,
    parameter = parameter,
    term = rownames(mat) %||% paste0("term", seq_along(estimate)),
    estimate = estimate,
    std_error = std_error,
    conf.low = estimate - z * std_error,
    conf.high = estimate + z * std_error,
    stringsAsFactors = FALSE
  )

  rownames(out) <- NULL

  out
}

.benchmark_normalize_coef_term <- function(term, time_var = NULL) {
  term <- as.character(term)

  term[term %in% c("(Intercept)", "Intercept", "mu.(Intercept)", "mu.Intercept")] <- "intercept"

  term <- sub("^mu\\.", "", term)

  if (!is.null(time_var)) {
    time_var <- as.character(time_var)[1L]

    if (!is.na(time_var) && nzchar(time_var)) {
      term <- gsub("(^|:)time_covariate($|:)", paste0("\\1", time_var, "\\2"), term, perl = TRUE)
    }
  }

  term[term %in% c("(Intercept)", "Intercept")] <- "intercept"

  term
}

.benchmark_coef_table_one <- function(fit, method, parameter = "mu", level = 0.95) {
  if (is.null(fit)) {
    return(data.frame())
  }

  out <- tryCatch(
    {
      if (inherits(fit, "gamlss.longitudinal")) {
        s <- summary(fit, include_vcov = TRUE, ci_level = level)

        tbl <- as.data.frame(s$coefficients, stringsAsFactors = FALSE)

        tbl <- tbl[tbl$parameter == parameter, , drop = FALSE]

        if (nrow(tbl) == 0L) {
          return(data.frame())
        }

        z <- stats::qnorm((1 + level) / 2)

        data.frame(
          method = method,
          parameter = parameter,
          term = sub(paste0("^", parameter, "\\."), "", tbl$term),
          estimate = tbl$estimate,
          std_error = tbl$std_error,
          conf.low = tbl$estimate - z * tbl$std_error,
          conf.high = tbl$estimate + z * tbl$std_error,
          stringsAsFactors = FALSE
        )
      } else if (inherits(fit, "gam")) {
        s <- summary(fit)

        .benchmark_coef_matrix_table(s$p.table, method = method, parameter = parameter, level = level)
      } else if (inherits(fit, "glmmTMB")) {
        s <- summary(fit)

        .benchmark_coef_matrix_table(s$coefficients$cond, method = method, parameter = parameter, level = level)
      } else {
        mat <- coef(summary(fit))

        .benchmark_coef_matrix_table(mat, method = method, parameter = parameter, level = level)
      }
    },
    error = function(e) data.frame()
  )

  out <- as.data.frame(out, stringsAsFactors = FALSE)

  if (nrow(out) == 0L) {
    return(out)
  }

  time_var <- if (inherits(fit, "gamlss.longitudinal")) fit$time_var else NULL

  out$term <- .benchmark_normalize_coef_term(out$term, time_var = time_var)

  out <- out[is.finite(out$estimate), , drop = FALSE]

  rownames(out) <- NULL

  out
}

.benchmark_wide_values <- function(x, value_col) {
  x <- as.data.frame(x, stringsAsFactors = FALSE)

  if (nrow(x) == 0L || !value_col %in% names(x)) {
    return(data.frame())
  }

  terms <- unique(x$term)

  methods <- unique(x$method)

  out <- data.frame(term = terms, stringsAsFactors = FALSE)

  for (method in methods) {
    values <- x[x$method == method, c("term", value_col), drop = FALSE]

    out[[method]] <- values[[value_col]][match(terms, values$term)]
  }

  out
}

# ---- benchmark-comparator-distribution-metrics.R ----

.benchmark_distribution_metric_empty <- function() {
  c(
    benchmark_mean_bias = NA_real_,
    benchmark_mean_mae = NA_real_,
    benchmark_mean_rmse = NA_real_,
    benchmark_q90_mae = NA_real_,
    benchmark_neg_log_score = NA_real_,
    benchmark_upper_tail_error_90 = NA_real_,
    benchmark_theta_time_abs_error = NA_real_,
    benchmark_interval_coverage_95 = NA_real_,
    benchmark_interval_width_95 = NA_real_,
    benchmark_pit_ks_p_value = NA_real_,
    benchmark_pit_mean_abs_error = NA_real_,
    benchmark_tail_error_lower_05 = NA_real_,
    benchmark_tail_error_upper_05 = NA_real_
  )
}

.benchmark_truth_family <- function(family, truth_family = NULL) {
  if (!is.null(truth_family)) {
    return(as.character(truth_family)[1L])
  }

  if (identical(family$family, "gaussian")) {
    return("NO")
  }

  if (identical(family$family, "Gamma")) {
    return("GA")
  }

  if (identical(family$family, "poisson")) {
    return("PO")
  }

  if (identical(family$family, "binomial")) {
    return("BI")
  }

  NA_character_
}

.benchmark_response_data <- function(data, formula) {
  response_var <- all.vars(stats::as.formula(formula))[1L]

  dat <- as.data.frame(data, stringsAsFactors = FALSE)

  dat$response <- dat[[response_var]]

  dat
}

.benchmark_distribution_summaries <- function(y, fitted, family, data, formula, truth_family = NULL, p = 0.9, interval_level = 0.95) {
  out <- .benchmark_distribution_metric_empty()

  dat <- .benchmark_response_data(data, formula)


  if ("true_mu" %in% names(dat)) {
    ok_truth <- is.finite(dat$true_mu) & is.finite(fitted)

    if (any(ok_truth)) {
      err_mu <- fitted[ok_truth] - dat$true_mu[ok_truth]

      out["benchmark_mean_bias"] <- mean(err_mu)

      out["benchmark_mean_mae"] <- mean(abs(err_mu))

      out["benchmark_mean_rmse"] <- sqrt(mean(err_mu^2))
    }
  }


  pred_dist <- .benchmark_predictive_distribution(y, fitted, family, p = p, interval_level = interval_level)

  truth_dist <- .coverage_true_margin_distribution(

    dat,
    family = .benchmark_truth_family(family, truth_family = truth_family),
    p = p
  )


  ok_q <- is.finite(pred_dist$q_p) & is.finite(truth_dist$q)

  if (any(ok_q)) {
    out["benchmark_q90_mae"] <- mean(abs(pred_dist$q_p[ok_q] - truth_dist$q[ok_q]), na.rm = TRUE)
  }

  ok_density <- is.finite(pred_dist$density) & pred_dist$density > 0

  if (any(ok_density)) {
    out["benchmark_neg_log_score"] <- mean(-log(pmax(pred_dist$density[ok_density], .Machine$double.xmin)), na.rm = TRUE)
  }


  comparator_at_truth <- .coverage_comparator_distribution(y, fitted, family, truth_dist$q, p = p)

  ok_tail <- is.finite(comparator_at_truth$cdf_at_q) & is.finite(truth_dist$cdf)

  if (any(ok_tail)) {
    out["benchmark_upper_tail_error_90"] <- mean(

      (1 - comparator_at_truth$cdf_at_q[ok_tail]) - (1 - truth_dist$cdf[ok_tail]),
      na.rm = TRUE
    )
  }


  ok_interval <- is.finite(y) & is.finite(pred_dist$lower) & is.finite(pred_dist$upper)

  if (any(ok_interval)) {
    out["benchmark_interval_coverage_95"] <- mean(y[ok_interval] >= pred_dist$lower[ok_interval] & y[ok_interval] <= pred_dist$upper[ok_interval])

    out["benchmark_interval_width_95"] <- mean(pred_dist$upper[ok_interval] - pred_dist$lower[ok_interval], na.rm = TRUE)
  }


  ok_pit <- is.finite(pred_dist$pit)

  if (sum(ok_pit) >= 3L) {
    pit <- pred_dist$pit[ok_pit]

    out["benchmark_pit_ks_p_value"] <- tryCatch(

      suppressWarnings(stats::ks.test(pit, "punif")$p.value),
      error = function(e) NA_real_
    )

    out["benchmark_pit_mean_abs_error"] <- abs(mean(pit, na.rm = TRUE) - 0.5)

    out["benchmark_tail_error_lower_05"] <- mean(pit <= 0.05, na.rm = TRUE) - 0.05

    out["benchmark_tail_error_upper_05"] <- mean(pit >= 0.95, na.rm = TRUE) - 0.05
  }


  out
}

# ---- benchmark-comparator-gamlss-metrics.R ----

.benchmark_predict_gamlss_vector <- function(fit, data, type, value_name = NULL, ...) {
  pred <- tryCatch(stats::predict(fit, newdata = data, type = type, ...), error = function(e) NULL)

  if (is.null(pred)) {
    return(rep(NA_real_, nrow(data)))
  }

  if (is.data.frame(pred)) {
    if (!is.null(value_name) && value_name %in% names(pred)) {
      pred <- pred[[value_name]]
    } else {
      pred <- pred[[ncol(pred)]]
    }
  }

  pred <- as.numeric(pred)

  if (length(pred) != nrow(data)) {
    return(rep(NA_real_, nrow(data)))
  }

  pred
}

.benchmark_gamlss_distribution_summaries <- function(fit, data, formula, truth_family = NULL, p = 0.9, interval_level = 0.95) {
  out <- .benchmark_distribution_metric_empty()

  y <- as.numeric(.benchmark_response(formula, data))

  fitted <- .benchmark_predict_response(fit, data)

  dat <- .benchmark_response_data(data, formula)


  if ("true_mu" %in% names(dat)) {
    ok_truth <- is.finite(dat$true_mu) & is.finite(fitted)

    if (any(ok_truth)) {
      err_mu <- fitted[ok_truth] - dat$true_mu[ok_truth]

      out["benchmark_mean_bias"] <- mean(err_mu)

      out["benchmark_mean_mae"] <- mean(abs(err_mu))

      out["benchmark_mean_rmse"] <- sqrt(mean(err_mu^2))
    }
  }


  alpha <- (1 - interval_level) / 2

  q_pred <- tryCatch(stats::predict(fit, newdata = data, type = "quantile", probs = c(alpha, 1 - alpha, p)), error = function(e) NULL)

  q_cols <- if (is.data.frame(q_pred) && ncol(q_pred) >= 6L) q_pred[(ncol(q_pred) - 2L):ncol(q_pred)] else NULL

  lower <- if (!is.null(q_cols)) as.numeric(q_cols[[1L]]) else rep(NA_real_, nrow(data))

  upper <- if (!is.null(q_cols)) as.numeric(q_cols[[2L]]) else rep(NA_real_, nrow(data))

  q_p <- if (!is.null(q_cols)) as.numeric(q_cols[[3L]]) else rep(NA_real_, nrow(data))


  pit <- .benchmark_predict_gamlss_vector(fit, data, type = "cdf", value_name = "cdf", q = y)

  density <- .benchmark_predict_gamlss_vector(fit, data, type = "density", value_name = "density", y = y)

  truth_dist <- .coverage_true_margin_distribution(

    dat,
    family = .benchmark_truth_family(stats::gaussian(), truth_family = truth_family %||% fit$margin_dist$family[1]),
    p = p
  )

  cdf_at_truth_q <- .benchmark_predict_gamlss_vector(fit, data, type = "cdf", value_name = "cdf", q = truth_dist$q)


  ok_q <- is.finite(q_p) & is.finite(truth_dist$q)

  if (any(ok_q)) {
    out["benchmark_q90_mae"] <- mean(abs(q_p[ok_q] - truth_dist$q[ok_q]), na.rm = TRUE)
  }

  ok_density <- is.finite(density) & density > 0

  if (any(ok_density)) {
    out["benchmark_neg_log_score"] <- mean(-log(pmax(density[ok_density], .Machine$double.xmin)), na.rm = TRUE)
  }


  ok_tail <- is.finite(cdf_at_truth_q) & is.finite(truth_dist$cdf)

  if (any(ok_tail)) {
    out["benchmark_upper_tail_error_90"] <- mean((1 - cdf_at_truth_q[ok_tail]) - (1 - truth_dist$cdf[ok_tail]), na.rm = TRUE)
  }


  ok_interval <- is.finite(y) & is.finite(lower) & is.finite(upper)

  if (any(ok_interval)) {
    out["benchmark_interval_coverage_95"] <- mean(y[ok_interval] >= lower[ok_interval] & y[ok_interval] <= upper[ok_interval])

    out["benchmark_interval_width_95"] <- mean(upper[ok_interval] - lower[ok_interval], na.rm = TRUE)
  }


  ok_pit <- is.finite(pit)

  if (sum(ok_pit) >= 3L) {
    pit <- pmin(pmax(pit[ok_pit], 0), 1)

    out["benchmark_pit_ks_p_value"] <- tryCatch(

      suppressWarnings(stats::ks.test(pit, "punif")$p.value),
      error = function(e) NA_real_
    )

    out["benchmark_pit_mean_abs_error"] <- abs(mean(pit, na.rm = TRUE) - 0.5)

    out["benchmark_tail_error_lower_05"] <- mean(pit <= 0.05, na.rm = TRUE) - 0.05

    out["benchmark_tail_error_upper_05"] <- mean(pit >= 0.95, na.rm = TRUE) - 0.05
  }


  out
}

# ---- benchmark-comparator-predictive-distribution.R ----

.benchmark_predictive_distribution <- function(y, fitted, family, p = 0.9, interval_level = 0.95) {
  n <- length(fitted)

  empty <- list(
    q_p = rep(NA_real_, n),
    lower = rep(NA_real_, n),
    upper = rep(NA_real_, n),
    pit = rep(NA_real_, n),
    density = rep(NA_real_, n)
  )

  ok <- is.finite(y) & is.finite(fitted)

  if (!any(ok)) {
    return(empty)
  }

  alpha <- (1 - interval_level) / 2


  if (identical(family$family, "gaussian")) {
    sigma_hat <- .coverage_comparator_dispersion(y, fitted, family)

    if (!is.finite(sigma_hat) || sigma_hat <= 0) {
      return(empty)
    }

    q_p <- fitted + stats::qnorm(p) * sigma_hat

    lower <- fitted + stats::qnorm(alpha) * sigma_hat

    upper <- fitted + stats::qnorm(1 - alpha) * sigma_hat

    pit <- stats::pnorm(y, mean = fitted, sd = sigma_hat)

    density <- stats::dnorm(y, mean = fitted, sd = sigma_hat)
  } else if (identical(family$family, "poisson")) {
    lambda <- pmax(fitted, .Machine$double.eps)

    q_p <- stats::qpois(p, lambda = lambda)

    lower <- stats::qpois(alpha, lambda = lambda)

    upper <- stats::qpois(1 - alpha, lambda = lambda)

    pit <- stats::ppois(y, lambda = lambda)

    density <- stats::dpois(y, lambda = lambda)
  } else if (identical(family$family, "binomial")) {
    prob <- pmin(pmax(fitted, .Machine$double.eps), 1 - .Machine$double.eps)

    q_p <- stats::qbinom(p, size = 1, prob = prob)

    lower <- stats::qbinom(alpha, size = 1, prob = prob)

    upper <- stats::qbinom(1 - alpha, size = 1, prob = prob)

    pit <- stats::pbinom(y, size = 1, prob = prob)

    density <- stats::dbinom(y, size = 1, prob = prob)
  } else if (identical(family$family, "Gamma")) {
    dispersion <- .coverage_comparator_dispersion(y, fitted, family)

    if (!is.finite(dispersion) || dispersion <= 0) {
      return(empty)
    }

    mu <- pmax(fitted, .Machine$double.eps)

    shape <- 1 / dispersion

    scale <- mu * dispersion

    q_p <- stats::qgamma(p, shape = shape, scale = scale)

    lower <- stats::qgamma(alpha, shape = shape, scale = scale)

    upper <- stats::qgamma(1 - alpha, shape = shape, scale = scale)

    pit <- stats::pgamma(y, shape = shape, scale = scale)

    density <- stats::dgamma(y, shape = shape, scale = scale)
  } else {
    return(empty)
  }


  list(
    q_p = as.numeric(q_p),
    lower = as.numeric(lower),
    upper = as.numeric(upper),
    pit = pmin(pmax(as.numeric(pit), 0), 1),
    density = as.numeric(density)
  )
}

# ---- benchmark-comparator-status.R ----

#' List available standard longitudinal benchmark comparators

#'

#' @return A data frame describing optional comparator backends.

#' @export

benchmark_comparator_status <- function() {
  packages <- c("stats", "geepack", "lme4", "mgcv", "mgcv", "glmmTMB")

  data.frame(
    comparator = c("glm", "gee", "glmm", "gam", "gamm", "glmmTMB"),
    comparator_class = c("glm", "gee", "glmm", "gam", "gamm", "glmm"),
    estimator = c("stats::glm", "geepack::geeglm", "lme4::lmer/glmer", "mgcv::gam", "mgcv::gam + s(subject, bs = 're')", "glmmTMB::glmmTMB"),
    package = packages,
    available = vapply(packages, requireNamespace, logical(1), quietly = TRUE),
    role = c(
      "independence mean baseline",
      "marginal mean baseline with working correlation",
      "random-intercept conditional mean baseline",
      "independence smooth mean baseline",
      "smooth mean baseline with subject random-effect smooth",
      "optional flexible GLMM baseline"
    ),
    stringsAsFactors = FALSE
  )
}

# ---- benchmark-comparator-utils.R ----

.benchmark_family <- function(family) {
  if (inherits(family, "family")) {
    return(family)
  }

  family <- match.arg(as.character(family)[1], c("gaussian", "poisson", "binomial", "gamma"))

  switch(family,
    gaussian = stats::gaussian(),
    poisson = stats::poisson(),
    binomial = stats::binomial(),
    gamma = stats::Gamma(link = "log")
  )
}

.benchmark_backtick <- function(x) {
  if (grepl("^[.A-Za-z][.A-Za-z0-9_]*$", x)) {
    return(x)
  }

  paste0("`", gsub("`", "\\\\`", x), "`")
}

.benchmark_formula_with_random_intercept <- function(formula, subject_var) {
  stats::as.formula(paste(deparse(formula), "+ (1 |", .benchmark_backtick(subject_var), ")"))
}

.benchmark_formula_with_subject_re_smooth <- function(formula, subject_var) {
  stats::as.formula(paste(deparse(formula), "+ s(", .benchmark_backtick(subject_var), ", bs = 're')"))
}

.benchmark_formula_has_smooth <- function(formula) {
  term_obj <- tryCatch(

    stats::terms(stats::as.formula(formula), specials = "s"),
    error = function(e) NULL
  )

  specials <- attr(term_obj, "specials")

  length(specials$s) > 0L
}

.benchmark_response <- function(formula, data) {
  response_var <- all.vars(stats::as.formula(formula))[1L]

  if (!response_var %in% names(data)) {
    stop("The response variable from 'formula' was not found in 'data'.", call. = FALSE)
  }

  data[[response_var]]
}

.benchmark_predict_response <- function(fit, data, glmm_re_form = NA) {
  pred_type <- if (inherits(fit, "gamlss.longitudinal")) "mean" else "response"

  predict_args <- list(
    object = fit,
    newdata = data,
    type = pred_type,
    allow.new.levels = TRUE
  )

  if (inherits(fit, c("merMod", "glmmTMB"))) {
    predict_args$re.form <- glmm_re_form
  }

  pred <- tryCatch(

    do.call(stats::predict, predict_args),
    error = function(e) {
      tryCatch(stats::predict(fit, newdata = data, type = pred_type), error = function(e2) rep(NA_real_, nrow(data)))
    }
  )

  as.numeric(pred)
}

.benchmark_score_predictions <- function(y, fitted) {
  ok <- is.finite(y) & is.finite(fitted)

  if (!any(ok)) {
    return(c(mae = NA_real_, rmse = NA_real_))
  }

  err <- fitted[ok] - y[ok]

  c(
    mae = mean(abs(err)),
    rmse = sqrt(mean(err^2))
  )
}

.benchmark_empty_result_row <- function(
    status_row,
    comparator,
    available,
    success,
    elapsed_sec,
    nobs,
    warning = NA_character_,
    error = NA_character_) {
  data.frame(
    method = comparator,
    comparator = comparator,
    comparator_class = status_row$comparator_class,
    estimator = status_row$estimator,
    package = status_row$package,
    available = available,
    success = success,
    elapsed_sec = elapsed_sec,
    nobs = nobs,
    logLik = NA_real_,
    AIC = NA_real_,
    mae = NA_real_,
    rmse = NA_real_,
    benchmark_mae = NA_real_,
    benchmark_rmse = NA_real_,
    as.list(.benchmark_distribution_metric_empty()),
    warning = warning,
    error = error,
    stringsAsFactors = FALSE
  )
}

# ---- benchmark-comparator-print.R ----

#' @export
print.gamlss_longitudinal_benchmark <- function(x, digits = max(3, getOption("digits") - 3), ...) {
  cat("\nStandard Longitudinal Comparator Benchmark\n")

  cat("-----------------------------------------\n")

  cat("Formula:", deparse(x$formula), "\n")

  cat("Subject:", x$subject_var, "\n")

  cat("Family:", x$family, "\n")

  cat(.benchmark_included_methods_line(x$results), "\n")

  cat("Estimand note: mean-fit metrics compare response-mean predictions; distributional metrics use fitted quantiles, densities, PIT values, intervals, and tail probabilities where available.\n\n")


  interpretation_table <- .benchmark_interpretation_table(x$interpretation)

  if (nrow(interpretation_table) > 0L) {
    cat("Benchmark Summary\n")

    cat("-----------------\n")

    .benchmark_print_by_area(interpretation_table, digits = digits)
  } else {
    lines <- .benchmark_interpretation_lines(x$interpretation)

    if (length(lines) > 0L) {
      cat("Benchmark Summary\n")

      cat("-----------------\n")

      cat(paste(lines, collapse = "\n"), "\n", sep = "")
    }
  }


  cat("\nBenchmark Details\n")

  cat("-----------------\n")


  successful_results <- x$results

  if ("success" %in% names(successful_results)) {
    successful_results <- successful_results[!is.na(successful_results$success) & successful_results$success == TRUE, , drop = FALSE]
  }

  metric_table <- .benchmark_metric_matrix(successful_results)

  if (nrow(metric_table) > 0L) {
    cat("Benchmark Metrics\n")

    cat("-----------------\n")

    print(.benchmark_print_method_columns(.benchmark_print_metric_labels(metric_table)), digits = digits, row.names = FALSE)
  }


  .benchmark_print_coefficient_comparison(x$coefficients, digits = digits)

  invisible(x)
}
