#' Select margin and copula families by joint longitudinal fit
#'
#' `select_joint_distribution()` screens marginal distribution and copula
#' combinations by fitting intercept-only [gamlss_longitudinal()] models and
#' ranking their joint likelihood, AIC, or BIC. Unlike [select_margin()] and
#' [select_copula()], this selector evaluates each candidate as a full joint
#' longitudinal model, so it can take substantially longer on large candidate
#' sets.
#'
#' @param data Long-format data frame.
#' @param response_var Response column name in `data`.
#' @param time_var,subject_var Time and subject identifier column names.
#' @param type Optional `gamlss::fitDist()` type passed to [select_margin()].
#' @param margin_families Optional marginal family names. When supplied, these
#'   are used directly as the marginal candidate set; otherwise candidates come
#'   from [select_margin()].
#' @param copula_families Candidate copula family codes. Supported values are
#'   `"N"`, `"C"`, `"F"`, `"G"`, `"J"`, and `"t"`.
#' @param criterion Ranking criterion, one of `"AIC"`, `"BIC"`, or `"logLik"`.
#' @param min_pairs Minimum number of complete adjacent response pairs required
#'   before fitting candidates.
#' @param time_intercepts Logical; if `TRUE`, pass through to
#'   [select_margin()] and use time-specific intercepts for each marginal
#'   distribution parameter in the joint screening fits.
#' @param copula_time_intercepts Logical; if `TRUE`, use time-specific
#'   intercepts for the copula dependence parameter in the joint screening fits.
#'   Time is treated as a factor, not as a linear trend.
#' @param try.gamlss Passed to [select_margin()].
#' @param trace Logical; passed to [select_margin()] and used to decide whether
#'   candidate [gamlss_longitudinal()] fit output is shown.
#' @param progress Logical; if `TRUE`, print candidate-level progress messages
#'   with elapsed time and an estimated remaining runtime.
#' @param fit_args Optional named list of arguments overriding the default
#'   [gamlss_longitudinal()] screening fit settings.
#' @param keep_fits Logical; if `TRUE`, attach successful fitted models in the
#'   `"fits"` attribute.
#'
#' @return A data frame with one row per margin-copula combination and class
#'   `joint_distribution_selection`. Failed fits are retained with an `error`
#'   message and missing fit metrics.
#' @export
select_joint_distribution <- function(
  data,
  response_var,
  time_var = "time",
  subject_var = "subject",
  type = NULL,
  margin_families = NULL,
  copula_families = c("N", "C", "F", "G", "J", "t"),
  criterion = c("AIC", "BIC", "logLik"),
  min_pairs = 10,
  time_intercepts = FALSE,
  copula_time_intercepts = FALSE,
  try.gamlss = FALSE,
  trace = FALSE,
  progress = TRUE,
  fit_args = list(),
  keep_fits = FALSE
) {
  criterion <- match.arg(criterion)
  if (!is.data.frame(data)) {
    data <- as.data.frame(data, stringsAsFactors = FALSE)
  } else {
    data <- as.data.frame(data, stringsAsFactors = FALSE)
  }
  .joint_selection_check_column(data, response_var, "response_var")
  .joint_selection_check_column(data, time_var, "time_var")
  .joint_selection_check_column(data, subject_var, "subject_var")
  if (!is.list(fit_args) ||
      (length(fit_args) > 0L && (is.null(names(fit_args)) || any(!nzchar(names(fit_args)))))) {
    stop("'fit_args' must be a named list.", call. = FALSE)
  }
  if (!is.logical(progress) || length(progress) != 1L || is.na(progress)) {
    stop("'progress' must be TRUE or FALSE.", call. = FALSE)
  }
  min_pairs <- as.integer(min_pairs)
  if (length(min_pairs) != 1L || !is.finite(min_pairs) || min_pairs < 1L) {
    stop("'min_pairs' must be a positive integer.", call. = FALSE)
  }
  n_pairs <- .joint_selection_count_pairs(
    data = data,
    response_var = response_var,
    time_var = time_var,
    subject_var = subject_var
  )
  if (n_pairs < min_pairs) {
    stop("At least ", min_pairs, " complete adjacent response pairs are required.", call. = FALSE)
  }

  copula_families <- unique(vapply(
    copula_families,
    .copula_family_code,
    character(1),
    USE.NAMES = FALSE
  ))
  if (length(copula_families) < 1L) {
    stop("'copula_families' must contain at least one supported family code.", call. = FALSE)
  }

  margin_selection <- .joint_selection_margin_candidates(
    data = data,
    response_var = response_var,
    time_var = time_var,
    type = type,
    margin_families = margin_families,
    time_intercepts = time_intercepts,
    try.gamlss = try.gamlss,
    trace = trace
  )
  margin_candidates <- as.data.frame(margin_selection)
  margin_candidates <- margin_candidates[
    isTRUE(length(margin_candidates$supported_by_longitudinal) > 0L) &
      margin_candidates$supported_by_longitudinal %in% TRUE,
    ,
    drop = FALSE
  ]
  if (nrow(margin_candidates) < 1L) {
    stop("No supported marginal families were retained for joint screening.", call. = FALSE)
  }

  combinations <- expand.grid(
    margin_family = margin_candidates$family,
    copula_family = copula_families,
    stringsAsFactors = FALSE
  )

  fit_store <- vector("list", nrow(combinations))
  rows <- vector("list", nrow(combinations))
  screen_start <- Sys.time()
  n_combinations <- nrow(combinations)
  if (isTRUE(progress)) {
    message(
      "Joint distribution screen: fitting ",
      n_combinations,
      " model(s) across ",
      nrow(margin_candidates),
      " margin candidate(s) and ",
      length(copula_families),
      " copula candidate(s)."
    )
  }
  for (ii in seq_len(nrow(combinations))) {
    margin_family <- combinations$margin_family[[ii]]
    copula_family <- combinations$copula_family[[ii]]
    candidate_label <- paste(margin_family, copula_family, sep = "+")
    if (isTRUE(progress)) {
      message("[", ii, "/", n_combinations, "] Fitting ", candidate_label, "...")
    }
    margin_dist <- .margin_family_object(margin_family)
    fit_one <- .joint_selection_fit_one(
      data = data,
      response_var = response_var,
      time_var = time_var,
      subject_var = subject_var,
      margin_family = margin_family,
      margin_dist = margin_dist,
      copula_family = copula_family,
      fit_args = fit_args,
      time_intercepts = time_intercepts,
      copula_time_intercepts = copula_time_intercepts,
      trace = trace
    )
    rows[[ii]] <- fit_one$row
    if (isTRUE(keep_fits)) {
      fit_store[[ii]] <- fit_one$fit
    }
    if (isTRUE(progress)) {
      elapsed_total <- as.numeric(difftime(Sys.time(), screen_start, units = "secs"))
      elapsed_fit <- as.numeric(fit_one$row$elapsed_sec[[1L]])
      remaining <- if (ii < n_combinations && elapsed_total > 0) {
        elapsed_total / ii * (n_combinations - ii)
      } else {
        0
      }
      status <- if (!is.na(fit_one$row$error[[1L]])) {
        paste0("failed: ", fit_one$row$error[[1L]])
      } else if (isTRUE(fit_one$row$converged[[1L]])) {
        "completed"
      } else {
        "completed without confirmed convergence"
      }
      message(
        "[", ii, "/", n_combinations, "] ",
        candidate_label,
        " ", status,
        " in ", .joint_selection_format_seconds(elapsed_fit),
        ". Elapsed ", .joint_selection_format_seconds(elapsed_total),
        "; ETA ", .joint_selection_format_seconds(remaining),
        "."
      )
    }
  }

  out <- do.call(rbind, rows)
  out$n_pairs <- n_pairs
  out$rank <- NA_integer_
  success <- is.na(out$error) & is.finite(out[[criterion]])
  ord_success <- order(out[success, criterion], decreasing = identical(criterion, "logLik"))
  success_idx <- which(success)[ord_success]
  failed_idx <- which(!success)
  out <- out[c(success_idx, failed_idx), , drop = FALSE]
  if (any(success)) {
    out$rank[seq_along(success_idx)] <- seq_along(success_idx)
  }
  out$delta_AIC <- if (any(is.finite(out$AIC))) out$AIC - min(out$AIC, na.rm = TRUE) else NA_real_
  rownames(out) <- NULL

  attr(out, "selected") <- if (length(success_idx) > 0L) {
    paste(out$margin_family[[1L]], out$copula_family[[1L]], sep = "+")
  } else {
    NA_character_
  }
  attr(out, "criterion") <- criterion
  attr(out, "margin_selection") <- margin_selection
  attr(out, "response_type") <- attr(margin_selection, "response_type")
  attr(out, "time_intercepts") <- isTRUE(time_intercepts)
  attr(out, "time_var") <- if (isTRUE(time_intercepts)) time_var else NULL
  attr(out, "copula_time_intercepts") <- isTRUE(copula_time_intercepts)
  attr(out, "copula_time_var") <- if (isTRUE(copula_time_intercepts)) time_var else NULL
  if (isTRUE(keep_fits)) {
    fit_store <- fit_store[c(success_idx, failed_idx)]
    attr(out, "fits") <- fit_store
  }
  class(out) <- c("joint_distribution_selection", "data.frame")
  out
}

.joint_selection_format_seconds <- function(seconds) {
  seconds <- as.numeric(seconds)[1L]
  if (!is.finite(seconds) || seconds < 0) {
    return("unknown")
  }
  if (seconds < 60) {
    return(sprintf("%.1fs", seconds))
  }
  minutes <- floor(seconds / 60)
  remaining <- seconds - minutes * 60
  if (minutes < 60) {
    return(sprintf("%dm %.0fs", minutes, remaining))
  }
  hours <- floor(minutes / 60)
  minutes <- minutes - hours * 60
  sprintf("%dh %dm", hours, minutes)
}

.joint_selection_margin_candidates <- function(
  data,
  response_var,
  time_var,
  type,
  margin_families,
  time_intercepts,
  try.gamlss,
  trace
) {
  if (!is.null(margin_families)) {
    if (!is.character(margin_families) || any(is.na(margin_families)) || any(!nzchar(margin_families))) {
      stop("'margin_families' must be a character vector of gamlss.dist family names.", call. = FALSE)
    }
    margin_families <- unique(margin_families)
    if (is.null(type)) {
      response <- as.numeric(data[[response_var]])
      response <- response[is.finite(response)]
      is_count <- all(response >= 0) && all(abs(response - round(response)) < .Machine$double.eps^0.5)
      type <- if (is_count) {
        "counts"
      } else if (all(response > 0)) {
        "realplus"
      } else {
        "realAll"
      }
    }
    out <- data.frame(
      family = margin_families,
      AIC = NA_real_,
      type = type,
      stringsAsFactors = FALSE
    )
    out$supported_by_longitudinal <- vapply(out$family, .joint_selection_margin_supported, logical(1))
    out$rank <- seq_len(nrow(out))
    out$delta_AIC <- NA_real_
    attr(out, "selected") <- if (nrow(out) > 0L) out$family[[1L]] else NA_character_
    attr(out, "response_type") <- type
    attr(out, "time_intercepts") <- isTRUE(time_intercepts)
    attr(out, "time_var") <- if (isTRUE(time_intercepts)) time_var else NULL
    class(out) <- c("margin_selection", "margin_screen", "data.frame")
    return(out)
  }

  margin_call <- quote(select_margin(
    data = data,
    response_var = response_var,
    time_var = if (isTRUE(time_intercepts)) time_var else NULL,
    type = type,
    families = NULL,
    time_intercepts = time_intercepts,
    try.gamlss = try.gamlss,
    trace = trace
  ))
  if (isTRUE(trace)) {
    eval(margin_call)
  } else {
    warnings <- character(0)
    withCallingHandlers(
      {
        utils::capture.output(ans <- eval(margin_call))
        ans
      },
      warning = function(w) {
        warnings <<- c(warnings, conditionMessage(w))
        invokeRestart("muffleWarning")
      }
    )
  }
}

.joint_selection_margin_supported <- function(family) {
  !is.null(.margin_family_object(family)) &&
    all(vapply(
      paste0(c("d", "p", "q"), family),
      exists,
      logical(1),
      envir = asNamespace("gamlss.dist"),
      inherits = FALSE
    ))
}

.joint_selection_check_column <- function(data, column, arg) {
  if (is.null(column) || !is.character(column) || length(column) != 1L || is.na(column) || !nzchar(column)) {
    stop("'", arg, "' must be a single column name.", call. = FALSE)
  }
  if (!column %in% names(data)) {
    stop(arg, "='", column, "' not found in 'data'.", call. = FALSE)
  }
  invisible(TRUE)
}

.joint_selection_count_pairs <- function(data, response_var, time_var, subject_var) {
  ord <- order(data[[subject_var]], data[[time_var]])
  data <- data[ord, , drop = FALSE]
  counts <- vapply(split(data, data[[subject_var]], drop = TRUE), function(subject_data) {
    subject_data <- subject_data[order(subject_data[[time_var]]), , drop = FALSE]
    if (nrow(subject_data) < 2L) {
      return(0L)
    }
    left <- seq_len(nrow(subject_data) - 1L)
    right <- left + 1L
    sum(is.finite(subject_data[[response_var]][left]) & is.finite(subject_data[[response_var]][right]))
  }, integer(1))
  sum(counts)
}

.joint_selection_fit_one <- function(
  data,
  response_var,
  time_var,
  subject_var,
  margin_family,
  margin_dist,
  copula_family,
  fit_args,
  time_intercepts,
  copula_time_intercepts,
  trace
) {
  start_time <- Sys.time()
  warnings <- character(0)
  fit <- NULL
  error <- NA_character_

  if (is.null(margin_dist)) {
    error <- paste0("Could not construct gamlss.dist family object for '", margin_family, "'.")
  } else {
    mu_formula <- if (isTRUE(time_intercepts)) {
      stats::as.formula(paste(.select_copula_formula_name(response_var), "~ factor(", .select_copula_formula_name(time_var), ")"))
    } else {
      stats::reformulate("1", response = response_var)
    }
    par_formula <- if (isTRUE(time_intercepts)) {
      stats::as.formula(paste("~ factor(", .select_copula_formula_name(time_var), ")"))
    } else {
      ~1
    }
    theta_formula <- if (isTRUE(copula_time_intercepts)) {
      stats::as.formula(paste("~ factor(", .select_copula_formula_name(time_var), ")"))
    } else {
      ~1
    }
    default_args <- list(
      dataset = data,
      margin_dist = margin_dist,
      copula_dist = copula_family,
      time_var = time_var,
      subject_var = subject_var,
      mu.formula = mu_formula,
      sigma.formula = par_formula,
      nu.formula = par_formula,
      tau.formula = par_formula,
      theta.formula = theta_formula,
      zeta.formula = ~1,
      include_dlcopdpar = TRUE,
      compute_vcov = FALSE,
      warm_start_joint = FALSE,
      verbose = 0
    )
    args <- utils::modifyList(default_args, fit_args)
    fit <- tryCatch(
      withCallingHandlers(
        {
          if (isTRUE(trace)) {
            do.call(gamlss_longitudinal, args)
          } else {
            utils::capture.output(ans <- do.call(gamlss_longitudinal, args))
            ans
          }
        },
        warning = function(w) {
          warnings <<- c(warnings, conditionMessage(w))
          invokeRestart("muffleWarning")
        }
      ),
      error = function(e) {
        error <<- conditionMessage(e)
        NULL
      }
    )
  }

  elapsed <- as.numeric(difftime(Sys.time(), start_time, units = "secs"))
  if (inherits(fit, "gamlss.longitudinal")) {
    summary_fit <- summary(fit, include_vcov = FALSE)
    fit_metrics <- summary_fit$fit
    invalid_reason <- .joint_selection_invalid_fit_reason(fit, fit_metrics)
    if (is.null(invalid_reason)) {
      row <- data.frame(
        margin_family = margin_family,
        copula_family = copula_family,
        logLik = as.numeric(fit_metrics$logLik),
        AIC = as.numeric(fit_metrics$AIC),
        BIC = as.numeric(fit_metrics$BIC),
        EDF = as.numeric(fit_metrics$model_selection["EDF", "joint"]),
        converged = isTRUE(fit$convergence$converged),
        hit_outer_limit = isTRUE(fit$convergence$hit_outer_limit),
        elapsed_sec = elapsed,
        warnings = paste(unique(warnings), collapse = "\n"),
        error = NA_character_,
        stringsAsFactors = FALSE
      )
    } else {
      row <- .joint_selection_failed_row(
        margin_family = margin_family,
        copula_family = copula_family,
        elapsed = elapsed,
        warnings = warnings,
        error = invalid_reason
      )
    }
  } else {
    row <- .joint_selection_failed_row(
      margin_family = margin_family,
      copula_family = copula_family,
      elapsed = elapsed,
      warnings = warnings,
      error = error
    )
  }
  list(row = row, fit = fit)
}

.joint_selection_failed_row <- function(margin_family, copula_family, elapsed, warnings, error) {
  data.frame(
    margin_family = margin_family,
    copula_family = copula_family,
    logLik = NA_real_,
    AIC = NA_real_,
    BIC = NA_real_,
    EDF = NA_real_,
    converged = FALSE,
    hit_outer_limit = NA,
    elapsed_sec = elapsed,
    warnings = paste(unique(warnings), collapse = "\n"),
    error = error,
    stringsAsFactors = FALSE
  )
}

.joint_selection_invalid_fit_reason <- function(fit, fit_metrics) {
  metric_values <- as.numeric(fit_metrics[c("logLik", "AIC", "BIC")])
  if (length(metric_values) != 3L || any(!is.finite(metric_values))) {
    return("Candidate fit did not provide finite likelihood criteria.")
  }

  model_selection <- fit_metrics$model_selection
  if (!is.null(model_selection) && all(c("marginal", "copula", "joint") %in% colnames(model_selection))) {
    final_loglik <- as.numeric(model_selection["LogLik", c("marginal", "copula", "joint")])
    history <- fit$log_lik_history
    if (length(final_loglik) == 3L &&
        all(is.finite(final_loglik)) &&
        all(abs(final_loglik) < .Machine$double.eps^0.5) &&
        !is.null(history) &&
        any(is.finite(history) & abs(history) > .Machine$double.eps^0.5)) {
      return("Candidate fit returned a zero final likelihood after non-zero likelihood history; treating fit as invalid.")
    }
  }

  NULL
}

#' @export
best_fit.joint_distribution_selection <- function(x, ...) {
  if (nrow(x) == 0L || !is.na(x$error[[1L]]) || !is.finite(x$AIC[[1L]])) {
    return(list(
      margin_family_name = NA_character_,
      margin_family = NULL,
      copula_family = NA_character_,
      criterion = attr(x, "criterion")
    ))
  }
  row <- as.list(as.data.frame(x)[1L, , drop = FALSE])
  margin_family_name <- row$margin_family
  copula_family <- row$copula_family
  row$margin_family <- NULL
  row$copula_family <- NULL
  c(
    list(
      margin_family_name = margin_family_name,
      margin_family = .margin_family_object(margin_family_name),
      copula_family = copula_family,
      criterion = attr(x, "criterion")
    ),
    row
  )
}

#' @export
best_fit_family.joint_distribution_selection <- function(x, ...) {
  best <- best_fit(x)
  list(
    margin_dist = best$margin_family,
    copula_dist = best$copula_family
  )
}

#' @export
`$.joint_distribution_selection` <- function(x, name) {
  if (identical(name, "best_fit")) {
    return(best_fit(x))
  }
  .subset2(as.data.frame(x), name, exact = FALSE)
}

#' @export
print.joint_distribution_selection <- function(x, ..., n = 10L) {
  cat("\nJoint Distribution Screen\n")
  cat("-------------------------\n")
  if (nrow(x) == 0L) {
    cat("No candidate combinations were retained.\n")
    return(invisible(x))
  }
  selected <- attr(x, "selected")
  cat("Selected:", if (length(selected) != 1L || is.na(selected)) "none" else selected, "\n")
  cat("Criterion:", attr(x, "criterion"), "\n\n")
  cols <- intersect(
    c(
      "rank", "margin_family", "copula_family", "logLik", "AIC", "BIC",
      "EDF", "converged", "elapsed_sec", "error"
    ),
    names(x)
  )
  display <- utils::head(as.data.frame(x)[cols], n = n)
  print(display, row.names = FALSE)
  if ("converged" %in% names(display) && any(!display$converged, na.rm = TRUE)) {
    cat("\nWarning: one or more printed joint candidate fits did not report convergence.\n")
  }
  invisible(x)
}
