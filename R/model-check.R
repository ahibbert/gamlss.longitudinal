.gl_residual_dependence_summary <- function(object, residual_lags = 1L) {

  out <- data.frame(lag = as.integer(residual_lags), normal_score_cor = NA_real_, n_pairs = 0L)

  if (length(residual_lags) == 0L) {

    return(data.frame(lag = integer(0), normal_score_cor = numeric(0), n_pairs = integer(0)))

  }


  copula_spec <- get_copula_dist(object$copula_dist)

  family_name <- .copula_family_code(copula_spec$copula_dist)

  family_num <- tryCatch(.copula_family_code(family_name), error = function(e) NA_character_)


  lag_summary <- tryCatch({

    fit_data <- .copula_v2_fit_data(object)

    rosenblatt_df <- .copula_v2_rosenblatt_series(fit_data, family_num)

    .copula_v2_rosenblatt_lag_summary(rosenblatt_df, lag_values = residual_lags)

  }, error = function(e) data.frame())


  if (nrow(lag_summary) == 0L) {

    return(out)

  }


  data.frame(

    lag = as.integer(lag_summary$lag),

    normal_score_cor = as.numeric(lag_summary$cor_z),

    n_pairs = as.integer(lag_summary$n_pairs),

    stringsAsFactors = FALSE

  )

}


.gl_pit_tail_summary <- function(pit, thresholds = c(0.05, 0.10)) {

  n <- sum(is.finite(pit))

  rows <- lapply(thresholds, function(threshold) {

    lower <- mean(pit <= threshold, na.rm = TRUE)

    upper <- mean(pit >= 1 - threshold, na.rm = TRUE)

    data.frame(

      threshold = threshold,

      lower = lower,

      upper = upper,

      expected = threshold,

      lower_ratio = lower / threshold,

      upper_ratio = upper / threshold,

      stringsAsFactors = FALSE

    )

  })

  out <- do.call(rbind, rows)

  if (n == 0L) {

    out[, c("lower", "upper", "lower_ratio", "upper_ratio")] <- NA_real_

  }

  tail_ratios <- c(out$lower_ratio, out$upper_ratio)

  attr(out, "tail_ratio_max") <- if (any(is.finite(tail_ratios))) {

    max(tail_ratios, na.rm = TRUE)

  } else {

    NA_real_

  }

  out

}


.gl_check_table <- function(summary_obj, scores, pit_stats, tail_stats, lag1_cor,

                            dependence_cor_cutoff = 0.25, vcov_method = NA_character_) {

  row <- function(area, quantity_checked, value, threshold_condition, default,

                  status, message, action) {

    severity <- switch(

      status,

      FAIL = "concern",

      REVIEW = "review",

      PASS = "ok",

      "ok"

    )

    data.frame(

      area = area,

      quantity_checked = quantity_checked,

      value = value,

      threshold_condition = threshold_condition,

      default = default,

      status = status,

      severity = severity,

      message = message,

      action = action,

      stringsAsFactors = FALSE

    )

  }


  converged <- isTRUE(summary_obj$convergence$converged)

  marginal_fail <- is.finite(pit_stats$ks_p_value) && pit_stats$ks_p_value < 0.05

  tail_fail <- is.finite(tail_stats$tail_ratio_max) && tail_stats$tail_ratio_max > 2

  copula_fail <- is.finite(lag1_cor) && abs(lag1_cor) > dependence_cor_cutoff

  variance_review <- identical(vcov_method, "numderiv")


  do.call(rbind, list(

    row(

      area = "Convergence",

      quantity_checked = "object$convergence$converged",

      value = if (converged) "TRUE" else "not TRUE",

      threshold_condition = "Not TRUE",

      default = "n/a",

      status = if (converged) "PASS" else "FAIL",

      message = if (converged) "Convergence was confirmed." else "Convergence was not confirmed.",

      action = if (converged) {

        "Continue with broader diagnostics."

      } else {

        "Refit with more iterations, different starts, or a simpler specification."

      }

    ),

    row(

      area = "Marginal fit",

      quantity_checked = "PIT Kolmogorov-Smirnov p-value vs Uniform(0, 1)",

      value = if (is.finite(pit_stats$ks_p_value)) formatC(pit_stats$ks_p_value, digits = 4, format = "fg") else NA_character_,

      threshold_condition = "ks_p_value < 0.05",

      default = "0.05",

      status = if (marginal_fail) "FAIL" else "PASS",

      message = if (marginal_fail) {

        "The marginal distribution is off by the PIT uniformity screen."

      } else {

        "The PIT uniformity screen did not flag marginal misfit."

      },

      action = if (marginal_fail) {

        "Inspect PIT, QQ, worm, and rootogram diagnostics; try a richer margin or covariate specification."

      } else {

        "Continue with visual marginal diagnostics."

      }

    ),

    row(

      area = "Tail fit",

      quantity_checked = "Maximum lower/upper PIT tail ratio over thresholds 0.05 and 0.10",

      value = if (is.finite(tail_stats$tail_ratio_max)) formatC(tail_stats$tail_ratio_max, digits = 4, format = "fg") else NA_character_,

      threshold_condition = "max(lower_ratio, upper_ratio) > 2",

      default = "2",

      status = if (tail_fail) "FAIL" else "PASS",

      message = if (tail_fail) {

        "Tail observations occur more often than the fitted margin expects."

      } else {

        "The basic PIT tail-ratio screen did not flag tail misfit."

      },

      action = if (tail_fail) {

        "Inspect lower/upper PIT tails and consider heavier-tailed or asymmetric margins."

      } else {

        "Continue with tail-sensitive diagnostics when tails are substantively important."

      }

    ),

    row(

      area = "Copula fit",

      quantity_checked = "Absolute lag-1 Rosenblatt normal-score residual correlation after fitted copula",

      value = if (is.finite(lag1_cor)) formatC(abs(lag1_cor), digits = 4, format = "fg") else NA_character_,

      threshold_condition = "abs(lag1_cor) > dependence_cor_cutoff",

      default = formatC(dependence_cor_cutoff, digits = 4, format = "fg"),

      status = if (copula_fail) "FAIL" else "PASS",

      message = if (copula_fail) {

        paste0(

          "Dependence remains after the copula in Rosenblatt normal-score residuals (|lag-1 cor| > ",

          dependence_cor_cutoff,

          ")."

        )

      } else {

        "The lag-1 Rosenblatt residual correlation screen did not flag residual dependence."

      },

      action = if (copula_fail) {

        "Consider a different copula family, time-varying dependence, richer serial structure, or a sensitivity refit before treating this as a failure."

      } else {

        "Continue with broader copula diagnostics."

      }

    ),

    row(

      area = "Variance calculation",

      quantity_checked = "Variance-covariance method from summary",

      value = if (is.na(vcov_method) || !nzchar(vcov_method)) NA_character_ else vcov_method,

      threshold_condition = 'vcov_method == "numderiv"',

      default = "n/a",

      status = if (variance_review) "REVIEW" else "PASS",

      message = if (variance_review) {

        "Variance-covariance inference used the numerical Hessian path."

      } else {

        "The variance-covariance method did not trigger the numerical-Hessian review screen."

      },

      action = if (variance_review) {

        "Cite intervals and tests as approximate numerical-Hessian inference."

      } else {

        "Continue with inference checks when reporting intervals or tests."

      }

    )

  ))

}


.gl_basic_checks <- function(checks) {

  checks[, c("area", "status", "value", "threshold_condition", "message"), drop = FALSE]

}


.gl_basic_checks_result <- function(checks) {

  if (any(checks$status == "FAIL", na.rm = TRUE)) {

    return("failed")

  }

  if (any(checks$status == "REVIEW", na.rm = TRUE)) {

    return("review")

  }

  "passed"

}


#' Check a fitted longitudinal GAMLSS-copula model

#'

#' `check_model()` turns diagnostics into a compact set of basic automated

#' checks for broad applied use. It does not replace visual inspection, but it

#' provides a stable first pass over convergence, marginal calibration, residual

#' dependence, and scoring summaries.

#'

#' @param object A fitted `gamlss.longitudinal` object.

#' @param include_vcov Logical; include variance-covariance inference metadata

#'   via [summary.gamlss.longitudinal()].

#' @param include_plots Logical; include standard diagnostic plot objects in

#'   `check$plots`. Visual review should usually use `plot(object)` or the

#'   explicit diagnostic helpers instead.

#' @param dependence_cor_cutoff Absolute lag-1 Rosenblatt normal-score residual

#'   correlation above which the dependence check is flagged. The default is a

#'   review threshold rather than a formal hypothesis test.

#' @param ... Passed to [summary.gamlss.longitudinal()] when `include_vcov` is

#'   `TRUE`.

#'

#' @return An object of class `gamlss_longitudinal_check`, including a compact

#'   `basic_checks` table, a full `checks` table, a `warnings` table containing

#'   failed checks, and overall `basic_checks_passed` and `basic_checks_result`

#'   fields.

#' @export

check_model <- function(

  object,

  include_vcov = FALSE,

  include_plots = FALSE,

  dependence_cor_cutoff = 0.25,

  ...

) {

  if (!inherits(object, "gamlss.longitudinal")) {

    stop("'object' must be a fitted 'gamlss.longitudinal' object.", call. = FALSE)

  }

  dependence_cor_cutoff <- as.numeric(dependence_cor_cutoff)

  if (length(dependence_cor_cutoff) != 1L || !is.finite(dependence_cor_cutoff) ||

      dependence_cor_cutoff <= 0 || dependence_cor_cutoff >= 1) {

    stop("'dependence_cor_cutoff' must be a single number between 0 and 1.", call. = FALSE)

  }


  s <- summary(object, include_vcov = include_vcov, ...)

  pit_out <- .gl_pit(object, randomize = FALSE)

  pit <- pmin(pmax(pit_out$pit, 0), 1)

  ks_p <- tryCatch(stats::ks.test(pit, "punif")$p.value, error = function(e) NA_real_)

  pit_stats <- data.frame(

    n = length(pit),

    mean = mean(pit, na.rm = TRUE),

    sd = stats::sd(pit, na.rm = TRUE),

    expected_sd = sqrt(1 / 12),

    ks_p_value = as.numeric(ks_p),

    stringsAsFactors = FALSE

  )

  tail_summary <- .gl_pit_tail_summary(pit)

  tail_stats <- data.frame(

    tail_ratio_max = attr(tail_summary, "tail_ratio_max"),

    stringsAsFactors = FALSE

  )


  scores <- as.data.frame(as.list(proscore(object, type = c("logs", "mae", "mse", "dss"))), stringsAsFactors = FALSE)

  residual_dependence <- .gl_residual_dependence_summary(object, residual_lags = 1L)

  lag1_cor <- residual_dependence$normal_score_cor[match(1L, residual_dependence$lag)]

  if (length(lag1_cor) == 0L) lag1_cor <- NA_real_

  copula_summary <- tryCatch(copula_time_summary(object), error = function(e) NULL)

  vcov_method <- s$vcov$method %||% s$fit$vcov_method %||% NA_character_

  checks <- .gl_check_table(

    summary_obj = list(fit = s$fit, convergence = object$convergence),

    scores = scores,

    pit_stats = pit_stats,

    tail_stats = tail_stats,

    lag1_cor = lag1_cor,

    dependence_cor_cutoff = dependence_cor_cutoff,

    vcov_method = vcov_method

  )

  warnings <- checks[checks$status == "FAIL", , drop = FALSE]

  basic_checks_result <- .gl_basic_checks_result(checks)

  if (nrow(warnings) > 0L) {

    warning(

      paste0(

        "Basic model checks failed: ",

        paste(warnings$area, collapse = ", "),

        ". Review check$warnings and broader diagnostics."

      ),

      call. = FALSE

    )

  }


  out <- list(

    model = s$model,

    fit = s$fit,

    convergence = object$convergence,

    scores = scores,

    pit = pit_stats,

    tail = tail_summary,

    residual_dependence = transform(residual_dependence, cutoff = dependence_cor_cutoff),

    copula = copula_summary,

    basic_checks = .gl_basic_checks(checks),

    basic_checks_passed = !any(checks$status == "FAIL", na.rm = TRUE),

    basic_checks_result = basic_checks_result,

    checks = checks,

    warnings = warnings,

    plots = if (isTRUE(include_plots)) {

      list(

        pithist = pithist(object, plot = TRUE),

        qqrplot = qqrplot(object, plot = TRUE),

        wormplot = wormplot(object, plot = TRUE),

        rootogram = rootogram(object, plot = TRUE)

      )

    } else {

      NULL

    }

  )

  class(out) <- "gamlss_longitudinal_check"

  out

}


#' @export

print.gamlss_longitudinal_check <- function(x, digits = max(3, getOption("digits") - 3), ...) {

  fmt_num <- function(value) {

    if (length(value) == 0L || is.null(value) || !is.finite(value)) {

      return("n/a")

    }

    formatC(value, digits = digits, format = "fg")

  }

  section <- function(title) {

    cat("\n", title, "\n", paste(rep("-", nchar(title)), collapse = ""), "\n", sep = "")

  }


  cat("\n")

  cat("GAMLSS Longitudinal Model Check\n")

  cat("===============================\n")

  cat("Margin: ", x$model$margin_dist, "    Copula: ", x$model$copula_dist, "\n", sep = "")

  cat("LogLik: ", fmt_num(x$fit$logLik), "    Converged: ", if (isTRUE(x$convergence$converged)) "yes" else "no", "\n", sep = "")


  if (!is.null(x$basic_checks)) {

    section("Basic Checks")

    basic_display <- x$basic_checks[, c("area", "status"), drop = FALSE]

    names(basic_display) <- c("Area", "Status")

    print(basic_display, row.names = FALSE, right = FALSE)

    cat("\nResult: ", toupper(x$basic_checks_result), "\n", sep = "")

    cat("Note: these are basic automated checks; broader model diagnostics should also be reviewed.\n")

  }


  section("Scores")

  print(x$scores, digits = digits, row.names = FALSE)


  section("PIT")

  print(x$pit, digits = digits, row.names = FALSE)


  if (!is.null(x$tail)) {

    section("Tail Calibration")

    print(x$tail, digits = digits, row.names = FALSE)

  }


  section("Residual Dependence")

  print(x$residual_dependence, digits = digits, row.names = FALSE)

  cat("\nUse check$checks for thresholds and check$warnings for failed-check details.\n")

  invisible(x)

}


