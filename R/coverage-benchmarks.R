# Consolidated support module for reviewer navigation.
# Source blocks are grouped mechanically from the files named below.

# ---- coverage-benchmark-setup.R ----

#' @keywords internal

#' @noRd

.coverage_standard_family <- function(family) {
  if (family %in% c("NO", "NO2")) {
    return(stats::gaussian())
  }

  if (family %in% c("PO", "ZIP", "ZIP2")) {
    return(stats::poisson())
  }

  if (family %in% c("GA", "EXP")) {
    return(stats::Gamma(link = "log"))
  }

  NULL
}

#' @keywords internal

#' @noRd

.coverage_standard_formula <- function(design, comparator) {
  if (design %in% c("intercept", "scale", "time_dependence")) {
    return(response ~ 1)
  }

  if (identical(design, "smooth") && identical(comparator, "gam")) {
    return(response ~ s(x, bs = "ps"))
  }

  response ~ x
}

# ---- coverage-benchmark-distributions.R ----

#' @keywords internal

#' @noRd

.coverage_true_margin_distribution <- function(dat, family, p = 0.9) {
  empty <- list(q = rep(NA_real_, nrow(dat)), cdf = rep(NA_real_, nrow(dat)))

  if (is.null(family)) {
    return(empty)
  }

  margin_dist <- tryCatch(do.call(get(family, envir = asNamespace("gamlss.dist")), list()), error = function(e) NULL)

  if (is.null(margin_dist) || is.null(margin_dist$family)) {
    return(empty)
  }

  family_name <- as.character(margin_dist$family[1])

  qfun <- tryCatch(get(paste0("q", family_name), envir = asNamespace("gamlss.dist"), inherits = FALSE), error = function(e) NULL)

  pfun <- tryCatch(get(paste0("p", family_name), envir = asNamespace("gamlss.dist"), inherits = FALSE), error = function(e) NULL)

  if (is.null(qfun) || is.null(pfun)) {
    return(empty)
  }

  par_names <- names(margin_dist$parameters)

  true_cols <- paste0("true_", par_names)

  if (!all(true_cols %in% names(dat))) {
    return(empty)
  }

  par_args <- stats::setNames(lapply(true_cols, function(nm) dat[[nm]]), par_names)

  q <- tryCatch(do.call(qfun, c(list(p = p), par_args)), error = function(e) rep(NA_real_, nrow(dat)))

  cdf <- tryCatch(do.call(pfun, c(list(q = q), par_args)), error = function(e) rep(NA_real_, nrow(dat)))

  list(q = as.numeric(q), cdf = as.numeric(cdf))
}

#' @keywords internal

#' @noRd

.coverage_comparator_dispersion <- function(y, fitted, family) {
  ok <- is.finite(y) & is.finite(fitted)

  if (sum(ok) < 3L) {
    return(NA_real_)
  }

  if (identical(family$family, "gaussian")) {
    return(stats::sd(y[ok] - fitted[ok]))
  }

  if (identical(family$family, "Gamma")) {
    mu <- pmax(fitted[ok], .Machine$double.eps)

    pearson <- (y[ok] - mu) / mu

    return(sum(pearson^2, na.rm = TRUE) / max(1L, sum(ok) - 1L))
  }

  NA_real_
}

#' @keywords internal

#' @noRd

.coverage_comparator_distribution <- function(y, fitted, family, q, p = 0.9) {
  n <- length(fitted)

  empty <- list(q = rep(NA_real_, n), cdf_at_q = rep(NA_real_, n))

  ok <- is.finite(y) & is.finite(fitted) & is.finite(q)

  if (!any(ok)) {
    return(empty)
  }

  if (identical(family$family, "gaussian")) {
    sigma_hat <- .coverage_comparator_dispersion(y, fitted, family)

    if (!is.finite(sigma_hat) || sigma_hat <= 0) {
      return(empty)
    }

    out_q <- fitted + stats::qnorm(p) * sigma_hat

    cdf <- stats::pnorm(q, mean = fitted, sd = sigma_hat)
  } else if (identical(family$family, "poisson")) {
    lambda <- pmax(fitted, .Machine$double.eps)

    out_q <- stats::qpois(p, lambda = lambda)

    cdf <- stats::ppois(q, lambda = lambda)
  } else if (identical(family$family, "Gamma")) {
    dispersion <- .coverage_comparator_dispersion(y, fitted, family)

    if (!is.finite(dispersion) || dispersion <= 0) {
      return(empty)
    }

    mu <- pmax(fitted, .Machine$double.eps)

    shape <- 1 / dispersion

    scale <- mu * dispersion

    out_q <- stats::qgamma(p, shape = shape, scale = scale)

    cdf <- stats::pgamma(q, shape = shape, scale = scale)
  } else {
    return(empty)
  }

  list(
    q = as.numeric(out_q),
    cdf_at_q = pmin(pmax(as.numeric(cdf), 0), 1)
  )
}

# ---- coverage-benchmark-truth-metrics.R ----

#' @keywords internal

#' @noRd

.coverage_benchmark_truth_metrics <- function(dat, fitted, family, gamlss_family = NULL) {
  n <- nrow(dat)

  empty <- c(
    benchmark_mean_bias = NA_real_,
    benchmark_mean_mae = NA_real_,
    benchmark_mean_rmse = NA_real_,
    benchmark_q90_mae = NA_real_,
    benchmark_neg_log_score = NA_real_,
    benchmark_upper_tail_error_90 = NA_real_,
    benchmark_interval_coverage_95 = NA_real_,
    benchmark_interval_width_95 = NA_real_,
    benchmark_pit_ks_p_value = NA_real_,
    benchmark_pit_mean_abs_error = NA_real_,
    benchmark_tail_error_lower_05 = NA_real_,
    benchmark_tail_error_upper_05 = NA_real_
  )

  if (length(fitted) != n) {
    return(empty)
  }

  out <- empty

  if ("true_mu" %in% names(dat)) {
    ok_truth <- is.finite(dat$true_mu) & is.finite(fitted)

    if (any(ok_truth)) {
      err_mu <- fitted[ok_truth] - dat$true_mu[ok_truth]

      out["benchmark_mean_bias"] <- mean(err_mu)

      out["benchmark_mean_mae"] <- mean(abs(err_mu))

      out["benchmark_mean_rmse"] <- sqrt(mean(err_mu^2))
    }
  }

  truth_dist <- .coverage_true_margin_distribution(dat, gamlss_family, p = 0.9)

  comparator_dist <- .coverage_comparator_distribution(dat$response, fitted, family, truth_dist$q, p = 0.9)

  ok_q90 <- is.finite(truth_dist$q) & is.finite(comparator_dist$q)

  if (any(ok_q90)) {
    out["benchmark_q90_mae"] <- mean(abs(comparator_dist$q[ok_q90] - truth_dist$q[ok_q90]), na.rm = TRUE)
  }

  pred_dist <- .benchmark_predictive_distribution(dat$response, fitted, family, p = 0.9)

  ok_density <- is.finite(pred_dist$density) & pred_dist$density > 0

  if (any(ok_density)) {
    out["benchmark_neg_log_score"] <- mean(-log(pmax(pred_dist$density[ok_density], .Machine$double.xmin)), na.rm = TRUE)
  }

  ok_tail <- is.finite(truth_dist$cdf) & is.finite(comparator_dist$cdf_at_q)

  if (any(ok_tail)) {
    true_upper_tail <- 1 - truth_dist$cdf[ok_tail]

    comparator_upper_tail <- 1 - comparator_dist$cdf_at_q[ok_tail]

    out["benchmark_upper_tail_error_90"] <- mean(comparator_upper_tail - true_upper_tail, na.rm = TRUE)
  }

  ok_obs <- is.finite(dat$response) & is.finite(fitted)

  if (!identical(family$family, "gaussian") || sum(ok_obs) < 3L) {
    return(out)
  }

  sigma_hat <- stats::sd(dat$response[ok_obs] - fitted[ok_obs])

  if (!is.finite(sigma_hat) || sigma_hat <= 0) {
    return(out)
  }

  lower <- fitted[ok_obs] + stats::qnorm(0.025) * sigma_hat

  upper <- fitted[ok_obs] + stats::qnorm(0.975) * sigma_hat

  out["benchmark_interval_coverage_95"] <- mean(dat$response[ok_obs] >= lower & dat$response[ok_obs] <= upper)

  out["benchmark_interval_width_95"] <- mean(upper - lower, na.rm = TRUE)

  pit <- stats::pnorm(dat$response[ok_obs], mean = fitted[ok_obs], sd = sigma_hat)

  pit <- pmin(pmax(pit, 0), 1)

  out["benchmark_pit_ks_p_value"] <- tryCatch(

    suppressWarnings(stats::ks.test(pit, "punif")$p.value),
    error = function(e) NA_real_
  )

  out["benchmark_pit_mean_abs_error"] <- abs(mean(pit, na.rm = TRUE) - 0.5)

  out["benchmark_tail_error_lower_05"] <- mean(pit <= 0.05, na.rm = TRUE) - 0.05

  out["benchmark_tail_error_upper_05"] <- mean(pit >= 0.95, na.rm = TRUE) - 0.05

  out
}

# ---- coverage-benchmark-gamlss-metrics.R ----

#' @keywords internal

#' @noRd

.coverage_benchmark_gamlss_metrics <- function(dat, fit, family) {
  n <- nrow(dat)

  empty <- c(
    benchmark_mae = NA_real_,
    benchmark_rmse = NA_real_,
    benchmark_mean_bias = NA_real_,
    benchmark_mean_mae = NA_real_,
    benchmark_mean_rmse = NA_real_,
    benchmark_q90_mae = NA_real_,
    benchmark_neg_log_score = NA_real_,
    benchmark_upper_tail_error_90 = NA_real_,
    benchmark_interval_coverage_95 = NA_real_,
    benchmark_interval_width_95 = NA_real_,
    benchmark_pit_ks_p_value = NA_real_,
    benchmark_pit_mean_abs_error = NA_real_,
    benchmark_tail_error_lower_05 = NA_real_,
    benchmark_tail_error_upper_05 = NA_real_
  )

  if (!inherits(fit, "gamlss.longitudinal")) {
    return(empty)
  }

  tryCatch(
    {
      diag_data <- tryCatch(.gl_fitted_distribution(fit, newdata = NULL, require_response = FALSE), error = function(e) NULL)

      if (is.null(diag_data) || length(diag_data$keep_index) == 0L) {
        return(empty)
      }

      idx <- diag_data$keep_index

      idx <- idx[idx >= 1L & idx <= n]

      if (length(idx) == 0L) {
        return(empty)
      }

      params <- diag_data$params

      mu_hat <- if ("mu" %in% names(params)) as.numeric(params$mu) else as.numeric(params[[1L]])

      fitted <- rep(NA_real_, n)

      fitted[idx] <- mu_hat[seq_along(idx)]

      out <- empty

      ok_obs <- is.finite(dat$response) & is.finite(fitted)

      if (any(ok_obs)) {
        err_obs <- fitted[ok_obs] - dat$response[ok_obs]

        out["benchmark_mae"] <- mean(abs(err_obs))

        out["benchmark_rmse"] <- sqrt(mean(err_obs^2))
      }

      if ("true_mu" %in% names(dat)) {
        ok_truth <- is.finite(dat$true_mu) & is.finite(fitted)

        if (any(ok_truth)) {
          err_mu <- fitted[ok_truth] - dat$true_mu[ok_truth]

          out["benchmark_mean_bias"] <- mean(err_mu)

          out["benchmark_mean_mae"] <- mean(abs(err_mu))

          out["benchmark_mean_rmse"] <- sqrt(mean(err_mu^2))
        }
      }

      truth_dist <- .coverage_true_margin_distribution(dat, family, p = 0.9)

      fitted_q90 <- rep(NA_real_, n)

      fitted_cdf_at_true_q90 <- rep(NA_real_, n)

      fitted_q025 <- rep(NA_real_, n)

      fitted_q975 <- rep(NA_real_, n)

      fitted_pit <- rep(NA_real_, n)

      fitted_density <- rep(NA_real_, n)

      fitted_q90[idx] <- tryCatch(

        as.numeric(.gl_call_family_fun("q", diag_data$family, 0.9, params))[seq_along(idx)],
        error = function(e) rep(NA_real_, length(idx))
      )

      fitted_cdf_at_true_q90[idx] <- tryCatch(

        as.numeric(.gl_call_family_fun("p", diag_data$family, truth_dist$q[idx], params))[seq_along(idx)],
        error = function(e) rep(NA_real_, length(idx))
      )

      fitted_q025[idx] <- tryCatch(

        as.numeric(.gl_call_family_fun("q", diag_data$family, 0.025, params))[seq_along(idx)],
        error = function(e) rep(NA_real_, length(idx))
      )

      fitted_q975[idx] <- tryCatch(

        as.numeric(.gl_call_family_fun("q", diag_data$family, 0.975, params))[seq_along(idx)],
        error = function(e) rep(NA_real_, length(idx))
      )

      fitted_pit[idx] <- tryCatch(

        as.numeric(.gl_call_family_fun("p", diag_data$family, dat$response[idx], params))[seq_along(idx)],
        error = function(e) rep(NA_real_, length(idx))
      )

      fitted_density[idx] <- tryCatch(

        as.numeric(.gl_call_family_fun("d", diag_data$family, dat$response[idx], params))[seq_along(idx)],
        error = function(e) rep(NA_real_, length(idx))
      )

      ok_q90 <- is.finite(truth_dist$q) & is.finite(fitted_q90)

      if (any(ok_q90)) {
        out["benchmark_q90_mae"] <- mean(abs(fitted_q90[ok_q90] - truth_dist$q[ok_q90]), na.rm = TRUE)
      }

      ok_density <- is.finite(fitted_density) & fitted_density > 0

      if (any(ok_density)) {
        out["benchmark_neg_log_score"] <- mean(-log(pmax(fitted_density[ok_density], .Machine$double.xmin)), na.rm = TRUE)
      }

      ok_tail <- is.finite(truth_dist$cdf) & is.finite(fitted_cdf_at_true_q90)

      if (any(ok_tail)) {
        true_upper_tail <- 1 - truth_dist$cdf[ok_tail]

        fitted_upper_tail <- 1 - fitted_cdf_at_true_q90[ok_tail]

        out["benchmark_upper_tail_error_90"] <- mean(fitted_upper_tail - true_upper_tail, na.rm = TRUE)
      }

      ok_interval <- is.finite(dat$response) & is.finite(fitted_q025) & is.finite(fitted_q975)

      if (any(ok_interval)) {
        out["benchmark_interval_coverage_95"] <- mean(
          dat$response[ok_interval] >= fitted_q025[ok_interval] &

            dat$response[ok_interval] <= fitted_q975[ok_interval]
        )

        out["benchmark_interval_width_95"] <- mean(fitted_q975[ok_interval] - fitted_q025[ok_interval], na.rm = TRUE)
      }

      ok_pit <- is.finite(fitted_pit)

      if (sum(ok_pit) >= 3L) {
        pit <- pmin(pmax(fitted_pit[ok_pit], 0), 1)

        out["benchmark_pit_ks_p_value"] <- tryCatch(

          suppressWarnings(stats::ks.test(pit, "punif")$p.value),
          error = function(e) NA_real_
        )

        out["benchmark_pit_mean_abs_error"] <- abs(mean(pit, na.rm = TRUE) - 0.5)

        out["benchmark_tail_error_lower_05"] <- mean(pit <= 0.05, na.rm = TRUE) - 0.05

        out["benchmark_tail_error_upper_05"] <- mean(pit >= 0.95, na.rm = TRUE) - 0.05
      }

      out
    },
    error = function(e) empty
  )
}
