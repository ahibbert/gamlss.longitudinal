.coverage_fit_gamlss <- function(dat, family, copula, design, max_inner_iter = 8, max_elapsed_sec = Inf) {

  .coverage_attach_namespace("gamlss")

  .coverage_attach_namespace("gamlss.dist")

  margin_dist <- do.call(get(family, envir = asNamespace("gamlss.dist")), list())

  formulas <- .coverage_fit_formulas(design)

  fml <- formulas$mu

  fit_vars <- unique(unlist(lapply(formulas[names(formulas) %in% names(margin_dist$parameters)], all.vars)))

  fit_vars <- fit_vars[fit_vars %in% names(dat)]

  fit_dat <- dat[, fit_vars, drop = FALSE]

  start <- Sys.time()

  captured <- .coverage_capture_conditions({

    fit <- NULL

    invisible(utils::capture.output({

      fit <- gamlss::gamlss(

        formula = fml,

        sigma.formula = formulas$sigma,

        nu.formula = formulas$nu,

        tau.formula = formulas$tau,

        data = fit_dat,

        family = margin_dist,

        control = gamlss::gamlss.control(n.cyc = max_inner_iter, trace = FALSE)

      )

    }))

    fit

  })

  elapsed <- as.numeric(difftime(Sys.time(), start, units = "secs"))

  success <- !inherits(captured$value, "error") &&

    is.finite(as.numeric(stats::logLik(captured$value)))

  converged <- success

  if (success && !is.null(captured$value$converged)) {

    converged <- isTRUE(captured$value$converged)

  }

  out <- data.frame(

    method = "gamlss",

    success = success,

    converged = converged,

    failure_type = .coverage_taxonomy(success, if (inherits(captured$value, "error")) conditionMessage(captured$value) else NA_character_, captured$warnings, elapsed, max_elapsed_sec),

    error = if (inherits(captured$value, "error")) conditionMessage(captured$value) else NA_character_,

    warnings = paste(unique(captured$warnings), collapse = " | "),

    marginal_loglik = if (success) as.numeric(stats::logLik(captured$value)) else NA_real_,

    copula_loglik = NA_real_,

    joint_loglik = NA_real_,

    elapsed_sec = elapsed,

    max_abs_param_error = NA_real_,

    max_rel_param_error = NA_real_,

    fitted_copula_tau = NA_real_,

    true_copula_tau = NA_real_,

    marginal_fit_method = "gamlss",

    marginal_fallback_used = FALSE,

    marginal_fallback_error = NA_character_,

    fit_attempt = "gamlss",

    fit_attempt_count = 1L,

    fit_attempt_trace = "gamlss",

    stringsAsFactors = FALSE

  )

  truth_eta <- .coverage_true_eta_coefficients(dat, family, copula, design)

  truth_eta <- truth_eta[truth_eta$parameter %in% names(margin_dist$parameters), , drop = FALSE]

  attr(out, "parameter_results") <- .coverage_parameter_results(

    .coverage_gamlss_eta_estimates(captured$value, margin_dist),

    truth_eta

  )

  out

}


#' @keywords internal

#' @noRd

.coverage_fit_gamlss2 <- function(dat, family, copula, design, max_elapsed_sec = Inf) {

  .coverage_attach_namespace("gamlss.dist")

  margin_dist <- do.call(get(family, envir = asNamespace("gamlss.dist")), list())

  fml <- if (design %in% c("intercept", "scale", "time_dependence")) response ~ 1 else response ~ x

  fit_vars <- unique(all.vars(fml))

  fit_vars <- fit_vars[fit_vars %in% names(dat)]

  fit_dat <- dat[, fit_vars, drop = FALSE]

  start <- Sys.time()

  fallback_used <- FALSE

  fallback_error <- NA_character_

  gamlss2_available <- requireNamespace("gamlss2", quietly = TRUE)

  captured <- if (gamlss2_available) {

    gamlss2_fit <- getExportedValue("gamlss2", "gamlss2")

    gamlss2_control <- getExportedValue("gamlss2", "gamlss2_control")

    .coverage_capture_conditions({

      fit <- gamlss2_fit(

        fml,

        data = fit_dat,

        family = margin_dist,

        control = gamlss2_control(trace = FALSE)

      )

      fit

    })

  } else {

    list(value = simpleError("Package 'gamlss2' is required for method = 'gamlss2'."), warnings = character(0))

  }


  success <- !inherits(captured$value, "error") &&

    is.finite(as.numeric(stats::logLik(captured$value)))

  if (!success) {

    fallback_error <- if (inherits(captured$value, "error")) conditionMessage(captured$value) else "gamlss2 returned non-finite logLik"

    fallback_used <- TRUE

    gamlss_captured <- .coverage_capture_conditions({

      fit <- NULL

      invisible(utils::capture.output({

        fit <- gamlss::gamlss(

          formula = fml,

          data = fit_dat,

          family = margin_dist,

          control = gamlss::gamlss.control(n.cyc = 8, trace = FALSE)

        )

      }))

      fit

    })

    if (!inherits(gamlss_captured$value, "error") && is.finite(as.numeric(stats::logLik(gamlss_captured$value)))) {

      captured <- gamlss_captured

      success <- TRUE

    } else {

      captured$warnings <- c(captured$warnings, gamlss_captured$warnings)

      captured$value <- if (inherits(gamlss_captured$value, "error")) gamlss_captured$value else captured$value

    }

  }

  elapsed <- as.numeric(difftime(Sys.time(), start, units = "secs"))

  success <- success && (!is.finite(max_elapsed_sec) || elapsed <= max_elapsed_sec)

  actual_method <- if (success && isTRUE(fallback_used)) "gamlss" else "gamlss2"

  out <- data.frame(

    method = "gamlss2",

    success = success,

    converged = success,

    failure_type = .coverage_taxonomy(success, if (inherits(captured$value, "error")) conditionMessage(captured$value) else NA_character_, captured$warnings, elapsed, max_elapsed_sec),

    error = if (inherits(captured$value, "error")) conditionMessage(captured$value) else NA_character_,

    warnings = paste(unique(captured$warnings), collapse = " | "),

    marginal_loglik = if (success) as.numeric(stats::logLik(captured$value)) else NA_real_,

    copula_loglik = NA_real_,

    joint_loglik = NA_real_,

    elapsed_sec = elapsed,

    max_abs_param_error = NA_real_,

    max_rel_param_error = NA_real_,

    fitted_copula_tau = NA_real_,

    true_copula_tau = NA_real_,

    marginal_fit_method = actual_method,

    marginal_fallback_used = fallback_used && success,

    marginal_fallback_error = fallback_error,

    fit_attempt = actual_method,

    fit_attempt_count = if (isTRUE(fallback_used)) 2L else 1L,

    fit_attempt_trace = if (isTRUE(fallback_used)) "gamlss2 > gamlss" else "gamlss2",

    stringsAsFactors = FALSE

  )

  truth_eta <- .coverage_true_eta_coefficients(dat, family, copula, design)

  truth_eta <- truth_eta[truth_eta$parameter %in% names(margin_dist$parameters), , drop = FALSE]

  estimates <- if (identical(actual_method, "gamlss")) {

    .coverage_gamlss_eta_estimates(captured$value, margin_dist)

  } else {

    .coverage_gamlss2_eta_estimates(captured$value)

  }

  attr(out, "parameter_results") <- .coverage_parameter_results(estimates, truth_eta)

  out

}


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


  tryCatch({

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

  }, error = function(e) empty)

}


#' @keywords internal

#' @noRd

.coverage_fit_standard_comparator <- function(dat, family, copula, design, method, max_elapsed_sec = Inf) {

  comparator_family <- .coverage_standard_family(family)

  start <- Sys.time()

  truth <- .coverage_truth_summary(dat, copula)

  true_tau <- truth["copula_tau"] %||% NA_real_

  empty_row <- function(success, failure_type, elapsed, error = NA_character_, warning = NA_character_) {

    data.frame(

      method = method,

      success = success,

      failure_type = failure_type,

      elapsed_sec = elapsed,

      max_abs_error = NA_real_,

      max_rel_error = NA_real_,

      fitted_copula_tau = NA_real_,

      true_copula_tau = unname(true_tau),

      marginal_loglik = NA_real_,

      copula_loglik = NA_real_,

      joint_loglik = NA_real_,

      marginal_fit_method = method,

      fit_attempt = method,

      fit_attempt_trace = method,

      benchmark_comparator = method,

      benchmark_class = NA_character_,

      benchmark_estimator = NA_character_,

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

      benchmark_tail_error_upper_05 = NA_real_,

      benchmark_error = error,

      benchmark_warning = warning,

      stringsAsFactors = FALSE

    )

  }


  if (is.null(comparator_family)) {

    elapsed <- as.numeric(difftime(Sys.time(), start, units = "secs"))

    return(empty_row(

      success = FALSE,

      failure_type = "unsupported comparator family",

      elapsed = elapsed,

      error = paste0("No standard comparator family mapping for GAMLSS family '", family, "'.")

    ))

  }


  dat_fit <- dat

  dat_fit$subject <- factor(dat_fit$subject)

  formula <- .coverage_standard_formula(design, method)

  captured <- .coverage_capture_conditions({

    benchmark_standard_models(

      data = dat_fit,

      formula = formula,

      subject_var = "subject",

      family = comparator_family,

      comparators = method,

      correlation = "exchangeable"

    )

  })

  elapsed <- as.numeric(difftime(Sys.time(), start, units = "secs"))

  if (inherits(captured$value, "error")) {

    return(empty_row(

      success = FALSE,

      failure_type = .coverage_taxonomy(FALSE, conditionMessage(captured$value), captured$warnings, elapsed, max_elapsed_sec),

      elapsed = elapsed,

      error = conditionMessage(captured$value),

      warning = if (length(captured$warnings)) paste(unique(captured$warnings), collapse = " | ") else NA_character_

    ))

  }


  row <- captured$value$results[1L, , drop = FALSE]

  success <- isTRUE(row$success) && (!is.finite(max_elapsed_sec) || elapsed <= max_elapsed_sec)

  failure_type <- if (success) {

    "ok"

  } else if (!isTRUE(row$available)) {

    "comparator unavailable"

  } else {

    .coverage_taxonomy(FALSE, row$error, row$warning, elapsed, max_elapsed_sec)

  }

  fit <- captured$value$fits[[method]]

  fitted <- if (success && !is.null(fit)) {

    .benchmark_predict_response(fit, dat_fit)

  } else {

    rep(NA_real_, nrow(dat_fit))

  }

  truth_metrics <- .coverage_benchmark_truth_metrics(dat_fit, fitted, comparator_family, gamlss_family = family)

  data.frame(

    method = method,

    success = success,

    failure_type = failure_type,

    elapsed_sec = elapsed,

    max_abs_error = row$mae,

    max_rel_error = NA_real_,

    fitted_copula_tau = NA_real_,

    true_copula_tau = unname(true_tau),

    marginal_loglik = row$logLik,

    copula_loglik = NA_real_,

    joint_loglik = row$logLik,

    marginal_fit_method = row$estimator,

    fit_attempt = method,

    fit_attempt_trace = row$estimator,

    benchmark_comparator = row$comparator,

    benchmark_class = row$comparator_class,

    benchmark_estimator = row$estimator,

    benchmark_mae = row$mae,

    benchmark_rmse = row$rmse,

    benchmark_mean_bias = unname(truth_metrics[["benchmark_mean_bias"]]),

    benchmark_mean_mae = unname(truth_metrics[["benchmark_mean_mae"]]),

    benchmark_mean_rmse = unname(truth_metrics[["benchmark_mean_rmse"]]),

    benchmark_q90_mae = unname(truth_metrics[["benchmark_q90_mae"]]),

    benchmark_neg_log_score = unname(truth_metrics[["benchmark_neg_log_score"]]),

    benchmark_upper_tail_error_90 = unname(truth_metrics[["benchmark_upper_tail_error_90"]]),

    benchmark_interval_coverage_95 = unname(truth_metrics[["benchmark_interval_coverage_95"]]),

    benchmark_interval_width_95 = unname(truth_metrics[["benchmark_interval_width_95"]]),

    benchmark_pit_ks_p_value = unname(truth_metrics[["benchmark_pit_ks_p_value"]]),

    benchmark_pit_mean_abs_error = unname(truth_metrics[["benchmark_pit_mean_abs_error"]]),

    benchmark_tail_error_lower_05 = unname(truth_metrics[["benchmark_tail_error_lower_05"]]),

    benchmark_tail_error_upper_05 = unname(truth_metrics[["benchmark_tail_error_upper_05"]]),

    benchmark_error = row$error,

    benchmark_warning = row$warning,

    stringsAsFactors = FALSE

  )

}


#' @keywords internal

#' @noRd

.coverage_longitudinal_start_fit <- function(

  dat,

  family,

  copula,

  design,

  max_outer_iter = 5,

  max_inner_iter = 8,

  max_elapsed_sec = 20

) {

  .coverage_attach_namespace("gamlss")

  .coverage_attach_namespace("gamlss.dist")

  margin_dist <- do.call(get(family, envir = asNamespace("gamlss.dist")), list())

  formulas <- .coverage_fit_formulas(design)

  fit <- NULL

  invisible(utils::capture.output({

    fit <- gamlss_longitudinal(

      dataset = dat,

      margin_dist = margin_dist,

      copula_dist = copula,

      time_var = "time",

      subject_var = "subject",

      mu.formula = formulas$mu,

      sigma.formula = formulas$sigma,

      nu.formula = formulas$nu,

      tau.formula = formulas$tau,

      theta.formula = formulas$theta,

      zeta.formula = formulas$zeta,

      include_dlcopdpar = FALSE,

      method = "RS",

      start_from = NA,

      warm_start_joint = FALSE,

      warm_start_joint_iter = 0L,

      max_outer_iter = max_outer_iter,

      max_inner_iter = max_inner_iter,

      outer_stop_crit = 1e-4,

      inner_stop_crit = 1e-4,

      max_elapsed_sec = max_elapsed_sec,

      compute_vcov = FALSE,

      verbose = 0

    )

  }))

  fit

}


#' @keywords internal

#' @noRd

.coverage_longitudinal_attempt <- function(

  dat,

  family,

  copula,

  method,

  design,

  max_outer_iter = 8,

  max_inner_iter = 8,

  max_elapsed_sec = 20,

  start_from = NA,

  warm_start_joint = TRUE,

  start_step_size = 0.5,

  cg_max_delta = 0.5,

  attempt_label = "default"

) {

  .coverage_attach_namespace("gamlss")

  .coverage_attach_namespace("gamlss.dist")

  margin_dist <- do.call(get(family, envir = asNamespace("gamlss.dist")), list())

  formulas <- .coverage_fit_formulas(design)

  fit_method <- if (identical(method, "cg")) "CG" else "RS"

  include_dlcopdpar <- identical(method, "rs_joint") || identical(method, "cg")

  start <- Sys.time()

  captured <- .coverage_capture_conditions({

    fit <- NULL

    invisible(utils::capture.output({

      fit <- gamlss_longitudinal(

        dataset = dat,

        margin_dist = margin_dist,

        copula_dist = copula,

        time_var = "time",

        subject_var = "subject",

        mu.formula = formulas$mu,

        sigma.formula = formulas$sigma,

        nu.formula = formulas$nu,

        tau.formula = formulas$tau,

        theta.formula = formulas$theta,

        zeta.formula = formulas$zeta,

        include_dlcopdpar = include_dlcopdpar,

        method = fit_method,

        start_from = start_from,

        warm_start_joint = isTRUE(warm_start_joint) && isTRUE(include_dlcopdpar) && all(is.na(start_from)),

        warm_start_joint_iter = min(5L, max_outer_iter),

        start_step_size = start_step_size,

        cg_max_delta = cg_max_delta,

        max_outer_iter = max_outer_iter,

        max_inner_iter = max_inner_iter,

        outer_stop_crit = 1e-4,

        inner_stop_crit = 1e-4,

        max_elapsed_sec = max_elapsed_sec,

        compute_vcov = FALSE,

        verbose = 0

      )

    }))

    fit

  })

  elapsed <- as.numeric(difftime(Sys.time(), start, units = "secs"))

  list(

    value = captured$value,

    warnings = captured$warnings,

    elapsed = elapsed,

    attempt_label = attempt_label

  )

}


#' @keywords internal

#' @noRd

.coverage_fit_longitudinal <- function(

  dat,

  family,

  copula,

  method,

  design,

  max_outer_iter = 8,

  max_inner_iter = 8,

  max_elapsed_sec = 20,

  start_from = NA

) {

  margin_dist <- do.call(get(family, envir = asNamespace("gamlss.dist")), list())

  attempts <- list()

  start_from_supplied <- !all(is.na(start_from))


  if (identical(method, "cg") && !start_from_supplied) {

    rs_start <- tryCatch(

      .coverage_longitudinal_start_fit(

        dat,

        family = family,

        copula = copula,

        design = design,

        max_outer_iter = min(5L, max_outer_iter),

        max_inner_iter = max_inner_iter,

        max_elapsed_sec = max_elapsed_sec

      )$par,

      error = function(e) NA

    )

    attempts[[length(attempts) + 1L]] <- .coverage_longitudinal_attempt(

      dat, family, copula, method, design,

      max_outer_iter = max_outer_iter,

      max_inner_iter = max_inner_iter,

      max_elapsed_sec = max_elapsed_sec,

      start_from = rs_start,

      warm_start_joint = FALSE,

      attempt_label = "rs_separate_start"

    )

  }


  attempts[[length(attempts) + 1L]] <- .coverage_longitudinal_attempt(

    dat, family, copula, method, design,

    max_outer_iter = max_outer_iter,

    max_inner_iter = max_inner_iter,

    max_elapsed_sec = max_elapsed_sec,

    start_from = start_from,

    warm_start_joint = TRUE,

    attempt_label = if (start_from_supplied) "explicit_start" else "default"

  )


  if (identical(method, "cg")) {

    attempts[[length(attempts) + 1L]] <- .coverage_longitudinal_attempt(

      dat, family, copula, method, design,

      max_outer_iter = max_outer_iter,

      max_inner_iter = max_inner_iter,

      max_elapsed_sec = max_elapsed_sec,

      start_from = NA,

      warm_start_joint = FALSE,

      attempt_label = "cold_start"

    )

    interior_start <- .coverage_truth_adjacent_start(dat, family, copula, design)

    attempts[[length(attempts) + 1L]] <- .coverage_longitudinal_attempt(

      dat, family, copula, method, design,

      max_outer_iter = max_outer_iter,

      max_inner_iter = max_inner_iter,

      max_elapsed_sec = max_elapsed_sec,

      start_from = interior_start,

      warm_start_joint = FALSE,

      start_step_size = 0.2,

      cg_max_delta = 0.2,

      attempt_label = "damped_interior_start"

    )

  }


  successful <- vapply(attempts, function(x) {

    fit <- x$value

    is_fit <- inherits(fit, "gamlss.longitudinal")

    if (!is_fit) return(FALSE)

    loglik <- fit$calc_lik_out_end$log_lik

    all(is.finite(as.numeric(loglik)))

  }, logical(1))

  chosen_idx <- if (any(successful)) which(successful)[1L] else length(attempts)

  chosen <- attempts[[chosen_idx]]

  elapsed <- sum(vapply(attempts[seq_len(chosen_idx)], `[[`, numeric(1), "elapsed"))

  fit <- chosen$value

  is_fit <- inherits(fit, "gamlss.longitudinal")

  loglik <- if (is_fit) fit$calc_lik_out_end$log_lik else c(marginal = NA_real_, copula = NA_real_, joint = NA_real_)

  estimates <- if (is_fit) .coverage_natural_estimates(fit) else stats::setNames(numeric(0), character(0))

  truth <- .coverage_truth_summary(dat, copula)

  common <- intersect(names(estimates), names(truth))

  abs_error <- if (length(common)) max(abs(estimates[common] - truth[common]), na.rm = TRUE) else NA_real_

  rel_error <- if (length(common)) {

    denom <- pmax(abs(truth[common]), 1e-8)

    max(abs(estimates[common] - truth[common]) / denom, na.rm = TRUE)

  } else {

    NA_real_

  }

  theta_tau <- if (is_fit && "theta" %in% names(estimates)) {

    zeta <- if ("zeta" %in% names(estimates)) estimates[["zeta"]] else 0

    .copula_par_to_tau(copula, estimates[["theta"]], zeta)

  } else {

    NA_real_

  }

  success <- is_fit &&

    all(is.finite(as.numeric(loglik))) &&

    !inherits(fit, "error") &&

    (!is.finite(max_elapsed_sec) || elapsed <= max_elapsed_sec)

  truth_metrics <- if (success) {

    .coverage_benchmark_gamlss_metrics(dat, fit, family)

  } else {

    .coverage_benchmark_gamlss_metrics(dat, NULL, family)

  }

  smooth_metrics <- if (success && identical(design, "smooth")) {

    .coverage_smooth_eta_recovery(fit, dat, copula)

  } else {

    c(smooth_eta_rmse = NA_real_, smooth_eta_max_abs_error = NA_real_, smooth_eta_n = 0)

  }

  out <- data.frame(

    method = method,

    success = success,

    converged = if (is_fit && !is.null(fit$convergence$converged)) isTRUE(fit$convergence$converged) else FALSE,

    failure_type = .coverage_taxonomy(success, if (inherits(fit, "error")) conditionMessage(fit) else NA_character_, chosen$warnings, elapsed, max_elapsed_sec),

    error = if (inherits(fit, "error")) conditionMessage(fit) else NA_character_,

    warnings = paste(unique(unlist(lapply(attempts[seq_len(chosen_idx)], `[[`, "warnings"))), collapse = " | "),

    marginal_loglik = as.numeric(loglik["marginal"]),

    copula_loglik = as.numeric(loglik["copula"]),

    joint_loglik = as.numeric(loglik["joint"]),

    elapsed_sec = elapsed,

    max_abs_param_error = abs_error,

    max_rel_param_error = rel_error,

    fitted_copula_tau = theta_tau,

    true_copula_tau = if ("copula_tau" %in% names(truth)) truth[["copula_tau"]] else NA_real_,

    marginal_fit_method = method,

    marginal_fallback_used = FALSE,

    marginal_fallback_error = NA_character_,

    fit_attempt = chosen$attempt_label,

    fit_attempt_count = chosen_idx,

    fit_attempt_trace = paste(vapply(attempts[seq_len(chosen_idx)], `[[`, character(1), "attempt_label"), collapse = " > "),

    smooth_eta_rmse = unname(smooth_metrics[["smooth_eta_rmse"]]),

    smooth_eta_max_abs_error = unname(smooth_metrics[["smooth_eta_max_abs_error"]]),

    smooth_eta_n = unname(smooth_metrics[["smooth_eta_n"]]),

    benchmark_comparator = "gamlss.longitudinal",

    benchmark_class = "gamlss_longitudinal",

    benchmark_estimator = paste0("gamlss.longitudinal::", method),

    benchmark_mae = unname(truth_metrics[["benchmark_mae"]]),

    benchmark_rmse = unname(truth_metrics[["benchmark_rmse"]]),

    benchmark_mean_bias = unname(truth_metrics[["benchmark_mean_bias"]]),

    benchmark_mean_mae = unname(truth_metrics[["benchmark_mean_mae"]]),

    benchmark_mean_rmse = unname(truth_metrics[["benchmark_mean_rmse"]]),

    benchmark_q90_mae = unname(truth_metrics[["benchmark_q90_mae"]]),

    benchmark_neg_log_score = unname(truth_metrics[["benchmark_neg_log_score"]]),

    benchmark_upper_tail_error_90 = unname(truth_metrics[["benchmark_upper_tail_error_90"]]),

    benchmark_interval_coverage_95 = unname(truth_metrics[["benchmark_interval_coverage_95"]]),

    benchmark_interval_width_95 = unname(truth_metrics[["benchmark_interval_width_95"]]),

    benchmark_pit_ks_p_value = unname(truth_metrics[["benchmark_pit_ks_p_value"]]),

    benchmark_pit_mean_abs_error = unname(truth_metrics[["benchmark_pit_mean_abs_error"]]),

    benchmark_tail_error_lower_05 = unname(truth_metrics[["benchmark_tail_error_lower_05"]]),

    benchmark_tail_error_upper_05 = unname(truth_metrics[["benchmark_tail_error_upper_05"]]),

    benchmark_error = NA_character_,

    benchmark_warning = NA_character_,

    stringsAsFactors = FALSE

  )

  parameter_results <- .coverage_parameter_results(

    .coverage_longitudinal_eta_estimates(fit),

    .coverage_true_eta_coefficients(dat, family, copula, design)

  )

  theta_time <- parameter_results[

    parameter_results$parameter == "theta" &

      parameter_results$term == "time_covariate",

    ,

    drop = FALSE

  ]

  out$benchmark_theta_time_abs_error <- if (nrow(theta_time) == 1L) {

    theta_time$abs_eta_error[[1L]]

  } else {

    NA_real_

  }

  attr(out, "parameter_results") <- parameter_results

  out

}


#' @keywords internal

#' @noRd
