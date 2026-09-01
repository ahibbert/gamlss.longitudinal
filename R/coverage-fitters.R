# Consolidated support module for reviewer navigation.
# Source blocks are grouped mechanically from the files named below.

# ---- coverage-fitters-longitudinal-start.R ----

#' @keywords internal

#' @noRd

.coverage_longitudinal_start_fit <- function(
    dat,
    family,
    copula,
    design,
    max_outer_iter = 5,
    max_inner_iter = 8,
    max_elapsed_sec = 20) {
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
      missingness = "segment",
      mu.formula = formulas$mu,
      sigma.formula = formulas$sigma,
      nu.formula = formulas$nu,
      tau.formula = formulas$tau,
      theta.formula = formulas$theta,
      zeta.formula = formulas$zeta,
      include_dlcopdpar = FALSE,
      method = "RS",
      start_from = NA,
      optimizer_control = gamlss_longitudinal_control(
        outer_tol = 1e-4,
        max_outer_iter = max_outer_iter,
        max_elapsed_sec = max_elapsed_sec,
        rs = list(
          inner_tol = 1e-4,
          max_inner_iter = max_inner_iter,
          warm_start_joint = FALSE,
          warm_start_joint_iter = 0L
        )
      ),
      compute_vcov = FALSE,
      verbose = 0
    )
  }))

  fit
}

# ---- coverage-fitters-longitudinal-attempt.R ----

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
    attempt_label = "default") {
  .coverage_attach_namespace("gamlss")

  .coverage_attach_namespace("gamlss.dist")

  margin_dist <- do.call(get(family, envir = asNamespace("gamlss.dist")), list())

  formulas <- .coverage_fit_formulas(design)

  fit_method <- if (identical(method, "cg")) "CG" else "RS"

  include_dlcopdpar <- identical(method, "rs_joint") || identical(method, "cg")
  optimizer_control <- if (identical(fit_method, "CG")) {
    gamlss_longitudinal_control(
      outer_tol = 1e-4,
      max_outer_iter = max_outer_iter,
      max_elapsed_sec = max_elapsed_sec,
      cg = list(max_delta = cg_max_delta)
    )
  } else {
    gamlss_longitudinal_control(
      outer_tol = 1e-4,
      max_outer_iter = max_outer_iter,
      max_elapsed_sec = max_elapsed_sec,
      rs = list(
        inner_tol = 1e-4,
        max_inner_iter = max_inner_iter,
        warm_start_joint = isTRUE(warm_start_joint) && isTRUE(include_dlcopdpar) && all(is.na(start_from)),
        warm_start_joint_iter = min(5L, max_outer_iter),
        start_step_size = start_step_size
      )
    )
  }

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
        missingness = "segment",
        mu.formula = formulas$mu,
        sigma.formula = formulas$sigma,
        nu.formula = formulas$nu,
        tau.formula = formulas$tau,
        theta.formula = formulas$theta,
        zeta.formula = formulas$zeta,
        include_dlcopdpar = include_dlcopdpar,
        method = fit_method,
        start_from = start_from,
        optimizer_control = optimizer_control,
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

# ---- coverage-fitters-longitudinal-selection.R ----

#' @keywords internal
#' @noRd
.coverage_longitudinal_attempt_success <- function(attempt) {
  fit <- attempt$value
  is_fit <- inherits(fit, "gamlss.longitudinal")
  if (!is_fit) {
    return(FALSE)
  }

  loglik <- fit$calc_lik_out_end$log_lik
  all(is.finite(as.numeric(loglik)))
}

#' @keywords internal
#' @noRd
.coverage_longitudinal_successful_attempts <- function(attempts) {
  vapply(attempts, .coverage_longitudinal_attempt_success, logical(1))
}

#' @keywords internal
#' @noRd
.coverage_longitudinal_chosen_index <- function(attempts) {
  successful <- .coverage_longitudinal_successful_attempts(attempts)
  if (any(successful)) {
    which(successful)[1L]
  } else {
    length(attempts)
  }
}

#' @keywords internal
#' @noRd
.coverage_longitudinal_elapsed <- function(attempts, chosen_idx) {
  sum(vapply(attempts[seq_len(chosen_idx)], `[[`, numeric(1), "elapsed"))
}

#' @keywords internal
#' @noRd
.coverage_longitudinal_attempt_warnings <- function(attempts, chosen_idx) {
  warnings <- unique(unlist(lapply(attempts[seq_len(chosen_idx)], `[[`, "warnings")))
  paste(warnings, collapse = " | ")
}

#' @keywords internal
#' @noRd
.coverage_longitudinal_attempt_trace <- function(attempts, chosen_idx) {
  paste(
    vapply(attempts[seq_len(chosen_idx)], `[[`, character(1), "attempt_label"),
    collapse = " > "
  )
}

# ---- coverage-fitters-longitudinal-result.R ----

#' @keywords internal
#' @noRd
.coverage_longitudinal_result_row <- function(
    method,
    success,
    is_fit,
    fit,
    chosen,
    chosen_idx,
    attempts,
    elapsed,
    max_elapsed_sec,
    loglik,
    abs_error,
    rel_error,
    theta_tau,
    truth,
    smooth_metrics,
    truth_metrics) {
  data.frame(
    method = method,
    success = success,
    converged = if (is_fit && !is.null(fit$convergence$converged)) isTRUE(fit$convergence$converged) else FALSE,
    failure_type = .coverage_taxonomy(success, if (inherits(fit, "error")) conditionMessage(fit) else NA_character_, chosen$warnings, elapsed, max_elapsed_sec),
    error = if (inherits(fit, "error")) conditionMessage(fit) else NA_character_,
    warnings = .coverage_longitudinal_attempt_warnings(attempts, chosen_idx),
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
    fit_attempt_trace = .coverage_longitudinal_attempt_trace(attempts, chosen_idx),
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
}

# ---- coverage-fitters-longitudinal.R ----

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
    start_from = NA) {
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

  chosen_idx <- .coverage_longitudinal_chosen_index(attempts)

  chosen <- attempts[[chosen_idx]]

  elapsed <- .coverage_longitudinal_elapsed(attempts, chosen_idx)

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

  out <- .coverage_longitudinal_result_row(
    method = method,
    success = success,
    is_fit = is_fit,
    fit = fit,
    chosen = chosen,
    chosen_idx = chosen_idx,
    attempts = attempts,
    elapsed = elapsed,
    max_elapsed_sec = max_elapsed_sec,
    loglik = loglik,
    abs_error = abs_error,
    rel_error = rel_error,
    theta_tau = theta_tau,
    truth = truth,
    smooth_metrics = smooth_metrics,
    truth_metrics = truth_metrics
  )

  parameter_results <- .coverage_parameter_results(
    .coverage_longitudinal_eta_estimates(fit),
    .coverage_true_eta_coefficients(dat, family, copula, design)
  )

  theta_time <- parameter_results[
    parameter_results$parameter == "theta" &

      parameter_results$term == "time_covariate", ,
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

# ---- coverage-fitters-gamlss.R ----

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

# ---- coverage-fitters-gamlss2.R ----

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

# ---- coverage-fitters-standard.R ----

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
