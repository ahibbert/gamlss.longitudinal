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
    trace) {
  start_time <- Sys.time()
  warnings <- character(0)
  fit <- NULL
  error <- NA_character_

  if (is.null(margin_dist)) {
    error <- paste0("Could not construct gamlss.dist family object for '", margin_family, "'.")
  } else {
    formulas <- .joint_selection_fit_formulas(
      response_var = response_var,
      time_var = time_var,
      time_intercepts = time_intercepts,
      copula_time_intercepts = copula_time_intercepts
    )
    default_args <- list(
      dataset = data,
      margin_dist = margin_dist,
      copula_dist = copula_family,
      time_var = time_var,
      subject_var = subject_var,
      mu.formula = formulas$mu_formula,
      sigma.formula = formulas$par_formula,
      nu.formula = formulas$par_formula,
      tau.formula = formulas$par_formula,
      theta.formula = formulas$theta_formula,
      zeta.formula = ~1,
      include_dlcopdpar = TRUE,
      compute_vcov = FALSE,
      optimizer_control = gamlss_longitudinal_control(
        rs = list(warm_start_joint = FALSE)
      ),
      verbose = 0
    )
    legacy_optimizer_names <- c(
      names(.gl_legacy_optimizer_map()),
      "use_backtracking", "backtracking_max_halves"
    )
    if (any(names(fit_args) %in% legacy_optimizer_names)) {
      default_args$optimizer_control <- NULL
      default_args$warm_start_joint <- FALSE
    }
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
      row <- .joint_selection_success_row(
        margin_family = margin_family,
        copula_family = copula_family,
        fit = fit,
        fit_metrics = fit_metrics,
        elapsed = elapsed,
        warnings = warnings
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
