.gl_bootstrap_refit <- function(object, response, fit_args) {
  if (is.null(object$dataset) || is.null(object$formulas_int)) {
    stop(

      "This fitted object does not store the internal dataset/formulas needed for refitting. Refit the model with the current package version.",
      call. = FALSE
    )
  }

  dat <- object$dataset

  if (!"response" %in% names(dat)) {
    stop("Stored model dataset does not contain a 'response' column.", call. = FALSE)
  }

  if (length(response) != nrow(dat)) {
    stop("Bootstrap response length does not match the stored model dataset.", call. = FALSE)
  }

  dat$response <- response

  args <- c(
    list(
      dataset = dat,
      margin_dist = object$margin_dist,
      copula_dist = object$copula_dist,
      time_var = "time",
      subject_var = "subject",
      mu.formula = object$formulas_int$mu,
      sigma.formula = object$formulas_int$sigma,
      nu.formula = object$formulas_int$nu,
      tau.formula = object$formulas_int$tau,
      theta.formula = object$formulas_int$theta,
      zeta.formula = object$formulas_int$zeta,
      start_from = object$par,
      include_dlcopdpar = object$include_dlcopdpar,
      compute_vcov = FALSE,
      verbose = 0
    ),
    fit_args
  )

  do.call(gamlss_longitudinal, args)
}

.gl_bootstrap_refit_dataset <- function(object, dat, fit_args) {
  if (is.null(object$formulas_int)) {
    stop(
      "This fitted object does not store the internal formulas needed for refitting. Refit the model with the current package version.",
      call. = FALSE
    )
  }

  args <- c(
    list(
      dataset = dat,
      margin_dist = object$margin_dist,
      copula_dist = object$copula_dist,
      time_var = "time",
      subject_var = "subject",
      mu.formula = object$formulas_int$mu,
      sigma.formula = object$formulas_int$sigma,
      nu.formula = object$formulas_int$nu,
      tau.formula = object$formulas_int$tau,
      theta.formula = object$formulas_int$theta,
      zeta.formula = object$formulas_int$zeta,
      start_from = object$par,
      include_dlcopdpar = object$include_dlcopdpar,
      compute_vcov = FALSE,
      verbose = 0
    ),
    fit_args
  )

  do.call(gamlss_longitudinal, args)
}

#' Bootstrap inference for fitted models

#'

#' `bootstrap_inference()` either simulates responses from a fitted

#' `gamlss.longitudinal` model, or resamples subjects with replacement,
#' refits the same model to each bootstrap dataset,

#' and summarizes the bootstrap distribution of selected fixed coefficients.

#' It is intended for opt-in applied uncertainty checks and should be run with

#' enough replicates outside CRAN-time tests for final reporting.

#'

#' @param object A fitted `gamlss.longitudinal` object.

#' @param R Number of bootstrap replicates.

#' @param terms Optional coefficient names, formula-term names such as

#'   `"mu.treatment"`, coefficient-name prefixes, or numeric indices to

#'   summarize.

#' @param level Confidence level for percentile intervals.

#' @param seed Optional random seed.

#' @param fit_args Optional named list of arguments passed to each refit, such

#'   as `max_outer_iter`, `max_inner_iter`, or convergence tolerances.

#' @param keep_fits Logical; keep successful refitted model objects.
#' @param type Character; `"parametric"` simulates from the fitted copula model,
#'   while `"cluster"` resamples subjects with replacement.

#' @param ... Additional arguments passed to [simulate.gamlss.longitudinal()].

#'

#' @return An object of class `gamlss_longitudinal_bootstrap`.

#' @export

bootstrap_inference <- function(
    object,
    R = 100,
    terms = NULL,
    level = 0.95,
    seed = NULL,
    fit_args = list(),
    keep_fits = FALSE,
    type = c("parametric", "cluster"),
    ...) {
  validated <- .gl_validate_bootstrap_args(object, R, fit_args)

  R <- validated$R

  fit_args <- validated$fit_args

  type <- match.arg(type)

  dots <- list(...)

  if (identical(type, "parametric")) {
    dots <- .gl_normalize_bootstrap_simulation_args(dots)
  } else if (length(dots) > 0L) {
    stop("'...' simulation arguments are only supported for type = 'parametric'.", call. = FALSE)
  }

  estimates <- stats::coef(object)

  idx <- .resolve_coefficient_terms(

    terms,
    names(estimates),
    arg = "terms",
    term_map = .fixed_term_coefficient_names(object)
  )

  terms_use <- names(estimates)[idx]

  seed_restore <- .gl_prepare_bootstrap_seed(seed)
  if (!is.null(seed_restore)) {
    on.exit(seed_restore(), add = TRUE)
  }

  sim <- if (identical(type, "parametric")) {
    do.call(simulate, c(list(object = object, nsim = R), dots))
  } else {
    NULL
  }

  boot_coef <- matrix(NA_real_, nrow = R, ncol = length(terms_use))

  colnames(boot_coef) <- terms_use

  errors <- rep(NA_character_, R)

  fits <- if (isTRUE(keep_fits)) vector("list", R) else NULL

  missing_mask <- is.na(object$response)

  for (b in seq_len(R)) {
    fit_b <- if (identical(type, "parametric")) {
      response_b <- sim[[b]]
      response_b[missing_mask] <- NA_real_

      tryCatch(
        suppressWarnings(.gl_bootstrap_refit(object, response_b, fit_args = fit_args)),
        error = function(e) e
      )
    } else {
      dat_b <- .gl_cluster_bootstrap_dataset(object)

      tryCatch(
        suppressWarnings(.gl_bootstrap_refit_dataset(object, dat_b, fit_args = fit_args)),
        error = function(e) e
      )
    }

    if (inherits(fit_b, "error")) {
      errors[[b]] <- conditionMessage(fit_b)

      next
    }

    coef_b <- stats::coef(fit_b)

    boot_coef[b, ] <- coef_b[terms_use]

    if (isTRUE(keep_fits)) {
      fits[[b]] <- fit_b
    }
  }

  summary <- .gl_bootstrap_summary(estimates, terms_use, boot_coef, level)

  .gl_bootstrap_result(
    summary = summary,
    boot_coef = boot_coef,
    errors = errors,
    R = R,
    level = level,
    fits = fits,
    simulation_type = type
  )
}

#' @export

print.gamlss_longitudinal_bootstrap <- function(x, digits = max(3, getOption("digits") - 3), ...) {
  title <- if (identical(x$simulation_type, "cluster")) {
    "Cluster Bootstrap for gamlss.longitudinal"
  } else {
    "Parametric Bootstrap for gamlss.longitudinal"
  }

  cat("\n", title, "\n", sep = "")

  cat(paste(rep("-", nchar(title)), collapse = ""), "\n", sep = "")

  cat("Replicates:", x$R, "\n")

  simulation_label <- if (identical(x$simulation_type, "cluster")) {
    "subject resampling"
  } else {
    "fitted copula model"
  }

  cat("Simulation:", simulation_label, "\n")

  cat("Failed refits:", x$failed_replicates, "\n\n")

  print(x$summary, digits = digits, row.names = FALSE)

  invisible(x)
}
