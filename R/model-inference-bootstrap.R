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
      data = dat,
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
      start = object$par,
      optimization = if (isTRUE(object$include_dlcopdpar)) "joint" else "separate",
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
      data = dat,
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
      start = object$par,
      optimization = if (isTRUE(object$include_dlcopdpar)) "joint" else "separate",
      compute_vcov = FALSE,
      verbose = 0
    ),
    fit_args
  )

  do.call(gamlss_longitudinal, args)
}

.gl_bootstrap_capture_attempt <- function(expr) {
  records <- list()
  record <- function(condition, type) {
    records[[length(records) + 1L]] <<- list(
      sequence = length(records) + 1L, type = type,
      condition = condition, classes = class(condition),
      message = conditionMessage(condition), call = conditionCall(condition)
    )
  }
  # Calling handlers observe conditions before the existing warning suppression.
  # Messages still propagate; duplicate warnings remain distinct records.
  value <- tryCatch(suppressWarnings(withCallingHandlers(force(expr),
    warning = function(w) record(w, "warning"),
    message = function(m) record(m, "message"),
    error = function(e) record(e, "error"))), error = function(e) e)
  list(value = value, conditions = records)
}

.gl_bootstrap_call_attempt_callback <- function(callback, record) {
  # Storage callbacks must not perturb later cluster draws or randomized refits.
  kind <- RNGkind()
  present <- exists(".Random.seed", envir = .GlobalEnv, inherits = FALSE)
  state <- if (present) get(".Random.seed", envir = .GlobalEnv, inherits = FALSE) else NULL
  on.exit({
    do.call(RNGkind, as.list(kind))
    if (present) assign(".Random.seed", state, envir = .GlobalEnv)
    else if (exists(".Random.seed", envir = .GlobalEnv, inherits = FALSE))
      rm(".Random.seed", envir = .GlobalEnv)
  }, add = TRUE)
  callback(record) # Deliberately outside refit error handling: storage errors abort.
  invisible(NULL)
}

#' Bootstrap inference for fitted models

#'

#' `bootstrap_inference()` either simulates responses from a fitted

#' `gamlss.longitudinal` model, or resamples subjects with replacement,
#' refits the same model to each bootstrap dataset,

#' and summarizes the bootstrap distribution of selected fixed coefficients.
#' Errors, explicit nonconvergence, and missing/non-finite target coefficients
#' are failed replicates; per-coefficient successful counts are reported.

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
#' @param attempt_callback Optional function of one record, called synchronously
#'   after every native refit return or error, before convergence/coefficient
#'   filtering. The record contains `attempt`, `R`, `type`, exact refit `data`,
#'   `seed`, `rng_kind`, `rng_before_refit`, `rng_after_refit`, `fit_args`,
#'   `coefficient_names`, complete `fit` (including failed eligibility), original
#'   `error` or `NULL`, and ordered `conditions` (original condition objects,
#'   types, classes, messages and calls). Conditions cover the refit only, not
#'   initial simulation or subsequent filtering. Callback return values are
#'   ignored; a callback error aborts before any later attempt. The default
#'   `NULL` retains existing behavior. This does not change `keep_fits`, retry
#'   failed fits, or certify inference. Callback RNG changes are restored; other
#'   external side effects are the callback author's responsibility. Use a
#'   callback that writes and verifies each record for durable attempt retention.

#'

#' @return An object of class `gamlss_longitudinal_bootstrap` with an exact
#'   fixed-coefficient target, replicate failures, successful counts, and
#'   adequacy status in `inference_contract`. Two successful replicates are a
#'   computational minimum only, not evidence of inferential adequacy.

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
    ..., attempt_callback = NULL) {
  .gl_require_converged_fit(object, "bootstrap inference")
  if (!is.null(attempt_callback) && !is.function(attempt_callback))
    stop("'attempt_callback' must be NULL or a function.", call. = FALSE)
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
    rng_before_refit <- NULL
    if (!is.null(attempt_callback)) {
      if (identical(type, "parametric")) {
        response_b <- sim[[b]]
        response_b[missing_mask] <- NA_real_
        dat_b <- object$dataset
        dat_b$response <- response_b
      } else {
        dat_b <- .gl_cluster_bootstrap_dataset(object)
      }
      rng_kind <- RNGkind()
      rng_before_refit <- if (exists(".Random.seed", .GlobalEnv, inherits = FALSE))
        get(".Random.seed", .GlobalEnv, inherits = FALSE) else NULL
      captured <- if (identical(type, "parametric")) {
        .gl_bootstrap_capture_attempt(.gl_bootstrap_refit(object, response_b, fit_args = fit_args))
      } else {
        .gl_bootstrap_capture_attempt(.gl_bootstrap_refit_dataset(object, dat_b, fit_args = fit_args))
      }
      fit_b <- captured$value
      .gl_bootstrap_call_attempt_callback(attempt_callback, list(
        attempt = b, R = R, type = type, data = dat_b, seed = seed,
        rng_kind = rng_kind, rng_before_refit = rng_before_refit,
        rng_after_refit = if (exists(".Random.seed", .GlobalEnv, inherits = FALSE))
          get(".Random.seed", .GlobalEnv, inherits = FALSE) else NULL,
        fit_args = fit_args, coefficient_names = terms_use,
        fit = if (inherits(fit_b, "error")) NULL else fit_b,
        error = if (inherits(fit_b, "error")) fit_b else NULL,
        conditions = captured$conditions
      ))
    } else {
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
    }

    if (inherits(fit_b, "error")) {
      errors[[b]] <- conditionMessage(fit_b)

      next
    }

    if (identical(fit_b$convergence$converged, FALSE)) {
      errors[[b]] <- "Bootstrap refit did not satisfy the model convergence contract."
      next
    }

    coef_b <- stats::coef(fit_b)

    if (!all(terms_use %in% names(coef_b)) || any(!is.finite(coef_b[terms_use]))) {
      errors[[b]] <- "Bootstrap refit returned missing or non-finite target fixed coefficients."
      next
    }

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
    "Cluster bootstrap"
  } else {
    "Parametric bootstrap"
  }

  cat("\n<", title, ">\n", sep = "")

  cat("Replicates:", x$R, "\n")

  simulation_label <- if (identical(x$simulation_type, "cluster")) {
    "subject resampling"
  } else {
    "fitted copula model"
  }

  cat("Simulation:", simulation_label, "\n")

  cat("Failed refits:", x$failed_replicates, "\n\n")

  cat("Target: selected fixed coefficients:", paste(x$inference_contract$coefficient_names, collapse = ", "), "\n")
  cat("Validity:", x$inference_contract$validity_status, "\n")
  cat("Adequacy:", x$inference_contract$inferential_adequacy, "\n")
  cat("Model-selection uncertainty is not included.\n\n")

  print(x$summary, digits = digits, row.names = FALSE)

  invisible(x)
}
