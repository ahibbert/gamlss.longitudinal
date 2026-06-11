.fixed_term_coefficient_names <- function(object) {

  if (!inherits(object, "gamlss.longitudinal") || is.null(object$model_matrix$x)) {

    return(NULL)

  }


  out <- list()

  for (parameter in names(object$model_matrix$x)) {

    X <- object$model_matrix$x[[parameter]]

    if (is.null(X)) next

    col_names <- colnames(X)

    if (length(col_names) == 0L) next


    assign <- attr(X, "assign")

    term_labels <- attr(X, "term.labels")

    if (length(assign) != length(col_names) || length(term_labels) == 0L) next


    for (term_idx in seq_along(term_labels)) {

      cols <- col_names[assign == term_idx]

      if (length(cols) == 0L) next

      out[[paste(parameter, term_labels[[term_idx]], sep = ".")]] <- paste(parameter, cols, sep = ".")

    }

  }

  out

}


.resolve_coefficient_terms <- function(terms, coefficient_names, arg = "terms", term_map = NULL) {

  if (is.null(terms)) {

    return(seq_along(coefficient_names))

  }


  if (!is.character(terms)) {

    idx <- terms

    if (length(idx) == 0L || any(is.na(idx)) || any(idx < 1L) || any(idx > length(coefficient_names))) {

      stop(sprintf("'%s' contains unknown coefficient names or indices.", arg), call. = FALSE)

    }

    return(as.integer(idx))

  }


  if (length(terms) == 0L || any(is.na(terms)) || any(!nzchar(terms))) {

    stop(sprintf("'%s' contains unknown coefficient names or indices.", arg), call. = FALSE)

  }


  idx_list <- lapply(terms, function(term) {

    exact <- match(term, coefficient_names)

    if (!is.na(exact)) {

      return(exact)

    }

    if (!is.null(term_map) && term %in% names(term_map)) {

      mapped_idx <- match(term_map[[term]], coefficient_names)

      mapped_idx <- mapped_idx[!is.na(mapped_idx)]

      if (length(mapped_idx) > 0L) {

        return(mapped_idx)

      }

    }

    prefix_idx <- which(startsWith(coefficient_names, term))

    if (length(prefix_idx) == 0L) {

      return(prefix_idx)

    }

    suffix <- substring(coefficient_names[prefix_idx], nchar(term) + 1L)

    main_idx <- prefix_idx[!grepl(":", suffix, fixed = TRUE)]

    if (length(main_idx) > 0L) {

      return(main_idx)

    }

    prefix_idx

  })

  missing <- vapply(idx_list, length, integer(1)) == 0L

  if (any(missing)) {

    stop(sprintf("'%s' contains unknown coefficient names or indices.", arg), call. = FALSE)

  }

  unique(as.integer(unlist(idx_list, use.names = FALSE)))

}


#' Wald tests for fixed coefficients

#'

#' `wald_test()` provides a small reporting-friendly hypothesis-test surface for

#' fitted `gamlss.longitudinal` models. It uses the same variance-covariance

#' route as [summary.gamlss.longitudinal()] and [confint.gamlss.longitudinal()],

#' so numerical-Hessian tests should be reported as approximate.

#'

#' @param object A fitted `gamlss.longitudinal` object.

#' @param terms Optional coefficient names, formula-term names such as

#'   `"mu.treatment"`, coefficient-name prefixes, or numeric indices. When `L`

#'   is `NULL`, these select coefficients for individual tests or a joint test.

#' @param L Optional contrast matrix. Columns must either be named with

#'   coefficient names or have one column per fixed coefficient in model order.

#' @param rhs Null-hypothesis value. Either a scalar or one value per tested row.

#' @param joint Logical; when `TRUE`, test selected `terms` jointly. Contrast

#'   matrices supplied through `L` are always tested jointly.

#' @param method Variance-covariance method passed to [vcov.gamlss.longitudinal()].

#' @param ... Additional arguments passed to [vcov.gamlss.longitudinal()].

#'

#' @return An object of class `gamlss_longitudinal_wald_test`.

#' @export

wald_test <- function(

  object,

  terms = NULL,

  L = NULL,

  rhs = 0,

  joint = FALSE,

  method = "analytical",

  ...

) {

  if (!inherits(object, "gamlss.longitudinal")) {

    stop("'object' must be a fitted 'gamlss.longitudinal' object.", call. = FALSE)

  }

  estimates <- stats::coef(object)

  if (length(estimates) == 0L) {

    stop("No fixed coefficients are available to test.", call. = FALSE)

  }

  vc <- .resolve_vcov(object, extra_args = list(method = method, ...))

  V <- vc$vcov$overall

  if (is.null(V) || !is.matrix(V)) {

    stop("A fixed-effect variance-covariance matrix is required for Wald tests.", call. = FALSE)

  }

  V_names <- colnames(V) %||% rownames(V)

  if (is.null(V_names)) {

    V_names <- names(estimates)

  }

  idx_v <- match(names(estimates), V_names)

  if (any(is.na(idx_v))) {

    stop("Variance-covariance matrix names do not match fitted coefficients.", call. = FALSE)

  }

  V <- V[idx_v, idx_v, drop = FALSE]

  rownames(V) <- colnames(V) <- names(estimates)


  if (!is.null(L)) {

    L <- as.matrix(L)

    if (nrow(L) == 0L || ncol(L) == 0L) {

      stop("'L' must contain at least one contrast row and one coefficient column.", call. = FALSE)

    }

    if (!is.null(colnames(L))) {

      idx <- match(colnames(L), names(estimates))

      if (any(is.na(idx))) {

        stop("Column names in 'L' must match coefficient names.", call. = FALSE)

      }

      L_full <- matrix(0, nrow = nrow(L), ncol = length(estimates))

      colnames(L_full) <- names(estimates)

      rownames(L_full) <- rownames(L)

      L_full[, idx] <- L

      L <- L_full

    } else if (ncol(L) != length(estimates)) {

      stop("Unnamed 'L' must have one column per fixed coefficient.", call. = FALSE)

    } else {

      colnames(L) <- names(estimates)

    }

    if (is.null(rownames(L))) {

      rownames(L) <- paste0("H", seq_len(nrow(L)))

    }

    joint <- TRUE

  } else {

    idx <- .resolve_coefficient_terms(

      terms,

      names(estimates),

      arg = "terms",

      term_map = .fixed_term_coefficient_names(object)

    )

    L <- diag(length(estimates))[idx, , drop = FALSE]

    colnames(L) <- names(estimates)

    rownames(L) <- names(estimates)[idx]

  }


  rhs <- rep(rhs, length.out = nrow(L))

  estimate <- as.numeric(L %*% estimates)

  diff <- estimate - rhs

  LVL <- L %*% V %*% t(L)


  if (isTRUE(joint)) {

    stat <- tryCatch(

      as.numeric(t(diff) %*% solve(LVL, diff)),

      error = function(e) NA_real_

    )

    out <- data.frame(

      hypothesis = paste(rownames(L), collapse = ", "),

      df = nrow(L),

      statistic = stat,

      p_value = stats::pchisq(stat, df = nrow(L), lower.tail = FALSE),

      method = vc$method %||% method,

      stringsAsFactors = FALSE

    )

  } else {

    se <- sqrt(pmax(0, diag(LVL)))

    z <- diff / se

    out <- data.frame(

      term = rownames(L),

      estimate = estimate,

      rhs = rhs,

      std_error = se,

      statistic = z,

      p_value = 2 * stats::pnorm(abs(z), lower.tail = FALSE),

      method = vc$method %||% method,

      stringsAsFactors = FALSE

    )

  }

  attr(out, "joint") <- isTRUE(joint)

  attr(out, "method_requested") <- method

  attr(out, "vcov_method") <- vc$method %||% method

  class(out) <- c("gamlss_longitudinal_wald_test", "data.frame")

  out

}


#' @export

print.gamlss_longitudinal_wald_test <- function(x, digits = max(3, getOption("digits") - 3), ...) {

  cat("\nWald Test for gamlss.longitudinal\n")

  cat("---------------------------------\n")

  cat("Test type:", if (isTRUE(attr(x, "joint"))) "joint" else "individual", "\n")

  cat("VCOV method:", attr(x, "vcov_method") %||% "unknown", "\n\n")

  print.data.frame(x, digits = digits, row.names = FALSE, ...)

  invisible(x)

}


.gl_model_edf <- function(object) {

  fixed <- length(object$par)

  smooth <- 0

  if (!is.null(object$df_s) && length(object$df_s) > 0L) {

    vals <- suppressWarnings(as.numeric(unlist(object$df_s, use.names = FALSE)))

    vals <- vals[is.finite(vals)]

    smooth <- sum(vals)

  }

  fixed + smooth

}


.gl_joint_loglik <- function(object) {

  ll <- NULL

  if (!is.null(object$calc_lik_out_end) && !is.null(object$calc_lik_out_end$log_lik)) {

    ll <- object$calc_lik_out_end$log_lik

  } else if (!is.null(object$calc_lik_out) && !is.null(object$calc_lik_out$log_lik)) {

    ll <- object$calc_lik_out$log_lik

  }

  if (!is.null(ll) && "joint" %in% names(ll)) {

    return(as.numeric(ll[["joint"]]))

  }

  NA_real_

}


#' Compare fitted models with likelihood-ratio summaries

#'

#' `likelihood_compare()` gives a compact sequential likelihood comparison for

#' nested or approximately nested `gamlss.longitudinal` models. It reports joint

#' log-likelihood, effective degrees of freedom, AIC, BIC, likelihood-ratio

#' increments, and chi-square reference p-values.

#'

#' @param ... Fitted `gamlss.longitudinal` objects, or a single list of fitted

#'   objects.

#' @param sort Logical; order models by effective degrees of freedom before

#'   computing sequential comparisons.

#'

#' @return An object of class `gamlss_longitudinal_likelihood_compare`.

#' @export

likelihood_compare <- function(..., sort = TRUE) {

  models <- list(...)

  if (length(models) == 1L && is.list(models[[1L]]) && !inherits(models[[1L]], "gamlss.longitudinal")) {

    models <- models[[1L]]

  }

  if (length(models) < 2L) {

    stop("At least two fitted models are required.", call. = FALSE)

  }

  ok <- vapply(models, inherits, logical(1), what = "gamlss.longitudinal")

  if (!all(ok)) {

    stop("All inputs must be fitted 'gamlss.longitudinal' objects.", call. = FALSE)

  }

  labels <- names(models)

  if (is.null(labels) || any(labels == "")) {

    labels <- paste0("model_", seq_along(models))

  }

  n_obs <- vapply(models, function(x) length(x$response), integer(1))

  if (length(unique(n_obs)) > 1L) {

    warning(

      "Models have different observation counts; likelihood-ratio comparisons may not be valid.",

      call. = FALSE

    )

  }

  df <- vapply(models, .gl_model_edf, numeric(1))

  loglik <- vapply(models, .gl_joint_loglik, numeric(1))

  if (isTRUE(sort)) {

    ord <- order(df, loglik)

    labels <- labels[ord]

    n_obs <- n_obs[ord]

    df <- df[ord]

    loglik <- loglik[ord]

  }

  aic <- -2 * loglik + 2 * df

  bic <- -2 * loglik + log(pmax(1, n_obs)) * df

  delta_df <- c(NA_real_, diff(df))

  lr <- c(NA_real_, 2 * diff(loglik))

  p_value <- rep(NA_real_, length(models))

  valid <- is.finite(lr) & is.finite(delta_df) & delta_df > 0

  p_value[valid] <- stats::pchisq(lr[valid], df = delta_df[valid], lower.tail = FALSE)


  out <- data.frame(

    model = labels,

    n_obs = as.integer(n_obs),

    df = as.numeric(df),

    logLik = as.numeric(loglik),

    AIC = as.numeric(aic),

    BIC = as.numeric(bic),

    delta_df = as.numeric(delta_df),

    LR_statistic = as.numeric(lr),

    p_value = as.numeric(p_value),

    stringsAsFactors = FALSE

  )

  class(out) <- c("gamlss_longitudinal_likelihood_compare", "data.frame")

  out

}


#' @export

print.gamlss_longitudinal_likelihood_compare <- function(x, digits = max(3, getOption("digits") - 3), ...) {

  cat("\nLikelihood Comparison for gamlss.longitudinal\n")

  cat("--------------------------------------------\n")

  cat("Sequential LR rows compare each model with the previous row.\n\n")

  print.data.frame(x, digits = digits, row.names = FALSE, ...)

  invisible(x)

}


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


#' Parametric bootstrap inference for fitted models

#'

#' `bootstrap_inference()` simulates responses from a fitted

#' `gamlss.longitudinal` model, refits the same model to each simulated response,

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

  ...

) {

  if (!inherits(object, "gamlss.longitudinal")) {

    stop("'object' must be a fitted 'gamlss.longitudinal' object.", call. = FALSE)

  }

  R <- as.integer(R)

  if (length(R) != 1L || is.na(R) || R < 1L) {

    stop("'R' must be a positive integer.", call. = FALSE)

  }

  if (!is.list(fit_args) || (length(fit_args) > 0L && (is.null(names(fit_args)) || any(names(fit_args) == "")))) {

    stop("'fit_args' must be a named list.", call. = FALSE)

  }

  dots <- list(...)

  if ("simulation_type" %in% names(dots)) {

    simulation_type <- dots$simulation_type

    if (!identical(simulation_type, "copula")) {

      stop("'simulation_type' is no longer supported; bootstrap_inference() simulates from the fitted copula model.", call. = FALSE)

    }

    dots$simulation_type <- NULL

  }

  estimates <- stats::coef(object)

  idx <- .resolve_coefficient_terms(

    terms,

    names(estimates),

    arg = "terms",

    term_map = .fixed_term_coefficient_names(object)

  )

  terms_use <- names(estimates)[idx]

  if (!is.null(seed)) {

    old_seed <- if (exists(".Random.seed", envir = .GlobalEnv, inherits = FALSE)) {

      get(".Random.seed", envir = .GlobalEnv, inherits = FALSE)

    } else {

      NULL

    }

    on.exit({

      if (is.null(old_seed)) {

        if (exists(".Random.seed", envir = .GlobalEnv, inherits = FALSE)) {

          rm(".Random.seed", envir = .GlobalEnv)

        }

      } else {

        assign(".Random.seed", old_seed, envir = .GlobalEnv)

      }

    }, add = TRUE)

    set.seed(seed)

  }


  sim <- do.call(simulate, c(list(object = object, nsim = R), dots))

  boot_coef <- matrix(NA_real_, nrow = R, ncol = length(terms_use))

  colnames(boot_coef) <- terms_use

  errors <- rep(NA_character_, R)

  fits <- if (isTRUE(keep_fits)) vector("list", R) else NULL

  missing_mask <- is.na(object$response)


  for (b in seq_len(R)) {

    response_b <- sim[[b]]

    response_b[missing_mask] <- NA_real_

    fit_b <- tryCatch(

      suppressWarnings(.gl_bootstrap_refit(object, response_b, fit_args = fit_args)),

      error = function(e) e

    )

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


  alpha <- 1 - level

  q_na <- function(x, prob) {

    x <- x[is.finite(x)]

    if (length(x) == 0L) return(NA_real_)

    as.numeric(stats::quantile(x, probs = prob, names = FALSE))

  }

  summary <- data.frame(

    term = terms_use,

    estimate = as.numeric(estimates[terms_use]),

    bootstrap_mean = colMeans(boot_coef, na.rm = TRUE),

    bootstrap_se = apply(boot_coef, 2, stats::sd, na.rm = TRUE),

    conf.low = apply(boot_coef, 2, q_na, prob = alpha / 2),

    conf.high = apply(boot_coef, 2, q_na, prob = 1 - alpha / 2),

    reps = colSums(is.finite(boot_coef)),

    stringsAsFactors = FALSE

  )


  out <- list(

    summary = summary,

    replicates = as.data.frame(boot_coef, stringsAsFactors = FALSE),

    errors = errors,

    R = R,

    successful_replicates = sum(stats::complete.cases(boot_coef)),

    failed_replicates = sum(!is.na(errors)),

    level = level,

    simulation_type = "copula",

    fits = fits

  )

  class(out) <- "gamlss_longitudinal_bootstrap"

  out

}


#' @export

print.gamlss_longitudinal_bootstrap <- function(x, digits = max(3, getOption("digits") - 3), ...) {

  cat("\nParametric Bootstrap for gamlss.longitudinal\n")

  cat("-------------------------------------------\n")

  cat("Replicates:", x$R, "\n")

  cat("Simulation:", "fitted copula model", "\n")

  cat("Failed refits:", x$failed_replicates, "\n\n")

  print(x$summary, digits = digits, row.names = FALSE)

  invisible(x)

}


