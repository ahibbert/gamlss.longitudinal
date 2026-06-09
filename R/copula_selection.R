#' Select a bivariate copula family from pseudo-observation pairs
#'
#' Screens supported copula families by maximising their native copula
#' log-likelihood on adjacent-time pseudo-observation pairs. Pseudo-observations
#' can be supplied directly as `u1`/`u2`, supplied as a row-aligned uniform
#' vector or column in `data`, or computed from a fitted
#' `gamlss.longitudinal` object.
#'
#' This helper is intended as a lightweight family-screening step. By default
#' it estimates a constant dependence parameter for each candidate family. With
#' `copula_time_intercepts = TRUE`, it estimates separate dependence parameters
#' for each observed adjacent time pair, treating time as a factor rather than
#' a linear trend. Richer covariate or smooth dependence structures should be
#' fitted afterwards with [gamlss_longitudinal()].
#'
#' @param data Optional long-format data frame.
#' @param object Optional fitted `gamlss.longitudinal` object.
#' @param u1,u2 Optional vectors of paired pseudo-observations.
#' @param u Optional row-aligned vector of pseudo-observations for `data`.
#' @param u_var Optional name of a pseudo-observation column in `data`.
#' @param response_var Optional response column name in `data`. Used to create
#'   pseudo-observations automatically when `u1`/`u2`, `u`, `u_var`, and
#'   `object` are not supplied.
#' @param margin_dist Optional `gamlss.dist` family object, or a
#'   [select_margin()] result, used to fit the temporary marginal model for
#'   automatic pseudo-observations. If omitted, [select_margin()] is called and
#'   a warning is issued.
#' @param mu.formula,sigma.formula,nu.formula,tau.formula Optional temporary
#'   marginal model formulas used when `select_copula()` creates
#'   pseudo-observations from `response_var`. If omitted, intercept-only
#'   formulas are used unless `time_intercepts = TRUE` or `margin_dist` is a
#'   time-intercept [select_margin()] result.
#' @param subject_var,time_var Subject and time column names used when building
#'   adjacent-time pairs from `data`.
#' @param families Candidate copula family codes. Supported values are `"N"`,
#'   `"C"`, `"F"`, `"G"`, `"J"`, and `"t"`.
#' @param lags Positive integer lag(s) used when forming adjacent pairs.
#' @param criterion Ranking criterion, one of `"AIC"`, `"BIC"`, or `"logLik"`.
#' @param t_df_grid Degrees-of-freedom grid used for the t-copula screen.
#' @param min_pairs Minimum number of complete pairs required.
#' @param time_intercepts Logical; when `select_copula()` creates
#'   pseudo-observations from `response_var`, use time-specific intercepts in
#'   the temporary marginal model and pass this setting through to
#'   [select_margin()] if the marginal family is auto-selected.
#' @param copula_time_intercepts Logical; if `TRUE`, screen each copula family
#'   with separate dependence parameters for each adjacent time-pair factor
#'   level. This does not impose a linear time trend.
#'
#' @return A data frame with one row per family and class
#'   `copula_selection`. The selected family is stored in the `selected`
#'   attribute.
#' @export
select_copula <- function(
  data = NULL,
  object = NULL,
  u1 = NULL,
  u2 = NULL,
  u = NULL,
  u_var = NULL,
  response_var = NULL,
  margin_dist = NULL,
  mu.formula = NULL,
  sigma.formula = NULL,
  nu.formula = NULL,
  tau.formula = NULL,
  subject_var = "subject",
  time_var = "time",
  families = c("N", "C", "F", "G", "J", "t"),
  lags = 1,
  criterion = c("AIC", "BIC", "logLik"),
  t_df_grid = c(3, 4, 6, 8, 12, 20, 30),
  min_pairs = 10,
  time_intercepts = FALSE,
  copula_time_intercepts = FALSE
) {
  criterion <- match.arg(criterion)
  families <- vapply(families, .copula_family_code, character(1), USE.NAMES = FALSE)
  lags <- as.integer(lags)
  if (length(lags) < 1L || any(!is.finite(lags)) || any(lags < 1L)) {
    stop("lags must contain positive integers.", call. = FALSE)
  }

  pairs <- .select_copula_pairs(
    data = data,
    object = object,
    u1 = u1,
    u2 = u2,
    u = u,
    u_var = u_var,
    response_var = response_var,
    margin_dist = margin_dist,
    mu.formula = mu.formula,
    sigma.formula = sigma.formula,
    nu.formula = nu.formula,
    tau.formula = tau.formula,
    subject_var = subject_var,
    time_var = time_var,
    lags = lags,
    time_intercepts = time_intercepts
  )
  keep <- is.finite(pairs$u1) & is.finite(pairs$u2)
  pairs <- pairs[keep, , drop = FALSE]
  if (nrow(pairs) < min_pairs) {
    stop("At least ", min_pairs, " complete pseudo-observation pairs are required.", call. = FALSE)
  }
  if (isTRUE(copula_time_intercepts) && !"copula_time" %in% names(pairs)) {
    stop(
      "'copula_time_intercepts = TRUE' requires data, u/u_var, response_var, or object input with time information.",
      call. = FALSE
    )
  }

  fits <- lapply(families, function(family) {
    .select_copula_fit_family(
      u1 = pairs$u1,
      u2 = pairs$u2,
      family = family,
      t_df_grid = t_df_grid,
      copula_time = if (isTRUE(copula_time_intercepts)) pairs$copula_time else NULL
    )
  })
  out <- do.call(rbind, fits)
  out$n_pairs <- nrow(pairs)
  out <- out[order(out[[criterion]], decreasing = identical(criterion, "logLik")), , drop = FALSE]
  rownames(out) <- NULL
  attr(out, "selected") <- out$family[1]
  attr(out, "criterion") <- criterion
  attr(out, "margin_selection") <- attr(pairs, "margin_selection")
  attr(out, "pseudo_observation_source") <- attr(pairs, "pseudo_observation_source")
  attr(out, "copula_time_intercepts") <- isTRUE(copula_time_intercepts)
  attr(out, "copula_time_levels") <- if (isTRUE(copula_time_intercepts)) unique(pairs$copula_time) else NULL
  class(out) <- c("copula_selection", "data.frame")
  out
}

#' @export
best_fit.copula_selection <- function(x, ...) {
  if (nrow(x) == 0L) {
    return(list(
      family = NA_character_,
      criterion = attr(x, "criterion")
    ))
  }
  row <- as.list(as.data.frame(x)[1L, , drop = FALSE])
  row$criterion <- attr(x, "criterion")
  row
}

#' @export
best_fit_family.copula_selection <- function(x, ...) {
  best_fit(x)$family
}

#' @export
`$.copula_selection` <- function(x, name) {
  if (identical(name, "best_fit")) {
    return(best_fit(x))
  }
  .subset2(as.data.frame(x), name, exact = FALSE)
}

.select_copula_pairs <- function(
  data,
  object,
  u1,
  u2,
  u,
  u_var,
  response_var,
  margin_dist,
  mu.formula,
  sigma.formula,
  nu.formula,
  tau.formula,
  subject_var,
  time_var,
  lags,
  time_intercepts = FALSE
) {
  if (!is.null(u1) || !is.null(u2)) {
    if (is.null(u1) || is.null(u2)) {
      stop("Both u1 and u2 must be supplied for direct pseudo-observation pairs.", call. = FALSE)
    }
    n <- max(length(u1), length(u2))
    return(data.frame(
      u1 = rep(.copula_clamp01(u1), length.out = n),
      u2 = rep(.copula_clamp01(u2), length.out = n)
    ))
  }

  if (!is.null(object)) {
    u <- .select_copula_u_from_fit(object)
    data <- data.frame(
      .subject = object$response_subject,
      .time = object$response_margin,
      .u = u
    )
    subject_var <- ".subject"
    time_var <- ".time"
    u_var <- ".u"
  } else {
    if (is.null(data)) {
      stop("Provide either object, u1/u2, or data with u/u_var.", call. = FALSE)
    }
    data <- as.data.frame(data)
    margin_selection_time_var <- .select_copula_margin_selection_time_var(margin_dist)
    if (!is.null(margin_selection_time_var) && !time_var %in% names(data) && margin_selection_time_var %in% names(data)) {
      time_var <- margin_selection_time_var
    }
    if (!is.null(u)) {
      if (length(u) != nrow(data)) {
        stop("u must have one value per row of data.", call. = FALSE)
      }
      data[[".u"]] <- u
      u_var <- ".u"
    }
    if (is.null(u_var)) {
      auto <- .select_copula_auto_u(
        data = data,
        response_var = response_var,
        margin_dist = margin_dist,
        mu.formula = mu.formula,
        sigma.formula = sigma.formula,
        nu.formula = nu.formula,
        tau.formula = tau.formula,
        time_var = time_var,
        time_intercepts = time_intercepts
      )
      data <- auto$data
      u_var <- auto$u_var
    }
  }

  if (!all(c(subject_var, time_var, u_var) %in% names(data))) {
    stop("data must contain subject_var, time_var, and u_var columns.", call. = FALSE)
  }
  pairs <- .select_copula_adjacent_pairs(data, subject_var = subject_var, time_var = time_var, u_var = u_var, lags = lags)
  if (exists("auto", inherits = FALSE)) {
    attr(pairs, "margin_selection") <- auto$margin_selection
    attr(pairs, "pseudo_observation_source") <- auto$source
  } else if (!is.null(object)) {
    attr(pairs, "pseudo_observation_source") <- "fitted_object"
  } else if (!is.null(u) || identical(u_var, ".u")) {
    attr(pairs, "pseudo_observation_source") <- "u"
  } else {
    attr(pairs, "pseudo_observation_source") <- u_var
  }
  pairs
}

.select_copula_auto_u <- function(
  data,
  response_var,
  margin_dist,
  mu.formula,
  sigma.formula,
  nu.formula,
  tau.formula,
  time_var,
  time_intercepts
) {
  if (is.null(response_var) || !is.character(response_var) || length(response_var) != 1L) {
    stop(
      "Provide u_var, u, u1/u2, object, or response_var so select_copula() can create pseudo-observations.",
      call. = FALSE
    )
  }
  if (!response_var %in% names(data)) {
    stop("response_var='", response_var, "' not found in 'data'.", call. = FALSE)
  }

  margin_selection <- NULL
  if (inherits(margin_dist, "margin_selection")) {
    margin_selection <- margin_dist
    if (isTRUE(attr(margin_selection, "time_intercepts"))) {
      time_intercepts <- TRUE
      selected_time_var <- attr(margin_selection, "time_var")
      if (!is.null(selected_time_var)) {
        time_var <- selected_time_var
      }
    }
    margin_dist <- best_fit_family(margin_selection)
    if (is.null(margin_dist)) {
      stop("The supplied select_margin() result did not contain a usable marginal family.", call. = FALSE)
    }
  } else if (is.null(margin_dist)) {
    warning(
      "'margin_dist' was not supplied; select_copula() is selecting a temporary marginal distribution with select_margin().",
      call. = FALSE
    )
    margin_selection <- suppressWarnings(suppressMessages(
      select_margin(
        data,
        response_var = response_var,
        time_var = if (isTRUE(time_intercepts)) time_var else NULL,
        time_intercepts = time_intercepts,
        trace = FALSE
      )
    ))
    margin_dist <- best_fit_family(margin_selection)
    if (is.null(margin_dist)) {
      stop("select_margin() did not return a usable marginal family.", call. = FALSE)
    }
  }

  if (isTRUE(time_intercepts)) {
    if (is.null(time_var) || !is.character(time_var) || length(time_var) != 1L) {
      stop("'time_var' must be a single column name when 'time_intercepts = TRUE'.", call. = FALSE)
    }
    if (!time_var %in% names(data)) {
      stop("time_var='", time_var, "' not found in 'data'.", call. = FALSE)
    }
    mu.formula <- .select_copula_time_response_formula(mu.formula, response_var, time_var)
    sigma.formula <- .select_copula_time_rhs_formula(sigma.formula, time_var)
    nu.formula <- .select_copula_time_rhs_formula(nu.formula, time_var)
    tau.formula <- .select_copula_time_rhs_formula(tau.formula, time_var)
  }

  pfun <- tryCatch(
    get(paste0("p", margin_dist$family[1]), envir = asNamespace("gamlss.dist"), mode = "function", inherits = FALSE),
    error = function(e) NULL
  )
  if (is.null(pfun)) {
    stop("'margin_dist' must provide a gamlss.dist CDF function.", call. = FALSE)
  }

  mu.formula <- .select_copula_response_formula(mu.formula, response_var)
  args <- list(
    formula = mu.formula,
    family = margin_dist,
    data = data,
    trace = FALSE
  )
  parameter_names <- names(margin_dist$parameters)
  if ("sigma" %in% parameter_names) args$sigma.fo <- .select_copula_rhs_formula(sigma.formula)
  if ("nu" %in% parameter_names) args$nu.fo <- .select_copula_rhs_formula(nu.formula)
  if ("tau" %in% parameter_names) args$tau.fo <- .select_copula_rhs_formula(tau.formula)

  formula_vars <- unique(unlist(
    lapply(
      c(list(args$formula), args[c("sigma.fo", "nu.fo", "tau.fo")]),
      function(formula) if (is.null(formula)) character(0) else all.vars(stats::as.formula(formula))
    ),
    use.names = FALSE
  ))
  missing_vars <- setdiff(formula_vars, names(data))
  if (length(missing_vars) > 0L) {
    stop(
      "Temporary marginal formula variable(s) not found in data: ",
      paste(missing_vars, collapse = ", "),
      call. = FALSE
    )
  }
  keep <- is.finite(data[[response_var]]) & stats::complete.cases(data[, formula_vars, drop = FALSE])
  if (sum(keep) < 3L) {
    stop("Need at least three complete response rows to create pseudo-observations.", call. = FALSE)
  }
  data_fit <- data[keep, , drop = FALSE]
  args$data <- data_fit[, formula_vars, drop = FALSE]

  margin_fit <- tryCatch(
    do.call(gamlss::gamlss, args),
    error = function(e) {
      stop("Temporary marginal model for copula selection failed: ", conditionMessage(e), call. = FALSE)
    }
  )

  cdf_args <- list(q = data_fit[[response_var]], y = data_fit[[response_var]], x = data_fit[[response_var]])
  for (parameter in parameter_names) {
    cdf_args[[parameter]] <- stats::fitted(margin_fit, what = parameter)
  }
  fixed_unlinked_values <- attr(margin_dist, "fixed_unlinked_values")
  if (length(fixed_unlinked_values) > 0L) {
    for (parameter in names(fixed_unlinked_values)) {
      value <- fixed_unlinked_values[[parameter]]
      if (is.numeric(value) && length(value) == 1L && is.finite(value)) {
        cdf_args[[parameter]] <- rep(value, nrow(data_fit))
      }
    }
  }
  cdf_args <- cdf_args[names(cdf_args) %in% formalArgs(pfun)]
  u <- tryCatch(
    do.call(pfun, cdf_args),
    error = function(e) {
      stop("Could not compute pseudo-observations from the temporary marginal model: ", conditionMessage(e), call. = FALSE)
    }
  )
  u_full <- rep(NA_real_, nrow(data))
  u_full[keep] <- .copula_clamp01(u)
  data[[".u_auto"]] <- u_full
  list(
    data = data,
    u_var = ".u_auto",
    margin_selection = margin_selection,
    source = if (is.null(margin_selection)) "margin_dist" else "select_margin"
  )
}

.select_copula_margin_selection_time_var <- function(margin_dist) {
  if (!inherits(margin_dist, "margin_selection") || !isTRUE(attr(margin_dist, "time_intercepts"))) {
    return(NULL)
  }
  time_var <- attr(margin_dist, "time_var")
  if (is.null(time_var) || !is.character(time_var) || length(time_var) != 1L || !nzchar(time_var)) {
    return(NULL)
  }
  time_var
}

.select_copula_time_response_formula <- function(formula, response_var, time_var) {
  if (!is.null(formula)) {
    return(formula)
  }
  stats::as.formula(paste(.select_copula_formula_name(response_var), "~ factor(", .select_copula_formula_name(time_var), ")"))
}

.select_copula_time_rhs_formula <- function(formula, time_var) {
  if (!is.null(formula)) {
    return(formula)
  }
  stats::as.formula(paste("~ factor(", .select_copula_formula_name(time_var), ")"))
}

.select_copula_formula_name <- function(name) {
  if (make.names(name) == name) {
    return(name)
  }
  paste0("`", gsub("`", "\\\\`", name), "`")
}

.select_copula_response_formula <- function(formula, response_var) {
  if (is.null(formula)) {
    return(stats::as.formula(paste(response_var, "~ 1")))
  }
  formula <- stats::as.formula(formula)
  if (length(formula) == 2L) {
    return(stats::as.formula(paste(response_var, deparse(formula), sep = " ")))
  }
  formula
}

.select_copula_rhs_formula <- function(formula) {
  if (is.null(formula)) {
    return(~1)
  }
  formula <- stats::as.formula(formula)
  if (length(formula) == 3L) {
    return(stats::as.formula(paste("~", deparse(formula[[3L]]))))
  }
  formula
}

.select_copula_u_from_fit <- function(object) {
  copula_link <- get_copula_dist(object$copula_dist)$copula_link
  eta_out <- calc_eta(
    object$par,
    object$model_matrix,
    object$margin_dist,
    copula_link,
    object$par_s
  )
  u <- calc_F_x(
    eta_out$eta_inv,
    object$model_matrix$x,
    object$margin_dist,
    object$response
  )
  .copula_clamp01(u)
}

.select_copula_adjacent_pairs <- function(data, subject_var, time_var, u_var, lags) {
  ord <- order(data[[subject_var]], data[[time_var]])
  data <- data[ord, , drop = FALSE]
  subjects <- unique(data[[subject_var]])
  out <- vector("list", length(subjects) * length(lags))
  k <- 0L

  for (subject in subjects) {
    subject_data <- data[data[[subject_var]] == subject, , drop = FALSE]
    subject_data <- subject_data[order(subject_data[[time_var]]), , drop = FALSE]
    n_time <- nrow(subject_data)
    for (lag in lags) {
      if (n_time <= lag) next
      left <- seq_len(n_time - lag)
      right <- left + lag
      k <- k + 1L
      out[[k]] <- data.frame(
        u1 = .copula_clamp01(subject_data[[u_var]][left]),
        u2 = .copula_clamp01(subject_data[[u_var]][right]),
        copula_time = paste(subject_data[[time_var]][left], subject_data[[time_var]][right], sep = "->")
      )
    }
  }

  if (k == 0L) {
    return(data.frame(u1 = numeric(), u2 = numeric()))
  }
  do.call(rbind, out[seq_len(k)])
}

.select_copula_fit_family <- function(u1, u2, family, t_df_grid, copula_time = NULL) {
  if (!is.null(copula_time)) {
    return(.select_copula_fit_family_by_time(u1, u2, family, t_df_grid, copula_time))
  }
  tau_start <- suppressWarnings(stats::cor(u1, u2, method = "kendall", use = "complete.obs"))
  if (!is.finite(tau_start)) tau_start <- 0

  fit <- switch(
    family,
    N = .select_copula_fit_one_par(u1, u2, family, lower = -0.95, upper = 0.95, start = .copula_tau_to_par("N", tau_start)),
    C = .select_copula_fit_one_par(u1, u2, family, lower = 1e-8, upper = 50, start = .copula_tau_to_par("C", pmax(tau_start, 0))),
    F = .select_copula_fit_frank(u1, u2, tau_start),
    G = .select_copula_fit_one_par(u1, u2, family, lower = 1 + 1e-8, upper = 50, start = .copula_tau_to_par("G", pmax(tau_start, 0))),
    J = .select_copula_fit_one_par(u1, u2, family, lower = 1 + 1e-8, upper = 50, start = .copula_tau_to_par("J", pmax(tau_start, 0))),
    t = .select_copula_fit_t(u1, u2, tau_start, t_df_grid),
    stop("Unsupported copula family: ", family, call. = FALSE)
  )

  k <- if (identical(family, "t")) 2 else 1
  data.frame(
    family = family,
    par = fit$par,
    par2 = fit$par2,
    tau = .copula_par_to_tau(family, fit$par, fit$par2),
    logLik = fit$logLik,
    AIC = -2 * fit$logLik + 2 * k,
    BIC = -2 * fit$logLik + log(length(u1)) * k,
    stringsAsFactors = FALSE
  )
}

.select_copula_fit_family_by_time <- function(u1, u2, family, t_df_grid, copula_time) {
  copula_time <- factor(copula_time, levels = unique(copula_time))
  if (length(copula_time) != length(u1)) {
    stop("'copula_time' must have one value per pseudo-observation pair.", call. = FALSE)
  }
  levels_time <- levels(copula_time)
  fits <- lapply(levels_time, function(level) {
    idx <- copula_time == level
    .select_copula_fit_family(u1[idx], u2[idx], family = family, t_df_grid = t_df_grid)
  })
  log_lik <- sum(vapply(fits, function(fit) fit$logLik[[1L]], numeric(1)))
  k_per_level <- if (identical(family, "t")) 2 else 1
  k_total <- length(levels_time) * k_per_level
  tau <- vapply(fits, function(fit) fit$tau[[1L]], numeric(1))
  data.frame(
    family = family,
    par = NA_real_,
    par2 = if (identical(family, "t")) NA_real_ else 0,
    tau = stats::weighted.mean(tau, w = as.numeric(table(copula_time)), na.rm = TRUE),
    logLik = log_lik,
    AIC = -2 * log_lik + 2 * k_total,
    BIC = -2 * log_lik + log(length(u1)) * k_total,
    n_copula_time_levels = length(levels_time),
    stringsAsFactors = FALSE
  )
}

.select_copula_fit_one_par <- function(u1, u2, family, lower, upper, start = NULL, par2 = 0) {
  lower <- as.numeric(lower)
  upper <- as.numeric(upper)
  if (!is.null(start) && is.finite(start)) {
    start <- pmin(pmax(as.numeric(start), lower), upper)
  }
  obj <- function(par) -.select_copula_loglik(u1, u2, family = family, par = par, par2 = par2)
  opt <- stats::optimize(obj, interval = c(lower, upper))
  candidates <- c(opt$minimum, lower, upper, start)
  candidates <- unique(candidates[is.finite(candidates)])
  ll <- vapply(candidates, function(par) -obj(par), numeric(1))
  best <- which.max(ll)
  list(par = candidates[best], par2 = par2, logLik = ll[best])
}

.select_copula_fit_frank <- function(u1, u2, tau_start) {
  if (abs(tau_start) < 1e-4) {
    tau_start <- 0.05
  }
  start <- .copula_tau_to_par("F", tau_start)
  .select_copula_fit_one_par(u1, u2, "F", lower = -50, upper = 50, start = start)
}

.select_copula_fit_t <- function(u1, u2, tau_start, t_df_grid) {
  t_df_grid <- unique(as.numeric(t_df_grid))
  t_df_grid <- t_df_grid[is.finite(t_df_grid) & t_df_grid > 2]
  if (length(t_df_grid) < 1L) {
    stop("t_df_grid must contain at least one finite value greater than 2.", call. = FALSE)
  }
  fits <- lapply(t_df_grid, function(df) {
    .select_copula_fit_one_par(
      u1,
      u2,
      family = "t",
      lower = -0.95,
      upper = 0.95,
      start = .copula_tau_to_par("t", tau_start),
      par2 = df
    )
  })
  ll <- vapply(fits, `[[`, numeric(1), "logLik")
  fits[[which.max(ll)]]
}

.select_copula_loglik <- function(u1, u2, family, par, par2 = 0) {
  dens <- .copula_pdf(u1, u2, family = family, par = par, par2 = par2)
  dens <- pmax(dens, .Machine$double.xmin)
  sum(log(dens[is.finite(dens)]))
}
