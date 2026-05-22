#' Select a bivariate copula family from pseudo-observation pairs
#'
#' Screens supported copula families by maximising their native copula
#' log-likelihood on adjacent-time pseudo-observation pairs. Pseudo-observations
#' can be supplied directly as `u1`/`u2`, supplied as a row-aligned uniform
#' vector or column in `data`, or computed from a fitted
#' `gamlss.longitudinal` object.
#'
#' This helper is intended as a lightweight family-screening step. It estimates
#' a constant dependence parameter for each candidate family; richer covariate
#' or smooth dependence structures should be fitted afterwards with
#' [gamlss.longitudinal()].
#'
#' @param data Optional long-format data frame.
#' @param object Optional fitted `gamlss.longitudinal` object.
#' @param u1,u2 Optional vectors of paired pseudo-observations.
#' @param u Optional row-aligned vector of pseudo-observations for `data`.
#' @param u_var Optional name of a pseudo-observation column in `data`.
#' @param subject_var,time_var Subject and time column names used when building
#'   adjacent-time pairs from `data`.
#' @param families Candidate copula family codes. Supported values are `"N"`,
#'   `"C"`, `"F"`, `"G"`, `"J"`, and `"t"`.
#' @param lags Positive integer lag(s) used when forming adjacent pairs.
#' @param criterion Ranking criterion, one of `"AIC"`, `"BIC"`, or `"logLik"`.
#' @param t_df_grid Degrees-of-freedom grid used for the t-copula screen.
#' @param min_pairs Minimum number of complete pairs required.
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
  subject_var = "subject",
  time_var = "time",
  families = c("N", "C", "F", "G", "J", "t"),
  lags = 1,
  criterion = c("AIC", "BIC", "logLik"),
  t_df_grid = c(3, 4, 6, 8, 12, 20, 30),
  min_pairs = 10
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
    subject_var = subject_var,
    time_var = time_var,
    lags = lags
  )
  keep <- is.finite(pairs$u1) & is.finite(pairs$u2)
  pairs <- pairs[keep, , drop = FALSE]
  if (nrow(pairs) < min_pairs) {
    stop("At least ", min_pairs, " complete pseudo-observation pairs are required.", call. = FALSE)
  }

  fits <- lapply(families, function(family) {
    .select_copula_fit_family(
      u1 = pairs$u1,
      u2 = pairs$u2,
      family = family,
      t_df_grid = t_df_grid
    )
  })
  out <- do.call(rbind, fits)
  out$n_pairs <- nrow(pairs)
  out <- out[order(out[[criterion]], decreasing = identical(criterion, "logLik")), , drop = FALSE]
  rownames(out) <- NULL
  attr(out, "selected") <- out$family[1]
  attr(out, "criterion") <- criterion
  class(out) <- c("copula_selection", "data.frame")
  out
}

.select_copula_pairs <- function(data, object, u1, u2, u, u_var, subject_var, time_var, lags) {
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
    if (!is.null(u)) {
      if (length(u) != nrow(data)) {
        stop("u must have one value per row of data.", call. = FALSE)
      }
      data[[".u"]] <- u
      u_var <- ".u"
    }
    if (is.null(u_var)) {
      stop("Provide u_var or u when selecting from data.", call. = FALSE)
    }
  }

  if (!all(c(subject_var, time_var, u_var) %in% names(data))) {
    stop("data must contain subject_var, time_var, and u_var columns.", call. = FALSE)
  }
  .select_copula_adjacent_pairs(data, subject_var = subject_var, time_var = time_var, u_var = u_var, lags = lags)
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
        u2 = .copula_clamp01(subject_data[[u_var]][right])
      )
    }
  }

  if (k == 0L) {
    return(data.frame(u1 = numeric(), u2 = numeric()))
  }
  do.call(rbind, out[seq_len(k)])
}

.select_copula_fit_family <- function(u1, u2, family, t_df_grid) {
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
