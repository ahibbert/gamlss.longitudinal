#' Evaluate marginal quantiles efficiently during simulation
#'
#' Continuous margins use their `gamlss.dist` quantile function directly.
#' Discrete margins use vectorised integer inversion of the corresponding CDF;
#' this avoids several `gamlss.dist` count-family quantile implementations that
#' scan the support separately for every simulated value.
#'
#' @noRd
.gl_simulation_quantiles <- function(family, p, args = list(), max_value = 1048575L) {
  family <- as.character(family)[1L]
  route <- .gl_capability_likelihood_route(family)

  if (identical(route, "exact_discrete_rectangle")) {
    return(.gl_simulation_discrete_quantiles(
      family = family,
      p = p,
      args = args,
      max_value = max_value
    ))
  }

  qfun <- get(
    paste0("q", family),
    envir = asNamespace("gamlss.dist"),
    inherits = FALSE
  )
  qargs <- c(list(p = p), args)
  qargs <- qargs[names(qargs) %in% formalArgs(qfun)]
  do.call(qfun, qargs)
}

.gl_simulation_discrete_quantiles <- function(family, p, args, max_value) {
  p <- as.numeric(p)
  if (any(!is.finite(p)) || any(p < 0 | p > 1)) {
    stop("Simulation probabilities must be finite and lie inside [0, 1].", call. = FALSE)
  }

  n <- length(p)
  if (!n) return(numeric())

  args <- lapply(args, rep, length.out = n)
  support_max <- if ("bd" %in% names(args)) {
    pmin(as.integer(args$bd), as.integer(max_value))
  } else {
    rep(as.integer(max_value), n)
  }
  pfun <- get(
    paste0("p", family),
    envir = asNamespace("gamlss.dist"),
    inherits = FALSE
  )
  cdf_below_zero <- .gl_simulation_discrete_cdf(pfun, rep(-1L, n), args)
  if (any(cdf_below_zero > 0)) {
    qfun <- get(
      paste0("q", family),
      envir = asNamespace("gamlss.dist"),
      inherits = FALSE
    )
    qargs <- c(list(p = p), args)
    qargs <- qargs[names(qargs) %in% formalArgs(qfun)]
    return(do.call(qfun, qargs))
  }
  lower <- rep(-1L, n)
  upper <- rep(0L, n)
  upper_cdf <- .gl_simulation_discrete_cdf(pfun, upper, args)
  unresolved <- upper_cdf < p

  while (any(unresolved)) {
    lower[unresolved] <- upper[unresolved]
    upper[unresolved] <- pmin(
      2L * upper[unresolved] + 1L,
      support_max[unresolved]
    )
    index <- which(unresolved)
    upper_cdf[index] <- .gl_simulation_discrete_cdf(
      pfun,
      upper[index],
      lapply(args, `[`, index)
    )
    unresolved <- upper_cdf < p
    if (any(unresolved & upper >= support_max)) {
      stop(
        "Could not bracket simulated ", family,
        " quantiles below the internal maximum of ", max_value, ".",
        call. = FALSE
      )
    }
  }

  active <- upper - lower > 1L
  while (any(active)) {
    index <- which(active)
    middle <- lower[index] + (upper[index] - lower[index]) %/% 2L
    middle_cdf <- .gl_simulation_discrete_cdf(
      pfun,
      middle,
      lapply(args, `[`, index)
    )
    move_upper <- middle_cdf >= p[index]
    upper[index[move_upper]] <- middle[move_upper]
    lower[index[!move_upper]] <- middle[!move_upper]
    active <- upper - lower > 1L
  }

  as.numeric(upper)
}

.gl_simulation_discrete_cdf <- function(pfun, q, args) {
  pargs <- c(list(q = q), args)
  pargs <- pargs[names(pargs) %in% formalArgs(pfun)]
  value <- do.call(pfun, pargs)
  if (any(!is.finite(value))) {
    stop("Non-finite marginal CDF encountered during simulation.", call. = FALSE)
  }
  as.numeric(value)
}
