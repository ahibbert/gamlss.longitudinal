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
  lower <- rep(-1L, n)
  upper <- rep(0L, n)
  upper_cdf <- .gl_simulation_discrete_cdf(family, upper, args, pfun)
  unresolved <- upper_cdf < p

  while (any(unresolved)) {
    lower[unresolved] <- upper[unresolved]
    upper[unresolved] <- pmin(
      2L * upper[unresolved] + 1L,
      support_max[unresolved]
    )
    index <- which(unresolved)
    upper_cdf[index] <- .gl_simulation_discrete_cdf(
      family, upper[index], lapply(args, `[`, index), pfun
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
      family, middle, lapply(args, `[`, index), pfun
    )
    move_upper <- middle_cdf >= p[index]
    upper[index[move_upper]] <- middle[move_upper]
    lower[index[!move_upper]] <- middle[!move_upper]
    active <- upper - lower > 1L
  }

  as.numeric(upper)
}

.gl_simulation_discrete_cdf <- function(family, q, args, pfun = NULL) {
  # Registered count margins have nonnegative support. Use the same CDF
  # evaluation (including count fast paths) as the rectangle likelihood.
  input <- c(list(q = q, x = q, y = q), args)
  .gl_evaluate_margin_cdf(
    input = input, margin_dist = list(family = family), response = q,
    FUN = pfun, use_fast = TRUE
  )
}
