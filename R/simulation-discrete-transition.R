#' Redraw the outgoing copula coordinate for an observed count
#'
#' Conditional on the count, a fresh uniform over its marginal CDF interval
#' integrates the inverse h-function against the copula rectangle transition.
#' A coordinate inherited from the incoming edge would retain latent history.
#' @noRd
.gl_simulation_redraw_discrete <- function(lower, upper, family, rows) {
  midpoint <- lower + (upper - lower) / 2
  invalid <- !is.finite(midpoint) | lower < 0 | upper > 1 |
    midpoint <= lower | midpoint >= upper
  if (any(invalid)) {
    stop(
      "Cannot resolve an interior ", family, " CDF interval for discrete ",
      "simulation at row(s) ",
      paste(utils::head(unique(rows[invalid]), 10L), collapse = ", "), ".",
      call. = FALSE
    )
  }
  u <- stats::runif(length(lower), min = lower, max = upper)
  # A positive-width interval near 0 or 1 can round an interior draw
  # to an endpoint. The midpoint is already known to be representable.
  rounded <- u <= 0 | u >= 1
  u[rounded] <- midpoint[rounded]
  u
}

.gl_simulation_discrete_previous_u <- function(u, family, args, rows) {
  count <- .gl_simulation_quantiles(family, u, args)
  pfun <- .gl_get_family_fun(family, "p")
  lower <- .gl_simulation_discrete_cdf(family, count - 1, args, pfun)
  upper <- .gl_simulation_discrete_cdf(family, count, args, pfun)
  .gl_simulation_redraw_discrete(lower, upper, family, rows)
}

.gl_simulation_discrete_inverse <- function(u1, target, family, par, par2) {
  values <- .copula_recycle(u1, target, par, par2)
  u1 <- values[[1L]]
  target <- values[[2L]]
  par <- values[[3L]]
  par2 <- values[[4L]]
  if (any(!is.finite(c(u1, target, par))) || any(u1 <= 0 | u1 >= 1) ||
      any(target <= 0 | target >= 1)) {
    stop("Discrete inverse h-function requires finite parameters and open-unit probabilities.",
         call. = FALSE)
  }

  # The continuous simulator clips to sqrt(eps). That can erase whole
  # count categories, so use the analytic inverses on open probabilities.
  out <- if (.copula_backend() == "native") {
    switch(family,
      N = .sim_invert_gaussian_hfunc1(u1, target, par),
      C = .gl_simulation_clayton_open_inverse(u1, target, par),
      F = .sim_invert_frank_hfunc1(u1, target, par),
      t = .sim_invert_t_hfunc1(u1, target, par, par2),
      G = .gl_simulation_archimedean_open_inverse(u1, target, family, par, par2),
      J = .gl_simulation_archimedean_open_inverse(u1, target, family, par, par2)
    )
  } else {
    .copula_require_vinecopula("the delegated VineCopula inverse h-function backend")
    VineCopula::BiCopHinv1(u1, target, family = .copula_family_number(family),
                          par = par, par2 = par2)
  }
  if (any(!is.finite(out) | out < 0 | out > 1)) {
    stop("Invalid copula inverse h-function result during discrete simulation.", call. = FALSE)
  }

  # Keep only floating-point endpoints away from zero and one.
  out[out == 0] <- .Machine$double.xmin
  out[out == 1] <- 1 - .Machine$double.eps / 2
  out
}

.gl_simulation_clayton_open_inverse <- function(u1, target, par) {
  theta <- .copula_clayton_theta(par)
  dependent <- theta > 1e-10
  if (!any(dependent)) return(target)
  th <- theta[dependent]
  power <- -th / (th + 1) * log(target[dependent])
  log_term <- -th * log(u1[dependent]) + power + log(-expm1(-power))
  log_sum <- pmax(log_term, 0) + log1p(exp(-abs(log_term)))
  target[dependent] <- exp(-log_sum / th)
  target
}

.gl_simulation_archimedean_open_inverse <- function(u1, target, family, par, par2) {
  eps <- sqrt(.Machine$double.eps)
  ordinary <- u1 >= 1e-12 & u1 <= 1 - 1e-12 & target >= eps & target <= 1 - eps
  out <- target
  out[ordinary] <- .sim_invert_archimedean_hfunc1(
    u1[ordinary], target[ordinary], family, par[ordinary], par2[ordinary]
  )
  tail <- par > 1 + 1e-10 & (!ordinary | out <= eps | out >= 1 - eps)
  if (!any(tail)) return(out)
  # Logit bisection resolves probabilities below the ordinary absolute
  # tolerance and avoids the clamped diagnostic h-functions in the tails.
  lower <- rep(stats::qlogis(.Machine$double.xmin), sum(tail))
  upper <- rep(stats::qlogis(1 - .Machine$double.eps / 2), sum(tail))
  for (i in seq_len(48L)) {
    mid <- lower + (upper - lower) / 2
    h <- .gl_simulation_archimedean_log_h(u1[tail], stats::plogis(mid), family, par[tail])
    below <- h < log(target[tail])
    lower[below] <- mid[below]
    upper[!below] <- mid[!below]
  }
  out[tail] <- stats::plogis(lower + (upper - lower) / 2)
  out
}

.gl_simulation_archimedean_log_h <- function(u, v, family, theta) {
  log_sum_exp <- function(a, b) pmax(a, b) + log1p(exp(-abs(a - b)))
  if (family == "G") {
    x <- -log(u)
    log_x <- log(x)
    log_s <- log_sum_exp(theta * log_x, theta * log(-log(v))) / theta
    return(-x * expm1(log_s - log_x) + (theta - 1) * (log_x - log_s))
  }
  log_a <- theta * log1p(-u)
  log_b <- theta * log1p(-v)
  log_s <- log_sum_exp(log_a, log_b + log(-expm1(log_a)))
  (theta - 1) * log1p(-u) + log(-expm1(log_b)) + (1 / theta - 1) * log_s
}
