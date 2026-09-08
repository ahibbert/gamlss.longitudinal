#' Stable copula-density average over a discrete probability rectangle
#'
#' The rectangle likelihood ratio is the average copula density over the two
#' marginal probability intervals. Fixed Gauss-Legendre quadrature evaluates
#' that average directly, avoiding cancellation in four-CDF subtraction. For a
#' Gaussian copula, upper-tail nodes are represented as log survival
#' probabilities and transformed directly to normal quantiles.
#'
#' @noRd
.gl_discrete_copula_log_ratio_quadrature <- function(
    family,
    par,
    par2,
    upper,
    lower,
    margin_log_d1,
    margin_log_d2,
    log_survival1 = NULL,
    log_survival2 = NULL) {
  family <- .copula_family_code(family)
  nodes_raw <- c(
    -0.9602898564975363, -0.7966664774136267,
    -0.5255324099163290, -0.1834346424956498,
     0.1834346424956498,  0.5255324099163290,
     0.7966664774136267,  0.9602898564975363
  )
  weights_raw <- c(
    0.1012285362903763, 0.2223810344533745,
    0.3137066458778873, 0.3626837833783620,
    0.3626837833783620, 0.3137066458778873,
    0.2223810344533745, 0.1012285362903763
  )
  nodes <- (nodes_raw + 1) / 2
  weights <- weights_raw / 2
  grid <- expand.grid(i = seq_along(nodes), j = seq_along(nodes))
  log_weights <- log(weights[grid$i]) + log(weights[grid$j])

  logspace_add <- function(a, b) {
    m <- pmax(a, b)
    m + log(exp(a - m) + exp(b - m))
  }
  log_sum_exp <- function(x) {
    m <- max(x)
    if (!is.finite(m)) return(m)
    m + log(sum(exp(x - m)))
  }
  node_probability <- function(k, node, log_survival, log_mass, lower_value, upper_value) {
    if (!is.null(log_survival) && is.finite(log_survival[[k]]) && is.finite(log_mass[[k]])) {
      log_s <- logspace_add(log_survival[[k]], log1p(-node) + log_mass[[k]])
      return(list(u = -expm1(log_s), log_survival = log_s))
    }
    list(
      u = lower_value + node * (upper_value - lower_value),
      log_survival = NA_real_
    )
  }

  vapply(seq_along(par), function(k) {
    log_density <- numeric(nrow(grid))
    for (g in seq_len(nrow(grid))) {
      p1 <- node_probability(
        k, nodes[grid$i[g]], log_survival1, margin_log_d1,
        lower[k, 1], upper[k, 1]
      )
      p2 <- node_probability(
        k, nodes[grid$j[g]], log_survival2, margin_log_d2,
        lower[k, 2], upper[k, 2]
      )
      if (identical(family, "N") && is.finite(p1$log_survival) && is.finite(p2$log_survival)) {
        z1 <- stats::qnorm(p1$log_survival, lower.tail = FALSE, log.p = TRUE)
        z2 <- stats::qnorm(p2$log_survival, lower.tail = FALSE, log.p = TRUE)
        rho <- .copula_gaussian_rho(par[[k]])
        denom <- 1 - rho^2
        log_density[g] <- (
          2 * rho * z1 * z2 - rho^2 * (z1^2 + z2^2)
        ) / (2 * denom) - 0.5 * log(denom)
      } else {
        log_density[g] <- .copula_logpdf(
          p1$u, p2$u, family = family, par = par[[k]], par2 = par2[[k]]
        )
      }
    }
    log_sum_exp(log_weights + log_density)
  }, numeric(1))
}
