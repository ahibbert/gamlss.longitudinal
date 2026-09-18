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
  n_rect <- length(par)
  n_grid <- nrow(grid)
  if (!n_rect) return(numeric(0))

  node_matrix <- function(node, lower_value, upper_value, log_survival, log_mass) {
    node <- matrix(node, nrow = n_rect, ncol = n_grid, byrow = TRUE)
    lower_value <- matrix(lower_value, nrow = n_rect, ncol = n_grid)
    upper_value <- matrix(upper_value, nrow = n_rect, ncol = n_grid)
    u <- lower_value + node * (upper_value - lower_value)
    log_s <- matrix(NA_real_, nrow = n_rect, ncol = n_grid)
    stable <- rep(FALSE, n_rect)
    if (!is.null(log_survival) && !is.null(log_mass)) {
      stable <- is.finite(log_survival) & is.finite(log_mass)
      if (any(stable)) {
        base <- matrix(log_survival[stable], nrow = sum(stable), ncol = n_grid)
        increment <- matrix(log_mass[stable], nrow = sum(stable), ncol = n_grid) +
          log1p(-node[stable, , drop = FALSE])
        maximum <- pmax(base, increment)
        stable_log_s <- maximum +
          log(exp(base - maximum) + exp(increment - maximum))
        log_s[stable, ] <- stable_log_s
        u[stable, ] <- -expm1(stable_log_s)
      }
    }
    list(u = u, log_survival = log_s, stable = stable)
  }

  p1 <- node_matrix(
    nodes[grid$i], lower[, 1], upper[, 1], log_survival1, margin_log_d1
  )
  p2 <- node_matrix(
    nodes[grid$j], lower[, 2], upper[, 2], log_survival2, margin_log_d2
  )
  log_density <- matrix(NA_real_, nrow = n_rect, ncol = n_grid)

  gaussian_stable <- identical(family, "N") & p1$stable & p2$stable
  if (any(gaussian_stable)) {
    z1 <- matrix(
      stats::qnorm(
        p1$log_survival[gaussian_stable, , drop = FALSE],
        lower.tail = FALSE,
        log.p = TRUE
      ),
      nrow = sum(gaussian_stable),
      ncol = n_grid
    )
    z2 <- matrix(
      stats::qnorm(
        p2$log_survival[gaussian_stable, , drop = FALSE],
        lower.tail = FALSE,
        log.p = TRUE
      ),
      nrow = sum(gaussian_stable),
      ncol = n_grid
    )
    rho <- .copula_gaussian_rho(par[gaussian_stable])
    denom <- 1 - rho^2
    log_density[gaussian_stable, ] <- (
      2 * rho * z1 * z2 - rho^2 * (z1^2 + z2^2)
    ) / (2 * denom) - 0.5 * log(denom)
  }

  generic <- !gaussian_stable
  if (any(generic)) {
    log_density[generic, ] <- matrix(
      .copula_logpdf(
        as.numeric(p1$u[generic, , drop = FALSE]),
        as.numeric(p2$u[generic, , drop = FALSE]),
        family = family,
        par = rep(par[generic], times = n_grid),
        par2 = rep(par2[generic], times = n_grid)
      ),
      nrow = sum(generic),
      ncol = n_grid
    )
  }

  weighted <- sweep(log_density, 2L, log_weights, "+")
  maximum <- apply(weighted, 1L, max)
  out <- maximum
  finite <- is.finite(maximum)
  if (any(finite)) {
    centered <- weighted[finite, , drop = FALSE] - maximum[finite]
    out[finite] <- maximum[finite] + log(rowSums(exp(centered)))
  }
  out
}

#' Direct Student-t copula rectangle probability by fixed quadrature
#'
#' Evaluating a discrete Student-t rectangle as four separately integrated CDF
#' corners is both expensive and vulnerable to cancellation. Integrate the
#' copula density over the probability rectangle directly instead. The same
#' fixed Gauss-Legendre rule is used by the stable log-ratio fallback above.
#'
#' @noRd
.copula_t_rectangle_prob_quadrature <- function(
    u1, u2, l1, l2, par, par2, floor_probability = TRUE) {
  vals <- .copula_recycle(
    as.numeric(u1), as.numeric(u2), as.numeric(l1), as.numeric(l2),
    .copula_gaussian_rho(par), .copula_t_df(par2)
  )
  n <- length(vals[[1L]])
  if (!n) return(numeric(0))

  upper <- cbind(vals[[1L]], vals[[2L]])
  lower <- cbind(vals[[3L]], vals[[4L]])
  width1 <- upper[, 1L] - lower[, 1L]
  width2 <- upper[, 2L] - lower[, 2L]
  valid <- is.finite(width1) & is.finite(width2) &
    is.finite(vals[[5L]]) & is.finite(vals[[6L]]) &
    lower[, 1L] >= 0 & upper[, 1L] <= 1 &
    lower[, 2L] >= 0 & upper[, 2L] <= 1 &
    width1 >= 0 & width2 >= 0

  out <- rep(NA_real_, n)
  zero <- valid & (width1 == 0 | width2 == 0)
  out[zero] <- 0
  positive <- valid & width1 > 0 & width2 > 0
  if (any(positive)) {
    log_width1 <- log(width1[positive])
    log_width2 <- log(width2[positive])
    log_ratio <- .gl_discrete_copula_log_ratio_quadrature(
      family = "t",
      par = vals[[5L]][positive],
      par2 = vals[[6L]][positive],
      upper = upper[positive, , drop = FALSE],
      lower = lower[positive, , drop = FALSE],
      margin_log_d1 = log_width1,
      margin_log_d2 = log_width2
    )
    probability <- exp(log_ratio + log_width1 + log_width2)
    probability <- pmin(pmax(probability, 0), pmin(width1[positive], width2[positive]))
    out[positive] <- probability
  }

  if (isTRUE(floor_probability)) pmax(out, 1e-300) else out
}
