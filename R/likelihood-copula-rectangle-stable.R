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

#' Stable Student-t conditional interval probability
#'
#' @noRd
.copula_t_conditional_interval <- function(upper, lower, x, rho, df) {
  scale <- sqrt((df + x^2) * (1 - rho^2) / (df + 1))
  upper_q <- stats::qt(upper, df = df)
  lower_q <- stats::qt(lower, df = df)
  z_upper <- (upper_q - rho * x) / scale
  z_lower <- (lower_q - rho * x) / scale

  # Subtract on the shorter tail. This avoids losing the count-cell mass when
  # both conditional CDF values round close to one in the ACTG upper tail.
  use_upper_tail <- z_lower > 0
  out <- numeric(length(z_upper))
  if (any(!use_upper_tail)) {
    ordinary <- !use_upper_tail
    out[ordinary] <- stats::pt(z_upper[ordinary], df = df[ordinary] + 1) -
      stats::pt(z_lower[ordinary], df = df[ordinary] + 1)
  }
  if (any(use_upper_tail)) {
    out[use_upper_tail] <-
      stats::pt(
        z_lower[use_upper_tail], df = df[use_upper_tail] + 1,
        lower.tail = FALSE
      ) -
      stats::pt(
        z_upper[use_upper_tail], df = df[use_upper_tail] + 1,
        lower.tail = FALSE
      )
  }
  pmin(pmax(out, 0), 1)
}

#' Direct Student-t copula rectangle probability by conditional quadrature
#'
#' A copula rectangle is the one-dimensional integral of a conditional CDF
#' difference. An eight-node Gauss-Legendre rule therefore needs eight rather
#' than 8 x 8 density evaluations while retaining the same t-copula rectangle
#' likelihood. Each pair is integrated over its narrower probability interval
#' to reduce curvature, and tail CDFs are subtracted on the stable side.
#'
#' @noRd
.copula_t_rectangle_prob_conditional <- function(
    u1, u2, l1, l2, par, par2, floor_probability = TRUE) {
  vals <- .copula_recycle(
    as.numeric(u1), as.numeric(u2), as.numeric(l1), as.numeric(l2),
    .copula_gaussian_rho(par), .copula_t_df(par2)
  )
  n <- length(vals[[1L]])
  if (!n) return(numeric(0))

  width1 <- vals[[1L]] - vals[[3L]]
  width2 <- vals[[2L]] - vals[[4L]]
  valid <- is.finite(width1) & is.finite(width2) &
    is.finite(vals[[5L]]) & is.finite(vals[[6L]]) &
    vals[[3L]] >= 0 & vals[[1L]] <= 1 &
    vals[[4L]] >= 0 & vals[[2L]] <= 1 &
    width1 >= 0 & width2 >= 0
  out <- rep(NA_real_, n)
  zero <- valid & (width1 == 0 | width2 == 0)
  out[zero] <- 0

  gaussian <- valid & !zero & vals[[6L]] >= .copula_t_gaussian_limit_df
  if (any(gaussian)) {
    out[gaussian] <- .copula_gaussian_rectangle_prob_cached(
      vals[[1L]][gaussian], vals[[2L]][gaussian],
      vals[[3L]][gaussian], vals[[4L]][gaussian], vals[[5L]][gaussian],
      floor_probability = FALSE
    )
  }

  ordinary <- valid & !zero & !gaussian
  if (any(ordinary)) {
    nodes <- (c(
      -0.9602898564975363, -0.7966664774136267,
      -0.5255324099163290, -0.1834346424956498,
       0.1834346424956498,  0.5255324099163290,
       0.7966664774136267,  0.9602898564975363
    ) + 1) / 2
    weights <- c(
      0.1012285362903763, 0.2223810344533745,
      0.3137066458778873, 0.3626837833783620,
      0.3626837833783620, 0.3137066458778873,
      0.2223810344533745, 0.1012285362903763
    ) / 2
    index <- which(ordinary)
    swap <- width2[index] < width1[index]
    integration_lower <- ifelse(swap, vals[[4L]][index], vals[[3L]][index])
    integration_width <- ifelse(swap, width2[index], width1[index])
    conditional_lower <- ifelse(swap, vals[[3L]][index], vals[[4L]][index])
    conditional_upper <- ifelse(swap, vals[[1L]][index], vals[[2L]][index])
    n_ordinary <- length(index)
    n_node <- length(nodes)
    integration_u <-
      matrix(integration_lower, nrow = n_ordinary, ncol = n_node) +
      matrix(integration_width, nrow = n_ordinary, ncol = n_node) *
      matrix(nodes, nrow = n_ordinary, ncol = n_node, byrow = TRUE)
    df <- rep(vals[[6L]][index], times = n_node)
    rho <- rep(vals[[5L]][index], times = n_node)
    conditional_mass <- .copula_t_conditional_interval(
      upper = rep(conditional_upper, times = n_node),
      lower = rep(conditional_lower, times = n_node),
      x = stats::qt(as.numeric(integration_u), df = df),
      rho = rho,
      df = df
    )
    conditional_mass <- matrix(
      conditional_mass, nrow = n_ordinary, ncol = n_node
    )
    probability <- integration_width * as.numeric(conditional_mass %*% weights)
    out[index] <- pmin(
      pmax(probability, 0),
      pmin(width1[index], width2[index])
    )
  }

  if (isTRUE(floor_probability)) pmax(out, 1e-300) else out
}

#' Direct Student-t copula rectangle probability by fixed 2D quadrature
#'
#' Evaluating a discrete Student-t rectangle as four separately integrated CDF
#' corners is both expensive and vulnerable to cancellation. Integrate the
#' copula density over the probability rectangle directly instead. The same
#' fixed Gauss-Legendre rule is used by the stable log-ratio fallback above.
#'
#' @noRd
.copula_t_rectangle_prob_quadrature_2d <- function(
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
