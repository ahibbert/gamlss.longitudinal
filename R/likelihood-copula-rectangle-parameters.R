.calc_copula_rectangle_par_derivative <- function(calc_lik, copula_dist, par_name, h = 1e-5) {
  row_id1 <- calc_lik$copula_row_id1
  n_pair <- length(row_id1)
  if (n_pair == 0L) {
    return(numeric(0))
  }

  par1 <- calc_lik$copula_par1
  par2 <- calc_lik$copula_par2
  hp <- h * pmax(1, abs(if (identical(par_name, "theta")) par1 else par2))
  par1_p <- par1_m <- par1
  par2_p <- par2_m <- par2
  if (identical(par_name, "theta")) {
    par1_p <- par1 + hp
    par1_m <- par1 - hp
    if (identical(copula_dist, "C")) par1_m <- pmax(par1_m, 1e-8)
  } else {
    par2_p <- par2 + hp
    par2_m <- pmax(par2 - hp, 1e-8)
  }

  u <- calc_lik$Fx_1_2
  l <- calc_lik$Fx_1_2_lower
  rect_p <- .copula_rectangle_prob(u[, 1], u[, 2], l[, 1], l[, 2], copula_dist, par1_p, par2_p)
  rect_m <- .copula_rectangle_prob(u[, 1], u[, 2], l[, 1], l[, 2], copula_dist, par1_m, par2_m)
  deriv <- (rect_p - rect_m) / (2 * hp)
  .gl_score_zero_unobserved_or_stop(
    deriv,
    calc_lik$pair_complete,
    paste0("Discrete rectangle derivative for copula parameter '", par_name, "'")
  )
}

#' Student-t rectangle derivatives from one shared quadrature grid
#'
#' The correlation derivative is integrated analytically. The degrees-of-
#' freedom derivative retains the existing central finite-difference semantics,
#' but reuses the probability nodes rather than rebuilding them for each side.
#'
#' @noRd
.calc_copula_t_rectangle_par_derivatives <- function(
    calc_lik, par_names = c("theta", "zeta"), h = 1e-5) {
  par_names <- intersect(as.character(par_names), c("theta", "zeta"))
  out <- setNames(vector("list", length(par_names)), par_names)
  n_pair <- length(calc_lik$copula_row_id1)
  if (!n_pair) {
    for (parameter in par_names) out[[parameter]] <- numeric(0)
    return(out)
  }

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
  n_grid <- nrow(grid)

  upper <- calc_lik$Fx_1_2
  lower <- calc_lik$Fx_1_2_lower
  width1 <- upper[, 1L] - lower[, 1L]
  width2 <- upper[, 2L] - lower[, 2L]
  rho <- .copula_gaussian_rho(calc_lik$copula_par1)
  df <- .copula_t_df(calc_lik$copula_par2)
  valid <- calc_lik$pair_complete & is.finite(width1) & is.finite(width2) &
    width1 > 0 & width2 > 0 & is.finite(rho) & is.finite(df)
  for (parameter in par_names) out[[parameter]] <- numeric(n_pair)
  if (!any(valid)) return(out)

  indices <- which(valid)
  n_valid <- length(indices)
  node1 <- matrix(nodes[grid$i], nrow = n_valid, ncol = n_grid, byrow = TRUE)
  node2 <- matrix(nodes[grid$j], nrow = n_valid, ncol = n_grid, byrow = TRUE)
  u1 <- matrix(lower[indices, 1L], nrow = n_valid, ncol = n_grid) +
    node1 * matrix(width1[indices], nrow = n_valid, ncol = n_grid)
  u2 <- matrix(lower[indices, 2L], nrow = n_valid, ncol = n_grid) +
    node2 * matrix(width2[indices], nrow = n_valid, ncol = n_grid)
  u1_vec <- as.numeric(u1)
  u2_vec <- as.numeric(u2)
  rho_vec <- rep(rho[indices], times = n_grid)
  df_vec <- rep(df[indices], times = n_grid)

  log_integral <- function(log_density) {
    weighted <- matrix(log_density, nrow = n_valid, ncol = n_grid)
    weighted <- sweep(weighted, 2L, log_weights, "+")
    maximum <- apply(weighted, 1L, max)
    result <- maximum
    finite <- is.finite(maximum)
    if (any(finite)) {
      result[finite] <- maximum[finite] + log(rowSums(exp(
        weighted[finite, , drop = FALSE] - maximum[finite]
      )))
    }
    result
  }

  base <- NULL
  if ("theta" %in% par_names) {
    base <- .copula_t_components(u1_vec, u2_vec, rho_vec, df_vec)
    dq_drho <- (
      -2 * base$x * base$y * base$rho_denom +
        2 * base$rho * base$numerator
    ) / base$rho_denom^2
    score_rho <- base$rho / base$rho_denom -
      (base$df + 2) * dq_drho / (2 * (base$df + base$q))
    weighted <- matrix(base$log_density, nrow = n_valid, ncol = n_grid)
    weighted <- sweep(weighted, 2L, log_weights, "+")
    maximum <- apply(weighted, 1L, max)
    normalized <- exp(weighted - maximum)
    normalized <- normalized / rowSums(normalized)
    rectangle_score <- rowSums(
      normalized * matrix(score_rho, nrow = n_valid, ncol = n_grid)
    )
    out$theta[indices] <- calc_lik$copula_rect_prob[indices] * rectangle_score
  }

  if ("zeta" %in% par_names) {
    hp <- h * pmax(1, abs(df[indices]))
    df_plus <- rep(df[indices] + hp, times = n_grid)
    df_minus <- rep(pmax(df[indices] - hp, 1e-8), times = n_grid)
    log_plus <- .copula_t_logpdf(u1_vec, u2_vec, rho_vec, df_plus)
    log_minus <- .copula_t_logpdf(u1_vec, u2_vec, rho_vec, df_minus)
    width_product <- width1[indices] * width2[indices]
    rectangle_plus <- width_product * exp(log_integral(log_plus))
    rectangle_minus <- width_product * exp(log_integral(log_minus))
    upper_bound <- pmin(width1[indices], width2[indices])
    rectangle_plus <- pmax(pmin(rectangle_plus, upper_bound), 1e-300)
    rectangle_minus <- pmax(pmin(rectangle_minus, upper_bound), 1e-300)
    out$zeta[indices] <- (rectangle_plus - rectangle_minus) / (2 * hp)
  }

  for (parameter in par_names) {
    out[[parameter]] <- .gl_score_zero_unobserved_or_stop(
      out[[parameter]],
      calc_lik$pair_complete,
      paste0("Discrete rectangle derivative for copula parameter '", parameter, "'")
    )
  }
  out
}
