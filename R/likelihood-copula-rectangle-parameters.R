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
  deriv[!is.finite(deriv)] <- 0
  deriv
}
