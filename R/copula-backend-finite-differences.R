.copula_one_sided_par_deriv <- function(pdf_fun, u1, u2, par0, h, log = FALSE) {
  dens0 <- rep(1, length(u1))
  dens1 <- pdf_fun(u1, u2, par0 + h)
  if (isTRUE(log)) {
    (log(dens1) - log(dens0)) / h
  } else {
    (dens1 - dens0) / h
  }
}

.copula_central_par_deriv <- function(pdf_fun, u1, u2, par0, h, log = FALSE) {
  dens_plus <- pdf_fun(u1, u2, par0 + h)
  dens_minus <- pdf_fun(u1, u2, par0 - h)
  if (isTRUE(log)) {
    (log(dens_plus) - log(dens_minus)) / (2 * h)
  } else {
    (dens_plus - dens_minus) / (2 * h)
  }
}

.copula_one_sided_par_deriv2 <- function(pdf_fun, u1, u2, par0, h) {
  (pdf_fun(u1, u2, par0 + 2 * h) - 2 * pdf_fun(u1, u2, par0 + h) + 1) / h^2
}

.copula_central_par_deriv2 <- function(pdf_fun, u1, u2, par0, h) {
  (pdf_fun(u1, u2, par0 + h) - 2 + pdf_fun(u1, u2, par0 - h)) / h^2
}
