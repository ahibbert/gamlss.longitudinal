.copula_frank_par_to_tau <- function(par) {
  theta <- .copula_frank_par(par)

  vapply(theta, function(th) {
    if (abs(th) <= 1e-8) {
      return(0)
    }

    integrand <- function(x) {
      ifelse(abs(x) < 1e-6, 1 - x / 2 + x^2 / 12, x / expm1(x))
    }

    deb <- stats::integrate(integrand, lower = 0, upper = th, rel.tol = 1e-8)$value / th

    1 - 4 / th + 4 * deb / th
  }, numeric(1), USE.NAMES = FALSE)
}

.copula_frank_tau_to_par <- function(tau) {
  tau <- pmin(pmax(as.numeric(tau), -0.999999), 0.999999)

  vapply(tau, function(tau_i) {
    if (abs(tau_i) <= 1e-8) {
      return(0)
    }

    lower <- if (tau_i < 0) -50 else 1e-8

    upper <- if (tau_i < 0) -1e-8 else 50

    stats::uniroot(

      function(th) .copula_frank_par_to_tau(th) - tau_i,
      lower = lower,
      upper = upper,
      tol = 1e-8
    )$root
  }, numeric(1), USE.NAMES = FALSE)
}
