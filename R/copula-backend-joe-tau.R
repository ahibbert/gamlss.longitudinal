.copula_joe_par_to_tau <- function(par) {
  theta <- .copula_joe_theta(par)

  vapply(theta, function(th) {
    if (th <= 1 + 1e-8) {
      return(0)
    }

    phi <- function(u) -log1p(-(1 - u)^th)

    dphi <- function(u) -th * (1 - u)^(th - 1) / (1 - (1 - u)^th)

    integrand <- function(u) phi(u) / dphi(u)

    1 + 4 * stats::integrate(integrand, lower = 0, upper = 1, rel.tol = 1e-8)$value
  }, numeric(1), USE.NAMES = FALSE)
}

.copula_joe_tau_to_par <- function(tau) {
  tau <- pmin(pmax(as.numeric(tau), 0), 0.999999)

  vapply(tau, function(tau_i) {
    if (tau_i <= 1e-8) {
      return(1)
    }

    stats::uniroot(

      function(th) .copula_joe_par_to_tau(th) - tau_i,
      lower = 1 + 1e-8,
      upper = 100,
      tol = 1e-8
    )$root
  }, numeric(1), USE.NAMES = FALSE)
}
