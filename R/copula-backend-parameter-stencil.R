#' Differentiate without shifting the evaluation point at a domain boundary
#' @noRd
.copula_bounded_first_difference <- function(fn, x, h, lower = -Inf, upper = Inf) {
  vals <- .copula_recycle(x, h, lower, upper)
  x <- vals[[1L]]
  h <- vals[[2L]]
  lower <- vals[[3L]]
  upper <- vals[[4L]]
  forward <- x - h < lower
  backward <- x + h > upper
  if (any(forward & backward)) stop("Copula derivative step exceeds its domain.", call. = FALSE)
  # Each row stays inside the backend domain; one-sided formulas evaluate
  # the derivative at x, unlike an asymmetric secant through a moved centre.
  x1 <- ifelse(backward, x - h, x + h)
  x2 <- ifelse(forward, x + 2 * h, x - h)
  x2[backward] <- x[backward] - 2 * h[backward]
  if (any(x1 < lower | x1 > upper | x2 < lower | x2 > upper)) {
    stop("Copula derivative nodes leave its domain.", call. = FALSE)
  }
  f0 <- fn(x)
  f1 <- fn(x1)
  f2 <- fn(x2)
  out <- (f1 - f2) / (2 * h)
  out[forward] <- (4 * (f1[forward] - f0[forward]) -
    (f2[forward] - f0[forward])) / (2 * h[forward])
  out[backward] <- (-4 * (f1[backward] - f0[backward]) +
    (f2[backward] - f0[backward])) / (2 * h[backward])
  out
}

#' Keep native df stencils within their Student-t or Gaussian branch
#' @noRd
.copula_t_df_derivative_domain <- function(df) {
  limit <- .copula_t_gaussian_limit_df
  if (any(df == limit)) stop("Student-t df is at the nondifferentiable Gaussian approximation join.", call. = FALSE)
  below_limit <- limit - .Machine$double.eps * 2^floor(log2(limit))
  list(lower = ifelse(df > limit, limit, 2.000001),
    # The native Gaussian branch includes the join itself. Keep all nodes on
    # the Student-t side strictly below that point, including endpoint nodes.
    upper = ifelse(df < limit, below_limit, Inf))
}

#' Argument/parameter mixed derivatives at the supplied copula parameters
#' @noRd
.copula_parameter_argument_mixed <- function(u1, u2, family, par, par2 = 0,
                                             parameter = "par") {
  family <- .copula_family_code(family)
  if (identical(.copula_backend(), "vinecopula")) {
    key <- if (parameter == "par") "par1" else "par2"
    return(list(
      u1 = .copula_deriv2(u1, u2, family, par, par2, deriv = paste0(key, "u1")),
      u2 = .copula_deriv2(u1, u2, family, par, par2, deriv = paste0(key, "u2"))))
  }
  vals <- .copula_recycle(u1, u2, par, par2)
  u1 <- vals[[1L]]; u2 <- vals[[2L]]
  par <- vals[[3L]]; par2 <- vals[[4L]]
  if (parameter == "par2") {
    x <- .copula_t_df(par2)
    domain <- .copula_t_df_derivative_domain(x)
    lower <- domain$lower
    upper <- domain$upper
    h <- 1e-5 * pmax(1, abs(x))
  } else {
    lower <- switch(family, N = -0.999999, t = -0.999999, C = 0, G = 1, J = 1, F = -Inf)
    upper <- if (family %in% c("N", "t")) 0.999999 else Inf
    x <- pmin(pmax(par, lower), upper)
    h <- 1e-5 * pmax(1, abs(x))
    if (family %in% c("N", "t")) h <- pmin(h, 1e-4 * (1 - abs(x)))
  }
  evaluate <- function(p, component) {
    p1 <- if (parameter == "par") p else par
    p2 <- if (parameter == "par2") p else par2
    .copula_density_argument_curvature(u1, u2, family, p1, p2)[[component]]
  }
  list(u1 = .copula_bounded_first_difference(function(p) evaluate(p, "du"), x, h, lower, upper),
    u2 = .copula_bounded_first_difference(function(p) evaluate(p, "dv"), x, h, lower, upper))
}
