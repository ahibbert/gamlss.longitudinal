#' First and second derivatives at an interior or one-sided boundary point
#' @noRd
.copula_bounded_log_derivatives <- function(fn, x, h, lower = -Inf, upper = Inf) {
  values <- .copula_recycle(x, h, lower, upper)
  x <- values[[1L]]; h <- values[[2L]]
  lower <- values[[3L]]; upper <- values[[4L]]
  forward <- x - 2 * h < lower
  backward <- x + 2 * h > upper
  if (any(forward & backward)) stop("Copula curvature step exceeds its domain.", call. = FALSE)
  offsets <- matrix(rep(c(-2, -1, 0, 1, 2), each = length(x)), nrow = length(x))
  offsets[forward, ] <- matrix(rep(0:4, each = sum(forward)), nrow = sum(forward))
  offsets[backward, ] <- matrix(rep(-(0:4), each = sum(backward)), nrow = sum(backward))
  nodes <- x + offsets * h
  if (any(nodes < lower | nodes > upper)) stop("Copula curvature node leaves its domain.", call. = FALSE)
  centre <- fn(x)
  delta <- vapply(seq_len(5L), function(i) fn(nodes[, i]) - centre, numeric(length(x)))
  dim(delta) <- c(length(x), 5L)
  w1 <- matrix(rep(c(1, -8, 0, 8, -1), each = length(x)), nrow = length(x))
  w2 <- matrix(rep(c(-1, 16, -30, 16, -1), each = length(x)), nrow = length(x))
  w1[forward, ] <- matrix(rep(c(-25, 48, -36, 16, -3), each = sum(forward)), nrow = sum(forward))
  w1[backward, ] <- matrix(rep(c(25, -48, 36, -16, 3), each = sum(backward)), nrow = sum(backward))
  w2[forward | backward, ] <- matrix(rep(c(35, -104, 114, -56, 11),
    each = sum(forward | backward)), nrow = sum(forward | backward))
  list(first = rowSums(delta * w1) / (12 * h),
    second = rowSums(delta * w2) / (12 * h^2))
}

#' Stable native Archimedean log density for parameter differentiation
#' @noRd
.copula_parameter_log_density <- function(u, v, family, theta, nu) {
  if (family == "C") {
    out <- numeric(length(theta))
    k <- theta != 0
    a <- log(u[k]); b <- log(v[k]); t <- theta[k]
    log_s <- log1p(expm1(-t * a) + expm1(-t * b))
    out[k] <- log1p(t) - (t + 1) * (a + b) - (2 + 1 / t) * log_s
    return(out)
  }
  if (family == "F") {
    out <- numeric(length(theta))
    k <- theta != 0
    t <- theta[k]
    den <- expm1(-t) + expm1(-t * u[k]) * expm1(-t * v[k])
    out[k] <- log(-t * expm1(-t)) - t * (u[k] + v[k]) - 2 * log(abs(den))
    return(out)
  }
  .copula_logpdf(u, v, family, theta, nu)
}

#' Native Student-t correlation log-density derivatives, including its limit
#' @noRd
.copula_t_rho_log_curvature <- function(u, v, rho, nu) {
  cmp <- .copula_t_components(u, v, rho, nu)
  d <- cmp$rho_denom
  b <- cmp$df * cmp$rho + cmp$x * cmp$y
  first <- -(cmp$df + 1) * cmp$rho / d + (cmp$df + 2) * b / cmp$denom
  second <- -(cmp$df + 1) * (1 + cmp$rho^2) / d^2 +
    (cmp$df + 2) * (cmp$df / cmp$denom + 2 * b^2 / cmp$denom^2)
  gaussian <- cmp$df >= .copula_t_gaussian_limit_df
  if (any(gaussian)) {
    k <- gaussian
    density <- .copula_gaussian_pdf(cmp$u1[k], cmp$u2[k], cmp$rho[k])
    first[k] <- .copula_gaussian_deriv(cmp$u1[k], cmp$u2[k], cmp$rho[k], "par", log = TRUE)
    second[k] <- .copula_gaussian_deriv2(cmp$u1[k], cmp$u2[k], cmp$rho[k], "par") / density - first[k]^2
  }
  list(first = first, second = second)
}

#' Parameter-only log curvature for the shared continuous Hessian
#' @noRd
.copula_hessian_parameter_log_curvature <- function(u, v, family, par, par2, density) {
  values <- .copula_recycle(.copula_clamp01(u), .copula_clamp01(v), par, par2, density)
  u <- values[[1L]]; v <- values[[2L]]
  par <- values[[3L]]; par2 <- values[[4L]]; density <- values[[5L]]
  native <- identical(.copula_backend(), "native")
  if (!native || family == "N") {
    first <- .copula_deriv(u, v, family, par, par2, "par", log = TRUE)
    second <- .copula_deriv2(u, v, family, par, par2, "par") / density - first^2
  } else if (family == "t") {
    rho <- .copula_t_rho_log_curvature(u, v, par, par2)
    first <- rho$first; second <- rho$second
  } else {
    lower <- switch(family, C = 0, F = -Inf, G = 1, J = 1)
    par <- pmax(par, lower)
    result <- .copula_bounded_log_derivatives(function(p)
      .copula_parameter_log_density(u, v, family, p, par2),
      par, h = 1e-4 * pmax(1, abs(par)), lower = lower)
    first <- result$first; second <- result$second
  }
  out <- list(theta = first, theta2 = second)
  if (family != "t") return(out)
  if (!native) {
    out$zeta <- .copula_deriv(u, v, family, par, par2, "par2", log = TRUE)
    out$zeta2 <- .copula_deriv2(u, v, family, par, par2, "par2") / density - out$zeta^2
    out$thetazeta <- .copula_deriv2(u, v, family, par, par2, "par1par2") / density - first * out$zeta
    return(out)
  }
  df <- .copula_t_df(par2)
  domain <- .copula_t_df_derivative_domain(df)
  lower <- domain$lower
  upper <- domain$upper
  result <- .copula_bounded_log_derivatives(function(p)
    .copula_t_logpdf(u, v, par, p), df, 1e-3 * pmax(1, df), lower, upper)
  cross <- .copula_bounded_log_derivatives(function(p)
    .copula_t_rho_log_curvature(u, v, par, p)$first,
    df, 1e-3 * pmax(1, df), lower, upper)
  out$zeta <- result$first; out$zeta2 <- result$second
  out$thetazeta <- cross$first
  out
}
