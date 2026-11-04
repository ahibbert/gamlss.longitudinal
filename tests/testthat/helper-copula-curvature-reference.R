# Independent symbolic reference. This file must not call package copula
# density/derivative helpers. D() differentiates the defining log density;
# marginal quantile Jacobians convert N/t natural coordinates back to u,v.
copula_curvature_reference <- local({
  expressions <- list(
    N = quote(-log(1 - theta^2) / 2 +
      (2 * theta * x * y - theta^2 * (x^2 + y^2)) /
        (2 * (1 - theta^2))),
    t = quote(-log(1 - theta^2) / 2 - (nu + 2) / 2 *
      log(1 + (x^2 - 2 * theta * x * y + y^2) /
        (nu * (1 - theta^2))) +
      (nu + 1) / 2 * (log(1 + x^2 / nu) + log(1 + y^2 / nu))),
    C = quote(log(1 + theta) - (theta + 1) * (log(x) + log(y)) -
      (2 + 1 / theta) * log(x^(-theta) + y^(-theta) - 1)),
    F = quote(log(theta * (1 - exp(-theta))) - theta * (x + y) -
      2 * log(abs(exp(-theta) - 1 +
        (exp(-theta * x) - 1) * (exp(-theta * y) - 1)))),
    G = quote(-((-log(x))^theta + (-log(y))^theta)^(1 / theta) +
      (theta - 1) * (log(-log(x)) + log(-log(y))) +
      (1 / theta - 2) * log((-log(x))^theta + (-log(y))^theta) +
      log(((-log(x))^theta + (-log(y))^theta)^(1 / theta) + theta - 1) -
      log(x) - log(y)),
    J = quote((theta - 1) * (log(1 - x) + log(1 - y)) +
      (1 / theta - 2) * log((1 - x)^theta + (1 - y)^theta -
        (1 - x)^theta * (1 - y)^theta) +
      log(theta * ((1 - x)^theta + (1 - y)^theta -
        (1 - x)^theta * (1 - y)^theta) +
        (theta - 1) * (1 - (1 - x)^theta) * (1 - (1 - y)^theta)))
  )
  # abs() has no D rule. Frank's denominator is negative for theta>0
  # and positive for theta<0; its log derivative is the derivative of log
  # either sign. Substitute the sign as a constant before differentiating.
  expressions$F <- quote(log(theta * (1 - exp(-theta))) - theta * (x + y) -
    2 * log(sign_den * (exp(-theta) - 1 +
      (exp(-theta * x) - 1) * (exp(-theta * y) - 1))))
  derivatives <- lapply(expressions, function(e) {
    list(logdensity = e, x = D(e, "x"), y = D(e, "y"),
      xy = D(D(e, "x"), "y"), xx = D(D(e, "x"), "x"),
      theta = D(e, "theta"), xtheta = D(D(e, "x"), "theta"),
      ytheta = D(D(e, "y"), "theta"), thetatheta = D(D(e, "theta"), "theta"))
  })
  function(u, v, family, theta, nu = 5) {
    stopifnot(length(u) == 1L, length(v) == 1L, u > 0, u < 1,
      v > 0, v < 1, family %in% names(expressions))
    if (family %in% c("C", "F") && theta == 0) {
      score <- if (family == "C") (1 + log(u)) * (1 + log(v)) else
        (1 - 2 * u) * (1 - 2 * v) / 2
      mixed_u <- if (family == "C") (1 + log(v)) / u else -(1 - 2 * v)
      mixed_v <- if (family == "C") (1 + log(u)) / v else -(1 - 2 * u)
      return(c(density = 1, log_u = 0, log_v = 0, log_uv = 0,
        density_uv = 0, density_uu = 0, density_theta = score,
        density_utheta = mixed_u, density_vtheta = mixed_v,
        density_thetatheta = NA_real_))
    }
    sign_den <- if (theta > 0) -1 else 1
    x <- u
    y <- v
    jacobian_x <- jacobian_y <- 1
    if (family == "N") {
      x <- stats::qnorm(u)
      y <- stats::qnorm(v)
      jacobian_x <- 1 / stats::dnorm(x)
      jacobian_y <- 1 / stats::dnorm(y)
    } else if (family == "t") {
      x <- stats::qt(u, nu)
      y <- stats::qt(v, nu)
      jacobian_x <- 1 / stats::dt(x, nu)
      jacobian_y <- 1 / stats::dt(y, nu)
    }
    reference_environment <- environment()
    values <- vapply(derivatives[[family]], eval, numeric(1),
      envir = reference_environment)
    if (family %in% c("G", "J") && theta == 1) {
      # Exact independence prevents subtracting large, equal symbolic terms
      # in the upper tail. Parameter derivatives retain their limiting values.
      values[c("logdensity", "x", "y", "xy", "xx")] <- 0
    }
    # The t expression omits constants independent of x,y; restore the
    # bivariate/univariate Student normalizing constants for density scale.
    log_density <- values[["logdensity"]]
    if (family == "t") {
      log_density <- log_density + lgamma((nu + 2) / 2) + lgamma(nu / 2) -
        2 * lgamma((nu + 1) / 2)
    }
    density <- exp(log_density)
    log_u <- values[["x"]] * jacobian_x
    log_v <- values[["y"]] * jacobian_y
    log_uv <- values[["xy"]] * jacobian_x * jacobian_y
    jacobian_x_derivative <- if (family == "N") x * jacobian_x^2 else
      if (family == "t") (nu + 1) * x / (nu + x^2) * jacobian_x^2 else 0
    c(density = density, log_u = log_u, log_v = log_v,
      log_uv = log_uv, density_uv = density * (log_uv + log_u * log_v),
      density_uu = density * (log_u^2 + values[["xx"]] * jacobian_x^2 +
        values[["x"]] * jacobian_x_derivative),
      density_theta = density * values[["theta"]],
      density_utheta = density * (values[["xtheta"]] * jacobian_x +
        log_u * values[["theta"]]),
      density_vtheta = density * (values[["ytheta"]] * jacobian_y +
        log_v * values[["theta"]]),
      density_thetatheta = density * (values[["thetatheta"]] + values[["theta"]]^2))
  }
})
