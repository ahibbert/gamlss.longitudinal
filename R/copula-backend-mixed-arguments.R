#' Mixed argument derivative of a continuous copula density
#'
#' Evaluate at the same probabilities as the density backend. A fixed-width
#' probability stencil is unsuitable in the tails: moving its centre changes
#' the point where curvature is evaluated. If L = log(c), then
#' c_uv = c * (L_uv + L_u * L_v).
#' @noRd
.copula_density_argument_curvature <- function(u1, u2, family, par, par2 = 0) {
  family <- .copula_family_code(family)
  native <- identical(.copula_backend(), "native")
  # VineCopula's density also clamps to [1e-12,1-1e-12]. Keep the
  # analytical curvature at that effective point for either backend.
  u1 <- .copula_clamp01(u1)
  u2 <- .copula_clamp01(u2)
  if (native) {
    par <- switch(family,
      N = .copula_gaussian_rho(par), t = .copula_gaussian_rho(par),
      C = .copula_clayton_theta(par), F = .copula_frank_par(par),
      G = .copula_gumbel_theta(par), J = .copula_joe_theta(par))
    if (family == "t") par2 <- .copula_t_df(par2)
  }
  values <- .copula_recycle(u1, u2, par, par2)
  u <- values[[1L]]
  v <- values[[2L]]
  theta <- values[[3L]]
  nu <- values[[4L]]
  density <- .copula_pdf(u, v, family, theta, nu)
  su <- sv <- suv <- suu <- svv <- rep(NA_real_, length(u))

  if (family == "N" || family == "t") {
    gaussian <- rep(family == "N", length(u))
    if (native && family == "t") gaussian <- nu >= .copula_t_gaussian_limit_df
    if (any(gaussian)) {
      k <- gaussian
      x <- stats::qnorm(u[k])
      y <- stats::qnorm(v[k])
      rho <- theta[k]
      den <- 1 - rho^2
      px <- stats::dnorm(x)
      py <- stats::dnorm(y)
      su[k] <- (rho * y - rho^2 * x) / (den * px)
      sv[k] <- (rho * x - rho^2 * y) / (den * py)
      suv[k] <- rho / (den * px * py)
      suu[k] <- (-rho^2 + (rho * y - rho^2 * x) * x) / (den * px^2)
      svv[k] <- (-rho^2 + (rho * x - rho^2 * y) * y) / (den * py^2)
    }
    if (any(!gaussian)) {
      k <- !gaussian
      df <- nu[k]
      x <- stats::qt(u[k], df)
      y <- stats::qt(v[k], df)
      rho <- theta[k]
      den <- df * (1 - rho^2) + x^2 - 2 * rho * x * y + y^2
      px <- stats::dt(x, df)
      py <- stats::dt(y, df)
      a <- x - rho * y
      b <- y - rho * x
      su[k] <- (-(df + 2) * a / den + (df + 1) * x / (df + x^2)) / px
      sv[k] <- (-(df + 2) * b / den + (df + 1) * y / (df + y^2)) / py
      suv[k] <- (df + 2) * (rho / den + 2 * a * b / den^2) / (px * py)
      lx <- su[k] * px
      ly <- sv[k] * py
      suu[k] <- (-(df + 2) * (1 / den - 2 * a^2 / den^2) +
        (df + 1) * (df - x^2) / (df + x^2)^2 +
        lx * (df + 1) * x / (df + x^2)) / px^2
      svv[k] <- (-(df + 2) * (1 / den - 2 * b^2 / den^2) +
        (df + 1) * (df - y^2) / (df + y^2)^2 +
        ly * (df + 1) * y / (df + y^2)) / py^2
    }
  } else if (family == "C") {
    a <- -theta * log(u)
    b <- -theta * log(v)
    m <- pmax(a, b)
    log_s <- m + log(exp(a - m) + exp(b - m) - exp(-m))
    au <- exp(a - log_s)
    bv <- exp(b - log_s)
    su <- (-(theta + 1) + (2 * theta + 1) * au) / u
    sv <- (-(theta + 1) + (2 * theta + 1) * bv) / v
    suv <- theta * (2 * theta + 1) * au * bv / (u * v)
    suu <- (theta + 1 - (2 * theta + 1) * (theta + 1) * au +
      theta * (2 * theta + 1) * au^2) / u^2
    svv <- (theta + 1 - (2 * theta + 1) * (theta + 1) * bv +
      theta * (2 * theta + 1) * bv^2) / v^2
    if (native) {
      # Match the backend density's explicit independence branch.
      independent <- theta <= 1e-10
      su[independent] <- sv[independent] <- suv[independent] <- 0
      suu[independent] <- svv[independent] <- 0
    }
  } else if (family == "F") {
    eu <- exp(-theta * u)
    ev <- exp(-theta * v)
    bu <- expm1(-theta * u)
    bv <- expm1(-theta * v)
    den <- expm1(-theta) + bu * bv
    du <- -theta * eu * bv
    dv <- -theta * ev * bu
    duv <- theta^2 * eu * ev
    su <- -theta - 2 * du / den
    sv <- -theta - 2 * dv / den
    suv <- -2 * (duv / den - (du / den) * (dv / den))
    suu <- 2 * theta * du / den + 2 * (du / den)^2
    svv <- 2 * theta * dv / den + 2 * (dv / den)^2
    independent <- if (native) abs(theta) <= 1e-8 else theta == 0
    su[independent] <- sv[independent] <- suv[independent] <- 0
    suu[independent] <- svv[independent] <- 0
  } else if (family == "G") {
    x <- -log(u)
    y <- -log(v)
    a <- x^theta + y^theta
    s <- a^(1 / theta)
    sx <- s * x^(theta - 1) / a
    sy <- s * y^(theta - 1) / a
    f <- -1 + (1 - 2 * theta) / s + 1 / (s + theta - 1)
    fp <- -(1 - 2 * theta) / s^2 - 1 / (s + theta - 1)^2
    su <- -(1 + (theta - 1) / x + sx * f) / u
    sv <- -(1 + (theta - 1) / y + sy * f) / v
    suv <- sx * sy * ((1 - theta) * f / s + fp) / (u * v)
    sxx <- (theta - 1) * sx / x + (1 - theta) * sx^2 / s
    syy <- (theta - 1) * sy / y + (1 - theta) * sy^2 / s
    suu <- (-(theta - 1) / x^2 + sxx * f + sx^2 * fp - su * u) / u^2
    svv <- (-(theta - 1) / y^2 + syy * f + sy^2 * fp - sv * v) / v^2
  } else if (family == "J") {
    a <- (1 - u)^theta
    b <- (1 - v)^theta
    s <- a + b - a * b
    au <- -theta * a / (1 - u)
    bv <- -theta * b / (1 - v)
    ds_u <- au * (1 - b)
    ds_v <- bv * (1 - a)
    ds_uv <- -au * bv
    # theta*s + (theta-1)*(1-a)*(1-b) = s+theta-1.
    k <- s + theta - 1
    su <- -(theta - 1) / (1 - u) + (1 / theta - 2) * ds_u / s + ds_u / k
    sv <- -(theta - 1) / (1 - v) + (1 / theta - 2) * ds_v / s + ds_v / k
    suv <- (1 / theta - 2) * (ds_uv / s - ds_u * ds_v / s^2) +
      ds_uv / k - ds_u * ds_v / k^2
    ds_uu <- theta * (theta - 1) * a / (1 - u)^2 * (1 - b)
    ds_vv <- theta * (theta - 1) * b / (1 - v)^2 * (1 - a)
    suu <- -(theta - 1) / (1 - u)^2 +
      (1 / theta - 2) * (ds_uu / s - ds_u^2 / s^2) + ds_uu / k - ds_u^2 / k^2
    svv <- -(theta - 1) / (1 - v)^2 +
      (1 / theta - 2) * (ds_vv / s - ds_v^2 / s^2) + ds_vv / k - ds_v^2 / k^2
  } else {
    stop("Unsupported continuous copula mixed derivative: ", family, call. = FALSE)
  }
  if (family %in% c("G", "J")) {
    independent <- theta == 1
    su[independent] <- sv[independent] <- suv[independent] <- 0
    suu[independent] <- svv[independent] <- 0
  }
  # Do not silently convert an unrepresentable derivative into zero here.
  # The caller decides how invalid or unavailable pairs are handled.
  list(density = density, du = density * su, dv = density * sv,
    duu = density * (suu + su^2), dvv = density * (svv + sv^2),
    duv = density * (suv + su * sv))
}

#' @noRd
.copula_density_mixed_arguments <- function(u1, u2, family, par, par2 = 0) {
  .copula_density_argument_curvature(u1, u2, family, par, par2)$duv
}
