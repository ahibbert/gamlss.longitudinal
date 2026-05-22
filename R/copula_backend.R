.copula_backend <- function() {
  backend <- getOption("gamlss.longitudinal.copula_backend", "native")
  backend <- tolower(as.character(backend)[1])
  if (!backend %in% c("vinecopula", "native")) {
    stop("Unknown copula backend: ", backend, call. = FALSE)
  }
  backend
}

.copula_require_vinecopula <- function(context = "this copula operation") {
  if (!requireNamespace("VineCopula", quietly = TRUE)) {
    stop(
      "VineCopula is required for ", context,
      ". Install VineCopula or use options(gamlss.longitudinal.copula_backend = 'native').",
      call. = FALSE
    )
  }
  invisible(TRUE)
}

.copula_family_code <- function(family) {
  if (!is.character(family) || length(family) != 1L || is.na(family)) {
    stop("Copula family must be a single character code.", call. = FALSE)
  }
  if (!family %in% c("N", "C", "F", "G", "J", "t")) {
    stop(
      "Unsupported copula family code '", family,
      "'. Use one of: N, C, F, G, J, t.",
      call. = FALSE
    )
  }
  family
}

.copula_family_number <- function(family) {
  family <- .copula_family_code(family)
  switch(
    family,
    N = 1,
    C = 3,
    F = 5,
    G = 4,
    J = 6,
    t = 2
  )
}

.copula_family_numbers <- function(family) {
  vapply(family, function(x) {
    if (identical(x, 0) || identical(x, "0")) {
      return(0)
    }
    .copula_family_number(as.character(x))
  }, numeric(1), USE.NAMES = FALSE)
}

.copula_gaussian_rho <- function(par) {
  pmin(pmax(as.numeric(par), -0.999999), 0.999999)
}

.copula_clamp01 <- function(u) {
  pmin(pmax(as.numeric(u), 1e-12), 1 - 1e-12)
}

.copula_recycle <- function(...) {
  args <- list(...)
  n <- max(vapply(args, length, integer(1)))
  lapply(args, rep, length.out = n)
}

.copula_indep_deriv <- function(u1, deriv, log = FALSE) {
  rep(0, length(u1))
}

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

.copula_gaussian_pdf <- function(u1, u2, par) {
  u1 <- .copula_clamp01(u1)
  u2 <- .copula_clamp01(u2)
  rho <- .copula_gaussian_rho(par)
  z1 <- stats::qnorm(u1)
  z2 <- stats::qnorm(u2)
  denom <- 1 - rho^2
  exp((2 * rho * z1 * z2 - rho^2 * (z1^2 + z2^2)) / (2 * denom)) / sqrt(denom)
}

.copula_gaussian_hfunc1 <- function(u1, u2, par) {
  u1 <- .copula_clamp01(u1)
  u2 <- .copula_clamp01(u2)
  rho <- .copula_gaussian_rho(par)
  stats::pnorm((stats::qnorm(u2) - rho * stats::qnorm(u1)) / sqrt(1 - rho^2))
}

.copula_gaussian_deriv <- function(u1, u2, par, deriv, log = FALSE) {
  u1 <- .copula_clamp01(u1)
  u2 <- .copula_clamp01(u2)
  rho <- .copula_gaussian_rho(par)
  n <- max(length(u1), length(u2), length(rho))
  u1 <- rep(u1, length.out = n)
  u2 <- rep(u2, length.out = n)
  rho <- rep(rho, length.out = n)

  z1 <- stats::qnorm(u1)
  z2 <- stats::qnorm(u2)
  denom <- 1 - rho^2
  density <- .copula_gaussian_pdf(u1, u2, rho)

  out <- switch(
    deriv,
    u1 = {
      dlog_du1 <- ((rho * z2 - rho^2 * z1) / denom) / stats::dnorm(z1)
      density * dlog_du1
    },
    u2 = {
      dlog_du2 <- ((rho * z1 - rho^2 * z2) / denom) / stats::dnorm(z2)
      density * dlog_du2
    },
    par = {
      sum_z2 <- z1^2 + z2^2
      numerator <- 2 * rho * z1 * z2 - rho^2 * sum_z2
      numerator_dr <- 2 * z1 * z2 - 2 * rho * sum_z2
      dlog_drho <- rho / denom +
        (numerator_dr * denom + 2 * rho * numerator) / (2 * denom^2)
      if (isTRUE(log)) dlog_drho else density * dlog_drho
    },
    stop("Unsupported Gaussian derivative: ", deriv, call. = FALSE)
  )

  out[!is.finite(out)] <- 0
  out
}

.copula_gaussian_deriv2 <- function(u1, u2, par, deriv) {
  u1 <- .copula_clamp01(u1)
  u2 <- .copula_clamp01(u2)
  rho <- .copula_gaussian_rho(par)
  n <- max(length(u1), length(u2), length(rho))
  u1 <- rep(u1, length.out = n)
  u2 <- rep(u2, length.out = n)
  rho <- rep(rho, length.out = n)

  z1 <- stats::qnorm(u1)
  z2 <- stats::qnorm(u2)
  phi1 <- stats::dnorm(z1)
  phi2 <- stats::dnorm(z2)
  denom <- 1 - rho^2
  density <- .copula_gaussian_pdf(u1, u2, rho)

  out <- switch(
    deriv,
    u1 = {
      dlog_du1 <- ((rho * z2 - rho^2 * z1) / denom) / phi1
      d2log_du1 <- (((rho * z2 - rho^2 * z1) / denom) * z1 -
        rho^2 / denom) / phi1^2
      density * (dlog_du1^2 + d2log_du1)
    },
    u2 = {
      dlog_du2 <- ((rho * z1 - rho^2 * z2) / denom) / phi2
      d2log_du2 <- (((rho * z1 - rho^2 * z2) / denom) * z2 -
        rho^2 / denom) / phi2^2
      density * (dlog_du2^2 + d2log_du2)
    },
    par = {
      sum_z2 <- z1^2 + z2^2
      numerator <- 2 * rho * z1 * z2 - rho^2 * sum_z2
      numerator_dr <- 2 * z1 * z2 - 2 * rho * sum_z2
      numerator_d2r <- -2 * sum_z2
      q <- numerator_dr * denom + 2 * rho * numerator
      q_dr <- numerator_d2r * denom + 2 * numerator
      dlog_drho <- rho / denom + q / (2 * denom^2)
      d2log_drho <- (1 + rho^2) / denom^2 +
        q_dr / (2 * denom^2) +
        2 * rho * q / denom^3
      density * (dlog_drho^2 + d2log_drho)
    },
    stop("Unsupported Gaussian second derivative: ", deriv, call. = FALSE)
  )

  out[!is.finite(out)] <- 0
  out
}

.copula_clayton_theta <- function(par) {
  pmax(as.numeric(par), 0)
}

.copula_clayton_s <- function(u1, u2, theta) {
  u1^(-theta) + u2^(-theta) - 1
}

.copula_clayton_cdf <- function(u1, u2, par) {
  u1 <- .copula_clamp01(u1)
  u2 <- .copula_clamp01(u2)
  theta <- .copula_clayton_theta(par)
  out <- rep(NA_real_, length.out = max(length(u1), length(u2), length(theta)))
  u1 <- rep(u1, length.out = length(out))
  u2 <- rep(u2, length.out = length(out))
  theta <- rep(theta, length.out = length(out))
  out <- u1 * u2
  dep <- theta > 1e-10
  if (any(dep)) {
    s <- .copula_clayton_s(u1[dep], u2[dep], theta[dep])
    out[dep] <- s^(-1 / theta[dep])
  }
  out
}

.copula_clayton_pdf <- function(u1, u2, par) {
  u1 <- .copula_clamp01(u1)
  u2 <- .copula_clamp01(u2)
  theta <- .copula_clayton_theta(par)
  out <- rep(1, length.out = max(length(u1), length(u2), length(theta)))
  u1 <- rep(u1, length.out = length(out))
  u2 <- rep(u2, length.out = length(out))
  theta <- rep(theta, length.out = length(out))
  dep <- theta > 1e-10
  if (any(dep)) {
    s <- .copula_clayton_s(u1[dep], u2[dep], theta[dep])
    out[dep] <- (theta[dep] + 1) *
      (u1[dep] * u2[dep])^(-theta[dep] - 1) *
      s^(-2 - 1 / theta[dep])
  }
  out
}

.copula_clayton_hfunc1 <- function(u1, u2, par) {
  u1 <- .copula_clamp01(u1)
  u2 <- .copula_clamp01(u2)
  theta <- .copula_clayton_theta(par)
  out <- rep(NA_real_, length.out = max(length(u1), length(u2), length(theta)))
  u1 <- rep(u1, length.out = length(out))
  u2 <- rep(u2, length.out = length(out))
  theta <- rep(theta, length.out = length(out))
  indep <- theta <= 1e-10
  out[indep] <- u2[indep]
  dep <- !indep
  if (any(dep)) {
    s <- .copula_clayton_s(u1[dep], u2[dep], theta[dep])
    out[dep] <- u1[dep]^(-theta[dep] - 1) * s^(-1 / theta[dep] - 1)
  }
  .copula_clamp01(out)
}

.copula_clayton_deriv <- function(u1, u2, par, deriv, log = FALSE) {
  u1 <- .copula_clamp01(u1)
  u2 <- .copula_clamp01(u2)
  theta <- .copula_clayton_theta(par)
  n <- max(length(u1), length(u2), length(theta))
  u1 <- rep(u1, length.out = n)
  u2 <- rep(u2, length.out = n)
  theta <- rep(theta, length.out = n)

  out <- numeric(n)
  indep <- theta <= 1e-8
  if (any(indep)) {
    out[indep] <- switch(
      deriv,
      u1 = .copula_indep_deriv(u1[indep], deriv, log),
      u2 = .copula_indep_deriv(u1[indep], deriv, log),
      par = .copula_one_sided_par_deriv(
        .copula_clayton_pdf,
        u1[indep],
        u2[indep],
        par0 = 0,
        h = 1e-5,
        log = log
      ),
      stop("Unsupported Clayton derivative: ", deriv, call. = FALSE)
    )
  }
  if (all(indep)) {
    out[!is.finite(out)] <- 0
    return(out)
  }

  dep <- !indep
  u1_dep <- u1[dep]
  u2_dep <- u2[dep]
  theta_dep <- theta[dep]

  s <- .copula_clayton_s(u1_dep, u2_dep, theta_dep)
  density <- .copula_clayton_pdf(u1_dep, u2_dep, theta_dep)

  out[dep] <- switch(
    deriv,
    u1 = {
      dlog_du1 <- -(theta_dep + 1) / u1_dep + (2 * theta_dep + 1) * u1_dep^(-theta_dep - 1) / s
      density * dlog_du1
    },
    u2 = {
      dlog_du2 <- -(theta_dep + 1) / u2_dep + (2 * theta_dep + 1) * u2_dep^(-theta_dep - 1) / s
      density * dlog_du2
    },
    par = {
      log_u1 <- log(u1_dep)
      log_u2 <- log(u2_dep)
      a <- -2 - 1 / theta_dep
      ds_dtheta <- -log_u1 * u1_dep^(-theta_dep) - log_u2 * u2_dep^(-theta_dep)
      dlog_dtheta <- 1 / (theta_dep + 1) -
        (log_u1 + log_u2) +
        log(s) / theta_dep^2 +
        a * ds_dtheta / s
      if (isTRUE(log)) dlog_dtheta else density * dlog_dtheta
    },
    stop("Unsupported Clayton derivative: ", deriv, call. = FALSE)
  )

  out[!is.finite(out)] <- 0
  out
}

.copula_clayton_deriv2 <- function(u1, u2, par, deriv) {
  u1 <- .copula_clamp01(u1)
  u2 <- .copula_clamp01(u2)
  theta <- .copula_clayton_theta(par)
  n <- max(length(u1), length(u2), length(theta))
  u1 <- rep(u1, length.out = n)
  u2 <- rep(u2, length.out = n)
  theta <- rep(theta, length.out = n)

  out <- numeric(n)
  indep <- theta <= 1e-8
  if (any(indep)) {
    out[indep] <- switch(
      deriv,
      u1 = rep(0, sum(indep)),
      u2 = rep(0, sum(indep)),
      par = .copula_one_sided_par_deriv2(
        .copula_clayton_pdf,
        u1[indep],
        u2[indep],
        par0 = 0,
        h = 1e-4
      ),
      stop("Unsupported Clayton second derivative: ", deriv, call. = FALSE)
    )
  }
  if (all(indep)) {
    out[!is.finite(out)] <- 0
    return(out)
  }

  dep <- !indep
  u1_dep <- u1[dep]
  u2_dep <- u2[dep]
  theta_dep <- theta[dep]
  s <- .copula_clayton_s(u1_dep, u2_dep, theta_dep)
  density <- .copula_clayton_pdf(u1_dep, u2_dep, theta_dep)

  out[dep] <- switch(
    deriv,
    u1 = {
      a1 <- u1_dep^(-theta_dep - 1)
      dlog_du1 <- -(theta_dep + 1) / u1_dep + (2 * theta_dep + 1) * a1 / s
      d2log_du1 <- (theta_dep + 1) / u1_dep^2 +
        (2 * theta_dep + 1) *
          (-(theta_dep + 1) * a1 / (u1_dep * s) + theta_dep * a1^2 / s^2)
      density * (dlog_du1^2 + d2log_du1)
    },
    u2 = {
      a2 <- u2_dep^(-theta_dep - 1)
      dlog_du2 <- -(theta_dep + 1) / u2_dep + (2 * theta_dep + 1) * a2 / s
      d2log_du2 <- (theta_dep + 1) / u2_dep^2 +
        (2 * theta_dep + 1) *
          (-(theta_dep + 1) * a2 / (u2_dep * s) + theta_dep * a2^2 / s^2)
      density * (dlog_du2^2 + d2log_du2)
    },
    par = {
      log_u1 <- log(u1_dep)
      log_u2 <- log(u2_dep)
      a <- -2 - 1 / theta_dep
      s_theta <- -log_u1 * u1_dep^(-theta_dep) - log_u2 * u2_dep^(-theta_dep)
      s_theta2 <- log_u1^2 * u1_dep^(-theta_dep) + log_u2^2 * u2_dep^(-theta_dep)
      dlog_dtheta <- 1 / (theta_dep + 1) -
        (log_u1 + log_u2) +
        log(s) / theta_dep^2 +
        a * s_theta / s
      d2log_dtheta <- -1 / (theta_dep + 1)^2 -
        2 * log(s) / theta_dep^3 +
        2 * s_theta / (theta_dep^2 * s) +
        a * (s_theta2 * s - s_theta^2) / s^2
      density * (dlog_dtheta^2 + d2log_dtheta)
    },
    stop("Unsupported Clayton second derivative: ", deriv, call. = FALSE)
  )

  out[!is.finite(out)] <- 0
  out
}

.copula_frank_par <- function(par) {
  as.numeric(par)
}

.copula_frank_cdf <- function(u1, u2, par) {
  vals <- .copula_recycle(.copula_clamp01(u1), .copula_clamp01(u2), .copula_frank_par(par))
  u1 <- vals[[1]]
  u2 <- vals[[2]]
  theta <- vals[[3]]
  out <- u1 * u2
  dep <- abs(theta) > 1e-8
  if (any(dep)) {
    et <- exp(-theta[dep])
    eu <- exp(-theta[dep] * u1[dep])
    ev <- exp(-theta[dep] * u2[dep])
    out[dep] <- -log1p((eu - 1) * (ev - 1) / (et - 1)) / theta[dep]
  }
  .copula_clamp01(out)
}

.copula_frank_pdf <- function(u1, u2, par) {
  vals <- .copula_recycle(.copula_clamp01(u1), .copula_clamp01(u2), .copula_frank_par(par))
  u1 <- vals[[1]]
  u2 <- vals[[2]]
  theta <- vals[[3]]
  out <- rep(1, length(theta))
  dep <- abs(theta) > 1e-8
  if (any(dep)) {
    et <- exp(-theta[dep])
    eu <- exp(-theta[dep] * u1[dep])
    ev <- exp(-theta[dep] * u2[dep])
    den <- et - 1 + (eu - 1) * (ev - 1)
    out[dep] <- -theta[dep] * (et - 1) * eu * ev / den^2
  }
  out[!is.finite(out)] <- 0
  out
}

.copula_frank_hfunc1 <- function(u1, u2, par) {
  vals <- .copula_recycle(.copula_clamp01(u1), .copula_clamp01(u2), .copula_frank_par(par))
  u1 <- vals[[1]]
  u2 <- vals[[2]]
  theta <- vals[[3]]
  out <- u2
  dep <- abs(theta) > 1e-8
  if (any(dep)) {
    et <- exp(-theta[dep])
    eu <- exp(-theta[dep] * u1[dep])
    ev <- exp(-theta[dep] * u2[dep])
    out[dep] <- eu * (ev - 1) / (et - 1 + (eu - 1) * (ev - 1))
  }
  .copula_clamp01(out)
}

.copula_frank_deriv <- function(u1, u2, par, deriv, log = FALSE) {
  vals <- .copula_recycle(.copula_clamp01(u1), .copula_clamp01(u2), .copula_frank_par(par))
  u1 <- vals[[1]]
  u2 <- vals[[2]]
  theta <- vals[[3]]

  indep <- abs(theta) <= 1e-8
  if (any(indep)) {
    out <- numeric(length(theta))
    out[indep] <- switch(
      deriv,
      u1 = rep(0, sum(indep)),
      u2 = rep(0, sum(indep)),
      par = .copula_central_par_deriv(
        .copula_frank_pdf,
        u1[indep],
        u2[indep],
        par0 = 0,
        h = 1e-5,
        log = log
      ),
      stop("Unsupported Frank derivative: ", deriv, call. = FALSE)
    )
    if (all(indep)) {
      out[!is.finite(out)] <- 0
      return(out)
    }
  } else {
    out <- numeric(length(theta))
  }

  dep <- !indep
  theta_dep <- theta[dep]
  u1_dep <- u1[dep]
  u2_dep <- u2[dep]
  et <- exp(-theta_dep)
  eu <- exp(-theta_dep * u1_dep)
  ev <- exp(-theta_dep * u2_dep)
  a <- et - 1
  den <- a + (eu - 1) * (ev - 1)
  density <- .copula_frank_pdf(u1_dep, u2_dep, theta_dep)

  score <- switch(
    deriv,
    u1 = -theta_dep + 2 * theta_dep * eu * (ev - 1) / den,
    u2 = -theta_dep + 2 * theta_dep * ev * (eu - 1) / den,
    par = {
      da <- -et
      dden <- da - u1_dep * eu * (ev - 1) - u2_dep * ev * (eu - 1)
      1 / theta_dep + da / a - (u1_dep + u2_dep) - 2 * dden / den
    },
    stop("Unsupported Frank derivative: ", deriv, call. = FALSE)
  )

  out[dep] <- if (isTRUE(log)) score else density * score
  out[!is.finite(out)] <- 0
  out
}

.copula_frank_deriv2 <- function(u1, u2, par, deriv) {
  vals <- .copula_recycle(.copula_clamp01(u1), .copula_clamp01(u2), .copula_frank_par(par))
  u1 <- vals[[1]]
  u2 <- vals[[2]]
  theta <- vals[[3]]

  indep <- abs(theta) <= 1e-6
  if (any(indep)) {
    out <- numeric(length(theta))
    out[indep] <- switch(
      deriv,
      u1 = rep(0, sum(indep)),
      u2 = rep(0, sum(indep)),
      par = .copula_central_par_deriv2(
        .copula_frank_pdf,
        u1[indep],
        u2[indep],
        par0 = 0,
        h = 1e-4
      ),
      stop("Unsupported Frank second derivative: ", deriv, call. = FALSE)
    )
    if (all(indep)) {
      out[!is.finite(out)] <- 0
      return(out)
    }
  } else {
    out <- numeric(length(theta))
  }

  dep <- !indep
  u1_dep <- u1[dep]
  u2_dep <- u2[dep]
  theta_dep <- theta[dep]

  out[dep] <- switch(
    deriv,
    u1 = {
      h <- pmin(1e-5, 0.25 * u1_dep, 0.25 * (1 - u1_dep))
      (
        .copula_frank_deriv(u1_dep + h, u2_dep, theta_dep, deriv = "u1") -
          .copula_frank_deriv(u1_dep - h, u2_dep, theta_dep, deriv = "u1")
      ) / (2 * h)
    },
    u2 = {
      h <- pmin(1e-5, 0.25 * u2_dep, 0.25 * (1 - u2_dep))
      (
        .copula_frank_deriv(u1_dep, u2_dep + h, theta_dep, deriv = "u2") -
          .copula_frank_deriv(u1_dep, u2_dep - h, theta_dep, deriv = "u2")
      ) / (2 * h)
    },
    par = {
      h <- pmin(1e-4, 0.25 * abs(theta_dep))
      (
        .copula_frank_deriv(u1_dep, u2_dep, theta_dep + h, deriv = "par") -
          .copula_frank_deriv(u1_dep, u2_dep, theta_dep - h, deriv = "par")
      ) / (2 * h)
    },
    stop("Unsupported Frank second derivative: ", deriv, call. = FALSE)
  )

  out[!is.finite(out)] <- 0
  out
}

.copula_frank_par_to_tau <- function(par) {
  theta <- .copula_frank_par(par)
  vapply(theta, function(th) {
    if (abs(th) <= 1e-8) return(0)
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
    if (abs(tau_i) <= 1e-8) return(0)
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

.copula_gumbel_theta <- function(par) {
  pmax(as.numeric(par), 1)
}

.copula_gumbel_parts <- function(u1, u2, par) {
  vals <- .copula_recycle(.copula_clamp01(u1), .copula_clamp01(u2), .copula_gumbel_theta(par))
  u1 <- vals[[1]]
  u2 <- vals[[2]]
  theta <- vals[[3]]
  x <- -log(u1)
  y <- -log(u2)
  a <- x^theta + y^theta
  s <- a^(1 / theta)
  list(u1 = u1, u2 = u2, theta = theta, x = x, y = y, a = a, s = s, cdf = exp(-s))
}

.copula_gumbel_cdf <- function(u1, u2, par) {
  .copula_gumbel_parts(u1, u2, par)$cdf
}

.copula_gumbel_pdf <- function(u1, u2, par) {
  p <- .copula_gumbel_parts(u1, u2, par)
  out <- p$cdf * (p$x * p$y)^(p$theta - 1) *
    p$s^(1 - 2 * p$theta) * (p$s + p$theta - 1) / (p$u1 * p$u2)
  out[!is.finite(out)] <- 0
  out
}

.copula_gumbel_hfunc1 <- function(u1, u2, par) {
  p <- .copula_gumbel_parts(u1, u2, par)
  out <- p$cdf * p$s^(1 - p$theta) * p$x^(p$theta - 1) / p$u1
  .copula_clamp01(out)
}

.copula_gumbel_deriv <- function(u1, u2, par, deriv, log = FALSE) {
  p <- .copula_gumbel_parts(u1, u2, par)

  indep <- p$theta <= 1 + 1e-8
  if (any(indep)) {
    out <- numeric(length(p$theta))
    out[indep] <- switch(
      deriv,
      u1 = rep(0, sum(indep)),
      u2 = rep(0, sum(indep)),
      par = .copula_one_sided_par_deriv(
        .copula_gumbel_pdf,
        p$u1[indep],
        p$u2[indep],
        par0 = 1,
        h = 1e-5,
        log = log
      ),
      stop("Unsupported Gumbel derivative: ", deriv, call. = FALSE)
    )
    if (all(indep)) {
      out[!is.finite(out)] <- 0
      return(out)
    }
  } else {
    out <- numeric(length(p$theta))
  }

  dep <- !indep
  p_dep <- lapply(p, function(x) x[dep])
  density <- .copula_gumbel_pdf(p_dep$u1, p_dep$u2, p_dep$theta)
  score <- switch(
    deriv,
    u1 = {
      ds_dx <- p_dep$s * p_dep$x^(p_dep$theta - 1) / p_dep$a
      dlog_dx <- 1 + (p_dep$theta - 1) / p_dep$x +
        ds_dx * (-1 + (1 - 2 * p_dep$theta) / p_dep$s + 1 / (p_dep$s + p_dep$theta - 1))
      -dlog_dx / p_dep$u1
    },
    u2 = {
      ds_dy <- p_dep$s * p_dep$y^(p_dep$theta - 1) / p_dep$a
      dlog_dy <- 1 + (p_dep$theta - 1) / p_dep$y +
        ds_dy * (-1 + (1 - 2 * p_dep$theta) / p_dep$s + 1 / (p_dep$s + p_dep$theta - 1))
      -dlog_dy / p_dep$u2
    },
    par = {
      log_x <- log(p_dep$x)
      log_y <- log(p_dep$y)
      da_dtheta <- p_dep$x^p_dep$theta * log_x + p_dep$y^p_dep$theta * log_y
      ds_dtheta <- p_dep$s * (da_dtheta / (p_dep$theta * p_dep$a) - log(p_dep$a) / p_dep$theta^2)
      -ds_dtheta +
        log_x + log_y -
        2 * log(p_dep$s) +
        (1 - 2 * p_dep$theta) * ds_dtheta / p_dep$s +
        (ds_dtheta + 1) / (p_dep$s + p_dep$theta - 1)
    },
    stop("Unsupported Gumbel derivative: ", deriv, call. = FALSE)
  )

  out[dep] <- if (isTRUE(log)) score else density * score
  out[!is.finite(out)] <- 0
  out
}

.copula_gumbel_deriv2 <- function(u1, u2, par, deriv) {
  p <- .copula_gumbel_parts(u1, u2, par)

  indep <- p$theta <= 1 + 1e-6
  if (any(indep)) {
    out <- numeric(length(p$theta))
    out[indep] <- switch(
      deriv,
      u1 = rep(0, sum(indep)),
      u2 = rep(0, sum(indep)),
      par = .copula_one_sided_par_deriv2(
        .copula_gumbel_pdf,
        p$u1[indep],
        p$u2[indep],
        par0 = 1,
        h = 1e-4
      ),
      stop("Unsupported Gumbel second derivative: ", deriv, call. = FALSE)
    )
    if (all(indep)) {
      out[!is.finite(out)] <- 0
      return(out)
    }
  } else {
    out <- numeric(length(p$theta))
  }

  dep <- !indep
  p_dep <- lapply(p, function(x) x[dep])
  out[dep] <- switch(
    deriv,
    u1 = {
      h <- pmin(1e-5, 0.25 * p_dep$u1, 0.25 * (1 - p_dep$u1))
      up <- p_dep$u1 + h
      um <- p_dep$u1 - h
      (
        .copula_gumbel_deriv(up, p_dep$u2, p_dep$theta, deriv = "u1") -
          .copula_gumbel_deriv(um, p_dep$u2, p_dep$theta, deriv = "u1")
      ) / (2 * h)
    },
    u2 = {
      h <- pmin(1e-5, 0.25 * p_dep$u2, 0.25 * (1 - p_dep$u2))
      up <- p_dep$u2 + h
      um <- p_dep$u2 - h
      (
        .copula_gumbel_deriv(p_dep$u1, up, p_dep$theta, deriv = "u2") -
          .copula_gumbel_deriv(p_dep$u1, um, p_dep$theta, deriv = "u2")
      ) / (2 * h)
    },
    par = {
      h <- pmin(1e-4, 0.25 * (p_dep$theta - 1))
      tp <- p_dep$theta + h
      tm <- p_dep$theta - h
      (
        .copula_gumbel_deriv(p_dep$u1, p_dep$u2, tp, deriv = "par") -
          .copula_gumbel_deriv(p_dep$u1, p_dep$u2, tm, deriv = "par")
      ) / (2 * h)
    },
    stop("Unsupported Gumbel second derivative: ", deriv, call. = FALSE)
  )

  out[!is.finite(out)] <- 0
  out
}

.copula_joe_theta <- function(par) {
  pmax(as.numeric(par), 1)
}

.copula_joe_parts <- function(u1, u2, par) {
  vals <- .copula_recycle(.copula_clamp01(u1), .copula_clamp01(u2), .copula_joe_theta(par))
  u1 <- vals[[1]]
  u2 <- vals[[2]]
  theta <- vals[[3]]
  a <- (1 - u1)^theta
  b <- (1 - u2)^theta
  s <- a + b - a * b
  list(u1 = u1, u2 = u2, theta = theta, a = a, b = b, s = s, cdf = 1 - s^(1 / theta))
}

.copula_joe_cdf <- function(u1, u2, par) {
  .copula_clamp01(.copula_joe_parts(u1, u2, par)$cdf)
}

.copula_joe_pdf <- function(u1, u2, par) {
  p <- .copula_joe_parts(u1, u2, par)
  au <- (1 - p$u1)^(p$theta - 1)
  bu <- (1 - p$u2)^(p$theta - 1)
  du <- 1 - p$a
  dv <- 1 - p$b
  m <- 1 / p$theta - 1
  out <- au * bu * p$s^(m - 1) * (p$theta * p$s + (p$theta - 1) * du * dv)
  out[!is.finite(out)] <- 0
  out
}

.copula_joe_hfunc1 <- function(u1, u2, par) {
  p <- .copula_joe_parts(u1, u2, par)
  out <- (1 - p$u1)^(p$theta - 1) * (1 - p$b) * p$s^(1 / p$theta - 1)
  .copula_clamp01(out)
}

.copula_joe_deriv <- function(u1, u2, par, deriv, log = FALSE) {
  p <- .copula_joe_parts(u1, u2, par)

  indep <- p$theta <= 1 + 1e-8
  if (any(indep)) {
    out <- numeric(length(p$theta))
    out[indep] <- switch(
      deriv,
      u1 = rep(0, sum(indep)),
      u2 = rep(0, sum(indep)),
      par = .copula_one_sided_par_deriv(
        .copula_joe_pdf,
        p$u1[indep],
        p$u2[indep],
        par0 = 1,
        h = 1e-5,
        log = log
      ),
      stop("Unsupported Joe derivative: ", deriv, call. = FALSE)
    )
    if (all(indep)) {
      out[!is.finite(out)] <- 0
      return(out)
    }
  } else {
    out <- numeric(length(p$theta))
  }

  dep <- !indep
  p_dep <- lapply(p, function(x) x[dep])
  r <- 1 - p_dep$u1
  q <- 1 - p_dep$u2
  log_r <- log(r)
  log_q <- log(q)
  du <- 1 - p_dep$a
  dv <- 1 - p_dep$b
  k <- p_dep$theta * p_dep$s + (p_dep$theta - 1) * du * dv
  density <- .copula_joe_pdf(p_dep$u1, p_dep$u2, p_dep$theta)

  score <- switch(
    deriv,
    u1 = {
      da_du <- -p_dep$theta * r^(p_dep$theta - 1)
      ds_du <- da_du * dv
      du_du <- -da_du
      dk_du <- p_dep$theta * ds_du + (p_dep$theta - 1) * du_du * dv
      -(p_dep$theta - 1) / r + (1 / p_dep$theta - 2) * ds_du / p_dep$s + dk_du / k
    },
    u2 = {
      db_du <- -p_dep$theta * q^(p_dep$theta - 1)
      ds_du <- db_du * du
      dv_du <- -db_du
      dk_du <- p_dep$theta * ds_du + (p_dep$theta - 1) * du * dv_du
      -(p_dep$theta - 1) / q + (1 / p_dep$theta - 2) * ds_du / p_dep$s + dk_du / k
    },
    par = {
      da_dt <- p_dep$a * log_r
      db_dt <- p_dep$b * log_q
      ds_dt <- da_dt * dv + db_dt * du
      du_dt <- -da_dt
      dv_dt <- -db_dt
      dk_dt <- p_dep$s + p_dep$theta * ds_dt + du * dv +
        (p_dep$theta - 1) * (du_dt * dv + du * dv_dt)
      log_r + log_q - log(p_dep$s) / p_dep$theta^2 +
        (1 / p_dep$theta - 2) * ds_dt / p_dep$s + dk_dt / k
    },
    stop("Unsupported Joe derivative: ", deriv, call. = FALSE)
  )

  out[dep] <- if (isTRUE(log)) score else density * score
  out[!is.finite(out)] <- 0
  out
}

.copula_joe_deriv2 <- function(u1, u2, par, deriv) {
  p <- .copula_joe_parts(u1, u2, par)

  indep <- p$theta <= 1 + 1e-6
  if (any(indep)) {
    out <- numeric(length(p$theta))
    out[indep] <- switch(
      deriv,
      u1 = rep(0, sum(indep)),
      u2 = rep(0, sum(indep)),
      par = .copula_one_sided_par_deriv2(
        .copula_joe_pdf,
        p$u1[indep],
        p$u2[indep],
        par0 = 1,
        h = 1e-4
      ),
      stop("Unsupported Joe second derivative: ", deriv, call. = FALSE)
    )
    if (all(indep)) {
      out[!is.finite(out)] <- 0
      return(out)
    }
  } else {
    out <- numeric(length(p$theta))
  }

  dep <- !indep
  p_dep <- lapply(p, function(x) x[dep])
  out[dep] <- switch(
    deriv,
    u1 = {
      h <- pmin(1e-5, 0.25 * p_dep$u1, 0.25 * (1 - p_dep$u1))
      (
        .copula_joe_deriv(p_dep$u1 + h, p_dep$u2, p_dep$theta, deriv = "u1") -
          .copula_joe_deriv(p_dep$u1 - h, p_dep$u2, p_dep$theta, deriv = "u1")
      ) / (2 * h)
    },
    u2 = {
      h <- pmin(1e-5, 0.25 * p_dep$u2, 0.25 * (1 - p_dep$u2))
      (
        .copula_joe_deriv(p_dep$u1, p_dep$u2 + h, p_dep$theta, deriv = "u2") -
          .copula_joe_deriv(p_dep$u1, p_dep$u2 - h, p_dep$theta, deriv = "u2")
      ) / (2 * h)
    },
    par = {
      h <- pmin(1e-4, 0.25 * (p_dep$theta - 1))
      (
        .copula_joe_deriv(p_dep$u1, p_dep$u2, p_dep$theta + h, deriv = "par") -
          .copula_joe_deriv(p_dep$u1, p_dep$u2, p_dep$theta - h, deriv = "par")
      ) / (2 * h)
    },
    stop("Unsupported Joe second derivative: ", deriv, call. = FALSE)
  )

  out[!is.finite(out)] <- 0
  out
}

.copula_joe_par_to_tau <- function(par) {
  theta <- .copula_joe_theta(par)
  vapply(theta, function(th) {
    if (th <= 1 + 1e-8) return(0)
    phi <- function(u) -log1p(-(1 - u)^th)
    dphi <- function(u) -th * (1 - u)^(th - 1) / (1 - (1 - u)^th)
    integrand <- function(u) phi(u) / dphi(u)
    1 + 4 * stats::integrate(integrand, lower = 0, upper = 1, rel.tol = 1e-8)$value
  }, numeric(1), USE.NAMES = FALSE)
}

.copula_joe_tau_to_par <- function(tau) {
  tau <- pmin(pmax(as.numeric(tau), 0), 0.999999)
  vapply(tau, function(tau_i) {
    if (tau_i <= 1e-8) return(1)
    stats::uniroot(
      function(th) .copula_joe_par_to_tau(th) - tau_i,
      lower = 1 + 1e-8,
      upper = 100,
      tol = 1e-8
    )$root
  }, numeric(1), USE.NAMES = FALSE)
}

.copula_t_df <- function(par2) {
  pmax(as.numeric(par2), 2.000001)
}

.copula_t_df_step <- function(df, rel = 1e-3) {
  pmin(pmax(rel * abs(df), 1e-4), 0.25 * (df - 2))
}

.copula_t_pdf <- function(u1, u2, par, par2) {
  vals <- .copula_recycle(.copula_clamp01(u1), .copula_clamp01(u2), .copula_gaussian_rho(par), .copula_t_df(par2))
  u1 <- vals[[1]]
  u2 <- vals[[2]]
  rho <- vals[[3]]
  df <- vals[[4]]
  x <- stats::qt(u1, df = df)
  y <- stats::qt(u2, df = df)
  q <- (x^2 - 2 * rho * x * y + y^2) / (1 - rho^2)
  log_biv <- lgamma((df + 2) / 2) - lgamma(df / 2) -
    log(df * pi) - 0.5 * log1p(-rho^2) -
    (df + 2) / 2 * log1p(q / df)
  out <- exp(log_biv - stats::dt(x, df = df, log = TRUE) - stats::dt(y, df = df, log = TRUE))
  out[!is.finite(out)] <- 0
  out
}

.copula_t_hfunc1 <- function(u1, u2, par, par2) {
  vals <- .copula_recycle(.copula_clamp01(u1), .copula_clamp01(u2), .copula_gaussian_rho(par), .copula_t_df(par2))
  u1 <- vals[[1]]
  u2 <- vals[[2]]
  rho <- vals[[3]]
  df <- vals[[4]]
  x <- stats::qt(u1, df = df)
  y <- stats::qt(u2, df = df)
  scale <- sqrt((df + x^2) * (1 - rho^2) / (df + 1))
  .copula_clamp01(stats::pt((y - rho * x) / scale, df = df + 1))
}

.copula_t_cdf_one <- function(u1, u2, rho, df) {
  x <- stats::qt(u1, df = df)
  y <- stats::qt(u2, df = df)
  integrand <- function(z) {
    scale <- sqrt((df + z^2) * (1 - rho^2) / (df + 1))
    stats::dt(z, df = df) * stats::pt((y - rho * z) / scale, df = df + 1)
  }
  stats::integrate(integrand, lower = -Inf, upper = x, rel.tol = 1e-7)$value
}

.copula_t_cdf <- function(u1, u2, par, par2) {
  vals <- .copula_recycle(.copula_clamp01(u1), .copula_clamp01(u2), .copula_gaussian_rho(par), .copula_t_df(par2))
  out <- mapply(
    .copula_t_cdf_one,
    vals[[1]],
    vals[[2]],
    vals[[3]],
    vals[[4]],
    SIMPLIFY = TRUE,
    USE.NAMES = FALSE
  )
  .copula_clamp01(out)
}

.copula_t_deriv <- function(u1, u2, par, par2, deriv, log = FALSE) {
  vals <- .copula_recycle(.copula_clamp01(u1), .copula_clamp01(u2), .copula_gaussian_rho(par), .copula_t_df(par2))
  u1 <- vals[[1]]
  u2 <- vals[[2]]
  rho <- vals[[3]]
  df <- vals[[4]]
  x <- stats::qt(u1, df = df)
  y <- stats::qt(u2, df = df)
  rho_denom <- 1 - rho^2
  numerator <- x^2 - 2 * rho * x * y + y^2
  denom <- df * rho_denom + numerator
  q <- numerator / rho_denom
  log_biv <- lgamma((df + 2) / 2) - lgamma(df / 2) -
    log(df * pi) - 0.5 * log1p(-rho^2) -
    (df + 2) / 2 * log1p(q / df)
  log_density <- log_biv - stats::dt(x, df = df, log = TRUE) - stats::dt(y, df = df, log = TRUE)
  density <- exp(log_density)
  density[!is.finite(density)] <- 0

  score <- switch(
    deriv,
    u1 = {
      dlog_dx <- -(df + 2) * (x - rho * y) / denom +
        (df + 1) * x / (df + x^2)
      dlog_dx / stats::dt(x, df = df)
    },
    u2 = {
      dlog_dy <- -(df + 2) * (y - rho * x) / denom +
        (df + 1) * y / (df + y^2)
      dlog_dy / stats::dt(y, df = df)
    },
    par = {
      dq_drho <- (-2 * x * y * rho_denom + 2 * rho * numerator) / rho_denom^2
      rho / rho_denom - (df + 2) * dq_drho / (2 * (df + q))
    },
    par2 = {
      h <- .copula_t_df_step(df)
      dcd_df <- (
        .copula_t_pdf(u1, u2, rho, df + h) -
          .copula_t_pdf(u1, u2, rho, df - h)
      ) / (2 * h)
      dcd_df / density
    },
    stop("Unsupported t derivative: ", deriv, call. = FALSE)
  )

  out <- if (isTRUE(log)) score else density * score
  out[!is.finite(out)] <- 0
  out
}

.copula_t_deriv2 <- function(u1, u2, par, par2, deriv) {
  vals <- .copula_recycle(.copula_clamp01(u1), .copula_clamp01(u2), .copula_gaussian_rho(par), .copula_t_df(par2))
  u1 <- vals[[1]]
  u2 <- vals[[2]]
  rho <- vals[[3]]
  df <- vals[[4]]

  out <- switch(
    deriv,
    u1 = {
      h <- pmin(1e-5, 0.25 * u1, 0.25 * (1 - u1))
      (
        .copula_t_deriv(u1 + h, u2, rho, df, deriv = "u1") -
          .copula_t_deriv(u1 - h, u2, rho, df, deriv = "u1")
      ) / (2 * h)
    },
    u2 = {
      h <- pmin(1e-5, 0.25 * u2, 0.25 * (1 - u2))
      (
        .copula_t_deriv(u1, u2 + h, rho, df, deriv = "u2") -
          .copula_t_deriv(u1, u2 - h, rho, df, deriv = "u2")
      ) / (2 * h)
    },
    par = {
      h <- pmin(1e-4, 0.25 * (1 - abs(rho)))
      (
        .copula_t_deriv(u1, u2, rho + h, df, deriv = "par") -
          .copula_t_deriv(u1, u2, rho - h, df, deriv = "par")
      ) / (2 * h)
    },
    par2 = {
      h <- .copula_t_df_step(df)
      (
        .copula_t_pdf(u1, u2, rho, df + h) -
          2 * .copula_t_pdf(u1, u2, rho, df) +
          .copula_t_pdf(u1, u2, rho, df - h)
      ) / (h^2)
    },
    par1par2 = {
      h <- .copula_t_df_step(df)
      (
        .copula_t_deriv(u1, u2, rho, df + h, deriv = "par") -
          .copula_t_deriv(u1, u2, rho, df - h, deriv = "par")
      ) / (2 * h)
    },
    stop("Unsupported t second derivative: ", deriv, call. = FALSE)
  )

  out[!is.finite(out)] <- 0
  out
}

.copula_pdf <- function(u1, u2, family, par, par2 = 0) {
  family <- .copula_family_code(family)
  if (.copula_backend() == "native" && family == "N") {
    return(.copula_gaussian_pdf(u1, u2, par))
  }
  if (.copula_backend() == "native" && family == "C") {
    return(.copula_clayton_pdf(u1, u2, par))
  }
  if (.copula_backend() == "native" && family == "F") {
    return(.copula_frank_pdf(u1, u2, par))
  }
  if (.copula_backend() == "native" && family == "G") {
    return(.copula_gumbel_pdf(u1, u2, par))
  }
  if (.copula_backend() == "native" && family == "J") {
    return(.copula_joe_pdf(u1, u2, par))
  }
  if (.copula_backend() == "native" && family == "t") {
    return(.copula_t_pdf(u1, u2, par, par2))
  }
  .copula_require_vinecopula("the delegated VineCopula PDF backend")
  VineCopula::BiCopPDF(
    u1,
    u2,
    family = .copula_family_number(family),
    par = par,
    par2 = par2
  )
}

.copula_deriv <- function(u1, u2, family, par, par2 = 0, deriv, log = FALSE) {
  family <- .copula_family_code(family)
  if (.copula_backend() == "native" && family == "N" && deriv %in% c("u1", "u2", "par")) {
    return(.copula_gaussian_deriv(u1, u2, par, deriv = deriv, log = log))
  }
  if (.copula_backend() == "native" && family == "C" && deriv %in% c("u1", "u2", "par")) {
    return(.copula_clayton_deriv(u1, u2, par, deriv = deriv, log = log))
  }
  if (.copula_backend() == "native" && family == "G" && deriv %in% c("u1", "u2", "par")) {
    return(.copula_gumbel_deriv(u1, u2, par, deriv = deriv, log = log))
  }
  if (.copula_backend() == "native" && family == "F" && deriv %in% c("u1", "u2", "par")) {
    return(.copula_frank_deriv(u1, u2, par, deriv = deriv, log = log))
  }
  if (.copula_backend() == "native" && family == "J" && deriv %in% c("u1", "u2", "par")) {
    return(.copula_joe_deriv(u1, u2, par, deriv = deriv, log = log))
  }
  if (.copula_backend() == "native" && family == "t" && deriv %in% c("u1", "u2", "par", "par2")) {
    return(.copula_t_deriv(u1, u2, par, par2, deriv = deriv, log = log))
  }
  .copula_require_vinecopula("the delegated VineCopula first-derivative backend")
  VineCopula::BiCopDeriv(
    u1,
    u2,
    family = .copula_family_number(family),
    par = par,
    par2 = par2,
    deriv = deriv,
    log = log
  )
}

.copula_deriv2 <- function(u1, u2, family, par, par2 = 0, deriv) {
  family <- .copula_family_code(family)
  if (.copula_backend() == "native" && family == "N" && deriv %in% c("u1", "u2", "par")) {
    return(.copula_gaussian_deriv2(u1, u2, par, deriv = deriv))
  }
  if (.copula_backend() == "native" && family == "C" && deriv %in% c("u1", "u2", "par")) {
    return(.copula_clayton_deriv2(u1, u2, par, deriv = deriv))
  }
  if (.copula_backend() == "native" && family == "G" && deriv %in% c("u1", "u2", "par")) {
    return(.copula_gumbel_deriv2(u1, u2, par, deriv = deriv))
  }
  if (.copula_backend() == "native" && family == "F" && deriv %in% c("u1", "u2", "par")) {
    return(.copula_frank_deriv2(u1, u2, par, deriv = deriv))
  }
  if (.copula_backend() == "native" && family == "J" && deriv %in% c("u1", "u2", "par")) {
    return(.copula_joe_deriv2(u1, u2, par, deriv = deriv))
  }
  if (.copula_backend() == "native" && family == "t" && deriv %in% c("u1", "u2", "par", "par2", "par1par2")) {
    return(.copula_t_deriv2(u1, u2, par, par2, deriv = deriv))
  }
  .copula_require_vinecopula("the delegated VineCopula second-derivative backend")
  VineCopula::BiCopDeriv2(
    u1,
    u2,
    family = .copula_family_number(family),
    par = par,
    par2 = par2,
    deriv = deriv
  )
}

.copula_tau_to_par <- function(family, tau) {
  family <- .copula_family_code(family)
  if (.copula_backend() == "native" && family == "N") {
    return(sin(pi * tau / 2))
  }
  if (.copula_backend() == "native" && family == "C") {
    tau <- pmin(pmax(as.numeric(tau), 0), 1 - 1e-12)
    return(2 * tau / (1 - tau))
  }
  if (.copula_backend() == "native" && family == "F") {
    return(.copula_frank_tau_to_par(tau))
  }
  if (.copula_backend() == "native" && family == "G") {
    tau <- pmin(pmax(as.numeric(tau), 0), 1 - 1e-12)
    return(1 / (1 - tau))
  }
  if (.copula_backend() == "native" && family == "J") {
    return(.copula_joe_tau_to_par(tau))
  }
  if (.copula_backend() == "native" && family == "t") {
    return(sin(pi * tau / 2))
  }
  .copula_require_vinecopula("the delegated VineCopula tau-to-parameter backend")
  VineCopula::BiCopTau2Par(
    family = .copula_family_number(family),
    tau = tau
  )
}

.copula_par_to_tau <- function(family, par, par2 = 0) {
  family <- .copula_family_code(family)
  if (.copula_backend() == "native" && family == "N") {
    return(2 / pi * asin(.copula_gaussian_rho(par)))
  }
  if (.copula_backend() == "native" && family == "C") {
    theta <- .copula_clayton_theta(par)
    return(theta / (theta + 2))
  }
  if (.copula_backend() == "native" && family == "F") {
    return(.copula_frank_par_to_tau(par))
  }
  if (.copula_backend() == "native" && family == "G") {
    theta <- .copula_gumbel_theta(par)
    return(1 - 1 / theta)
  }
  if (.copula_backend() == "native" && family == "J") {
    return(.copula_joe_par_to_tau(par))
  }
  if (.copula_backend() == "native" && family == "t") {
    return(2 / pi * asin(.copula_gaussian_rho(par)))
  }
  .copula_require_vinecopula("the delegated VineCopula parameter-to-tau backend")
  VineCopula::BiCopPar2Tau(
    family = .copula_family_number(family),
    par = par,
    par2 = par2
  )
}

.copula_cdf <- function(u1, u2, family, par, par2 = 0) {
  family <- .copula_family_code(family)
  if (.copula_backend() == "native" && family == "N" && all(abs(.copula_gaussian_rho(par)) <= 1e-12)) {
    vals <- .copula_recycle(.copula_clamp01(u1), .copula_clamp01(u2), par)
    return(vals[[1]] * vals[[2]])
  }
  if (.copula_backend() == "native" && family == "C") {
    return(.copula_clayton_cdf(u1, u2, par))
  }
  if (.copula_backend() == "native" && family == "F") {
    return(.copula_frank_cdf(u1, u2, par))
  }
  if (.copula_backend() == "native" && family == "G") {
    return(.copula_gumbel_cdf(u1, u2, par))
  }
  if (.copula_backend() == "native" && family == "J") {
    return(.copula_joe_cdf(u1, u2, par))
  }
  if (.copula_backend() == "native" && family == "t") {
    return(.copula_t_cdf(u1, u2, par, par2))
  }
  .copula_require_vinecopula("the delegated VineCopula CDF backend")
  VineCopula::BiCopCDF(
    u1,
    u2,
    family = .copula_family_number(family),
    par = par,
    par2 = par2
  )
}

.copula_hfunc1 <- function(u1, u2, family, par, par2 = 0) {
  family <- .copula_family_code(family)
  if (.copula_backend() == "native" && family == "N") {
    return(.copula_gaussian_hfunc1(u1, u2, par))
  }
  if (.copula_backend() == "native" && family == "C") {
    return(.copula_clayton_hfunc1(u1, u2, par))
  }
  if (.copula_backend() == "native" && family == "F") {
    return(.copula_frank_hfunc1(u1, u2, par))
  }
  if (.copula_backend() == "native" && family == "G") {
    return(.copula_gumbel_hfunc1(u1, u2, par))
  }
  if (.copula_backend() == "native" && family == "J") {
    return(.copula_joe_hfunc1(u1, u2, par))
  }
  if (.copula_backend() == "native" && family == "t") {
    return(.copula_t_hfunc1(u1, u2, par, par2))
  }
  .copula_require_vinecopula("the delegated VineCopula h-function backend")
  VineCopula::BiCopHfunc1(
    u1,
    u2,
    family = .copula_family_number(family),
    par = par,
    par2 = par2
  )
}

.copula_dvine <- function(order, family, par, par2 = 0) {
  .copula_require_vinecopula("D-vine simulation setup via .copula_dvine()")
  VineCopula::D2RVine(order, .copula_family_numbers(family), par, par2)
}

.copula_rvine_sim <- function(n, rvm) {
  .copula_require_vinecopula("R-vine simulation via .copula_rvine_sim()")
  VineCopula::RVineSim(n, rvm)
}
