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


