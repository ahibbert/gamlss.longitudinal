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


