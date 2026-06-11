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


