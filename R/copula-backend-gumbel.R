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


