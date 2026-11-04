.copula_gumbel_deriv <- function(u1, u2, par, deriv, log = FALSE) {
  p <- .copula_gumbel_parts(u1, u2, par)

  indep <- p$theta <= 1 + 1e-8

  if (any(indep)) {
    out <- numeric(length(p$theta))

    out[indep] <- switch(deriv,
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

  score <- switch(deriv,
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
