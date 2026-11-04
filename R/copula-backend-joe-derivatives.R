.copula_joe_deriv <- function(u1, u2, par, deriv, log = FALSE) {
  p <- .copula_joe_parts(u1, u2, par)

  indep <- p$theta <= 1 + 1e-8

  if (any(indep)) {
    out <- numeric(length(p$theta))

    out[indep] <- switch(deriv,
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

  score <- switch(deriv,
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
