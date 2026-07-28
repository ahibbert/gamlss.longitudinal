.copula_v2_clamp01 <- function(x) {
  pmin(pmax(x, 0.001), 0.999)
}

.copula_v2_tau_from_par <- function(family_num, par, par2 = NA_real_) {
  # Handle NA inputs immediately

  if (!is.finite(par)) {
    return(NA_real_)
  }

  tau <- tryCatch(
    {
      if (is.finite(par2)) {
        suppressWarnings(.copula_par_to_tau(family = family_num, par = par, par2 = par2))
      } else {
        suppressWarnings(.copula_par_to_tau(family = family_num, par = par, par2 = 0))
      }
    },
    error = function(e) NA_real_
  )

  if (is.finite(tau)) {
    return(as.numeric(tau))
  }

  # Fallback for Gaussian copulas using formula.

  if (identical(family_num, "N") && is.finite(par)) {
    return(2 / pi * asin(max(min(par, 0.999999), -0.999999)))
  }

  NA_real_
}

.copula_v2_bicop_cdf <- function(u1, u2, family_num, par, par2 = NA_real_) {
  n <- max(length(u1), length(u2), length(par), length(par2))

  if (!is.character(family_num) || length(family_num) != 1L || n < 1) {
    return(rep(NA_real_, length(u1)))
  }

  u1 <- rep(u1, length.out = n)

  u2 <- rep(u2, length.out = n)

  par <- rep(par, length.out = n)

  par2 <- rep(par2, length.out = n)

  vapply(seq_len(n), function(i) {
    if (!is.finite(u1[i]) || !is.finite(u2[i]) || !is.finite(par[i])) {
      return(NA_real_)
    }

    tryCatch(
      {
        .copula_cdf(

          u1[i],
          u2[i],
          family = family_num,
          par = par[i],
          par2 = if (is.finite(par2[i])) par2[i] else 0
        )
      },
      error = function(e) NA_real_
    )
  }, numeric(1), USE.NAMES = FALSE)
}

.copula_v2_bicop_cond_u2_given_u1 <- function(u1, u2, family_num, par, par2 = NA_real_) {
  n <- max(length(u1), length(u2), length(par), length(par2))

  if (!is.character(family_num) || length(family_num) != 1L || n < 1) {
    return(rep(NA_real_, length(u1)))
  }

  u1 <- rep(u1, length.out = n)

  u2 <- rep(u2, length.out = n)

  par <- rep(par, length.out = n)

  par2 <- rep(par2, length.out = n)

  out <- vapply(seq_len(n), function(i) {
    if (!is.finite(u1[i]) || !is.finite(u2[i]) || !is.finite(par[i])) {
      return(NA_real_)
    }

    tryCatch(
      {
        # BiCopHfunc1 gives dC(u1, u2) / du1, i.e. F(U2 <= u2 | U1 = u1).

        .copula_hfunc1(

          u1[i],
          u2[i],
          family = family_num,
          par = par[i],
          par2 = if (is.finite(par2[i])) par2[i] else 0
        )
      },
      error = function(e) NA_real_
    )
  }, numeric(1), USE.NAMES = FALSE)

  .copula_v2_clamp01(as.numeric(out))
}

.copula_v2_message_plot <- function(title, subtitle, message) {
  ggplot2::ggplot(data.frame(x = 0, y = 0), ggplot2::aes(x = x, y = y)) +
    ggplot2::geom_text(label = message, size = 4) +
    ggplot2::xlim(-1, 1) +
    ggplot2::ylim(-1, 1) +
    ggplot2::labs(title = title, subtitle = subtitle) +
    ggplot2::theme_void()
}
