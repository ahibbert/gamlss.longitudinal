#' Finite-difference derivatives of marginal CDF bounds
#'
#' Computes observation-level derivatives of CDF bounds with respect to natural
#' margin parameters. These feed the discrete rectangle likelihood derivatives.
#'
#' @noRd
.calc_zibnb_F_bounds_derivatives <- function(
    eta_inv, response, par_names = c("mu", "sigma", "nu", "tau")) {
  required <- c("mu", "sigma", "nu", "tau")
  if (!all(required %in% names(eta_inv))) return(NULL)

  n_obs <- length(response)
  mu <- rep(as.numeric(eta_inv$mu), length.out = n_obs)
  sigma <- rep(as.numeric(eta_inv$sigma), length.out = n_obs)
  nu <- rep(as.numeric(eta_inv$nu), length.out = n_obs)
  tau <- rep(as.numeric(eta_inv$tau), length.out = n_obs)
  observed <- is.finite(response)
  valid <- observed & response >= 0 & response == floor(response) &
    is.finite(mu) & mu > 0 & is.finite(sigma) & sigma > 0 &
    is.finite(nu) & nu > 0 & is.finite(tau) & tau > 0 & tau < 1
  if (any(observed & !valid)) return(NULL)

  q <- integer(n_obs)
  q[valid] <- as.integer(response[valid])
  m <- 1 / sigma + 1
  shape_n <- mu * nu / sigma
  shape_k <- 1 / nu

  probability <- numeric(n_obs)
  score_n <- numeric(n_obs)
  score_m <- numeric(n_obs)
  score_k <- numeric(n_obs)
  mass <- numeric(n_obs)
  mass_score_n <- numeric(n_obs)
  mass_score_m <- numeric(n_obs)
  mass_score_k <- numeric(n_obs)

  active <- which(valid)
  if (length(active)) {
    p_x <- numeric(n_obs)
    s_n <- numeric(n_obs)
    s_m <- numeric(n_obs)
    s_k <- numeric(n_obs)
    p_x[active] <- exp(
      lbeta(shape_n[active], m[active] + shape_k[active]) -
        lbeta(shape_n[active], m[active])
    )
    s_n[active] <- digamma(shape_n[active] + m[active]) -
      digamma(shape_n[active] + m[active] + shape_k[active])
    s_m[active] <- digamma(m[active] + shape_k[active]) -
      digamma(shape_n[active] + m[active] + shape_k[active]) -
      digamma(m[active]) + digamma(shape_n[active] + m[active])
    s_k[active] <- digamma(m[active] + shape_k[active]) -
      digamma(shape_n[active] + m[active] + shape_k[active])

    probability[active] <- p_x[active]
    score_n[active] <- p_x[active] * s_n[active]
    score_m[active] <- p_x[active] * s_m[active]
    score_k[active] <- p_x[active] * s_k[active]
    at_mass <- active[q[active] == 0L]
    if (length(at_mass)) {
      mass[at_mass] <- p_x[at_mass]
      mass_score_n[at_mass] <- p_x[at_mass] * s_n[at_mass]
      mass_score_m[at_mass] <- p_x[at_mass] * s_m[at_mass]
      mass_score_k[at_mass] <- p_x[at_mass] * s_k[at_mass]
    }

    max_q <- max(q[active])
    if (max_q >= 1L) {
      for (x in 0:(max_q - 1L)) {
        next_active <- active[q[active] >= x + 1L]
        if (!length(next_active)) break
        denominator <- x + shape_n[next_active] + m[next_active] + shape_k[next_active]
        ratio <- ((x + shape_n[next_active]) / denominator) *
          ((x + shape_k[next_active]) / (x + 1))
        p_x[next_active] <- p_x[next_active] * ratio
        s_n[next_active] <- s_n[next_active] +
          1 / (x + shape_n[next_active]) - 1 / denominator
        s_m[next_active] <- s_m[next_active] - 1 / denominator
        s_k[next_active] <- s_k[next_active] - 1 / denominator +
          1 / (x + shape_k[next_active])

        probability[next_active] <- probability[next_active] + p_x[next_active]
        score_n[next_active] <- score_n[next_active] +
          p_x[next_active] * s_n[next_active]
        score_m[next_active] <- score_m[next_active] +
          p_x[next_active] * s_m[next_active]
        score_k[next_active] <- score_k[next_active] +
          p_x[next_active] * s_k[next_active]
        at_mass <- next_active[q[next_active] == x + 1L]
        if (length(at_mass)) {
          mass[at_mass] <- p_x[at_mass]
          mass_score_n[at_mass] <- p_x[at_mass] * s_n[at_mass]
          mass_score_m[at_mass] <- p_x[at_mass] * s_m[at_mass]
          mass_score_k[at_mass] <- p_x[at_mass] * s_k[at_mass]
        }
      }
    }
  }

  d_n <- list(upper = score_n, lower = score_n - mass_score_n)
  d_m <- list(upper = score_m, lower = score_m - mass_score_m)
  d_k <- list(upper = score_k, lower = score_k - mass_score_k)
  positive_lower <- valid & q >= 1L
  zero_lower <- valid & !positive_lower

  derivatives <- list(
    mu = list(
      upper = (1 - tau) * d_n$upper * nu / sigma,
      lower = (1 - tau) * d_n$lower * nu / sigma
    ),
    sigma = list(
      upper = (1 - tau) * (
        d_n$upper * (-shape_n / sigma) + d_m$upper * (-1 / sigma^2)
      ),
      lower = (1 - tau) * (
        d_n$lower * (-shape_n / sigma) + d_m$lower * (-1 / sigma^2)
      )
    ),
    nu = list(
      upper = (1 - tau) * (
        d_n$upper * mu / sigma + d_k$upper * (-1 / nu^2)
      ),
      lower = (1 - tau) * (
        d_n$lower * mu / sigma + d_k$lower * (-1 / nu^2)
      )
    ),
    tau = list(
      upper = 1 - probability,
      lower = 1 - (probability - mass)
    )
  )
  for (parameter in names(derivatives)) {
    derivatives[[parameter]]$upper[!valid] <- 0
    derivatives[[parameter]]$lower[!positive_lower] <- 0
    derivatives[[parameter]]$lower[zero_lower] <- 0
  }
  derivatives[intersect(par_names, names(derivatives))]
}

#' @noRd
.calc_F_bounds_derivatives <- function(
    eta_inv,
    mm,
    margin_dist,
    response,
    par_names = NULL,
    h = 1e-4,
    cdf_fn = calc_F_x) {
  if (is.list(mm) && all(c("x", "s") %in% names(mm))) {
    mm <- mm$x
  }

  margin_par_names <- names(eta_inv)[names(eta_inv) %in% c("mu", "sigma", "nu", "tau")]
  if (!is.null(par_names)) {
    margin_par_names <- intersect(margin_par_names, par_names)
  }

  if (identical(as.character(margin_dist$family[1L]), "ZIBNB") &&
      isTRUE(getOption("gamlss.longitudinal.fast_zibnb_cdf_scores", TRUE))) {
    fast <- .calc_zibnb_F_bounds_derivatives(
      eta_inv = eta_inv,
      response = response,
      par_names = margin_par_names
    )
    if (!is.null(fast)) return(fast)
  }

  out <- setNames(vector("list", length(margin_par_names)), margin_par_names)
  for (pn in margin_par_names) {
    hp <- h * pmax(1, abs(as.numeric(eta_inv[[pn]])))
    eta_plus <- eta_minus <- eta_inv
    eta_plus[[pn]] <- eta_plus[[pn]] + hp
    eta_minus[[pn]] <- eta_minus[[pn]] - hp
    if (pn %in% c("mu", "sigma", "tau", "zeta")) {
      eta_plus[[pn]] <- pmax(eta_plus[[pn]], 1e-8)
      eta_minus[[pn]] <- pmax(eta_minus[[pn]], 1e-8)
    }
    if (pn == "nu" && identical(as.character(margin_dist$family[1]), "DEL")) {
      eta_plus[[pn]] <- pmin(pmax(eta_plus[[pn]], 1e-8), 1 - 1e-8)
      eta_minus[[pn]] <- pmin(pmax(eta_minus[[pn]], 1e-8), 1 - 1e-8)
    }
    denom <- eta_plus[[pn]] - eta_minus[[pn]]
    denom[!is.finite(denom) | abs(denom) < .Machine$double.eps] <- NA_real_

    # A family CDF can fail on only one side of an otherwise valid numerical
    # perturbation (for example, pZISICHEL at an extreme but finite sigma).
    # Preserve the central difference wherever both sides are evaluable and
    # use the valid base likelihood state for a one-sided derivative otherwise.
    # Both invalid sides remain NA and are rejected by the validator below.
    Fu_plus <- cdf_fn(eta_plus, mm, margin_dist, response, strict = FALSE)
    Fu_minus <- cdf_fn(eta_minus, mm, margin_dist, response, strict = FALSE)
    Fl_plus <- cdf_fn(eta_plus, mm, margin_dist, response - 1, strict = FALSE)
    Fl_minus <- cdf_fn(eta_minus, mm, margin_dist, response - 1, strict = FALSE)

    upper <- (Fu_plus - Fu_minus) / denom
    lower <- (Fl_plus - Fl_minus) / denom
    needs_base <-
      xor(is.finite(Fu_plus), is.finite(Fu_minus)) |
      xor(is.finite(Fl_plus), is.finite(Fl_minus))
    if (any(needs_base)) {
      Fu_base <- cdf_fn(eta_inv, mm, margin_dist, response)
      Fl_base <- cdf_fn(eta_inv, mm, margin_dist, response - 1)
      plus_denom <- eta_plus[[pn]] - eta_inv[[pn]]
      minus_denom <- eta_inv[[pn]] - eta_minus[[pn]]

      upper_plus_only <- is.finite(Fu_plus) & !is.finite(Fu_minus)
      upper_minus_only <- !is.finite(Fu_plus) & is.finite(Fu_minus)
      lower_plus_only <- is.finite(Fl_plus) & !is.finite(Fl_minus)
      lower_minus_only <- !is.finite(Fl_plus) & is.finite(Fl_minus)
      upper[upper_plus_only] <-
        (Fu_plus[upper_plus_only] - Fu_base[upper_plus_only]) /
        plus_denom[upper_plus_only]
      upper[upper_minus_only] <-
        (Fu_base[upper_minus_only] - Fu_minus[upper_minus_only]) /
        minus_denom[upper_minus_only]
      lower[lower_plus_only] <-
        (Fl_plus[lower_plus_only] - Fl_base[lower_plus_only]) /
        plus_denom[lower_plus_only]
      lower[lower_minus_only] <-
        (Fl_base[lower_minus_only] - Fl_minus[lower_minus_only]) /
        minus_denom[lower_minus_only]
    }

    out[[pn]] <- list(
      upper = .gl_validate_margin_cdf_derivative(
        upper,
        response = response,
        margin_dist = margin_dist,
        context = "Upper discrete rectangle-bound derivative",
        parameter = pn
      ),
      lower = .gl_validate_margin_cdf_derivative(
        lower,
        response = response,
        margin_dist = margin_dist,
        context = "Lower discrete rectangle-bound derivative",
        parameter = pn
      )
    )
  }
  out
}
