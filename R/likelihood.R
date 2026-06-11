#' First derivative of a copula CDF with respect to the first margin
#'
#' Wrapper around the h-function used by likelihood and derivative assembly.
#'
#' @noRd
.copula_cdf_du1 <- function(u1, u2, family, par, par2 = 0) {
  .copula_hfunc1(u1, u2, family = family, par = par, par2 = par2)
}

#' First derivative of a copula CDF with respect to the second margin
#'
#' Uses symmetry of the h-function implementation.
#'
#' @noRd
.copula_cdf_du2 <- function(u1, u2, family, par, par2 = 0) {
  .copula_hfunc1(u2, u1, family = family, par = par, par2 = par2)
}

#' Cached Gaussian copula CDF
#'
#' Evaluates the Gaussian copula CDF with boundary shortcuts and duplicated-row
#' caching. Uses `mvtnorm` when available, otherwise falls back to `VineCopula`.
#'
#' @noRd
.copula_gaussian_cdf_cached <- function(u1, u2, par) {
  vals <- .copula_recycle(as.numeric(u1), as.numeric(u2), .copula_gaussian_rho(par))
  u1 <- vals[[1]]
  u2 <- vals[[2]]
  rho <- vals[[3]]
  out <- rep(NA_real_, length(u1))

  finite <- is.finite(u1) & is.finite(u2) & is.finite(rho)
  if (!any(finite)) {
    return(out)
  }

  u1 <- pmin(pmax(u1, 0), 1)
  u2 <- pmin(pmax(u2, 0), 1)

  zero <- finite & (u1 <= 0 | u2 <= 0)
  out[zero] <- 0

  one1 <- finite & !zero & u1 >= 1
  out[one1] <- u2[one1]

  one2 <- finite & !zero & !one1 & u2 >= 1
  out[one2] <- u1[one2]

  independent <- finite & !zero & !one1 & !one2 & abs(rho) <= 1e-12
  out[independent] <- u1[independent] * u2[independent]

  remaining <- finite & !zero & !one1 & !one2 & !independent
  if (any(remaining)) {
    key <- paste(signif(u1[remaining], 15), signif(u2[remaining], 15), signif(rho[remaining], 15), sep = "\r")
    unique_key <- !duplicated(key)
    u1_unique <- u1[remaining][unique_key]
    u2_unique <- u2[remaining][unique_key]
    rho_unique <- rho[remaining][unique_key]
    if (requireNamespace("mvtnorm", quietly = TRUE)) {
      z1 <- stats::qnorm(u1_unique)
      z2 <- stats::qnorm(u2_unique)
      unique_vals <- vapply(seq_along(rho_unique), function(ii) {
        as.numeric(mvtnorm::pmvnorm(
          upper = c(z1[[ii]], z2[[ii]]),
          corr = matrix(c(1, rho_unique[[ii]], rho_unique[[ii]], 1), 2, 2)
        ))
      }, numeric(1))
    } else {
      .copula_require_vinecopula("the Gaussian copula CDF backend")
      unique_vals <- VineCopula::BiCopCDF(
        u1_unique,
        u2_unique,
        family = 1,
        par = rho_unique
      )
    }
    out[remaining] <- unique_vals[match(key, key[unique_key])]
  }

  pmin(pmax(as.numeric(out), 0), 1)
}

#' Cached Gaussian rectangle probability
#'
#' Computes discrete/count likelihood rectangle masses from four cached
#' Gaussian CDF evaluations.
#'
#' @noRd
.copula_gaussian_rectangle_prob_cached <- function(u1, u2, l1, l2, par) {
  vals <- .copula_recycle(as.numeric(u1), as.numeric(u2), as.numeric(l1), as.numeric(l2), .copula_gaussian_rho(par))
  n <- length(vals[[1]])
  if (n == 0L) {
    return(numeric(0))
  }
  cdf <- .copula_gaussian_cdf_cached(
    c(vals[[1]], vals[[3]], vals[[1]], vals[[3]]),
    c(vals[[2]], vals[[2]], vals[[4]], vals[[4]]),
    rep(vals[[5]], 4L)
  )
  rect <- cdf[seq_len(n)] -
    cdf[n + seq_len(n)] -
    cdf[2L * n + seq_len(n)] +
    cdf[3L * n + seq_len(n)]
  pmax(as.numeric(rect), 1e-300)
}

#' Copula rectangle probability
#'
#' Returns positive lower-bounded rectangle probabilities for discrete margins,
#' using the cached Gaussian path where possible and generic inclusion-exclusion
#' for other copulas.
#'
#' @noRd
.copula_rectangle_prob <- function(u1, u2, l1, l2, family, par, par2 = 0) {
  family <- .copula_family_code(family)
  if (identical(family, "N")) {
    return(.copula_gaussian_rectangle_prob_cached(u1, u2, l1, l2, par))
  }
  rect <- .copula_cdf(u1, u2, family = family, par = par, par2 = par2) -
    .copula_cdf(l1, u2, family = family, par = par, par2 = par2) -
    .copula_cdf(u1, l2, family = family, par = par, par2 = par2) +
    .copula_cdf(l1, l2, family = family, par = par, par2 = par2)
  pmax(as.numeric(rect), 1e-300)
}

#' Finite-difference derivatives of marginal CDF bounds
#'
#' Computes observation-level derivatives of CDF bounds with respect to natural
#' margin parameters. These feed the discrete rectangle likelihood derivatives.
#'
#' @noRd
.calc_F_bounds_derivatives <- function(eta_inv, mm, margin_dist, response, par_names = NULL, h = 1e-4) {
  if (is.list(mm) && all(c("x", "s") %in% names(mm))) {
    mm <- mm$x
  }

  margin_par_names <- names(eta_inv)[names(eta_inv) %in% c("mu", "sigma", "nu", "tau")]
  if (!is.null(par_names)) {
    margin_par_names <- intersect(margin_par_names, par_names)
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

    Fu_plus <- calc_F_x(eta_plus, mm, margin_dist, response)
    Fu_minus <- calc_F_x(eta_minus, mm, margin_dist, response)
    Fl_plus <- calc_F_x(eta_plus, mm, margin_dist, response - 1)
    Fl_minus <- calc_F_x(eta_minus, mm, margin_dist, response - 1)

    out[[pn]] <- list(
      upper = as.numeric((Fu_plus - Fu_minus) / denom),
      lower = as.numeric((Fl_plus - Fl_minus) / denom)
    )
    out[[pn]]$upper[!is.finite(out[[pn]]$upper)] <- 0
    out[[pn]]$lower[!is.finite(out[[pn]]$lower)] <- 0
  }
  out
}

.calc_copula_rectangle_par_derivative <- function(calc_lik, copula_dist, par_name, h = 1e-5) {
  row_id1 <- calc_lik$copula_row_id1
  n_pair <- length(row_id1)
  if (n_pair == 0L) return(numeric(0))

  par1 <- calc_lik$copula_par1
  par2 <- calc_lik$copula_par2
  hp <- h * pmax(1, abs(if (identical(par_name, "theta")) par1 else par2))
  par1_p <- par1_m <- par1
  par2_p <- par2_m <- par2
  if (identical(par_name, "theta")) {
    par1_p <- par1 + hp
    par1_m <- par1 - hp
    if (identical(copula_dist, "C")) par1_m <- pmax(par1_m, 1e-8)
  } else {
    par2_p <- par2 + hp
    par2_m <- pmax(par2 - hp, 1e-8)
  }

  u <- calc_lik$Fx_1_2
  l <- calc_lik$Fx_1_2_lower
  rect_p <- .copula_rectangle_prob(u[, 1], u[, 2], l[, 1], l[, 2], copula_dist, par1_p, par2_p)
  rect_m <- .copula_rectangle_prob(u[, 1], u[, 2], l[, 1], l[, 2], copula_dist, par1_m, par2_m)
  deriv <- (rect_p - rect_m) / (2 * hp)
  deriv[!is.finite(deriv)] <- 0
  deriv
}

.calc_discrete_rectangle_scores <- function(
  eta_inv,
  mm,
  margin_dist,
  copula_dist,
  response,
  response_margin,
  response_subject,
  pair_cache = NULL,
  calc_lik = NULL,
  method = c("analytical", "finite"),
  h = 1e-5
) {
  method <- match.arg(method)
  if (is.null(pair_cache)) {
    pair_cache <- build_copula_pair_cache(response, response_margin, response_subject)
  }
  if (is.null(calc_lik)) {
    calc_lik <- calc_likelihood_minimal(
      eta_inv, mm = mm, margin_dist, copula_dist, response = response,
      response_margin = response_margin, response_subject = response_subject,
      pair_cache = pair_cache
    )
  }

  score_names <- intersect(names(mm), c("mu", "sigma", "nu", "tau", "theta", "zeta"))
  if (method == "finite") {
    base <- calc_lik$log_lik["joint"]
    out <- setNames(vector("list", length(score_names)), score_names)
    for (pn in score_names) {
      n <- length(eta_inv[[pn]])
      out[[pn]] <- numeric(n)
      for (ii in seq_len(n)) {
        hp <- h * max(1, abs(eta_inv[[pn]][ii]))
        plus <- minus <- eta_inv
        plus[[pn]][ii] <- plus[[pn]][ii] + hp
        minus[[pn]][ii] <- minus[[pn]][ii] - hp
        if (pn %in% c("mu", "sigma", "tau", "zeta")) {
          plus[[pn]][ii] <- max(plus[[pn]][ii], 1e-8)
          minus[[pn]][ii] <- max(minus[[pn]][ii], 1e-8)
        }
        if (pn == "nu" && identical(as.character(margin_dist$family[1]), "DEL")) {
          plus[[pn]][ii] <- min(max(plus[[pn]][ii], 1e-8), 1 - 1e-8)
          minus[[pn]][ii] <- min(max(minus[[pn]][ii], 1e-8), 1 - 1e-8)
        }
        if (pn == "theta" && identical(copula_dist, "C")) {
          minus[[pn]][ii] <- max(minus[[pn]][ii], 1e-8)
        }
        lp <- calc_likelihood_minimal(
          plus, mm = mm, margin_dist, copula_dist, response = response,
          response_margin = response_margin, response_subject = response_subject,
          pair_cache = pair_cache, calc_margin_deriv = FALSE
        )$log_lik["joint"]
        lm <- calc_likelihood_minimal(
          minus, mm = mm, margin_dist, copula_dist, response = response,
          response_margin = response_margin, response_subject = response_subject,
          pair_cache = pair_cache, calc_margin_deriv = FALSE
        )$log_lik["joint"]
        if (is.finite(lp) && is.finite(lm)) {
          out[[pn]][ii] <- (lp - lm) / (plus[[pn]][ii] - minus[[pn]][ii])
        } else if (is.finite(lp) && is.finite(base)) {
          out[[pn]][ii] <- (lp - base) / hp
        } else if (is.finite(lm) && is.finite(base)) {
          out[[pn]][ii] <- (base - lm) / hp
        }
      }
    }
    return(out)
  }

  margin_par <- intersect(score_names, c("mu", "sigma", "nu", "tau"))
  copula_par <- intersect(score_names, c("theta", "zeta"))
  out <- setNames(vector("list", length(score_names)), score_names)

  margin_deriv_subnames <- c(mu = "m", sigma = "d", nu = "v", tau = "t")
  for (pn in margin_par) {
    d_name <- paste0("dld", margin_deriv_subnames[pn])
    hit <- grep(paste0("^", d_name, "$"), names(calc_lik$margin_deriv))
    if (length(hit) == 0) hit <- grep(d_name, names(calc_lik$margin_deriv))
    out[[pn]] <- if (length(hit) == 0) rep(0, length(response)) else as.numeric(calc_lik$margin_deriv[[hit[1]]])
    out[[pn]][!is.finite(out[[pn]])] <- 0
  }

  row_id1 <- calc_lik$copula_row_id1
  row_id2 <- calc_lik$copula_row_id2
  pair_ok <- calc_lik$pair_complete
  rect <- calc_lik$copula_rect_prob
  u <- calc_lik$Fx_1_2
  l <- calc_lik$Fx_1_2_lower
  par1 <- calc_lik$copula_par1
  par2 <- calc_lik$copula_par2

  if (length(row_id1) > 0L && length(margin_par) > 0L) {
    F_deriv <- .calc_F_bounds_derivatives(eta_inv, mm, margin_dist, response, par_names = margin_par, h = h)
    du1_uu <- .copula_cdf_du1(u[, 1], u[, 2], copula_dist, par1, par2)
    du1_ul <- .copula_cdf_du1(u[, 1], l[, 2], copula_dist, par1, par2)
    dl1_lu <- .copula_cdf_du1(l[, 1], u[, 2], copula_dist, par1, par2)
    dl1_ll <- .copula_cdf_du1(l[, 1], l[, 2], copula_dist, par1, par2)
    du2_uu <- .copula_cdf_du2(u[, 1], u[, 2], copula_dist, par1, par2)
    du2_lu <- .copula_cdf_du2(l[, 1], u[, 2], copula_dist, par1, par2)
    dl2_ul <- .copula_cdf_du2(u[, 1], l[, 2], copula_dist, par1, par2)
    dl2_ll <- .copula_cdf_du2(l[, 1], l[, 2], copula_dist, par1, par2)

    for (pn in margin_par) {
      dU <- F_deriv[[pn]]$upper
      dL <- F_deriv[[pn]]$lower
      drect1 <- (du1_uu - du1_ul) * dU[row_id1] + (-dl1_lu + dl1_ll) * dL[row_id1]
      drect2 <- (du2_uu - du2_lu) * dU[row_id2] + (-dl2_ul + dl2_ll) * dL[row_id2]
      contrib1 <- drect1 / rect - out[[pn]][row_id1]
      contrib2 <- drect2 / rect - out[[pn]][row_id2]
      contrib1[!pair_ok | !is.finite(contrib1)] <- 0
      contrib2[!pair_ok | !is.finite(contrib2)] <- 0
      sum1 <- rowsum(contrib1, row_id1, reorder = FALSE)
      sum2 <- rowsum(contrib2, row_id2, reorder = FALSE)
      out[[pn]][as.integer(rownames(sum1))] <- out[[pn]][as.integer(rownames(sum1))] + sum1[, 1]
      out[[pn]][as.integer(rownames(sum2))] <- out[[pn]][as.integer(rownames(sum2))] + sum2[, 1]
    }
  }

  for (pn in copula_par) {
    n_par <- length(eta_inv[[pn]])
    out[[pn]] <- numeric(n_par)
    if (length(row_id1) == 0L) next
    drect <- .calc_copula_rectangle_par_derivative(calc_lik, copula_dist, pn, h = h)
    contrib <- drect / rect
    contrib[!pair_ok | !is.finite(contrib)] <- 0
    if (n_par == length(response)) {
      par_idx <- row_id1
    } else {
      par_idx <- calc_lik$copula_theta_index_map[row_id1]
    }
    valid <- is.finite(par_idx) & par_idx >= 1 & par_idx <= n_par
    if (any(valid)) {
      summed <- rowsum(contrib[valid], par_idx[valid], reorder = FALSE)
      out[[pn]][as.integer(rownames(summed))] <- summed[, 1]
    }
  }

  out
}

.score_list_to_beta_gradient <- function(score_list, eta_out, mm_cg, beta_names) {
  grad <- rep(0, length(beta_names))
  names(grad) <- beta_names
  for (pn in names(score_list)) {
    if (!pn %in% names(mm_cg$x) || !pn %in% names(eta_out$eta_dr)) next
    score_eta <- as.numeric(score_list[[pn]]) * as.numeric(eta_out$eta_dr[[pn]])
    score_eta[!is.finite(score_eta)] <- 0
    X <- as.matrix(mm_cg$x[[pn]])
    par_grad <- as.numeric(crossprod(X, score_eta))
    x_names <- colnames(X)
    names(par_grad) <- ifelse(
      startsWith(x_names, paste0(pn, ".")),
      x_names,
      paste(pn, x_names, sep = ".")
    )
    idx <- match(names(par_grad), names(grad))
    valid <- !is.na(idx)
    grad[idx[valid]] <- par_grad[valid]
  }
  grad
}

.cg_analytical_gradient <- function(
  beta_vec,
  mm_cg,
  eta_out,
  calc_lik,
  margin_dist,
  copula_dist,
  include_dlcopdpar,
  response,
  response_margin,
  response_subject
) {
  eta <- eta_out$eta
  eta_inv <- eta_out$eta_inv
  eta_dr <- eta_out$eta_dr

  grad <- rep(0, length(beta_vec))
  names(grad) <- names(beta_vec)

  margin_par <- intersect(names(mm_cg$x), c("mu", "sigma", "nu", "tau"))
  copula_par <- intersect(names(mm_cg$x), c("theta", "zeta"))

  if (identical(calc_lik$likelihood_type, "discrete_rectangle")) {
    score_list <- .calc_discrete_rectangle_scores(
      eta_inv,
      mm_cg$x,
      margin_dist,
      copula_dist,
      response,
      response_margin,
      response_subject,
      calc_lik = calc_lik,
      method = "analytical"
    )
    return(.score_list_to_beta_gradient(score_list, eta_out, mm_cg, names(beta_vec)))
  }

  copula_derivatives <- calc_copula_derivatives(
    eta_inv,
    calc_lik$Fx_1_2,
    copula_dist,
    par1 = calc_lik$copula_par1,
    par2 = calc_lik$copula_par2,
    pair_complete = calc_lik$pair_complete
  )

  margin_score_natural <- list()
  if(length(margin_par) > 0) {
    margin_deriv_subnames <- c("m", "d", "v", "t")
    names(margin_deriv_subnames) <- c("mu", "sigma", "nu", "tau")

    for(par_name in margin_par) {
      d_name <- paste0("dld", margin_deriv_subnames[par_name])
      hit <- grep(paste0("^", d_name, "$"), names(calc_lik$margin_deriv))
      if(length(hit) == 0) {
        hit <- grep(d_name, names(calc_lik$margin_deriv))
      }
      if(length(hit) == 0) {
        margin_score_natural[[par_name]] <- rep(0, length(eta[[par_name]]))
      } else {
        margin_score_natural[[par_name]] <- as.numeric(calc_lik$margin_deriv[[hit[1]]])
      }
    }

    if(isTRUE(include_dlcopdpar)) {
      nd_impact_F <- calc_Fx_derivatives(eta_inv, mm_cg$x, margin_dist, response = response)

      order_margin <- data.frame(time = response_margin, subject = response_subject)
      margin_deriv_1 <- matrix(0, ncol = length(margin_par), nrow = length(response))
      colnames(margin_deriv_1) <- paste0("dld", margin_par)
      for(par_name in margin_par) {
        margin_deriv_1[, paste0("dld", par_name)] <- margin_score_natural[[par_name]]
      }

      margin_components_base <- cbind(
        order_margin,
        response = response,
        margin_p = calc_lik$margin_p,
        margin_d = calc_lik$margin_d,
        margin_deriv_1,
        mu = eta_inv[["mu"]]
      )
      names(margin_components_base)[seq_len(ncol(order_margin))] <- c("time", "subject")

      copula_components <- cbind(
        calc_lik$order_copula,
        row_id1 = calc_lik$copula_row_id1,
        row_id2 = calc_lik$copula_row_id2,
        dcdu1 = copula_derivatives$dcdu1,
        dcdu2 = copula_derivatives$dcdu2,
        copula_d = calc_lik$copula_d
      )

      for(par_name in margin_par) {
        margin_components <- cbind(
          margin_components_base,
          F_nd = nd_impact_F[[par_name]]
        )
        margin_components_Ft_plus <- margin_components
        margin_components_Ft_plus$time <- normalize_lag_time(margin_components_Ft_plus$time)
        margin_plus <- merge(
          margin_components,
          margin_components_Ft_plus,
          by = c("time", "subject"),
          all.x = TRUE
        )
        copula_merged <- merge(
          copula_components,
          margin_plus,
          by.x = c("time1", "subject1"),
          by.y = c("time", "subject"),
          all.x = TRUE
        )
        d1_cop <- calc_deriv_copula_wrt_margin(
          copula_merged,
          margin_par,
          par_name,
          calc_d2 = FALSE
        )[, which(margin_par == par_name)]
        n_score <- length(margin_score_natural[[par_name]])
        if(length(d1_cop) >= n_score) {
          margin_score_natural[[par_name]] <- margin_score_natural[[par_name]] + d1_cop[seq_len(n_score)]
        }
      }
    }
  }

  for(par_name in margin_par) {
    score_eta <- as.numeric(margin_score_natural[[par_name]]) * as.numeric(eta_dr[[par_name]])
    score_eta[!is.finite(score_eta)] <- 0
    X <- as.matrix(mm_cg$x[[par_name]])
    par_grad <- as.numeric(crossprod(X, score_eta))
    x_names <- colnames(X)
    names(par_grad) <- ifelse(
      startsWith(x_names, paste0(par_name, ".")),
      x_names,
      paste(par_name, x_names, sep = ".")
    )
    idx <- match(names(par_grad), names(grad))
    valid <- !is.na(idx)
    grad[idx[valid]] <- par_grad[valid]
  }

  if("theta" %in% copula_par) {
    n_par <- length(eta[["theta"]])
    d1_full <- rep(0, n_par)
    row_id1 <- calc_lik$copula_row_id1
    if(length(row_id1) > 0) {
      if(n_par == length(response)) {
        par_idx <- row_id1
      } else {
        par_idx <- calc_lik$copula_theta_index_map[row_id1]
      }
      valid_idx <- is.finite(par_idx) & par_idx >= 1 & par_idx <= n_par
      if(any(valid_idx)) {
        d1_sum <- rowsum(copula_derivatives$dldth[valid_idx], par_idx[valid_idx], reorder = FALSE)
        d1_full[as.integer(rownames(d1_sum))] <- d1_sum[, 1]
      }
    }
    score_eta <- d1_full * as.numeric(eta_dr[["theta"]])
    score_eta[!is.finite(score_eta)] <- 0
    par_grad <- as.numeric(crossprod(as.matrix(mm_cg$x[["theta"]]), score_eta))
    x_names <- colnames(mm_cg$x[["theta"]])
    names(par_grad) <- ifelse(
      startsWith(x_names, "theta."),
      x_names,
      paste("theta", x_names, sep = ".")
    )
    idx <- match(names(par_grad), names(grad))
    valid <- !is.na(idx)
    grad[idx[valid]] <- par_grad[valid]
  }

  if("zeta" %in% copula_par && "dldz" %in% names(copula_derivatives)) {
    n_par <- length(eta[["zeta"]])
    d1_full <- rep(0, n_par)
    row_id1 <- calc_lik$copula_row_id1
    if(length(row_id1) > 0) {
      if(n_par == length(response)) {
        par_idx <- row_id1
      } else {
        par_idx <- calc_lik$copula_theta_index_map[row_id1]
      }
      valid_idx <- is.finite(par_idx) & par_idx >= 1 & par_idx <= n_par
      if(any(valid_idx)) {
        d1_sum <- rowsum(copula_derivatives$dldz[valid_idx], par_idx[valid_idx], reorder = FALSE)
        d1_full[as.integer(rownames(d1_sum))] <- d1_sum[, 1]
      }
    }
    score_eta <- d1_full * as.numeric(eta_dr[["zeta"]])
    score_eta[!is.finite(score_eta)] <- 0
    par_grad <- as.numeric(crossprod(as.matrix(mm_cg$x[["zeta"]]), score_eta))
    x_names <- colnames(mm_cg$x[["zeta"]])
    names(par_grad) <- ifelse(
      startsWith(x_names, "zeta."),
      x_names,
      paste("zeta", x_names, sep = ".")
    )
    idx <- match(names(par_grad), names(grad))
    valid <- !is.na(idx)
    grad[idx[valid]] <- par_grad[valid]
  }

  grad
}

#' Calculate the likelihood components for the joint model
#'
#' This function calculates the marginal and copula log likelihoods and components
#' for the joint model by organizing response data by margin and subject for
#' efficient pair-based copula calculations.
#'
#' @param response A numeric vector of response values.
#' @param response_margin A numeric vector indicating the margin (time) for each response.
#' @param response_subject A numeric vector indicating the subject for each response.
#' @return A list containing:
#' \item{log_lik}{A named vector with marginal, copula, and joint log-likelihoods.}
#' \item{margin_d}{A numeric vector of marginal densities.}
#' \item{copula_d}{A numeric vector of copula densities.}
#' \item{margin_p}{A numeric vector of marginal distribution function values.}
#' \item{Fx_1_2}{A matrix of marginal distribution function values for pairs of margins.}
#' \item{order_copula}{A matrix indicating the order of margins and subjects for copula calculations.}
#' \item{margin_deriv}{A list of marginal derivatives.}
#' \item{order_copula}{A matrix indicating the order of margins and subjects for copula calculations.}
#'
#' @keywords internal
#' @noRd
build_copula_pair_cache <- function(response, response_margin, response_subject) {
  margin_names=sort(unique(response_margin))
  num_margins=length(margin_names)
  n_obs=length(response)
  obs_response=!is.na(response)

  base_df=data.frame(
    row_id=seq_len(n_obs),
    time=response_margin,
    subject=response_subject,
    observed=obs_response,
    stringsAsFactors = FALSE
  )

  pair_df_all=list()
  if(num_margins>1) {
    for (i in seq_len(num_margins-1)) {
      t1=margin_names[i]
      t2=margin_names[i+1]

      left=base_df[base_df$time==t1,c("row_id","subject","time","observed")]
      right=base_df[base_df$time==t2,c("row_id","subject","time","observed")]
      names(left)=c("row_id1","subject","time1","observed1")
      names(right)=c("row_id2","subject","time2","observed2")

      pair_i=merge(left,right,by="subject",all=FALSE)
      if(nrow(pair_i)>0) {
        pair_df_all[[length(pair_df_all)+1]]=pair_i
      }
    }
  }

  if(length(pair_df_all)==0) {
    pair_df=data.frame(
      subject=response_subject[0],
      row_id1=integer(0),
      time1=response_margin[0],
      observed1=logical(0),
      row_id2=integer(0),
      time2=response_margin[0],
      observed2=logical(0)
    )
  } else {
    pair_df=do.call(rbind,pair_df_all)
  }

  order_copula=as.matrix(pair_df[,c("time1","subject","time2","subject")])
  colnames(order_copula)=c("time1","subject1","time2","subject2")

  observed_pair_base=rep(FALSE,nrow(pair_df))
  if(nrow(pair_df)>0) {
    observed_pair_base=pair_df$observed1 & pair_df$observed2
  }

  Fx_1_2_template=matrix(NA_real_,nrow=nrow(pair_df),ncol=2)
  colnames(Fx_1_2_template)=c("u1","u2")

  theta_rows=which(response_margin %in% margin_names[seq_len(max(1, num_margins-1))])
  theta_index_map=rep(NA_integer_,n_obs)
  theta_index_map[theta_rows]=seq_along(theta_rows)

  cache=list(
    row_id1=pair_df$row_id1,
    row_id2=pair_df$row_id2,
    Fx_1_2_template=Fx_1_2_template,
    order_copula=order_copula,
    observed_pair_base=observed_pair_base,
    theta_index_map=theta_index_map,
    margin_names=margin_names,
    num_margins=num_margins,
    n_obs=n_obs
  )
  cache
}

.build_margin_eval_cache <- function(margin_dist, calc_d2=FALSE) {
  if(calc_d2==TRUE) {
    to_include=grepl("dld",names(margin_dist))|grepl("d2ld",names(margin_dist))
  } else {
    to_include=grepl("dld",names(margin_dist))
  }

  margin_deriv_names=names(margin_dist)[to_include]
  margin_deriv_cache=lapply(margin_deriv_names, function(deriv_name) {
    FUN=margin_dist[[deriv_name]]
    list(name=deriv_name,FUN=FUN,args=formalArgs(FUN))
  })

  margin_pFUN <- get(
    paste("p", margin_dist$family[1], sep = ""),
    envir = asNamespace("gamlss.dist"),
    mode = "function",
    inherits = FALSE
  )
  margin_dFUN <- get(
    paste("d", margin_dist$family[1], sep = ""),
    envir = asNamespace("gamlss.dist"),
    mode = "function",
    inherits = FALSE
  )

  list(
    calc_d2=calc_d2,
    margin_deriv_cache=margin_deriv_cache,
    margin_pFUN=margin_pFUN,
    margin_p_args=formalArgs(margin_pFUN),
    margin_dFUN=margin_dFUN,
    margin_d_args=formalArgs(margin_dFUN),
    family_call_cache = new.env(parent = emptyenv())
  )
}

.call_margin_family_cached <- function(FUN, args, FUN_args, cacheable = FALSE, cache_env = NULL, cache_prefix = "") {
  call_args <- args[FUN_args]
  safe_call <- function(call_args, fallback_n = NULL) {
    if (is.null(fallback_n)) {
      fallback_n <- max(c(1L, vapply(call_args, length, integer(1))), na.rm = TRUE)
    }
    tryCatch(
      do.call(FUN, args = call_args),
      error = function(e) rep(NA_real_, fallback_n)
    )
  }
  if (!isTRUE(cacheable) || length(call_args) == 0L) {
    return(safe_call(call_args))
  }

  arg_lengths <- vapply(call_args, length, integer(1))
  n <- max(arg_lengths)
  if (n <= 1L) {
    return(safe_call(call_args, fallback_n = n))
  }

  if (any(!arg_lengths %in% c(1L, n))) {
    return(safe_call(call_args, fallback_n = n))
  }
  if (!all(vapply(call_args, function(x) is.atomic(x) && !is.character(x), logical(1)))) {
    return(safe_call(call_args, fallback_n = n))
  }

  expanded <- lapply(call_args, rep, length.out = n)
  key_parts <- lapply(expanded, function(x) {
    if (is.numeric(x) || is.integer(x)) {
      format(signif(as.numeric(x), 15), scientific = FALSE, trim = TRUE)
    } else {
      as.character(x)
    }
  })
  key <- do.call(paste, c(key_parts, sep = "\r"))
  unique_idx <- !duplicated(key)
  if (sum(unique_idx) > 0.8 * n) {
    return(safe_call(call_args, fallback_n = n))
  }

  if (!is.null(cache_env)) {
    cached <- rep(NA_real_, n)
    hit <- logical(n)
    lookup_key <- paste0(cache_prefix, "\r", key)
    for (ii in seq_len(n)) {
      if (exists(lookup_key[[ii]], envir = cache_env, inherits = FALSE)) {
        cached[[ii]] <- get(lookup_key[[ii]], envir = cache_env, inherits = FALSE)
        hit[[ii]] <- TRUE
      }
    }
    if (all(hit)) {
      return(cached)
    }
    need <- !hit
    need_unique <- !duplicated(key[need])
    need_args <- lapply(expanded, function(x) x[need][need_unique])
    names(need_args) <- names(call_args)
    need_val <- safe_call(need_args, fallback_n = sum(need_unique))
    if (length(need_val) != sum(need_unique)) {
      return(safe_call(call_args, fallback_n = n))
    }
    need_lookup <- paste0(cache_prefix, "\r", key[need][need_unique])
    for (jj in seq_along(need_lookup)) {
      assign(need_lookup[[jj]], need_val[[jj]], envir = cache_env)
    }
    cached[need] <- need_val[match(key[need], key[need][need_unique])]
    return(cached)
  }

  unique_args <- lapply(expanded, `[`, unique_idx)
  names(unique_args) <- names(call_args)
  unique_val <- safe_call(unique_args, fallback_n = sum(unique_idx))
  if (length(unique_val) != sum(unique_idx)) {
    return(safe_call(call_args, fallback_n = n))
  }
  unique_val[match(key, key[unique_idx])]
}

.call_fast_count_family <- function(prefix, family_name, args) {
  if (!identical(family_name, "PO") || !"mu" %in% names(args)) {
    return(NULL)
  }
  x <- args$x %||% args$q %||% args$y
  mu <- args$mu
  if (is.null(x) || is.null(mu)) {
    return(NULL)
  }
  if (identical(prefix, "d")) {
    return(stats::dpois(x, lambda = mu))
  }
  if (identical(prefix, "p")) {
    return(stats::ppois(x, lambda = mu))
  }
  NULL
}

#' Calculate minimal joint likelihood components
#'
#' @param pair_cache Optional cache built by build_copula_pair_cache to reuse pair indexing across repeated likelihood calls.
#'
#' @noRd
calc_likelihood_minimal <- function(eta_inv,mm,margin_dist,copula_dist,calc_d2=FALSE,response,response_margin,response_subject,penalize_smooth=FALSE,par_s=NA,pair_cache=NULL,margin_eval_cache=NULL,calc_margin_deriv=TRUE) {
  #Setup input matrix of response and parameters
  #response=dataset$response; response_subject=dataset$subject; response_margin=dataset$time; dataset=NA
  if(is.null(pair_cache)) {
    pair_cache=build_copula_pair_cache(response,response_margin,response_subject)
  }
  if(is.null(margin_eval_cache) || !identical(margin_eval_cache$calc_d2, calc_d2)) {
    margin_eval_cache=.build_margin_eval_cache(margin_dist, calc_d2=calc_d2)
  }

  margin_names=pair_cache$margin_names
  num_margins=pair_cache$num_margins
  n_obs=pair_cache$n_obs
  discrete_margin <- .is_discrete_margin(margin_dist)

  obs_response=!is.na(response)

  order_margin=cbind(response_margin,response_subject)
  colnames(order_margin)=c("time","subject")

  margin_deriv_input=list()
  margin_deriv_input[["y"]]=response
  margin_deriv_input[["q"]]=response
  margin_deriv_input[["x"]]=response
  for (par_name in names(mm)) {
    if (par_name %in% c("mu","sigma","nu","tau")) {
      margin_deriv_input[[par_name]]=eta_inv[[par_name]]
    }
  }
  fixed_unlinked_values <- attr(margin_dist, "fixed_unlinked_values")
  if (length(fixed_unlinked_values) > 0L) {
    for (par_name in names(fixed_unlinked_values)) {
      value <- fixed_unlinked_values[[par_name]]
      if (is.numeric(value) && length(value) == 1L && is.finite(value)) {
        margin_deriv_input[[par_name]] <- rep(value, n_obs)
      }
    }
  }

  ################## MARGIN DERIVATIVES
  margin_deriv=list()
  if(isTRUE(calc_margin_deriv)) {
    for (deriv_info in margin_eval_cache$margin_deriv_cache) {
      FUN_args=names(margin_deriv_input)[names(margin_deriv_input)%in%deriv_info$args]
      deriv_val <- tryCatch(
        do.call(deriv_info$FUN, args = margin_deriv_input[FUN_args]),
        error = function(e) rep(0, n_obs)
      )
      if(length(deriv_val)==n_obs) {
        deriv_val[!obs_response]=0
        deriv_val[!is.finite(deriv_val)]=0
      }
      margin_deriv[[deriv_info$name]]=deriv_val
    }
  }

  FUN_args=names(margin_deriv_input)[names(margin_deriv_input)%in%margin_eval_cache$margin_p_args]
  margin_p <- .call_fast_count_family("p", margin_dist$family[1], margin_deriv_input)
  if (is.null(margin_p)) {
    margin_p=.call_margin_family_cached(
    margin_eval_cache$margin_pFUN,
    margin_deriv_input,
    FUN_args,
    cacheable = discrete_margin,
    cache_env = margin_eval_cache$family_call_cache,
    cache_prefix = paste0(margin_dist$family[1], ":p:upper")
    )
  }
  margin_p[!obs_response]=NA
  margin_p[!is.finite(margin_p)]=NA

  margin_p_lower <- NULL
  likelihood_type <- "continuous_density"
  if (discrete_margin) {
    margin_deriv_input_lower <- margin_deriv_input
    lower_response <- response - 1
    margin_deriv_input_lower[["q"]] <- lower_response
    margin_deriv_input_lower[["x"]] <- lower_response
    margin_p_lower <- rep(NA_real_, length(response))
    negative_lower <- obs_response & is.finite(lower_response) & lower_response < 0
    margin_p_lower[negative_lower] <- 0
    valid_lower <- obs_response & is.finite(lower_response) & lower_response >= 0
    if (any(valid_lower)) {
      lower_call_args <- lapply(margin_deriv_input_lower, function(value) {
        if (length(value) == n_obs) value[valid_lower] else value
      })
      lower_eval <- .call_fast_count_family("p", margin_dist$family[1], lower_call_args)
      if (is.null(lower_eval)) {
        lower_eval <- .call_margin_family_cached(
          margin_eval_cache$margin_pFUN,
          lower_call_args,
          names(lower_call_args)[names(lower_call_args) %in% margin_eval_cache$margin_p_args],
          cacheable = TRUE,
          cache_env = margin_eval_cache$family_call_cache,
          cache_prefix = paste0(margin_dist$family[1], ":p:lower")
        )
      }
      margin_p_lower[valid_lower] <- lower_eval
    }
    margin_p_lower[!obs_response] <- NA
    margin_p_lower[!is.finite(margin_p_lower)] <- NA
    likelihood_type <- "discrete_rectangle"
  }

  FUN_args=names(margin_deriv_input)[names(margin_deriv_input)%in%margin_eval_cache$margin_d_args]
  margin_d <- .call_fast_count_family("d", margin_dist$family[1], margin_deriv_input)
  if (is.null(margin_d)) {
    margin_d=.call_margin_family_cached(
    margin_eval_cache$margin_dFUN,
    margin_deriv_input,
    FUN_args,
    cacheable = discrete_margin,
    cache_env = margin_eval_cache$family_call_cache,
    cache_prefix = paste0(margin_dist$family[1], ":d")
    )
  }
  margin_d[!obs_response]=NA
  margin_d[!is.finite(margin_d) | margin_d<=0]=NA

  ################COPULA DERIVATIVES
  #First calculate margin F(x1), F(x2) as inputs to copula

  row_id1=pair_cache$row_id1
  row_id2=pair_cache$row_id2
  order_copula=pair_cache$order_copula

  Fx_1_2=pair_cache$Fx_1_2_template
  if(is.null(Fx_1_2)) {
    Fx_1_2=matrix(NA_real_,nrow=length(row_id1),ncol=2)
    colnames(Fx_1_2)=c("u1","u2")
  }
  if(length(row_id1)>0) {
    Fx_1_2[,1]=margin_p[row_id1]
    Fx_1_2[,2]=margin_p[row_id2]
  }

  Fx_1_2_lower <- NULL
  if (identical(likelihood_type, "discrete_rectangle")) {
    Fx_1_2_lower <- Fx_1_2
    if(length(row_id1)>0) {
      Fx_1_2_lower[,1]=margin_p_lower[row_id1]
      Fx_1_2_lower[,2]=margin_p_lower[row_id2]
    }
  }

  pair_complete=pair_cache$observed_pair_base & is.finite(Fx_1_2[,1]) & is.finite(Fx_1_2[,2])
  if (identical(likelihood_type, "discrete_rectangle")) {
    pair_complete <- pair_complete & is.finite(Fx_1_2_lower[,1]) & is.finite(Fx_1_2_lower[,2]) &
      is.finite(margin_d[row_id1]) & is.finite(margin_d[row_id2]) &
      margin_d[row_id1] > 0 & margin_d[row_id2] > 0
  }

  par1=rep(NA_real_,length(row_id1))
  par2=rep(NA_real_,length(row_id1))
  if(length(row_id1)>0) {
    theta_len <- length(eta_inv[["theta"]])
    if(theta_len == n_obs) {
      theta_idx <- row_id1
    } else {
      theta_idx=pair_cache$theta_index_map[row_id1]
    }

    par1=eta_inv[["theta"]][theta_idx]
    if("zeta" %in% names(eta_inv)) {
      par2=eta_inv[["zeta"]][theta_idx]
    } else {
      par2=rep(0,length(par1))
    }
  }

  pair_complete=pair_complete & is.finite(par1) & is.finite(par2)

  Fx_eval=Fx_1_2
  if(nrow(Fx_eval)>0) {
    Fx_eval[!is.finite(Fx_eval)]=0.5
    Fx_eval[Fx_eval>1]=1
    Fx_eval[Fx_eval<0]=0
  }

  par1_eval=par1
  par2_eval=par2
  par1_eval[!is.finite(par1_eval)]=0
  par2_eval[!is.finite(par2_eval)]=0

  if(copula_dist=="C") {
    par1_eval[par1_eval>=28]=27.9
  }

  copula_rect_prob <- NULL
  if(length(par1_eval)==0) {
    copula_d=numeric(0)
    copula_rect_prob <- numeric(0)
  } else if (identical(likelihood_type, "discrete_rectangle")) {
    Fx_lower_eval <- Fx_1_2_lower
    Fx_lower_eval[!is.finite(Fx_lower_eval)] <- 0.5
    Fx_lower_eval[Fx_lower_eval > 1] <- 1
    Fx_lower_eval[Fx_lower_eval < 0] <- 0
    copula_rect_prob <- .copula_rectangle_prob(
      Fx_eval[,1], Fx_eval[,2], Fx_lower_eval[,1], Fx_lower_eval[,2],
      family = copula_dist, par = par1_eval, par2 = par2_eval
    )
    denom <- margin_d[row_id1] * margin_d[row_id2]
    copula_d <- copula_rect_prob / denom
  } else {
    copula_d=.copula_pdf(Fx_eval[,1],Fx_eval[,2],family = copula_dist,par=par1_eval,par2=par2_eval)
  }
  if(length(copula_d)>0) {
    copula_d[!is.finite(copula_d) | copula_d<=0]=1
    copula_d[!pair_complete]=1
  }

  ########COMBINE MARGINS AND COPULA DERVIATIVES

  margin_loglik_terms=log(margin_d[!is.na(margin_d)])
  margin_loglik_terms=margin_loglik_terms[is.finite(margin_loglik_terms)]
  copula_loglik_terms=log(copula_d[pair_complete])
  copula_loglik_terms=copula_loglik_terms[is.finite(copula_loglik_terms)]

  log_lik=c(sum(margin_loglik_terms),sum(copula_loglik_terms),sum(margin_loglik_terms)+sum(copula_loglik_terms))
  names(log_lik)=c("marginal","copula","joint")

  copula_p=rep(NA_real_,length(copula_d))

  return_list=list(log_lik,margin_d,copula_d,margin_p,copula_p,Fx_1_2,order_copula,margin_deriv,pair_complete,par1,par2,row_id1,row_id2,pair_cache$theta_index_map,likelihood_type,margin_p_lower,Fx_1_2_lower,copula_rect_prob)
  names(return_list)=c("log_lik","margin_d","copula_d","margin_p","copula_p","Fx_1_2","order_copula","margin_deriv","pair_complete","copula_par1","copula_par2","copula_row_id1","copula_row_id2","copula_theta_index_map","likelihood_type","margin_p_lower","Fx_1_2_lower","copula_rect_prob")
  return(return_list)
}

.calc_likelihood_update_copula <- function(eta_inv, base_lik, copula_dist, pair_cache) {
  row_id1 <- pair_cache$row_id1
  row_id2 <- pair_cache$row_id2
  Fx_1_2 <- base_lik$Fx_1_2
  n_obs <- pair_cache$n_obs

  pair_complete <- pair_cache$observed_pair_base & is.finite(Fx_1_2[, 1]) & is.finite(Fx_1_2[, 2])
  if (identical(base_lik$likelihood_type, "discrete_rectangle")) {
    pair_complete <- pair_complete &
      is.finite(base_lik$Fx_1_2_lower[, 1]) & is.finite(base_lik$Fx_1_2_lower[, 2]) &
      is.finite(base_lik$margin_d[row_id1]) & is.finite(base_lik$margin_d[row_id2]) &
      base_lik$margin_d[row_id1] > 0 & base_lik$margin_d[row_id2] > 0
  }

  par1 <- rep(NA_real_, length(row_id1))
  par2 <- rep(NA_real_, length(row_id1))
  if(length(row_id1) > 0) {
    theta_len <- length(eta_inv[["theta"]])
    if(theta_len == n_obs) {
      theta_idx <- row_id1
    } else {
      theta_idx <- pair_cache$theta_index_map[row_id1]
    }

    par1 <- eta_inv[["theta"]][theta_idx]
    if("zeta" %in% names(eta_inv)) {
      par2 <- eta_inv[["zeta"]][theta_idx]
    } else {
      par2 <- rep(0, length(par1))
    }
  }

  pair_complete <- pair_complete & is.finite(par1) & is.finite(par2)

  Fx_eval <- Fx_1_2
  if(nrow(Fx_eval) > 0) {
    Fx_eval[!is.finite(Fx_eval)] <- 0.5
    Fx_eval[Fx_eval > 1] <- 1
    Fx_eval[Fx_eval < 0] <- 0
  }

  par1_eval <- par1
  par2_eval <- par2
  par1_eval[!is.finite(par1_eval)] <- 0
  par2_eval[!is.finite(par2_eval)] <- 0

  if(copula_dist == "C") {
    par1_eval[par1_eval >= 28] <- 27.9
  }

  copula_rect_prob <- base_lik$copula_rect_prob
  if(length(par1_eval) == 0) {
    copula_d <- numeric(0)
    copula_rect_prob <- numeric(0)
  } else if (identical(base_lik$likelihood_type, "discrete_rectangle")) {
    Fx_lower_eval <- base_lik$Fx_1_2_lower
    Fx_lower_eval[!is.finite(Fx_lower_eval)] <- 0.5
    Fx_lower_eval[Fx_lower_eval > 1] <- 1
    Fx_lower_eval[Fx_lower_eval < 0] <- 0
    copula_rect_prob <- .copula_rectangle_prob(
      Fx_eval[, 1], Fx_eval[, 2], Fx_lower_eval[, 1], Fx_lower_eval[, 2],
      family = copula_dist, par = par1_eval, par2 = par2_eval
    )
    copula_d <- copula_rect_prob / (base_lik$margin_d[row_id1] * base_lik$margin_d[row_id2])
  } else {
    copula_d <- .copula_pdf(Fx_eval[, 1], Fx_eval[, 2], family = copula_dist, par = par1_eval, par2 = par2_eval)
  }
  if(length(copula_d) > 0) {
    copula_d[!is.finite(copula_d) | copula_d <= 0] <- 1
    copula_d[!pair_complete] <- 1
  }

  copula_loglik_terms <- log(copula_d[pair_complete])
  copula_loglik_terms <- copula_loglik_terms[is.finite(copula_loglik_terms)]
  copula_loglik <- sum(copula_loglik_terms)
  marginal_loglik <- as.numeric(base_lik$log_lik["marginal"])
  log_lik <- c(marginal_loglik, copula_loglik, marginal_loglik + copula_loglik)
  names(log_lik) <- c("marginal", "copula", "joint")

  base_lik$log_lik <- log_lik
  base_lik$copula_d <- copula_d
  base_lik$copula_rect_prob <- copula_rect_prob
  base_lik$pair_complete <- pair_complete
  base_lik$copula_par1 <- par1
  base_lik$copula_par2 <- par2
  base_lik
}

#' @keywords internal
#' @noRd
score_function_v2 <- function(eta,dldpar,d2ldpar,dpardeta,response=NA,phi=1,step_size=1,verbose=FALSE,crit_wk=0.0000001) {

  u_k=dldeta = dldpar * dpardeta
  f_k=d2ldpar
  w_k=-f_k*(dpardeta*dpardeta)

  #Stop if weights are too small
  w_k[abs(w_k)<crit_wk]=1
  u_k[abs(w_k)<crit_wk]=0

  w_k[abs(u_k)<crit_wk]=1
  u_k[abs(u_k)<crit_wk]=0

  z_k=(1-phi)*eta+phi*(eta+step_size*(u_k/w_k))

  if(verbose==TRUE) {
    steps_mean=round(rbind(colMeans(as.matrix(eta))
                           ,colMeans(as.matrix(dldpar-dlcopdpar))
                           ,colMeans(as.matrix(dlcopdpar))
                           ,colMeans(as.matrix(dpardeta))
                           ,colMeans(as.matrix(dpardeta*dpardeta))
                           ,colMeans(as.matrix(f_k))
                           ,colMeans(as.matrix(w_k))
                           ,colMeans(as.matrix(u_k))
                           ,colMeans(as.matrix(u_k/w_k))
                           ,colMeans(as.matrix(z_k))
    ),8)
    rownames(steps_mean)=c("eta","dldpar","dlcopdpar","dpardeta","dpardeta2","f_k","w_k","u_k","(1/w_k)*u_k","z_k")
    print(steps_mean)
  }
  return_list=list(colMeans(as.matrix(z_k)),as.matrix(u_k),as.matrix(f_k),as.matrix(w_k),as.matrix(z_k))
  names(return_list)=c("par","u_k","f_k","w_k","z_k")
  return(return_list)
}

calc_deriv_copula_wrt_margin = function(input,margin_par,par_name,calc_d2=FALSE) {

    #Calculate copula derivative with respect to marginal parameters
    #input=copula_merged
    num_col <- function(nm) {
      val <- if (is.data.frame(input)) input[[nm]] else input[, nm]
      suppressWarnings(as.numeric(val))
    }

    if(calc_d2==FALSE) {

      dlcopdpar_row1=matrix(0,nrow=nrow(input),ncol=length(margin_par))
      dlcopdpar_row2=matrix(0,nrow=nrow(input),ncol=length(margin_par))
      i=1
      for (inner_par_name in margin_par) {

        if(inner_par_name==par_name) {
          #Take parameters from input for clarity
          dc_tplus_du_t=num_col("dcdu1")
          dc_tplus_du_tplus=num_col("dcdu2")
          l_t=num_col(paste(paste("dld",inner_par_name,sep=""),".x",sep=""))
          l_t_plus=num_col(paste(paste("dld",inner_par_name,sep=""),".y",sep=""))
          x_t=num_col("response.x")
          x_t_plus=num_col("response.y")
          f_t=num_col("margin_d.x")
          f_t_plus=num_col("margin_d.y")
          c_tplus=num_col("copula_d")
          mu_t=num_col("mu.x")
          mu_t_plus=num_col("mu.y")

          F_nd_t=num_col("F_nd.x")
          F_nd_t_plus=num_col("F_nd.y")

          du_t_dmu=F_nd_t
          du_t_plus_dmu=F_nd_t_plus

          # Exact endpoint attribution for pair log-copula derivative:
          # row_id1 gets (dc/du1)*(du1/dpar)/c, row_id2 gets (dc/du2)*(du2/dpar)/c.
          dlogc_row1=(dc_tplus_du_t * du_t_dmu) / c_tplus
          dlogc_row2=(dc_tplus_du_tplus * du_t_plus_dmu) / c_tplus
          dlogc_row1[!is.finite(dlogc_row1)] = 0
          dlogc_row2[!is.finite(dlogc_row2)] = 0

          dlcopdpar_row1[,i]=dlogc_row1
          dlcopdpar_row2[,i]=dlogc_row2

        }
        i=i+1
      }
      colnames(dlcopdpar_row1)=paste("dlcopd",margin_par,sep="")
      colnames(dlcopdpar_row2)=paste("dlcopd",margin_par,sep="")

      # Prefer explicit row-index accumulation to avoid merge-order instability.
      if(all(c("row_id1","row_id2") %in% colnames(input))) {
        n_obs <- suppressWarnings(max(c(num_col("row_id1"), num_col("row_id2")), na.rm = TRUE))
        if(!is.finite(n_obs) || n_obs < 1) {
          stop("Invalid row ids in copula-to-margin derivative assembly.")
        }
        d1_cop <- matrix(0, nrow = as.integer(n_obs), ncol = length(margin_par))
        colnames(d1_cop) <- margin_par

        row_id1 <- as.integer(num_col("row_id1"))
        row_id2 <- as.integer(num_col("row_id2"))
        for(j in seq_along(margin_par)) {
          contrib1 <- as.numeric(dlcopdpar_row1[, j])
          contrib2 <- as.numeric(dlcopdpar_row2[, j])
          valid1 <- is.finite(contrib1) & is.finite(row_id1) & row_id1 >= 1 & row_id1 <= n_obs
          valid2 <- is.finite(contrib2) & is.finite(row_id2) & row_id2 >= 1 & row_id2 <= n_obs
          d1_cop[row_id1[valid1], j] <- d1_cop[row_id1[valid1], j] + contrib1[valid1]
          d1_cop[row_id2[valid2], j] <- d1_cop[row_id2[valid2], j] + contrib2[valid2]
        }
        return(d1_cop)
      }

      dlcopdpar <- dlcopdpar_row1 + dlcopdpar_row2

      par_dlcopdpar=dlcopdpar[,paste("dlcopd",margin_par,sep="")]
      merged_dlcopdpar=merge(cbind(input[,c("time1","time2","subject1","subject2")],par_dlcopdpar),cbind(input[,c("time1","time2","subject1","subject2")],par_dlcopdpar),by.x=c("time2","subject2"),by.y=c("time1","subject1"),all=TRUE)
      merged_dlcopdpar[is.na(merged_dlcopdpar)]=0

      x_comp=grepl("dlcopd",colnames(merged_dlcopdpar))&grepl(".x",colnames(merged_dlcopdpar))
      y_comp=grepl("dlcopd",colnames(merged_dlcopdpar))&grepl(".y",colnames(merged_dlcopdpar))

      d1_cop=0.5*(merged_dlcopdpar[,x_comp]+merged_dlcopdpar[,y_comp])

      return(d1_cop)

    } else {

      d2lcopdpar2=matrix(0,nrow=nrow(input),ncol=length(margin_par))
      i=1
      for (inner_par_name in margin_par) {

        if(inner_par_name==par_name) {
          #Take parameters from input for clarity
          dc_tplus_du_t=num_col("dcdu1")
          dc_tplus_du_tplus=num_col("dcdu2")
          #l_t=input[,paste(paste("dld",inner_par_name,sep=""),".x",sep="")]
          #l_t_plus=input[,paste(paste("dld",inner_par_name,sep=""),".y",sep="")]
          #x_t=input[,"response.x"]
          #x_t_plus=input[,"response.y"]
          #f_t=input[,"margin_d.x"]
          #f_t_plus=input[,"margin_d.y"]
          c_tplus=num_col("copula_d")
          mu_t=num_col("mu.x")
          mu_t_plus=num_col("mu.y")

          F_nd_t=num_col("F_nd.x")
          F_nd_t_plus=num_col("F_nd.y")

          du_t_dmu=F_nd_t
          du_t_plus_dmu=F_nd_t_plus

          dc_plus_dt_dmu=dc_tplus_du_t * du_t_dmu
          dc_plus_dt_plus_dmu=dc_tplus_du_tplus * du_t_plus_dmu
          dc_plus_dt_dmu[is.nan(dc_plus_dt_dmu)]=0
          dc_plus_dt_plus_dmu[is.nan(dc_plus_dt_plus_dmu)]=0
          dcdmu_tplus=((dc_plus_dt_dmu + dc_plus_dt_plus_dmu) / c_tplus)
          dcdmu_tplus[is.nan(dcdmu_tplus)|is.na(dcdmu_tplus)]=0

          #dlcopdpar[,i]=dcdmu_tplus

          #######NOW FOR SECOND DERIVATIVE OF COPULA TERM

          F_nd2=num_col("F_nd2.x")
          F_nd2_plus=num_col("F_nd2.y")

          d2u_t_dmu2=F_nd2
          d2u_t_plus_dmu2=F_nd2_plus

          d2cdu_t2=num_col("d2cdu12")
          d2cdu_t_plus2=num_col("d2cdu22")
          d2cdu_t2[is.nan(d2cdu_t2)]=0
          d2cdu_t_plus2[is.nan(d2cdu_t_plus2)]=0

          d2cdmu2=  d2cdu_t2*du_t_dmu^2 +
                    dc_tplus_du_t * d2u_t_dmu2 +
                    d2cdu_t_plus2*du_t_plus_dmu^2 +
                    dc_tplus_du_tplus * d2u_t_plus_dmu2

          d2lcdmu2=as.matrix((d2cdmu2*c_tplus-(dcdmu_tplus^2))/(c_tplus^2))
          d2lcdmu2=num_col("c_nd2")

          d2lcopdpar2[,i]=d2lcdmu2

        }
        i=i+1
      }
      colnames(d2lcopdpar2)=paste("d2lcopd",margin_par,sep="")

      par_d2lcopdpar=d2lcopdpar2[,paste("d2lcopd",margin_par,sep="")]
      merged_d2lcopdpar=merge(cbind(input[,c("time1","time2","subject1","subject2")],par_d2lcopdpar)
                              ,cbind(input[,c("time1","time2","subject1","subject2")],par_d2lcopdpar)
                              ,by.x=c("time2","subject2"),by.y=c("time1","subject1"),all=TRUE)
      merged_d2lcopdpar[is.na(merged_d2lcopdpar)]=0

      x_comp=grepl("d2lcopd",colnames(merged_d2lcopdpar))&grepl(".x",colnames(merged_d2lcopdpar))
      y_comp=grepl("d2lcopd",colnames(merged_d2lcopdpar))&grepl(".y",colnames(merged_d2lcopdpar))

      d2_cop=0.5*(merged_d2lcopdpar[,x_comp]+merged_d2lcopdpar[,y_comp])

      #plot(d2lcopdpar2[,paste("d2lcopd",par_name,sep="")],input[,"c_nd2"],main="d2",ylab="numerical")

      return(d2_cop)
    }


}

.calc_dlcopdpar_indexed <- function(
  row_id1,
  row_id2,
  dcdu1,
  dcdu2,
  copula_d,
  F_nd,
  n_obs,
  pair_complete = NULL
) {
  row_id1 <- as.integer(row_id1)
  row_id2 <- as.integer(row_id2)
  n_pair <- length(row_id1)

  if (length(row_id2) != n_pair || length(dcdu1) != n_pair || length(dcdu2) != n_pair ||
      length(copula_d) != n_pair) {
    stop("Copula derivative inputs have inconsistent pair lengths.", call. = FALSE)
  }
  if (length(F_nd) != n_obs) {
    stop("F derivative length does not match the number of observations.", call. = FALSE)
  }

  if (is.null(pair_complete)) {
    pair_complete <- rep(TRUE, n_pair)
  } else {
    pair_complete <- as.logical(pair_complete)
    if (length(pair_complete) != n_pair) {
      stop("pair_complete length does not match copula pair length.", call. = FALSE)
    }
  }

  dlogc_row1 <- (as.numeric(dcdu1) * as.numeric(F_nd[row_id1])) / as.numeric(copula_d)
  dlogc_row2 <- (as.numeric(dcdu2) * as.numeric(F_nd[row_id2])) / as.numeric(copula_d)
  dlogc_row1[!pair_complete | !is.finite(dlogc_row1)] <- 0
  dlogc_row2[!pair_complete | !is.finite(dlogc_row2)] <- 0

  out <- numeric(n_obs)
  valid1 <- is.finite(row_id1) & row_id1 >= 1L & row_id1 <= n_obs
  valid2 <- is.finite(row_id2) & row_id2 >= 1L & row_id2 <= n_obs
  if (any(valid1)) {
    sum1 <- rowsum(dlogc_row1[valid1], row_id1[valid1], reorder = FALSE)
    out[as.integer(rownames(sum1))] <- out[as.integer(rownames(sum1))] + sum1[, 1]
  }
  if (any(valid2)) {
    sum2 <- rowsum(dlogc_row2[valid2], row_id2[valid2], reorder = FALSE)
    out[as.integer(rownames(sum2))] <- out[as.integer(rownames(sum2))] + sum2[, 1]
  }
  out
}

check_dlcopdpar_gradient_margin_score <- function(
  eta,
  eta_inv,
  par_name,
  margin_dist,
  copula_dist,
  dataset,
  mm,
  pair_cache,
  d1,
  base_loglik,
  verbose=FALSE
) {
  if (!par_name %in% c("mu", "sigma", "nu", "tau")) {
    return(list(warned = FALSE, message = NULL))
  }

  eta_vec <- as.numeric(eta[[par_name]])
  score_vec <- as.numeric(d1)
  finite_idx <- which(is.finite(eta_vec) & is.finite(score_vec))
  if (length(finite_idx) == 0) {
    return(list(warned = FALSE, message = NULL))
  }

  probe_candidates <- unique(round(seq(1, length(finite_idx), length.out = min(3, length(finite_idx)))))
  probe_idx <- finite_idx[probe_candidates]
  eps <- 1e-6
  tolerance <- 1e-3
  diffs <- numeric(length(probe_idx))

  linkinv_fun <- margin_dist[[paste(par_name, ".linkinv", sep = "")]]
  if (is.null(linkinv_fun)) {
    return(list(warned = FALSE, message = NULL))
  }

  for (k in seq_along(probe_idx)) {
    idx <- probe_idx[k]
    eta_pert <- eta
    eta_pert[[par_name]][idx] <- eta_pert[[par_name]][idx] + eps

    eta_inv_pert <- eta_inv
    eta_inv_pert[[par_name]] <- linkinv_fun(eta_pert[[par_name]])

    lik_pert <- calc_likelihood_minimal(
      eta_inv_pert,
      mm = mm,
      margin_dist = margin_dist,
      copula_dist = copula_dist,
      calc_d2 = FALSE,
      response = dataset$response,
      response_margin = dataset$time,
      response_subject = dataset$subject,
      pair_cache = pair_cache,
      calc_margin_deriv = FALSE
    )$log_lik["joint"]

    finite_diff <- (lik_pert - base_loglik) / eps
    diffs[k] <- finite_diff - score_vec[idx]
  }

  max_abs_diff <- max(abs(diffs), na.rm = TRUE)
  if (!is.finite(max_abs_diff)) {
    return(list(warned = FALSE, message = NULL))
  }

  message_text <- paste0(
    "DLCOPDGRAD check for ", par_name,
    ": max abs(score - finite_diff) = ", signif(max_abs_diff, 4),
    " over ", length(probe_idx), " probe row(s)."
  )

  if (max_abs_diff > tolerance) {
    return(list(warned = TRUE, message = paste0(message_text, " Potential score mismatch detected.")))
  }

  if (!is.null(verbose) && verbose > 1) {
    message(message_text)
  }

  return(list(warned = FALSE, message = message_text))
}

calc_copula_derivatives = function(eta_inv, Fx_1_2, copula_dist, calc_d2=FALSE, calc_d2_marginal=FALSE, par1=NULL, par2=NULL, pair_complete=NULL) {

  if(is.null(par1)) {
    par1=eta_inv[["theta"]]
  }

  if(is.null(par2)) {
    if("zeta" %in% names(eta_inv)) {
      par2=eta_inv[["zeta"]]
    } else {
      par2=eta_inv[["theta"]]*0
    }
  }

  if(is.null(pair_complete)) {
    pair_complete=rep(TRUE,length(par1))
  }

  if(length(par1)==0) {
    if("zeta" %in% names(eta_inv)) {
      if(calc_d2==TRUE) {
        return(list(dldth=numeric(0),dcdth=numeric(0),dldz=numeric(0),dcdz=numeric(0),dcdu1=numeric(0),dcdu2=numeric(0),d2ldth2=numeric(0),d2ldz2=numeric(0),d2ldthdz=numeric(0),d2cdu12=numeric(0),d2cdu22=numeric(0)))
      }
      return(list(dldth=numeric(0),dcdth=numeric(0),dldz=numeric(0),dcdz=numeric(0),dcdu1=numeric(0),dcdu2=numeric(0)))
    }
    if(calc_d2==TRUE) {
      return(list(dldth=numeric(0),dcdth=numeric(0),dcdu1=numeric(0),dcdu2=numeric(0),d2ldth2=numeric(0),d2cdu12=numeric(0),d2cdu22=numeric(0)))
    }
    return(list(dldth=numeric(0),dcdth=numeric(0),dcdu1=numeric(0),dcdu2=numeric(0)))
  }

  Fx_eval=as.matrix(Fx_1_2)
  Fx_eval[!is.finite(Fx_eval)]=0.5
  Fx_eval[Fx_eval>1]=1
  Fx_eval[Fx_eval<0]=0

  par1_eval=par1
  par2_eval=par2
  par1_eval[!is.finite(par1_eval)]=0
  par2_eval[!is.finite(par2_eval)]=0

  if(copula_dist=="C") {
    par1_eval[par1_eval>=28]=27.9
  }

  copula_d=.copula_pdf(Fx_eval[,1],Fx_eval[,2],family = copula_dist,par=par1_eval,par2=par2_eval)
  copula_d[!is.finite(copula_d) | copula_d<=0]=1
  copula_d[!pair_complete]=1

  dldth=.copula_deriv(Fx_eval[,1],Fx_eval[,2],family = copula_dist,par=par1_eval,par2=par2_eval,deriv="par",log=TRUE)
  dcdth=copula_d*dldth

  if(calc_d2==TRUE) {
    d2cdth=.copula_deriv2(Fx_eval[,1],Fx_eval[,2],family = copula_dist,par=par1_eval,par2=par2_eval,deriv="par")
    d2ldth2=(1/(copula_d^2))*(copula_d*d2cdth-dcdth^2)
  }

  if("zeta" %in% names(eta_inv)) {
    dldz=.copula_deriv(Fx_eval[,1],Fx_eval[,2],family = copula_dist,par=par1_eval,par2=par2_eval,deriv="par2",log=TRUE)
    dcdz=copula_d*dldz

    if(calc_d2==TRUE) {
      d2cdz=.copula_deriv2(Fx_eval[,1],Fx_eval[,2],family = copula_dist,par=par1_eval,par2=par2_eval,deriv="par2")
      d2ldz2=(1/(copula_d^2))*(copula_d*d2cdz-dcdz^2)

      d2cdthdz=.copula_deriv2(Fx_eval[,1],Fx_eval[,2],family = copula_dist,par=par1_eval,par2=par2_eval,deriv="par1par2")
      d2ldthdz=(d2cdthdz*copula_d-dcdth*dcdz)/(copula_d^2)
    }

  }
  dcdu1=.copula_deriv(Fx_eval[,1],Fx_eval[,2],family = copula_dist,par=par1_eval,par2=par2_eval,deriv="u1",log=FALSE)
  dcdu2=.copula_deriv(Fx_eval[,1],Fx_eval[,2],family = copula_dist,par=par1_eval,par2=par2_eval,deriv="u2",log=FALSE)

  dldth[!is.finite(dldth)]=0; if(calc_d2==TRUE) {d2ldth2[!is.finite(d2ldth2)]=0  }
  dldth[!pair_complete]=0
  dcdth[!pair_complete]=0
  dcdu1[!pair_complete]=0
  dcdu2[!pair_complete]=0

  if(calc_d2==TRUE) {
    d2cdu12=.copula_deriv2(Fx_eval[,1],Fx_eval[,2],family = copula_dist,par=par1_eval,par2=par2_eval,deriv="u1")
    d2cdu22=.copula_deriv2(Fx_eval[,1],Fx_eval[,2],family = copula_dist,par=par1_eval,par2=par2_eval,deriv="u2")
    d2cdu12[!is.finite(d2cdu12)]=0
    d2cdu22[!is.finite(d2cdu22)]=0
    d2cdu12[!pair_complete]=0
    d2cdu22[!pair_complete]=0
  }

  if("zeta" %in% names(eta_inv)) {
    dldz[!is.finite(dldz)]=0
    dcdz[!is.finite(dcdz)]=0
    dldz[!pair_complete]=0
    dcdz[!pair_complete]=0
    if(calc_d2==TRUE) {
      d2ldz2[!is.finite(d2ldz2)]=0
      d2ldthdz[!is.finite(d2ldthdz)]=0
      d2ldz2[!pair_complete]=0
      d2ldthdz[!pair_complete]=0
    }
  }

  if(calc_d2==TRUE) {
    d2ldth2[!pair_complete]=0
  }

  ############# RETURN LIST

  if("zeta" %in% names(eta_inv)) {
    if(calc_d2==TRUE) {
      return_list=list(dldth,dcdth,dldz,dcdz,dcdu1,dcdu2,d2ldth2,d2ldz2,d2ldthdz, d2cdu12, d2cdu22)
      names(return_list)=c("dldth","dcdth","dldz","dcdz","dcdu1","dcdu2","d2ldth2","d2ldz2","d2ldthdz","d2cdu12","d2cdu22")
    } else {
      return_list=list(dldth,dcdth,dldz,dcdz,dcdu1,dcdu2)
      names(return_list)=c("dldth","dcdth","dldz","dcdz","dcdu1","dcdu2")
    }
  } else {
    if(calc_d2==TRUE) {
      return_list=list(dldth,dcdth,dcdu1,dcdu2,d2ldth2, d2cdu12, d2cdu22)
      names(return_list)=c("dldth","dcdth","dcdu1","dcdu2","d2ldth2","d2cdu12","d2cdu22")
    } else {
      return_list=list(dldth,dcdth,dcdu1,dcdu2)
      names(return_list)=c("dldth","dcdth","dcdu1","dcdu2")
    }
  }

  return(return_list)
}

calc_Fx_derivatives = function(eta_inv, mm, margin_dist,response,par_names=NULL) {
  # Allow callers to pass full model matrix object; we only need fixed-effect blocks.
  if (is.list(mm) && all(c("x", "s") %in% names(mm))) {
    mm = mm$x
  }

  nd_impact=nd_impact_m=nd_impact_c=rep(0,length(names(eta_inv)))
  nd_impact_F=list() ###### WE DON"T NEED TO BE CALCULATING THIS FOR ALL PARAMETERS EVERY TIME

  margin_par_names=names(eta_inv)[names(eta_inv) %in% c("mu","sigma","nu","tau")]
  if(!is.null(par_names)) {
    margin_par_names=intersect(margin_par_names, par_names)
  }

  for (eta_par_names_nd in margin_par_names) {
    adj_fac=.0001
    change=change_m=change_c=c(0,0)
    change_F=matrix(0,nrow=length(eta_inv[[1]]),ncol=2)
    i=1
    for (adj in c(-1*adj_fac,adj_fac)) {
      eta_inv_adj=eta_inv
      eta_inv_adj[[eta_par_names_nd]]=eta_inv_adj[[eta_par_names_nd]]+adj
      #change[i]=calc_likelihood_minimal(eta_inv_adj,mm,margin_dist,copula_dist)$log_lik["joint"]
      #change_m[i]=calc_likelihood_minimal(eta_inv_adj,mm,margin_dist,copula_dist)$log_lik["marginal"]
      #change_c[i]=calc_likelihood_minimal(eta_inv_adj,mm,margin_dist,copula_dist)$log_lik["copula"]
      change_F[,i]=calc_F_x(eta_inv_adj,mm,margin_dist,response)#calc_likelihood_minimal(eta_inv_adj,mm,margin_dist,copula_dist)$margin_p
      #print(change_F[,i]==calc_F_x(eta_inv_adj,mm,margin_dist))
      i=i+1
    }
    #nd_impact[eta_par_names_nd]=(change[2]-change[1])/(2*adj_fac)
    #nd_impact_m[eta_par_names_nd]=(change_m[2]-change_m[1])/(2*adj_fac)
    #nd_impact_c[eta_par_names_nd]=(change_c[2]-change_c[1])/(2*adj_fac)
    nd_impact_F[[eta_par_names_nd]]=(change_F[,2]-change_F[,1])/(2*adj_fac)
  }
  return(nd_impact_F)
}

calc_Fx2_derivatives = function(eta_inv, mm, margin_dist,response,testing=FALSE,response_margin=NA,response_subject=NA) {
  nd_impact=nd_impact_m=nd_impact_c=rep(0,length(names(eta_inv)))
  nd_impact_F=nd_impact_c=list() ###### WE DON"T NEED TO BE CALCULATING THIS FOR ALL PARAMETERS EVERY TIME

  margin_par_names=names(eta_inv)[names(eta_inv) %in% c("mu","sigma","nu","tau")]

  for (eta_par_names_nd in margin_par_names) {
    adj_fac=.0001
    change_F=matrix(0,nrow=length(eta_inv[[1]]),ncol=3)
    change_c=matrix(0,nrow=length(eta_inv[["theta"]]),ncol=3)
    i=1
    for (adj in c(-1*adj_fac,adj_fac)) {
      eta_inv_adj=eta_inv
      eta_inv_adj[[eta_par_names_nd]]=eta_inv_adj[[eta_par_names_nd]]+adj
      #change[i]=calc_likelihood_minimal(eta_inv_adj,mm,margin_dist,copula_dist)$log_lik["joint"]
      #change_m[i]=calc_likelihood_minimal(eta_inv_adj,mm,margin_dist,copula_dist)$log_lik["marginal"]
      if(testing==TRUE) {
        change_c[,i]=log(calc_likelihood_minimal(eta_inv_adj,mm,margin_dist,copula_dist,calc_d2=FALSE,response,response_margin,response_subject)$copula_d)
      }
      change_F[,i]=calc_F_x(eta_inv_adj,mm,margin_dist,response)#calc_likelihood_minimal(eta_inv_adj,mm,margin_dist,copula_dist)$margin_p
      #print(change_F[,i]==calc_F_x(eta_inv_adj,mm,margin_dist))
      i=i+1
    }
    change_F[,3]=calc_F_x(eta_inv,mm,margin_dist,response)
    if(testing==TRUE) {
      change_c[,3]=log(calc_likelihood_minimal(eta_inv,mm,margin_dist,copula_dist,calc_d2=FALSE,response,response_margin,response_subject)$copula_d)
    }
    #nd_impact[eta_par_names_nd]=(change[2]-change[1])/(2*adj_fac)
    #nd_impact_m[eta_par_names_nd]=(change_m[2]-change_m[1])/(2*adj_fac)
    nd_impact_c[[eta_par_names_nd]]=(change_c[,2]+change_c[,1]-2*change_c[,3])/(adj_fac^2)
    nd_impact_F[[eta_par_names_nd]]=(change_F[,2]+change_F[,1]-2*change_F[,3])/(adj_fac^2)
  }
  if(testing==FALSE) {
    return(nd_impact_F)
  } else {
    return(list(nd_impact_F,nd_impact_c))
  }
}

calc_true_SE_numderiv_only = function(eta_inv, mm, margin_dist,response,testing=FALSE,response_margin=NA,response_subject=NA) {

  adj_fac=.001
  nd_impact=rep(0,length(names(eta_inv)))
  names(nd_impact)=margin_par_names=names(eta_inv)#[names(eta_inv) %in% c("mu","sigma","nu","tau")]

  for (eta_par_names_nd in margin_par_names) {

    change=rep(0,length(names(eta_inv)))
    i=1
    for (adj in c(-1*adj_fac,adj_fac)) {
      eta_inv_adj=eta_inv
      eta_inv_adj[[eta_par_names_nd]]=eta_inv_adj[[eta_par_names_nd]]+adj
      change[i]=calc_likelihood_minimal(eta_inv_adj,mm,margin_dist,copula_dist,calc_d2 = FALSE,response,response_margin,response_subject)$log_lik["joint"]
      i=i+1
    }
    change[3]=calc_likelihood_minimal(eta_inv,mm,margin_dist,copula_dist,calc_d2 = FALSE,response,response_margin,response_subject)$log_lik["joint"]
    nd_impact[eta_par_names_nd]=(change[2]+change[1]-2*change[3])/(adj_fac^2)

    #print(c(change,nd_impact[eta_par_names_nd]))
  }

  nd_cross=matrix(0,nrow=length(names(eta_inv)),ncol=length(names(eta_inv)))
  colnames(nd_cross)=rownames(nd_cross)=names(eta_inv)
  for (name1 in margin_par_names) {
    for (name2 in margin_par_names) {
      for(adj1 in c(-1*adj_fac,adj_fac)) {
        for(adj2 in c(-1*adj_fac,adj_fac)) {
          if(name1!=name2) {
            eta_inv_adj=eta_inv
            eta_inv_adj[[name1]]=eta_inv_adj[[name1]]+adj1
            eta_inv_adj[[name2]]=eta_inv_adj[[name2]]+adj2
            change=calc_likelihood_minimal(eta_inv_adj,mm,margin_dist,copula_dist,calc_d2 = FALSE,response,response_margin,response_subject)$log_lik["joint"]
            nd_cross[name1,name2]=nd_cross[name1,name2]+change*if(adj1==adj2){1} else {-1}
          }
        }
      }
    }
  }
  nd_cross=nd_cross/(4*(adj_fac^2))

  nd2=(diag(nd_impact)+nd_cross)

  return(nd2)
}

calc_true_SE_numderiv_only_rowwise = function(eta_inv, mm, margin_dist,response,testing=FALSE,response_margin=NA,response_subject=NA) {

  adj_fac=.00001
  nd_impact_m=nd_impact_c=list()

  for (eta_par_names_nd in margin_par_names) {

    change_m=change_c=list()
    i=1
    for (adj in c(-1*adj_fac,adj_fac)) {
      eta_inv_adj=eta_inv
      eta_inv_adj[[eta_par_names_nd]]=eta_inv_adj[[eta_par_names_nd]]+adj

      calc_lik_temp=calc_likelihood_minimal(eta_inv_adj,mm,margin_dist,copula_dist,calc_d2 = FALSE,response,response_margin,response_subject)

      change_m[[i]]=calc_lik_temp$margin_d
      change_c[[i]]=calc_lik_temp$copula_d
      i=i+1
    }
    calc_lik_temp=calc_likelihood_minimal(eta_inv,mm,margin_dist,copula_dist,calc_d2 = FALSE,response,response_margin,response_subject)
    change_m[[3]]=calc_lik_temp$margin_d
    change_c[[3]]=calc_lik_temp$copula_d
    nd_impact_m[[eta_par_names_nd]]=(change_m[[3]]+change_m[[1]]-2*change_m[[3]])/(adj_fac^2)
    nd_impact_c[[eta_par_names_nd]]=(change_c[[3]]+change_c[[1]]-2*change_c[[3]])/(adj_fac^2)
  }
  names(nd_impact_m)=names(nd_impact_c)=margin_par_names=names(eta_inv)

  nd_cross_m
  colnames(nd_cross)=rownames(nd_cross)=names(eta_inv)
  for (name1 in margin_par_names) {
    for (name2 in margin_par_names) {
      for(adj1 in c(-1*adj_fac,adj_fac)) {
        for(adj2 in c(-1*adj_fac,adj_fac)) {
          if(name1!=name2) {
            eta_inv_adj=eta_inv
            eta_inv_adj[[name1]]=eta_inv_adj[[name1]]+adj1
            eta_inv_adj[[name2]]=eta_inv_adj[[name2]]+adj2
            change=calc_likelihood_minimal(eta_inv_adj,mm,margin_dist,copula_dist,calc_d2 = FALSE,response,response_margin,response_subject)$log_lik["joint"]
            nd_cross[name1,name2]=nd_cross[name1,name2]+change*if(adj1==adj2){1} else {-1}
          }
        }
      }
    }
  }
  nd_cross=nd_cross/(4*(adj_fac^2))

  nd2=(diag(nd_impact)+nd_cross)

  return(nd2)
}

calc_true_SE_numderiv_only_covariates = function(object, par, mm, margin_dist,response,testing=FALSE,response_margin=NA,response_subject=NA,h=.0001, progress=interactive()) {

  #object=fit; par=NA; numderiv=TRUE; sep_d2=TRUE

  response=object$response
  response_margin=object$response_margin
  response_subject=object$response_subject

  #margin_names=unique(object$response_margin)
  #num_margins=length(margin_names)

  #se_out=object$par*0;
  margin_dist=object$margin_dist; copula_dist=object$copula_dist; copula_link=get_copula_dist(copula_dist)$copula_link
  mm=object$model_matrix

  par_cov=object$par
  par_s=object$par_s

  input_par=par_cov
  progress=isTRUE(progress)

  adj_fac=h
  par_names=names(input_par)
  nd_impact=rep(0,length(par_names))
  names(nd_impact)=par_names#[names(eta_inv) %in% c("mu","sigma","nu","tau")]

  # Reuse fixed copula pairing metadata across all numerical derivative evaluations.
  pair_cache=build_copula_pair_cache(response,response_margin,response_subject)
  eval_joint_loglik <- function(par_vec) {
    eta_out=calc_eta(par_cov=par_vec,mm=mm,margin_dist,copula_link,par_s)
    eta_inv=eta_out$eta_inv
    calc_likelihood_minimal(
      eta_inv,
      mm=mm$x,
      margin_dist,
      copula_dist,
      calc_d2=FALSE,
      response=response,
      response_margin=response_margin,
      response_subject=response_subject,
      pair_cache=pair_cache,
      calc_margin_deriv=FALSE
    )$log_lik["joint"]
  }

  base_loglik=eval_joint_loglik(input_par)

  if(progress) cat("Calculating numerical first derivates for Hessian matrix...\n")
  pb_first <- NULL
  if(progress) {
    pb_first <- utils::txtProgressBar(min=0, max=length(par_names), style=3)
    on.exit(close(pb_first), add=TRUE)
  }
  first_counter <- 0L
  for (par_names_nd in par_names) {

    #print(par_names_nd)
    change=rep(0,3)
    i=1
    for (adj in c(-1*adj_fac,adj_fac)) {

      par_cov=input_par
      par_cov[[par_names_nd]]=par_cov[[par_names_nd]]+adj

      change[i]=eval_joint_loglik(par_cov)
      i=i+1
    }
    change[3]=base_loglik
    nd_impact[par_names_nd]=(change[2]+change[1]-2*change[3])/(adj_fac^2)

    first_counter <- first_counter + 1L
    if(progress) {
      utils::setTxtProgressBar(pb_first, first_counter)
    }

    #print(c(change,nd_impact[eta_par_names_nd]))
  }
  if(progress) cat("\n")

  if(progress) cat("Calculating numerical second derivates for Hessian matrix... this may take a while\n")
  p=length(par_names)
  nd_cross=matrix(0,nrow=p,ncol=p)
  colnames(nd_cross)=rownames(nd_cross)=par_names
  second_total <- if(p > 1) choose(p,2) else 0
  pb_second <- NULL
  if(progress && second_total > 0) {
    pb_second <- utils::txtProgressBar(min=0, max=second_total, style=3)
    on.exit(close(pb_second), add=TRUE)
  }
  second_counter <- 0L
  if(p > 1) {
    for (i in 1:(p-1)) {
      name1=par_names[i]
      for (j in (i+1):p) {
        name2=par_names[j]
        cross_sum=0
        for(adj1 in c(-1*adj_fac,adj_fac)) {
          for(adj2 in c(-1*adj_fac,adj_fac)) {
            par=input_par
            par[[name1]]=par[[name1]]+adj1
            par[[name2]]=par[[name2]]+adj2

            change=eval_joint_loglik(par)
            cross_sum=cross_sum+change*if(adj1==adj2){1} else {-1}
          }
        }
        nd_cross[name1,name2]=cross_sum
        nd_cross[name2,name1]=cross_sum

        second_counter <- second_counter + 1L
        if(progress && !is.null(pb_second)) {
          utils::setTxtProgressBar(pb_second, second_counter)
        }
      }
    }
  }
  if(progress && !is.null(pb_second)) cat("\n")
  nd_cross=nd_cross/(4*(adj_fac^2))

  nd2=(diag(nd_impact)+nd_cross)

  return(nd2)
}

#' Log-likelihood for a fitted longitudinal GAMLSS-copula model
#'
#' @param object A fitted `gamlss.longitudinal` object.
#' @param ... Additional arguments, currently unused.
#'
#' @return Named log-likelihood components.
#' @export
