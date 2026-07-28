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
    h = 1e-5) {
  method <- match.arg(method)
  if (is.null(pair_cache)) {
    pair_cache <- build_copula_pair_cache(response, response_margin, response_subject)
  }
  if (is.null(calc_lik)) {
    calc_lik <- calc_likelihood_minimal(
      eta_inv,
      mm = mm, margin_dist, copula_dist, response = response,
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
          plus,
          mm = mm, margin_dist, copula_dist, response = response,
          response_margin = response_margin, response_subject = response_subject,
          pair_cache = pair_cache, calc_margin_deriv = FALSE
        )$log_lik["joint"]
        lm <- calc_likelihood_minimal(
          minus,
          mm = mm, margin_dist, copula_dist, response = response,
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
