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

.cg_copula_parameter_gradient <- function(
    par_name,
    derivative,
    eta,
    eta_dr,
    mm_cg,
    calc_lik,
    response) {
  n_par <- length(eta[[par_name]])
  d1_full <- rep(0, n_par)
  row_id1 <- calc_lik$copula_row_id1
  if (length(row_id1) > 0) {
    if (n_par == length(response)) {
      par_idx <- row_id1
    } else {
      par_idx <- calc_lik$copula_theta_index_map[row_id1]
    }
    valid_idx <- is.finite(par_idx) & par_idx >= 1 & par_idx <= n_par
    if (any(valid_idx)) {
      d1_sum <- rowsum(derivative[valid_idx], par_idx[valid_idx], reorder = FALSE)
      d1_full[as.integer(rownames(d1_sum))] <- d1_sum[, 1]
    }
  }
  score_eta <- d1_full * as.numeric(eta_dr[[par_name]])
  score_eta[!is.finite(score_eta)] <- 0
  par_grad <- as.numeric(crossprod(as.matrix(mm_cg$x[[par_name]]), score_eta))
  x_names <- colnames(mm_cg$x[[par_name]])
  names(par_grad) <- ifelse(
    startsWith(x_names, paste0(par_name, ".")),
    x_names,
    paste(par_name, x_names, sep = ".")
  )
  par_grad
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
    response_subject) {
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

  margin_score_natural <- .cg_margin_natural_scores(
    margin_par = margin_par,
    eta = eta,
    eta_inv = eta_inv,
    mm_cg = mm_cg,
    calc_lik = calc_lik,
    margin_dist = margin_dist,
    copula_derivatives = copula_derivatives,
    include_dlcopdpar = include_dlcopdpar,
    response = response,
    response_margin = response_margin,
    response_subject = response_subject
  )

  for (par_name in margin_par) {
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

  if ("theta" %in% copula_par) {
    par_grad <- .cg_copula_parameter_gradient(
      "theta",
      copula_derivatives$dldth,
      eta,
      eta_dr,
      mm_cg,
      calc_lik,
      response
    )
    idx <- match(names(par_grad), names(grad))
    valid <- !is.na(idx)
    grad[idx[valid]] <- par_grad[valid]
  }

  if ("zeta" %in% copula_par && "dldz" %in% names(copula_derivatives)) {
    par_grad <- .cg_copula_parameter_gradient(
      "zeta",
      copula_derivatives$dldz,
      eta,
      eta_dr,
      mm_cg,
      calc_lik,
      response
    )
    idx <- match(names(par_grad), names(grad))
    valid <- !is.na(idx)
    grad[idx[valid]] <- par_grad[valid]
  }

  grad
}
