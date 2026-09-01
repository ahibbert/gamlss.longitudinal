#' Compute analytical (semi-analytical) Hessian for gamlss.longitudinal
#'
#' This is the main entry point called by vcov.gamlss.longitudinal when
#' method = "analytical".
#'
#' @param object A fitted gamlss.longitudinal object.
#' @param progress Logical; show a progress indicator.
#' @param h Step size for CDF finite differences (default 1e-4).
#'
#' @return A square matrix of dimension length(object$par) x length(object$par)
#'   representing the Hessian of the joint log-likelihood evaluated at the MLE.
#'
#' @keywords internal
#' @noRd
calc_analytical_hessian <- function(object, progress = interactive(), h = 1e-4) {
  progress <- isTRUE(progress)

  response <- object$response

  response_margin <- object$response_margin

  response_subject <- object$response_subject

  margin_dist <- object$margin_dist

  copula_dist <- object$copula_dist

  copula_link <- get_copula_dist(copula_dist)$copula_link

  mm <- object$model_matrix

  par_cov <- object$par

  par_s <- object$par_s

  if (progress) cat("Analytical Hessian: computing eta and likelihood components...\n")

  eta_out <- calc_eta(par_cov, mm, margin_dist, copula_link, par_s = par_s)

  eta_inv <- eta_out$eta_inv

  eta_dr <- eta_out$eta_dr

  eta_d2 <- .calc_eta_d2_linkinv(eta_out$eta, margin_dist, copula_link)

  .warn_gg_near_zero_nu(margin_dist, eta_inv)

  pair_cache <- build_copula_pair_cache(response, response_margin, response_subject)

  # Full likelihood call to get margin_deriv (analytical d1 and d2 from gamlss family)

  calc_lik <- calc_likelihood_minimal(

    eta_inv,
    mm = mm$x, margin_dist, copula_dist,
    calc_d2 = TRUE,
    response = response, response_margin = response_margin,
    response_subject = response_subject,
    pair_cache = pair_cache
  )

  if (identical(calc_lik$likelihood_type, "discrete_rectangle")) {
    if (progress) cat("Analytical Hessian: assembling discrete rectangle contributions...\n")

    margin_d1l <- .get_margin_d1l(calc_lik$margin_deriv, mm)

    margin_d2l <- .calc_margin_d2l_fd(eta_inv, mm, margin_dist, response, h = h)

    cop_hess <- .calc_discrete_rectangle_hessian_contributions(
      eta_inv = eta_inv,
      eta_dr = eta_dr,
      eta_d2 = eta_d2,
      pair_cache = pair_cache,
      margin_dist = margin_dist,
      copula_dist = copula_dist,
      response = response,
      calc_lik = calc_lik,
      mm = mm,
      h = h
    )

    H <- .assemble_covariate_hessian(object, margin_d1l, margin_d2l, cop_hess, eta_dr, eta_d2, mm, pair_cache)

    if (progress) cat("Analytical Hessian: done.\n")

    return(H)
  }

  # Attach margin_p to eta_inv for use inside copula hessian function

  eta_inv[["margin_p_cache"]] <- calc_lik$margin_p

  if (progress) cat("Analytical Hessian: computing CDF derivatives...\n")

  # Evaluate only on observed (non-NA) response rows; set to 0 elsewhere

  dF <- .calc_dFdpar(eta_inv, mm, margin_dist, response, h = h)

  d2F <- .calc_d2Fdpar2(eta_inv, mm, margin_dist, response, h = h)

  d2F_x <- .calc_d2Fdpar_cross(eta_inv, mm, margin_dist, response, h = h)

  if (progress) cat("Analytical Hessian: assembling copula Hessian contributions...\n")

  cop_hess <- .calc_copula_hessian_contributions(
    eta_inv, pair_cache, copula_dist, dF, d2F, d2F_x
  )

  if (progress) cat("Analytical Hessian: extracting margin second derivatives...\n")

  margin_d1l <- .get_margin_d1l(calc_lik$margin_deriv, mm)

  margin_d2l <- .calc_margin_d2l_fd(eta_inv, mm, margin_dist, response, h = h)

  if (progress) cat("Analytical Hessian: assembling covariate Hessian matrix...\n")

  H <- .assemble_covariate_hessian(object, margin_d1l, margin_d2l, cop_hess, eta_dr, eta_d2, mm, pair_cache)

  if (progress) cat("Analytical Hessian: done.\n")

  H
}
