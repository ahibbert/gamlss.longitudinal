#' Variance-covariance matrix for a fitted longitudinal GAMLSS-copula model
#'
#' @param object A fitted `gamlss.longitudinal` object.
#' @param par Optional parameter list for evaluating uncertainty away from the
#'   fitted coefficients.
#' @param sep_d2 Logical legacy argument retained for compatibility.
#' @param numderiv Logical; use the numerical Hessian path.
#' @param method Character; variance-covariance method to use. `"analytical"` (default)
#'   uses the semi-analytical Hessian from `R/hessian-analytical.R`.
#'   `"numderiv"` uses full finite-difference numerical second derivatives as a
#'   slower reference path. `"sandwich"` uses a cluster-robust sandwich estimator
#'   for fixed coefficients, with clusters defaulting to subjects. The legacy
#'   `numderiv` logical argument is still accepted and maps to `method = "numderiv"`
#'   when `TRUE`.
#' @param progress Logical; show progress bars for slow Hessian calculations.
#' @param h Numeric finite-difference step used by the analytical Hessian
#'   helper.
#' @param cluster Optional cluster labels for `method = "sandwich"`. Defaults to
#'   the fitted subject identifiers.
#' @param sandwich_h Numeric finite-difference step used for sandwich cluster
#'   score contributions.
#' @param sandwich_adjust Logical; apply a finite-sample cluster correction to
#'   sandwich covariance estimates.
#' @param sandwich_bread_method Character; Hessian method used for the sandwich
#'   bread.
#' @param ... Additional arguments, currently unused.
#'
#' @return A list containing variance-covariance matrices and standard errors.
#' @export
vcov.gamlss.longitudinal <- function(object, par = NA, sep_d2 = TRUE, numderiv = FALSE,
                                     method = c("analytical", "numderiv", "analytical_only", "sandwich"),
                                     progress = interactive(), h = 1e-4,
                                     cluster = NULL,
                                     sandwich_h = 1e-5,
                                     sandwich_adjust = TRUE,
                                     sandwich_bread_method = c("analytical", "numderiv", "analytical_only"),
                                     ...) {
  # object=fit; par=NA; numderiv=TRUE; sep_d2=TRUE

  vcov_setup <- .gl_prepare_vcov_evaluation(
    object = object,
    par = par,
    numderiv = numderiv,
    method = method,
    progress = progress
  )
  method <- vcov_setup$method
  method_requested <- vcov_setup$method_requested
  method_used <- vcov_setup$method_used
  progress <- vcov_setup$progress
  include_dlcopdpar <- vcov_setup$include_dlcopdpar
  response <- vcov_setup$response
  response_margin <- vcov_setup$response_margin
  response_subject <- vcov_setup$response_subject
  margin_names <- vcov_setup$margin_names
  num_margins <- vcov_setup$num_margins
  margin_dist <- vcov_setup$margin_dist
  copula_dist <- vcov_setup$copula_dist
  copula_link <- vcov_setup$copula_link
  mm <- vcov_setup$mm
  par_cov <- vcov_setup$par_cov
  par_s <- vcov_setup$par_s
  eta_out <- vcov_setup$eta_out
  eta_inv <- vcov_setup$eta_inv
  eta_dr <- vcov_setup$eta_dr
  eta <- vcov_setup$eta

  if (identical(method, "sandwich")) {
    sandwich_path <- .gl_vcov_compute_sandwich(
      object = object,
      par_cov = par_cov,
      par_s = par_s,
      mm = mm,
      margin_dist = margin_dist,
      response = response,
      response_margin = response_margin,
      response_subject = response_subject,
      cluster = cluster,
      score_h = sandwich_h,
      bread_h = h,
      adjust = sandwich_adjust,
      bread_method = sandwich_bread_method,
      progress = progress
    )

    return(.gl_vcov_build_result(
      object = object,
      eta_inv = eta_inv,
      response = response,
      vcov_final = sandwich_path$vcov_final,
      se_final = sandwich_path$se_final,
      method_used = sandwich_path$method_used,
      method_requested = method_requested,
      hessian_diagnostics = sandwich_path$hessian_diagnostics
    ))
  }

  method_info <- .gl_vcov_apply_margin_preflight(method, method_used, margin_dist, eta_inv)
  method <- method_info$method
  method_used <- method_info$method_used

  # if(!all(is.na(par))) {response=eta_inv[["mu"]]}

  calc_lik_out <- calc_likelihood_minimal(eta_inv,
    mm = mm$x, margin_dist, copula_dist, calc_d2 = TRUE,
    response = response, response_margin = response_margin, response_subject = response_subject
  )

  Fx_1_2 <- calc_lik_out$Fx_1_2
  margin_p <- calc_lik_out$margin_p
  margin_d <- calc_lik_out$margin_d
  copula_d <- calc_lik_out$copula_d

  method_info <- .gl_vcov_apply_likelihood_preflight(method, method_used, calc_lik_out, response)
  method <- method_info$method
  method_used <- method_info$method_used

  ### Calculate derivaties: margin and copula d1 and d2

  margin_derivatives <- calc_lik_out$margin_deriv

  copula_derivatives <- calc_copula_derivatives(eta_inv, Fx_1_2, copula_dist, calc_d2 = TRUE)

  # Calculate numerical derivatives for F(x) - also used as a reference for overall likelihood derivative

  nd_impact_F <- calc_Fx_derivatives(eta_inv, mm$x, margin_dist, response)

  nd_impact_F2 <- calc_Fx2_derivatives(eta_inv, mm$x, margin_dist, response)

  vcov_path <- .gl_vcov_compute_primary(
    object = object,
    par_cov = par_cov,
    mm = mm,
    margin_dist = margin_dist,
    response = response,
    response_margin = response_margin,
    response_subject = response_subject,
    method = method,
    progress = progress,
    h = h
  )
  method_used <- vcov_path$method_used

  d2_mat <- NULL
  if (!(method %in% c("numderiv", "analytical", "analytical_only"))) {
    d2_mat <- .gl_vcov_legacy_derivative_matrix(
      object = object,
      eta_inv = eta_inv,
      mm = mm,
      margin_dist = margin_dist,
      response = response,
      response_margin = response_margin,
      response_subject = response_subject,
      margin_names = margin_names,
      num_margins = num_margins,
      margin_p = margin_p,
      margin_d = margin_d,
      copula_d = copula_d,
      copula_derivatives = copula_derivatives,
      margin_derivatives = margin_derivatives,
      nd_impact_F = nd_impact_F,
      nd_impact_F2 = nd_impact_F2,
      eta = eta,
      include_dlcopdpar = include_dlcopdpar,
      sep_d2 = sep_d2
    )
  }

  vcov_solved <- .gl_vcov_solve_if_needed(vcov_path, method, d2_mat, response)

  return(.gl_vcov_build_result(
    object = object,
    eta_inv = eta_inv,
    response = response,
    vcov_final = vcov_solved$vcov_final,
    se_final = vcov_solved$se_final,
    method_used = method_used,
    method_requested = method_requested,
    hessian_diagnostics = vcov_solved$hessian_diagnostics
  ))
}
