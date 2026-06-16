#' Validate imported fitting controls, and update to be in the right format if needed or error.
#' 
#' This function is called internally by `fit_model()` to validate and transform the controls for model fitting. 
#' It checks that inputted controls are the correct type, in acceptable ranges, and consistent with each other. 
#' It also sets defaults for any controls that were not explicitly provided by the user. 
#' The final controls are returned as a list that can be used in the fitting process.
#'
#' @noRd
.gl_normalize_fit_controls <- function(
    method,
    start_from,
    warm_start_joint,
    warm_start_joint_iter,
    backtracking_max_halves,
    cg_max_stall,
    cg_max_delta,
    cg_lambda_update_every,
    cg_max_lambda_updates,
    cg_raw_loglik_drop_tol,
    cg_line_search,
    cg_max_line_search_evals,
    cg_gradient_method,
    discrete_score_method,
    cg_zeta_hessian,
    cg_hessian_method,
    vcov_method,
    vcov_numderiv,
    rs_smooth_trust_radius) {
  backtracking_max_halves <- .gl_normalize_backtracking_halves(backtracking_max_halves)
  rs_smooth_trust_radius <- .gl_validate_rs_smooth_trust_radius(rs_smooth_trust_radius)

  method <- toupper(as.character(method)[1])
  cg_line_search <- match.arg(as.character(cg_line_search)[1], c("first", "best"))
  cg_gradient_method <- match.arg(as.character(cg_gradient_method)[1], c("analytical", "forward", "central"))
  discrete_score_method <- match.arg(as.character(discrete_score_method)[1], c("analytical", "finite"))
  cg_zeta_hessian <- match.arg(as.character(cg_zeta_hessian)[1], c("analytical", "finite"))
  cg_hessian_method <- match.arg(as.character(cg_hessian_method)[1], c("analytical", "finite", "auto"))

  cg_max_line_search_evals <- .gl_normalize_cg_line_search_evals(cg_max_line_search_evals)

  if (!method %in% c("RS", "CG")) {
    stop("ERROR: method must be one of 'RS' or 'CG'.")
  }

  user_supplied_start <- !all(is.na(start_from))

  warm_start_controls <- .gl_normalize_warm_start_controls(warm_start_joint, warm_start_joint_iter)
  warm_start_joint <- warm_start_controls$warm_start_joint
  warm_start_joint_iter <- warm_start_controls$warm_start_joint_iter

  vcov_controls <- .gl_normalize_vcov_controls(vcov_method, vcov_numderiv)
  vcov_method <- vcov_controls$vcov_method
  vcov_numderiv <- vcov_controls$vcov_numderiv

  cg_lambda_controls <- .gl_normalize_cg_lambda_controls(
    cg_lambda_update_every = cg_lambda_update_every,
    cg_max_lambda_updates = cg_max_lambda_updates,
    cg_raw_loglik_drop_tol = cg_raw_loglik_drop_tol
  )
  cg_lambda_update_every <- cg_lambda_controls$cg_lambda_update_every
  cg_max_lambda_updates <- cg_lambda_controls$cg_max_lambda_updates
  cg_raw_loglik_drop_tol <- cg_lambda_controls$cg_raw_loglik_drop_tol

  cg_fallback_controls <- .gl_normalize_cg_fallback_controls(cg_max_stall, cg_max_delta)
  cg_max_stall <- cg_fallback_controls$cg_max_stall
  cg_max_delta <- cg_fallback_controls$cg_max_delta

  list(
    method = method,
    user_supplied_start = user_supplied_start,
    warm_start_joint = warm_start_joint,
    warm_start_joint_iter = warm_start_joint_iter,
    backtracking_max_halves = backtracking_max_halves,
    cg_max_stall = cg_max_stall,
    cg_max_delta = cg_max_delta,
    cg_lambda_update_every = cg_lambda_update_every,
    cg_max_lambda_updates = cg_max_lambda_updates,
    cg_raw_loglik_drop_tol = cg_raw_loglik_drop_tol,
    cg_line_search = cg_line_search,
    cg_max_line_search_evals = cg_max_line_search_evals,
    cg_gradient_method = cg_gradient_method,
    discrete_score_method = discrete_score_method,
    cg_zeta_hessian = cg_zeta_hessian,
    cg_hessian_method = cg_hessian_method,
    vcov_method = vcov_method,
    vcov_numderiv = vcov_numderiv,
    rs_smooth_trust_radius = rs_smooth_trust_radius
  )
}
