#' Compute the effective degrees of freedom for one RS smooth
#'
#' @noRd
.gl_rs_smooth_term_edf <- function(B, weights, penalty, lambda) {
  if (!is.matrix(B) || !nrow(B) || !ncol(B)) {
    stop("RS smooth EDF requires a non-empty smooth design matrix.", call. = FALSE)
  }
  weights <- as.numeric(weights)
  if (length(weights) != nrow(B) || any(!is.finite(weights))) {
    stop("RS smooth EDF requires one finite working weight per design row.", call. = FALSE)
  }
  if (!is.matrix(penalty) || !identical(dim(penalty), c(ncol(B), ncol(B))) ||
      any(!is.finite(penalty))) {
    stop("RS smooth EDF requires a finite square penalty matrix matching the smooth basis.", call. = FALSE)
  }
  lambda <- as.numeric(lambda)[1L]
  if (!is.finite(lambda) || lambda < 0) {
    stop("RS smooth EDF requires a finite non-negative smoothing parameter.", call. = FALSE)
  }

  information <- crossprod(B, B * weights)
  edf <- sum(.solve_linear_system(information + lambda * penalty) * information)
  if (!is.finite(edf) || edf < -sqrt(.Machine$double.eps)) {
    stop("RS smooth EDF calculation returned an invalid value.", call. = FALSE)
  }
  max(0, as.numeric(edf))
}

#' Refresh all RS smooth EDFs from the final optimizer state
#'
#' An RS block can legitimately require no accepted coefficient update in the
#' last fit. Its working weights and fitted smooth remain valid, but the
#' iteration-local EDF value may then still be the zero used at initialization.
#' Information criteria must instead use EDFs recomputed from the final stored
#' weights, penalties, and smoothing parameters for every smooth block.
#'
#' @noRd
.gl_refresh_rs_smooth_edf <- function(df_s, mm_s, lambda_s, weights_final) {
  if (is.null(mm_s) || !length(mm_s)) return(df_s)
  if (is.null(df_s)) df_s <- list()

  for (par_name in names(mm_s)) {
    smooths <- mm_s[[par_name]]
    if (is.null(smooths) || !length(smooths)) next
    if (is.null(weights_final[[par_name]])) {
      stop("Final RS working weights are missing for smooth parameter '", par_name, "'.", call. = FALSE)
    }
    if (is.null(df_s[[par_name]])) df_s[[par_name]] <- list()

    for (smooth_name in names(smooths)) {
      B <- smooths[[smooth_name]]
      penalty <- attr(B, "penalty")
      lambda <- lambda_s[[par_name]][[smooth_name]]
      if (is.null(lambda)) {
        stop(
          "Final RS smoothing parameter is missing for '", par_name, ": ",
          smooth_name, "'.", call. = FALSE
        )
      }
      value <- .gl_rs_smooth_term_edf(
        B = B,
        weights = weights_final[[par_name]],
        penalty = penalty,
        lambda = lambda
      )
      df_s[[par_name]][[smooth_name]] <- stats::setNames(value, smooth_name)
    }
  }

  df_s
}

# Re-evaluate all working weights at one common final coefficient state. This
# makes no WLS update and performs no smoothing search or convergence test.
.gl_refresh_rs_final_weights <- function(
    par_cov, par_s, mm, rs_calc_eta, margin_dist, copula_dist, dataset,
    pair_cache, margin_eval_cache, copula_link, discrete_score_method,
    include_dlcopdpar, check_elapsed_budget,
    likelihood_fn = calc_likelihood_minimal,
    score_state_fn = .gl_evaluate_rs_parameter_score_state) {
  eta <- rs_calc_eta(par_cov_current = par_cov, par_s_current = par_s)
  likelihood <- likelihood_fn(
    eta$eta_inv, mm = mm$x, margin_dist = margin_dist, copula_dist = copula_dist,
    calc_d2 = FALSE, response = dataset$response, response_margin = dataset$time,
    response_subject = dataset$subject, pair_cache = pair_cache,
    margin_eval_cache = margin_eval_cache)
  if (identical(likelihood$valid, FALSE) || !is.finite(likelihood$log_lik["joint"])) {
    stop("Invalid likelihood while refreshing final RS working weights.", call. = FALSE)
  }
  weights <- list()
  for (pn in names(mm$x)) {
    check_elapsed_budget("RS final working-weight block")
    state <- score_state_fn(
      par_name = pn, eta = eta$eta, eta_inv = eta$eta_inv,
      Fx_1_2 = likelihood$Fx_1_2, copula_dist = copula_dist,
      calc_lik_out = likelihood, mm = mm, margin_dist = margin_dist, dataset = dataset,
      pair_cache = pair_cache, discrete_score_method = discrete_score_method,
      include_dlcopdpar = include_dlcopdpar, margin_deriv = likelihood$margin_deriv,
      copula_d = likelihood$copula_d, log_lik = likelihood$log_lik,
      check_dlcopdpar_gradient = FALSE, outer_only_run_counter = 0L,
      verbose = 0, timer = numeric(), timer_start = Sys.time())
    inputs <- .gl_prepare_rs_score_inputs(eta$eta, eta$eta_dr, state$d1, pn)
    score <- score_function_v2(inputs$eta, inputs$d1, -(inputs$d1 * inputs$d1), inputs$eta_dr)
    weights[[pn]] <- as.numeric(score$w_k)
  }
  list(weights = weights, likelihood = likelihood)
}
