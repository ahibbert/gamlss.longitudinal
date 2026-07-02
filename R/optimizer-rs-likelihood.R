#' Evaluate the RS likelihood state for one inner iteration
#'
#' @noRd
.gl_evaluate_rs_iteration_likelihood_state <- function(
    rs_calc_eta,
    par_cov,
    par_s,
    mm,
    margin_dist,
    copula_dist,
    dataset,
    pair_cache,
    margin_eval_cache,
    current_calc_lik_out = NULL,
    par_name = NULL,
    first_outer_run,
    outer_start_log_lik,
    timer,
    timer_start,
    log_lik_history,
    par_history,
    validate_eta_fn = .gl_validate_rs_eta_lengths,
    likelihood_context_fn = .gl_rs_likelihood_context,
    outer_start_fn = .gl_update_rs_outer_start_state,
    timer_fn = .gl_record_rs_timer_step,
    history_fn = .gl_append_optimizer_history,
    state_builder_fn = .gl_build_rs_iteration_likelihood_state) {
  eta_out <- rs_calc_eta(par_cov_current = par_cov, par_s_current = par_s)
  eta <- eta_out$eta
  eta_dr <- eta_out$eta_dr
  eta_inv <- eta_out$eta_inv

  validate_eta_fn(
    eta_inv = eta_inv,
    mm = mm,
    response = dataset$response
  )

  likelihood_context_args <- list(
    eta_inv = eta_inv,
    mm = mm,
    margin_dist = margin_dist,
    copula_dist = copula_dist,
    dataset = dataset,
    pair_cache = pair_cache,
    margin_eval_cache = margin_eval_cache
  )
  likelihood_context_formals <- names(formals(likelihood_context_fn))
  if ("..." %in% likelihood_context_formals || "par_name" %in% likelihood_context_formals) {
    likelihood_context_args$par_name <- par_name
  }
  if ("..." %in% likelihood_context_formals || "current_calc_lik_out" %in% likelihood_context_formals) {
    likelihood_context_args$current_calc_lik_out <- current_calc_lik_out
  }
  likelihood_context <- do.call(likelihood_context_fn, likelihood_context_args)

  outer_start_state <- outer_start_fn(
    first_outer_run = first_outer_run,
    outer_start_log_lik = outer_start_log_lik,
    log_lik = likelihood_context$log_lik
  )

  timer <- timer_fn(timer, timer_start, "Calc Lik")
  timer <- timer_fn(timer, timer_start, "Numerical Derivatives")

  history_state <- history_fn(
    log_lik_history = log_lik_history,
    par_history = par_history,
    calc_lik_out = likelihood_context$calc_lik_out,
    par_cov = par_cov
  )

  state_builder_fn(
    eta_out = eta_out,
    eta = eta,
    eta_dr = eta_dr,
    eta_inv = eta_inv,
    likelihood_context = likelihood_context,
    outer_start_state = outer_start_state,
    timer = timer,
    history_state = history_state
  )
}

#' Evaluate RS derivative and parameter-score state for one inner iteration
#'
#' @noRd
.gl_evaluate_rs_parameter_score_state <- function(
    par_name,
    eta,
    eta_inv,
    Fx_1_2,
    copula_dist,
    calc_lik_out,
    mm,
    margin_dist,
    dataset,
    pair_cache,
    discrete_score_method,
    include_dlcopdpar,
    margin_deriv,
    copula_d,
    log_lik,
    check_dlcopdpar_gradient,
    outer_only_run_counter,
    verbose,
    timer,
    timer_start,
    copula_context_fn = .gl_rs_copula_derivative_context,
    timer_fn = .gl_record_rs_timer_step,
    discrete_scores_fn = .gl_rs_discrete_scores,
    parameter_score_fn = .gl_rs_parameter_score,
    print_fn = print,
    state_builder_fn = .gl_build_rs_parameter_score_state) {
  copula_context_args <- list(
    eta_inv = eta_inv,
    Fx_1_2 = Fx_1_2,
    copula_dist = copula_dist,
    calc_lik_out = calc_lik_out
  )
  copula_context_formals <- names(formals(copula_context_fn))
  if ("..." %in% copula_context_formals || "par_name" %in% copula_context_formals) {
    copula_context_args$par_name <- par_name
  }
  if ("..." %in% copula_context_formals || "include_dlcopdpar" %in% copula_context_formals) {
    copula_context_args$include_dlcopdpar <- include_dlcopdpar
  }
  copula_context <- do.call(copula_context_fn, copula_context_args)

  timer <- timer_fn(timer, timer_start, "Copula Derivatives")

  discrete_scores <- discrete_scores_fn(
    calc_lik_out = calc_lik_out,
    eta_inv = eta_inv,
    mm = mm,
    margin_dist = margin_dist,
    copula_dist = copula_dist,
    dataset = dataset,
    pair_cache = pair_cache,
    discrete_score_method = discrete_score_method
  )

  score_assembly <- parameter_score_fn(
    par_name = par_name,
    discrete_scores = discrete_scores,
    include_dlcopdpar = include_dlcopdpar,
    eta = eta,
    eta_inv = eta_inv,
    response = dataset$response,
    margin_deriv = margin_deriv,
    margin_dist = margin_dist,
    copula_dist = copula_dist,
    dataset = dataset,
    mm = mm,
    calc_lik_out = calc_lik_out,
    copula_derivatives = copula_context$copula_derivatives,
    dcdu1 = copula_context$dcdu1,
    dcdu2 = copula_context$dcdu2,
    copula_d = copula_d,
    log_lik = log_lik,
    pair_cache = pair_cache,
    check_dlcopdpar_gradient = check_dlcopdpar_gradient,
    outer_only_run_counter = outer_only_run_counter,
    verbose = verbose
  )

  if (identical(score_assembly$path, "margin")) {
    timer <- timer_fn(timer, timer_start, "Margin Derivatives")

    if (verbose >= 4) {
      print_fn(timer)
    }
  }

  state_builder_fn(
    copula_context = copula_context,
    discrete_scores = discrete_scores,
    score_assembly = score_assembly,
    timer = timer
  )
}
