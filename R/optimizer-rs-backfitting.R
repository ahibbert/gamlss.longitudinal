#' Run one RS inner-loop backfitting block update
#'
#' @noRd
.gl_run_rs_inner_backfitting_step <- function(
    eta,
    eta_dr,
    d1,
    par_name,
    rs_design_cache,
    par_cov,
    par_s,
    timer,
    timer_start,
    rs_smooth_trust_radius,
    rs_calc_eta,
    calc_lik_out,
    pair_cache,
    margin_eval_cache,
    lambda_penalty_K,
    lambda_s,
    margin_dist,
    copula_dist,
    dataset,
    mm,
    copula_link,
    df_s,
    step_size,
    rs_update_lambda,
    inner_run_counter,
    outer_run_counter,
    outer_only_run_counter,
    verbose,
    rs_block_trace,
    run_counter,
    use_backtracking,
    backtracking_max_halves,
    plot_results,
    log_lik_history,
    par_history,
    true_val,
    backfitting_inputs_fn = .gl_rs_backfitting_inputs,
    timer_fn = .gl_record_rs_timer_step,
    runner_fn = .gl_build_rs_backfitting_runner,
    acceptance_fn = .gl_run_rs_backfitting_acceptance) {
  backfitting_inputs <- backfitting_inputs_fn(
    eta = eta,
    eta_dr = eta_dr,
    d1 = d1,
    par_name = par_name,
    rs_design_cache = rs_design_cache,
    par_cov = par_cov,
    par_s = par_s
  )

  timer <- timer_fn(timer, timer_start, "Backfitting")

  backfitting_iteration <- runner_fn(
    design_info = backfitting_inputs$design_info,
    w_k_vec = backfitting_inputs$w_k_vec,
    z_k = backfitting_inputs$z_k,
    rs_smooth_trust_radius = rs_smooth_trust_radius,
    rs_calc_eta = rs_calc_eta,
    calc_lik_out = calc_lik_out,
    pair_cache = pair_cache,
    margin_eval_cache = margin_eval_cache
  )

  K <- lambda_penalty_K

  rs_step <- acceptance_fn(
    lambda_s = lambda_s,
    par_name = par_name,
    par_s = par_s,
    par_cov = par_cov,
    beta_start = backfitting_inputs$beta_start,
    K = K,
    margin_dist = margin_dist,
    copula_dist = copula_dist,
    dataset = dataset,
    mm = mm,
    copula_link = copula_link,
    df_s = df_s,
    step_size = step_size,
    rs_update_lambda = rs_update_lambda,
    inner_run_counter = inner_run_counter,
    outer_run_counter = outer_run_counter,
    outer_only_run_counter = outer_only_run_counter,
    verbose = verbose,
    backfitting_fn = backfitting_iteration,
    calc_lik_out = calc_lik_out,
    rs_block_trace = rs_block_trace,
    run_counter = run_counter,
    timer_start = timer_start,
    use_backtracking = use_backtracking,
    backtracking_max_halves = backtracking_max_halves,
    plot_results = plot_results,
    log_lik_history = log_lik_history,
    par_history = par_history,
    true_val = true_val
  )

  timer <- timer_fn(
    timer,
    timer_start,
    rs_step$timer_label,
    elapsed_sec = rs_step$elapsed_sec
  )

  list(
    backfitting_inputs = backfitting_inputs,
    score = backfitting_inputs$score,
    d1 = backfitting_inputs$d1,
    eta_dr_vec = backfitting_inputs$eta_dr_vec,
    lambda_s = rs_step$lambda_s,
    par_cov = rs_step$par_cov,
    par_s = rs_step$par_s,
    calc_lik_out_end = rs_step$calc_lik_out_end,
    df_s = rs_step$df_s,
    rs_block_trace = rs_step$rs_block_trace,
    change_log_lik = rs_step$change_log_lik,
    run_counter = rs_step$run_counter,
    outer_run_counter = rs_step$outer_run_counter,
    inner_run_counter = rs_step$inner_run_counter,
    rs_step = rs_step,
    timer = timer
  )
}
