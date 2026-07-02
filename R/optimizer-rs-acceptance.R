#' Run one RS proposal, acceptance, and diagnostics iteration
#' 
#' Prepare proposal with .gl_rs_prepare_backfitting_proposal (/R/optimizer-rs-proposal.R), 
#' run backfitting with .gl_run_backfitting (/R/optimizer-backfitting.R), 
#' and accept/reject with .gl_rs_accept_backfitting_step (/R/optimizer-rs-acceptance.R).
#' 
#' @noRd
.gl_run_rs_backfitting_acceptance <- function(
    lambda_s,
    par_name,
    par_s,
    par_cov,
    beta_start,
    K,
    margin_dist,
    copula_dist,
    dataset,
    mm,
    copula_link,
    df_s,
    step_size,
    rs_update_lambda,
    inner_run_counter,
    outer_only_run_counter,
    outer_run_counter,
    verbose,
    backfitting_fn,
    calc_lik_out,
    rs_block_trace,
    run_counter,
    timer_start,
    use_backtracking,
    backtracking_max_halves,
    plot_results,
    log_lik_history,
    par_history,
    true_val,
    prepare_proposal_fn = .gl_rs_prepare_backfitting_proposal,
    accept_step_fn = .gl_rs_accept_backfitting_step,
    apply_acceptance_fn = .gl_apply_rs_acceptance_state,
    report_acceptance_fn = .gl_report_rs_acceptance,
    plot_progress_fn = .gl_plot_rs_progress,
    difftime_fn = difftime,
    sys_time_fn = Sys.time) {
  proposal <- prepare_proposal_fn(
    lambda_s = lambda_s,
    par_name = par_name,
    par_s = par_s,
    par_cov = par_cov,
    beta_start = beta_start,
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
    outer_only_run_counter = outer_only_run_counter,
    verbose = verbose,
    backfitting_fn = backfitting_fn
  )
  lambda_s <- proposal$lambda_s
  backfitting_iteration_results <- proposal$proposed_results

  current_results <- list(
    par_cov = par_cov,
    par_s = par_s,
    calc_lik_out_end = calc_lik_out,
    GAIC_lambda_k = NA_real_,
    df_s = df_s
  )

  rs_acceptance <- accept_step_fn(
    proposed_results = backfitting_iteration_results,
    current_results = current_results,
    nominal_step_size = step_size,
    use_backtracking = use_backtracking,
    backtracking_max_halves = backtracking_max_halves,
    proposal_fn = function(trial_step) {
      backfitting_fn(
        par_s = par_s,
        par_cov = par_cov,
        beta_start = beta_start,
        lambda_s = lambda_s,
        first_inner_run = FALSE,
        K = K,
        margin_dist = margin_dist,
        copula_dist = copula_dist,
        dataset = dataset,
        mm = mm,
        copula_link = copula_link,
        df_s = df_s,
        step_size = trial_step,
        par_name = par_name
      )
    },
    outer_iteration = outer_only_run_counter,
    inner_iteration = inner_run_counter,
    global_inner_iteration = outer_run_counter,
    parameter = par_name,
    elapsed_sec = as.numeric(difftime_fn(sys_time_fn(), timer_start, units = "secs"))
  )

  acceptance_state <- apply_acceptance_fn(
    rs_acceptance = rs_acceptance,
    rs_block_trace = rs_block_trace,
    calc_lik_out = calc_lik_out,
    run_counter = run_counter,
    outer_run_counter = outer_run_counter,
    inner_run_counter = inner_run_counter
  )

  report_acceptance_fn(
    par_name = par_name,
    step_size = step_size,
    use_backtracking = use_backtracking,
    rs_acceptance = rs_acceptance,
    calc_lik_out_end = acceptance_state$calc_lik_out_end,
    verbose = verbose
  )

  if (isTRUE(plot_results)) {
    plot_progress_fn(
      log_lik_history = log_lik_history,
      par_history = par_history,
      par_count = length(acceptance_state$par_cov),
      true_val = true_val
    )
  }

  list(
    lambda_s = lambda_s,
    par_cov = acceptance_state$par_cov,
    par_s = acceptance_state$par_s,
    calc_lik_out_end = acceptance_state$calc_lik_out_end,
    df_s = acceptance_state$df_s,
    rs_block_trace = acceptance_state$rs_block_trace,
    change_log_lik = acceptance_state$change_log_lik,
    run_counter = acceptance_state$run_counter,
    outer_run_counter = acceptance_state$outer_run_counter,
    inner_run_counter = acceptance_state$inner_run_counter,
    rs_acceptance = rs_acceptance,
    timer_label = "Plotting",
    elapsed_sec = as.numeric(difftime_fn(sys_time_fn(), timer_start, units = "secs"))
  )
}
