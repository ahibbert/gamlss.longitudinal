#' Initialise optimizer starting values
#'
#' Creates the counters, history matrices, diagnostic placeholders, and likelihood caches shared by the RS and CG optimizer paths.
#'
#' @noRd
.gl_initialize_optimizer_state <- function(
    par_cov,
    start_step_size,
    dataset,
    margin_dist) {
  pair_cache <- build_copula_pair_cache(
    response = dataset$response,
    response_margin = dataset$time,
    response_subject = dataset$subject
  )

  list(
    first_outer_run = TRUE,
    outer_log_lik_change = 0,
    outer_start_log_lik = 0,
    outer_end_log_lik = 0,
    log_lik_history = matrix(ncol = 3, nrow = 0),
    par_history = {
      out <- matrix(ncol = length(par_cov), nrow = 0)
      colnames(out) <- names(par_cov)
      out
    },
    cg_stop_reason = NA_character_,
    cg_last_grad_inf = NA_real_,
    cg_last_step_l2 = NA_real_,
    cg_best_raw_loglik = -Inf,
    cg_best_iteration = NA_integer_,
    cg_raw_loglik_drop_from_best = NA_real_,
    rs_block_trace = list(),
    outer_run_counter = 1,
    outer_only_run_counter = 1,
    outer_negative_streak = 0,
    step_size = start_step_size,
    weights_final = list(),
    pair_cache = pair_cache,
    margin_eval_cache = .build_margin_eval_cache(margin_dist, calc_d2 = FALSE)
  )
}

#' Appends optimizer likelihood and parameter histories through each iteration for tracing and diagnostics
#'
#' @noRd
.gl_append_optimizer_history <- function(
    log_lik_history,
    par_history,
    calc_lik_out,
    par_cov) {
  list(
    log_lik_history = rbind(log_lik_history, calc_lik_out$log_lik),
    par_history = rbind(par_history, par_cov[colnames(par_history)])
  )
}

#' Initialize unit final weights for optimizer paths without IRLS weights (i.e. CG)
#'
#' @noRd
.gl_initialize_unit_weights <- function(mm) {
  weights_final <- list()
  for (pn in names(mm$x)) {
    weights_final[[pn]] <- rep(1, nrow(mm$x[[pn]]))
  }
  weights_final
}

#' Initialize RS per-parameter inner-loop state
#'
#' @noRd
.gl_initialize_rs_parameter_state <- function() {
  list(
    first_inner_run = TRUE,
    change_log_lik = 0,
    beta_change_inner = 99,
    run_counter = 1,
    inner_run_counter = 1
  )
}

#' Initialize RS outer-iteration starting log-likelihood when needed
#'
#' @noRd
.gl_update_rs_outer_start_state <- function(
    first_outer_run,
    outer_start_log_lik,
    log_lik) {
  if (isTRUE(first_outer_run)) {
    outer_start_log_lik <- log_lik["joint"]
    first_outer_run <- FALSE
  }

  list(
    first_outer_run = first_outer_run,
    outer_start_log_lik = outer_start_log_lik
  )
}
