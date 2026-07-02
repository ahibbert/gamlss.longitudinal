test_that("RS optimizer runner wires loop helpers and returns updated state", {
  calls <- character()
  outer_calls <- 0L
  inner_calls <- 0L

  record <- function(label) {
    calls <<- c(calls, label)
  }

  mm <- list(x = list(mu = matrix(1, nrow = 2, ncol = 1)))
  colnames(mm$x$mu) <- "mu.(Intercept)"

  out <- gamlss.longitudinal:::.gl_run_rs_optimizer(
    first_outer_run = TRUE,
    outer_log_lik_change = 0,
    outer_stop_crit = 1e-4,
    outer_only_run_counter = 1L,
    max_outer_iter = 10L,
    check_elapsed_budget = function(stage) record(paste("budget", stage)),
    verbose = 0,
    mm = mm,
    rs_calc_eta = function(...) list(),
    par_cov = c(`mu.(Intercept)` = 0),
    par_s = list(mu = list()),
    margin_dist = "NO",
    copula_dist = "Gaussian",
    dataset = data.frame(y = 1:2),
    pair_cache = list(pair = TRUE),
    margin_eval_cache = list(margin = TRUE),
    outer_start_log_lik = 0,
    log_lik_history = matrix(ncol = 3, nrow = 0),
    par_history = matrix(ncol = 1, nrow = 0, dimnames = list(NULL, "mu.(Intercept)")),
    copula_link = "identity",
    discrete_score_method = "analytic",
    include_dlcopdpar = FALSE,
    check_dlcopdpar_gradient = FALSE,
    rs_design_cache = list(mu = list()),
    rs_smooth_trust_radius = Inf,
    lambda_penalty_K = list(),
    lambda_s = list(mu = 1),
    df_s = list(mu = 1),
    step_size = 1,
    rs_update_lambda = FALSE,
    outer_run_counter = 1L,
    rs_block_trace = list(),
    use_backtracking = FALSE,
    backtracking_max_halves = 0L,
    plot_results = FALSE,
    true_val = NULL,
    inner_stop_crit = 1e-4,
    max_inner_iter = 10L,
    outer_negative_streak = 0L,
    step_adjustment = 0.5,
    max_steps = 5L,
    start_step_size = 1,
    max_negative_outer_streak = 5L,
    weights_final = list(),
    outer_continue_fn = function(...) {
      outer_calls <<- outer_calls + 1L
      record(paste0("outer-continue-", outer_calls))
      outer_calls == 1L
    },
    inner_continue_fn = function(...) {
      inner_calls <<- inner_calls + 1L
      record(paste0("inner-continue-", inner_calls))
      inner_calls == 1L
    },
    initialize_parameter_state_fn = function() {
      record("initialize-parameter")
      list(
        first_inner_run = TRUE,
        change_log_lik = 0,
        beta_change_inner = 99,
        run_counter = 1L,
        inner_run_counter = 1L
      )
    },
    likelihood_state_fn = function(
      rs_calc_eta,
      par_cov,
      par_s,
      mm,
      margin_dist,
      copula_dist,
      dataset,
      pair_cache,
      margin_eval_cache,
      first_outer_run,
      outer_start_log_lik,
      timer,
      timer_start,
      log_lik_history,
      par_history
    ) {
      record("likelihood-state")
      list(
        eta_out = list(path = "stub"),
        eta = list(mu = c(0, 0)),
        eta_dr = list(mu = c(1, 1)),
        eta_inv = list(mu = c(0, 0)),
        calc_lik_out = list(log_lik = c(marginal = -4, copula = -1, joint = -5)),
        log_lik = c(marginal = -4, copula = -1, joint = -5),
        margin_d = list(),
        margin_p = list(),
        margin_deriv = list(),
        copula_d = list(),
        copula_p = list(),
        Fx_1_2 = matrix(0.5, nrow = 1, ncol = 2),
        order_copula = 1L,
        first_outer_run = FALSE,
        outer_start_log_lik = -5,
        timer = c(timer, "likelihood"),
        log_lik_history = rbind(log_lik_history, c(marginal = -4, copula = -1, joint = -5)),
        par_history = rbind(par_history, par_cov)
      )
    },
    score_state_fn = function(
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
      timer_start
    ) {
      record("score-state")
      list(
        Fx_1_2 = Fx_1_2,
        copula_derivatives = list(),
        dldth = NULL,
        dcdth = NULL,
        dcdu1 = NULL,
        dcdu2 = NULL,
        discrete_scores = list(),
        score_assembly = list(),
        d1 = list(mu = c(1, 1)),
        d1_m = list(),
        d1_cop = list(),
        timer = c(timer, "score")
      )
    },
    backfitting_step_fn = function(
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
      true_val
    ) {
      record("backfitting-step")
      list(
        score = list(w_k = c(0.5, 0.75)),
        d1 = d1,
        eta_dr_vec = c(1, 1),
        lambda_s = list(mu = 2),
        par_cov = c(`mu.(Intercept)` = 1.25),
        par_s = list(mu = list(updated = TRUE)),
        calc_lik_out_end = list(log_lik = c(marginal = -2, copula = -1, joint = -3)),
        df_s = list(mu = 1.5),
        rs_block_trace = c(rs_block_trace, list(list(parameter = par_name))),
        change_log_lik = 0,
        run_counter = 2L,
        outer_run_counter = 2L,
        inner_run_counter = 2L,
        timer = c(timer, "backfitting")
      )
    },
    outer_iteration_state_fn = function(
      calc_lik_out_end,
      outer_start_log_lik,
      outer_only_run_counter,
      outer_negative_streak,
      step_adjustment,
      max_steps,
      start_step_size,
      max_negative_outer_streak,
      outer_stop_crit,
      verbose
    ) {
      record("outer-state")
      list(
        step_size = 0.5,
        outer_only_run_counter = 2L,
        outer_end_log_lik = calc_lik_out_end$log_lik["joint"],
        outer_log_lik_change = calc_lik_out_end$log_lik["joint"] - outer_start_log_lik,
        outer_negative_streak = 0L
      )
    },
    time_fn = function() as.POSIXct("2024-01-01 00:00:00", tz = "UTC"),
    result_fn = function(...) {
      record("result")
      .gl_build_rs_optimizer_result(...)
    }
  )

  expect_equal(
    calls,
    c(
      "outer-continue-1",
      "budget RS outer iteration",
      "initialize-parameter",
      "inner-continue-1",
      "likelihood-state",
      "score-state",
      "backfitting-step",
      "inner-continue-2",
      "outer-state",
      "outer-continue-2",
      "result"
    )
  )
  expect_equal(out$par_cov, c(`mu.(Intercept)` = 1.25))
  expect_equal(out$lambda_s, list(mu = 2))
  expect_equal(out$df_s, list(mu = 1.5))
  expect_equal(out$weights_final$mu, c(0.5, 0.75))
  expect_equal(out$calc_lik_out_end$log_lik["joint"], c(joint = -3))
  expect_equal(out$outer_log_lik_change, c(joint = 2))
  expect_equal(out$outer_only_run_counter, 2L)
  expect_false(out$first_outer_run)
  expect_equal(nrow(out$log_lik_history), 1)
  expect_equal(nrow(out$par_history), 1)
})

test_that("RS optimizer result helper preserves returned-state contract", {
  result <- .gl_build_rs_optimizer_result(
    par_cov = c(mu = 1),
    par_s = list(mu = list()),
    lambda_s = list(mu = 2),
    df_s = list(mu = 1.5),
    calc_lik_out_end = list(log_lik = c(joint = -2)),
    rs_block_trace = list(list(parameter = "mu")),
    outer_run_counter = 4L,
    outer_only_run_counter = 3L,
    outer_end_log_lik = -2,
    outer_log_lik_change = 1,
    outer_negative_streak = 0L,
    step_size = 0.5,
    weights_final = list(mu = c(1, 1)),
    log_lik_history = matrix(-2, nrow = 1),
    par_history = matrix(1, nrow = 1),
    first_outer_run = FALSE,
    outer_start_log_lik = -3
  )

  expect_named(
    result,
    c(
      "par_cov",
      "par_s",
      "lambda_s",
      "df_s",
      "calc_lik_out_end",
      "rs_block_trace",
      "outer_run_counter",
      "outer_only_run_counter",
      "outer_end_log_lik",
      "outer_log_lik_change",
      "outer_negative_streak",
      "step_size",
      "weights_final",
      "log_lik_history",
      "par_history",
      "first_outer_run",
      "outer_start_log_lik"
    )
  )
  expect_equal(result$outer_log_lik_change, 1)
  expect_equal(result$weights_final, list(mu = c(1, 1)))
  expect_false(result$first_outer_run)
})
