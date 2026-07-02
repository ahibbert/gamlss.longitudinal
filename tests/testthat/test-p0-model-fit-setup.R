test_that("fit optimizer context sequences setup helpers and returns optimizer fields", {
  calls <- character()
  captured <- list()

  mm <- list(x = list(mu = matrix(1, nrow = 2)), s = list(mu = list()))
  dataset <- data.frame(response = c(1, 2), time = c(1, 2), subject = c("a", "a"))
  par_cov <- c(`mu.(Intercept)` = 1)
  par_s <- list(mu = list())
  df_s <- list(mu = list())
  lambda_s <- list(mu = list())
  rs_design_cache <- list(cache = TRUE)
  rs_calc_eta <- function(...) list(eta = TRUE)
  pair_cache <- list(pair = TRUE)
  margin_eval_cache <- list(margin = TRUE)

  out <- gamlss.longitudinal:::.gl_initialize_fit_optimizer_context(
    start_from = NA,
    warm_start_par_s = NULL,
    mm = mm,
    margin_dist = "NO",
    copula_dist = "N",
    dataset = dataset,
    lambda_start = 7,
    start_step_size = 0.5,
    copula_link = "identity",
    inner_stop_crit = NA,
    outer_stop_crit = NA,
    cg_grad_tol = NA,
    cg_step_tol = NA,
    method = "RS",
    verbose = 0,
    initial_state_fn = function(start_from, warm_start_par_s, mm, margin_dist,
                                copula_dist, dataset, lambda_start) {
      calls <<- c(calls, "initial")
      captured$initial <<- list(
        start_from = start_from,
        warm_start_par_s = warm_start_par_s,
        mm = mm,
        margin_dist = margin_dist,
        copula_dist = copula_dist,
        dataset = dataset,
        lambda_start = lambda_start
      )
      list(par_cov = par_cov, par_s = par_s, df_s = df_s, lambda_s = lambda_s)
    },
    rs_design_cache_fn = function(mm, par_s) {
      calls <<- c(calls, "design")
      captured$design <<- list(mm = mm, par_s = par_s)
      rs_design_cache
    },
    optimizer_state_fn = function(par_cov, start_step_size, dataset, margin_dist) {
      calls <<- c(calls, "optimizer")
      captured$optimizer <<- list(
        par_cov = par_cov,
        start_step_size = start_step_size,
        dataset = dataset,
        margin_dist = margin_dist
      )
      list(
        first_outer_run = FALSE,
        outer_log_lik_change = 1,
        outer_start_log_lik = 2,
        outer_end_log_lik = 3,
        log_lik_history = matrix(4),
        par_history = matrix(5),
        cg_stop_reason = "not-run",
        cg_last_grad_inf = 6,
        cg_last_step_l2 = 7,
        cg_best_raw_loglik = 8,
        cg_best_iteration = 9L,
        cg_raw_loglik_drop_from_best = 10,
        rs_block_trace = data.frame(),
        outer_run_counter = 11L,
        outer_only_run_counter = 12L,
        outer_negative_streak = 13L,
        step_size = 0.5,
        weights_final = list(),
        pair_cache = pair_cache,
        margin_eval_cache = margin_eval_cache
      )
    },
    rs_eta_calculator_fn = function(rs_design_cache, mm, margin_dist, copula_link) {
      calls <<- c(calls, "eta")
      captured$eta <<- list(
        rs_design_cache = rs_design_cache,
        mm = mm,
        margin_dist = margin_dist,
        copula_link = copula_link
      )
      rs_calc_eta
    },
    stop_criteria_fn = function(inner_stop_crit, outer_stop_crit, cg_grad_tol,
                                cg_step_tol, method, par_cov, par_s, mm,
                                margin_dist, copula_dist, copula_link, dataset,
                                pair_cache, margin_eval_cache, verbose) {
      calls <<- c(calls, "stopping")
      captured$stopping <<- list(
        inner_stop_crit = inner_stop_crit,
        outer_stop_crit = outer_stop_crit,
        cg_grad_tol = cg_grad_tol,
        cg_step_tol = cg_step_tol,
        method = method,
        par_cov = par_cov,
        par_s = par_s,
        mm = mm,
        margin_dist = margin_dist,
        copula_dist = copula_dist,
        copula_link = copula_link,
        dataset = dataset,
        pair_cache = pair_cache,
        margin_eval_cache = margin_eval_cache,
        verbose = verbose
      )
      list(
        inner_stop_crit = 0.1,
        outer_stop_crit = 0.2,
        cg_grad_tol_eff = 0.3,
        cg_step_tol_eff = 0.4
      )
    }
  )

  expect_equal(calls, c("initial", "design", "optimizer", "eta", "stopping"))
  expect_identical(captured$initial$mm, mm)
  expect_identical(captured$design$par_s, par_s)
  expect_equal(captured$optimizer$par_cov, par_cov)
  expect_identical(captured$eta$rs_design_cache, rs_design_cache)
  expect_identical(captured$stopping$pair_cache, pair_cache)
  expect_identical(captured$stopping$margin_eval_cache, margin_eval_cache)

  expect_equal(out$par_cov, par_cov)
  expect_identical(out$par_s, par_s)
  expect_identical(out$rs_design_cache, rs_design_cache)
  expect_identical(out$rs_calc_eta, rs_calc_eta)
  expect_equal(out$outer_only_run_counter, 12L)
  expect_equal(out$inner_stop_crit, 0.1)
  expect_equal(out$outer_stop_crit, 0.2)
  expect_equal(out$cg_grad_tol_eff, 0.3)
  expect_equal(out$cg_step_tol_eff, 0.4)
})

test_that("fit optimizer context bundle preserves initialized state fields", {
  rs_calc_eta <- function(...) list()
  initial_state <- list(
    par_cov = c(`mu.(Intercept)` = 1),
    par_s = list(mu = list()),
    df_s = list(mu = 1),
    lambda_s = list(mu = 2)
  )
  optimizer_state <- list(
    first_outer_run = TRUE,
    outer_log_lik_change = 0.1,
    outer_start_log_lik = -5,
    outer_end_log_lik = -4,
    log_lik_history = matrix(-4, nrow = 1),
    par_history = matrix(1, nrow = 1),
    cg_stop_reason = "not-run",
    cg_last_grad_inf = 0.2,
    cg_last_step_l2 = 0.3,
    cg_best_raw_loglik = -4,
    cg_best_iteration = 2L,
    cg_raw_loglik_drop_from_best = 0,
    rs_block_trace = list(data.frame(parameter = "mu")),
    outer_run_counter = 3L,
    outer_only_run_counter = 4L,
    outer_negative_streak = 0L,
    step_size = 0.5,
    weights_final = list(mu = c(1, 1)),
    pair_cache = list(pair = TRUE),
    margin_eval_cache = list(margin = TRUE)
  )
  stop_criteria <- list(
    inner_stop_crit = 0.01,
    outer_stop_crit = 0.02,
    cg_grad_tol_eff = 0.03,
    cg_step_tol_eff = 0.04
  )

  out <- gamlss.longitudinal:::.gl_fit_optimizer_context_bundle(
    initial_state = initial_state,
    rs_design_cache = list(cache = TRUE),
    rs_calc_eta = rs_calc_eta,
    optimizer_state = optimizer_state,
    stop_criteria = stop_criteria
  )

  expect_equal(out$par_cov, initial_state$par_cov)
  expect_identical(out$par_s, initial_state$par_s)
  expect_identical(out$rs_design_cache, list(cache = TRUE))
  expect_identical(out$rs_calc_eta, rs_calc_eta)
  expect_equal(out$outer_only_run_counter, 4L)
  expect_identical(out$weights_final, optimizer_state$weights_final)
  expect_identical(out$pair_cache, optimizer_state$pair_cache)
  expect_identical(out$margin_eval_cache, optimizer_state$margin_eval_cache)
  expect_equal(out$inner_stop_crit, 0.01)
  expect_equal(out$outer_stop_crit, 0.02)
  expect_equal(out$cg_grad_tol_eff, 0.03)
  expect_equal(out$cg_step_tol_eff, 0.04)
})
