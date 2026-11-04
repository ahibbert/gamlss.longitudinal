test_that("optimizer state initializes counters histories and caches", {
  dataset <- data.frame(
    response = c(1, 2, 1.5, 2.5),
    time = c(1, 2, 1, 2),
    subject = c("a", "a", "b", "b")
  )
  par_cov <- c(`mu.(Intercept)` = 1, `sigma.(Intercept)` = 0, `theta.(Intercept)` = 0)

  state <- gamlss.longitudinal:::.gl_initialize_optimizer_state(
    par_cov = par_cov,
    start_step_size = 0.25,
    dataset = dataset,
    margin_dist = gamlss.dist::NO()
  )

  expect_true(state$first_outer_run)
  expect_equal(state$outer_log_lik_change, 0)
  expect_equal(state$outer_start_log_lik, 0)
  expect_equal(state$outer_end_log_lik, 0)
  expect_equal(dim(state$log_lik_history), c(0L, 3L))
  expect_equal(dim(state$par_history), c(0L, length(par_cov)))
  expect_equal(colnames(state$par_history), names(par_cov))
  expect_identical(state$cg_stop_reason, NA_character_)
  expect_equal(state$cg_best_raw_loglik, -Inf)
  expect_identical(state$cg_best_iteration, NA_integer_)
  expect_equal(state$outer_run_counter, 1)
  expect_equal(state$outer_only_run_counter, 1)
  expect_equal(state$outer_negative_streak, 0)
  expect_equal(state$step_size, 0.25)
  expect_equal(state$weights_final, list())
  expect_equal(state$rs_block_trace, list())
  expect_true(is.list(state$pair_cache))
  expect_true(is.list(state$margin_eval_cache))
})

test_that("optimizer history helper appends likelihood and parameters in history order", {
  log_lik_history <- matrix(
    c(-3, -2, -5),
    nrow = 1,
    dimnames = list(NULL, c("margin", "copula", "joint"))
  )
  par_history <- matrix(
    c(1, 2),
    nrow = 1,
    dimnames = list(NULL, c("mu.(Intercept)", "sigma.(Intercept)"))
  )
  calc_lik_out <- list(log_lik = c(margin = -4, copula = -3, joint = -7))
  par_cov <- c(`sigma.(Intercept)` = 4, `mu.(Intercept)` = 3, unused = 99)

  out <- gamlss.longitudinal:::.gl_append_optimizer_history(
    log_lik_history = log_lik_history,
    par_history = par_history,
    calc_lik_out = calc_lik_out,
    par_cov = par_cov
  )

  expect_equal(nrow(out$log_lik_history), 2)
  expect_equal(colnames(out$log_lik_history), c("margin", "copula", "joint"))
  expect_equal(out$log_lik_history[2, ], calc_lik_out$log_lik)
  expect_equal(nrow(out$par_history), 2)
  expect_equal(colnames(out$par_history), colnames(par_history))
  expect_equal(out$par_history[2, ], c(`mu.(Intercept)` = 3, `sigma.(Intercept)` = 4))
})

test_that("unit weight helper initializes one final weight per observation and parameter", {
  mm <- list(
    x = list(
      mu = matrix(1, nrow = 3, ncol = 1),
      sigma = matrix(1, nrow = 2, ncol = 1)
    )
  )

  out <- gamlss.longitudinal:::.gl_initialize_unit_weights(mm)

  expect_named(out, c("mu", "sigma"))
  expect_equal(out$mu, c(1, 1, 1))
  expect_equal(out$sigma, c(1, 1))
})

test_that("RS parameter state helper initializes inner-loop bookkeeping", {
  out <- gamlss.longitudinal:::.gl_initialize_rs_parameter_state()

  expect_true(out$first_inner_run)
  expect_equal(out$change_log_lik, 0)
  expect_equal(out$beta_change_inner, 99)
  expect_equal(out$run_counter, 1)
  expect_equal(out$inner_run_counter, 1)
})

test_that("RS outer start state captures the first joint log-likelihood", {
  out <- gamlss.longitudinal:::.gl_update_rs_outer_start_state(
    first_outer_run = TRUE,
    outer_start_log_lik = 0,
    log_lik = c(margin = -2, copula = -3, joint = -5)
  )

  expect_false(out$first_outer_run)
  expect_equal(out$outer_start_log_lik, c(joint = -5))
})

test_that("RS outer start state preserves established starts", {
  out <- gamlss.longitudinal:::.gl_update_rs_outer_start_state(
    first_outer_run = FALSE,
    outer_start_log_lik = c(joint = -10),
    log_lik = c(margin = -2, copula = -3, joint = -5)
  )

  expect_false(out$first_outer_run)
  expect_equal(out$outer_start_log_lik, c(joint = -10))
})
