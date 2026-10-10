test_that("RS BFGS refinement trigger is one-shot and scheduled", {
  due <- gamlss.longitudinal:::.gl_should_run_rs_bfgs_refinement
  expect_false(due(9L, TRUE, 10L, FALSE))
  expect_true(due(10L, TRUE, 10L, FALSE))
  expect_true(due(11L, TRUE, 10L, FALSE))
  expect_false(due(10L, FALSE, 10L, FALSE))
  expect_false(due(10L, TRUE, 10L, TRUE))
  expect_true(due(
    12L, TRUE, 2L, FALSE,
    refinements_completed = 1L,
    max_refinements = 2L,
    lambda_update_scheduled = TRUE,
    objective_change = 2e-4,
    outer_tol = 1e-4
  ))
  expect_false(due(
    12L, TRUE, 2L, FALSE,
    refinements_completed = 1L,
    max_refinements = 2L,
    lambda_update_scheduled = TRUE,
    objective_change = 5e-5,
    outer_tol = 1e-4
  ))
  expect_true(due(
    12L, TRUE, 2L, FALSE,
    refinements_completed = 2L,
    max_refinements = 3L,
    lambda_update_scheduled = TRUE,
    objective_change = 5e-5,
    outer_tol = 1e-4,
    prior_refinement_incomplete = TRUE
  ))
  expect_false(due(
    11L, TRUE, 2L, FALSE,
    refinements_completed = 1L,
    max_refinements = 2L,
    lambda_update_scheduled = FALSE,
    objective_change = 2e-4,
    outer_tol = 1e-4
  ))
})

test_that("RS BFGS refinement improves a valid penalized state", {
  budgets <- character()
  mm <- list(
    x = list(mu = matrix(1, nrow = 1, dimnames = list(NULL, "(Intercept)"))),
    s = list(mu = list())
  )
  par_cov <- c(`mu.(Intercept)` = 0)

  out <- gamlss.longitudinal:::.gl_run_rs_bfgs_refinement(
    mm = mm,
    par_cov = par_cov,
    par_s = list(mu = list()),
    lambda_s = list(mu = list()),
    margin_dist = list(),
    copula_link = list(),
    copula_dist = "N",
    dataset = data.frame(response = 1, time = 1, subject = 1),
    pair_cache = list(),
    margin_eval_cache = list(),
    include_dlcopdpar = FALSE,
    max_iter = 20L,
    reltol = 1e-12,
    check_elapsed_budget = function(stage) budgets <<- c(budgets, stage),
    build_model_fn = function(mm, par_cov, par_s) {
      list(mm = mm, beta = par_cov)
    },
    build_penalty_fn = function(beta_names, lambda_current, par_s, mm) {
      matrix(0, length(beta_names), length(beta_names),
        dimnames = list(beta_names, beta_names))
    },
    evaluate_fn = function(beta_vec, ...) {
      list(
        loglik = -(unname(beta_vec) - 2)^2,
        eta_out = list(),
        calc_lik = list()
      )
    },
    gradient_fn = function(beta, ...) {
      stats::setNames(-2 * (unname(beta) - 2), names(beta))
    }
  )

  expect_true(out$diagnostics$accepted)
  expect_equal(unname(out$par_cov), 2, tolerance = 1e-6)
  expect_gt(out$diagnostics$improvement, 0)
  expect_true(length(budgets) > 0L)
  expect_true(all(budgets == "RS BFGS refinement"))
})

test_that("RS BFGS refinement scales sparse large coefficients", {
  captured_control <- NULL
  par_cov <- c(a = 0.2, b = -7, c = 100)
  out <- gamlss.longitudinal:::.gl_run_rs_bfgs_refinement(
    mm = list(x = list(mu = matrix(1, 1, 1)), s = list(mu = list())),
    par_cov = par_cov,
    par_s = list(mu = list()),
    lambda_s = list(mu = list()),
    margin_dist = list(),
    copula_link = list(),
    copula_dist = "N",
    dataset = data.frame(response = 1, time = 1, subject = 1),
    pair_cache = list(),
    margin_eval_cache = list(),
    include_dlcopdpar = FALSE,
    max_iter = 1L,
    reltol = 1e-12,
    parscale_max = 5,
    check_elapsed_budget = function(stage) NULL,
    build_model_fn = function(mm, par_cov, par_s) list(mm = mm, beta = par_cov),
    build_penalty_fn = function(beta_names, ...) matrix(0, 3, 3),
    evaluate_fn = function(beta_vec, ...) {
      list(loglik = -sum(beta_vec^2), eta_out = list(), calc_lik = list())
    },
    gradient_fn = function(beta, ...) -2 * beta,
    optim_fn = function(par, fn, gr, method, control) {
      captured_control <<- control
      list(par = par, value = fn(par), counts = c("function" = 1L, gradient = 1L),
           convergence = 0L, message = NULL)
    }
  )
  expect_equal(captured_control$parscale, c(1, 5, 5))
  expect_equal(out$diagnostics$parscale_min, 1)
  expect_equal(out$diagnostics$parscale_max, 5)
  expect_identical(out$diagnostics$evaluation_calls, 3L)
  expect_identical(out$diagnostics$evaluation_cache_hits, 2L)
  expect_identical(out$diagnostics$evaluation_cache_misses, 1L)
})

test_that("RS BFGS refinement can restrict optimization to named parameter blocks", {
  captured_start <- NULL
  captured_gradient <- NULL
  par_cov <- c(`mu.(Intercept)` = 1, `theta.(Intercept)` = -2)

  out <- gamlss.longitudinal:::.gl_run_rs_bfgs_refinement(
    mm = list(
      x = list(mu = matrix(1, 1, 1), theta = matrix(1, 1, 1)),
      s = list(mu = list(), theta = list())
    ),
    par_cov = par_cov,
    par_s = list(mu = list(), theta = list()),
    lambda_s = list(mu = list(), theta = list()),
    margin_dist = list(),
    copula_link = list(),
    copula_dist = "N",
    dataset = data.frame(response = 1, time = 1, subject = 1),
    pair_cache = list(),
    margin_eval_cache = list(),
    include_dlcopdpar = FALSE,
    max_iter = 1L,
    reltol = 1e-12,
    active_parameters = "theta",
    check_elapsed_budget = function(stage) NULL,
    build_model_fn = function(mm, par_cov, par_s) list(mm = mm, beta = par_cov),
    build_penalty_fn = function(beta_names, ...) {
      matrix(0, length(beta_names), length(beta_names))
    },
    evaluate_fn = function(beta_vec, ...) {
      list(
        loglik = -sum((beta_vec - c(1, 3))^2),
        eta_out = list(),
        calc_lik = list()
      )
    },
    gradient_fn = function(beta, ...) -2 * (beta - c(1, 3)),
    optim_fn = function(par, fn, gr, method, control) {
      captured_start <<- par
      captured_gradient <<- gr(par)
      list(
        par = stats::setNames(3, names(par)),
        value = fn(stats::setNames(3, names(par))),
        counts = c("function" = 1L, gradient = 1L),
        convergence = 0L,
        message = NULL
      )
    }
  )

  expect_named(captured_start, "theta.(Intercept)")
  expect_named(captured_gradient, "theta.(Intercept)")
  expect_equal(out$par_cov[["mu.(Intercept)"]], 1)
  expect_equal(out$par_cov[["theta.(Intercept)"]], 3)
  expect_identical(out$diagnostics$active_parameters, "theta")
  expect_identical(out$diagnostics$active_coefficient_count, 1L)
})

test_that("RS coordinate rescue improves only selected large fixed effects", {
  checked <- character()
  out <- gamlss.longitudinal:::.gl_run_rs_coordinate_rescue(
    beta = c(a = 5, b = 0.2, smooth = 8),
    fixed_names = c("a", "b"),
    value_fn = function(beta) -(beta[["a"]] - 1)^2 - (beta[["b"]] - 2)^2,
    threshold = 1.5,
    radius_floor = 4,
    max_coordinates = 2L,
    tolerance = 1e-8,
    check_elapsed_budget = function(stage) checked <<- c(checked, stage)
  )

  expect_true(out$used)
  expect_true(out$improved)
  expect_equal(out$beta[["a"]], 1, tolerance = 1e-5)
  expect_equal(out$beta[["b"]], 0.2)
  expect_equal(out$beta[["smooth"]], 8)
  expect_identical(out$trace$coefficient, "a")
  expect_true(out$trace$accepted)
  expect_identical(checked, "RS BFGS coordinate rescue")
})

test_that("RS BFGS refinement safely retains state after optimizer errors", {
  mm <- list(
    x = list(mu = matrix(1, nrow = 1, dimnames = list(NULL, "(Intercept)"))),
    s = list(mu = list())
  )
  par_cov <- c(`mu.(Intercept)` = 1)
  out <- gamlss.longitudinal:::.gl_run_rs_bfgs_refinement(
    mm = mm,
    par_cov = par_cov,
    par_s = list(mu = list()),
    lambda_s = list(mu = list()),
    margin_dist = list(),
    copula_link = list(),
    copula_dist = "N",
    dataset = data.frame(response = 1, time = 1, subject = 1),
    pair_cache = list(),
    margin_eval_cache = list(),
    include_dlcopdpar = FALSE,
    max_iter = 20L,
    reltol = 1e-12,
    check_elapsed_budget = function(stage) NULL,
    build_model_fn = function(mm, par_cov, par_s) list(mm = mm, beta = par_cov),
    build_penalty_fn = function(beta_names, ...) matrix(0, 1, 1),
    evaluate_fn = function(beta_vec, ...) {
      list(loglik = -1, eta_out = list(), calc_lik = list())
    },
    optim_fn = function(...) stop("example optimizer failure")
  )

  expect_false(out$diagnostics$accepted)
  expect_identical(out$diagnostics$reason, "optim_error")
  expect_equal(out$par_cov, par_cov)
})

test_that("RS BFGS refinement falls back only for classed analytical boundary errors", {
  boundary_error <- structure(
    list(message = "gradient boundary", call = NULL),
    class = c("gamlss_longitudinal_margin_cdf_error", "error", "condition")
  )
  fallback_gradient <- NULL
  out <- gamlss.longitudinal:::.gl_run_rs_bfgs_refinement(
    mm = list(x = list(mu = matrix(1, 1, 1)), s = list(mu = list())),
    par_cov = c(`mu.(Intercept)` = 0),
    par_s = list(mu = list()),
    lambda_s = list(mu = list()),
    margin_dist = list(),
    copula_link = list(),
    copula_dist = "N",
    dataset = data.frame(response = 1, time = 1, subject = 1),
    pair_cache = list(),
    margin_eval_cache = list(),
    include_dlcopdpar = FALSE,
    max_iter = 1L,
    reltol = 1e-12,
    check_elapsed_budget = function(stage) NULL,
    build_model_fn = function(mm, par_cov, par_s) list(mm = mm, beta = par_cov),
    build_penalty_fn = function(beta_names, ...) matrix(0, 1, 1),
    evaluate_fn = function(beta_vec, ...) {
      list(loglik = -(unname(beta_vec) - 2)^2, eta_out = list(), calc_lik = list())
    },
    gradient_fn = function(beta, ...) {
      if (unname(beta) > 0.5) stop(boundary_error)
      stats::setNames(-2 * (unname(beta) - 2), names(beta))
    },
    optim_fn = function(par, fn, gr, method, control) {
      expect_true(is.finite(fn(par)))
      expect_true(is.finite(gr(par)))
      candidate <- stats::setNames(1, names(par))
      expect_true(is.finite(fn(candidate)))
      fallback_gradient <<- gr(candidate)
      list(
        par = par,
        value = fn(par),
        counts = c("function" = 3L, gradient = 2L),
        convergence = 0L,
        message = NULL
      )
    }
  )

  expect_equal(unname(fallback_gradient), -2, tolerance = 1e-5)
  expect_identical(out$diagnostics$analytical_gradient_failures, 1L)
  expect_identical(out$diagnostics$finite_gradient_fallbacks, 1L)
  expect_identical(out$diagnostics$invalid_gradient_barriers, 0L)
})

test_that("RS BFGS refinement directs invalid-gradient requests toward the last valid point", {
  barrier_gradient <- NULL
  out <- gamlss.longitudinal:::.gl_run_rs_bfgs_refinement(
    mm = list(x = list(mu = matrix(1, 1, 1)), s = list(mu = list())),
    par_cov = c(`mu.(Intercept)` = 0),
    par_s = list(mu = list()),
    lambda_s = list(mu = list()),
    margin_dist = list(),
    copula_link = list(),
    copula_dist = "N",
    dataset = data.frame(response = 1, time = 1, subject = 1),
    pair_cache = list(),
    margin_eval_cache = list(),
    include_dlcopdpar = FALSE,
    max_iter = 1L,
    reltol = 1e-12,
    check_elapsed_budget = function(stage) NULL,
    build_model_fn = function(mm, par_cov, par_s) list(mm = mm, beta = par_cov),
    build_penalty_fn = function(beta_names, ...) matrix(0, 1, 1),
    evaluate_fn = function(beta_vec, ...) {
      if (unname(beta_vec) > 0.5) return(NULL)
      list(loglik = -(unname(beta_vec) - 2)^2, eta_out = list(), calc_lik = list())
    },
    gradient_fn = function(beta, ...) {
      stats::setNames(-2 * (unname(beta) - 2), names(beta))
    },
    optim_fn = function(par, fn, gr, method, control) {
      expect_true(is.finite(gr(par)))
      candidate <- stats::setNames(1, names(par))
      expect_equal(fn(candidate), 1e100)
      barrier_gradient <<- gr(candidate)
      list(
        par = par,
        value = fn(par),
        counts = c("function" = 2L, gradient = 2L),
        convergence = 0L,
        message = NULL
      )
    }
  )

  expect_gt(unname(barrier_gradient), 0)
  expect_identical(out$diagnostics$analytical_gradient_failures, 0L)
  expect_identical(out$diagnostics$finite_gradient_fallbacks, 0L)
  expect_identical(out$diagnostics$invalid_gradient_barriers, 1L)
})

test_that("RS BFGS refinement propagates elapsed-time budget errors", {
  budget_error <- structure(
    list(message = "budget exhausted", call = NULL),
    class = c("gamlss.longitudinal_time_limit_error", "error", "condition")
  )

  expect_error(
    gamlss.longitudinal:::.gl_run_rs_bfgs_refinement(
      mm = list(x = list(mu = matrix(1, 1, 1)), s = list(mu = list())),
      par_cov = c(`mu.(Intercept)` = 0),
      par_s = list(mu = list()),
      lambda_s = list(mu = list()),
      margin_dist = list(),
      copula_link = list(),
      copula_dist = "N",
      dataset = data.frame(response = 1, time = 1, subject = 1),
      pair_cache = list(),
      margin_eval_cache = list(),
      include_dlcopdpar = FALSE,
      max_iter = 1L,
      reltol = 1e-12,
      check_elapsed_budget = function(stage) stop(budget_error),
      build_model_fn = function(mm, par_cov, par_s) list(mm = mm, beta = par_cov),
      build_penalty_fn = function(beta_names, ...) matrix(0, 1, 1),
      evaluate_fn = function(beta_vec, ...) list(loglik = 0, eta_out = list(), calc_lik = list())
    ),
    class = "gamlss.longitudinal_time_limit_error"
  )
})
