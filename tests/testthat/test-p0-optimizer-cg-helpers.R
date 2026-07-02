test_that("CG model helper augments smooth basis into fixed design", {
  B <- matrix(c(1, 0, 0, 1, 1, 1), nrow = 3)
  mm <- list(
    x = list(mu = cbind(`(Intercept)` = 1, x = c(0, 1, 2))),
    s = list(mu = list(`s(x)` = B))
  )
  par_cov <- c(`mu.(Intercept)` = 1, mu.x = 2)
  par_s <- list(mu = list(`s(x)` = c(`mu.s(x).1` = 3, `mu.s(x).2` = 4)))

  out <- gamlss.longitudinal:::.gl_build_cg_model(mm, par_cov, par_s)

  expect_equal(out$beta, c(par_cov, par_s$mu$`s(x)`))
  expect_equal(ncol(out$mm$x$mu), 4)
  expect_equal(colnames(out$mm$x$mu), c("(Intercept)", "x", "s(x).1", "s(x).2"))
  expect_equal(out$mm$s$mu, list())
})

test_that("CG beta unpack helper restores fixed and smooth state", {
  par_cov <- c(`mu.(Intercept)` = 1, mu.x = 2)
  par_s <- list(mu = list(`s(x)` = c(`mu.s(x).1` = 0, `mu.s(x).2` = 0)))
  beta <- c(`mu.(Intercept)` = 10, mu.x = 20, `mu.s(x).1` = 30, `mu.s(x).2` = 40)

  out <- gamlss.longitudinal:::.gl_unpack_cg_beta(
    beta,
    par_cov_template = par_cov,
    par_s_template = par_s
  )

  expect_equal(out$par_cov, beta[names(par_cov)])
  expect_equal(out$par_s$mu$`s(x)`, beta[names(par_s$mu$`s(x)`)])
})

test_that("CG optimizer state helper initializes augmented model penalty and counters", {
  mm <- list(x = list(mu = matrix(1, nrow = 2, ncol = 1)), s = list(mu = list()))
  par_cov <- c(`mu.(Intercept)` = 1)
  par_s <- list(mu = list(`s(x)` = c(`mu.s(x).1` = 2)))
  lambda_s <- list(mu = list(`s(x)` = 0.5))
  captured <- list()

  out <- gamlss.longitudinal:::.gl_initialize_cg_optimizer_state(
    mm = mm,
    par_cov = par_cov,
    par_s = par_s,
    lambda_s = lambda_s,
    cg_max_delta = 0.25,
    build_model_fn = function(mm, par_cov, par_s) {
      captured$model <<- list(mm = mm, par_cov = par_cov, par_s = par_s)
      list(
        mm = list(x = list(mu = cbind(`(Intercept)` = c(1, 1), sx = c(0, 1)))),
        beta = c(`mu.(Intercept)` = 1, `mu.s(x).1` = 2)
      )
    },
    build_penalty_fn = function(beta_names, lambda_current) {
      captured$penalty <<- list(beta_names = beta_names, lambda_current = lambda_current)
      diag(c(0, lambda_current$mu$`s(x)`), nrow = 2)
    }
  )

  expect_identical(captured$model$mm, mm)
  expect_equal(captured$model$par_cov, par_cov)
  expect_identical(captured$model$par_s, par_s)
  expect_equal(captured$penalty$beta_names, c("mu.(Intercept)", "mu.s(x).1"))
  expect_identical(captured$penalty$lambda_current, lambda_s)
  expect_equal(out$beta_all, c(`mu.(Intercept)` = 1, `mu.s(x).1` = 2))
  expect_equal(out$penalty_mat, diag(c(0, 0.5), nrow = 2))
  expect_equal(out$cg_trust_radius, 0.25)
  expect_equal(out$cg_stall_count, 0L)
  expect_false(out$cg_converged)
  expect_equal(out$cg_lambda_update_count, 0L)
  expect_true(out$cg_has_smooths)
  expect_s3_class(out$cg_lambda_trace, "data.frame")
  expect_equal(nrow(out$cg_lambda_trace), 0L)
  expect_equal(out$cg_step_trace, list())
})

test_that("CG optimizer state helper detects models without smooth lambdas", {
  out <- gamlss.longitudinal:::.gl_initialize_cg_optimizer_state(
    mm = list(x = list(mu = matrix(1, nrow = 2, ncol = 1))),
    par_cov = c(`mu.(Intercept)` = 1),
    par_s = list(mu = list()),
    lambda_s = list(mu = list()),
    cg_max_delta = 1,
    build_model_fn = function(...) {
      list(
        mm = list(x = list(mu = matrix(1, nrow = 2, ncol = 1))),
        beta = c(`mu.(Intercept)` = 1)
      )
    },
    build_penalty_fn = function(beta_names, lambda_current) matrix(0, nrow = length(beta_names), ncol = length(beta_names))
  )

  expect_false(out$cg_has_smooths)
  expect_equal(out$cg_trust_radius, 1)
  expect_equal(out$penalty_mat, matrix(0, nrow = 1, ncol = 1))
})

test_that("CG iteration start helper evaluates likelihood history best objective and gradient", {
  calls <- character()
  captured <- list()
  beta <- c(`mu.(Intercept)` = 1)
  mm_cg <- list(x = list(mu = matrix(1, nrow = 2, ncol = 1)))
  penalty <- matrix(0, nrow = 1, ncol = 1)
  log_lik_history <- matrix(ncol = 3, nrow = 0)
  par_history <- matrix(ncol = 1, nrow = 0, dimnames = list(NULL, "mu.(Intercept)"))
  dataset <- data.frame(response = c(1, 2), time = c(1, 2), subject = c("a", "a"))
  eval_start <- list(
    loglik = -3,
    calc_lik = list(log_lik = c(margin = -1, copula = -2, joint = -3)),
    eta_out = list(eta_inv = list(mu = c(1, 2))),
    par_cov = beta
  )

  out <- gamlss.longitudinal:::.gl_evaluate_cg_iteration_start(
    beta_vec = beta,
    mm_cg = mm_cg,
    penalty_current = penalty,
    log_lik_history = log_lik_history,
    par_history = par_history,
    best_raw_loglik = -Inf,
    best_iteration = NA_integer_,
    current_iteration = 4L,
    margin_dist = "NO",
    copula_dist = "N",
    include_dlcopdpar = TRUE,
    dataset = dataset,
    cg_gradient_method = "central",
    eval_fn = function(beta_vec, mm_cg) {
      calls <<- c(calls, "eval")
      captured$eval_beta <<- beta_vec
      captured$eval_mm_cg <<- mm_cg
      eval_start
    },
    objective_fn = function(beta_vec, loglik, penalty_current) {
      calls <<- c(calls, "objective")
      captured$objective_beta <<- beta_vec
      captured$objective_loglik <<- loglik
      captured$objective_penalty <<- penalty_current
      -2
    },
    finite_gradient_fn = function(...) {
      stop("finite gradient should be passed through but not called directly", call. = FALSE)
    },
    history_fn = function(log_lik_history, par_history, calc_lik_out, par_cov) {
      calls <<- c(calls, "history")
      captured$history_calc_lik_out <<- calc_lik_out
      captured$history_par_cov <<- par_cov
      list(
        log_lik_history = rbind(log_lik_history, calc_lik_out$log_lik),
        par_history = rbind(par_history, par_cov[colnames(par_history)])
      )
    },
    best_fn = function(candidate_loglik, best_raw_loglik, best_iteration, current_iteration) {
      calls <<- c(calls, "best")
      captured$best_candidate <<- candidate_loglik
      captured$best_previous <<- best_raw_loglik
      captured$best_iteration_previous <<- best_iteration
      captured$best_current_iteration <<- current_iteration
      list(best_raw_loglik = candidate_loglik, best_iteration = current_iteration)
    },
    gradient_fn = function(beta_vec, mm_cg, eval_start, margin_dist, copula_dist,
                           include_dlcopdpar, dataset, cg_gradient_method,
                           finite_gradient_fn) {
      calls <<- c(calls, "gradient")
      captured$gradient_beta <<- beta_vec
      captured$gradient_mm_cg <<- mm_cg
      captured$gradient_eval_start <<- eval_start
      captured$gradient_margin_dist <<- margin_dist
      captured$gradient_copula_dist <<- copula_dist
      captured$gradient_include_dlcopdpar <<- include_dlcopdpar
      captured$gradient_dataset <<- dataset
      captured$gradient_method <<- cg_gradient_method
      captured$gradient_finite_fn <<- finite_gradient_fn
      c(`mu.(Intercept)` = 0.5)
    }
  )

  expect_equal(calls, c("eval", "history", "best", "objective", "gradient"))
  expect_equal(captured$eval_beta, beta)
  expect_identical(captured$eval_mm_cg, mm_cg)
  expect_identical(captured$history_calc_lik_out, eval_start$calc_lik)
  expect_equal(captured$history_par_cov, beta)
  expect_equal(captured$best_candidate, -3)
  expect_equal(captured$best_current_iteration, 4L)
  expect_equal(captured$objective_beta, beta)
  expect_equal(captured$objective_loglik, -3)
  expect_equal(captured$objective_penalty, penalty)
  expect_identical(captured$gradient_eval_start, eval_start)
  expect_true(captured$gradient_include_dlcopdpar)
  expect_equal(captured$gradient_method, "central")
  expect_type(captured$gradient_finite_fn, "closure")
  expect_identical(out$eval_start, eval_start)
  expect_equal(out$log_lik_history[1, ], eval_start$calc_lik$log_lik)
  expect_equal(out$par_history[1, ], beta)
  expect_equal(out$outer_start_log_lik, -3)
  expect_equal(out$best_raw_loglik, -3)
  expect_equal(out$best_iteration, 4L)
  expect_equal(out$obj_start, -2)
  expect_equal(out$grad, c(`mu.(Intercept)` = 0.5))
})

test_that("CG iteration start helper preserves non-finite likelihood error", {
  expect_error(
    gamlss.longitudinal:::.gl_evaluate_cg_iteration_start(
      beta_vec = c(x = 1),
      mm_cg = list(),
      penalty_current = matrix(0, nrow = 1, ncol = 1),
      log_lik_history = matrix(ncol = 3, nrow = 0),
      par_history = matrix(ncol = 1, nrow = 0),
      best_raw_loglik = -Inf,
      best_iteration = NA_integer_,
      current_iteration = 1L,
      margin_dist = "NO",
      copula_dist = "N",
      include_dlcopdpar = FALSE,
      dataset = data.frame(response = 1, time = 1, subject = "a"),
      cg_gradient_method = "forward",
      eval_fn = function(...) list(loglik = Inf),
      objective_fn = function(...) 0,
      finite_gradient_fn = function(...) 0
    ),
    "CG failed: current likelihood is not finite.",
    fixed = TRUE
  )
})

test_that("CG curvature and line-search helper sequences Hessian lambda EDF candidates and line search", {
  calls <- character()
  captured <- list()
  dataset <- data.frame(response = c(1, 2), time = c(1, 2), subject = c("a", "a"))
  mm_cg <- list(x = list(mu = matrix(1, nrow = 2, ncol = 1)))
  beta <- c(`mu.(Intercept)` = 1)
  grad <- c(`mu.(Intercept)` = 0.4)
  penalty <- matrix(0.1, nrow = 1, ncol = 1)
  H_initial <- matrix(-2, nrow = 1, ncol = 1)
  H_adjusted <- matrix(-3, nrow = 1, ncol = 1)
  penalty_updated <- matrix(0.2, nrow = 1, ncol = 1)
  lambda_current <- list(mu = list())
  lambda_updated <- list(mu = list(`s(x)` = 2))
  lambda_trace <- data.frame()
  lambda_trace_updated <- data.frame(parameter = "mu")
  candidate_steps <- list(c(`mu.(Intercept)` = 0.1))
  line_best <- list(step = c(`mu.(Intercept)` = 0.1), step_l2 = 0.1)

  out <- gamlss.longitudinal:::.gl_prepare_cg_curvature_line_search_state(
    dataset = dataset,
    margin_dist = "NO",
    copula_dist = "N",
    mm_cg = mm_cg,
    beta_vec = beta,
    grad_vec = grad,
    lambda_current = lambda_current,
    lambda_trace = lambda_trace,
    penalty_current = penalty,
    lambda_update_count = 0L,
    update_lambda = TRUE,
    max_lambda_updates = 2L,
    lambda_update_every = 1L,
    outer_iteration = 3L,
    trust_radius = 0.5,
    max_delta = 1,
    step_tol = 0.01,
    build_penalty_fn = function(...) NULL,
    eval_fn = function(...) NULL,
    edf_fn = function(H_obs_current, penalty_current, beta_names) {
      calls <<- c(calls, "edf")
      captured$edf_H <<- H_obs_current
      captured$edf_penalty <<- penalty_current
      captured$edf_beta_names <<- beta_names
      list(mu = list(`s(x)` = 1.5))
    },
    objective_fn = function(...) NULL,
    lambda_penalty_K = 2,
    cg_zeta_hessian = "finite",
    finite_hessian_fn = function(...) NULL,
    obj_start = -10,
    armijo_c1 = 1e-4,
    line_search = "best",
    max_line_search_evals = 5L,
    use_backtracking = TRUE,
    backtracking_max_halves = 4L,
    verbose = 2,
    observed_hessian_fn = function(tmp_obj, beta_vec, mm_cg, context) {
      calls <<- c(calls, "observed_hessian")
      captured$observed_tmp_obj <<- tmp_obj
      captured$observed_beta <<- beta_vec
      captured$observed_mm_cg <<- mm_cg
      captured$observed_context <<- context
      H_initial
    },
    hessian_object_fn = function(dataset, margin_dist, copula_dist, mm_cg, beta_vec) {
      calls <<- c(calls, "hessian_object")
      captured$hessian_dataset <<- dataset
      captured$hessian_margin_dist <<- margin_dist
      captured$hessian_copula_dist <<- copula_dist
      captured$hessian_mm_cg <<- mm_cg
      captured$hessian_beta <<- beta_vec
      list(tmp = TRUE)
    },
    zeta_hessian_fn = function(H_obs, beta_vec, mm_cg, zeta_hessian, finite_hessian_fn, verbose) {
      calls <<- c(calls, "zeta")
      captured$zeta_H <<- H_obs
      captured$zeta_beta <<- beta_vec
      captured$zeta_method <<- zeta_hessian
      captured$zeta_verbose <<- verbose
      list(H_obs = H_adjusted, H_zeta_fd = matrix(9, nrow = 1, ncol = 1))
    },
    lambda_schedule_fn = function(lambda_current, lambda_trace, penalty_current, lambda_update_count,
                                  update_lambda, max_lambda_updates, lambda_update_every,
                                  outer_iteration, H_obs_current, beta_vec, grad_vec, mm_cg,
                                  trust_radius, max_delta, step_tol, build_penalty_fn,
                                  eval_fn, edf_fn, objective_fn, lambda_penalty_K, verbose) {
      calls <<- c(calls, "lambda")
      captured$lambda_current <<- lambda_current
      captured$lambda_trace <<- lambda_trace
      captured$lambda_penalty_current <<- penalty_current
      captured$lambda_count <<- lambda_update_count
      captured$lambda_update <<- update_lambda
      captured$lambda_outer_iteration <<- outer_iteration
      captured$lambda_H <<- H_obs_current
      captured$lambda_beta <<- beta_vec
      captured$lambda_grad <<- grad_vec
      captured$lambda_trust_radius <<- trust_radius
      captured$lambda_penalty_K <<- lambda_penalty_K
      list(
        lambda = lambda_updated,
        lambda_trace = lambda_trace_updated,
        penalty_mat = penalty_updated,
        trust_radius = 0.25,
        lambda_update_count = 1L,
        lambda_changed = TRUE
      )
    },
    candidate_steps_fn = function(g_pen, H_pen, trust_radius) {
      calls <<- c(calls, "candidates")
      captured$candidate_g_pen <<- g_pen
      captured$candidate_H_pen <<- H_pen
      captured$candidate_trust_radius <<- trust_radius
      candidate_steps
    },
    line_search_fn = function(candidate_steps, beta_vec, mm_cg, penalty_current, obj_start,
                              trust_radius, max_delta, armijo_c1, line_search,
                              max_line_search_evals, use_backtracking,
                              backtracking_max_halves, eval_fn, objective_fn, verbose) {
      calls <<- c(calls, "line_search")
      captured$line_candidate_steps <<- candidate_steps
      captured$line_beta <<- beta_vec
      captured$line_penalty <<- penalty_current
      captured$line_obj_start <<- obj_start
      captured$line_trust_radius <<- trust_radius
      captured$line_max_delta <<- max_delta
      captured$line_armijo_c1 <<- armijo_c1
      captured$line_search <<- line_search
      captured$line_max_evals <<- max_line_search_evals
      captured$line_backtracking <<- use_backtracking
      captured$line_halves <<- backtracking_max_halves
      captured$line_verbose <<- verbose
      list(best = line_best, line_eval_count = 4L)
    }
  )

  expect_equal(calls, c("hessian_object", "observed_hessian", "zeta", "lambda", "edf", "candidates", "line_search"))
  expect_identical(captured$hessian_dataset, dataset)
  expect_equal(captured$observed_context, "outer iteration 3")
  expect_equal(captured$zeta_H, H_initial)
  expect_equal(captured$zeta_method, "finite")
  expect_equal(captured$lambda_H, H_adjusted)
  expect_equal(captured$lambda_grad, grad)
  expect_true(captured$lambda_update)
  expect_equal(captured$edf_H, H_adjusted)
  expect_equal(captured$edf_penalty, penalty_updated)
  expect_equal(captured$candidate_g_pen, grad - as.numeric(penalty_updated %*% beta))
  expect_equal(captured$candidate_H_pen, H_adjusted - penalty_updated)
  expect_equal(captured$candidate_trust_radius, 0.25)
  expect_identical(captured$line_candidate_steps, candidate_steps)
  expect_equal(captured$line_obj_start, -10)
  expect_equal(captured$line_trust_radius, 0.25)
  expect_equal(captured$line_max_evals, 5L)
  expect_true(captured$line_backtracking)
  expect_equal(out$tmp_obj, list(tmp = TRUE))
  expect_equal(out$H_obs, H_adjusted)
  expect_equal(out$H_zeta_fd, matrix(9, nrow = 1, ncol = 1))
  expect_identical(out$lambda_s, lambda_updated)
  expect_equal(out$cg_lambda_trace, lambda_trace_updated)
  expect_equal(out$penalty_mat, penalty_updated)
  expect_equal(out$cg_trust_radius, 0.25)
  expect_equal(out$cg_lambda_update_count, 1L)
  expect_true(out$lambda_changed)
  expect_equal(out$df_s, list(mu = list(`s(x)` = 1.5)))
  expect_equal(out$g_pen, grad - as.numeric(penalty_updated %*% beta))
  expect_equal(out$H_pen, H_adjusted - penalty_updated)
  expect_identical(out$candidate_steps, candidate_steps)
  expect_identical(out$line_search_out, list(best = line_best, line_eval_count = 4L))
  expect_identical(out$best, line_best)
  expect_equal(out$line_eval_count, 4L)
})

test_that("CG beta evaluator returns likelihood state for valid coefficients", {
  dataset <- data.frame(response = c(1, 2), time = c(1, 2), subject = c("a", "a"))
  mm_cg <- list(x = list(mu = cbind(`(Intercept)` = c(1, 1))), s = list(mu = list()))
  par_cov <- c(`mu.(Intercept)` = 1)

  out <- gamlss.longitudinal:::.gl_evaluate_cg_beta(
    beta_vec = par_cov,
    mm_cg = mm_cg,
    par_cov_template = par_cov,
    par_s_template = list(mu = list()),
    margin_dist = gamlss.dist::NO(),
    copula_link = "identity",
    copula_dist = "N",
    dataset = dataset,
    pair_cache = list(pair = TRUE),
    margin_eval_cache = list(margin = TRUE),
    calc_eta_fn = function(par, model_matrix, margin_dist, copula_link, par_s) {
      expect_equal(par, par_cov)
      expect_equal(par_s, list(mu = list()))
      list(eta_inv = list(mu = c(1, 2), sigma = c(1, 1), theta = c(0.2, 0.3)))
    },
    likelihood_fn = function(eta_inv, mm, margin_dist, copula_dist, calc_d2,
                             response, response_margin, response_subject,
                             pair_cache, margin_eval_cache) {
      expect_false(calc_d2)
      expect_equal(mm, mm_cg$x)
      expect_equal(copula_dist, "N")
      expect_equal(response, dataset$response)
      expect_equal(response_margin, dataset$time)
      expect_equal(response_subject, dataset$subject)
      expect_equal(pair_cache, list(pair = TRUE))
      expect_equal(margin_eval_cache, list(margin = TRUE))
      list(log_lik = c(margin = -1, copula = -2, joint = -3))
    }
  )

  expect_equal(out$loglik, -3)
  expect_equal(unname(out$calc_lik$log_lik["joint"]), -3)
  expect_equal(out$eta_out$eta_inv$theta, c(0.2, 0.3))
  expect_equal(out$par_cov, par_cov)
  expect_equal(out$par_s, list(mu = list()))
})

test_that("CG runtime helper bundle binds fit state and injected evaluators", {
  dataset <- data.frame(response = 1, time = 1, subject = "a")
  mm <- list(
    x = list(mu = cbind(`(Intercept)` = 1)),
    s = list(mu = list())
  )
  par_cov <- c(`mu.(Intercept)` = 1)
  par_s <- list(mu = list())

  runtime <- gamlss.longitudinal:::.gl_build_cg_runtime_helpers(
    par_cov = par_cov,
    par_s = par_s,
    mm = mm,
    margin_dist = gamlss.dist::NO(),
    copula_link = "identity",
    copula_dist = "N",
    dataset = dataset,
    pair_cache = list(pair = TRUE),
    margin_eval_cache = list(margin = TRUE),
    cg_gradient_method = "central",
    cg_hessian_method = "finite",
    verbose = 0,
    calc_eta_fn = function(par, model_matrix, margin_dist, copula_link, par_s) {
      list(eta_inv = list(mu = unname(par), sigma = 1, theta = 0))
    },
    likelihood_fn = function(eta_inv, mm, margin_dist, copula_dist, calc_d2,
                             response, response_margin, response_subject,
                             pair_cache, margin_eval_cache) {
      list(log_lik = c(margin = 0, copula = 0, joint = -(eta_inv$mu - 2)^2))
    },
    analytical_hessian_fn = function(obj, progress = FALSE) {
      stop("finite Hessian path should not call analytical Hessian")
    }
  )

  expect_named(
    runtime,
    c(
      "build_model", "build_penalty", "evaluate", "objective", "gradient",
      "finite_hessian_block", "observed_hessian", "smooth_edf_list"
    )
  )

  eval_out <- runtime$evaluate(par_cov, mm)
  expect_equal(eval_out$loglik, -1)

  expect_equal(runtime$build_penalty(names(par_cov), list()), matrix(0, 1, 1, dimnames = list(names(par_cov), names(par_cov))))
  expect_equal(runtime$objective(par_cov, eval_out$loglik, matrix(0, 1, 1)), -1)
  expect_equal(unname(runtime$gradient(par_cov, eval_out$loglik, mm)), 2, tolerance = 1e-5)

  H <- runtime$observed_hessian(list(), par_cov, mm, context = "unit test")
  expect_equal(as.numeric(H), -2, tolerance = 1e-4)
  expect_equal(runtime$smooth_edf_list(H, matrix(0, 1, 1), names(par_cov)), list(mu = list()))
})

test_that("CG gradient computation delegates to finite gradient helper", {
  beta <- c(a = 1, b = 2)
  mm_cg <- list(x = TRUE)
  eval_start <- list(loglik = -10, eta_out = list(), calc_lik = list())
  dataset <- data.frame(response = 1, time = 1, subject = "a")

  out <- gamlss.longitudinal:::.gl_compute_cg_gradient(
    beta_vec = beta,
    mm_cg = mm_cg,
    eval_start = eval_start,
    margin_dist = gamlss.dist::NO(),
    copula_dist = "N",
    include_dlcopdpar = FALSE,
    dataset = dataset,
    cg_gradient_method = "central",
    finite_gradient_fn = function(beta_vec, base_ll, mm_cg_arg) {
      expect_equal(beta_vec, beta)
      expect_equal(base_ll, -10)
      expect_equal(mm_cg_arg, mm_cg)
      c(a = 3, b = 4)
    },
    analytical_gradient_fn = function(...) {
      stop("finite gradient path should not call analytical gradient")
    }
  )

  expect_equal(out, c(a = 3, b = 4))
})

test_that("CG gradient computation forwards analytical gradient inputs", {
  beta <- c(a = 1)
  mm_cg <- list(x = TRUE)
  eta_out <- list(eta_inv = list(mu = 1))
  calc_lik <- list(log_lik = c(joint = -1))
  eval_start <- list(loglik = -1, eta_out = eta_out, calc_lik = calc_lik)
  dataset <- data.frame(response = 7, time = 2, subject = "id1")
  margin_dist <- gamlss.dist::NO()

  out <- gamlss.longitudinal:::.gl_compute_cg_gradient(
    beta_vec = beta,
    mm_cg = mm_cg,
    eval_start = eval_start,
    margin_dist = margin_dist,
    copula_dist = "N",
    include_dlcopdpar = TRUE,
    dataset = dataset,
    cg_gradient_method = "analytical",
    finite_gradient_fn = function(...) {
      stop("analytical gradient path should not call finite gradient")
    },
    analytical_gradient_fn = function(beta_arg, mm_arg, eta_arg, lik_arg,
                                      margin_arg, copula_arg, include_arg,
                                      response_arg, time_arg, subject_arg) {
      expect_equal(beta_arg, beta)
      expect_equal(mm_arg, mm_cg)
      expect_equal(eta_arg, eta_out)
      expect_equal(lik_arg, calc_lik)
      expect_equal(margin_arg, margin_dist)
      expect_equal(copula_arg, "N")
      expect_true(include_arg)
      expect_equal(response_arg, dataset$response)
      expect_equal(time_arg, dataset$time)
      expect_equal(subject_arg, dataset$subject)
      c(a = 5)
    }
  )

  expect_equal(out, c(a = 5))
})

test_that("CG beta evaluator rejects non-finite eta values", {
  dataset <- data.frame(response = 1, time = 1, subject = "a")
  mm_cg <- list(x = list(mu = matrix(1, nrow = 1)), s = list(mu = list()))

  out <- gamlss.longitudinal:::.gl_evaluate_cg_beta(
    beta_vec = c(`mu.(Intercept)` = 1),
    mm_cg = mm_cg,
    par_cov_template = c(`mu.(Intercept)` = 1),
    par_s_template = list(mu = list()),
    margin_dist = gamlss.dist::NO(),
    copula_link = "identity",
    copula_dist = "N",
    dataset = dataset,
    pair_cache = list(),
    margin_eval_cache = list(),
    calc_eta_fn = function(...) list(eta_inv = list(mu = Inf, sigma = 1, theta = 0)),
    likelihood_fn = function(...) stop("invalid eta should skip likelihood")
  )

  expect_null(out)
})

test_that("CG beta evaluator rejects Gaussian and t copula boundary theta", {
  dataset <- data.frame(response = 1, time = 1, subject = "a")
  mm_cg <- list(x = list(mu = matrix(1, nrow = 1)), s = list(mu = list()))

  for(copula in c("N", "t")) {
    out <- gamlss.longitudinal:::.gl_evaluate_cg_beta(
      beta_vec = c(`mu.(Intercept)` = 1),
      mm_cg = mm_cg,
      par_cov_template = c(`mu.(Intercept)` = 1),
      par_s_template = list(mu = list()),
      margin_dist = gamlss.dist::NO(),
      copula_link = "identity",
      copula_dist = copula,
      dataset = dataset,
      pair_cache = list(),
      margin_eval_cache = list(),
      calc_eta_fn = function(...) list(eta_inv = list(mu = 1, sigma = 1, theta = 0.999)),
      likelihood_fn = function(...) stop("boundary theta should skip likelihood")
    )

    expect_null(out)
  }
})

test_that("CG beta evaluator rejects invalid positive distribution parameters", {
  dataset <- data.frame(response = 1, time = 1, subject = "a")
  mm_cg <- list(x = list(mu = matrix(1, nrow = 1)), s = list(mu = list()))

  out <- gamlss.longitudinal:::.gl_evaluate_cg_beta(
    beta_vec = c(`mu.(Intercept)` = 1),
    mm_cg = mm_cg,
    par_cov_template = c(`mu.(Intercept)` = 1),
    par_s_template = list(mu = list()),
    margin_dist = gamlss.dist::NO(),
    copula_link = "identity",
    copula_dist = "N",
    dataset = dataset,
    pair_cache = list(),
    margin_eval_cache = list(),
    calc_eta_fn = function(...) list(eta_inv = list(mu = 1, sigma = 0, theta = 0)),
    likelihood_fn = function(...) stop("non-positive sigma should skip likelihood")
  )

  expect_null(out)
})

test_that("CG beta evaluator converts likelihood errors to invalid candidate", {
  dataset <- data.frame(response = 1, time = 1, subject = "a")
  mm_cg <- list(x = list(mu = matrix(1, nrow = 1)), s = list(mu = list()))

  out <- gamlss.longitudinal:::.gl_evaluate_cg_beta(
    beta_vec = c(`mu.(Intercept)` = 1),
    mm_cg = mm_cg,
    par_cov_template = c(`mu.(Intercept)` = 1),
    par_s_template = list(mu = list()),
    margin_dist = gamlss.dist::NO(),
    copula_link = "identity",
    copula_dist = "N",
    dataset = dataset,
    pair_cache = list(),
    margin_eval_cache = list(),
    calc_eta_fn = function(...) list(eta_inv = list(mu = 1, sigma = 1, theta = 0)),
    likelihood_fn = function(...) stop("bad candidate")
  )

  expect_null(out)
})

test_that("CG Hessian object helper preserves analytical Hessian input schema", {
  dataset <- data.frame(response = c(1, 2), time = c(1, 2), subject = c("a", "a"))
  mm_cg <- list(
    x = list(mu = matrix(1, nrow = 2), sigma = matrix(2, nrow = 2)),
    s = list(mu = list(), sigma = list())
  )
  beta <- c(`mu.(Intercept)` = 1, `sigma.(Intercept)` = 0)

  out <- gamlss.longitudinal:::.gl_build_cg_hessian_object(
    dataset = dataset,
    margin_dist = gamlss.dist::NO(),
    copula_dist = "N",
    mm_cg = mm_cg,
    beta_vec = beta
  )

  expect_equal(
    names(out),
    c(
      "response", "response_margin", "response_subject", "margin_dist",
      "copula_dist", "model_matrix", "par", "par_s"
    )
  )
  expect_equal(out$response, dataset$response)
  expect_equal(out$response_margin, dataset$time)
  expect_equal(out$response_subject, dataset$subject)
  expect_equal(out$copula_dist, "N")
  expect_equal(out$model_matrix, mm_cg)
  expect_equal(out$par, beta)
  expect_equal(out$par_s, list(mu = list(), sigma = list()))
})

test_that("CG best raw log-likelihood helper updates only on finite improvement", {
  improved <- gamlss.longitudinal:::.gl_update_cg_best_loglik(
    candidate_loglik = 12,
    best_raw_loglik = 10,
    best_iteration = 1L,
    current_iteration = 2L
  )
  expect_equal(improved$best_raw_loglik, 12)
  expect_equal(improved$best_iteration, 2L)

  unchanged <- gamlss.longitudinal:::.gl_update_cg_best_loglik(
    candidate_loglik = 9,
    best_raw_loglik = 10,
    best_iteration = 1L,
    current_iteration = 2L
  )
  expect_equal(unchanged$best_raw_loglik, 10)
  expect_equal(unchanged$best_iteration, 1L)

  nonfinite <- gamlss.longitudinal:::.gl_update_cg_best_loglik(
    candidate_loglik = NA_real_,
    best_raw_loglik = 10,
    best_iteration = 1L,
    current_iteration = 2L
  )
  expect_equal(nonfinite$best_raw_loglik, 10)
  expect_equal(nonfinite$best_iteration, 1L)
})

test_that("CG raw log-likelihood drop helper includes prevented deterioration", {
  expect_equal(
    gamlss.longitudinal:::.gl_cg_raw_loglik_drop(
      best_raw_loglik = 10,
      current_loglik = 8
    ),
    2
  )
  expect_equal(
    gamlss.longitudinal:::.gl_cg_raw_loglik_drop(
      best_raw_loglik = 10,
      current_loglik = 8,
      prevented_deterioration = TRUE,
      prevented_raw_loglik_drop = 3
    ),
    3
  )
  expect_equal(
    gamlss.longitudinal:::.gl_cg_raw_loglik_drop(
      best_raw_loglik = 10,
      current_loglik = 8,
      prevented_deterioration = TRUE,
      prevented_raw_loglik_drop = NA_real_
    ),
    2
  )
})

test_that("CG penalty and objective helpers preserve penalty calculations", {
  B <- matrix(c(1, 0, 0, 1, 1, 1), nrow = 3)
  penalty <- matrix(c(2, 0.5, 0.5, 3), nrow = 2)
  attr(B, "penalty") <- penalty
  mm <- list(
    x = list(mu = cbind(`(Intercept)` = 1, x = c(0, 1, 2))),
    s = list(mu = list(`s(x)` = B))
  )
  par_s <- list(mu = list(`s(x)` = c(`mu.s(x).1` = 0, `mu.s(x).2` = 0)))
  lambda_s <- list(mu = list(`s(x)` = 4))
  beta_names <- c("mu.(Intercept)", "mu.x", "mu.s(x).1", "mu.s(x).2")

  P <- gamlss.longitudinal:::.gl_build_cg_penalty(beta_names, lambda_s, par_s = par_s, mm = mm)
  expected <- matrix(0, 4, 4, dimnames = list(beta_names, beta_names))
  expected[3:4, 3:4] <- 4 * penalty

  expect_equal(P, expected)

  beta <- c(`mu.(Intercept)` = 1, mu.x = 2, `mu.s(x).1` = 3, `mu.s(x).2` = 4)
  expect_equal(
    gamlss.longitudinal:::.gl_cg_objective(beta, loglik = 10, penalty_current = P),
    10 - 0.5 * sum(as.numeric(beta) * as.numeric(P %*% beta))
  )
})

test_that("CG lambda candidate helper normalizes base values and bounds candidates", {
  expect_equal(gamlss.longitudinal:::.gl_cg_lambda_base(0), 1)
  expect_equal(gamlss.longitudinal:::.gl_cg_lambda_base(NA_real_), 1)
  expect_equal(gamlss.longitudinal:::.gl_cg_lambda_base(2), 2)

  expect_equal(
    gamlss.longitudinal:::.gl_cg_lambda_candidates(2),
    c(0.2, 0.5, 1, 2, 4, 8, 20)
  )
  expect_true(all(gamlss.longitudinal:::.gl_cg_lambda_candidates(1e9) <= 1e6))
  expect_true(all(gamlss.longitudinal:::.gl_cg_lambda_candidates(1e-6) >= 0.01))
})

test_that("CG step limiter applies trust radius before max-coordinate cap", {
  delta <- c(x = 3, y = 4)

  trust_limited <- gamlss.longitudinal:::.gl_limit_cg_step(
    delta,
    trust_radius = 2.5,
    max_delta = 10
  )
  expect_equal(sqrt(sum(trust_limited^2)), 2.5, tolerance = 1e-12)

  coordinate_limited <- gamlss.longitudinal:::.gl_limit_cg_step(
    delta,
    trust_radius = 10,
    max_delta = 2
  )
  expect_equal(max(abs(coordinate_limited)), 2, tolerance = 1e-12)

  both_limited <- gamlss.longitudinal:::.gl_limit_cg_step(
    delta,
    trust_radius = 2.5,
    max_delta = 1
  )
  expect_equal(both_limited, c(x = 0.75, y = 1), tolerance = 1e-12)
})

test_that("CG trust-radius shrink helper respects step tolerance floor", {
  expect_equal(
    gamlss.longitudinal:::.gl_shrink_cg_trust_radius(trust_radius = 10, step_tol = 0.1),
    5
  )
  expect_equal(
    gamlss.longitudinal:::.gl_shrink_cg_trust_radius(trust_radius = 0.15, step_tol = 0.1),
    0.1
  )
})

test_that("CG trust-radius expansion helper expands only boundary-sized accepted steps", {
  expect_equal(
    gamlss.longitudinal:::.gl_expand_cg_trust_radius(
      trust_radius = 2,
      step_l2 = 1.5,
      step_tol = 0.1,
      max_delta = 10
    ),
    2
  )
  expect_equal(
    gamlss.longitudinal:::.gl_expand_cg_trust_radius(
      trust_radius = 2,
      step_l2 = 1.6,
      step_tol = 0.1,
      max_delta = 10
    ),
    3
  )
  expect_equal(
    gamlss.longitudinal:::.gl_expand_cg_trust_radius(
      trust_radius = 2,
      step_l2 = 2,
      step_tol = 0.1,
      max_delta = 2.5
    ),
    2.5
  )
})

test_that("CG Hessian source paths use the expected package files", {
  paths <- gamlss.longitudinal:::.gl_cg_hessian_source_paths("pkgroot")

  expect_equal(
    basename(paths),
    c(
      "hessian-linkinv-derivatives.R",
      "hessian-fd-step.R",
      "hessian-warnings.R",
      "hessian-margin-cdf.R",
      "hessian-copula.R",
      "hessian-assembly.R",
      "hessian-analytical.R"
    )
  )
  expect_true(all(grepl("pkgroot", paths, fixed = TRUE)))
})

test_that("CG Hessian availability helper is a no-op when already loaded", {
  sourced <- character(0)

  expect_invisible(gamlss.longitudinal:::.gl_ensure_cg_hessian_available(
    package_root = "pkgroot",
    exists_fn = function(name, mode) TRUE,
    source_fn = function(file, local) sourced <<- c(sourced, file),
    file_exists_fn = function(paths) rep(TRUE, length(paths))
  ))
  expect_length(sourced, 0)
})

test_that("CG Hessian availability helper sources all Hessian files when available", {
  sourced <- character(0)
  exists_calls <- 0L
  exists_after_source <- function(name, mode) {
    exists_calls <<- exists_calls + 1L
    exists_calls > 1L
  }

  expect_invisible(gamlss.longitudinal:::.gl_ensure_cg_hessian_available(
    package_root = "pkgroot",
    exists_fn = exists_after_source,
    source_fn = function(file, local) {
      expect_false(local)
      sourced <<- c(sourced, file)
    },
    file_exists_fn = function(paths) rep(TRUE, length(paths))
  ))

  expect_equal(basename(sourced), basename(gamlss.longitudinal:::.gl_cg_hessian_source_paths()))
})

test_that("CG Hessian availability helper preserves missing-helper error", {
  expect_error(
    gamlss.longitudinal:::.gl_ensure_cg_hessian_available(
      package_root = "pkgroot",
      exists_fn = function(name, mode) FALSE,
      source_fn = function(file, local) stop("should not source"),
      file_exists_fn = function(paths) rep(FALSE, length(paths))
    ),
    "CG requires calc_analytical_hessian\\(\\); source the R/hessian-\\*.R files first\\."
  )
})

test_that("CG trust-radius expansion helper leaves non-finite inputs unchanged", {
  expect_equal(
    gamlss.longitudinal:::.gl_expand_cg_trust_radius(
      trust_radius = 2,
      step_l2 = NA_real_,
      step_tol = 0.1,
      max_delta = 10
    ),
    2
  )
})

test_that("CG candidate-step builder starts with normalized gradient step", {
  g_pen <- c(x = 3, y = 4)
  H_pen <- diag(c(-2, -3))

  steps <- gamlss.longitudinal:::.gl_build_cg_candidate_steps(
    g_pen = g_pen,
    H_pen = H_pen,
    trust_radius = 2.5
  )

  expect_true(length(steps) >= 1L)
  expect_equal(steps[[1]], c(1.5, 2), tolerance = 1e-12, ignore_attr = TRUE)
})

test_that("CG candidate-step builder appends paired ridge Newton steps", {
  g_pen <- c(x = 2)
  H_pen <- matrix(-4, nrow = 1)

  steps <- gamlss.longitudinal:::.gl_build_cg_candidate_steps(
    g_pen = g_pen,
    H_pen = H_pen,
    trust_radius = 1
  )

  ridge_values <- c(0, 1e-8, 1e-6, 1e-4, 1e-2, 1, 10, 100)
  expected_newton <- 2 / (4 + ridge_values)

  expect_equal(steps[[1]], c(1), tolerance = 1e-12, ignore_attr = TRUE)
  expect_equal(
    vapply(steps[seq(2, length(steps), by = 2)], function(x) x[[1]], numeric(1)),
    expected_newton,
    tolerance = 1e-12
  )
  expect_equal(
    vapply(steps[seq(3, length(steps), by = 2)], function(x) x[[1]], numeric(1)),
    -expected_newton,
    tolerance = 1e-12
  )
})

test_that("CG candidate-step builder skips non-finite directions", {
  steps <- gamlss.longitudinal:::.gl_build_cg_candidate_steps(
    g_pen = c(x = NA_real_),
    H_pen = matrix(0, nrow = 1),
    trust_radius = 1
  )

  expect_equal(steps, list())
})

test_that("CG lambda scorer evaluates candidates and reports best penalized objective", {
  candidates <- c(0.5, 1)
  lambda_current <- list(mu = list(`s(x)` = 1))
  beta <- c(x = 0)
  grad <- c(x = 1)
  H_obs <- matrix(-2, nrow = 1, dimnames = list("x", "x"))
  build_penalty <- function(beta_names, lambda_try) {
    matrix(lambda_try$mu$`s(x)`, nrow = 1, dimnames = list(beta_names, beta_names))
  }
  eval_candidate <- function(beta_try, mm_cg) {
    list(loglik = -as.numeric((beta_try["x"] - 0.5)^2))
  }
  edf_constant <- function(H_obs_current, penalty_current, beta_names) {
    list(mu = list(`s(x)` = 0.5))
  }

  out <- gamlss.longitudinal:::.gl_score_cg_lambda_candidates(
    candidates = candidates,
    parameter = "mu",
    smooth = "s(x)",
    lambda_current = lambda_current,
    H_obs_current = H_obs,
    beta_vec = beta,
    grad_vec = grad,
    mm_cg = NULL,
    trust_radius = Inf,
    max_delta = Inf,
    build_penalty_fn = build_penalty,
    eval_fn = eval_candidate,
    edf_fn = edf_constant,
    objective_fn = gamlss.longitudinal:::.gl_cg_objective,
    lambda_penalty_K = 2
  )

  beta_try <- 1 / (2 + candidates)
  raw_loglik <- -((beta_try - 0.5)^2)
  penalty <- beta_try^2 * candidates
  penalized <- raw_loglik - 0.5 * penalty
  gaic <- -2 * raw_loglik + 2 * 0.5

  expect_equal(out$raw_loglik, raw_loglik, tolerance = 1e-12)
  expect_equal(out$penalty_value, penalty, tolerance = 1e-12)
  expect_equal(out$penalized_loglik, penalized, tolerance = 1e-12)
  expect_equal(out$edf_values, rep(0.5, length(candidates)))
  expect_equal(out$gaic_score, gaic, tolerance = 1e-12)
  expect_equal(out$best, which.max(penalized))
})

test_that("CG lambda scorer preserves invalid-candidate sentinels", {
  candidates <- c(0.5, 1)
  lambda_current <- list(mu = list(`s(x)` = 1))
  beta <- c(x = 0)
  grad <- c(x = 1)
  H_obs <- matrix(-2, nrow = 1, dimnames = list("x", "x"))
  build_penalty <- function(beta_names, lambda_try) {
    matrix(lambda_try$mu$`s(x)`, nrow = 1, dimnames = list(beta_names, beta_names))
  }
  eval_invalid <- function(beta_try, mm_cg) NULL
  edf_constant <- function(H_obs_current, penalty_current, beta_names) {
    list(mu = list(`s(x)` = 0.5))
  }

  out <- gamlss.longitudinal:::.gl_score_cg_lambda_candidates(
    candidates = candidates,
    parameter = "mu",
    smooth = "s(x)",
    lambda_current = lambda_current,
    H_obs_current = H_obs,
    beta_vec = beta,
    grad_vec = grad,
    mm_cg = NULL,
    trust_radius = Inf,
    max_delta = Inf,
    build_penalty_fn = build_penalty,
    eval_fn = eval_invalid,
    edf_fn = edf_constant,
    objective_fn = gamlss.longitudinal:::.gl_cg_objective,
    lambda_penalty_K = 2
  )

  expect_length(out$best, 0)
  expect_true(all(is.na(out$raw_loglik)))
  expect_true(all(is.na(out$penalty_value)))
  expect_true(all(is.na(out$penalized_loglik)))
  expect_true(all(is.na(out$edf_values)))
  expect_true(all(is.infinite(out$gaic_score)))
})

test_that("CG lambda updater returns unchanged state when disabled", {
  lambda_current <- list(mu = list(`s(x)` = 1))
  lambda_trace <- data.frame(existing = 1)

  out <- gamlss.longitudinal:::.gl_update_cg_lambdas(
    lambda_current = lambda_current,
    lambda_trace = lambda_trace,
    H_obs_current = matrix(-2, nrow = 1, dimnames = list("x", "x")),
    beta_vec = c(x = 0),
    grad_vec = c(x = 1),
    mm_cg = NULL,
    trust_radius = Inf,
    outer_iteration = 1,
    update_lambda = FALSE,
    max_delta = Inf,
    build_penalty_fn = function(beta_names, lambda_try) {
      stop("disabled path should not score candidates")
    },
    eval_fn = function(beta_try, mm_cg) NULL,
    edf_fn = function(H_obs_current, penalty_current, beta_names) list(),
    objective_fn = gamlss.longitudinal:::.gl_cg_objective,
    lambda_penalty_K = 2,
    verbose = 0
  )

  expect_equal(out$lambda, lambda_current)
  expect_equal(out$lambda_trace, lambda_trace)
})

test_that("CG lambda updater applies best candidate and appends trace rows", {
  lambda_current <- list(mu = list(`s(x)` = 1))
  beta <- c(x = 0)
  grad <- c(x = 1)
  H_obs <- matrix(-2, nrow = 1, dimnames = list("x", "x"))
  build_penalty <- function(beta_names, lambda_try) {
    matrix(lambda_try$mu$`s(x)`, nrow = 1, dimnames = list(beta_names, beta_names))
  }
  eval_candidate <- function(beta_try, mm_cg) {
    list(loglik = -as.numeric((beta_try["x"] - 0.5)^2))
  }
  edf_constant <- function(H_obs_current, penalty_current, beta_names) {
    list(mu = list(`s(x)` = 0.5))
  }

  out <- gamlss.longitudinal:::.gl_update_cg_lambdas(
    lambda_current = lambda_current,
    lambda_trace = data.frame(),
    H_obs_current = H_obs,
    beta_vec = beta,
    grad_vec = grad,
    mm_cg = NULL,
    trust_radius = Inf,
    outer_iteration = 3,
    update_lambda = TRUE,
    max_delta = Inf,
    build_penalty_fn = build_penalty,
    eval_fn = eval_candidate,
    edf_fn = edf_constant,
    objective_fn = gamlss.longitudinal:::.gl_cg_objective,
    lambda_penalty_K = 2,
    verbose = 0
  )

  candidates <- gamlss.longitudinal:::.gl_cg_lambda_candidates(1)
  beta_try <- 1 / (2 + candidates)
  raw_loglik <- -((beta_try - 0.5)^2)
  penalty <- beta_try^2 * candidates
  penalized <- raw_loglik - 0.5 * penalty
  best <- which.max(penalized)

  expect_equal(out$lambda$mu$`s(x)`, candidates[best])
  expect_equal(nrow(out$lambda_trace), length(candidates))
  expect_equal(out$lambda_trace$outer_iteration, rep(3L, length(candidates)))
  expect_equal(out$lambda_trace$chosen, seq_along(candidates) == best)
})

test_that("CG lambda-change predicate compares flattened lambda values with tolerance", {
  lambda_before <- list(mu = list(`s(x)` = 1), sigma = list(`s(z)` = 2))

  expect_false(gamlss.longitudinal:::.gl_cg_lambdas_changed(
    lambda_before = lambda_before,
    lambda_after = list(mu = list(`s(x)` = 1), sigma = list(`s(z)` = 2))
  ))

  expect_false(gamlss.longitudinal:::.gl_cg_lambdas_changed(
    lambda_before = lambda_before,
    lambda_after = list(mu = list(`s(x)` = 1 + 1e-14), sigma = list(`s(z)` = 2))
  ))

  expect_true(gamlss.longitudinal:::.gl_cg_lambdas_changed(
    lambda_before = lambda_before,
    lambda_after = list(mu = list(`s(x)` = 1.1), sigma = list(`s(z)` = 2))
  ))
})

test_that("CG scheduled lambda update returns unchanged state when not scheduled", {
  lambda_current <- list(mu = list(`s(x)` = 1))
  lambda_trace <- data.frame(existing = 1)
  penalty_current <- matrix(1, nrow = 1, dimnames = list("x", "x"))

  out <- gamlss.longitudinal:::.gl_maybe_update_cg_lambdas_on_schedule(
    lambda_current = lambda_current,
    lambda_trace = lambda_trace,
    penalty_current = penalty_current,
    lambda_update_count = 0L,
    update_lambda = TRUE,
    max_lambda_updates = 5L,
    lambda_update_every = 2L,
    outer_iteration = 1L,
    H_obs_current = matrix(-1, nrow = 1),
    beta_vec = c(x = 1),
    grad_vec = c(x = 1),
    mm_cg = NULL,
    trust_radius = 10,
    max_delta = 10,
    step_tol = 1e-6,
    build_penalty_fn = function(...) stop("unscheduled update should not rebuild penalty"),
    eval_fn = function(...) NULL,
    edf_fn = function(...) list(),
    objective_fn = gamlss.longitudinal:::.gl_cg_objective,
    lambda_penalty_K = 2,
    verbose = 0,
    update_lambdas_fn = function(...) stop("unscheduled update should not call updater")
  )

  expect_equal(out$lambda, lambda_current)
  expect_equal(out$lambda_trace, lambda_trace)
  expect_equal(out$penalty_mat, penalty_current)
  expect_equal(out$trust_radius, 10)
  expect_equal(out$lambda_update_count, 0L)
  expect_false(out$lambda_changed)
})

test_that("CG scheduled lambda update increments count without shrink when lambda is unchanged", {
  lambda_current <- list(mu = list(`s(x)` = 1))
  penalty_current <- matrix(1, nrow = 1, dimnames = list("x", "x"))

  out <- gamlss.longitudinal:::.gl_maybe_update_cg_lambdas_on_schedule(
    lambda_current = lambda_current,
    lambda_trace = data.frame(),
    penalty_current = penalty_current,
    lambda_update_count = 0L,
    update_lambda = TRUE,
    max_lambda_updates = 5L,
    lambda_update_every = 2L,
    outer_iteration = 2L,
    H_obs_current = matrix(-1, nrow = 1),
    beta_vec = c(x = 1),
    grad_vec = c(x = 1),
    mm_cg = NULL,
    trust_radius = 10,
    max_delta = 10,
    step_tol = 1e-6,
    build_penalty_fn = function(beta_names, lambda_try) {
      matrix(lambda_try$mu$`s(x)`, nrow = 1, dimnames = list(beta_names, beta_names))
    },
    eval_fn = function(...) NULL,
    edf_fn = function(...) list(),
    objective_fn = gamlss.longitudinal:::.gl_cg_objective,
    lambda_penalty_K = 2,
    verbose = 0,
    update_lambdas_fn = function(...) list(lambda = lambda_current, lambda_trace = data.frame(updated = TRUE)),
    shrink_trust_radius_fn = function(...) stop("unchanged lambda should not shrink trust radius")
  )

  expect_equal(out$lambda, lambda_current)
  expect_equal(out$penalty_mat, penalty_current)
  expect_equal(out$trust_radius, 10)
  expect_equal(out$lambda_update_count, 1L)
  expect_false(out$lambda_changed)
  expect_true(out$lambda_trace$updated)
})

test_that("CG scheduled lambda update shrinks trust radius when lambda changes", {
  lambda_current <- list(mu = list(`s(x)` = 1))
  lambda_updated <- list(mu = list(`s(x)` = 2))

  out <- gamlss.longitudinal:::.gl_maybe_update_cg_lambdas_on_schedule(
    lambda_current = lambda_current,
    lambda_trace = data.frame(),
    penalty_current = matrix(1, nrow = 1, dimnames = list("x", "x")),
    lambda_update_count = 1L,
    update_lambda = TRUE,
    max_lambda_updates = 5L,
    lambda_update_every = 2L,
    outer_iteration = 4L,
    H_obs_current = matrix(-1, nrow = 1),
    beta_vec = c(x = 1),
    grad_vec = c(x = 1),
    mm_cg = NULL,
    trust_radius = 10,
    max_delta = 10,
    step_tol = 0.5,
    build_penalty_fn = function(beta_names, lambda_try) {
      matrix(lambda_try$mu$`s(x)`, nrow = 1, dimnames = list(beta_names, beta_names))
    },
    eval_fn = function(...) NULL,
    edf_fn = function(...) list(),
    objective_fn = gamlss.longitudinal:::.gl_cg_objective,
    lambda_penalty_K = 2,
    verbose = 0,
    update_lambdas_fn = function(...) list(lambda = lambda_updated, lambda_trace = data.frame(updated = TRUE)),
    shrink_trust_radius_fn = function(trust_radius, step_tol) {
      expect_equal(trust_radius, 10)
      expect_equal(step_tol, 0.5)
      5
    }
  )

  expect_equal(out$lambda, lambda_updated)
  expect_equal(out$penalty_mat, matrix(2, nrow = 1, dimnames = list("x", "x")))
  expect_equal(out$trust_radius, 5)
  expect_equal(out$lambda_update_count, 2L)
  expect_true(out$lambda_changed)
  expect_true(out$lambda_trace$updated)
})

test_that("CG scheduled lambda update delegates lambda-change decision", {
  lambda_current <- list(mu = list(`s(x)` = 1))
  lambda_updated <- list(mu = list(`s(x)` = 2))
  calls <- character()

  out <- gamlss.longitudinal:::.gl_maybe_update_cg_lambdas_on_schedule(
    lambda_current = lambda_current,
    lambda_trace = data.frame(),
    penalty_current = matrix(1, nrow = 1, dimnames = list("x", "x")),
    lambda_update_count = 1L,
    update_lambda = TRUE,
    max_lambda_updates = 5L,
    lambda_update_every = 2L,
    outer_iteration = 4L,
    H_obs_current = matrix(-1, nrow = 1),
    beta_vec = c(x = 1),
    grad_vec = c(x = 1),
    mm_cg = NULL,
    trust_radius = 10,
    max_delta = 10,
    step_tol = 0.5,
    build_penalty_fn = function(beta_names, lambda_try) {
      matrix(lambda_try$mu$`s(x)`, nrow = 1, dimnames = list(beta_names, beta_names))
    },
    eval_fn = function(...) NULL,
    edf_fn = function(...) list(),
    objective_fn = gamlss.longitudinal:::.gl_cg_objective,
    lambda_penalty_K = 2,
    verbose = 0,
    update_lambdas_fn = function(...) list(lambda = lambda_updated, lambda_trace = data.frame(updated = TRUE)),
    shrink_trust_radius_fn = function(trust_radius, step_tol) {
      calls <<- c(calls, "shrink")
      trust_radius / 2
    },
    lambdas_changed_fn = function(lambda_before, lambda_after) {
      calls <<- c(calls, "changed")
      expect_equal(lambda_before, lambda_current)
      expect_equal(lambda_after, lambda_updated)
      FALSE
    }
  )

  expect_equal(calls, "changed")
  expect_equal(out$lambda, lambda_updated)
  expect_equal(out$trust_radius, 10)
  expect_false(out$lambda_changed)
})

test_that("CG line-search evaluator chooses best improving step", {
  candidate_steps <- list(c(x = 0.2), c(x = 0.8))
  beta <- c(x = 0)
  penalty <- matrix(0, nrow = 1, dimnames = list("x", "x"))
  eval_candidate <- function(beta_try, mm_cg) {
    list(loglik = -as.numeric((beta_try["x"] - 0.7)^2))
  }

  out <- gamlss.longitudinal:::.gl_evaluate_cg_line_search(
    candidate_steps = candidate_steps,
    beta_vec = beta,
    mm_cg = NULL,
    penalty_current = penalty,
    obj_start = -1,
    trust_radius = Inf,
    max_delta = Inf,
    armijo_c1 = 0,
    line_search = "best",
    max_line_search_evals = 10,
    max_backtrack = 0,
    eval_fn = eval_candidate,
    objective_fn = gamlss.longitudinal:::.gl_cg_objective
  )

  expect_equal(out$line_eval_count, 2L)
  expect_equal(out$best$beta, c(x = 0.8))
  expect_true(out$best$improvement > 0)
})

test_that("CG line-search candidate helper evaluates one limited step", {
  penalty <- matrix(0, nrow = 2, ncol = 2, dimnames = list(c("x", "y"), c("x", "y")))
  captured <- list()

  out <- gamlss.longitudinal:::.gl_evaluate_cg_line_search_candidate(
    delta0 = c(x = 3, y = 4),
    backtrack_index = 2L,
    beta_vec = c(x = 1, y = 2),
    mm_cg = list(mm = TRUE),
    penalty_current = penalty,
    obj_start = -10,
    trust_radius = 2.5,
    max_delta = Inf,
    eval_fn = function(beta_try, mm_cg) {
      captured$beta <<- beta_try
      captured$mm <<- mm_cg
      list(loglik = -4)
    },
    objective_fn = gamlss.longitudinal:::.gl_cg_objective
  )

  expect_equal(captured$mm, list(mm = TRUE))
  expect_equal(out$beta, c(x = 2.5, y = 4))
  expect_equal(out$eval, list(loglik = -4))
  expect_equal(out$improvement, 6)
  expect_equal(out$step_l2, 2.5, tolerance = 1e-12)
})

test_that("CG line-search candidate helper rejects invalid likelihood states", {
  penalty <- matrix(0, nrow = 1, ncol = 1, dimnames = list("x", "x"))

  expect_null(gamlss.longitudinal:::.gl_evaluate_cg_line_search_candidate(
    delta0 = c(x = 1),
    backtrack_index = 1L,
    beta_vec = c(x = 0),
    mm_cg = NULL,
    penalty_current = penalty,
    obj_start = 0,
    trust_radius = Inf,
    max_delta = Inf,
    eval_fn = function(...) NULL,
    objective_fn = gamlss.longitudinal:::.gl_cg_objective
  ))

  expect_null(gamlss.longitudinal:::.gl_evaluate_cg_line_search_candidate(
    delta0 = c(x = 1),
    backtrack_index = 1L,
    beta_vec = c(x = 0),
    mm_cg = NULL,
    penalty_current = penalty,
    obj_start = 0,
    trust_radius = Inf,
    max_delta = Inf,
    eval_fn = function(...) list(loglik = NA_real_),
    objective_fn = gamlss.longitudinal:::.gl_cg_objective
  ))
})

test_that("CG line-search evaluator delegates candidate evaluation", {
  calls <- list()

  out <- gamlss.longitudinal:::.gl_evaluate_cg_line_search(
    candidate_steps = list(c(x = 1)),
    beta_vec = c(x = 0),
    mm_cg = list(mm = TRUE),
    penalty_current = matrix(0, nrow = 1, dimnames = list("x", "x")),
    obj_start = -1,
    trust_radius = Inf,
    max_delta = Inf,
    armijo_c1 = 0,
    line_search = "best",
    max_line_search_evals = 2L,
    max_backtrack = 1L,
    eval_fn = function(...) NULL,
    objective_fn = gamlss.longitudinal:::.gl_cg_objective,
    candidate_eval_fn = function(...) {
      calls[[length(calls) + 1L]] <<- list(...)
      list(beta = c(x = 1), eval = list(loglik = 0), improvement = 1, step_l2 = 1)
    }
  )

  expect_equal(out$line_eval_count, 2L)
  expect_equal(length(calls), 2L)
  expect_equal(calls[[1]]$backtrack_index, 1L)
  expect_equal(calls[[2]]$backtrack_index, 2L)
  expect_equal(calls[[1]]$mm_cg, list(mm = TRUE))
  expect_equal(out$best$beta, c(x = 1))
})

test_that("CG line-search evaluator can stop at first improving step", {
  candidate_steps <- list(c(x = 0.2), c(x = 0.8))
  beta <- c(x = 0)
  penalty <- matrix(0, nrow = 1, dimnames = list("x", "x"))
  eval_candidate <- function(beta_try, mm_cg) {
    list(loglik = -as.numeric((beta_try["x"] - 0.7)^2))
  }

  out <- gamlss.longitudinal:::.gl_evaluate_cg_line_search(
    candidate_steps = candidate_steps,
    beta_vec = beta,
    mm_cg = NULL,
    penalty_current = penalty,
    obj_start = -1,
    trust_radius = Inf,
    max_delta = Inf,
    armijo_c1 = 0,
    line_search = "first",
    max_line_search_evals = 10,
    max_backtrack = 0,
    eval_fn = eval_candidate,
    objective_fn = gamlss.longitudinal:::.gl_cg_objective
  )

  expect_equal(out$line_eval_count, 1L)
  expect_equal(out$best$beta, c(x = 0.2))
})

test_that("CG line-search evaluator respects backtracking and evaluation limit", {
  candidate_steps <- list(c(x = 4))
  beta <- c(x = 0)
  penalty <- matrix(0, nrow = 1, dimnames = list("x", "x"))
  eval_candidate <- function(beta_try, mm_cg) {
    list(loglik = -as.numeric((beta_try["x"] - 1)^2))
  }

  out <- gamlss.longitudinal:::.gl_evaluate_cg_line_search(
    candidate_steps = candidate_steps,
    beta_vec = beta,
    mm_cg = NULL,
    penalty_current = penalty,
    obj_start = -10,
    trust_radius = Inf,
    max_delta = Inf,
    armijo_c1 = 0,
    line_search = "best",
    max_line_search_evals = 2,
    max_backtrack = 4,
    eval_fn = eval_candidate,
    objective_fn = gamlss.longitudinal:::.gl_cg_objective
  )

  expect_equal(out$line_eval_count, 2L)
  expect_equal(out$best$beta, c(x = 2))
})

test_that("CG line-search evaluator returns NULL best for invalid candidates", {
  candidate_steps <- list(c(x = 1), c(x = 2))
  beta <- c(x = 0)
  penalty <- matrix(0, nrow = 1, dimnames = list("x", "x"))
  eval_invalid <- function(beta_try, mm_cg) NULL

  out <- gamlss.longitudinal:::.gl_evaluate_cg_line_search(
    candidate_steps = candidate_steps,
    beta_vec = beta,
    mm_cg = NULL,
    penalty_current = penalty,
    obj_start = 0,
    trust_radius = Inf,
    max_delta = Inf,
    armijo_c1 = 0,
    line_search = "best",
    max_line_search_evals = 10,
    max_backtrack = 0,
    eval_fn = eval_invalid,
    objective_fn = gamlss.longitudinal:::.gl_cg_objective
  )

  expect_null(out$best)
  expect_equal(out$line_eval_count, 2L)
})

test_that("CG line-search runner forwards backtracking controls", {
  out <- gamlss.longitudinal:::.gl_run_cg_line_search(
    candidate_steps = list(c(x = 1)),
    beta_vec = c(x = 0),
    mm_cg = list(mm = TRUE),
    penalty_current = matrix(0, nrow = 1),
    obj_start = -1,
    trust_radius = 2,
    max_delta = 3,
    armijo_c1 = 0.1,
    line_search = "best",
    max_line_search_evals = 4L,
    use_backtracking = TRUE,
    backtracking_max_halves = 5L,
    eval_fn = function(...) NULL,
    objective_fn = gamlss.longitudinal:::.gl_cg_objective,
    verbose = 0,
    line_search_fn = function(candidate_steps, beta_vec, mm_cg, penalty_current,
                              obj_start, trust_radius, max_delta, armijo_c1,
                              line_search, max_line_search_evals, max_backtrack,
                              eval_fn, objective_fn) {
      expect_equal(candidate_steps, list(c(x = 1)))
      expect_equal(beta_vec, c(x = 0))
      expect_equal(mm_cg, list(mm = TRUE))
      expect_equal(obj_start, -1)
      expect_equal(trust_radius, 2)
      expect_equal(max_delta, 3)
      expect_equal(armijo_c1, 0.1)
      expect_equal(line_search, "best")
      expect_equal(max_line_search_evals, 4L)
      expect_equal(max_backtrack, 5L)
      list(best = list(step_l2 = 1), line_eval_count = 9L)
    }
  )

  expect_equal(out$line_eval_count, 9L)
})

test_that("CG line-search runner disables backtracking when requested", {
  out <- gamlss.longitudinal:::.gl_run_cg_line_search(
    candidate_steps = list(c(x = 0.8)),
    beta_vec = c(x = 0),
    mm_cg = NULL,
    penalty_current = matrix(0, nrow = 1, dimnames = list("x", "x")),
    obj_start = -1,
    trust_radius = Inf,
    max_delta = Inf,
    armijo_c1 = 0,
    line_search = "best",
    max_line_search_evals = 10L,
    use_backtracking = FALSE,
    backtracking_max_halves = 5L,
    eval_fn = function(beta_try, mm_cg) {
      list(loglik = -as.numeric((beta_try["x"] - 0.7)^2))
    },
    objective_fn = gamlss.longitudinal:::.gl_cg_objective,
    verbose = 0
  )

  expect_equal(out$line_eval_count, 1L)
  expect_equal(out$best$beta, c(x = 0.8))
})

test_that("CG step acceptance helper handles missing best step", {
  out <- gamlss.longitudinal:::.gl_assess_cg_step_acceptance(
    best = NULL,
    best_raw_loglik = 10,
    outer_start_loglik = 9,
    raw_loglik_drop_tol = 1,
    lambda_update_count = 1L
  )

  expect_false(out$has_step)
  expect_false(out$accept)
  expect_false(out$prevented_deterioration)
  expect_true(is.na(out$prevented_raw_loglik_drop))
})

test_that("CG step acceptance helper accepts improving or tolerated steps", {
  best <- list(eval = list(loglik = 9.5))

  out <- gamlss.longitudinal:::.gl_assess_cg_step_acceptance(
    best = best,
    best_raw_loglik = 10,
    outer_start_loglik = 9,
    raw_loglik_drop_tol = 1,
    lambda_update_count = 1L
  )

  expect_true(out$has_step)
  expect_true(out$accept)
  expect_false(out$prevented_deterioration)
  expect_true(is.na(out$prevented_raw_loglik_drop))
})

test_that("CG step acceptance helper prevents raw log-likelihood deterioration after lambda updates", {
  best <- list(eval = list(loglik = 8))

  out <- gamlss.longitudinal:::.gl_assess_cg_step_acceptance(
    best = best,
    best_raw_loglik = 10,
    outer_start_loglik = 9,
    raw_loglik_drop_tol = 1,
    lambda_update_count = 1L
  )

  expect_true(out$has_step)
  expect_false(out$accept)
  expect_true(out$prevented_deterioration)
  expect_equal(out$prevented_raw_loglik_drop, 2)
})

test_that("CG step acceptance helper disables deterioration guard before lambda updates", {
  best <- list(eval = list(loglik = 8))

  out <- gamlss.longitudinal:::.gl_assess_cg_step_acceptance(
    best = best,
    best_raw_loglik = 10,
    outer_start_loglik = 9,
    raw_loglik_drop_tol = 1,
    lambda_update_count = 0L
  )

  expect_true(out$accept)
  expect_false(out$prevented_deterioration)
})

test_that("CG step application handles missing line-search step", {
  eval_start <- list(calc_lik = list(log_lik = c(joint = -10)))
  beta <- c(`mu.(Intercept)` = 1)
  par_s <- list(mu = list())
  step_acceptance <- list(
    has_step = FALSE,
    accept = FALSE,
    prevented_deterioration = FALSE,
    prevented_raw_loglik_drop = NA_real_
  )

  out <- gamlss.longitudinal:::.gl_apply_cg_step_acceptance(
    step_acceptance = step_acceptance,
    best = NULL,
    eval_start = eval_start,
    beta_vec = beta,
    par_cov_template = beta,
    par_s_template = par_s,
    stall_count = 2L,
    trust_radius = 4,
    step_tol = 0.5,
    max_delta = 10,
    max_stall = 5L,
    verbose = 0
  )

  expect_equal(out$beta, beta)
  expect_equal(out$par_cov, beta)
  expect_equal(out$par_s, par_s)
  expect_equal(out$calc_lik_out_end, eval_start$calc_lik)
  expect_equal(out$stall_count, 3L)
  expect_equal(out$trust_radius, 2)
  expect_true(is.na(out$accepted_improvement))
  expect_null(out$best)
  expect_false(out$prevented_deterioration)
})

test_that("CG step application preserves current state when deterioration is prevented", {
  eval_start <- list(calc_lik = list(log_lik = c(joint = -10)))
  beta <- c(`mu.(Intercept)` = 1)
  par_s <- list(mu = list())
  best <- list(
    beta = c(`mu.(Intercept)` = 2),
    eval = list(calc_lik = list(log_lik = c(joint = -20))),
    improvement = 1,
    step_l2 = 1
  )
  step_acceptance <- list(
    has_step = TRUE,
    accept = FALSE,
    prevented_deterioration = TRUE,
    prevented_raw_loglik_drop = 4
  )

  out <- gamlss.longitudinal:::.gl_apply_cg_step_acceptance(
    step_acceptance = step_acceptance,
    best = best,
    eval_start = eval_start,
    beta_vec = beta,
    par_cov_template = beta,
    par_s_template = par_s,
    stall_count = 2L,
    trust_radius = 4,
    step_tol = 0.5,
    max_delta = 10,
    max_stall = 5L,
    verbose = 0
  )

  expect_equal(out$beta, beta)
  expect_equal(out$par_cov, beta)
  expect_equal(out$par_s, par_s)
  expect_equal(out$calc_lik_out_end, eval_start$calc_lik)
  expect_equal(out$stall_count, 2L)
  expect_equal(out$trust_radius, 4)
  expect_true(is.na(out$accepted_improvement))
  expect_null(out$best)
  expect_true(out$prevented_deterioration)
  expect_equal(out$prevented_raw_loglik_drop, 4)
})

test_that("CG step application accepts step and unpacks coefficient state", {
  eval_start <- list(calc_lik = list(log_lik = c(joint = -10)))
  beta <- c(`mu.(Intercept)` = 1, `mu.s(x).1` = 0)
  par_cov <- c(`mu.(Intercept)` = 1)
  par_s <- list(mu = list(`s(x)` = c(`mu.s(x).1` = 0)))
  best <- list(
    beta = c(`mu.(Intercept)` = 2, `mu.s(x).1` = 3),
    eval = list(calc_lik = list(log_lik = c(joint = -5))),
    improvement = 1.5,
    step_l2 = 4
  )
  step_acceptance <- list(
    has_step = TRUE,
    accept = TRUE,
    prevented_deterioration = FALSE,
    prevented_raw_loglik_drop = NA_real_
  )

  out <- gamlss.longitudinal:::.gl_apply_cg_step_acceptance(
    step_acceptance = step_acceptance,
    best = best,
    eval_start = eval_start,
    beta_vec = beta,
    par_cov_template = par_cov,
    par_s_template = par_s,
    stall_count = 2L,
    trust_radius = 4,
    step_tol = 0.5,
    max_delta = 10,
    max_stall = 5L,
    verbose = 0
  )

  expect_equal(out$beta, best$beta)
  expect_equal(out$par_cov, best$beta[names(par_cov)])
  expect_equal(out$par_s$mu$`s(x)`, best$beta[names(par_s$mu$`s(x)`)])
  expect_equal(out$calc_lik_out_end, best$eval$calc_lik)
  expect_equal(out$stall_count, 0L)
  expect_equal(out$trust_radius, 6)
  expect_equal(out$accepted_improvement, 1.5)
  expect_equal(out$best, best)
  expect_false(out$prevented_deterioration)
})

test_that("CG line-search diagnostics helper sequences acceptance apply and diagnostics", {
  calls <- character()
  captured <- list()
  best <- list(id = "best")
  eval_start <- list(calc_lik = list(log_lik = c(joint = -10)))
  beta <- c(`mu.(Intercept)` = 1)
  par_s <- list(mu = list())
  step_acceptance <- list(
    has_step = TRUE,
    accept = TRUE,
    prevented_deterioration = FALSE,
    prevented_raw_loglik_drop = NA_real_
  )
  step_state <- list(
    beta = c(`mu.(Intercept)` = 2),
    par_cov = c(`mu.(Intercept)` = 2),
    par_s = par_s,
    calc_lik_out_end = list(log_lik = c(joint = -4)),
    stall_count = 0L,
    trust_radius = 0.75,
    accepted_improvement = 1.25,
    best = list(step_l2 = 0.2),
    prevented_deterioration = FALSE,
    prevented_raw_loglik_drop = NA_real_
  )
  diagnostics <- list(
    outer_log_lik_change = 6,
    best_raw_loglik = -4,
    best_iteration = 5L,
    raw_loglik_drop_from_best = 0,
    grad_inf = 0.05,
    step_l2 = 0.2,
    stopping = list(
      tolerance_met = TRUE,
      max_stall_hit = FALSE,
      deterioration_hit = FALSE,
      stop_requested = TRUE
    ),
    step_trace = list(data.frame(iter = 5L))
  )

  out <- gamlss.longitudinal:::.gl_apply_cg_line_search_diagnostics_state(
    best = best,
    eval_start = eval_start,
    beta_vec = beta,
    par_cov_template = beta,
    par_s_template = par_s,
    stall_count = 1L,
    trust_radius = 0.5,
    step_tol = 0.01,
    max_delta = 1,
    max_stall = 3L,
    verbose = 2,
    best_raw_loglik = -5,
    best_iteration = 4L,
    outer_iteration = 5L,
    outer_start_log_lik = -10,
    obj_start = -12,
    raw_loglik_drop_tol = 2,
    lambda_update_count = 1L,
    step_trace = list(),
    g_pen = c(`mu.(Intercept)` = 0.1),
    trust_radius_start = 0.5,
    line_eval_count = 4L,
    lambda_changed = TRUE,
    outer_stop_crit = 0.1,
    grad_tol = 0.2,
    assess_step_fn = function(best, best_raw_loglik, outer_start_loglik,
                              raw_loglik_drop_tol, lambda_update_count) {
      calls <<- c(calls, "assess")
      captured$assess_best <<- best
      captured$assess_best_raw_loglik <<- best_raw_loglik
      captured$assess_outer_start <<- outer_start_loglik
      captured$assess_drop_tol <<- raw_loglik_drop_tol
      captured$assess_lambda_count <<- lambda_update_count
      step_acceptance
    },
    apply_step_fn = function(step_acceptance, best, eval_start, beta_vec,
                             par_cov_template, par_s_template, stall_count,
                             trust_radius, step_tol, max_delta, max_stall,
                             verbose) {
      calls <<- c(calls, "apply")
      captured$apply_acceptance <<- step_acceptance
      captured$apply_best <<- best
      captured$apply_eval_start <<- eval_start
      captured$apply_beta <<- beta_vec
      captured$apply_stall_count <<- stall_count
      captured$apply_trust_radius <<- trust_radius
      captured$apply_step_tol <<- step_tol
      captured$apply_max_delta <<- max_delta
      captured$apply_max_stall <<- max_stall
      captured$apply_verbose <<- verbose
      step_state
    },
    diagnostics_fn = function(step_trace, g_pen, best, best_raw_loglik,
                              best_iteration, outer_iteration,
                              outer_start_log_lik, outer_end_log_lik,
                              obj_start, accepted_improvement,
                              trust_radius_start, trust_radius_end,
                              line_eval_count, lambda_update_count,
                              lambda_changed, stall_count,
                              prevented_deterioration,
                              prevented_raw_loglik_drop, max_stall,
                              raw_loglik_drop_tol, outer_stop_crit,
                              grad_tol, step_tol, verbose) {
      calls <<- c(calls, "diagnostics")
      captured$diagnostics_best <<- best
      captured$diagnostics_outer_end <<- outer_end_log_lik
      captured$diagnostics_accepted_improvement <<- accepted_improvement
      captured$diagnostics_trust_radius_end <<- trust_radius_end
      captured$diagnostics_line_eval_count <<- line_eval_count
      captured$diagnostics_lambda_changed <<- lambda_changed
      captured$diagnostics_stall_count <<- stall_count
      captured$diagnostics_prevented <<- prevented_deterioration
      captured$diagnostics_grad_tol <<- grad_tol
      captured$diagnostics_step_tol <<- step_tol
      diagnostics
    },
    state_builder_fn = function(...) {
      calls <<- c(calls, "state")
      gamlss.longitudinal:::.gl_build_cg_line_search_diagnostics_state(...)
    }
  )

  expect_equal(calls, c("assess", "apply", "diagnostics", "state"))
  expect_identical(captured$assess_best, best)
  expect_equal(captured$assess_best_raw_loglik, -5)
  expect_equal(captured$assess_outer_start, -10)
  expect_equal(captured$assess_lambda_count, 1L)
  expect_identical(captured$apply_acceptance, step_acceptance)
  expect_identical(captured$apply_best, best)
  expect_equal(captured$apply_beta, beta)
  expect_equal(captured$apply_trust_radius, 0.5)
  expect_equal(captured$diagnostics_best, step_state$best)
  expect_equal(captured$diagnostics_outer_end, -4)
  expect_equal(captured$diagnostics_accepted_improvement, 1.25)
  expect_equal(captured$diagnostics_trust_radius_end, 0.75)
  expect_true(captured$diagnostics_lambda_changed)
  expect_equal(captured$diagnostics_grad_tol, 0.2)
  expect_equal(captured$diagnostics_step_tol, 0.01)
  expect_identical(out$step_acceptance, step_acceptance)
  expect_identical(out$step_state, step_state)
  expect_equal(out$beta, step_state$beta)
  expect_equal(out$par_cov, step_state$par_cov)
  expect_identical(out$calc_lik_out_end, step_state$calc_lik_out_end)
  expect_equal(out$outer_end_log_lik, -4)
  expect_identical(out$iteration_diagnostics, diagnostics)
  expect_equal(out$outer_log_lik_change, 6)
  expect_equal(out$best_raw_loglik, -4)
  expect_equal(out$grad_inf, 0.05)
  expect_true(out$tolerance_met)
  expect_false(out$max_stall_hit)
  expect_false(out$deterioration_hit)
  expect_true(out$stop_requested)
  expect_identical(out$step_trace, diagnostics$step_trace)
})

test_that("CG line-search diagnostics state builder preserves returned-state contract", {
  step_acceptance <- list(
    has_step = TRUE,
    accept = TRUE,
    prevented_deterioration = FALSE,
    prevented_raw_loglik_drop = NA_real_
  )
  step_state <- list(
    beta = c(`mu.(Intercept)` = 2),
    par_cov = c(`mu.(Intercept)` = 2),
    par_s = list(mu = list()),
    calc_lik_out_end = list(log_lik = c(joint = -4)),
    stall_count = 0L,
    trust_radius = 0.75,
    accepted_improvement = 1.25,
    best = list(step_l2 = 0.2),
    prevented_deterioration = FALSE,
    prevented_raw_loglik_drop = NA_real_
  )
  diagnostics <- list(
    outer_log_lik_change = 6,
    best_raw_loglik = -4,
    best_iteration = 5L,
    raw_loglik_drop_from_best = 0,
    grad_inf = 0.05,
    step_l2 = 0.2,
    stopping = list(
      tolerance_met = TRUE,
      max_stall_hit = FALSE,
      deterioration_hit = FALSE,
      stop_requested = TRUE
    ),
    step_trace = list(data.frame(iter = 5L))
  )

  out <- gamlss.longitudinal:::.gl_build_cg_line_search_diagnostics_state(
    step_acceptance = step_acceptance,
    step_state = step_state,
    outer_end_log_lik = -4,
    iteration_diagnostics = diagnostics
  )

  expect_identical(out$step_acceptance, step_acceptance)
  expect_identical(out$step_state, step_state)
  expect_equal(out$beta, step_state$beta)
  expect_equal(out$par_cov, step_state$par_cov)
  expect_identical(out$calc_lik_out_end, step_state$calc_lik_out_end)
  expect_equal(out$outer_end_log_lik, -4)
  expect_identical(out$iteration_diagnostics, diagnostics)
  expect_equal(out$outer_log_lik_change, 6)
  expect_equal(out$best_raw_loglik, -4)
  expect_equal(out$best_iteration, 5L)
  expect_equal(out$raw_loglik_drop_from_best, 0)
  expect_equal(out$grad_inf, 0.05)
  expect_equal(out$step_l2, 0.2)
  expect_identical(out$stopping, diagnostics$stopping)
  expect_true(out$tolerance_met)
  expect_false(out$max_stall_hit)
  expect_false(out$deterioration_hit)
  expect_true(out$stop_requested)
  expect_identical(out$step_trace, diagnostics$step_trace)
})

test_that("CG stopping helper detects tolerance convergence", {
  out <- gamlss.longitudinal:::.gl_assess_cg_stopping(
    outer_log_lik_change = 1e-7,
    grad_inf = 1e-7,
    step_l2 = 1e-7,
    stall_count = 0L,
    max_stall = 3L,
    raw_loglik_drop_from_best = 0,
    raw_loglik_drop_tol = Inf,
    lambda_update_count = 0L,
    prevented_deterioration = FALSE,
    outer_stop_crit = 1e-6,
    grad_tol = 1e-6,
    step_tol = 1e-6
  )

  expect_true(out$tolerance_met)
  expect_false(out$max_stall_hit)
  expect_false(out$deterioration_hit)
  expect_true(out$stop_requested)
})

test_that("CG stopping helper detects max stall", {
  out <- gamlss.longitudinal:::.gl_assess_cg_stopping(
    outer_log_lik_change = 1,
    grad_inf = 1,
    step_l2 = 1,
    stall_count = 3L,
    max_stall = 3L,
    raw_loglik_drop_from_best = 0,
    raw_loglik_drop_tol = Inf,
    lambda_update_count = 0L,
    prevented_deterioration = FALSE,
    outer_stop_crit = 1e-6,
    grad_tol = 1e-6,
    step_tol = 1e-6
  )

  expect_false(out$tolerance_met)
  expect_true(out$max_stall_hit)
  expect_false(out$deterioration_hit)
  expect_true(out$stop_requested)
})

test_that("CG stopping helper detects deterioration after lambda updates", {
  out <- gamlss.longitudinal:::.gl_assess_cg_stopping(
    outer_log_lik_change = 1,
    grad_inf = 1,
    step_l2 = 1,
    stall_count = 0L,
    max_stall = 3L,
    raw_loglik_drop_from_best = 2,
    raw_loglik_drop_tol = 1,
    lambda_update_count = 1L,
    prevented_deterioration = FALSE,
    outer_stop_crit = 1e-6,
    grad_tol = 1e-6,
    step_tol = 1e-6
  )

  expect_false(out$tolerance_met)
  expect_false(out$max_stall_hit)
  expect_true(out$deterioration_hit)
  expect_true(out$stop_requested)
})

test_that("CG stopping helper carries prevented deterioration into stop decision", {
  out <- gamlss.longitudinal:::.gl_assess_cg_stopping(
    outer_log_lik_change = 1,
    grad_inf = 1,
    step_l2 = 1,
    stall_count = 0L,
    max_stall = 3L,
    raw_loglik_drop_from_best = 0,
    raw_loglik_drop_tol = 1,
    lambda_update_count = 0L,
    prevented_deterioration = TRUE,
    outer_stop_crit = 1e-6,
    grad_tol = 1e-6,
    step_tol = 1e-6
  )

  expect_true(out$deterioration_hit)
  expect_true(out$stop_requested)
})

test_that("CG stopping helper keeps iterating when no stop criterion is hit", {
  out <- gamlss.longitudinal:::.gl_assess_cg_stopping(
    outer_log_lik_change = 1,
    grad_inf = 1,
    step_l2 = 1,
    stall_count = 0L,
    max_stall = 3L,
    raw_loglik_drop_from_best = 0,
    raw_loglik_drop_tol = 1,
    lambda_update_count = 1L,
    prevented_deterioration = FALSE,
    outer_stop_crit = 1e-6,
    grad_tol = 1e-6,
    step_tol = 1e-6
  )

  expect_false(out$tolerance_met)
  expect_false(out$max_stall_hit)
  expect_false(out$deterioration_hit)
  expect_false(out$stop_requested)
})

test_that("CG iteration diagnostics update best likelihood and append trace", {
  out <- gamlss.longitudinal:::.gl_update_cg_iteration_diagnostics(
    step_trace = list(),
    g_pen = c(a = 0.01, b = -0.02),
    best = list(step_l2 = 0.03),
    best_raw_loglik = 10,
    best_iteration = 1L,
    outer_iteration = 2L,
    outer_start_log_lik = 10,
    outer_end_log_lik = 11,
    obj_start = 9,
    accepted_improvement = 0.5,
    trust_radius_start = 1,
    trust_radius_end = 1.2,
    line_eval_count = 3L,
    lambda_update_count = 0L,
    lambda_changed = FALSE,
    stall_count = 0L,
    prevented_deterioration = FALSE,
    prevented_raw_loglik_drop = NA_real_,
    max_stall = 3L,
    raw_loglik_drop_tol = Inf,
    outer_stop_crit = 1e-6,
    grad_tol = 1e-6,
    step_tol = 1e-6,
    verbose = 0
  )

  expect_equal(out$outer_log_lik_change, 1)
  expect_equal(out$best_raw_loglik, 11)
  expect_equal(out$best_iteration, 2L)
  expect_equal(out$raw_loglik_drop_from_best, 0)
  expect_equal(out$grad_inf, 0.02)
  expect_equal(out$step_l2, 0.03)
  expect_false(out$stopping$stop_requested)
  expect_length(out$step_trace, 1L)
  expect_equal(out$step_trace[[1]]$outer_iteration, 2L)
  expect_true(out$step_trace[[1]]$accepted_step)
  expect_equal(out$step_trace[[1]]$line_search_evals, 3L)
})

test_that("CG iteration diagnostics carry prevented deterioration into stopping", {
  out <- gamlss.longitudinal:::.gl_update_cg_iteration_diagnostics(
    step_trace = list(data.frame(existing = TRUE)),
    g_pen = c(a = 1),
    best = NULL,
    best_raw_loglik = 12,
    best_iteration = 1L,
    outer_iteration = 4L,
    outer_start_log_lik = 11,
    outer_end_log_lik = 10,
    obj_start = 10,
    accepted_improvement = NA_real_,
    trust_radius_start = 1,
    trust_radius_end = 0.5,
    line_eval_count = 0L,
    lambda_update_count = 0L,
    lambda_changed = FALSE,
    stall_count = 1L,
    prevented_deterioration = TRUE,
    prevented_raw_loglik_drop = 2,
    max_stall = 3L,
    raw_loglik_drop_tol = 1,
    outer_stop_crit = 1e-6,
    grad_tol = 1e-6,
    step_tol = 1e-6,
    verbose = 0
  )

  expect_equal(out$outer_log_lik_change, -1)
  expect_equal(out$best_raw_loglik, 12)
  expect_equal(out$best_iteration, 1L)
  expect_equal(out$raw_loglik_drop_from_best, 2)
  expect_equal(out$step_l2, 0)
  expect_true(out$stopping$deterioration_hit)
  expect_true(out$stopping$stop_requested)
  expect_length(out$step_trace, 2L)
  expect_false(out$step_trace[[2]]$accepted_step)
  expect_true(out$step_trace[[2]]$raw_deterioration_hit)
})

test_that("CG step trace row builder preserves reviewer-facing schema", {
  row <- gamlss.longitudinal:::.gl_build_cg_step_trace_row(
    outer_iteration = 2,
    start_loglik = 10.5,
    end_loglik = 11.5,
    raw_loglik_change = 1,
    start_penalized_loglik = 9,
    accepted_penalized_improvement = 0.5,
    grad_inf = 0.01,
    step_l2 = 0.02,
    trust_radius_start = 0.5,
    trust_radius_end = 0.75,
    line_search_evals = 4,
    accepted_step = TRUE,
    lambda_update_count = 1,
    lambda_changed = TRUE,
    stall_count = 0,
    tolerance_met = FALSE,
    max_stall_hit = FALSE,
    raw_deterioration_hit = FALSE,
    raw_loglik_drop_from_best = 0
  )

  expect_equal(nrow(row), 1L)
  expect_equal(
    names(row),
    c(
      "outer_iteration", "start_logLik", "end_logLik", "raw_logLik_change",
      "start_penalized_logLik", "accepted_penalized_improvement",
      "grad_inf", "step_l2", "trust_radius_start", "trust_radius_end",
      "line_search_evals", "accepted_step", "lambda_update_count",
      "lambda_changed", "stall_count", "tolerance_met", "max_stall_hit",
      "raw_deterioration_hit", "raw_loglik_drop_from_best"
    )
  )
  expect_type(row$outer_iteration, "integer")
  expect_type(row$line_search_evals, "integer")
  expect_type(row$accepted_step, "logical")
  expect_true(row$accepted_step)
  expect_true(row$lambda_changed)
})

test_that("CG step trace row builder normalizes logical flags", {
  row <- gamlss.longitudinal:::.gl_build_cg_step_trace_row(
    outer_iteration = 1,
    start_loglik = 1,
    end_loglik = 1,
    raw_loglik_change = 0,
    start_penalized_loglik = 1,
    accepted_penalized_improvement = NA_real_,
    grad_inf = 1,
    step_l2 = 0,
    trust_radius_start = 1,
    trust_radius_end = 0.5,
    line_search_evals = 0,
    accepted_step = NULL,
    lambda_update_count = 0,
    lambda_changed = NULL,
    stall_count = 1,
    tolerance_met = NULL,
    max_stall_hit = TRUE,
    raw_deterioration_hit = NULL,
    raw_loglik_drop_from_best = 0
  )

  expect_false(row$accepted_step)
  expect_false(row$lambda_changed)
  expect_false(row$tolerance_met)
  expect_true(row$max_stall_hit)
  expect_false(row$raw_deterioration_hit)
})

test_that("CG lambda trace row builder preserves candidate schema", {
  scores <- list(
    raw_loglik = c(1, 2, 3),
    penalty_value = c(0.1, 0.2, 0.3),
    penalized_loglik = c(0.9, 1.8, 2.7),
    edf_values = c(1, 1.5, 2),
    gaic_score = c(10, 9, 8)
  )

  rows <- gamlss.longitudinal:::.gl_build_cg_lambda_trace_rows(
    outer_iteration = 4,
    parameter = "mu",
    smooth = "s(x)",
    lambda_before = 2,
    candidates = c(1, 2, 4),
    lambda_scores = scores,
    best = 3L
  )

  expect_equal(nrow(rows), 3L)
  expect_equal(
    names(rows),
    c(
      "outer_iteration", "parameter", "smooth", "lambda_before",
      "lambda_candidate", "raw_logLik_after_step", "smooth_penalty_after_step",
      "penalized_logLik_after_step", "edf_after_step", "gaic_score", "chosen"
    )
  )
  expect_type(rows$outer_iteration, "integer")
  expect_equal(rows$parameter, rep("mu", 3))
  expect_equal(rows$smooth, rep("s(x)", 3))
  expect_equal(rows$lambda_before, rep(2, 3))
  expect_equal(rows$lambda_candidate, c(1, 2, 4))
  expect_equal(rows$chosen, c(FALSE, FALSE, TRUE))
})

test_that("CG stop reason helper preserves reason priority", {
  expect_equal(
    gamlss.longitudinal:::.gl_cg_stop_reason(
      tolerance_met = TRUE,
      deterioration_hit = TRUE
    ),
    "tolerance"
  )
  expect_equal(
    gamlss.longitudinal:::.gl_cg_stop_reason(
      tolerance_met = FALSE,
      deterioration_hit = TRUE
    ),
    "raw_loglik_deterioration"
  )
  expect_equal(
    gamlss.longitudinal:::.gl_cg_stop_reason(
      tolerance_met = FALSE,
      deterioration_hit = FALSE
    ),
    "max_stall"
  )
})

test_that("CG convergence delay helper requires first lambda update opportunity", {
  expect_true(
    gamlss.longitudinal:::.gl_should_delay_cg_convergence_for_lambda_update(
      update_lambda = TRUE,
      has_smooths = TRUE,
      lambda_update_count = 0L
    )
  )
  expect_false(
    gamlss.longitudinal:::.gl_should_delay_cg_convergence_for_lambda_update(
      update_lambda = FALSE,
      has_smooths = TRUE,
      lambda_update_count = 0L
    )
  )
  expect_false(
    gamlss.longitudinal:::.gl_should_delay_cg_convergence_for_lambda_update(
      update_lambda = TRUE,
      has_smooths = FALSE,
      lambda_update_count = 0L
    )
  )
  expect_false(
    gamlss.longitudinal:::.gl_should_delay_cg_convergence_for_lambda_update(
      update_lambda = TRUE,
      has_smooths = TRUE,
      lambda_update_count = 1L
    )
  )
})

test_that("CG stop handler delays convergence for first smoother lambda update", {
  lambda_current <- list(mu = list(`s(x)` = 1))
  lambda_trace <- data.frame()
  beta <- c(`mu.(Intercept)` = 1)
  H_obs <- matrix(-2, nrow = 1, dimnames = list(names(beta), names(beta)))
  grad <- c(`mu.(Intercept)` = 0.1)
  penalty <- matrix(2, nrow = 1, dimnames = list(names(beta), names(beta)))

  out <- gamlss.longitudinal:::.gl_handle_cg_stop_request(
    update_lambda = TRUE,
    has_smooths = TRUE,
    lambda_update_count = 0L,
    lambda_current = lambda_current,
    lambda_trace = lambda_trace,
    H_obs_current = H_obs,
    beta_vec = beta,
    grad_vec = grad,
    mm_cg = NULL,
    trust_radius = 1,
    outer_iteration = 3,
    max_delta = 10,
    build_penalty_fn = function(beta_names, lambda_new) {
      expect_equal(lambda_new$mu$`s(x)`, 2)
      penalty
    },
    eval_fn = function(...) NULL,
    edf_fn = function(H_obs_current, penalty_mat, beta_names) {
      expect_equal(penalty_mat, penalty)
      list(mu = list(`s(x)` = 0.5))
    },
    objective_fn = gamlss.longitudinal:::.gl_cg_objective,
    lambda_penalty_K = 2,
    tolerance_met = TRUE,
    deterioration_hit = FALSE,
    raw_loglik_drop_from_best = 0,
    verbose = 0,
    update_lambdas_fn = function(lambda_current, lambda_trace, ...) {
      list(lambda = list(mu = list(`s(x)` = 2)), lambda_trace = data.frame(updated = TRUE))
    }
  )

  expect_equal(out$lambda$mu$`s(x)`, 2)
  expect_equal(out$lambda_trace, data.frame(updated = TRUE))
  expect_equal(out$penalty_mat, penalty)
  expect_equal(out$df_s, list(mu = list(`s(x)` = 0.5)))
  expect_equal(out$lambda_update_count, 1L)
  expect_true(out$stall_count_reset)
  expect_identical(out$stop_reason, NA_character_)
  expect_false(out$converged)
})

test_that("CG stop handler marks final convergence reason when no delay is needed", {
  lambda_current <- list(mu = list())
  lambda_trace <- data.frame()
  beta <- c(`mu.(Intercept)` = 1)
  H_obs <- matrix(-2, nrow = 1, dimnames = list(names(beta), names(beta)))
  grad <- c(`mu.(Intercept)` = 0.1)

  out <- gamlss.longitudinal:::.gl_handle_cg_stop_request(
    update_lambda = FALSE,
    has_smooths = FALSE,
    lambda_update_count = 2L,
    lambda_current = lambda_current,
    lambda_trace = lambda_trace,
    H_obs_current = H_obs,
    beta_vec = beta,
    grad_vec = grad,
    mm_cg = NULL,
    trust_radius = 1,
    outer_iteration = 3,
    max_delta = 10,
    build_penalty_fn = function(...) stop("no delayed update should occur"),
    eval_fn = function(...) NULL,
    edf_fn = function(...) NULL,
    objective_fn = gamlss.longitudinal:::.gl_cg_objective,
    lambda_penalty_K = 2,
    tolerance_met = FALSE,
    deterioration_hit = TRUE,
    raw_loglik_drop_from_best = 4,
    verbose = 0
  )

  expect_equal(out$lambda, lambda_current)
  expect_equal(out$lambda_trace, lambda_trace)
  expect_null(out$penalty_mat)
  expect_null(out$df_s)
  expect_equal(out$lambda_update_count, 2L)
  expect_false(out$stall_count_reset)
  expect_equal(out$stop_reason, "raw_loglik_deterioration")
  expect_true(out$converged)
})

test_that("CG stop state helper applies delayed-update penalty EDF and stall reset", {
  penalty_current <- matrix(0, nrow = 1)
  penalty_new <- matrix(2, nrow = 1)
  df_current <- list(mu = list(old = 1))
  df_new <- list(mu = list(new = 2))
  lambda_new <- list(mu = list(`s(x)` = 3))
  lambda_trace_new <- data.frame(updated = TRUE)
  captured <- NULL

  out <- gamlss.longitudinal:::.gl_apply_cg_stop_request_state(
    update_lambda = TRUE,
    has_smooths = TRUE,
    lambda_update_count = 0L,
    lambda_current = list(mu = list(`s(x)` = 1)),
    lambda_trace = data.frame(),
    penalty_current = penalty_current,
    df_s_current = df_current,
    stall_count = 4L,
    H_obs_current = matrix(-2, nrow = 1),
    beta_vec = c(x = 1),
    grad_vec = c(x = 0.1),
    mm_cg = list(),
    trust_radius = 1,
    outer_iteration = 2L,
    max_delta = 5,
    build_penalty_fn = function(...) NULL,
    eval_fn = function(...) NULL,
    edf_fn = function(...) NULL,
    objective_fn = function(...) NULL,
    lambda_penalty_K = 2,
    tolerance_met = TRUE,
    deterioration_hit = FALSE,
    raw_loglik_drop_from_best = 0,
    verbose = 0,
    stop_request_fn = function(update_lambda, has_smooths, lambda_update_count,
                               lambda_current, lambda_trace, H_obs_current,
                               beta_vec, grad_vec, mm_cg, trust_radius,
                               outer_iteration, max_delta, build_penalty_fn,
                               eval_fn, edf_fn, objective_fn,
                               lambda_penalty_K, tolerance_met,
                               deterioration_hit, raw_loglik_drop_from_best,
                               verbose) {
      captured <<- list(
        update_lambda = update_lambda,
        has_smooths = has_smooths,
        lambda_update_count = lambda_update_count,
        lambda_current = lambda_current,
        outer_iteration = outer_iteration,
        tolerance_met = tolerance_met
      )
      list(
        lambda = lambda_new,
        lambda_trace = lambda_trace_new,
        penalty_mat = penalty_new,
        df_s = df_new,
        lambda_update_count = 1L,
        stall_count_reset = TRUE,
        stop_reason = NA_character_,
        converged = FALSE
      )
    }
  )

  expect_true(captured$update_lambda)
  expect_true(captured$has_smooths)
  expect_equal(captured$lambda_update_count, 0L)
  expect_equal(captured$outer_iteration, 2L)
  expect_identical(out$lambda, lambda_new)
  expect_equal(out$lambda_trace, lambda_trace_new)
  expect_equal(out$lambda_update_count, 1L)
  expect_equal(out$penalty_mat, penalty_new)
  expect_equal(out$df_s, df_new)
  expect_equal(out$stall_count, 0L)
  expect_identical(out$stop_reason, NA_character_)
  expect_false(out$converged)
})

test_that("CG stop state helper preserves current penalty EDF and stall count for final stop", {
  penalty_current <- matrix(0, nrow = 1)
  df_current <- list(mu = list(old = 1))
  lambda_current <- list(mu = list())
  lambda_trace <- data.frame()

  out <- gamlss.longitudinal:::.gl_apply_cg_stop_request_state(
    update_lambda = FALSE,
    has_smooths = FALSE,
    lambda_update_count = 3L,
    lambda_current = lambda_current,
    lambda_trace = lambda_trace,
    penalty_current = penalty_current,
    df_s_current = df_current,
    stall_count = 4L,
    H_obs_current = matrix(-2, nrow = 1),
    beta_vec = c(x = 1),
    grad_vec = c(x = 0.1),
    mm_cg = list(),
    trust_radius = 1,
    outer_iteration = 2L,
    max_delta = 5,
    build_penalty_fn = function(...) NULL,
    eval_fn = function(...) NULL,
    edf_fn = function(...) NULL,
    objective_fn = function(...) NULL,
    lambda_penalty_K = 2,
    tolerance_met = FALSE,
    deterioration_hit = TRUE,
    raw_loglik_drop_from_best = 4,
    verbose = 0,
    stop_request_fn = function(...) {
      list(
        lambda = lambda_current,
        lambda_trace = lambda_trace,
        penalty_mat = NULL,
        df_s = NULL,
        lambda_update_count = 3L,
        stall_count_reset = FALSE,
        stop_reason = "raw_loglik_deterioration",
        converged = TRUE
      )
    }
  )

  expect_identical(out$lambda, lambda_current)
  expect_equal(out$lambda_trace, lambda_trace)
  expect_equal(out$lambda_update_count, 3L)
  expect_equal(out$penalty_mat, penalty_current)
  expect_equal(out$df_s, df_current)
  expect_equal(out$stall_count, 4L)
  expect_equal(out$stop_reason, "raw_loglik_deterioration")
  expect_true(out$converged)
})

test_that("CG optional stop-state wrapper preserves state when no stop is requested", {
  lambda_s <- list(mu = list())
  lambda_trace <- data.frame()
  penalty_mat <- matrix(0, nrow = 1)
  df_s <- list(mu = list())

  out <- gamlss.longitudinal:::.gl_maybe_apply_cg_stop_request_state(
    stop_requested = FALSE,
    cg_update_lambda = TRUE,
    cg_has_smooths = TRUE,
    cg_lambda_update_count = 2L,
    lambda_s = lambda_s,
    cg_lambda_trace = lambda_trace,
    penalty_mat = penalty_mat,
    df_s = df_s,
    cg_stall_count = 4L,
    H_obs = matrix(-2, nrow = 1),
    beta_all = c(x = 1),
    grad = c(x = 0.1),
    mm_cg = list(),
    cg_trust_radius = 1,
    outer_only_run_counter = 3L,
    cg_max_delta = 5,
    build_cg_penalty = function(...) stop("stop helper should not be called"),
    cg_eval = function(...) NULL,
    cg_smooth_edf_list = function(...) NULL,
    cg_objective = function(...) NULL,
    lambda_penalty_K = 2,
    cg_tolerance_met = FALSE,
    cg_deterioration_hit = FALSE,
    cg_raw_loglik_drop_from_best = 0,
    verbose = 0,
    stop_request_state_fn = function(...) stop("stop helper should not be called")
  )

  expect_identical(out$lambda_s, lambda_s)
  expect_equal(out$cg_lambda_trace, lambda_trace)
  expect_equal(out$cg_lambda_update_count, 2L)
  expect_equal(out$penalty_mat, penalty_mat)
  expect_equal(out$df_s, df_s)
  expect_equal(out$cg_stall_count, 4L)
  expect_identical(out$cg_stop_reason, NA_character_)
  expect_false(out$cg_converged)
})

test_that("CG optional stop-state wrapper delegates requested stop", {
  captured <- NULL

  out <- gamlss.longitudinal:::.gl_maybe_apply_cg_stop_request_state(
    stop_requested = TRUE,
    cg_update_lambda = TRUE,
    cg_has_smooths = TRUE,
    cg_lambda_update_count = 0L,
    lambda_s = list(mu = list(`s(x)` = 1)),
    cg_lambda_trace = data.frame(),
    penalty_mat = matrix(0, nrow = 1),
    df_s = list(mu = list(old = 1)),
    cg_stall_count = 4L,
    H_obs = matrix(-2, nrow = 1),
    beta_all = c(x = 1),
    grad = c(x = 0.1),
    mm_cg = list(marker = TRUE),
    cg_trust_radius = 1,
    outer_only_run_counter = 3L,
    cg_max_delta = 5,
    build_cg_penalty = function(...) NULL,
    cg_eval = function(...) NULL,
    cg_smooth_edf_list = function(...) NULL,
    cg_objective = function(...) NULL,
    lambda_penalty_K = 2,
    cg_tolerance_met = TRUE,
    cg_deterioration_hit = FALSE,
    cg_raw_loglik_drop_from_best = 0,
    verbose = 0,
    stop_request_state_fn = function(...) {
      captured <<- list(...)
      list(
        lambda = list(mu = list(`s(x)` = 2)),
        lambda_trace = data.frame(updated = TRUE),
        lambda_update_count = 1L,
        penalty_mat = matrix(2, nrow = 1),
        df_s = list(mu = list(new = 2)),
        stall_count = 0L,
        stop_reason = NA_character_,
        converged = FALSE
      )
    }
  )

  expect_true(captured$update_lambda)
  expect_true(captured$has_smooths)
  expect_equal(captured$lambda_update_count, 0L)
  expect_equal(captured$outer_iteration, 3L)
  expect_equal(captured$tolerance_met, TRUE)
  expect_equal(out$lambda_s$mu$`s(x)`, 2)
  expect_equal(out$cg_lambda_trace, data.frame(updated = TRUE))
  expect_equal(out$cg_lambda_update_count, 1L)
  expect_equal(out$penalty_mat, matrix(2, nrow = 1))
  expect_equal(out$df_s, list(mu = list(new = 2)))
  expect_equal(out$cg_stall_count, 0L)
  expect_identical(out$cg_stop_reason, NA_character_)
  expect_false(out$cg_converged)
})

test_that("CG final smooth EDF refresh updates penalty and EDF when Hessian is available", {
  dataset <- data.frame(response = c(1, 2), time = c(1, 2), subject = c("a", "a"))
  beta <- c(`mu.(Intercept)` = 1)
  mm_cg <- list(x = list(mu = matrix(1, nrow = 2)), s = list(mu = list()))
  old_penalty <- matrix(0, nrow = 1, dimnames = list(names(beta), names(beta)))
  new_penalty <- matrix(2, nrow = 1, dimnames = list(names(beta), names(beta)))
  final_H <- matrix(-3, nrow = 1, dimnames = list(names(beta), names(beta)))

  out <- gamlss.longitudinal:::.gl_refresh_final_cg_smooth_edf(
    dataset = dataset,
    margin_dist = gamlss.dist::NO(),
    copula_dist = "N",
    mm_cg = mm_cg,
    beta_vec = beta,
    lambda_current = list(mu = list()),
    penalty_current = old_penalty,
    df_s_current = list(mu = list()),
    observed_hessian_fn = function(obj, beta_vec, mm_cg_arg, context) {
      expect_equal(context, "final smooth EDF update")
      expect_equal(obj$response, dataset$response)
      expect_equal(beta_vec, beta)
      expect_equal(mm_cg_arg, mm_cg)
      final_H
    },
    build_penalty_fn = function(beta_names, lambda_current) {
      expect_equal(beta_names, names(beta))
      new_penalty
    },
    edf_fn = function(H_obs_current, penalty_current, beta_names) {
      expect_equal(H_obs_current, final_H)
      expect_equal(penalty_current, new_penalty)
      list(mu = list(`s(x)` = 0.75))
    }
  )

  expect_equal(out$final_H, final_H)
  expect_equal(out$penalty_mat, new_penalty)
  expect_equal(out$df_s, list(mu = list(`s(x)` = 0.75)))
})

test_that("CG final smooth EDF refresh preserves current state when Hessian fails", {
  dataset <- data.frame(response = c(1, 2), time = c(1, 2), subject = c("a", "a"))
  beta <- c(`mu.(Intercept)` = 1)
  mm_cg <- list(x = list(mu = matrix(1, nrow = 2)), s = list(mu = list()))
  old_penalty <- matrix(0, nrow = 1, dimnames = list(names(beta), names(beta)))
  old_df <- list(mu = list(`s(x)` = 0.25))

  out <- gamlss.longitudinal:::.gl_refresh_final_cg_smooth_edf(
    dataset = dataset,
    margin_dist = gamlss.dist::NO(),
    copula_dist = "N",
    mm_cg = mm_cg,
    beta_vec = beta,
    lambda_current = list(mu = list()),
    penalty_current = old_penalty,
    df_s_current = old_df,
    observed_hessian_fn = function(...) stop("no final Hessian"),
    build_penalty_fn = function(...) stop("penalty should not be rebuilt"),
    edf_fn = function(...) stop("EDF should not be recomputed")
  )

  expect_null(out$final_H)
  expect_equal(out$penalty_mat, old_penalty)
  expect_equal(out$df_s, old_df)
})

test_that("CG optimizer finalization helper refreshes EDF and initializes weights", {
  captured <- list()
  dataset <- data.frame(response = c(1, 2), time = c(1, 2), subject = c("a", "a"))
  beta <- c(`mu.(Intercept)` = 1)
  mm <- list(x = list(mu = matrix(1, nrow = 2), sigma = matrix(1, nrow = 3)))
  mm_cg <- list(x = list(mu = matrix(1, nrow = 2)))
  old_penalty <- matrix(0, nrow = 1, ncol = 1)
  new_penalty <- matrix(2, nrow = 1, ncol = 1)
  old_df <- list(mu = list(old = 1))
  new_df <- list(mu = list(new = 2))
  final_edf <- list(final_H = matrix(-1, nrow = 1), penalty_mat = new_penalty, df_s = new_df)
  weights <- list(mu = c(1, 1), sigma = c(1, 1, 1))

  out <- gamlss.longitudinal:::.gl_finalize_cg_optimizer_state(
    dataset = dataset,
    margin_dist = "NO",
    copula_dist = "N",
    mm = mm,
    mm_cg = mm_cg,
    beta_vec = beta,
    lambda_current = list(mu = list()),
    penalty_current = old_penalty,
    df_s_current = old_df,
    observed_hessian_fn = function(...) NULL,
    build_penalty_fn = function(...) NULL,
    edf_fn = function(...) NULL,
    final_edf_fn = function(dataset, margin_dist, copula_dist, mm_cg, beta_vec,
                            lambda_current, penalty_current, df_s_current,
                            observed_hessian_fn, build_penalty_fn, edf_fn) {
      captured$final_edf <<- list(
        dataset = dataset,
        margin_dist = margin_dist,
        copula_dist = copula_dist,
        mm_cg = mm_cg,
        beta_vec = beta_vec,
        lambda_current = lambda_current,
        penalty_current = penalty_current,
        df_s_current = df_s_current,
        observed_hessian_fn = observed_hessian_fn,
        build_penalty_fn = build_penalty_fn,
        edf_fn = edf_fn
      )
      final_edf
    },
    weights_fn = function(mm) {
      captured$weights_mm <<- mm
      weights
    }
  )

  expect_identical(captured$final_edf$dataset, dataset)
  expect_identical(captured$final_edf$margin_dist, "NO")
  expect_identical(captured$final_edf$copula_dist, "N")
  expect_identical(captured$final_edf$mm_cg, mm_cg)
  expect_equal(captured$final_edf$beta_vec, beta)
  expect_equal(captured$final_edf$penalty_current, old_penalty)
  expect_equal(captured$final_edf$df_s_current, old_df)
  expect_identical(captured$weights_mm, mm)
  expect_identical(out$final_edf, final_edf)
  expect_equal(out$penalty_mat, new_penalty)
  expect_equal(out$df_s, new_df)
  expect_identical(out$weights_final, weights)
})

test_that("CG zeta Hessian override is inactive unless finite mode and zeta terms are present", {
  H_obs <- diag(2)
  dimnames(H_obs) <- list(c("mu.x", "theta.x"), c("mu.x", "theta.x"))
  beta <- c(mu.x = 1, theta.x = 2)

  analytical <- gamlss.longitudinal:::.gl_apply_cg_zeta_hessian_override(
    H_obs = H_obs,
    beta_vec = beta,
    mm_cg = NULL,
    zeta_hessian = "analytical",
    finite_hessian_fn = function(...) stop("finite path should not run"),
    verbose = 0
  )
  expect_equal(analytical$H_obs, H_obs)
  expect_null(analytical$H_zeta_fd)

  no_zeta <- gamlss.longitudinal:::.gl_apply_cg_zeta_hessian_override(
    H_obs = H_obs,
    beta_vec = beta,
    mm_cg = NULL,
    zeta_hessian = "finite",
    finite_hessian_fn = function(...) stop("finite path should not run without zeta"),
    verbose = 0
  )
  expect_equal(no_zeta$H_obs, H_obs)
  expect_null(no_zeta$H_zeta_fd)
})

test_that("CG zeta Hessian override replaces finite zeta block", {
  beta <- c(mu.x = 1, zeta.a = 2, zeta.b = 3)
  H_obs <- diag(3)
  dimnames(H_obs) <- list(names(beta), names(beta))
  H_zeta <- matrix(
    c(-2, -0.5, -1, -4),
    nrow = 2,
    dimnames = list(c("zeta.a", "zeta.b"), c("zeta.a", "zeta.b"))
  )

  out <- gamlss.longitudinal:::.gl_apply_cg_zeta_hessian_override(
    H_obs = H_obs,
    beta_vec = beta,
    mm_cg = list(marker = TRUE),
    zeta_hessian = "finite",
    finite_hessian_fn = function(beta_vec, block_names, mm_cg) {
      expect_equal(beta_vec, beta)
      expect_equal(block_names, c("zeta.a", "zeta.b"))
      expect_equal(mm_cg, list(marker = TRUE))
      H_zeta
    },
    verbose = 0
  )

  expected <- H_obs
  expected[c("zeta.a", "zeta.b"), c("zeta.a", "zeta.b")] <- 0.5 * (H_zeta + t(H_zeta))
  expect_equal(out$H_obs, expected)
  expect_equal(out$H_zeta_fd, H_zeta)
})

test_that("CG zeta Hessian override preserves analytical block when finite block is invalid", {
  beta <- c(mu.x = 1, zeta.a = 2)
  H_obs <- diag(2)
  dimnames(H_obs) <- list(names(beta), names(beta))
  H_zeta <- matrix(NA_real_, nrow = 1, dimnames = list("zeta.a", "zeta.a"))

  out <- gamlss.longitudinal:::.gl_apply_cg_zeta_hessian_override(
    H_obs = H_obs,
    beta_vec = beta,
    mm_cg = NULL,
    zeta_hessian = "finite",
    finite_hessian_fn = function(...) H_zeta,
    verbose = 0
  )

  expect_equal(out$H_obs, H_obs)
  expect_equal(out$H_zeta_fd, H_zeta)
})

test_that("CG finite gradient helper approximates quadratic score", {
  eval_quadratic <- function(beta_vec, mm_cg) {
    list(loglik = -(beta_vec["x"]^2 + 2 * beta_vec["y"]^2 + beta_vec["x"] * beta_vec["y"]))
  }
  beta <- c(x = 1.5, y = -0.75)
  base_ll <- eval_quadratic(beta, NULL)$loglik
  expected <- c(
    x = -2 * unname(beta["x"]) - unname(beta["y"]),
    y = -4 * unname(beta["y"]) - unname(beta["x"])
  )

  central <- gamlss.longitudinal:::.gl_cg_finite_gradient(
    beta_vec = beta,
    base_ll = base_ll,
    mm_cg = NULL,
    eval_fn = eval_quadratic,
    gradient_method = "central"
  )
  forward <- gamlss.longitudinal:::.gl_cg_finite_gradient(
    beta_vec = beta,
    base_ll = base_ll,
    mm_cg = NULL,
    eval_fn = eval_quadratic,
    gradient_method = "forward"
  )

  expect_equal(central, expected, tolerance = 1e-7)
  expect_equal(forward, expected, tolerance = 1e-4)
})

test_that("CG finite Hessian block helper approximates quadratic curvature", {
  eval_quadratic <- function(beta_vec, mm_cg) {
    list(loglik = -(beta_vec["x"]^2 + 2 * beta_vec["y"]^2 + beta_vec["x"] * beta_vec["y"]))
  }
  beta <- c(x = 1.5, y = -0.75, z = 2)
  expected <- matrix(
    c(-2, -1, -1, -4),
    nrow = 2,
    dimnames = list(c("x", "y"), c("x", "y"))
  )

  H <- gamlss.longitudinal:::.gl_cg_finite_hessian_block(
    beta_vec = beta,
    block_names = c("x", "y"),
    mm_cg = NULL,
    eval_fn = eval_quadratic
  )

  expect_equal(H, expected, tolerance = 1e-5)
})

test_that("CG Hessian checker validates shape and finite entries", {
  H <- diag(2)

  expect_true(gamlss.longitudinal:::.gl_cg_hessian_ok(H, c("x", "y")))
  expect_false(gamlss.longitudinal:::.gl_cg_hessian_ok(H, c("x", "y", "z")))
  H[1, 1] <- NA_real_
  expect_false(gamlss.longitudinal:::.gl_cg_hessian_ok(H, c("x", "y")))
})

test_that("CG observed Hessian helper uses valid analytical Hessian", {
  beta <- c(x = 1, y = 2)
  H <- matrix(c(-2, -0.5, -1, -4), nrow = 2, dimnames = list(names(beta), names(beta)))
  finite_fn <- function(beta_vec, beta_names, mm_cg) stop("finite path should not be used")

  out <- gamlss.longitudinal:::.gl_cg_observed_hessian(
    tmp_obj = list(),
    beta_vec = beta,
    mm_cg = NULL,
    analytical_fn = function(obj) H,
    finite_fn = finite_fn,
    hessian_method = "analytical",
    verbose = 0
  )

  expect_equal(out, 0.5 * (H + t(H)))
})

test_that("CG observed Hessian helper falls back to finite Hessian", {
  beta <- c(x = 1, y = 2)
  H_fd <- matrix(c(-2, -1, -1, -4), nrow = 2, dimnames = list(names(beta), names(beta)))
  finite_fn <- function(beta_vec, beta_names, mm_cg) {
    expect_equal(beta_names, names(beta))
    H_fd
  }

  expect_warning(
    out <- gamlss.longitudinal:::.gl_cg_observed_hessian(
      tmp_obj = list(),
      beta_vec = beta,
      mm_cg = NULL,
      analytical_fn = function(obj) matrix(NA_real_, 2, 2),
      finite_fn = finite_fn,
      hessian_method = "analytical",
      verbose = 1,
      context = "unit test"
    ),
    "falling back to finite-difference Hessian"
  )
  expect_equal(out, H_fd)
})

test_that("CG observed Hessian helper supports explicit finite mode", {
  beta <- c(x = 1, y = 2)
  H_fd <- matrix(c(-2, -1, -1, -4), nrow = 2, dimnames = list(names(beta), names(beta)))
  finite_fn <- function(beta_vec, beta_names, mm_cg) H_fd

  out <- gamlss.longitudinal:::.gl_cg_observed_hessian(
    tmp_obj = list(),
    beta_vec = beta,
    mm_cg = NULL,
    analytical_fn = function(obj) stop("analytical path should not be used"),
    finite_fn = finite_fn,
    hessian_method = "finite",
    verbose = 1
  )

  expect_equal(out, H_fd)
})

test_that("CG observed Hessian helper errors when finite fallback is unusable", {
  beta <- c(x = 1, y = 2)
  finite_fn <- function(beta_vec, beta_names, mm_cg) matrix(NA_real_, 2, 2)

  expect_error(
    gamlss.longitudinal:::.gl_cg_observed_hessian(
      tmp_obj = list(),
      beta_vec = beta,
      mm_cg = NULL,
      analytical_fn = function(obj) stop("bad analytical path"),
      finite_fn = finite_fn,
      hessian_method = "analytical",
      verbose = 0,
      context = "unit test"
    ),
    "CG failed to construct a usable Hessian during unit test"
  )
})

test_that("CG smooth EDF helper computes bounded trace contribution", {
  beta_names <- c("mu.s(x).1", "mu.s(x).2")
  par_s <- list(mu = list(`s(x)` = c(`mu.s(x).1` = 0, `mu.s(x).2` = 0)))
  H_obs <- matrix(
    c(-2, 0, 0, -2),
    nrow = 2,
    dimnames = list(beta_names, beta_names)
  )
  penalty <- diag(1, 2)
  dimnames(penalty) <- list(beta_names, beta_names)

  out <- gamlss.longitudinal:::.gl_cg_smooth_edf_list(
    H_obs_current = H_obs,
    penalty_current = penalty,
    beta_names = beta_names,
    par_s = par_s
  )

  expect_named(out, "mu")
  expect_named(out$mu, "s(x)")
  expect_equal(out$mu$`s(x)`, 2 * 2 / (2 + 1 + 2e-8), tolerance = 1e-8)
})

test_that("CG smooth EDF helper falls back to block size when solve is unusable", {
  beta_names <- c("mu.s(x).1", "mu.s(x).2")
  par_s <- list(mu = list(`s(x)` = c(`mu.s(x).1` = 0, `mu.s(x).2` = 0)))
  H_obs <- matrix(Inf, nrow = 2, ncol = 2, dimnames = list(beta_names, beta_names))
  penalty <- matrix(0, nrow = 2, ncol = 2, dimnames = list(beta_names, beta_names))

  out <- gamlss.longitudinal:::.gl_cg_smooth_edf_list(
    H_obs_current = H_obs,
    penalty_current = penalty,
    beta_names = beta_names,
    par_s = par_s
  )

  expect_equal(out$mu$`s(x)`, length(beta_names))
})
