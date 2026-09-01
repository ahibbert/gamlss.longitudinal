test_that("convergence metadata preserves RS and CG stopping semantics", {
  rs_info <- gamlss.longitudinal:::.gl_build_convergence_info(
    method = "RS",
    outer_log_lik_change = 1e-8,
    outer_stop_crit = 1e-6,
    outer_only_run_counter = 4L,
    max_outer_iter = 10L,
    cg_stop_reason = NA_character_,
    cg_last_grad_inf = NA_real_,
    cg_last_step_l2 = NA_real_,
    cg_best_raw_loglik = -Inf,
    cg_best_iteration = NA_integer_,
    cg_raw_loglik_drop_from_best = NA_real_,
    cg_raw_loglik_drop_tol = NA_real_,
    cg_gradient_method = "forward",
    cg_zeta_hessian = "analytical",
    cg_hessian_method = "analytical",
    objective = -10
  )

  expect_true(rs_info$converged)
  expect_false(rs_info$hit_outer_limit)
  expect_identical(rs_info$stop_reason, "converged")
  expect_identical(rs_info$method, "RS")
  expect_true(is.na(rs_info$cg_gradient_method))

  cg_info <- gamlss.longitudinal:::.gl_build_convergence_info(
    method = "CG",
    outer_log_lik_change = 0,
    outer_stop_crit = 1e-6,
    outer_only_run_counter = 11L,
    max_outer_iter = 10L,
    cg_stop_reason = "max_stall",
    cg_last_grad_inf = 0.2,
    cg_last_step_l2 = 0.01,
    cg_best_raw_loglik = -12,
    cg_best_iteration = 3L,
    cg_raw_loglik_drop_from_best = 0.5,
    cg_raw_loglik_drop_tol = 10,
    cg_gradient_method = "central",
    cg_zeta_hessian = "finite",
    cg_hessian_method = "auto",
    objective = -10,
    cg_grad_tol = 0.1,
    cg_step_tol = 0.001
  )

  expect_false(cg_info$converged)
  expect_true(cg_info$hit_outer_limit)
  expect_true(cg_info$hit_max_stall)
  expect_identical(cg_info$stop_reason, "max_stall")
  expect_equal(cg_info$grad_inf, 0.2)
  expect_identical(cg_info$cg_gradient_method, "central")
})

test_that("information criteria partition marginal and copula effective degrees of freedom", {
  par_cov <- c(`mu.(Intercept)` = 1, `theta.(Intercept)` = 0.2, `zeta.(Intercept)` = 5)
  par_s <- list(
    mu = list(`s(x)` = c(`mu.s(x).1` = 0, `mu.s(x).2` = 0)),
    theta = list(`s(time)` = c(`theta.s(time).1` = 0)),
    zeta = list()
  )
  df_s <- list(
    mu = list(`s(x)` = 1.5),
    theta = list(`s(time)` = 0.5),
    zeta = list()
  )
  calc_lik <- list(log_lik = c(marginal = -10, copula = -5, joint = -15))
  dataset <- data.frame(response = c(1, 2, 3, 4), time = c(1, 1, 2, 2), subject = c(1, 2, 1, 2))

  aics <- gamlss.longitudinal:::.gl_fit_information_criteria(
    par_cov = par_cov,
    par_s = par_s,
    df_s = df_s,
    calc_lik_out_end = calc_lik,
    dataset = dataset
  )

  expect_identical(rownames(aics), c("LogLik", "AIC", "BIC", "EDF"))
  expect_equal(unname(aics["LogLik", ]), c(-10, -5, -15))
  expect_equal(unname(aics["AIC", ]), c(25, 15, 40))
  expect_equal(unname(aics["EDF", ]), c(2.5, 2.5, 5))
})

test_that("core fit object fields are assembled with stable names", {
  dataset <- data.frame(response = c(1, 2), time = c(1, 2), subject = c(1, 1))
  mm <- list(x = list(mu = matrix(1, nrow = 2, ncol = 1)), s = list(mu = list()))
  calc_lik <- list(log_lik = c(marginal = -2, copula = -1, joint = -3))
  convergence_info <- list(converged = TRUE)

  out <- gamlss.longitudinal:::.gl_build_fit_object_core(
    par_cov = c(`mu.(Intercept)` = 1),
    log_lik_history = matrix(-3),
    par_history = matrix(1),
    calc_lik_out_end = calc_lik,
    mm = mm,
    margin_dist = "NO",
    copula_dist = "N",
    include_dlcopdpar = TRUE,
    dataset = dataset,
    dataset_original = dataset,
    response_var = "response",
    time_var = "time",
    subject_var = "subject",
    formulas = list(mu = response ~ 1),
    formulas_int = list(mu = response ~ 1),
    var_map = list(response = "response"),
    par_s = list(mu = list()),
    lambda_s = list(mu = list()),
    df_s = list(mu = list()),
    weights_final = list(mu = c(1, 1)),
    method = "RS",
    warm_start_info = list(used = FALSE),
    convergence_info = convergence_info
  )

  expect_named(out, c(
    "par", "log_lik_history", "par_history", "calc_lik_out_end",
    "model_matrix", "margin_dist", "copula_dist", "include_dlcopdpar",
    "response", "response_margin", "response_subject", "par_s",
    "lambda_s", "df_s", "weights", "dataset", "dataset_original",
    "response_var", "time_var", "subject_var", "formulas", "formulas_int",
    "var_map", "optim_method", "capability_registry_version",
    "capability_route", "warm_start_joint", "convergence"
  ))
  expect_equal(out$par, c(`mu.(Intercept)` = 1))
  expect_identical(out$response, dataset$response)
  expect_identical(out$response_margin, dataset$time)
  expect_identical(out$response_subject, dataset$subject)
  expect_identical(out$convergence, convergence_info)
})

test_that("fit optimizer traces are attached by optimizer method", {
  base <- list(ok = TRUE)

  rs_out <- gamlss.longitudinal:::.gl_attach_fit_optimizer_traces(
    return_list = base,
    method = "RS",
    rs_block_trace = list(data.frame(parameter = "mu"), data.frame(parameter = "sigma")),
    cg_lambda_trace = data.frame(iteration = 1L),
    cg_step_trace = list(data.frame(iteration = 1L))
  )

  cg_out <- gamlss.longitudinal:::.gl_attach_fit_optimizer_traces(
    return_list = base,
    method = "CG",
    rs_block_trace = list(data.frame(parameter = "mu")),
    cg_lambda_trace = data.frame(iteration = 1L),
    cg_step_trace = list(data.frame(iteration = 1L), data.frame(iteration = 2L))
  )

  empty_rs <- gamlss.longitudinal:::.gl_attach_fit_optimizer_traces(
    return_list = base,
    method = "RS",
    rs_block_trace = list(),
    cg_lambda_trace = data.frame(),
    cg_step_trace = list()
  )

  expect_identical(rs_out$rs_block_trace$parameter, c("mu", "sigma"))
  expect_false("cg_step_trace" %in% names(rs_out))
  expect_identical(cg_out$cg_lambda_trace$iteration, 1L)
  expect_identical(cg_out$cg_step_trace$iteration, 1:2)
  expect_false("rs_block_trace" %in% names(cg_out))
  expect_true(is.data.frame(empty_rs$rs_block_trace))
  expect_identical(nrow(empty_rs$rs_block_trace), 0L)
})

test_that("final fit object assembly preserves fitted object fields and traces", {
  dataset <- data.frame(response = c(1, 2), time = c(1, 2), subject = c(1, 1), x = c(0, 1))
  mm <- list(
    x = list(mu = matrix(1, nrow = 2, ncol = 1, dimnames = list(NULL, "mu.(Intercept)"))),
    s = list(mu = list()),
    formulae = list(mu = response ~ 1)
  )
  calc_lik <- list(log_lik = c(marginal = -2, copula = -1, joint = -3))
  convergence_info <- gamlss.longitudinal:::.gl_build_convergence_info(
    method = "RS",
    outer_log_lik_change = 0,
    outer_stop_crit = 1e-6,
    outer_only_run_counter = 2L,
    max_outer_iter = 10L,
    cg_stop_reason = NA_character_,
    cg_last_grad_inf = NA_real_,
    cg_last_step_l2 = NA_real_,
    cg_best_raw_loglik = -Inf,
    cg_best_iteration = NA_integer_,
    cg_raw_loglik_drop_from_best = NA_real_,
    cg_raw_loglik_drop_tol = NA_real_,
    cg_gradient_method = "forward",
    cg_zeta_hessian = "analytical",
    cg_hessian_method = "analytical"
  )

  fit <- gamlss.longitudinal:::.gl_finalize_fit_object(
    par_cov = c(`mu.(Intercept)` = 1),
    log_lik_history = matrix(c(-2, -1, -3), nrow = 1, dimnames = list(NULL, c("marginal", "copula", "joint"))),
    par_history = matrix(1, nrow = 1, dimnames = list(NULL, "mu.(Intercept)")),
    calc_lik_out_end = calc_lik,
    mm = mm,
    margin_dist = list(family = c("NO", "Normal")),
    copula_dist = "N",
    include_dlcopdpar = TRUE,
    dataset = dataset,
    dataset_original = dataset,
    response_var = "response",
    time_var = "time",
    subject_var = "subject",
    formulas = list(mu = response ~ 1),
    formulas_int = list(mu = response ~ 1),
    var_map = list(response = "response"),
    par_s = list(mu = list()),
    lambda_s = list(mu = list()),
    df_s = list(mu = list()),
    weights_final = list(mu = c(1, 1)),
    method = "RS",
    warm_start_info = list(used = FALSE),
    rs_block_trace = list(data.frame(parameter = "mu", accepted = TRUE)),
    cg_lambda_trace = data.frame(),
    cg_step_trace = list(),
    convergence_info = convergence_info,
    fit_start_time = Sys.time(),
    compute_vcov = FALSE,
    vcov_numderiv = FALSE,
    vcov_method = "analytical",
    verbose = 0
  )

  expect_s3_class(fit, "gamlss.longitudinal")
  expect_equal(fit$par, c(`mu.(Intercept)` = 1))
  expect_identical(fit$calc_lik_out_end, calc_lik)
  expect_identical(fit$model_matrix, mm)
  expect_identical(fit$response_var, "response")
  expect_true(is.data.frame(fit$rs_block_trace))
  expect_equal(fit$rs_block_trace$parameter, "mu")
  expect_false(fit$vcov_meta$precomputed)
  expect_identical(fit$vcov_meta$method, "analytical")
})

test_that("final fit workflow builds convergence metadata before final object assembly", {
  captured <- list()
  optimizer_state <- list(
    par_cov = c(`mu.(Intercept)` = 1),
    log_lik_history = matrix(c(-2, -1, -3), nrow = 1, dimnames = list(NULL, c("marginal", "copula", "joint"))),
    par_history = matrix(1, nrow = 1, dimnames = list(NULL, "mu.(Intercept)")),
    calc_lik_out_end = list(log_lik = c(marginal = -2, copula = -1, joint = -3)),
    par_s = list(mu = list()),
    lambda_s = list(mu = list()),
    df_s = list(mu = list()),
    weights_final = list(mu = c(1, 1)),
    rs_block_trace = list(data.frame(parameter = "mu")),
    cg_lambda_trace = data.frame(),
    cg_step_trace = list(),
    outer_log_lik_change = 0,
    outer_stop_crit = 1e-6,
    outer_only_run_counter = 2L,
    cg_stop_reason = NA_character_,
    cg_last_grad_inf = NA_real_,
    cg_last_step_l2 = NA_real_,
    cg_best_raw_loglik = -Inf,
    cg_best_iteration = NA_integer_,
    cg_raw_loglik_drop_from_best = NA_real_
  )
  fit_data <- list(
    dataset = data.frame(response = c(1, 2), time = c(1, 2), subject = c(1, 1)),
    dataset_original = data.frame(y = c(1, 2), t = c(1, 2), id = c(1, 1)),
    response_var = "y",
    formulas_int = list(mu = response ~ 1),
    var_map = list(y = "response")
  )
  mm <- list(x = list(mu = matrix(1, nrow = 2, ncol = 1)), s = list(mu = list()))
  matrix_bundle <- list(mm = mm, copula_link = "identity")
  convergence_info <- list(converged = TRUE, hit_outer_limit = FALSE)
  fit_out <- structure(list(ok = TRUE), class = "gamlss.longitudinal")

  out <- gamlss.longitudinal:::.gl_finalize_fit_workflow(
    optimizer_state = optimizer_state,
    fit_data = fit_data,
    matrix_bundle = matrix_bundle,
    margin_dist = list(family = c("NO", "Normal")),
    copula_dist = "N",
    include_dlcopdpar = TRUE,
    method = "RS",
    original_formulas = list(mu = y ~ 1),
    time_var = "t",
    subject_var = "id",
    warm_start_info = list(used = FALSE),
    fit_start_time = Sys.time(),
    compute_vcov = FALSE,
    vcov_numderiv = FALSE,
    vcov_method = "analytical",
    verbose = 0,
    max_outer_iter = 10L,
    cg_raw_loglik_drop_tol = 10,
    cg_gradient_method = "forward",
    cg_zeta_hessian = "analytical",
    cg_hessian_method = "analytical",
    convergence_fn = function(...) {
      captured$convergence <<- list(...)
      convergence_info
    },
    finalize_fn = function(...) {
      captured$finalize <<- list(...)
      fit_out
    }
  )

  expect_s3_class(out, "gamlss.longitudinal")
  expect_equal(captured$convergence$outer_log_lik_change, optimizer_state$outer_log_lik_change)
  expect_equal(captured$convergence$outer_stop_crit, optimizer_state$outer_stop_crit)
  expect_equal(captured$convergence$cg_raw_loglik_drop_tol, 10)
  expect_equal(captured$finalize$par_cov, optimizer_state$par_cov)
  expect_identical(captured$finalize$mm, mm)
  expect_identical(captured$finalize$dataset, fit_data$dataset)
  expect_identical(captured$finalize$dataset_original, fit_data$dataset_original)
  expect_identical(captured$finalize$formulas, list(mu = y ~ 1))
  expect_identical(captured$finalize$formulas_int, fit_data$formulas_int)
  expect_identical(captured$finalize$convergence_info, convergence_info)
})

test_that("prepared fit finalization bridge forwards normalized workflow state", {
  captured <- list()
  optimizer_state <- list(par_cov = c(`mu.(Intercept)` = 1))
  fit_data <- list(dataset = data.frame(response = c(1, 2)))
  matrix_bundle <- list(mm = list(x = list(), s = list()))
  warm_start_info <- list(used = TRUE, outer_iter = 2L)
  fit_start_time <- Sys.time()
  original_formulas <- list(mu = y ~ 1, sigma = ~ 1)

  fit_workflow <- list(
    controls = list(
      method = "CG",
      vcov_numderiv = TRUE,
      vcov_method = "numderiv",
      cg_raw_loglik_drop_tol = 7,
      cg_gradient_method = "central",
      cg_zeta_hessian = "finite",
      cg_hessian_method = "auto"
    ),
    fit_data = fit_data,
    matrix_bundle = matrix_bundle,
    warm_start = list(warm_start_info = warm_start_info)
  )

  out <- gamlss.longitudinal:::.gl_finalize_prepared_fit(
    optimizer_state = optimizer_state,
    fit_workflow = fit_workflow,
    margin_dist = "NO",
    copula_dist = "N",
    include_dlcopdpar = TRUE,
    original_formulas = original_formulas,
    time_var = "time",
    subject_var = "subject",
    fit_start_time = fit_start_time,
    compute_vcov = FALSE,
    verbose = 2,
    max_outer_iter = 5L,
    finalize_fn = function(...) {
      captured$args <<- list(...)
      list(ok = TRUE)
    }
  )

  expect_identical(out, list(ok = TRUE))
  expect_identical(captured$args$optimizer_state, optimizer_state)
  expect_identical(captured$args$fit_data, fit_data)
  expect_identical(captured$args$matrix_bundle, matrix_bundle)
  expect_equal(captured$args$method, "CG")
  expect_identical(captured$args$original_formulas, original_formulas)
  expect_identical(captured$args$warm_start_info, warm_start_info)
  expect_true(captured$args$vcov_numderiv)
  expect_equal(captured$args$vcov_method, "numderiv")
  expect_equal(captured$args$cg_raw_loglik_drop_tol, 7)
  expect_equal(captured$args$cg_gradient_method, "central")
  expect_equal(captured$args$cg_zeta_hessian, "finite")
  expect_equal(captured$args$cg_hessian_method, "auto")
})
