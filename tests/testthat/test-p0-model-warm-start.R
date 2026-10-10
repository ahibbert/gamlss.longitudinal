make_warm_start_args <- function(...) {
  defaults <- list(
    start_from = NA,
    method = "RS",
    include_dlcopdpar = TRUE,
    warm_start_joint = TRUE,
    warm_start_joint_iter = 2L,
    user_supplied_start = FALSE,
    dataset_original = data.frame(response = 1:2, time = 1:2, subject = 1),
    margin_dist = list(family = c("NO", "Normal")),
    copula_dist = "N",
    time_var = "time",
    subject_var = "subject",
    missingness = "error",
    mu.formula = response ~ 1,
    sigma.formula = ~ 1,
    nu.formula = ~ 1,
    tau.formula = ~ 1,
    theta.formula = ~ 1,
    zeta.formula = ~ 1,
    optimizer_control = gamlss_longitudinal_control(
      outer_tol = 1e-6,
      rs = list(
        inner_tol = 1e-6,
        start_step_size = 0.5,
        step_adjustment = 1,
        max_steps = 5L,
        max_inner_iter = 3L,
        initial_inner_iter = 2L,
        max_negative_outer_streak = 10L,
        adaptive_step_size = TRUE,
        lambda_update_every = 7L,
        bfgs_refine = TRUE,
        bfgs_refine_after = 4L,
        bfgs_max_iter = 25L,
        bfgs_reltol = 1e-10,
        bfgs_coordinate_rescue = TRUE,
        bfgs_coordinate_threshold = 2,
        bfgs_coordinate_radius = 5,
        bfgs_coordinate_max = 4L,
        bfgs_coordinate_tol = 1e-4,
        bfgs_repolish_iter = 12L
      )
    ),
    true_val = NA,
    vcov_method = "analytical",
    vcov_numderiv = FALSE,
    use_Rcpp = FALSE,
    lambda_start = NA,
    lambda_penalty_K = 2,
    verbose = 0,
    fit_fn = function(...) {
      list(
        par = c(`mu.(Intercept)` = 1),
        par_s = list(mu = list(`s(x)` = c(`mu.s(x).1` = 0.1))),
        calc_lik_out_end = list(log_lik = c(marginal = -2, copula = -1, joint = -3))
      )
    }
  )
  utils::modifyList(defaults, list(...))
}

run_warm_start <- function(...) {
  do.call(gamlss.longitudinal:::.gl_run_joint_warm_start, make_warm_start_args(...))
}

run_warm_start_fit_call <- function(...) {
  args <- make_warm_start_args(...)
  formals_needed <- names(formals(gamlss.longitudinal:::.gl_call_joint_warm_start_fit))
  args <- args[intersect(names(args), formals_needed)]
  do.call(gamlss.longitudinal:::.gl_call_joint_warm_start_fit, args)
}

test_that("warm-start failures teach the structured optimizer control", {
  err <- tryCatch(
    gamlss.longitudinal:::.gl_validate_joint_warm_start_fit(
      NULL,
      simpleError("example failure")
    ),
    error = identity
  )
  expect_match(conditionMessage(err), "optimizer_control =", fixed = TRUE)
  expect_match(conditionMessage(err), "rs = list(warm_start_joint = FALSE)", fixed = TRUE)
  expect_false(grepl("Set warm_start_joint = FALSE", conditionMessage(err), fixed = TRUE))
})

test_that("warm-start policy preserves the previous joint RS conditions", {
  expect_true(gamlss.longitudinal:::.gl_should_run_joint_warm_start(
    method = "RS",
    include_dlcopdpar = TRUE,
    warm_start_joint = TRUE,
    warm_start_joint_iter = 1L,
    user_supplied_start = FALSE
  ))
  expect_false(gamlss.longitudinal:::.gl_should_run_joint_warm_start("CG", TRUE, TRUE, 1L, FALSE))
  expect_false(gamlss.longitudinal:::.gl_should_run_joint_warm_start("RS", FALSE, TRUE, 1L, FALSE))
  expect_false(gamlss.longitudinal:::.gl_should_run_joint_warm_start("RS", TRUE, FALSE, 1L, FALSE))
  expect_false(gamlss.longitudinal:::.gl_should_run_joint_warm_start("RS", TRUE, TRUE, 0L, FALSE))
  expect_false(gamlss.longitudinal:::.gl_should_run_joint_warm_start("RS", TRUE, TRUE, 1L, TRUE))
})

test_that("warm-start helper returns default metadata when skipped", {
  out <- run_warm_start(
    user_supplied_start = TRUE,
    start_from = c(`mu.(Intercept)` = 2),
    fit_fn = function(...) stop("should not run")
  )

  expect_equal(out$start_from, c(`mu.(Intercept)` = 2))
  expect_null(out$warm_start_par_s)
  expect_false(out$warm_start_info$used)
  expect_equal(out$warm_start_info$outer_iter, 0L)
  expect_null(out$warm_start_info$log_lik)
})

test_that("warm-start helper captures fit output warnings and starts", {
  fit_args <- NULL
  warm_fit <- function(...) {
    fit_args <<- list(...)
    cat("warm output line\n")
    warning("warm warning")
    list(
      par = c(`mu.(Intercept)` = 1.2),
      par_s = list(mu = list(`s(x)` = c(`mu.s(x).1` = 0.4))),
      calc_lik_out_end = list(log_lik = c(marginal = -4, copula = -2, joint = -6))
    )
  }

  top_output <- capture.output(
    out <- run_warm_start(verbose = 2, fit_fn = warm_fit)
  )

  expect_equal(out$start_from, c(`mu.(Intercept)` = 1.2))
  expect_equal(out$warm_start_par_s$mu$`s(x)`, c(`mu.s(x).1` = 0.4))
  expect_true(out$warm_start_info$used)
  expect_equal(out$warm_start_info$outer_iter, 2L)
  expect_equal(out$warm_start_info$log_lik, c(marginal = -4, copula = -2, joint = -6))
  expect_true(out$warm_start_info$carries_smooth)
  expect_false(out$warm_start_info$carries_lambda)
  expect_match(out$warm_start_info$captured_output, "warm output line")
  expect_equal(out$warm_start_info$captured_warnings, "warm warning")
  expect_match(paste(top_output, collapse = "\n"), "Running separate RS warm-start phase", fixed = TRUE)
  expect_match(paste(top_output, collapse = "\n"), "warm output line", fixed = TRUE)

  expect_identical(fit_args$optimization, "separate")
  expect_false("check_dlcopdpar_gradient" %in% names(fit_args))
  expect_false(fit_args$optimizer_control$rs$warm_start_joint)
  expect_equal(fit_args$optimizer_control$rs$warm_start_joint_iter, 0L)
  expect_equal(fit_args$optimizer_control$shared$max_outer_iter, 2L)
  expect_identical(fit_args$optimizer_control$rs$lambda_update_every, 7L)
  expect_identical(fit_args$optimizer_control$rs$initial_inner_iter, 2L)
  expect_true(fit_args$optimizer_control$rs$adaptive_step_size)
  expect_false(fit_args$optimizer_control$rs$bfgs_refine)
  expect_identical(fit_args$optimizer_control$rs$bfgs_refine_after, 4L)
  expect_identical(fit_args$optimizer_control$rs$bfgs_max_iter, 25L)
  expect_equal(fit_args$optimizer_control$rs$bfgs_reltol, 1e-10)
  expect_true(fit_args$optimizer_control$rs$bfgs_coordinate_rescue)
  expect_equal(fit_args$optimizer_control$rs$bfgs_coordinate_threshold, 2)
  expect_equal(fit_args$optimizer_control$rs$bfgs_coordinate_radius, 5)
  expect_identical(fit_args$optimizer_control$rs$bfgs_coordinate_max, 4L)
  expect_equal(fit_args$optimizer_control$rs$bfgs_coordinate_tol, 1e-4)
  expect_identical(fit_args$optimizer_control$rs$bfgs_repolish_iter, 12L)
  expect_false(fit_args$compute_vcov)
})

test_that("warm-start result helper records accepted fit metadata", {
  warm_fit <- list(
    par = c(`mu.(Intercept)` = 1.2),
    par_s = list(mu = list(`s(x)` = c(`mu.s(x).1` = 0.4))),
    lambda_s = list(mu = list(`s(x)` = 3)),
    calc_lik_out_end = list(log_lik = c(marginal = -4, copula = -2, joint = -6))
  )

  expect_invisible(gamlss.longitudinal:::.gl_validate_joint_warm_start_fit(warm_fit, NULL))
  out <- gamlss.longitudinal:::.gl_joint_warm_start_result(
    warm_fit = warm_fit,
    warm_output = "warm output",
    warm_warnings = "warm warning",
    warm_start_joint_iter = 2L
  )

  expect_equal(out$start_from, c(`mu.(Intercept)` = 1.2))
  expect_equal(out$warm_start_par_s, warm_fit$par_s)
  expect_equal(out$warm_start_lambda_s, warm_fit$lambda_s)
  expect_true(out$warm_start_info$carries_lambda)
  expect_true(out$warm_start_info$used)
  expect_equal(out$warm_start_info$outer_iter, 2L)
  expect_false(out$warm_start_info$include_dlcopdpar)
  expect_equal(out$warm_start_info$log_lik, c(marginal = -4, copula = -2, joint = -6))
  expect_true(out$warm_start_info$carries_smooth)
  expect_equal(out$warm_start_info$captured_output, "warm output")
  expect_equal(out$warm_start_info$captured_warnings, "warm warning")
})

test_that("warm-start fit-call helper captures output warnings and cold-start arguments", {
  fit_args <- NULL

  out <- run_warm_start_fit_call(
    fit_fn = function(...) {
      fit_args <<- list(...)
      cat("fit call output\n")
      warning("fit call warning")
      list(
        par = c(`mu.(Intercept)` = 2),
        par_s = list(mu = list()),
        calc_lik_out_end = list(log_lik = c(joint = -5))
      )
    }
  )

  expect_equal(out$fit$par, c(`mu.(Intercept)` = 2))
  expect_match(out$output, "fit call output")
  expect_equal(out$warnings, "fit call warning")
  expect_null(out$error)
  expect_identical(fit_args$optimization, "separate")
  expect_false("check_dlcopdpar_gradient" %in% names(fit_args))
  expect_false(fit_args$optimizer_control$rs$warm_start_joint)
  expect_equal(fit_args$optimizer_control$rs$warm_start_joint_iter, 0L)
  expect_equal(fit_args$optimizer_control$shared$max_outer_iter, 2L)
  expect_identical(fit_args$optimizer_control$rs$lambda_update_every, 7L)
  expect_identical(fit_args$optimizer_control$rs$initial_inner_iter, 2L)
  expect_true(fit_args$optimizer_control$rs$adaptive_step_size)
  expect_false(fit_args$optimizer_control$rs$bfgs_refine)
  expect_identical(fit_args$optimizer_control$rs$bfgs_refine_after, 4L)
  expect_identical(fit_args$optimizer_control$rs$bfgs_max_iter, 25L)
  expect_equal(fit_args$optimizer_control$rs$bfgs_reltol, 1e-10)
  expect_true(fit_args$optimizer_control$rs$bfgs_coordinate_rescue)
  expect_equal(fit_args$optimizer_control$rs$bfgs_coordinate_threshold, 2)
  expect_equal(fit_args$optimizer_control$rs$bfgs_coordinate_radius, 5)
  expect_identical(fit_args$optimizer_control$rs$bfgs_coordinate_max, 4L)
  expect_equal(fit_args$optimizer_control$rs$bfgs_coordinate_tol, 1e-4)
  expect_identical(fit_args$optimizer_control$rs$bfgs_repolish_iter, 12L)
  expect_false(fit_args$compute_vcov)
})

test_that("warm-start runner delegates separate fit call", {
  captured <- NULL

  out <- run_warm_start(
    warm_start_fit_call_fn = function(...) {
      captured <<- list(...)
      list(
        fit = list(
          par = c(`mu.(Intercept)` = 3),
          par_s = list(mu = list()),
          calc_lik_out_end = list(log_lik = c(marginal = -1, copula = -2, joint = -3))
        ),
        output = "delegated output",
        warnings = "delegated warning",
        error = NULL
      )
    }
  )

  expect_equal(captured$warm_start_joint_iter, 2L)
  expect_equal(captured$method, "RS")
  expect_equal(out$start_from, c(`mu.(Intercept)` = 3))
  expect_equal(out$warm_start_info$captured_output, "delegated output")
  expect_equal(out$warm_start_info$captured_warnings, "delegated warning")
})

test_that("warm-start helper preserves failure messages", {
  expect_error(
    run_warm_start(fit_fn = function(...) stop("boom")),
    "Separate RS warm-start phase failed: boom",
    fixed = TRUE
  )
  expect_error(
    run_warm_start(fit_fn = function(...) list(foo = list(), calc_lik_out_end = list(log_lik = c(joint = -1)))),
    "Separate RS warm-start phase did not return coefficient starting values",
    fixed = TRUE
  )
})

test_that("a fitted model start carries fixed smooth and smoothing-penalty state", {
  prior <- structure(list(
    par = c(
      `mu.(Intercept)` = 1,
      `mu.x` = 2,
      `sigma.(Intercept)` = -0.5
    ),
    par_s = list(
      mu = list(`s(z)` = c(`mu.s(z).1` = 0.2, `mu.s(z).2` = -0.1)),
      sigma = list()
    ),
    lambda_s = list(mu = list(`s(z)` = 7), sigma = list())
  ), class = "gamlss.longitudinal")

  warm <- run_warm_start(
    start_from = prior,
    user_supplied_start = TRUE,
    fit_fn = function(...) stop("a supplied fitted state must skip generic warm start")
  )
  expect_equal(warm$start_from, prior$par)
  expect_equal(warm$warm_start_par_s, prior$par_s)
  expect_equal(warm$warm_start_lambda_s, prior$lambda_s)
  expect_false(warm$warm_start_info$used)

  mm <- list(
    x = list(
      mu = structure(matrix(1, 3, 1), dimnames = list(NULL, "(Intercept)")),
      sigma = structure(matrix(1, 3, 1), dimnames = list(NULL, "(Intercept)"))
    ),
    s = list(
      mu = list(`s(z)` = structure(matrix(1, 3, 2), penalty = diag(2))),
      sigma = list()
    )
  )
  initial <- gamlss.longitudinal:::.gl_build_initial_parameter_state(
    start_from = warm$start_from,
    warm_start_par_s = warm$warm_start_par_s,
    warm_start_lambda_s = warm$warm_start_lambda_s,
    mm = mm,
    margin_dist = list(),
    copula_dist = "G",
    dataset = data.frame(response = 1:3),
    lambda_start = NA
  )

  expect_equal(
    initial$par_cov,
    c(`mu.(Intercept)` = 1, `sigma.(Intercept)` = -0.5)
  )
  expect_equal(initial$par_s$mu$`s(z)`, prior$par_s$mu$`s(z)`)
  expect_equal(unname(initial$lambda_s$mu$`s(z)`), 7)
})

test_that("reduced-formula warm starts keep mapped coefficients with projection fallback", {
  prior <- structure(list(
    par = c(`sigma.intercept` = 10, `sigma.week2` = -10),
    par_s = list(sigma = list()),
    lambda_s = list(sigma = list()),
    model_matrix = list(
      x = list(sigma = data.frame(
        intercept = rep(1, 4),
        week2 = c(0, 1, 0, 1)
      )),
      s = list(sigma = list())
    )
  ), class = "gamlss.longitudinal")
  target_mm <- list(
    x = list(sigma = structure(
      matrix(1, nrow = 4, ncol = 1),
      dimnames = list(NULL, "intercept")
    )),
    s = list(sigma = list())
  )

  projected <- gamlss.longitudinal:::.gl_project_fitted_start(
    prior, target_mm
  )

  expect_equal(unname(projected$par[["sigma.intercept"]]), 5)
  expect_identical(projected$projected, "sigma")
  expect_equal(unname(projected$rmse[["sigma"]]), 5)

  warm <- run_warm_start(
    start_from = prior,
    user_supplied_start = TRUE,
    target_mm = target_mm,
    fit_fn = function(...) stop("a supplied fitted state must skip generic warm start")
  )
  expect_equal(unname(warm$start_from[["sigma.intercept"]]), 10)
  expect_null(warm$mapped_start_fallback)
  expect_equal(
    unname(warm$projected_start_fallback$par[["sigma.intercept"]]),
    5
  )
  expect_identical(warm$warm_start_info$projected_parameters, "sigma")
  expect_false(warm$warm_start_info$projection_used)
})
