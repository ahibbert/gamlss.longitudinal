make_control_args <- function(...) {
  defaults <- list(
    method = "rs",
    start_from = NA,
    warm_start_joint = TRUE,
    warm_start_joint_iter = 5,
    backtracking_max_halves = 50,
    cg_max_stall = 5,
    cg_max_delta = 0.5,
    cg_lambda_update_every = 10,
    cg_max_lambda_updates = NA,
    cg_raw_loglik_drop_tol = 10,
    cg_line_search = "best",
    cg_max_line_search_evals = NA,
    cg_gradient_method = "forward",
    discrete_score_method = c("analytical", "finite"),
    cg_zeta_hessian = "analytical",
    cg_hessian_method = c("analytical", "finite", "auto"),
    vcov_method = c("analytical", "numderiv"),
    vcov_numderiv = FALSE,
    rs_smooth_trust_radius = Inf
  )
  utils::modifyList(defaults, list(...))
}

normalize_fit_controls <- function(...) {
  do.call(gamlss.longitudinal:::.gl_normalize_fit_controls, make_control_args(...))
}

test_that("fit control normalization preserves existing defaults", {
  controls <- normalize_fit_controls()

  expect_equal(controls$method, "RS")
  expect_false(controls$user_supplied_start)
  expect_equal(controls$warm_start_joint_iter, 5L)
  expect_equal(controls$backtracking_max_halves, 50L)
  expect_equal(controls$cg_max_line_search_evals, Inf)
  expect_equal(controls$cg_max_lambda_updates, Inf)
  expect_equal(controls$cg_gradient_method, "forward")
  expect_equal(controls$discrete_score_method, "analytical")
  expect_equal(controls$cg_hessian_method, "analytical")
  expect_equal(controls$vcov_method, "analytical")
  expect_false(controls$vcov_numderiv)
})

test_that("fit control normalization preserves compatibility switches", {
  controls <- normalize_fit_controls(
    method = "cg",
    start_from = c(mu = 1),
    cg_line_search = "first",
    cg_max_line_search_evals = 3,
    cg_max_lambda_updates = 2,
    vcov_numderiv = TRUE
  )

  expect_equal(controls$method, "CG")
  expect_true(controls$user_supplied_start)
  expect_equal(controls$cg_line_search, "first")
  expect_equal(controls$cg_max_line_search_evals, 3L)
  expect_equal(controls$cg_max_lambda_updates, 2L)
  expect_equal(controls$vcov_method, "numderiv")
  expect_true(controls$vcov_numderiv)
})

test_that("fit control normalization preserves existing validation messages", {
  expect_error(
    normalize_fit_controls(backtracking_max_halves = -1),
    "backtracking_max_halves must be a single non-negative integer",
    fixed = TRUE
  )
  expect_error(
    normalize_fit_controls(method = "bad"),
    "method must be one of 'RS' or 'CG'",
    fixed = TRUE
  )
  expect_error(
    normalize_fit_controls(warm_start_joint = NA),
    "warm_start_joint must be TRUE or FALSE",
    fixed = TRUE
  )
})

test_that("fit-control scalar normalizers preserve validation contracts", {
  expect_error(
    gamlss.longitudinal:::.gl_normalize_backtracking_halves(2.9),
    "backtracking_max_halves must be a single non-negative integer",
    fixed = TRUE
  )
  expect_error(
    gamlss.longitudinal:::.gl_normalize_backtracking_halves(NA_real_),
    "backtracking_max_halves must be a single non-negative integer",
    fixed = TRUE
  )

  expect_equal(gamlss.longitudinal:::.gl_validate_rs_smooth_trust_radius(Inf), Inf)
  expect_error(
    gamlss.longitudinal:::.gl_validate_rs_smooth_trust_radius(0),
    "rs_smooth_trust_radius must be a single positive numeric value or Inf",
    fixed = TRUE
  )

  warm <- gamlss.longitudinal:::.gl_normalize_warm_start_controls(TRUE, 3.8)
  expect_true(warm$warm_start_joint)
  expect_equal(warm$warm_start_joint_iter, 3L)

  vcov <- gamlss.longitudinal:::.gl_normalize_vcov_controls("analytical", TRUE)
  expect_equal(vcov$vcov_method, "numderiv")
  expect_true(vcov$vcov_numderiv)
})

test_that("CG fit-control normalizers preserve NA-as-Inf and fallback policies", {
  expect_equal(gamlss.longitudinal:::.gl_normalize_cg_line_search_evals(NA), Inf)
  expect_error(
    gamlss.longitudinal:::.gl_normalize_cg_line_search_evals(4.9),
    "cg_max_line_search_evals must be a single non-negative integer or NA",
    fixed = TRUE
  )
  expect_error(
    gamlss.longitudinal:::.gl_normalize_cg_line_search_evals(-1),
    "cg_max_line_search_evals must be a single non-negative integer or NA",
    fixed = TRUE
  )

  expect_error(
    gamlss.longitudinal:::.gl_normalize_cg_lambda_controls(
      cg_lambda_update_every = 2.8,
      cg_max_lambda_updates = NA,
      cg_raw_loglik_drop_tol = NA
    ),
    "cg_lambda_update_every must be a positive integer",
    fixed = TRUE
  )
  lambda <- gamlss.longitudinal:::.gl_normalize_cg_lambda_controls(
    cg_lambda_update_every = 2L,
    cg_max_lambda_updates = NA,
    cg_raw_loglik_drop_tol = NA
  )
  expect_equal(lambda$cg_lambda_update_every, 2L)
  expect_equal(lambda$cg_max_lambda_updates, Inf)
  expect_equal(lambda$cg_raw_loglik_drop_tol, Inf)

  fallback <- gamlss.longitudinal:::.gl_normalize_cg_fallback_controls(
    cg_max_stall = 0,
    cg_max_delta = -1
  )
  expect_equal(fallback$cg_max_stall, 5L)
  expect_equal(fallback$cg_max_delta, 0.5)
})

test_that("step control normalization preserves defaults and automatic messages", {
  joint_msg <- capture.output(
    joint_controls <- gamlss.longitudinal:::.gl_normalize_step_controls(
      method = "RS",
      include_dlcopdpar = TRUE,
      start_step_size = 0.5,
      max_steps = 5,
      step_adjustment = NA,
      verbose = 1
    )
  )
  expect_equal(joint_controls$start_step_size, 0.5)
  expect_equal(joint_controls$max_steps, 5L)
  expect_equal(joint_controls$step_adjustment, 1)
  expect_match(paste(joint_msg, collapse = "\n"), "automatic step_adjustment=1 for joint RS", fixed = TRUE)

  separate_msg <- capture.output(
    separate_controls <- gamlss.longitudinal:::.gl_normalize_step_controls(
      method = "RS",
      include_dlcopdpar = FALSE,
      start_step_size = 0.25,
      max_steps = 2,
      step_adjustment = NA,
      verbose = 1
    )
  )
  expect_equal(separate_controls$step_adjustment, 1)
  expect_match(paste(separate_msg, collapse = "\n"), "automatic step_adjustment=1 for separate RS", fixed = TRUE)

  cg_msg <- capture.output(
    cg_controls <- gamlss.longitudinal:::.gl_normalize_step_controls(
      method = "CG",
      include_dlcopdpar = TRUE,
      start_step_size = 0.5,
      max_steps = 5,
      step_adjustment = NA,
      verbose = 1
    )
  )
  expect_equal(cg_controls$step_adjustment, 1)
  expect_match(paste(cg_msg, collapse = "\n"), "automatic step_adjustment=1 for CG", fixed = TRUE)
})

test_that("step control normalization preserves explicit values and validation messages", {
  controls <- gamlss.longitudinal:::.gl_normalize_step_controls(
    method = "RS",
    include_dlcopdpar = TRUE,
    start_step_size = 0.2,
    max_steps = 3.9,
    step_adjustment = 0.8,
    verbose = 0
  )
  expect_equal(controls$start_step_size, 0.2)
  expect_equal(controls$max_steps, 3L)
  expect_equal(controls$step_adjustment, 0.8)

  expect_error(
    gamlss.longitudinal:::.gl_normalize_step_controls(
      method = "RS",
      include_dlcopdpar = TRUE,
      start_step_size = 0,
      max_steps = 5,
      step_adjustment = NA,
      verbose = 0
    ),
    "start_step_size must be a single positive finite numeric value",
    fixed = TRUE
  )
  expect_error(
    gamlss.longitudinal:::.gl_normalize_step_controls(
      method = "RS",
      include_dlcopdpar = TRUE,
      start_step_size = 0.5,
      max_steps = -1,
      step_adjustment = NA,
      verbose = 0
    ),
    "max_steps must be a single non-negative integer",
    fixed = TRUE
  )
  expect_error(
    gamlss.longitudinal:::.gl_normalize_step_controls(
      method = "RS",
      include_dlcopdpar = TRUE,
      start_step_size = 0.5,
      max_steps = 5,
      step_adjustment = 0,
      verbose = 0
    ),
    "step_adjustment must be a single positive numeric value, or NA for the method-specific default",
    fixed = TRUE
  )
})

test_that("elapsed budget helper preserves timing policy and error message", {
  start <- as.POSIXct("2026-01-01 00:00:00", tz = "UTC")

  expect_invisible(gamlss.longitudinal:::.gl_check_elapsed_budget(
    fit_start_time = start,
    max_elapsed_sec = Inf,
    stage = "setup",
    now = start + 100
  ))
  expect_invisible(gamlss.longitudinal:::.gl_check_elapsed_budget(
    fit_start_time = start,
    max_elapsed_sec = 0,
    stage = "setup",
    now = start + 100
  ))
  expect_invisible(gamlss.longitudinal:::.gl_check_elapsed_budget(
    fit_start_time = start,
    max_elapsed_sec = 10,
    stage = "setup",
    now = start + 9
  ))

  expect_error(
    gamlss.longitudinal:::.gl_check_elapsed_budget(
      fit_start_time = start,
      max_elapsed_sec = 10,
      stage = "CG outer iteration",
      now = start + 11
    ),
    "Model exceeded max_elapsed_sec during CG outer iteration \\(elapsed 11.0 sec > 10.0 sec\\)."
  )
})

test_that("elapsed budget checker factory binds fit timing state", {
  start <- as.POSIXct("2026-01-01 00:00:00", tz = "UTC")
  captured <- list()

  checker <- gamlss.longitudinal:::.gl_build_elapsed_budget_checker(
    fit_start_time = start,
    max_elapsed_sec = 12,
    check_fn = function(...) {
      captured$args <<- list(...)
      invisible(TRUE)
    }
  )

  expect_invisible(checker())
  expect_identical(captured$args$fit_start_time, start)
  expect_equal(captured$args$max_elapsed_sec, 12)
  expect_equal(captured$args$stage, "optimisation")

  expect_invisible(checker("RS inner loop"))
  expect_equal(captured$args$stage, "RS inner loop")
})
