named_hessian <- function(x, names = c("mu.x", "sigma.x")) {
  dimnames(x) <- list(names, names)
  x
}

expect_inference_failure <- function(expr, code) {
  err <- tryCatch(expr, error = identity)
  expect_s3_class(err, "gamlss_longitudinal_inference_unavailable")
  expect_true(code %in% err$diagnostics$failure_codes)
  invisible(err)
}

test_that("inference control exposes stable profiles and records overrides", {
  standard <- inference_control("standard")
  strict <- inference_control("strict")
  overridden <- inference_control("standard", condition_max = 123)

  expect_s3_class(standard, "gamlss_longitudinal_inference_control")
  expect_false(standard$check_agreement)
  expect_true(strict$check_agreement)
  expect_lt(strict$symmetry_tol, standard$symmetry_tol)
  expect_equal(overridden$condition_max, 123)
  expect_equal(overridden$expert_override, list(condition_max = 123))
  expect_error(inference_control("standard", rank_tol = 0), "positive finite")
})

test_that("valid curvature returns signed covariance and complete diagnostics", {
  H <- named_hessian(matrix(c(-4, -0.2, -0.2, -2), 2))
  gradient <- list(
    gradient = c("mu.x" = 1e-5, "sigma.x" = -1e-5),
    steps = c(1e-4, 1e-4),
    scaled_max = 1e-5
  )
  out <- gamlss.longitudinal:::.gl_validate_hessian_inference(
    H, parameter_names = colnames(H), gradient = gradient,
    source = "analytical", fallback = list(from = "none")
  )

  expect_true(all(diag(out$vcov) > 0))
  expect_equal(out$se, sqrt(diag(out$vcov)))
  expect_identical(out$hessian_diagnostics$status, "available")
  expect_identical(out$hessian_diagnostics$failure_codes, character())
  expect_identical(out$hessian_diagnostics$hessian_source, "analytical")
  expect_identical(out$hessian_diagnostics$validation_profile, "standard")
  expect_true(all(c(
    "relative_symmetry_error", "min_information_eigenvalue",
    "max_information_eigenvalue", "numerical_rank",
    "scaled_condition_number", "scaled_gradient_max",
    "covariance_diagonal"
  ) %in% names(out$hessian_diagnostics)))
})

test_that("hard curvature boundaries are classed inference failures", {
  expect_inference_failure(
    gamlss.longitudinal:::.gl_solve_hessian_vcov(
      named_hessian(diag(c(-1, NA_real_)))
    ),
    "hessian_nonfinite"
  )
  expect_inference_failure(
    gamlss.longitudinal:::.gl_solve_hessian_vcov(
      named_hessian(matrix(c(-1, 0, 0, 0), 2))
    ),
    "information_not_positive_definite"
  )
  expect_inference_failure(
    gamlss.longitudinal:::.gl_solve_hessian_vcov(
      named_hessian(matrix(c(-1, 0, 0, 1), 2))
    ),
    "information_not_positive_definite"
  )
  expect_inference_failure(
    gamlss.longitudinal:::.gl_solve_hessian_vcov(
      named_hessian(matrix(c(-1, 0.1, 0.2, -1), 2))
    ),
    "hessian_asymmetric"
  )
  expect_inference_failure(
    gamlss.longitudinal:::.gl_solve_hessian_vcov(
      named_hessian(-matrix(c(1, 1 - 1e-13, 1 - 1e-13, 1), 2))
    ),
    "information_ill_conditioned"
  )
})

test_that("scaled conditioning and fitted-score checks are invariant to units", {
  H1 <- named_hessian(-matrix(c(4, 1, 1, 9), 2))
  D <- diag(c(100, 0.01))
  H2 <- named_hessian(D %*% H1 %*% D)
  g1 <- c(0.01, -0.02)
  g2 <- as.numeric(D %*% g1)
  scaled1 <- max(abs(g1) / sqrt(diag(-H1)))
  scaled2 <- max(abs(g2) / sqrt(diag(-H2)))

  out1 <- gamlss.longitudinal:::.gl_validate_hessian_inference(
    H1, parameter_names = colnames(H1),
    gradient = list(gradient = g1, steps = c(1, 1), scaled_max = scaled1)
  )
  out2 <- gamlss.longitudinal:::.gl_validate_hessian_inference(
    H2, parameter_names = colnames(H2),
    gradient = list(gradient = g2, steps = c(1, 1), scaled_max = scaled2)
  )

  expect_equal(out1$hessian_diagnostics$scaled_condition_number,
               out2$hessian_diagnostics$scaled_condition_number,
               tolerance = 1e-10)
  expect_equal(out1$hessian_diagnostics$scaled_gradient_max,
               out2$hessian_diagnostics$scaled_gradient_max,
               tolerance = 1e-12)
})

test_that("high fitted score and Hessian disagreement are explicit failures", {
  H <- named_hessian(-diag(c(2, 3)))
  expect_inference_failure(
    gamlss.longitudinal:::.gl_validate_hessian_inference(
      H, parameter_names = colnames(H),
      gradient = list(gradient = c(1, 1), steps = c(1, 1), scaled_max = 1)
    ),
    "fitted_gradient_too_large"
  )

  strict <- inference_control("strict")
  expect_inference_failure(
    gamlss.longitudinal:::.gl_validate_hessian_inference(
      H, parameter_names = colnames(H), control = strict,
      source = "analytical", reference_hessian = H * 2
    ),
    "hessian_method_disagreement"
  )
})

test_that("summary refuses nonpositive covariance rather than masking it", {
  object <- structure(
    list(par = c("mu.x" = 1, "sigma.x" = 0)),
    class = "gamlss.longitudinal"
  )
  bad_vcov <- list(
    vcov = list(overall = named_hessian(diag(c(1, -1)))),
    se = NULL,
    hessian_diagnostics = list(status = "available")
  )
  expect_inference_failure(
    gamlss.longitudinal:::.gl_summary_coefficient_table(object, bad_vcov),
    "derived_variance_nonpositive"
  )
})

test_that("representative BCPE/t and NBI/Clayton paths pass standard validation", {
  skip_on_cran()
  skip_if_not_installed("gamlss")
  skip_if_not_installed("gamlss.dist")

  cases <- list(
    list(
      seed = 812L,
      margin = gamlss.dist::BCPE(),
      copula = "t",
      margin_params = list(mu = 10, sigma = 0.2, nu = 1, tau = 2),
      copula_params = list(theta = 0.3, zeta = 8),
      formulas = list(nu.formula = "~ 1", tau.formula = "~ 1", zeta.formula = "~ 1")
    ),
    list(
      seed = 813L,
      margin = gamlss.dist::NBI(),
      copula = "C",
      margin_params = list(mu = 3, sigma = 0.4),
      copula_params = list(theta = 1),
      formulas = list(nu.formula = "~ 1", tau.formula = "~ 1", zeta.formula = "~ 1")
    )
  )

  for (case in cases) {
    dat <- simulate_longitudinal_dataset(
      n = 40, times = 1:3,
      margin_dist = case$margin,
      copula_dist = case$copula,
      margin_params = case$margin_params,
      copula_params = case$copula_params,
      seed = case$seed
    )
    fit_args <- c(list(
      dataset = dat,
      margin_dist = case$margin,
      copula_dist = case$copula,
      time_var = "time",
      subject_var = "subject",
      mu.formula = "response ~ 1",
      sigma.formula = "~ 1",
      theta.formula = "~ 1",
      compute_vcov = FALSE,
      warm_start_joint = FALSE,
      max_outer_iter = 8L,
      max_inner_iter = 8L,
      outer_stop_crit = 0.05,
      inner_stop_crit = 0.05,
      verbose = 0
    ), case$formulas)
    invisible(utils::capture.output(
      fit <- suppressWarnings(do.call(gamlss_longitudinal, fit_args))
    ))
    vc <- vcov(fit, method = "numderiv", progress = FALSE)

    expect_identical(vc$hessian_diagnostics$status, "available")
    expect_identical(vc$hessian_diagnostics$validation_profile, "standard")
    expect_true(all(diag(vc$vcov$overall) > 0))
    expect_true(is.finite(vc$hessian_diagnostics$scaled_gradient_max))
  }
})
