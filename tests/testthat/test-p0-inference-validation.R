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
  expect_identical(standard$information_eigenvalue_min, 0)
  expect_true(strict$check_agreement)
  expect_lt(strict$symmetry_tol, standard$symmetry_tol)
  expect_equal(overridden$condition_max, 123)
  expect_equal(overridden$expert_override, list(condition_max = 123))
  expect_identical(standard$defaults_version, "0.1.0-provisional")
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

test_that("signed information eigenvalue threshold is registered and strict", {
  control <- inference_control("standard")
  expect_identical(control$information_eigenvalue_min, 0)

  H_positive <- matrix(-.Machine$double.eps, 1L, 1L,
                       dimnames = list("mu.x", "mu.x"))
  available <- gamlss.longitudinal:::.gl_validate_hessian_inference(
    H_positive,
    parameter_names = "mu.x",
    control = control,
    gradient = list(gradient = c("mu.x" = 0), steps = 1e-4, scaled_max = 0)
  )
  expect_identical(available$hessian_diagnostics$status, "available")
  expect_identical(
    available$hessian_diagnostics$effective_thresholds$information_eigenvalue_min,
    0
  )

  H_boundary <- matrix(0, 1L, 1L, dimnames = list("mu.x", "mu.x"))
  expect_inference_failure(
    gamlss.longitudinal:::.gl_validate_hessian_inference(
      H_boundary,
      parameter_names = "mu.x",
      control = control,
      gradient = list(gradient = c("mu.x" = 0), steps = 1e-4, scaled_max = 0)
    ),
    "information_not_positive_definite"
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

test_that("inference validation has tested inclusive threshold boundaries", {
  H <- named_hessian(-diag(2))
  ctl <- inference_control("standard", condition_max = 1, gradient_tol = 0.25)
  at_boundary <- gamlss.longitudinal:::.gl_validate_hessian_inference(
    H,
    parameter_names = colnames(H),
    control = ctl,
    gradient = list(
      gradient = c(0.25, 0),
      steps = c(1e-4, 1e-4),
      scaled_max = 0.25
    )
  )
  expect_identical(at_boundary$hessian_diagnostics$status, "available")
  expect_equal(at_boundary$hessian_diagnostics$scaled_condition_number, 1)

  H_condition <- named_hessian(-matrix(c(1, 0.5, 0.5, 1), 2))
  condition_boundary <- kappa(-H_condition, exact = TRUE)
  at_condition <- gamlss.longitudinal:::.gl_validate_hessian_inference(
    H_condition,
    parameter_names = colnames(H_condition),
    control = inference_control("standard", condition_max = condition_boundary),
    gradient = list(
      gradient = c(0, 0), steps = c(1e-4, 1e-4), scaled_max = 0
    )
  )
  expect_equal(
    at_condition$hessian_diagnostics$scaled_condition_number,
    condition_boundary
  )
  expect_inference_failure(
    gamlss.longitudinal:::.gl_validate_hessian_inference(
      H_condition,
      parameter_names = colnames(H_condition),
      control = inference_control(
        "standard", condition_max = condition_boundary - 1e-10
      ),
      gradient = list(
        gradient = c(0, 0), steps = c(1e-4, 1e-4), scaled_max = 0
      )
    ),
    "information_ill_conditioned"
  )

  expect_inference_failure(
    gamlss.longitudinal:::.gl_validate_hessian_inference(
      H,
      parameter_names = colnames(H),
      control = ctl,
      gradient = list(
        gradient = c(0.250001, 0),
        steps = c(1e-4, 1e-4),
        scaled_max = 0.250001
      )
    ),
    "fitted_gradient_too_large"
  )
  expect_inference_failure(
    gamlss.longitudinal:::.gl_validate_hessian_inference(
      H,
      parameter_names = colnames(H),
      control = ctl,
      gradient = NULL
    ),
    "fitted_gradient_not_checked"
  )
})

test_that("symmetry, rank, and Hessian-agreement boundaries are inclusive", {
  gradient <- list(gradient = c(0, 0), steps = c(1e-4, 1e-4), scaled_max = 0)

  symmetry_tol <- 0.01
  H_at_symmetry <- named_hessian(matrix(c(-1, 0, symmetry_tol, -1), 2))
  at_symmetry <- gamlss.longitudinal:::.gl_validate_hessian_inference(
    H_at_symmetry, parameter_names = colnames(H_at_symmetry), gradient = gradient,
    control = inference_control("standard", symmetry_tol = symmetry_tol)
  )
  expect_equal(at_symmetry$hessian_diagnostics$relative_symmetry_error, symmetry_tol)
  expect_inference_failure(
    gamlss.longitudinal:::.gl_validate_hessian_inference(
      named_hessian(matrix(c(-1, 0, symmetry_tol + 1e-6, -1), 2)),
      parameter_names = colnames(H_at_symmetry), gradient = gradient,
      control = inference_control("standard", symmetry_tol = symmetry_tol)
    ),
    "hessian_asymmetric"
  )

  H_rank <- named_hessian(-matrix(c(1, 0.5, 0.5, 1), 2))
  at_rank <- gamlss.longitudinal:::.gl_validate_hessian_inference(
    H_rank, parameter_names = colnames(H_rank), gradient = gradient,
    control = inference_control("standard", rank_tol = 0.6, condition_max = 10)
  )
  expect_identical(at_rank$hessian_diagnostics$numerical_rank, 2L)
  expect_inference_failure(
    gamlss.longitudinal:::.gl_validate_hessian_inference(
      H_rank, parameter_names = colnames(H_rank), gradient = gradient,
      control = inference_control("standard", rank_tol = 0.600001, condition_max = 10)
    ),
    "information_rank_deficient"
  )

  agreement_tol <- 0.25
  H <- named_hessian(-diag(2))
  at_agreement <- gamlss.longitudinal:::.gl_validate_hessian_inference(
    H, parameter_names = colnames(H), gradient = gradient,
    source = "analytical", reference_hessian = H * (1 + agreement_tol),
    control = inference_control("strict", agreement_tol = agreement_tol)
  )
  expect_equal(
    at_agreement$hessian_diagnostics$hessian_agreement_relative,
    agreement_tol
  )
  expect_inference_failure(
    gamlss.longitudinal:::.gl_validate_hessian_inference(
      H, parameter_names = colnames(H), gradient = gradient,
      source = "analytical", reference_hessian = H * (1 + agreement_tol + 1e-6),
      control = inference_control("strict", agreement_tol = agreement_tol)
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

test_that("classed inference failure propagates through all inference consumers", {
  dat <- make_fixture_factor_time_interaction(n_subject = 8L)
  fit <- fit_fixture_model(
    dat,
    include_dlcopdpar = TRUE,
    mu_formula = "y ~ age",
    sigma_formula = "~ 1",
    theta_formula = "~ 1",
    max_outer_iter = 2L,
    max_inner_iter = 2L,
    compute_vcov = FALSE
  )
  unavailable <- list(
    status = "unavailable",
    failure_codes = "information_not_positive_definite",
    validation_profile = "standard",
    validation_defaults_version = "0.1.0-provisional"
  )
  fit$vcov <- NULL
  fit$vcov_meta <- list(
    precomputed = FALSE,
    method = "analytical",
    inference_status = "unavailable",
    hessian_diagnostics = unavailable
  )

  consumers <- list(
    confint = function() confint(fit),
    wald_test = function() wald_test(fit, terms = names(fit$par)[1]),
    prediction_se = function() predict(fit, type = "mu", se.fit = TRUE),
    prediction_interval = function() predict(
      fit, type = "mu", interval = "confidence"
    ),
    summary = function() summary(fit, include_vcov = TRUE),
    publication_table = function() publication_table(
      fit, table = "coefficients", output = "data.frame"
    ),
    term_plot = function() plot_terms(fit, data = dat),
    marginal_effects = function() marginal_effects(
      fit,
      newdata = dat,
      variable = "age",
      values = stats::quantile(dat$age, c(0.25, 0.75)),
      se.fit = TRUE
    )
  )

  for (consumer_name in names(consumers)) {
    err <- tryCatch(
      suppressMessages(suppressWarnings(
        invisible(utils::capture.output(consumers[[consumer_name]]()))
      )),
      error = identity
    )
    expect_s3_class(
      err,
      "gamlss_longitudinal_inference_unavailable"
    )
    expect_identical(
      err$diagnostics$failure_codes,
      "information_not_positive_definite"
    )
  }
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
