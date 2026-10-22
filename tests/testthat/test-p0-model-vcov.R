make_minimal_vcov_object <- function() {
  mm <- list(
    x = list(
      mu = matrix(1, nrow = 2, ncol = 1, dimnames = list(NULL, "mu.(Intercept)")),
      sigma = matrix(1, nrow = 2, ncol = 1, dimnames = list(NULL, "sigma.(Intercept)")),
      theta = matrix(1, nrow = 1, ncol = 1, dimnames = list(NULL, "theta.(Intercept)"))
    ),
    s = list(mu = list(), sigma = list(), theta = list())
  )

  list(
    response = c(1, 2),
    response_margin = c(1, 2),
    response_subject = c(1, 1),
    margin_dist = gamlss.dist::NO(),
    copula_dist = "N",
    model_matrix = mm,
    par = c(`mu.(Intercept)` = 0, `sigma.(Intercept)` = 0, `theta.(Intercept)` = 0),
    par_s = list(mu = NA, sigma = NA, theta = NA)
  )
}

test_that("vcov setup prepares fitted-object state and eta values", {
  object <- make_minimal_vcov_object()

  out <- gamlss.longitudinal:::.gl_prepare_vcov_evaluation(
    object = object,
    par = NA,
    numderiv = FALSE,
    method = c("analytical", "numderiv", "analytical_only"),
    progress = TRUE
  )

  expect_equal(out$method, "analytical")
  expect_equal(out$method_requested, "analytical")
  expect_equal(out$method_used, "analytical")
  expect_true(out$progress)
  expect_true(out$include_dlcopdpar)
  expect_identical(out$response, object$response)
  expect_identical(out$response_margin, object$response_margin)
  expect_identical(out$response_subject, object$response_subject)
  expect_equal(out$margin_names, c(1, 2))
  expect_equal(out$num_margins, 2L)
  expect_identical(out$mm, object$model_matrix)
  expect_equal(out$par_cov, object$par)
  expect_named(out$eta, c("mu", "sigma", "theta"))
  expect_equal(unname(out$eta$mu), c(0, 0))
  expect_equal(unname(out$eta$sigma), c(0, 0))
  expect_length(out$eta_inv$theta, 1L)
})

test_that("cluster log-likelihood and scores aggregate fitted likelihood", {
  dat <- make_fixture_factor_time_interaction(n_subject = 8L)
  fit <- fit_fixture_model(
    dat,
    include_dlcopdpar = TRUE,
    mu_formula = "y ~ age",
    sigma_formula = "~ 1",
    theta_formula = "~ 1",
    max_outer_iter = 2,
    max_inner_iter = 2
  )

  contrib <- gamlss.longitudinal:::.gl_cluster_joint_loglik_contributions(fit)
  expect_equal(
    sum(contrib),
    unname(fit$calc_lik_out_end$log_lik["joint"]),
    tolerance = 1e-6
  )

  scores <- gamlss.longitudinal:::.gl_cluster_score_matrix(fit, h = 1e-4)
  expect_equal(nrow(scores), length(unique(fit$response_subject)))
  expect_equal(colnames(scores), names(fit$par))
  expect_true(all(is.finite(scores)))

  step <- 1e-4 * max(1, abs(fit$par[[1]]))
  plus <- minus <- fit$par
  plus[[1]] <- plus[[1]] + step
  minus[[1]] <- minus[[1]] - step
  total_score <- (
    sum(gamlss.longitudinal:::.gl_cluster_joint_loglik_contributions(fit, par_cov = plus)) -
      sum(gamlss.longitudinal:::.gl_cluster_joint_loglik_contributions(fit, par_cov = minus))
  ) / (2 * step)
  expect_equal(sum(scores[, 1]), total_score, tolerance = 1e-8)
})

test_that("sandwich vcov is available through inference helpers", {
  dat <- make_fixture_factor_time_interaction(n_subject = 8L)
  fit <- fit_fixture_model(
    dat,
    include_dlcopdpar = TRUE,
    mu_formula = "y ~ age",
    sigma_formula = "~ 1",
    theta_formula = "~ 1",
    max_outer_iter = 2,
    max_inner_iter = 2
  )

  vc <- vcov(
    fit,
    method = "sandwich",
    sandwich_h = 1e-4,
    sandwich_bread_method = "analytical"
  )
  expect_equal(vc$method, "sandwich_cluster")
  expect_equal(rownames(vc$vcov$overall), names(fit$par))
  expect_true(all(is.finite(vc$se$overall)))
  expect_equal(vc$vcov$cluster_robust, vc$vcov$overall)
  expect_equal(vc$se$cluster_robust, vc$se$overall)
  expect_equal(rownames(vc$vcov$model_based), names(fit$par))
  expect_equal(vc$se$model_based, sqrt(diag(vc$vcov$model_based)))
  expect_equal(vc$hessian_diagnostics$n_clusters, length(unique(fit$response_subject)))

  s <- summary(fit, include_vcov = TRUE, vcov_method = "sandwich", sandwich_h = 1e-4)
  expect_true(all(c(
    "std_error", "model_std_error", "robust_std_error", "p_value"
  ) %in% names(s$coefficients)))
  expect_true(any(is.finite(s$coefficients$std_error)))
  expect_equal(s$coefficients$std_error, s$coefficients$robust_std_error)
  expect_true(all(is.finite(s$coefficients$model_std_error)))
  expect_equal(
    s$coefficients$p_value,
    2 * stats::pnorm(abs(s$coefficients$estimate / s$coefficients$robust_std_error),
      lower.tail = FALSE
    )
  )

  ci <- confint(fit, method = "sandwich", sandwich_h = 1e-4)
  expect_equal(nrow(ci), length(fit$par))
  expect_true(all(is.finite(ci)))

  wt <- wald_test(fit, terms = names(fit$par)[1], method = "sandwich", sandwich_h = 1e-4)
  expect_s3_class(wt, "gamlss_longitudinal_wald_test")
  expect_equal(wt$method, "sandwich_cluster")

  fit$likelihood_contract <- list(objective = "segmented")
  fit$missingness <- list(
    objective = "segmented",
    n_subjects_with_gaps = 1L,
    n_segments = length(unique(fit$response_subject)) + 1L
  )
  warning_seen <- FALSE
  segmented_summary <- withCallingHandlers(
    summary(fit, sandwich_h = 1e-4),
    gamlss_longitudinal_segmented_inference_warning = function(w) {
      warning_seen <<- TRUE
      invokeRestart("muffleWarning")
    }
  )
  expect_true(warning_seen)
  expect_s3_class(segmented_summary, "summary.gamlss.longitudinal")
  expect_identical(segmented_summary$fit$vcov_method, "sandwich_cluster")
  expect_identical(segmented_summary$fit$vcov_method_requested, "analytical")
  expect_true(all(is.finite(segmented_summary$coefficients$std_error)))
  segmented_print <- paste(capture.output(print(segmented_summary)), collapse = "\n")
  expect_match(
    segmented_print,
    "subject-clustered sandwich covariance",
    fixed = TRUE
  )
  expect_match(segmented_print, "model_se", fixed = TRUE)
  expect_match(segmented_print, "robust_se", fixed = TRUE)
  expect_match(
    segmented_print,
    "Tests and p-values use subject-clustered robust standard errors",
    fixed = TRUE
  )
})

test_that("summary reuses fit-time sandwich covariance for a segmented fit", {
  dat <- make_fixture_factor_time_interaction(n_subject = 40L)
  gap_ids <- seq(5L, 40L, by = 10L)
  keep <- !(dat$id %in% gap_ids & dat$time_raw == levels(dat$time_raw)[2L])
  dat <- dat[keep, , drop = FALSE]

  fit <- suppressWarnings(fit_fixture_model(
    dat,
    missingness = "segment",
    include_dlcopdpar = FALSE,
    mu_formula = "y ~ 1",
    sigma_formula = "~ 1",
    theta_formula = "~ 1",
    compute_vcov = TRUE,
    max_outer_iter = 4L,
    max_inner_iter = 3L,
    outer_stop_crit = 10,
    inner_stop_crit = 10
  ))

  expect_true(fit$convergence$converged)
  expect_true(fit$vcov_meta$precomputed)
  expect_identical(fit$vcov_meta$method, "analytical")
  expect_identical(fit$vcov_meta$method_used, "sandwich_cluster")

  testthat::local_mocked_bindings(
    vcov.gamlss.longitudinal = function(...) {
      stop("summary recomputed covariance instead of using the fit-time cache")
    },
    .package = "gamlss.longitudinal"
  )
  summary_out <- summary(fit)
  explicit_summary_out <- summary(fit, vcov_method = "sandwich")

  expect_s3_class(summary_out, "summary.gamlss.longitudinal")
  expect_s3_class(explicit_summary_out, "summary.gamlss.longitudinal")
  expect_identical(summary_out$fit$vcov_method, "sandwich_cluster")
  expect_identical(explicit_summary_out$fit$vcov_method, "sandwich_cluster")
  expect_true(all(is.finite(summary_out$coefficients$std_error)))
  expect_equal(
    summary_out$coefficients$std_error,
    summary_out$coefficients$robust_std_error
  )
  expect_true(all(is.finite(summary_out$coefficients$model_std_error)))
})

test_that("vcov setup preserves legacy numderiv and parameter override behavior", {
  object <- make_minimal_vcov_object()
  override <- list(
    par = c(`mu.(Intercept)` = 1, `sigma.(Intercept)` = 0.5, `theta.(Intercept)` = 0.25),
    par_s = list(mu = NA, sigma = NA, theta = NA)
  )

  out <- gamlss.longitudinal:::.gl_prepare_vcov_evaluation(
    object = object,
    par = override,
    numderiv = TRUE,
    method = c("analytical", "numderiv", "analytical_only"),
    progress = FALSE
  )

  expect_equal(out$method, "numderiv")
  expect_equal(out$method_requested, "numderiv")
  expect_equal(out$method_used, "numderiv")
  expect_false(out$progress)
  expect_equal(out$par_cov, override$par)
  expect_equal(unname(out$eta$mu), c(1, 1))
  expect_equal(unname(out$eta$sigma), c(0.5, 0.5))
  expect_equal(unname(out$eta$theta), 0.25)
})

test_that("current vcov solve methods require valid signed likelihood curvature", {
  hessian <- -diag(c(2, 4))
  dimnames(hessian) <- list(c("a", "b"), c("a", "b"))

  out <- gamlss.longitudinal:::.gl_vcov_solve_if_needed(
    vcov_path = list(
      hessian_nd = hessian,
      vcov_final = NULL,
      se_final = NULL,
      hessian_diagnostics = NULL
    ),
    method = "numderiv",
    d2_mat = NULL,
    response = c(1, 2, 3),
    gradient = list(
      gradient = c(a = 0, b = 0),
      steps = c(1e-4, 1e-4),
      scaled_max = 0
    )
  )

  expect_equal(out$vcov_final, solve(-hessian))
  expect_equal(out$se_final, sqrt(diag(solve(-hessian))))
  expect_true(is.list(out$hessian_diagnostics))
})

test_that("analytical vcov fallback uses a guarded larger numerical Hessian step", {
  parameter_names <- c("mu.(Intercept)", "sigma.(Intercept)")
  calls <- numeric()
  hessian <- -diag(c(2, 4))
  dimnames(hessian) <- list(parameter_names, parameter_names)
  object <- list(par = stats::setNames(c(0, 0), parameter_names))

  testthat::local_mocked_bindings(
    .gl_source_analytical_hessian_helpers = function() invisible(TRUE),
    calc_analytical_hessian = function(...) stop("analytical path failed"),
    calc_true_SE_numderiv_only_covariates = function(..., h = NULL, progress = NULL) {
      calls <<- c(calls, h)
      hessian
    },
    .gl_vcov_fitted_gradient = function(...) {
      list(
        gradient = stats::setNames(c(0, 0), parameter_names),
        steps = c(1e-4, 1e-4), scaled_max = 0
      )
    },
    .package = "gamlss.longitudinal"
  )

  expect_warning(
    out <- gamlss.longitudinal:::.gl_vcov_compute_primary(
      object = object,
      par_cov = object$par,
      mm = list(x = list()),
      margin_dist = NULL,
      response = numeric(),
      response_margin = numeric(),
      response_subject = numeric(),
      method = "analytical",
      progress = FALSE,
      h = 1e-4,
      inference = inference_control("standard")
    ),
    "Analytical Hessian vcov failed",
    fixed = TRUE
  )

  expect_equal(calls, 1e-3)
  expect_identical(out$method_used, "numderiv")
  expect_equal(out$hessian_fallback$numerical_hessian_step, 1e-3)
  expect_equal(out$vcov_final, solve(-hessian))
  expect_identical(out$hessian_diagnostics$status, "available")
})

test_that("vcov details retain guarded numerical-Hessian fallback provenance", {
  testthat::local_mocked_bindings(
    .gl_compute_smooth_vcov = function(...) list(smooth_vcov = list(), smooth_se = list()),
    .gl_enrich_vcov_contract = function(out, object) out,
    .package = "gamlss.longitudinal"
  )
  fallback <- list(
    from = "analytical", to = "numderiv", numerical_hessian_step = 1e-3
  )
  out <- gamlss.longitudinal:::.gl_vcov_build_result(
    object = list(), eta_inv = list(), response = numeric(),
    vcov_final = diag(2), se_final = c(1, 1), method_used = "numderiv",
    method_requested = "analytical", hessian_diagnostics = list(status = "available"),
    hessian_fallback = fallback
  )
  expect_identical(out$hessian_fallback, fallback)
})
