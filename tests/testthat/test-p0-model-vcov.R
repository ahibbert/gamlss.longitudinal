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
  expect_equal(vc$hessian_diagnostics$n_clusters, length(unique(fit$response_subject)))

  s <- summary(fit, include_vcov = TRUE, vcov_method = "sandwich", sandwich_h = 1e-4)
  expect_true(all(c("std_error", "p_value") %in% names(s$coefficients)))
  expect_true(any(is.finite(s$coefficients$std_error)))

  ci <- confint(fit, method = "sandwich", sandwich_h = 1e-4)
  expect_equal(nrow(ci), length(fit$par))
  expect_true(all(is.finite(ci)))

  wt <- wald_test(fit, terms = names(fit$par)[1], method = "sandwich", sandwich_h = 1e-4)
  expect_s3_class(wt, "gamlss_longitudinal_wald_test")
  expect_equal(wt$method, "sandwich_cluster")
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
