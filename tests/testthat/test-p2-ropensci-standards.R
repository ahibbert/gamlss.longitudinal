test_that("standard regression accessors expose fitted model components", {
  skip_if_not_installed("gamlss")
  skip_if_not_installed("gamlss.dist")

  dat <- make_fixture_with_structural_missing_rows()
  fit <- fit_fixture_model(
    dat,
    include_dlcopdpar = FALSE,
    max_outer_iter = 2,
    theta_formula = "~ 1"
  )

  expect_s3_class(stats::formula(fit), "formula")
  expect_s3_class(stats::terms(fit), "terms")
  expect_s3_class(stats::formula(fit, parameter = "theta"), "formula")
  expect_s3_class(stats::formula(fit, internal = TRUE), "formula")

  expect_equal(stats::nobs(fit), sum(is.finite(fit$response)))
  expect_equal(stats::nobs(fit, type = "expanded"), nrow(fit$dataset))
  expect_equal(stats::nobs(fit, type = "submitted"), nrow(fit$dataset_original))

  mf_expanded <- stats::model.frame(fit, type = "expanded")
  mf_observed <- stats::model.frame(fit, type = "observed")
  expect_true(all(c("subject", "time", "time_covariate", "response") %in% names(mf_expanded)))
  expect_equal(nrow(mf_expanded), stats::nobs(fit, type = "expanded"))
  expect_equal(nrow(mf_observed), stats::nobs(fit))
  expect_equal(rownames(mf_expanded), as.character(seq_len(nrow(mf_expanded))))

  fit_mu <- stats::fitted(fit)
  res_response <- stats::residuals(fit, type = "response")
  res_pearson <- stats::residuals(fit, type = "pearson")
  res_quantile <- stats::residuals(fit, type = "quantile")

  expect_equal(length(fit_mu), length(fit$response))
  expect_true(all(is.finite(stats::fitted(fit, finite = TRUE))))
  expect_equal(length(res_response), stats::nobs(fit))
  expect_equal(length(res_pearson), stats::nobs(fit))
  expect_equal(length(res_quantile), stats::nobs(fit))
  expect_true(all(is.finite(res_response)))
  expect_true(all(is.finite(res_pearson)))
  expect_true(all(is.finite(res_quantile)))

  expect_true(is.list(fit$convergence))
  expect_true(all(c("converged", "stop_reason", "outer_iterations", "method") %in% names(fit$convergence)))
})

test_that("known-truth example data are shipped for user verification", {
  path <- system.file("extdata", "known-truth-longitudinal.csv", package = "gamlss.longitudinal")
  if (!nzchar(path)) {
    path <- file.path("inst", "extdata", "known-truth-longitudinal.csv")
  }
  expect_true(file.exists(path))

  dat <- utils::read.csv(path)
  expect_true(all(c("subject", "time", "x", "response", "true_mu", "true_sigma") %in% names(dat)))
  expect_equal(dat$true_mu, 1 + 0.5 * dat$x)
  expect_equal(unique(dat$true_sigma), 0.2)
})

test_that("rank-deficient fixed-effect matrices are detected", {
  skip_if_not_installed("gamlss")
  skip_if_not_installed("gamlss.dist")

  dat <- make_fixture_factor_time_interaction(n_subject = 10L)
  dat$age_copy <- dat$age

  expect_warning(
    fit_fixture_model(
      dat,
      include_dlcopdpar = FALSE,
      mu_formula = "y ~ age + age_copy",
      sigma_formula = "~ 1",
      theta_formula = "~ 1",
      max_outer_iter = 2
    ),
    "rank deficient"
  )
})

test_that("benchmark report paths receive markdown suffixes", {
  results <- data.frame(
    family = "NO",
    copula = "N",
    design = "intercept",
    n_subject = 3L,
    n_time = 3L,
    dependence = "moderate",
    missingness = "none",
    start_mode = "default",
    method = c("rs_separate", "gee"),
    benchmark_mean_rmse = c(0.1, 0.3)
  )
  out_base <- file.path(tempdir(), "standards-report")
  out_path <- gamlss.longitudinal::write_benchmark_report(
    results,
    path = out_base,
    metrics = "benchmark_mean_rmse"
  )
  expect_match(out_path, "\\.md$")
  expect_true(file.exists(out_path))
  expect_equal(basename(out_path), "standards-report.md")
})
