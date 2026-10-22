test_that("complete Hessians match an independent likelihood for all six copulas", {
  old <- options(gamlss.longitudinal.copula_backend = "native")
  on.exit(options(old), add = TRUE)
  for (tail in c(FALSE, TRUE)) for (family in c("N", "C", "F", "G", "J", "t")) {
    object <- curvature_model_fixture(family, tail)
    original <- object
    ref <- curvature_reference_hessian(object)
    for (backend in c("native", "vinecopula")) {
      options(gamlss.longitudinal.copula_backend = backend)
      actual <- calc_analytical_hessian(object, progress = FALSE)
      error <- sqrt(sum((actual - ref)^2)) / sqrt(sum(ref^2))
      if (exists("validation_details", inherits = TRUE)) {
        validation_details$hessians[[paste(backend, family, tail)]] <-
          list(backend = backend, family = family, tail = tail,
            relative_error = error, actual = actual, reference = ref)
      }
      expect_true(is.finite(error) && error < 2e-4,
        info = paste(backend, family, tail, error))
      expect_identical(object, original)
    }
  }
})

test_that("stationary Gaussian model and sandwich covariance match independent algebra", {
  old <- options(gamlss.longitudinal.copula_backend = "native")
  on.exit(options(old), add = TRUE)
  object <- curvature_model_fixture()
  before <- object
  ref <- curvature_reference_covariance(object)
  expect_lt(max(abs(colSums(ref$scores))), 1e-5)
  for (backend in c("native", "vinecopula")) {
    options(gamlss.longitudinal.copula_backend = backend)
    for (method in c("analytical_only", "sandwich")) {
      actual <- stats::vcov(object, method = method, details = TRUE,
        sandwich_bread_method = "analytical_only", progress = FALSE)
      expected <- if (method == "sandwich") ref$sandwich else ref$model
      if (exists("validation_details", inherits = TRUE)) {
        validation_details$covariance[[paste(backend, method)]] <-
          list(backend = backend, method = method, actual = actual$vcov$overall,
            reference = expected, diagnostics = actual$hessian_diagnostics)
      }
      expect_identical(dim(actual$vcov$overall), dim(expected))
      expect_identical(rownames(actual$vcov$overall), names(object$par))
      expect_equal(as.numeric(actual$vcov$overall), as.numeric(expected), tolerance = 2e-5)
      expect_equal(unname(actual$se$overall), unname(sqrt(diag(expected))), tolerance = 2e-5)
      expect_identical(actual$continuous_curvature,
        .gl_vcov_continuous_curvature_marker(object))
      expect_identical(object, before)
    }
  }
})
