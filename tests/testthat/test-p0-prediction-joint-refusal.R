test_that("joint prediction uncertainty refuses before model or covariance work", {
  touched <- 0L
  forbidden <- function(...) { touched <<- touched + 1L; stop("unexpected work") }
  testthat::local_mocked_bindings(
    .gl_require_converged_fit = forbidden,
    .gl_fitted_distribution = forbidden,
    .resolve_vcov = forbidden,
    .package = "gamlss.longitudinal"
  )
  object <- structure(list(), class = "gamlss.longitudinal")
  for (type in c("mu", "mean", "response", "median", "parameters",
                 "quantile", "cdf", "density", "probability")) {
    for (request in list(list(se.fit = TRUE), list(interval = "confidence"))) {
      for (scope in c("joint", "j")) {
        err <- expect_error(do.call(stats::predict, c(
          list(object = object, type = type, coefficients = scope), request
        )), class = "gamlss_longitudinal_inference_unavailable")
        expect_identical(err$diagnostics$failure_codes,
                         "prediction_joint_uncertainty_unsupported")
      }
    }
  }
  expect_error(stats::predict(object, se.fit = TRUE, coef = "joint"),
               class = "gamlss_longitudinal_inference_unavailable")
  expect_identical(touched, 0L)
})

test_that("default and explicit fixed prediction retain the fixed gradient calculation", {
  X <- cbind(mu.intercept = 1, mu.x = c(0, 1))
  V <- diag(c(4, 9))
  dimnames(V) <- list(colnames(X), colnames(X))
  scopes <- list()
  object <- structure(list(par = c(mu.intercept = 2, mu.x = 1),
    model_matrix = list(x = list(mu = X)), copula_dist = "N", par_s = list()),
    class = "gamlss.longitudinal")
  testthat::local_mocked_bindings(
    .gl_require_converged_fit = function(...) invisible(NULL),
    .gl_fitted_distribution = function(...) list(params = list(mu = c(2, 3)), family = NULL),
    .gl_prediction_frame = function(...) data.frame(subject = c("a", "b"),
      time = c(1, 1), response = c(NA_real_, NA_real_)),
    .resolve_vcov = function(object, extra_args) {
      scopes <<- c(scopes, list(extra_args$coefficients))
      list(vcov = list(overall = V), inference_contract = list(
        validity_status = "available", method_used = "analytical"))
    },
    .gl_require_available_inference = function(...) invisible(NULL),
    calc_eta = function(...) list(eta_dr = list(mu = c(1, 2))),
    .package = "gamlss.longitudinal"
  )
  ordinary <- stats::predict(object, type = "mu", interval = "confidence")
  explicit <- stats::predict(object, type = "mu", interval = "confidence", coefficients = "fixed")
  expect_identical(ordinary, explicit)
  expect_identical(scopes, list(NULL, "fixed"))
  expected_se <- c(2, 2 * sqrt(13))
  expect_equal(as.numeric(ordinary$se.fit), expected_se)
  expect_equal(ordinary$fit, c(2, 3))
  expect_equal(as.numeric(ordinary$conf.low), c(2, 3) - stats::qnorm(.975) * expected_se)
  expect_equal(as.numeric(ordinary$conf.high), c(2, 3) + stats::qnorm(.975) * expected_se)
  expect_identical(attr(ordinary, "inference_contract")$contract_id, "prediction_mu_delta")
})

test_that("point-only joint scope leaves point predictions unchanged", {
  touched <- 0L
  testthat::local_mocked_bindings(
    .gl_warn_nonconverged_prediction = function(...) invisible(NULL),
    .gl_fitted_distribution = function(...) list(params = list(mu = c(2, 3)), family = NULL),
    .resolve_vcov = function(...) { touched <<- touched + 1L; stop("unexpected covariance") },
    .package = "gamlss.longitudinal"
  )
  object <- structure(list(), class = "gamlss.longitudinal")
  ordinary <- stats::predict(object, type = "mu")
  explicit <- stats::predict(object, type = "mu", coefficients = "joint")
  expect_identical(explicit, ordinary)
  expect_equal(explicit, c(2, 3))
  expect_identical(touched, 0L)
})

test_that("omitted prediction scope reuses a valid current cache without recomputation", {
  X <- cbind(mu.intercept = 1, mu.x = c(0, 1))
  V <- diag(c(4, 9))
  dimnames(V) <- list(colnames(X), colnames(X))
  object <- structure(list(par = c(mu.intercept = 2, mu.x = 1),
    margin_dist = gamlss.dist::NO(), copula_dist = "N", par_s = list(),
    model_matrix = list(x = list(mu = X)),
    vcov = list(method = "analytical", vcov = list(overall = V),
      hessian_diagnostics = list(status = "available")),
    vcov_meta = list(method = "analytical", numderiv = FALSE,
      cache_version = gamlss.longitudinal:::.gl_vcov_cache_version())),
    class = "gamlss.longitudinal")
  object$vcov$continuous_curvature <-
    gamlss.longitudinal:::.gl_vcov_continuous_curvature_marker(object)
  original <- object
  recomputations <- 0L
  testthat::local_mocked_bindings(
    .gl_require_converged_fit = function(...) invisible(NULL),
    .gl_fitted_distribution = function(...) list(params = list(mu = c(2, 3)), family = NULL),
    .gl_prediction_frame = function(...) data.frame(subject = c("a", "b"),
      time = c(1, 1), response = c(NA_real_, NA_real_)),
    calc_eta = function(...) list(eta_dr = list(mu = c(1, 2))),
    vcov.gamlss.longitudinal = function(...) {
      recomputations <<- recomputations + 1L
      stop("fixture refuses covariance recomputation", call. = FALSE)
    },
    .package = "gamlss.longitudinal"
  )
  expect_true(gamlss.longitudinal:::.can_use_cached_vcov(object,
    method = "analytical", extra_args = list(method = "analytical")))
  expect_false(gamlss.longitudinal:::.can_use_cached_vcov(object,
    method = "analytical", extra_args = list(method = "analytical", coefficients = "fixed")))
  # The real resolver, cache predicate, contract and SE calculation run here.
  for (request in list(list(se.fit = TRUE), list(interval = "confidence"))) {
    result <- do.call(stats::predict, c(list(object = object, type = "mu"), request))
    expect_equal(as.numeric(result$se.fit), c(2, 2 * sqrt(13)))
    expect_identical(recomputations, 0L)
  }
  expect_identical(object, original)
  # Explicit fixed already bypassed the cache before this change; preserve it.
  expect_error(stats::predict(object, type = "mu", se.fit = TRUE, coefficients = "fixed"),
    "fixture refuses covariance recomputation", fixed = TRUE)
  expect_identical(recomputations, 1L)
  expect_error(stats::predict(object, type = "mu", se.fit = TRUE, coefficients = "joint"),
    class = "gamlss_longitudinal_inference_unavailable")
  expect_identical(recomputations, 1L)
  expect_identical(object, original)
})
