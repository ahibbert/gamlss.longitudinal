test_that("all covariance Hessian compatibility routes use the shared correction", {
  old <- options(gamlss.longitudinal.copula_backend = "native")
  on.exit(options(old), add = TRUE)
  for (family in c("N", "C", "F", "G", "J", "t")) {
    theta <- if (family %in% c("N", "t")) 0.6 else 2
    args <- list(c(1e-8, 0.2), c(0.37, 0.7), family, theta, 5,
      c(TRUE, FALSE), family == "t")
    expect_identical(do.call(.gl_covariance_gaussian_pair_derivatives, args),
      do.call(.copula_hessian_pair_derivatives, args))
  }
  delegate <- .gl_covariance_analytical_hessian
  fake <- new.env(parent = environment(delegate))
  fake$calc_analytical_hessian <- function(object, progress, h)
    list(object = object, progress = progress, h = h)
  environment(delegate) <- fake
  for (family in c("N", "C", "F", "G", "J", "t")) {
    object <- list(copula_dist = family)
    expect_identical(delegate(object, FALSE, 2e-4),
      fake$calc_analytical_hessian(object, FALSE, 2e-4))
  }
})

test_that("new markers reject affected caches and preserve unrelated inference", {
  old <- options(gamlss.longitudinal.copula_backend = "native")
  on.exit(options(old), add = TRUE)
  covariance <- list(vcov = list(overall = diag(2)),
    hessian_diagnostics = list(status = "available"), method = "analytical")
  fit <- structure(list(copula_dist = "N", margin_dist = gamlss.dist::NO(),
    vcov = covariance, vcov_meta = list(method = "analytical", numderiv = FALSE,
      cache_version = .gl_vcov_cache_version())), class = "gamlss.longitudinal")
  for (family in c("N", "C", "F", "G", "J", "t")) {
    fit$copula_dist <- family
    fit$vcov$continuous_curvature <- NULL
    expect_false(.can_use_cached_vcov(fit, method = "analytical"))
    fit$vcov$gaussian_mixed_curvature_version <- 1L
    expect_false(.can_use_cached_vcov(fit, method = "analytical"))
    fit$vcov$continuous_curvature <- list(version = 1L, backend = "native", copula = family)
    expect_false(.can_use_cached_vcov(fit, method = "analytical"))
    fit$vcov$continuous_curvature <- list(version = 2L, backend = "native", copula = family)
    expect_false(.can_use_cached_vcov(fit, method = "analytical"))
    fit$vcov$continuous_curvature <- .gl_vcov_continuous_curvature_marker(fit)
    expect_true(.can_use_cached_vcov(fit, method = "analytical"))
    fit$vcov$continuous_curvature$boundary_version <- NULL
    expect_false(.can_use_cached_vcov(fit, method = "analytical"))
    fit$vcov$continuous_curvature <- .gl_vcov_continuous_curvature_marker(fit)
    fit$vcov$continuous_curvature$boundary_version <- 1L
    expect_false(.can_use_cached_vcov(fit, method = "analytical"))
    fit$vcov$continuous_curvature <- .gl_vcov_continuous_curvature_marker(fit)
    options(gamlss.longitudinal.copula_backend = "vinecopula")
    expect_false(.can_use_cached_vcov(fit, method = "analytical"))
    options(gamlss.longitudinal.copula_backend = "native")
  }
  fit$margin_dist <- gamlss.dist::NBI()
  fit$vcov$continuous_curvature <- NULL
  expect_true(.can_use_cached_vcov(fit, method = "analytical"))
  fit$margin_dist <- gamlss.dist::NO()
  fit$vcov$method <- fit$vcov_meta$method <- "numderiv"
  fit$vcov_meta$numderiv <- TRUE
  expect_true(.can_use_cached_vcov(fit, numderiv = TRUE, method = "numderiv"))
  fit$vcov_meta$numderiv <- FALSE
  fit$vcov$method <- fit$vcov_meta$method <- "sandwich_cluster"
  fit$vcov$hessian_diagnostics$bread_method <- "numderiv"
  expect_true(.can_use_cached_vcov(fit, method = "sandwich"))
  fit$vcov$hessian_diagnostics$bread_method <- "analytical"
  expect_false(.can_use_cached_vcov(fit, method = "sandwich"))
  fit$vcov$continuous_curvature <- .gl_vcov_continuous_curvature_marker(fit)
  expect_true(.can_use_cached_vcov(fit, method = "sandwich"))
  fit$vcov$hessian_diagnostics$status <- "unavailable"
  expect_false(.can_use_cached_vcov(fit, method = "sandwich"))
  expect_identical(.gl_vcov_cache_version(), 2L)
})

test_that("new covariance results carry only the applicable curvature marker", {
  build <- .gl_vcov_build_result
  fake <- new.env(parent = environment(build))
  fake$.gl_compute_smooth_vcov <- function(...) list(smooth_vcov = list(), smooth_se = list())
  fake$.gl_enrich_vcov_contract <- function(out, object) out
  environment(build) <- fake
  object <- list(copula_dist = "t", margin_dist = gamlss.dist::NO())
  before <- object
  result <- function(method, diagnostics = list(status = "available"))
    build(object, list(), numeric(), diag(2), c(1, 1), method, method, diagnostics)
  expect_identical(result("analytical")$continuous_curvature,
    .gl_vcov_continuous_curvature_marker(object))
  expect_null(result("numderiv")$continuous_curvature)
  expect_null(result("sandwich_cluster", list(bread_method = "numderiv"))$continuous_curvature)
  expect_identical(result("sandwich_cluster", list(bread_method = "analytical"))$continuous_curvature,
    .gl_vcov_continuous_curvature_marker(object))
  expect_identical(object, before)
  object$margin_dist <- gamlss.dist::NBI()
  expect_null(result("analytical")$continuous_curvature)
})

test_that("stale caches request recomputation without mutating the saved object", {
  resolve <- .resolve_vcov
  fake <- new.env(parent = environment(resolve))
  calls <- 0L
  fake$vcov.gamlss.longitudinal <- function(object, numderiv, details, method) {
    calls <<- calls + 1L
    list(synthetic_recomputation = TRUE, method = method)
  }
  environment(resolve) <- fake
  object <- structure(list(copula_dist = "G", margin_dist = gamlss.dist::NO(),
    vcov = list(vcov = list(overall = diag(2)), method = "analytical"),
    vcov_meta = list(method = "analytical", numderiv = FALSE,
      cache_version = .gl_vcov_cache_version())), class = "gamlss.longitudinal")
  before <- object
  expect_true(resolve(object)$synthetic_recomputation)
  expect_identical(calls, 1L)
  expect_identical(object, before)
})
