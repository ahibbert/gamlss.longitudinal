# Fit-free cache routing and metadata tests. These do not validate the Hessian
# correction numerically; the independent chain-rule suite does that separately.
chain_epoch_cache_fixture <- function(discrete = FALSE, method = "analytical",
                                      bread = "analytical", version = 3L) {
  fit <- structure(list(
    copula_dist = "N",
    margin_dist = if (discrete) gamlss.dist::NBI() else gamlss.dist::NO(),
    vcov = list(vcov = list(overall = diag(2)), method = method,
      hessian_diagnostics = list(status = "available", bread_method = bread)),
    vcov_meta = list(method = method, numderiv = identical(method, "numderiv"),
      cache_version = version)), class = "gamlss.longitudinal")
  if (!discrete && .gl_vcov_uses_analytical_curvature(method, fit$vcov$hessian_diagnostics)) {
    fit$vcov$continuous_curvature <- .gl_vcov_continuous_curvature_marker(fit)
  }
  fit
}

test_that("epoch 2 analytical caches are rejected for all margins and copulas", {
  old <- options(gamlss.longitudinal.copula_backend = "native")
  on.exit(options(old), add = TRUE)
  for (backend in c("native", "vinecopula")) {
    options(gamlss.longitudinal.copula_backend = backend)
    for (discrete in c(FALSE, TRUE)) for (family in c("N", "C", "F", "G", "J", "t")) {
      for (method in c("analytical", "analytical_only", "sandwich_cluster")) {
        fit <- chain_epoch_cache_fixture(discrete, method, version = 2L)
        fit$copula_dist <- family
        if (!discrete) fit$vcov$continuous_curvature <- .gl_vcov_continuous_curvature_marker(fit)
        before <- fit
        expect_false(.can_use_cached_vcov(fit, method = method))
        expect_identical(fit, before)
        fit$vcov_meta$cache_version <- 3L
        expect_true(.can_use_cached_vcov(fit, method = method))
      }
    }
  }
})

test_that("numerical caches regenerate conservatively without an analytical claim", {
  for (discrete in c(FALSE, TRUE)) {
    for (method in c("numderiv", "sandwich_cluster")) {
      fit <- chain_epoch_cache_fixture(discrete, method, bread = "numderiv", version = 2L)
      fit$vcov$continuous_curvature <- NULL
      numeric_flag <- identical(method, "numderiv")
      expect_false(.can_use_cached_vcov(fit, numderiv = numeric_flag, method = method))
      fit$vcov_meta$cache_version <- 3L
      expect_true(.can_use_cached_vcov(fit, numderiv = numeric_flag, method = method))
      stamped <- .gl_vcov_stamp_epoch(fit$vcov)
      expect_identical(stamped$inference_epoch, 3L)
      expect_null(stamped$analytical_margin_chain_rule_version)
      expect_identical(stamped$vcov, fit$vcov$vcov)
      expect_identical(stamped$method, method)
    }
  }
})

test_that("missing or malformed metadata cannot bypass epoch or method matching", {
  fit <- chain_epoch_cache_fixture(TRUE)
  expect_true(.can_use_cached_vcov(fit, method = "analytical"))
  for (bad in list(NULL, 3L, matrix(3L), list(),
                  list(cache_version = 3L),
                  list(cache_version = 3L, method = "analytical"),
                  list(cache_version = 3L, numderiv = FALSE))) {
    changed <- fit
    changed$vcov_meta <- bad
    expect_false(.can_use_cached_vcov(changed, method = "analytical"))
  }
  for (bad in list(NA, 0L, c(FALSE, FALSE), "FALSE", logical())) {
    changed <- fit
    changed$vcov_meta$numderiv <- bad
    expect_false(.can_use_cached_vcov(changed, method = "analytical"))
  }
  for (field in c("method", "method_used")) for (bad in list(NA_character_, "", character(), c("analytical", "numderiv"))) {
    changed <- fit
    changed$vcov_meta[[field]] <- bad
    expect_false(.can_use_cached_vcov(changed, method = "analytical"))
  }
  for (bad in list(NULL, diag(2), 1)) {
    changed <- fit
    changed$vcov <- bad
    expect_false(.can_use_cached_vcov(changed, method = "analytical"))
  }
  expect_false(.can_use_cached_vcov(fit, numderiv = TRUE, method = "analytical"))
  expect_false(.can_use_cached_vcov(fit, method = "numderiv"))
  expect_false(.can_use_cached_vcov(fit, method = "analytical", extra_args = list(h = 1e-3)))
  fit$vcov$hessian_diagnostics$status <- "unavailable"
  expect_false(.can_use_cached_vcov(fit, method = "analytical"))
})

test_that("fixed resolver recomputes old caches and preserves current cache behavior", {
  resolve <- .resolve_vcov
  local <- new.env(parent = environment(resolve))
  calls <- 0L
  local$vcov.gamlss.longitudinal <- function(object, numderiv, details, method) {
    calls <<- calls + 1L
    expect_false(numderiv)
    expect_true(details)
    expect_identical(method, "analytical")
    list(recomputed = TRUE, inference_epoch = 3L)
  }
  local$.gl_enrich_vcov_contract <- function(value, object) value
  environment(resolve) <- local
  for (discrete in c(FALSE, TRUE)) {
    fit <- chain_epoch_cache_fixture(discrete, version = 2L)
    before <- fit
    expect_true(resolve(fit)$recomputed)
    expect_identical(fit, before)
    fit$vcov_meta$cache_version <- 3L
    expect_identical(resolve(fit), fit$vcov)
  }
  expect_identical(calls, 2L)
})

test_that("fit cache writer records the current epoch without changing estimates", {
  attach <- .gl_attach_fit_vcov
  local <- new.env(parent = environment(attach))
  local$vcov.gamlss.longitudinal <- function(object, numderiv, method, details, progress) {
    .gl_vcov_stamp_epoch(list(vcov = list(overall = diag(2)), method = method,
      method_requested = method, hessian_diagnostics = list(status = "available")))
  }
  environment(attach) <- local
  fit <- chain_epoch_cache_fixture(TRUE, version = 2L)
  fit$par <- c(mu_x = 0.2, mu_z = 0.4)
  before <- fit
  after <- attach(fit, TRUE, FALSE, "analytical", 0)
  expect_identical(after$vcov_meta$cache_version, 3L)
  expect_true(after$vcov_meta$precomputed)
  expect_identical(after$vcov$inference_epoch, 3L)
  expect_identical(after$vcov$analytical_margin_chain_rule_version, 1L)
  expect_identical(after$par, before$par)
  expect_identical(fit, before)
})

test_that("epoch invalidation does not loosen the stored unavailable inference gate", {
  resolve <- .resolve_vcov
  local <- new.env(parent = environment(resolve))
  local$.gl_inference_unavailable <- function(diagnostics) simpleError("stored inference unavailable")
  local$vcov.gamlss.longitudinal <- function(...) list(explicitly_recomputed = TRUE)
  environment(resolve) <- local
  fit <- chain_epoch_cache_fixture(TRUE, version = 2L)
  fit$vcov_meta$inference_status <- "unavailable"
  fit$vcov_meta$hessian_diagnostics <- list(failure_codes = character())
  expect_error(resolve(fit), "stored inference unavailable", fixed = TRUE)
  expect_true(resolve(fit, extra_args = list(inference = list()))$explicitly_recomputed)
})

test_that("fixed assembler stamps used curvature rather than requested curvature", {
  build <- .gl_vcov_build_result
  local <- new.env(parent = environment(build))
  local$.gl_compute_smooth_vcov <- function(...) list(smooth_vcov = list(), smooth_se = list())
  local$.gl_enrich_vcov_contract <- function(value, object) value
  environment(build) <- local
  for (discrete in c(FALSE, TRUE)) {
    fit <- chain_epoch_cache_fixture(discrete)
    for (used in c("analytical", "numderiv", "sandwich_cluster")) {
      for (bread in c("analytical", "numderiv")) {
        result <- build(fit, list(), numeric(), diag(2), c(1, 1), used, "analytical",
                        list(status = "available", bread_method = bread))
        expect_identical(result$inference_epoch, 3L)
        analytical <- used == "analytical" || (used == "sandwich_cluster" && bread == "analytical")
        expect_identical(result$analytical_margin_chain_rule_version, if (analytical) 1L else NULL)
        expect_identical(result$method, used)
        expect_identical(result$method_requested, "analytical")
        matrix_result <- .gl_vcov_format_result(result)
        expect_identical(attr(matrix_result, "gamlss_longitudinal_details")$inference_epoch, 3L)
      }
    }
  }
})

test_that("joint results never enter fixed caches and joint requests always recompute", {
  fit <- chain_epoch_cache_fixture(TRUE)
  expect_false(.can_use_cached_vcov(fit, extra_args = list(coefficients = "joint")))
  fit$vcov$coefficients <- "joint"
  expect_false(.can_use_cached_vcov(fit, method = "analytical"))
  resolve <- .resolve_vcov
  local <- new.env(parent = environment(resolve))
  calls <- 0L
  local$vcov.gamlss.longitudinal <- function(object, numderiv, details, method, coefficients) {
    calls <<- calls + 1L
    expect_identical(coefficients, "joint")
    .gl_vcov_stamp_epoch(list(coefficients = coefficients, method = method,
      vcov = list(overall = diag(2))))
  }
  environment(resolve) <- local
  before <- fit
  first <- resolve(fit, extra_args = list(coefficients = "joint"))
  second <- resolve(fit, extra_args = list(coefficients = "joint"))
  expect_identical(calls, 2L)
  expect_identical(first, second)
  expect_identical(first$inference_epoch, 3L)
  expect_identical(fit, before)
})
