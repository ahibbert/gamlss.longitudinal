make_standard_generic_object <- function(objective = "ordinary", converged = TRUE) {
  structure(list(
    par = c(`mu.(Intercept)` = 1, `theta.(Intercept)` = 0.2),
    df_s = list(mu = list(smooth = 0.5)),
    response = c(1, 2, NA, 4),
    response_subject = c(1L, 1L, 2L, 2L),
    calc_lik_out_end = list(log_lik = c(marginal = -8, copula = -2, joint = -10)),
    convergence = list(converged = converged, stop_reason = if (converged) "tolerance" else "iteration_limit"),
    likelihood_contract = list(
      objective = objective,
      between_segment_assumption = if (objective == "segmented") "independent" else "not_applicable"
    ),
    missingness = list(objective = objective)
  ), class = "gamlss.longitudinal")
}

test_that("standard likelihood and information-criterion generics are scalar", {
  fit <- make_standard_generic_object()
  ll <- stats::logLik(fit)

  expect_s3_class(ll, "logLik")
  expect_length(ll, 1L)
  expect_equal(as.numeric(ll), -10)
  expect_equal(attr(ll, "df"), 2.5)
  expect_identical(attr(ll, "nobs"), 2L)
  expect_identical(attr(ll, "nresponses"), 3L)
  expect_length(stats::AIC(fit), 1L)
  expect_equal(stats::AIC(fit), 25)
  expect_equal(stats::BIC(fit), 20 + log(2) * 2.5)
})

test_that("information criteria reject nonconverged fits", {
  fit <- make_standard_generic_object(converged = FALSE)
  expect_error(stats::AIC(fit), class = "gamlss.longitudinal_nonconvergence_error")
  expect_error(stats::BIC(fit), class = "gamlss.longitudinal_nonconvergence_error")
})

test_that("segmented likelihood criteria remain available under their assumption", {
  fit <- make_standard_generic_object(objective = "segmented")
  expect_equal(stats::AIC(fit), 25)
  expect_identical(attr(stats::logLik(fit), "likelihood_contract")$objective, "segmented")
})

test_that("standard vcov output is a coefficient-aligned matrix", {
  coefficient_names <- c("mu.(Intercept)", "theta.(Intercept)")
  mat <- diag(c(0.2, 0.3))
  dimnames(mat) <- list(coefficient_names, coefficient_names)
  details <- list(
    vcov = list(overall = mat, smooth_vcov = list(), smooth_se = list()),
    se = list(overall = sqrt(diag(mat)), smooth_se = list()),
    method = "analytical",
    method_requested = "analytical",
    hessian_diagnostics = list(status = "available")
  )

  out <- gamlss.longitudinal:::.gl_vcov_format_result(details)
  expect_true(is.matrix(out))
  expect_identical(rownames(out), coefficient_names)
  expect_identical(colnames(out), coefficient_names)
  expect_equal(unclass(out)[seq_along(out)], unclass(mat)[seq_along(mat)])
  expect_identical(out$method, "analytical")
  expect_identical(out$vcov$overall, mat)
})

test_that("model-Hessian inference routes to sandwich for segmented objectives", {
  fit <- make_standard_generic_object(objective = "segmented")
  warning_seen <- FALSE
  method <- withCallingHandlers(
    gamlss.longitudinal:::.gl_require_supported_missingness_inference(fit, "analytical"),
    gamlss_longitudinal_segmented_inference_warning = function(w) {
      warning_seen <<- TRUE
      invokeRestart("muffleWarning")
    }
  )
  expect_true(warning_seen)
  expect_identical(method, "sandwich")
  expect_identical(
    gamlss.longitudinal:::.gl_require_supported_missingness_inference(fit, "sandwich"),
    "sandwich"
  )
})

test_that("legacy unversioned vcov caches are not reused", {
  fit <- make_standard_generic_object()
  fit$vcov <- list(vcov = list(overall = diag(2)), hessian_diagnostics = list(status = "available"))
  fit$vcov_meta <- list(method = "analytical", numderiv = FALSE)
  expect_false(gamlss.longitudinal:::.can_use_cached_vcov(fit, method = "analytical"))

  fit$vcov_meta$cache_version <- gamlss.longitudinal:::.gl_vcov_cache_version()
  expect_true(gamlss.longitudinal:::.can_use_cached_vcov(fit, method = "analytical"))

  fit$vcov$method <- "sandwich_cluster"
  fit$vcov_meta$method_used <- "sandwich_cluster"
  expect_true(gamlss.longitudinal:::.can_use_cached_vcov(fit, method = "sandwich"))
})

test_that("fit finalization caches routed segmented inference for summary reuse", {
  fit <- make_standard_generic_object(objective = "segmented")
  fit$par_s <- list()
  coefficient_names <- names(fit$par)
  covariance <- diag(c(0.2, 0.3))
  dimnames(covariance) <- list(coefficient_names, coefficient_names)
  calls <- 0L

  testthat::local_mocked_bindings(
    vcov.gamlss.longitudinal = function(object, numderiv, method, details, progress, ...) {
      calls <<- calls + 1L
      warning(structure(
        list(message = "routed", call = NULL),
        class = c(
          "gamlss_longitudinal_segmented_inference_warning",
          "warning",
          "condition"
        )
      ))
      list(
        vcov = list(overall = covariance, smooth_vcov = list(), smooth_se = list()),
        se = list(overall = sqrt(diag(covariance)), smooth_se = list()),
        method = "sandwich_cluster",
        method_requested = "analytical",
        hessian_diagnostics = list(status = "available")
      )
    },
    .package = "gamlss.longitudinal"
  )

  expect_warning(
    fit <- gamlss.longitudinal:::.gl_attach_fit_vcov(
      return_list = fit,
      compute_vcov = TRUE,
      vcov_numderiv = FALSE,
      vcov_method = "analytical",
      verbose = 0
    ),
    NA
  )
  expect_equal(calls, 1L)
  expect_true(fit$vcov_meta$precomputed)
  expect_identical(fit$vcov_meta$method, "analytical")
  expect_identical(fit$vcov_meta$method_used, "sandwich_cluster")

  resolved <- gamlss.longitudinal:::.resolve_vcov(fit)
  expect_equal(calls, 1L)
  expect_equal(as.numeric(resolved$vcov$overall), as.numeric(covariance))
  expect_identical(dimnames(resolved$vcov$overall), dimnames(covariance))
  expect_identical(resolved$method, "sandwich_cluster")
})

test_that("stored, summary, and standard criteria agree on a fitted object", {
  dat <- make_fixture_factor_time_interaction(n_subject = 8L)
  fit <- fit_fixture_model(
    dat,
    include_dlcopdpar = FALSE,
    mu_formula = "y ~ age",
    sigma_formula = "~ 1",
    theta_formula = "~ 1",
    max_outer_iter = 3L,
    max_inner_iter = 3L,
    outer_stop_crit = 10,
    inner_stop_crit = 10,
    compute_vcov = FALSE
  )
  expect_true(fit$convergence$converged)

  summary_out <- summary(fit, include_vcov = FALSE)
  expect_equal(as.numeric(logLik(fit)), fit$model_selection["LogLik", "joint"])
  expect_equal(AIC(fit), fit$model_selection["AIC", "joint"])
  expect_equal(BIC(fit), fit$model_selection["BIC", "joint"])
  expect_equal(summary_out$fit$AIC, AIC(fit))
  expect_equal(summary_out$fit$BIC, BIC(fit))
  expect_identical(attr(logLik(fit), "nobs"), length(unique(fit$response_subject)))
  expect_identical(attr(logLik(fit), "nresponses"), nobs(fit))

  fit2 <- fit
  comparison <- stats::AIC(fit, fit2)
  expect_s3_class(comparison, "data.frame")
  expect_equal(comparison$AIC, rep(AIC(fit), 2L))
})
