test_that("public inference consumers use corrected fresh and cached uncertainty", {
  old <- options(gamlss.longitudinal.copula_backend = "native")
  on.exit(options(old), add = TRUE)
  base <- curvature_model_fixture()
  reference <- curvature_reference_covariance(base)
  z <- stats::qnorm(0.975)
  for (backend in c("native", "vinecopula")) {
    options(gamlss.longitudinal.copula_backend = backend)
    for (method in c("analytical_only", "sandwich")) {
      expected <- if (method == "sandwich") reference$sandwich else reference$model
      se <- sqrt(diag(expected))
      object <- .gl_attach_fit_vcov(base, TRUE, FALSE, method, 0)
      expect_equal(as.numeric(object$vcov$vcov$overall), as.numeric(expected), tolerance = 2e-5)
      expect_true(.can_use_cached_vcov(object, method = method))
      before <- object
      initial_calls <- getOption("curvature.validation.counters")$covariance_calls
      check_summary <- function(fit) {
        result <- summary(fit, vcov_method = method)
        row <- match(names(fit$par), result$coefficients$term)
        expect_equal(result$coefficients$std_error[row], as.numeric(se), tolerance = 2e-5)
        if (method == "sandwich") {
          expect_equal(result$coefficients$model_std_error[row],
            as.numeric(sqrt(diag(reference$model))), tolerance = 2e-5)
          expect_equal(result$coefficients$robust_std_error[row], as.numeric(se), tolerance = 2e-5)
        }
      }
      check_intervals <- function(fit) {
        result <- stats::confint(fit, method = method)
        expect_equal(as.numeric(result[, 1L]), as.numeric(fit$par - z * se), tolerance = 2e-5)
        expect_equal(as.numeric(result[, 2L]), as.numeric(fit$par + z * se), tolerance = 2e-5)
      }
      check_wald <- function(fit) {
        result <- wald_test(fit, terms = "mu.x", method = method)
        expect_equal(result$std_error, unname(se[["mu.x"]]), tolerance = 2e-5)
        expect_equal(result$statistic, unname(fit$par[["mu.x"]] / se[["mu.x"]]), tolerance = 2e-5)
      }
      check_prediction <- function(fit) {
        result <- stats::predict(fit, type = "mu", se.fit = TRUE,
          interval = "confidence", vcov_method = method)
        X <- fit$model_matrix$x$mu
        V <- expected[c("mu.intercept", "mu.x"), c("mu.intercept", "mu.x")]
        predicted_se <- sqrt(rowSums((X %*% V) * X))
        mu <- as.numeric(X %*% fit$par[c("mu.intercept", "mu.x")])
        expect_equal(as.numeric(result$se.fit), predicted_se, tolerance = 2e-5)
        expect_equal(as.numeric(result$conf.low), mu - z * predicted_se, tolerance = 2e-5)
        expect_equal(as.numeric(result$conf.high), mu + z * predicted_se, tolerance = 2e-5)
        expect_identical(attr(result, "inference_contract")$contract_id, "prediction_mu_delta")
        expect_identical(attr(result, "inference_contract")$method_used,
          if (method == "sandwich") "sandwich_cluster" else "analytical")
      }
      check_summary(object); check_intervals(object); check_wald(object); check_prediction(object)
      expect_identical(getOption("curvature.validation.counters")$covariance_calls, initial_calls)
      stale <- object
      stale$vcov$continuous_curvature <- list(version = 1L, backend = backend, copula = "N")
      expect_false(.can_use_cached_vcov(stale, method = method))
      stale$vcov$vcov$overall[,] <- 0
      diag(stale$vcov$vcov$overall) <- 100
      stale$vcov$se$overall[] <- 10
      original_stale <- stale
      check_summary(stale); check_intervals(stale); check_wald(stale); check_prediction(stale)
      expect_identical(stale, original_stale)
      expect_identical(object, before)
      expect_identical(object$vcov$smooth_inference_contract$contract_id, "smooth_penalized_conditional")
      device_path <- file.path(validation_output_dir, paste(backend, method, "fixed-terms.pdf", sep = "-"))
      grDevices::pdf(device_path)
      plotted <- tryCatch(plot_fixed_terms(object, vcov_obj = object$vcov),
        finally = grDevices::dev.off())
      expect_equal(plotted$mu$x$se, abs(plotted$mu$x$x) * se[["mu.x"]], tolerance = 2e-5)
      expect_equal(plotted$mu$x$ci_lower,
        plotted$mu$x$fitted - z * plotted$mu$x$se, tolerance = 2e-5)
      unavailable <- object
      unavailable$vcov_meta$inference_status <- "unavailable"
      unavailable$vcov_meta$hessian_diagnostics <- list(failure_codes = "covariance_nonfinite")
      expect_error(stats::confint(unavailable, method = method),
        class = "gamlss_longitudinal_inference_unavailable")
    }
  }
})
