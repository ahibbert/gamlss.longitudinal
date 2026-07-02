make_summary_helper_object <- function() {
  structure(
    list(
      response = c(1, 2, NA),
      response_subject = c("a", "a", "b"),
      response_margin = c(1, 2, 1),
      margin_dist = gamlss.dist::NO(),
      copula_dist = "N",
      par = c(
        theta.intercept = 0.2,
        mu.intercept = 1,
        sigma.intercept = -0.1
      ),
      par_s = list(
        mu = list(`s(age)` = c(0.1, 0.2)),
        theta = list()
      ),
      df_s = list(
        mu = list(`s(age)` = 1.5),
        theta = list()
      ),
      calc_lik_out_end = list(
        log_lik = c(marginal = -10, copula = -2, joint = -12)
      ),
      vcov_meta = list()
    ),
    class = "gamlss.longitudinal"
  )
}

test_that("summary helpers build model metadata and criteria", {
  object <- make_summary_helper_object()

  model_info <- gamlss.longitudinal:::.gl_summary_model_info(object)
  fit_criteria <- gamlss.longitudinal:::.gl_summary_fit_criteria(object, n_obs = model_info$n_obs)

  expect_equal(model_info$n_obs, 3L)
  expect_equal(model_info$n_subjects, 2L)
  expect_equal(model_info$n_timepoints, 2L)
  expect_equal(model_info$n_fixed, 3L)
  expect_equal(model_info$n_smooth_terms, 1L)
  expect_equal(model_info$edf_smooth, 1.5)

  expect_equal(fit_criteria$loglik_joint, -12)
  expect_equal(fit_criteria$edf_vec, c(marginal = 3.5, copula = 1, joint = 4.5))
  expect_equal(as.numeric(fit_criteria$aic_vec["joint"]), 33)
  expect_true(all(c("LogLik", "AIC", "BIC", "EDF") %in% rownames(fit_criteria$model_selection)))
})

test_that("summary helpers build coefficient and smooth-term tables", {
  object <- make_summary_helper_object()
  V <- diag(c(0.04, 0.09, 0.16))
  dimnames(V) <- list(names(object$par), names(object$par))
  vcov_out <- list(vcov = list(overall = V), se = NULL)

  coef_tbl <- gamlss.longitudinal:::.gl_summary_coefficient_table(object, vcov_out = vcov_out)
  smooth_tbl <- gamlss.longitudinal:::.gl_summary_smooth_terms(object)

  expect_equal(coef_tbl$term, c("mu.intercept", "sigma.intercept", "theta.intercept"))
  expect_equal(coef_tbl$parameter, c("mu", "sigma", "theta"))
  expect_equal(coef_tbl$std_error, c(0.3, 0.4, 0.2))
  expect_true(all(is.finite(coef_tbl$p_value)))
  expect_false(any(c(".original_order", ".param_rank") %in% names(coef_tbl)))

  expect_equal(nrow(smooth_tbl), 1L)
  expect_equal(smooth_tbl$parameter, "mu")
  expect_equal(smooth_tbl$smooth_term, "s(age)")
  expect_equal(smooth_tbl$edf, 1.5)
})

test_that("summary object builder preserves summary return contract", {
  object <- make_summary_helper_object()
  model_info <- gamlss.longitudinal:::.gl_summary_model_info(object)
  fit_criteria <- gamlss.longitudinal:::.gl_summary_fit_criteria(object, n_obs = model_info$n_obs)
  coef_tbl <- gamlss.longitudinal:::.gl_summary_coefficient_table(object, vcov_out = NULL)
  smooth_tbl <- gamlss.longitudinal:::.gl_summary_smooth_terms(object)
  vcov_out <- list(
    method = "analytical",
    method_requested = "analytical",
    hessian_diagnostics = list(status = "ok")
  )

  out <- gamlss.longitudinal:::.gl_build_summary_object(
    object = object,
    model_info = model_info,
    fit_criteria = fit_criteria,
    coef_tbl = coef_tbl,
    smooth_terms = smooth_tbl,
    vcov_out = vcov_out,
    include_vcov = TRUE,
    numderiv = FALSE,
    ci_level = 0.9
  )

  expect_s3_class(out, "summary.gamlss.longitudinal")
  expect_identical(out$model, model_info)
  expect_equal(out$fit$logLik, -12)
  expect_equal(out$fit$AIC, as.numeric(fit_criteria$aic_vec["joint"]))
  expect_equal(out$fit$BIC, as.numeric(fit_criteria$bic_vec["joint"]))
  expect_equal(out$fit$ci_level, 0.9)
  expect_true(out$fit$vcov_included)
  expect_false(out$fit$vcov_numderiv)
  expect_equal(out$fit$vcov_method, "analytical")
  expect_equal(out$fit$vcov_method_requested, "analytical")
  expect_equal(out$fit$hessian_diagnostics, list(status = "ok"))
  expect_identical(out$fit$model_selection, fit_criteria$model_selection)
  expect_identical(out$smooth_terms, smooth_tbl)
  expect_identical(out$coefficients, coef_tbl)
  expect_identical(out$vcov, vcov_out)
})

test_that("summary coefficient display preserves print formatting rules", {
  coef_tbl <- data.frame(
    term = c("mu.intercept", "theta.zeta", "sigma.scale"),
    estimate = c(1.23456, NA_real_, -0.12345),
    std_error = c(0.33333, 0.44444, NA_real_),
    p_value = c(0.000001, NA_real_, 0.04994),
    signif = c("***", NA_character_, "*"),
    parameter = c("mu", "theta", "sigma"),
    stringsAsFactors = FALSE
  )

  out <- gamlss.longitudinal:::.gl_summary_coefficient_display(coef_tbl, digits = 3)

  expect_equal(out$estimate, c("1.235", "NA", "-0.123"))
  expect_equal(out$std_error, c("0.333", "0.444", "NA"))
  expect_equal(out$p_value, c("<0.0001", "NA", "0.0499"))
  expect_equal(out$signif, c("***", "", "*"))
  expect_equal(out$parameter, c("mu", "theta", "sigma"))
})

test_that("summary coefficient printer groups parameters in model order", {
  coef_disp <- data.frame(
    term = c("theta.zeta", "mu.intercept", "sigma.scale"),
    estimate = c("0.200", "1.000", "-0.100"),
    std_error = c("0.100", "0.200", "0.300"),
    p_value = c("0.0500", "0.0100", "0.2000"),
    signif = c(".", "*", ""),
    parameter = c("theta", "mu", "sigma"),
    stringsAsFactors = FALSE
  )

  out <- capture.output(
    gamlss.longitudinal:::.gl_print_summary_coefficient_blocks(coef_disp)
  )

  param_headers <- grep("^  \\[", out, value = TRUE)
  expect_equal(param_headers, c("  [mu]", "  [sigma]", "  [theta]"))
  expect_true(any(grepl("mu.intercept", out, fixed = TRUE)))
  expect_true(any(grepl("sigma.scale", out, fixed = TRUE)))
  expect_true(any(grepl("theta.zeta", out, fixed = TRUE)))
  returned <- NULL
  invisible(
    capture.output(
      returned <- gamlss.longitudinal:::.gl_print_summary_coefficient_blocks(coef_disp)
    )
  )
  expect_identical(returned, coef_disp)
})

test_that("summary smooth-term section handles populated and empty tables", {
  smooth_terms <- data.frame(
    parameter = "mu",
    smooth_term = "s(age)",
    edf = 1.23456,
    n_coef = 2L,
    stringsAsFactors = FALSE
  )

  returned <- NULL
  populated <- capture.output(
    returned <- gamlss.longitudinal:::.gl_print_summary_smooth_terms(smooth_terms, digits = 2)
  )
  empty <- capture.output(
    gamlss.longitudinal:::.gl_print_summary_smooth_terms(NULL, digits = 2)
  )

  expect_identical(returned, smooth_terms)
  expect_true(any(grepl("Smooth terms:", populated, fixed = TRUE)))
  expect_true(any(grepl("s(age)", populated, fixed = TRUE)))
  expect_true(any(grepl("1.23", populated, fixed = TRUE)))
  expect_true(any(grepl("Use plot(object)", populated, fixed = TRUE)))
  expect_true(any(grepl("None", empty, fixed = TRUE)))
})

test_that("summary model-selection section prints stored and fallback criteria", {
  fit <- list(
    model_selection = matrix(
      c(-10, 25, 30, 2),
      nrow = 1,
      dimnames = list("joint", c("LogLik", "AIC", "BIC", "EDF"))
    ),
    logLik = -10,
    AIC = 25,
    BIC = 30
  )

  returned <- NULL
  stored <- capture.output(
    returned <- gamlss.longitudinal:::.gl_print_summary_model_selection(fit, digits = 2)
  )
  fallback <- capture.output(
    gamlss.longitudinal:::.gl_print_summary_model_selection(
      list(model_selection = NULL, logLik = -10, AIC = 25, BIC = 30),
      digits = 2
    )
  )

  expect_identical(returned, fit)
  expect_true(any(grepl("Model Selection Criteria:", stored, fixed = TRUE)))
  expect_true(any(grepl("LogLik", stored, fixed = TRUE)))
  expect_true(any(grepl("logLik", fallback, fixed = TRUE)))
  expect_true(any(grepl("BIC", fallback, fixed = TRUE)))
})

test_that("summary method delegates to summary helpers", {
  object <- make_summary_helper_object()

  out <- summary(object, include_vcov = FALSE)

  expect_s3_class(out, "summary.gamlss.longitudinal")
  expect_equal(out$model$n_obs, 3L)
  expect_equal(out$fit$logLik, -12)
  expect_equal(out$coefficients$term, c("mu.intercept", "sigma.intercept", "theta.intercept"))
  expect_equal(nrow(out$smooth_terms), 1L)
})

test_that("summary print method uses grouped coefficient display", {
  object <- make_summary_helper_object()
  out <- summary(object, include_vcov = FALSE)

  printed <- capture.output(print(out))

  expect_true(any(grepl("GAMLSS Longitudinal Model Summary", printed, fixed = TRUE)))
  expect_true(any(grepl("Fixed coefficients:", printed, fixed = TRUE)))
  expect_true(any(grepl("  [mu]", printed, fixed = TRUE)))
  expect_true(any(grepl("  [sigma]", printed, fixed = TRUE)))
  expect_true(any(grepl("  [theta]", printed, fixed = TRUE)))
  expect_true(any(grepl("Model Selection Criteria:", printed, fixed = TRUE)))
})
