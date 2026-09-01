test_that("inference registry is complete, unique, and explicit about scope", {
  contracts <- inference_contracts()

  expect_s3_class(contracts, "data.frame")
  expect_named(contracts, c(
    "contract_version", "contract_id", "producer", "method", "estimand",
    "coefficient_blocks", "conditioning", "omitted_uncertainty", "assumptions",
    "failure_states", "validation_status", "approximation"
  ))
  expect_identical(anyDuplicated(contracts$contract_id), 0L)
  expect_true(all(nzchar(unlist(contracts, use.names = FALSE))))
  expect_setequal(contracts$contract_id, c(
    "fixed_hessian_analytical", "fixed_hessian_numerical",
    "fixed_sandwich_cluster", "smooth_penalized_conditional",
    "bootstrap_parametric_fixed", "bootstrap_cluster_fixed", "wald_fixed",
    "likelihood_ratio_nested", "confint_fixed", "prediction_mu_delta",
    "fitted_distribution_plugin", "marginal_effect_mu_delta",
    "fixed_term_pointwise", "smooth_term_pointwise",
    "publication_coefficients", "publication_predictions_point"
  ))

  fixed <- contracts[contracts$contract_id == "fixed_hessian_analytical", ]
  expect_match(fixed$conditioning, "smooth")
  expect_match(fixed$omitted_uncertainty, "fixed-smooth")
  expect_match(fixed$omitted_uncertainty, "smoothing-parameter")

  prediction <- contracts[contracts$contract_id == "prediction_mu_delta", ]
  expect_match(prediction$coefficient_blocks, "mu fixed")
  expect_match(prediction$approximation, "not a prediction interval")
})

test_that("vcov contract records exact fixed and smooth coefficient blocks", {
  object <- list(
    par = c("mu.x" = 1, "sigma.x" = 0, "theta.x" = 0.2),
    par_s = list(mu = list("s(time)" = c(0.1, 0.2)), sigma = list())
  )
  V <- diag(3)
  dimnames(V) <- list(names(object$par), names(object$par))
  vc <- list(
    method = "analytical",
    vcov = list(
      overall = V,
      smooth_vcov = list(mu = list("s(time)" = diag(2))),
      smooth_se = list(mu = list("s(time)" = c(1, 1)))
    ),
    hessian_diagnostics = list(status = "available", failure_codes = character())
  )

  out <- gamlss.longitudinal:::.gl_enrich_vcov_contract(vc, object)

  expect_identical(out$inference_contract$contract_id, "fixed_hessian_analytical")
  expect_identical(out$inference_contract$coefficient_names, names(object$par))
  expect_identical(
    out$smooth_inference_contract$coefficient_names,
    c("mu:s(time)[1]", "mu:s(time)[2]")
  )
  expect_identical(
    attr(out$vcov$overall, "inference_contract")$contract_id,
    "fixed_hessian_analytical"
  )
  expect_identical(
    attr(out$vcov$smooth_vcov$mu[["s(time)"]], "inference_contract")$coefficient_names,
    c("mu:s(time)[1]", "mu:s(time)[2]")
  )

  vc$method <- "sandwich_cluster"
  sand <- gamlss.longitudinal:::.gl_enrich_vcov_contract(vc, object)
  expect_identical(sand$inference_contract$contract_id, "fixed_sandwich_cluster")
})

make_contract_fit <- function() {
  p <- c("mu.x" = 1, "sigma.x" = 0)
  V <- diag(c(0.04, 0.09))
  dimnames(V) <- list(names(p), names(p))
  structure(list(
    par = p,
    par_s = list(mu = list(), sigma = list()),
    vcov = list(
      method = "analytical",
      vcov = list(overall = V),
      se = list(overall = setNames(sqrt(diag(V)), names(p))),
      hessian_diagnostics = list(status = "available", failure_codes = character())
    ),
    vcov_meta = list(
      numderiv = FALSE, method = "analytical", inference_status = "available",
      cache_version = gamlss.longitudinal:::.gl_vcov_cache_version()
    )
  ), class = "gamlss.longitudinal")
}

test_that("confint and Wald consumers preserve estimates and carry target metadata", {
  fit <- make_contract_fit()

  ci <- confint(fit)
  expect_equal(unname(rowMeans(ci)), unname(fit$par))
  expect_identical(attr(ci, "inference_contract")$contract_id, "confint_fixed")
  expect_identical(attr(ci, "inference_contract")$coefficient_names, names(fit$par))
  expect_identical(
    attr(ci, "inference_contract")$covariance_contract$contract_id,
    "fixed_hessian_analytical"
  )

  wt <- wald_test(fit, terms = "mu.x")
  expect_identical(attr(wt, "inference_contract")$contract_id, "wald_fixed")
  expect_identical(attr(wt, "inference_contract")$coefficient_names, "mu.x")
})

make_contract_prediction_object <- function() {
  X <- matrix(1, nrow = 3, ncol = 1, dimnames = list(NULL, "mu.intercept"))
  X_sigma <- matrix(1, nrow = 3, ncol = 1, dimnames = list(NULL, "sigma.intercept"))
  structure(list(
    response = c(1, 2, NA), response_margin = 1:3,
    response_subject = c("a", "a", "b"), margin_dist = gamlss.dist::NO(),
    copula_dist = "N", par = c(mu.intercept = 10, sigma.intercept = 0),
    par_s = list(mu = NA, sigma = NA),
    model_matrix = list(
      x = list(mu = X, sigma = X_sigma),
      s = list(mu = list(), sigma = list())
    )
  ), class = "gamlss.longitudinal")
}

test_that("prediction and distribution consumers distinguish confidence from variation", {
  object <- make_contract_prediction_object()
  out <- gamlss.longitudinal:::.gl_prediction_interval_frame(
    object = object,
    newdata = NULL,
    fit_values = rep(10, 3),
    interval = "confidence",
    level = 0.95,
    vcov_method = "analytical",
    se_fn = function(...) rep(0.1, 3)
  )
  expect_identical(attr(out, "inference_contract")$contract_id, "prediction_mu_delta")
  expect_identical(attr(out, "inference_contract")$coefficient_names, "mu.intercept")

  diag_data <- gamlss.longitudinal:::.gl_fitted_distribution(object, require_response = FALSE)
  pred <- gamlss.longitudinal:::.gl_prediction_frame(object, require_response = FALSE)
  quant <- gamlss.longitudinal:::.gl_prediction_quantile_frame(
    pred, diag_data$params, diag_data$family, probs = c(0.1, 0.9)
  )
  expect_identical(
    attr(quant, "inference_contract")$contract_id,
    "fitted_distribution_plugin"
  )
})

test_that("bootstrap contracts identify selected targets and failed refits", {
  boot_coef <- matrix(c(1, 2, NA, NA), nrow = 2, byrow = TRUE)
  colnames(boot_coef) <- c("mu.x", "theta.x")
  summary <- data.frame(term = colnames(boot_coef), reps = c(1L, 1L))
  errors <- c(NA_character_, "explicit nonconvergence")

  parametric <- gamlss.longitudinal:::.gl_bootstrap_result(
    summary, boot_coef, errors, R = 2L, level = 0.95,
    simulation_type = "copula"
  )
  expect_equal(parametric$successful_replicates, 1L)
  expect_equal(parametric$failed_replicates, 1L)
  expect_identical(
    parametric$inference_contract$contract_id,
    "bootstrap_parametric_fixed"
  )
  expect_identical(parametric$inference_contract$coefficient_names, colnames(boot_coef))
  expect_match(parametric$inference_contract$observed_failures, "nonconvergence")

  cluster <- gamlss.longitudinal:::.gl_bootstrap_result(
    summary, boot_coef, errors, R = 2L, level = 0.95,
    simulation_type = "cluster"
  )
  expect_identical(cluster$inference_contract$contract_id, "bootstrap_cluster_fixed")
})

test_that("likelihood comparison records conditional reference assumptions", {
  make_model <- function(ll, npar) structure(list(
    par = setNames(seq_len(npar), paste0("b", seq_len(npar))),
    df_s = list(), response = 1:20,
    response_subject = rep(1:10, each = 2), response_margin = 1:20,
    margin_dist = list(family = "NO"), copula_dist = "N",
    model_matrix = list(
      x = list(joint = matrix(
        1, nrow = 20, ncol = npar,
        dimnames = list(NULL, paste0("b", seq_len(npar)))
      )),
      s = list()
    ),
    calc_lik_out_end = list(log_lik = c(joint = ll))
  ), class = "gamlss.longitudinal")

  out <- gamlss.longitudinal:::.gl_likelihood_compare_table(
    list(reduced = make_model(-20, 1), full = make_model(-18, 2)),
    sort = FALSE
  )
  expect_identical(
    attr(out, "inference_contract")$contract_id,
    "likelihood_ratio_nested"
  )
  expect_true(is.finite(out$p_value[[2]]))
})

test_that("user-facing inference documentation names material omissions", {
  root <- normalizePath(testthat::test_path("..", ".."), winslash = "/")
  files <- c(
    "R/inference-contracts.R", "R/model-predict.R", "R/model-vcov.R",
    "R/model-effects.R", "vignettes/inference-uncertainty.Rmd"
  )
  text <- paste(unlist(lapply(file.path(root, files), readLines, warn = FALSE)), collapse = "\n")

  expect_match(text, "fixed-smooth")
  expect_match(text, "smoothing-parameter uncertainty")
  expect_match(text, "not (a )?prediction interval")
  expect_match(text, "inference_contracts")
})
