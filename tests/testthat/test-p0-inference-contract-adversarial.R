make_adversarial_contract_fit <- function() {
  p <- c("mu.x" = 1, "sigma.x" = 0)
  V <- diag(c(0.04, 0.09))
  dimnames(V) <- list(names(p), names(p))
  structure(list(
    par = p,
    par_s = list(mu = list(), sigma = list()),
    vcov = list(
      method = "analytical", method_requested = "analytical",
      vcov = list(overall = V),
      se = list(overall = stats::setNames(sqrt(diag(V)), names(p))),
      hessian_diagnostics = list(status = "available", failure_codes = character())
    ),
    vcov_meta = list(numderiv = FALSE, method = "analytical", inference_status = "available")
  ), class = "gamlss.longitudinal")
}

make_adversarial_prediction_object <- function() {
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

predict.fake_inference_effect <- function(object, newdata, type, ...) {
  out <- data.frame(
    subject = seq_len(nrow(newdata)), time = seq_len(nrow(newdata)),
    response = NA_real_, fit = rep(2, nrow(newdata)),
    se.fit = rep(0.2, nrow(newdata))
  )
  attr(out, "inference_contract") <- object$prediction_contract
  out
}
registerS3method(
  "predict", "fake_inference_effect", predict.fake_inference_effect,
  envir = asNamespace("stats")
)

test_that("duplicate joint Wald contrasts fail with a classed unavailable condition", {
  fit <- make_adversarial_contract_fit()
  L <- rbind(first = c(1, 0), duplicate = c(1, 0))
  colnames(L) <- names(fit$par)

  err <- expect_error(
    wald_test(fit, L = L),
    class = "gamlss_longitudinal_inference_unavailable"
  )
  expect_identical(err$diagnostics$status, "unavailable")
  expect_identical(err$diagnostics$failure_codes, "wald_contrast_covariance_singular")
})

test_that("likelihood references remain unverified without nesting evidence", {
  model <- function(ll, p) structure(list(
    response = 1:4,
    par = stats::setNames(seq_len(p), paste0("b", seq_len(p))),
    df_s = list(),
    calc_lik_out_end = list(log_lik = c(joint = ll))
  ), class = "gamlss.longitudinal")

  out <- gamlss.longitudinal:::.gl_likelihood_compare_table(
    list(reduced = model(-10, 1), full = model(-8, 2))
  )
  expect_identical(out$reference_status[[2L]], "unverified")
  expect_true(is.na(out$p_value[[2L]]))
  expect_match(out$reference_failure[[2L]], "unverified")

  err <- expect_error(
    gamlss.longitudinal:::.gl_likelihood_compare_table(
      list(reduced = model(-8, 1), full = model(-10, 2))
    ),
    class = "gamlss_longitudinal_inference_unavailable"
  )
  expect_identical(err$diagnostics$failure_codes, "negative_likelihood_ratio")
})

test_that("prediction metadata propagates analytical-to-numerical fallback", {
  object <- make_adversarial_prediction_object()
  covariance_contract <- gamlss.longitudinal:::.gl_fixed_inference_contract(
    list(
      method_requested = "analytical",
      method = "numderiv_fallback",
      hessian_diagnostics = list(
        status = "available", failure_codes = character(),
        fallback_reason = "analytical_derivative_disagreement"
      )
    ),
    coefficient_names = "mu.intercept"
  )
  se_fn <- function(...) {
    out <- rep(0.1, 3)
    attr(out, "inference_contract") <- covariance_contract
    out
  }
  out <- gamlss.longitudinal:::.gl_prediction_interval_frame(
    object, NULL, rep(10, 3), "confidence", 0.95, "analytical", se_fn = se_fn
  )
  contract <- attr(out, "inference_contract")
  expect_identical(contract$method_requested, "analytical")
  expect_identical(contract$method_used, "numderiv_fallback")
  expect_true(contract$fallback_used)
  expect_identical(
    contract$diagnostics$fallback_reason,
    "analytical_derivative_disagreement"
  )

  expect_error(
    gamlss.longitudinal:::.gl_prediction_interval_frame(
      object, NULL, rep(10, 3), "confidence", 0.95, "analytical",
      se_fn = function(...) rep(NA_real_, 3)
    ),
    class = "gamlss_longitudinal_inference_unavailable"
  )

  unmatched <- object
  V <- matrix(1, 1, 1, dimnames = list("sigma.intercept", "sigma.intercept"))
  unmatched$vcov <- list(
    method = "analytical", method_requested = "analytical",
    vcov = list(overall = V),
    hessian_diagnostics = list(status = "available", failure_codes = character())
  )
  unmatched$vcov_meta <- list(numderiv = FALSE, method = "analytical")
  expect_error(
    gamlss.longitudinal:::.gl_predict_response_se(unmatched),
    class = "gamlss_longitudinal_inference_unavailable"
  )
})

test_that("marginal effects retain the nested prediction covariance contract", {
  covariance_contract <- gamlss.longitudinal:::.gl_fixed_inference_contract(
    list(
      method_requested = "analytical", method = "numderiv_fallback",
      hessian_diagnostics = list(status = "available", failure_codes = character())
    ), coefficient_names = "mu.x"
  )
  prediction_contract <- gamlss.longitudinal:::.gl_inference_contract(
    "prediction_mu_delta", coefficient_names = "mu.x",
    method = "numderiv_fallback", validity_status = "available"
  )
  prediction_contract$covariance_contract <- covariance_contract
  prediction_contract$method_requested <- "analytical"
  prediction_contract$method_used <- "numderiv_fallback"
  prediction_contract$fallback_used <- TRUE
  object <- structure(
    list(par = c(mu.x = 1), prediction_contract = prediction_contract),
    class = c("fake_inference_effect", "gamlss.longitudinal")
  )
  out <- marginal_effects(
    object, newdata = data.frame(x = 1:3), variable = "x",
    values = c(1, 2), se.fit = TRUE
  )
  contract <- attr(out, "inference_contract")
  expect_identical(contract$prediction_contract$contract_id, "prediction_mu_delta")
  expect_identical(contract$covariance_contract$contract_id, "fixed_hessian_numerical")
  expect_identical(contract$method_used, "numderiv_fallback")
  expect_true(contract$fallback_used)
})

test_that("smooth covariance contracts distinguish complete partial and unavailable blocks", {
  object <- list(
    par = c(mu.intercept = 1),
    par_s = list(mu = list(a = c(0, 0), b = c(0, 0), c = c(0, 0)))
  )
  V <- matrix(1, 1, 1, dimnames = list("mu.intercept", "mu.intercept"))
  vc <- list(
    method = "analytical",
    vcov = list(
      overall = V,
      smooth_vcov = list(mu = list(a = diag(2))),
      smooth_se = list(mu = list(b = c(0.1, 0.2)))
    ),
    hessian_diagnostics = list(status = "available", failure_codes = character())
  )
  out <- gamlss.longitudinal:::.gl_enrich_vcov_contract(vc, object)
  status <- out$smooth_inference_contract$block_status
  expect_identical(out$smooth_inference_contract$validity_status, "partial_approximate")
  expect_identical(status$status, c("complete", "diagonal_only", "unavailable"))
})

test_that("smooth plot metadata reports unavailable bands", {
  B <- cbind(a = 1:3, b = c(0, 1, 0))
  object <- structure(list(
    par = c(mu.intercept = 1),
    par_s = list(mu = list("s(time)" = c(0.1, 0.2))),
    model_matrix = list(s = list(mu = list("s(time)" = B))),
    response = 1:3
  ), class = "gamlss.longitudinal")
  out <- suppressWarnings(plot_smooth_terms(
    object,
    vcov_obj = list(vcov = list(smooth_vcov = NULL, smooth_se = NULL)),
    setup_mfrow = FALSE, show_legend = FALSE, even_grid = FALSE
  ))
  expect_identical(attr(out, "inference_contract")$validity_status, "unavailable")
  expect_identical(
    attr(out$mu[["s(time)"]], "inference_contract")$validity_status,
    "unavailable"
  )
  expect_s3_class(out$mu[["s(time)"]]$plot, "ggplot")
})

test_that("sandwich cluster contributions reject invalid included densities", {
  object <- list(
    response = c(1, 2), response_margin = 1:2,
    response_subject = c(1, 1), margin_dist = list(family = "NO"),
    copula_dist = "N", par = c(mu = 0), par_s = list(),
    model_matrix = list(x = list(), s = list())
  )
  testthat::local_mocked_bindings(
    calc_eta = function(...) list(eta_inv = list()),
    build_copula_pair_cache = function(...) list(),
    get_copula_dist = function(...) list(copula_link = "identity"),
    calc_likelihood_minimal = function(...) list(
      margin_d = c(1, NA_real_),
      copula_row_id1 = integer(), copula_row_id2 = integer(),
      pair_complete = logical(), copula_d = numeric(),
      log_lik = c(joint = 0)
    ),
    .package = "gamlss.longitudinal"
  )

  err <- expect_error(
    gamlss.longitudinal:::.gl_cluster_joint_loglik_contributions(object),
    class = "gamlss_longitudinal_inference_unavailable"
  )
  expect_identical(err$diagnostics$failure_codes, "sandwich_invalid_margin_contribution")
})

test_that("bootstrap computational minimum is not labelled inferential adequacy", {
  boot_coef <- matrix(c(1, 2, 1.1, 2.1), nrow = 2, byrow = TRUE)
  colnames(boot_coef) <- c("mu.x", "theta.x")
  summary <- data.frame(term = colnames(boot_coef), reps = c(2L, 2L))
  out <- gamlss.longitudinal:::.gl_bootstrap_result(
    summary, boot_coef, rep(NA_character_, 2), R = 2L, level = 0.95
  )
  expect_identical(out$inference_contract$computational_minimum, 2L)
  expect_identical(
    out$inference_contract$inferential_adequacy,
    "not_assessed_from_replicate_count"
  )
  expect_match(out$inference_contract$adequacy_note, "not evidence")
})
