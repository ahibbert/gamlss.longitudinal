make_likelihood_compare_fit <- function(loglik, n = 4L, par = c(mu = 1), df_s = NULL) {
  structure(
    list(
      response = seq_len(n),
      response_subject = rep(seq_len(ceiling(n / 2)), each = 2, length.out = n),
      response_margin = seq_len(n),
      par = par,
      df_s = df_s,
      margin_dist = list(family = "NO"),
      copula_dist = "N",
      model_matrix = list(
        x = list(joint = matrix(1, nrow = n, ncol = length(par),
                                dimnames = list(NULL, names(par)))),
        s = list()
      ),
      calc_lik_out_end = list(log_lik = c(joint = loglik))
    ),
    class = "gamlss.longitudinal"
  )
}

test_that("likelihood comparison input helper validates and labels models", {
  reduced <- make_likelihood_compare_fit(-10)
  full <- make_likelihood_compare_fit(-8)

  out <- gamlss.longitudinal:::.gl_likelihood_compare_models(list(reduced, full))
  expect_equal(names(out), c("model_1", "model_2"))
  expect_s3_class(out[[1]], "gamlss.longitudinal")

  labelled <- gamlss.longitudinal:::.gl_likelihood_compare_models(
    list(reduced, full),
    labels = c("fit_zeta", "fit_full")
  )
  expect_equal(names(labelled), c("fit_zeta", "fit_full"))

  named <- gamlss.longitudinal:::.gl_likelihood_compare_models(list(reduced = reduced, full = full))
  expect_equal(names(named), c("reduced", "full"))

  nested <- gamlss.longitudinal:::.gl_likelihood_compare_models(list(list(a = reduced, b = full)))
  expect_equal(names(nested), c("a", "b"))

  expect_error(
    gamlss.longitudinal:::.gl_likelihood_compare_models(list(reduced)),
    "At least two fitted models"
  )
  expect_error(
    gamlss.longitudinal:::.gl_likelihood_compare_models(list(reduced, list())),
    "All inputs must be fitted"
  )
})

test_that("likelihood comparison captures call labels", {
  labels <- gamlss.longitudinal:::.gl_likelihood_compare_call_labels(
    quote(list(fit_zeta, chosen = fit_full))
  )
  expect_equal(labels, c("fit_zeta", "chosen"))

  inline_list_labels <- gamlss.longitudinal:::.gl_likelihood_compare_call_labels(
    quote(list(list(fit_zeta, fit_full)))
  )
  expect_equal(inline_list_labels, c("fit_zeta", "fit_full"))
})

test_that("likelihood comparison refuses AIC sorting of sequential LR rows", {
  reduced <- make_likelihood_compare_fit(-10, par = c(mu = 1))
  full <- make_likelihood_compare_fit(-6, par = c(mu = 1, sigma = 2))
  models <- list(full = full, reduced = reduced)

  expect_error(
    gamlss.longitudinal:::.gl_likelihood_compare_table(models, sort = TRUE),
    "not supported"
  )
})

test_that("likelihood comparison preserves order and withholds reference for sample mismatch", {
  a <- make_likelihood_compare_fit(-10, n = 4L, par = c(mu = 1))
  b <- make_likelihood_compare_fit(-6, n = 5L, par = c(mu = 1, sigma = 2))
  models <- list(a = a, b = b)

  out <- gamlss.longitudinal:::.gl_likelihood_compare_table(models, sort = FALSE)

  expect_equal(out$model, c("a", "b"))
  expect_equal(out$n_obs, c(4L, 5L))
  expect_identical(out$reference_status[[2L]], "unavailable")
  expect_match(out$reference_failure[[2L]], "observed_sample_mismatch")
  expect_true(is.na(out$p_value[[2L]]))
})

test_that("likelihood_compare delegates to helper output", {
  fit_reduced <- make_likelihood_compare_fit(-10, par = c(mu = 1))
  fit_full <- make_likelihood_compare_fit(-6, par = c(mu = 1, sigma = 2))

  out <- likelihood_compare(fit_reduced, fit_full)

  expect_s3_class(out, "gamlss_longitudinal_likelihood_compare")
  expect_equal(out$model, c("fit_reduced", "fit_full"))
  expect_equal(out$AIC, c(22, 16))
  expect_identical(out$reference_status[[2L]], "reference_available")
  expect_true(is.finite(out$p_value[[2L]]))
})
