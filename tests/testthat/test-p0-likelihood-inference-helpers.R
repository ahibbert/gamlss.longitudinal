make_likelihood_compare_fit <- function(loglik, n = 4L, par = c(mu = 1), df_s = NULL) {
  structure(
    list(
      response = seq_len(n),
      par = par,
      df_s = df_s,
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

test_that("likelihood comparison table computes sorted sequential LR metrics", {
  reduced <- make_likelihood_compare_fit(-10, par = c(mu = 1))
  full <- make_likelihood_compare_fit(-6, par = c(mu = 1, sigma = 2))
  models <- list(full = full, reduced = reduced)

  out <- gamlss.longitudinal:::.gl_likelihood_compare_table(models, sort = TRUE)

  expect_s3_class(out, "gamlss_longitudinal_likelihood_compare")
  expect_equal(out$model, c("reduced", "full"))
  expect_equal(out$df, c(1, 2))
  expect_equal(out$logLik, c(-10, -6))
  expect_equal(out$AIC, c(22, 16))
  expect_equal(out$delta_df, c(NA_real_, 1))
  expect_equal(out$LR_statistic, c(NA_real_, 8))
  expect_equal(out$p_value[2], stats::pchisq(8, df = 1, lower.tail = FALSE))
})

test_that("likelihood comparison table can preserve input order and warns for n mismatch", {
  a <- make_likelihood_compare_fit(-10, n = 4L, par = c(mu = 1))
  b <- make_likelihood_compare_fit(-6, n = 5L, par = c(mu = 1, sigma = 2))
  models <- list(b = b, a = a)

  expect_warning(
    out <- gamlss.longitudinal:::.gl_likelihood_compare_table(models, sort = FALSE),
    "different observation counts"
  )

  expect_equal(out$model, c("b", "a"))
  expect_equal(out$n_obs, c(5L, 4L))
})

test_that("likelihood_compare delegates to helper output", {
  reduced <- make_likelihood_compare_fit(-10, par = c(mu = 1))
  full <- make_likelihood_compare_fit(-6, par = c(mu = 1, sigma = 2))

  out <- likelihood_compare(reduced = reduced, full = full)

  expect_s3_class(out, "gamlss_longitudinal_likelihood_compare")
  expect_equal(out$model, c("reduced", "full"))
  expect_equal(out$LR_statistic[2], 8)
})
