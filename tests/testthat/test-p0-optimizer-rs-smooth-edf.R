test_that("RS smooth EDF refresh replaces stale initialization values", {
  B_mu <- cbind(c(-1, 0, 1), c(1, -2, 1))
  B_sigma <- cbind(c(-0.5, 0, 0.5))
  attr(B_mu, "penalty") <- diag(c(1, 2))
  attr(B_sigma, "penalty") <- matrix(1, 1, 1)

  mm_s <- list(
    mu = list(`s(x)` = B_mu),
    sigma = list(`s(z)` = B_sigma),
    nu = list()
  )
  lambda_s <- list(
    mu = list(`s(x)` = 3),
    sigma = list(`s(z)` = 2),
    nu = list()
  )
  weights <- list(mu = c(1, 2, 1), sigma = c(2, 1, 2), nu = rep(1, 3))
  stale <- list(mu = list(`s(x)` = 0), sigma = list(`s(z)` = 0), nu = list())

  out <- gamlss.longitudinal:::.gl_refresh_rs_smooth_edf(
    df_s = stale,
    mm_s = mm_s,
    lambda_s = lambda_s,
    weights_final = weights
  )

  expected <- function(B, w, S, lambda) {
    information <- crossprod(B, B * w)
    sum(solve(information + lambda * S) * information)
  }
  expect_equal(unname(out$mu$`s(x)`), expected(B_mu, weights$mu, attr(B_mu, "penalty"), 3))
  expect_equal(unname(out$sigma$`s(z)`), expected(B_sigma, weights$sigma, attr(B_sigma, "penalty"), 2))
  expect_identical(out$nu, list())
  expect_gt(unname(out$mu$`s(x)`), 0)
  expect_gt(unname(out$sigma$`s(z)`), 0)
})

test_that("RS smooth EDF validates final-state inputs", {
  B <- matrix(c(-0.5, 0.5), ncol = 1)
  attr(B, "penalty") <- matrix(1, 1, 1)

  expect_error(
    gamlss.longitudinal:::.gl_rs_smooth_term_edf(B, 1, attr(B, "penalty"), 2),
    "one finite working weight"
  )
  expect_error(
    gamlss.longitudinal:::.gl_rs_smooth_term_edf(B, c(1, 1), attr(B, "penalty"), -1),
    "non-negative smoothing parameter"
  )
})
