test_that("margin CDF setup resolves p-function arguments for finite differences", {
  eta_inv <- list(mu = c(1, 2), sigma = c(0.5, 0.8), unused = c(9, 9))

  out <- .hessian_margin_cdf_setup(
    eta_inv = eta_inv,
    margin_dist = gamlss.dist::NO(),
    response = c(0, 1)
  )

  expect_equal(out$margin_pars, c("mu", "sigma"))
  expect_true(is.function(out$pfun))
  expect_named(out$args_base, c("q", "mu", "sigma"))
  expect_equal(out$args_base$q, c(0, 1))
  expect_equal(out$args_base$mu, c(1, 2))
  expect_equal(out$args_base$sigma, c(0.5, 0.8))
})

test_that("margin CDF finite-difference helpers preserve expected nested shapes", {
  eta_inv <- list(mu = c(0, 1), sigma = c(1, 1.5))
  response <- c(-0.2, 0.4)
  mm <- list(x = list(mu = matrix(1, 2, 1), sigma = matrix(1, 2, 1)))

  d1 <- .calc_dFdpar(eta_inv, mm, gamlss.dist::NO(), response)
  d2_diag <- .calc_d2Fdpar2(eta_inv, mm, gamlss.dist::NO(), response)
  d2_cross <- .calc_d2Fdpar_cross(eta_inv, mm, gamlss.dist::NO(), response)

  expect_named(d1, c("mu", "sigma"))
  expect_named(d2_diag, c("mu", "sigma"))
  expect_named(d2_cross, c("mu", "sigma"))
  expect_equal(length(d1$mu), length(response))
  expect_equal(length(d2_diag$sigma), length(response))
  expect_equal(length(d2_cross$mu$sigma), length(response))
  expect_equal(d2_cross$mu$sigma, d2_cross$sigma$mu)
  expect_true(all(is.finite(unlist(d1))))
  expect_true(all(is.finite(unlist(d2_diag))))
  expect_true(all(is.finite(unlist(d2_cross))))
})

test_that("margin density setup resolves d-function arguments for finite differences", {
  eta_inv <- list(mu = c(1, 2), sigma = c(0.5, 0.8))
  mm <- list(x = list(mu = matrix(1, 2, 1), sigma = matrix(1, 2, 1)))

  out <- .hessian_margin_density_setup(
    eta_inv = eta_inv,
    mm = mm,
    margin_dist = gamlss.dist::NO(),
    response = c(0, 1)
  )

  expect_equal(out$margin_pars, c("mu", "sigma"))
  expect_true(is.function(out$dfun))
  expect_named(out$base_args, c("x", "log", "mu", "sigma"))
  expect_true(out$base_args$log)
})

test_that("margin log-density finite-difference Hessian returns symmetric nested blocks", {
  eta_inv <- list(mu = c(0, 1), sigma = c(1, 1.5))
  response <- c(-0.2, 0.4)
  mm <- list(x = list(mu = matrix(1, 2, 1), sigma = matrix(1, 2, 1)))

  out <- .calc_margin_d2l_fd(eta_inv, mm, gamlss.dist::NO(), response)

  expect_named(out, c("mu", "sigma"))
  expect_named(out$mu, c("mu", "sigma"))
  expect_equal(length(out$mu$mu), length(response))
  expect_equal(length(out$sigma$sigma), length(response))
  expect_equal(out$mu$sigma, out$sigma$mu)
  expect_true(all(is.finite(unlist(out))))
})
