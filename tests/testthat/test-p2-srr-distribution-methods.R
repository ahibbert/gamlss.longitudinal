test_that("srr copula analytical derivatives agree with finite differences", {
  u1 <- c(0.2, 0.5, 0.8)
  u2 <- c(0.3, 0.6, 0.7)
  rho <- 0.35
  h <- 1e-5

  analytic <- gamlss.longitudinal:::.copula_gaussian_deriv(
    u1,
    u2,
    par = rho,
    deriv = "par",
    log = FALSE
  )
  finite <- (
    gamlss.longitudinal:::.copula_gaussian_pdf(u1, u2, rho + h) -
      gamlss.longitudinal:::.copula_gaussian_pdf(u1, u2, rho - h)
  ) / (2 * h)

  expect_equal(analytic, finite, tolerance = 1e-5)
})

test_that("srr t-copula numerical integration returns stable probabilities", {
  cdf <- gamlss.longitudinal:::.copula_t_cdf(
    u1 = c(0.1, 0.5, 0.9),
    u2 = c(0.2, 0.5, 0.8),
    par = c(-0.2, 0, 0.4),
    par2 = c(4, 8, 12)
  )

  expect_true(all(is.finite(cdf)))
  expect_true(all(cdf > 0 & cdf < 1))
})

test_that("srr count distribution tails use finite p/q/d family operations", {
  skip_if_not_installed("gamlss.dist")

  probs <- c(0.1, 0.5, 0.9, 0.99)
  mu <- rep(3, length(probs))
  q <- gamlss.dist::qPO(probs, mu = mu)
  cdf_at_q <- gamlss.dist::pPO(q, mu = mu)
  cdf_before_q <- gamlss.dist::pPO(pmax(q - 1, 0), mu = mu)
  density <- gamlss.dist::dPO(q, mu = mu)

  expect_true(all(is.finite(q)))
  expect_true(all(is.finite(cdf_at_q)))
  expect_true(all(is.finite(density)))
  expect_true(all(cdf_at_q + .Machine$double.eps^0.5 >= probs))
  expect_true(all(cdf_before_q <= probs | q == 0))
  expect_true(all(density >= 0))
})
