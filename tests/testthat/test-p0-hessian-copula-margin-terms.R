test_that("copula Hessian margin accumulators preserve expected nested shapes", {
  out <- .copula_hessian_margin_accumulators(
    margin_pars = c("mu", "sigma"),
    n_obs = 3L,
    n_pairs = 2L,
    has_zeta = TRUE
  )

  expect_named(
    out,
    c(
      "cop_d1l_margin",
      "cop_d2l_margin",
      "cop_d2l_margin_theta_u1",
      "cop_d2l_margin_theta_u2",
      "cop_d2l_margin_zeta_u1",
      "cop_d2l_margin_zeta_u2"
    )
  )
  expect_named(out$cop_d1l_margin, c("mu", "sigma"))
  expect_named(out$cop_d2l_margin$mu, c("mu", "sigma"))
  expect_equal(out$cop_d1l_margin$mu, numeric(3L))
  expect_equal(out$cop_d2l_margin$mu$sigma, numeric(3L))
  expect_equal(out$cop_d2l_margin_theta_u1$sigma, numeric(2L))
  expect_equal(out$cop_d2l_margin_zeta_u2$mu, numeric(2L))

  no_zeta <- .copula_hessian_margin_accumulators(
    margin_pars = "mu",
    n_obs = 2L,
    n_pairs = 1L,
    has_zeta = FALSE
  )
  expect_null(no_zeta$cop_d2l_margin_zeta_u1)
  expect_null(no_zeta$cop_d2l_margin_zeta_u2)
})

test_that("copula Hessian margin pair helpers expose derivative formulas", {
  derivs <- list(
    c_val = c(2, 4),
    dcdu1 = c(3, 5),
    dcdu2 = c(7, 11),
    d2cdu1_2 = c(13, 17),
    d2cdu2_2 = c(19, 23),
    dcdth = c(29, 31),
    d2cdthu1 = c(37, 41),
    d2cdthu2 = c(43, 47),
    zeta = list(
      dcdz = c(53, 59),
      d2cdzu1 = c(61, 67),
      d2cdzu2 = c(71, 73)
    )
  )
  k <- 2L
  cv <- derivs$c_val[k]

  d1 <- .copula_hessian_margin_d1_pair(derivs, k, dFi1 = 0.2, dFi2 = 0.3, cv = cv)
  diag <- .copula_hessian_margin_diag_pair(
    derivs, k, dFi1 = 0.2, dFi2 = 0.3, d2Fi1 = 0.4, d2Fi2 = 0.5, cv = cv
  )
  cross <- .copula_hessian_margin_cross_pair(
    derivs, k, dFi1 = 0.2, dFi2 = 0.3, dFj1 = 0.6, dFj2 = 0.7,
    d2Fij_cross1 = 0.8, d2Fij_cross2 = 0.9, cv = cv
  )
  theta <- .copula_hessian_margin_theta_pair(derivs, k, dFi1 = 0.2, dFi2 = 0.3, cv = cv)
  zeta <- .copula_hessian_margin_zeta_pair(derivs, k, dFi1 = 0.2, dFi2 = 0.3, cv = cv)

  expect_equal(d1$i1, derivs$dcdu1[k] * 0.2 / cv)
  expect_equal(d1$i2, derivs$dcdu2[k] * 0.3 / cv)
  expect_equal(
    diag$i1,
    (derivs$d2cdu1_2[k] * 0.2^2 + derivs$dcdu1[k] * 0.4) / cv -
      (derivs$dcdu1[k] * 0.2 / cv)^2
  )
  expect_equal(
    diag$i2,
    (derivs$d2cdu2_2[k] * 0.3^2 + derivs$dcdu2[k] * 0.5) / cv -
      (derivs$dcdu2[k] * 0.3 / cv)^2
  )
  expect_equal(
    cross$i1,
    (derivs$d2cdu1_2[k] * 0.2 * 0.6 + derivs$dcdu1[k] * 0.8) / cv -
      (derivs$dcdu1[k] * 0.2 / cv) * (derivs$dcdu1[k] * 0.6 / cv)
  )
  expect_equal(
    cross$i2,
    (derivs$d2cdu2_2[k] * 0.3 * 0.7 + derivs$dcdu2[k] * 0.9) / cv -
      (derivs$dcdu2[k] * 0.3 / cv) * (derivs$dcdu2[k] * 0.7 / cv)
  )
  expect_equal(theta$i1, derivs$d2cdthu1[k] * 0.2 / cv - (derivs$dcdu1[k] * 0.2 / cv) * (derivs$dcdth[k] / cv))
  expect_equal(theta$i2, derivs$d2cdthu2[k] * 0.3 / cv - (derivs$dcdu2[k] * 0.3 / cv) * (derivs$dcdth[k] / cv))
  expect_equal(zeta$i1, derivs$zeta$d2cdzu1[k] * 0.2 / cv - (derivs$dcdu1[k] * 0.2 / cv) * (derivs$zeta$dcdz[k] / cv))
  expect_equal(zeta$i2, derivs$zeta$d2cdzu2[k] * 0.3 / cv - (derivs$dcdu2[k] * 0.3 / cv) * (derivs$zeta$dcdz[k] / cv))

  expect_equal(.copula_hessian_finite_or_zero(c(1, NaN, Inf, -Inf)), c(1, 0, 0, 0))
})
