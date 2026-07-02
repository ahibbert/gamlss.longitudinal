test_that("Hessian assembly helper maps coefficient names to parameter blocks", {
  expect_equal(.hessian_assembly_param_of("mu.age"), "mu")
  expect_equal(.hessian_assembly_param_of("sigma.(Intercept)"), "sigma")
  expect_equal(.hessian_assembly_param_of("theta"), "theta")
})

test_that("Hessian assembly helper recovers fixed intercept and smooth basis columns", {
  X <- matrix(
    c(
      1, 20, 0.1,
      1, 30, 0.2
    ),
    nrow = 2,
    byrow = TRUE
  )
  colnames(X) <- c("(Intercept)", "age", "1")
  mm_x <- list(mu = X)
  param_block <- c(mu = "mu", mu.age = "mu", `mu.s(age).1` = "mu", sigma.age = "sigma")

  expect_equal(.hessian_assembly_x_col("mu", param_block, mm_x), c(1, 1))
  expect_equal(.hessian_assembly_x_col("mu.age", param_block, mm_x), c(20, 30))
  expect_equal(.hessian_assembly_x_col("mu.s(age).1", param_block, mm_x), c(0.1, 0.2))
  expect_null(.hessian_assembly_x_col("sigma.age", param_block, mm_x))
})

test_that("Hessian assembly margin derivative helpers combine margin and copula pieces", {
  margin_d1l <- list(mu = c(1, 2), sigma = c(3, 4))
  margin_d2l <- list(
    mu = list(mu = c(0.1, 0.2), sigma = c(0.3, 0.4)),
    sigma = list(mu = c(0.3, 0.4), sigma = c(0.5, 0.6))
  )
  copula_hess <- list(
    cop_d1l_margin = list(mu = c(10, 20)),
    cop_d2l_margin = list(
      mu = list(mu = c(1, 1), sigma = c(2, 2)),
      sigma = list(mu = c(2, 2), sigma = c(3, 3))
    )
  )
  eta_dr <- list(mu = c(1, 1), sigma = c(1, 1), theta = c(1, 1))

  expect_equal(
    .hessian_assembly_margin_d1_obs("mu", margin_d1l, copula_hess, eta_dr),
    c(11, 22)
  )
  expect_equal(
    .hessian_assembly_margin_d1_obs("sigma", margin_d1l, copula_hess, eta_dr),
    c(3, 4)
  )
  expect_equal(
    .hessian_assembly_margin_d2_obs("mu", "sigma", margin_d2l, copula_hess, eta_dr),
    c(2.3, 2.4)
  )
  expect_equal(
    .hessian_assembly_margin_d2_obs("mu", "theta", margin_d2l, copula_hess, eta_dr),
    c(0, 0)
  )
})

test_that("Hessian assembly margin-copula helper uses pair-indexed theta rows", {
  pair_cache <- list(theta_index_map = c(1L, 2L, NA_integer_, 3L))
  row_id1 <- c(1L, 2L, 3L, 4L)
  row_id2 <- c(2L, 3L, 4L, 1L)
  pair_ok <- c(TRUE, TRUE, TRUE, TRUE)
  id1_ok <- row_id1[pair_ok]
  id2_ok <- row_id2[pair_ok]
  copula_hess <- list(
    cop_d2l_margin_theta_u1 = list(mu = c(0.1, 0.2, 99, 0.3)),
    cop_d2l_margin_theta_u2 = list(mu = c(0.4, 0.5, 99, 0.6)),
    cop_d2l_margin_zeta_u1 = list(mu = c(1.1, 1.2, 99, 1.3)),
    cop_d2l_margin_zeta_u2 = list(mu = c(1.4, 1.5, 99, 1.6))
  )
  xa <- c(10, 20, 30, 40)
  dra <- c(1, 2, 3, 4)
  xb <- c(100, 200, 300)
  drb <- c(0.5, 0.6, 0.7)

  theta_out <- .hessian_assembly_margin_copula_block(
    pa = "mu",
    pb = "theta",
    target = "theta",
    xa = xa,
    xb = xb,
    dra = dra,
    drb = drb,
    margin_d2l = list(mu = list()),
    copula_hess = copula_hess,
    pair_cache = pair_cache,
    row_id1 = row_id1,
    pair_ok = pair_ok,
    id1_ok = id1_ok,
    id2_ok = id2_ok
  )
  zeta_out <- .hessian_assembly_margin_copula_block(
    pa = "zeta",
    pb = "mu",
    target = "zeta",
    xa = xb,
    xb = xa,
    dra = drb,
    drb = dra,
    margin_d2l = list(mu = list()),
    copula_hess = copula_hess,
    pair_cache = pair_cache,
    row_id1 = row_id1,
    pair_ok = pair_ok,
    id1_ok = id1_ok,
    id2_ok = id2_ok
  )

  expect_equal(
    theta_out,
    10 * 1 * 0.1 * 0.5 * 100 + 20 * 2 * 0.4 * 0.5 * 100 +
      20 * 2 * 0.2 * 0.6 * 200 + 30 * 3 * 0.5 * 0.6 * 200 +
      40 * 4 * 0.3 * 0.7 * 300 + 10 * 1 * 0.6 * 0.7 * 300
  )
  expect_equal(
    zeta_out,
    10 * 1 * 1.1 * 0.5 * 100 + 20 * 2 * 1.4 * 0.5 * 100 +
      20 * 2 * 1.2 * 0.6 * 200 + 30 * 3 * 1.5 * 0.6 * 200 +
      40 * 4 * 1.3 * 0.7 * 300 + 10 * 1 * 1.6 * 0.7 * 300
  )
})

test_that("Hessian assembly margin-margin helper combines observation and cross-pair terms", {
  margin_d1l <- list(mu = c(1, 2, 3))
  margin_d2l <- list(mu = list(mu = c(0.1, 0.2, 0.3)))
  copula_hess <- list(
    cop_d1l_margin = list(mu = c(10, 20, 30)),
    cop_d2l_margin = list(mu = list(mu = c(1, 2, 3))),
    cross_pair_contribs = list(mu = list(mu = c(0.5, 0.6)))
  )
  eta_dr <- list(mu = c(2, 3, 4))
  eta_d2 <- list(mu = c(0.01, 0.02, 0.03))
  xa <- c(1, 2, 3)
  xb <- c(4, 5, 6)
  pair_ok <- c(TRUE, TRUE)
  id1_ok <- c(1L, 2L)
  id2_ok <- c(2L, 3L)

  out <- .hessian_assembly_margin_margin_block(
    pa = "mu",
    pb = "mu",
    xa = xa,
    xb = xb,
    dra = eta_dr$mu,
    drb = eta_dr$mu,
    margin_d1l = margin_d1l,
    margin_d2l = margin_d2l,
    copula_hess = copula_hess,
    eta_dr = eta_dr,
    eta_d2 = eta_d2,
    pair_ok = pair_ok,
    id1_ok = id1_ok,
    id2_ok = id2_ok
  )

  obs_terms <- (c(0.1, 0.2, 0.3) + c(1, 2, 3)) * eta_dr$mu * eta_dr$mu +
    (c(1, 2, 3) + c(10, 20, 30)) * eta_d2$mu
  obs_expected <- sum(xa * xb * obs_terms)
  cross_expected <-
    xa[1] * eta_dr$mu[1] * 0.5 * xb[2] * eta_dr$mu[2] +
    xa[2] * eta_dr$mu[2] * 0.5 * xb[1] * eta_dr$mu[1] +
    xa[2] * eta_dr$mu[2] * 0.6 * xb[3] * eta_dr$mu[3] +
    xa[3] * eta_dr$mu[3] * 0.6 * xb[2] * eta_dr$mu[2]

  expect_equal(out, obs_expected + cross_expected)
})
