test_that("Hessian assembly copula parameter blocks use theta-indexed rows", {
  pair_cache <- list(theta_index_map = c(1L, 2L, 2L, 3L))
  row_id1 <- c(1L, 2L, 3L, 4L)
  pair_ok <- c(TRUE, TRUE, FALSE, TRUE)
  copula_hess <- list(
    cop_d2l_theta = c(1, 2, 99, 3),
    cop_d2l_zeta = c(4, 5, 99, 6),
    cop_d2l_thetazeta = c(7, 8, 99, 9)
  )
  xa <- c(10, 20, 30)
  xb <- c(1, 2, 3)
  dra <- c(1, 10, 100)
  drb <- c(2, 20, 200)

  expect_equal(
    .hessian_assembly_copula_parameter_block(
      "theta", "theta", xa, xb, dra, drb,
      copula_hess, pair_cache, row_id1, pair_ok
    ),
    10 * 1 * 1 + 20 * 2 * 2 + 30 * 3 * 3
  )
  expect_equal(
    .hessian_assembly_copula_parameter_block(
      "zeta", "zeta", xa, xb, dra, drb,
      copula_hess, pair_cache, row_id1, pair_ok
    ),
    10 * 1 * 4 + 20 * 2 * 5 + 30 * 3 * 6
  )
  expect_equal(
    .hessian_assembly_copula_parameter_block(
      "theta", "zeta", xa, xb, dra, drb,
      copula_hess, pair_cache, row_id1, pair_ok
    ),
    10 * 1 * 7 * 1 * 2 + 20 * 2 * 8 * 10 * 20 + 30 * 3 * 9 * 100 * 200
  )
  expect_null(
    .hessian_assembly_copula_parameter_block(
      "mu", "sigma", xa, xb, dra, drb,
      copula_hess, pair_cache, row_id1, pair_ok
    )
  )
})
