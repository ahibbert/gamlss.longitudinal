test_that("copula-to-margin first derivative helper accumulates endpoint rows", {
  input <- data.frame(
    row_id1 = c(1, 2),
    row_id2 = c(2, 3),
    time1 = c(1, 2),
    time2 = c(2, 3),
    subject1 = c(1, 1),
    subject2 = c(1, 1),
    dcdu1 = c(2, 4),
    dcdu2 = c(3, 5),
    dldmu.x = c(0, 0),
    dldmu.y = c(0, 0),
    response.x = c(1, 2),
    response.y = c(2, 3),
    margin_d.x = c(1, 1),
    margin_d.y = c(1, 1),
    copula_d = c(2, 4),
    mu.x = c(0, 0),
    mu.y = c(0, 0),
    F_nd.x = c(0.5, 0.1),
    F_nd.y = c(0.25, 0.2)
  )

  out <- .calc_deriv_copula_wrt_margin_d1(input, margin_par = c("mu", "sigma"), par_name = "mu")

  expect_equal(colnames(out), c("mu", "sigma"))
  expect_equal(out[, "mu"], c(0.5, 0.375 + 0.1, 0.25), tolerance = 1e-12)
  expect_equal(out[, "sigma"], c(0, 0, 0))
  expect_equal(
    calc_deriv_copula_wrt_margin(input, c("mu", "sigma"), "mu", calc_d2 = FALSE),
    out
  )
})

test_that("copula-to-margin first derivative helper validates row ids", {
  input <- data.frame(
    row_id1 = NA_real_,
    row_id2 = NA_real_,
    dcdu1 = 1,
    dcdu2 = 1,
    dldmu.x = 0,
    dldmu.y = 0,
    response.x = 1,
    response.y = 2,
    margin_d.x = 1,
    margin_d.y = 1,
    copula_d = 1,
    mu.x = 0,
    mu.y = 0,
    F_nd.x = 1,
    F_nd.y = 1
  )

  expect_error(
    .calc_deriv_copula_wrt_margin_d1(input, margin_par = "mu", par_name = "mu"),
    "Invalid row ids",
    fixed = TRUE
  )
})
