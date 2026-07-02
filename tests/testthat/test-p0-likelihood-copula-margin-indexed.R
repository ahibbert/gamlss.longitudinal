test_that("indexed dlcopdpar accumulator sums endpoint contributions by observation", {
  out <- .calc_dlcopdpar_indexed(
    row_id1 = c(1, 2, 2),
    row_id2 = c(2, 3, 4),
    dcdu1 = c(2, 4, 8),
    dcdu2 = c(3, 6, 9),
    copula_d = c(2, 2, 2),
    F_nd = c(10, 20, 30, 40),
    n_obs = 4,
    pair_complete = c(TRUE, FALSE, TRUE)
  )

  expect_equal(out, c(10, 80 + 30, 0, 180))
})

test_that("indexed dlcopdpar accumulator validates input lengths", {
  expect_error(
    .calc_dlcopdpar_indexed(
      row_id1 = c(1, 2),
      row_id2 = c(2, 3),
      dcdu1 = 1,
      dcdu2 = c(1, 1),
      copula_d = c(1, 1),
      F_nd = c(1, 1, 1),
      n_obs = 3
    ),
    "inconsistent pair lengths"
  )

  expect_error(
    .calc_dlcopdpar_indexed(
      row_id1 = c(1, 2),
      row_id2 = c(2, 3),
      dcdu1 = c(1, 1),
      dcdu2 = c(1, 1),
      copula_d = c(1, 1),
      F_nd = c(1, 1),
      n_obs = 3
    ),
    "F derivative length"
  )

  expect_error(
    .calc_dlcopdpar_indexed(
      row_id1 = c(1, 2),
      row_id2 = c(2, 3),
      dcdu1 = c(1, 1),
      dcdu2 = c(1, 1),
      copula_d = c(1, 1),
      F_nd = c(1, 1, 1),
      n_obs = 3,
      pair_complete = TRUE
    ),
    "pair_complete length"
  )
})
