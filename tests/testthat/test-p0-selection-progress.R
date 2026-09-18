test_that("select_margin reports stages by default and verbose zero is quiet", {
  skip_if_not_installed("gamlss")

  y <- seq(-1, 1, length.out = 30)
  default_output <- capture.output(invisible(select_margin(
    y,
    type = "realAll",
    families = "NO",
    trace = FALSE
  )))
  quiet_output <- capture.output(invisible(select_margin(
    y,
    type = "realAll",
    families = "NO",
    trace = FALSE,
    verbose = 0
  )))

  expect_match(paste(default_output, collapse = "\n"), "[select_margin] Screening", fixed = TRUE)
  expect_match(paste(default_output, collapse = "\n"), "[select_margin] Complete", fixed = TRUE)
  expect_identical(quiet_output, character())
})

test_that("select_copula reports each candidate by default and verbose zero is quiet", {
  u1 <- seq(0.05, 0.95, length.out = 20)
  u2 <- rev(u1)

  default_output <- capture.output(invisible(select_copula(
    u1 = u1,
    u2 = u2,
    families = c("N", "F")
  )))
  quiet_output <- capture.output(invisible(select_copula(
    u1 = u1,
    u2 = u2,
    families = c("N", "F"),
    verbose = 0
  )))

  expect_match(paste(default_output, collapse = "\n"), "[1/2] Fitting", fixed = TRUE)
  expect_match(paste(default_output, collapse = "\n"), "[2/2] Fitting", fixed = TRUE)
  expect_identical(quiet_output, character())
})

test_that("selection verbosity uses the shared validation contract", {
  expect_error(
    select_margin(1:5, families = "PO", verbose = 4),
    "integer from 0 to 3",
    fixed = TRUE
  )
  expect_error(
    select_copula(u1 = 1:10 / 11, u2 = 1:10 / 11, verbose = NA),
    "integer from 0 to 3",
    fixed = TRUE
  )
})
