test_that("diagnostic distribution lookup uses gamlss.dist namespace", {
  skip_if_not_installed("gamlss.dist")

  pLOGNO2 <- function(...) NA_real_
  p <- gamlss.longitudinal:::.gl_call_family_fun(
    "p",
    "LOGNO2",
    1,
    list(mu = 1, sigma = 0.5)
  )

  expect_equal(p, gamlss.dist::pLOGNO2(1, mu = 1, sigma = 0.5))
})
