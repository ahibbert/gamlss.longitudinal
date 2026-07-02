test_that("copula diagnostics pair data builds ordered lag pairs and left-row parameters", {
  fit_data <- data.frame(
    subject = c(1, 1, 1, 2, 2),
    time = c(1, 2, 3, 1, 3),
    u = c(0.1, 0.2, 0.3, 0.4, 0.5),
    theta_fit = c(2, 3, 4, 5, 6),
    zeta_fit = c(NA, NA, NA, NA, NA),
    tau_fit = c(0.2, 0.3, 0.4, 0.5, 0.6)
  )

  lag1 <- .copula_v2_pair_data(fit_data, lags = 1)
  lag2 <- .copula_v2_pair_data(fit_data, lags = 2)

  expect_equal(nrow(lag1), 2L)
  expect_equal(lag1$subject, c(1, 1))
  expect_equal(lag1$time_pair, c("T1 vs T2", "T2 vs T3"))
  expect_equal(lag1$theta_pair, c(2, 3))
  expect_equal(lag1$tau_fit, c(0.2, 0.3))
  expect_equal(nrow(lag2), 1L)
  expect_equal(lag2$time_pair, "T1 vs T3")
})

test_that("copula diagnostics pair data validates usable time structure", {
  one_time <- data.frame(
    subject = c(1, 2),
    time = c(1, 1),
    u = c(0.1, 0.2),
    theta_fit = c(1, 1),
    zeta_fit = c(NA, NA),
    tau_fit = c(0.1, 0.1)
  )

  expect_error(.copula_v2_pair_data(one_time), "Need at least two time points")
})

test_that("copula diagnostics grouping supports defaults pair columns and supplied data", {
  pair_data <- data.frame(
    subject = c(1, 2),
    time_left = c("1", "1"),
    time_pair = c("T1 vs T2", "T1 vs T2"),
    lag = c(1, 2)
  )
  object <- list(var_map = list(id = "subject", visit = "time", arm_original = "arm"))
  dat <- data.frame(id = c(1, 2), visit = c(1, 1), arm_original = c("A", "B"))

  default_group <- .copula_v2_attach_group(pair_data, object, by = NULL)
  lag_group <- .copula_v2_attach_group(pair_data, object, by = "lag")
  arm_group <- .copula_v2_attach_group(pair_data, object, by = "arm_original", data = dat)

  expect_equal(as.character(default_group$split_group), c("T1 vs T2", "T1 vs T2"))
  expect_equal(as.character(lag_group$split_group), c("1", "2"))
  expect_equal(as.character(arm_group$split_group), c("A", "B"))
  expect_error(.copula_v2_attach_group(pair_data, object, by = c("a", "b")), "single column name")
})

test_that("copula diagnostics transform maps uniform values to normal scale", {
  pair_data <- data.frame(u1 = c(0, 0.5, 1), u2 = c(0.25, 0.75, 1))

  transformed <- .copula_v2_transform_data(pair_data, transform = "normal")
  unchanged <- .copula_v2_transform_data(pair_data, transform = "uniform")

  expect_true(all(is.finite(transformed$u1)))
  expect_true(all(is.finite(transformed$u2)))
  expect_equal(transformed$u1[2], 0, tolerance = 1e-12)
  expect_equal(unchanged, pair_data)
})
