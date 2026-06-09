test_that("select_joint_distribution ranks joint fits and exposes best-fit helpers", {
  skip_if_not_installed("gamlss")

  dat <- simulate_longitudinal_dataset(
    n = 12,
    times = 1:3,
    margin_dist = gamlss.dist::NO(),
    margin_params = list(mu = 0, sigma = 1),
    copula_dist = "N",
    copula_params = list(theta = 0.3),
    seed = 1001
  )

  selected <- select_joint_distribution(
    data = dat,
    response_var = "response",
    time_var = "time",
    subject_var = "subject",
    type = "realAll",
    margin_families = "NO",
    copula_families = c("N", "C"),
    progress = FALSE,
    fit_args = list(
      max_outer_iter = 3,
      max_inner_iter = 3,
      outer_stop_crit = 1,
      inner_stop_crit = 1,
      backtracking_max_halves = 0
    )
  )

  expect_s3_class(selected, "joint_distribution_selection")
  expect_equal(sort(unique(selected$margin_family)), "NO")
  expect_equal(sort(unique(selected$copula_family)), c("C", "N"))
  expect_true(any(is.finite(selected$AIC)))
  expect_true(all(is.na(selected$error[is.finite(selected$AIC)])))
  expect_equal(selected$rank[is.finite(selected$AIC)], seq_len(sum(is.finite(selected$AIC))))
  expect_equal(attr(selected, "selected"), paste(selected$margin_family[[1L]], selected$copula_family[[1L]], sep = "+"))

  best <- best_fit(selected)
  best_family <- best_fit_family(selected)
  expect_equal(best$margin_family_name, "NO")
  expect_equal(best$margin_family$family[[1L]], "NO")
  expect_true(best$copula_family %in% c("N", "C"))
  expect_equal(best_family$margin_dist$family[[1L]], "NO")
  expect_equal(best_family$copula_dist, best$copula_family)
  expect_equal(selected$best_fit$copula_family, best$copula_family)
})

test_that("select_joint_distribution can screen with time-intercept margins", {
  skip_if_not_installed("gamlss")

  dat <- simulate_longitudinal_dataset(
    n = 10,
    times = 1:3,
    margin_dist = gamlss.dist::NO(),
    margin_params = list(mu = 0, sigma = 1),
    copula_dist = "N",
    copula_params = list(theta = 0.25),
    seed = 1002
  )
  dat$response <- dat$response + rep(c(-0.4, 0.1, 0.5), times = length(unique(dat$subject)))

  selected <- select_joint_distribution(
    data = dat,
    response_var = "response",
    time_var = "time",
    subject_var = "subject",
    type = "realAll",
    margin_families = "NO",
    copula_families = "N",
    time_intercepts = TRUE,
    progress = FALSE,
    fit_args = list(
      max_outer_iter = 2,
      max_inner_iter = 2,
      outer_stop_crit = 1,
      inner_stop_crit = 1,
      backtracking_max_halves = 0
    )
  )

  expect_s3_class(selected, "joint_distribution_selection")
  expect_true(isTRUE(attr(selected, "time_intercepts")))
  expect_equal(attr(selected, "time_var"), "time")
  expect_true(isTRUE(attr(attr(selected, "margin_selection"), "time_intercepts")))
  expect_equal(attr(attr(selected, "margin_selection"), "time_var"), "time")
  expect_equal(selected$margin_family, "NO")
  expect_equal(selected$copula_family, "N")
})

test_that("select_joint_distribution can screen with factor copula time intercepts", {
  skip_if_not_installed("gamlss")

  dat <- simulate_longitudinal_dataset(
    n = 10,
    times = 1:3,
    margin_dist = gamlss.dist::NO(),
    margin_params = list(mu = 0, sigma = 1),
    copula_dist = "N",
    copula_params = list(theta = 0.25),
    seed = 1003
  )

  selected <- select_joint_distribution(
    data = dat,
    response_var = "response",
    time_var = "time",
    subject_var = "subject",
    type = "realAll",
    margin_families = "NO",
    copula_families = "N",
    copula_time_intercepts = TRUE,
    progress = FALSE,
    fit_args = list(
      max_outer_iter = 2,
      max_inner_iter = 2,
      outer_stop_crit = 1,
      inner_stop_crit = 1,
      backtracking_max_halves = 0
    )
  )

  expect_s3_class(selected, "joint_distribution_selection")
  expect_true(isTRUE(attr(selected, "copula_time_intercepts")))
  expect_equal(attr(selected, "copula_time_var"), "time")
  expect_equal(selected$margin_family, "NO")
  expect_equal(selected$copula_family, "N")
})

test_that("select_joint_distribution retains failures and includes t candidates", {
  skip_if_not_installed("gamlss")

  dat <- simulate_longitudinal_dataset(
    n = 8,
    times = 1:3,
    margin_dist = gamlss.dist::NO(),
    margin_params = list(mu = 0, sigma = 1),
    copula_dist = "N",
    copula_params = list(theta = 0.2),
    seed = 1004
  )

  messages <- character(0)
  selected <- withCallingHandlers(
    select_joint_distribution(
      data = dat,
      response_var = "response",
      time_var = "time",
      subject_var = "subject",
      type = "realAll",
      margin_families = "NO",
      copula_families = c("N", "t"),
      progress = TRUE,
      fit_args = list(max_elapsed_sec = 1e-9)
    ),
    message = function(m) {
      messages <<- c(messages, conditionMessage(m))
      invokeRestart("muffleMessage")
    }
  )

  expect_true(any(grepl("Joint distribution screen", messages, fixed = TRUE)))
  expect_s3_class(selected, "joint_distribution_selection")
  expect_equal(sort(unique(selected$copula_family)), c("N", "t"))
  expect_true(all(!is.na(selected$error)))
  expect_true(all(!is.finite(selected$AIC)))
  expect_true(is.na(attr(selected, "selected")))
})
