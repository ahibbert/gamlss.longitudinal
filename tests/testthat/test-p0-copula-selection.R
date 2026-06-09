test_that("select_copula recovers the simulated family from native pseudo-observations", {
  cases <- list(
    N = list(theta = 0.6),
    C = list(theta = 4),
    F = list(theta = 9),
    G = list(theta = 2.4),
    J = list(theta = 3),
    t = list(theta = 0.6, zeta = 4)
  )

  for (family in names(cases)) {
    dat <- simulate_longitudinal_dataset(
      n = 500,
      times = 1:5,
      margin_dist = gamlss.dist::NO(),
      margin_params = list(mu = 0, sigma = 1),
      copula_dist = family,
      copula_params = cases[[family]],
      seed = 100 + match(family, names(cases))
    )

    selected <- select_copula(
      data = dat,
      u_var = "u",
      families = c("N", "C", "F", "G", "J", "t"),
      t_df_grid = c(3, 4, 6, 10, 20)
    )

    expect_s3_class(selected, "copula_selection")
    expect_equal(attr(selected, "selected"), family)
    expect_equal(selected$family[1], family)
    expect_equal(selected$best_fit$family, family)
    expect_equal(best_fit(selected)$family, family)
    expect_equal(best_fit_family(selected), family)
    expect_true(all(is.finite(selected$logLik)))
    expect_true(all(is.finite(selected$AIC)))
  }
})

test_that("select_copula accepts direct pseudo-observation pairs", {
  dat <- simulate_longitudinal_dataset(
    n = 200,
    times = 1:4,
    margin_dist = gamlss.dist::NO(),
    margin_params = list(mu = 0, sigma = 1),
    copula_dist = "C",
    copula_params = list(theta = 3),
    seed = 321
  )
  pairs <- dat[order(dat$subject, dat$time), c("subject", "time", "u")]
  pair_data <- do.call(rbind, lapply(split(pairs, pairs$subject), function(x) {
    data.frame(u1 = x$u[-nrow(x)], u2 = x$u[-1])
  }))

  selected <- select_copula(
    u1 = pair_data$u1,
    u2 = pair_data$u2,
    families = c("N", "C", "F", "G", "J", "t")
  )

  expect_equal(attr(selected, "selected"), "C")
})

test_that("select_copula creates pseudo-observations from a supplied margin", {
  dat <- simulate_longitudinal_dataset(
    n = 80,
    times = 1:4,
    margin_dist = gamlss.dist::NO(),
    margin_params = list(mu = 0, sigma = 1),
    copula_dist = "N",
    copula_params = list(theta = 0.5),
    seed = 456
  )

  selected <- select_copula(
    data = dat,
    response_var = "response",
    margin_dist = gamlss.dist::NO(),
    families = c("N", "C"),
    min_pairs = 10
  )

  expect_s3_class(selected, "copula_selection")
  expect_true(best_fit_family(selected) %in% c("N", "C"))
  expect_equal(attr(selected, "pseudo_observation_source"), "margin_dist")
  expect_null(attr(selected, "margin_selection"))
  expect_true(all(is.finite(selected$AIC)))
})

test_that("select_copula can auto-select the temporary margin with a warning", {
  skip_if_not_installed("gamlss")

  dat <- simulate_longitudinal_dataset(
    n = 30,
    times = 1:4,
    margin_dist = gamlss.dist::NO(),
    margin_params = list(mu = 0, sigma = 1),
    copula_dist = "N",
    copula_params = list(theta = 0.4),
    seed = 789
  )

  captured <- capture_warnings(suppressMessages(
    select_copula(
      data = dat,
      response_var = "response",
      families = "N",
      min_pairs = 10
    )
  ))
  selected <- captured$value

  expect_true(any(grepl("margin_dist", captured$warnings, fixed = TRUE)))
  expect_s3_class(selected, "copula_selection")
  expect_equal(best_fit_family(selected), "N")
  expect_equal(attr(selected, "pseudo_observation_source"), "select_margin")
  expect_s3_class(attr(selected, "margin_selection"), "margin_selection")
})

test_that("select_copula reuses time-intercept margin selections", {
  skip_if_not_installed("gamlss")

  dat <- simulate_longitudinal_dataset(
    n = 35,
    times = 1:4,
    margin_dist = gamlss.dist::NO(),
    margin_params = list(mu = 0, sigma = 1),
    copula_dist = "N",
    copula_params = list(theta = 0.35),
    seed = 790
  )
  dat$response <- dat$response + rep(c(-0.8, -0.2, 0.4, 1.0), times = length(unique(dat$subject)))

  margin_selection <- suppressWarnings(suppressMessages(select_margin(
    dat,
    response_var = "response",
    time_var = "time",
    time_intercepts = TRUE,
    type = "realAll",
    families = "NO",
    trace = FALSE
  )))

  selected <- select_copula(
    data = dat,
    response_var = "response",
    margin_dist = margin_selection,
    families = "N",
    min_pairs = 10
  )

  expect_s3_class(selected, "copula_selection")
  expect_equal(best_fit_family(selected), "N")
  expect_equal(attr(selected, "pseudo_observation_source"), "select_margin")
  expect_s3_class(attr(selected, "margin_selection"), "margin_selection")
  expect_true(isTRUE(attr(attr(selected, "margin_selection"), "time_intercepts")))
  expect_true(all(is.finite(selected$AIC)))
})

test_that("select_copula can screen copulas with factor time-pair intercepts", {
  dat <- simulate_longitudinal_dataset(
    n = 60,
    times = 1:4,
    margin_dist = gamlss.dist::NO(),
    margin_params = list(mu = 0, sigma = 1),
    copula_dist = "N",
    copula_params = list(theta = 0.4),
    seed = 791
  )

  selected <- select_copula(
    data = dat,
    u_var = "u",
    families = "N",
    copula_time_intercepts = TRUE,
    min_pairs = 10
  )

  expect_s3_class(selected, "copula_selection")
  expect_true(isTRUE(attr(selected, "copula_time_intercepts")))
  expect_equal(attr(selected, "copula_time_levels"), c("1->2", "2->3", "3->4"))
  expect_equal(selected$n_copula_time_levels, 3)
  expect_true(is.finite(selected$AIC))

  expect_error(
    select_copula(
      u1 = dat$u[seq_len(20)],
      u2 = dat$u[seq_len(20) + 1L],
      families = "N",
      copula_time_intercepts = TRUE
    ),
    "time information"
  )
})

test_that("select_copula accepts fitted longitudinal objects", {
  dat <- make_fixture_factor_time_interaction(n_subject = 10L)
  fit <- fit_fixture_model(dat, include_dlcopdpar = TRUE)

  selected <- select_copula(
    object = fit,
    families = "N",
    min_pairs = 10
  )

  expect_s3_class(selected, "copula_selection")
  expect_equal(best_fit_family(selected), "N")
  expect_equal(attr(selected, "pseudo_observation_source"), "fitted_object")
})
