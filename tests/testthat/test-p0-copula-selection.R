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
