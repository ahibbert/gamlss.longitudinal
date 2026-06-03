test_that("srr extended multi-seed recovery and benchmark checks are opt-in", {
  skip_if_not(identical(Sys.getenv("GAMLSS_LONGITUDINAL_EXTENDED_TESTS"), "true"))
  skip_if_not_installed("gamlss")
  skip_if_not_installed("gamlss.dist")

  grid <- gamlss.longitudinal:::.coverage_make_case_grid(
    families = "NO",
    copulas = "N",
    methods = "rs_separate",
    designs = "intercept"
  )

  results <- do.call(rbind, lapply(c(2026L, 2027L), function(seed) {
    out <- gamlss.longitudinal:::.coverage_run_grid(
      grid,
      n = 24L,
      times = 1:3,
      seed = seed,
      max_outer_iter = 3L,
      max_inner_iter = 3L,
      max_elapsed_sec = 15,
      dependence = "moderate"
    )
    out$seed_id <- seed
    out
  }))

  expect_equal(length(unique(results$seed_id)), 2L)
  expect_true(all(results$success))
  expect_true(all(is.finite(results$joint_loglik)))
  expect_true(all(is.finite(results$benchmark_mean_bias)))
})

test_that("srr extended stress scenarios are opt-in", {
  skip_if_not(identical(Sys.getenv("GAMLSS_LONGITUDINAL_EXTENDED_TESTS"), "true"))
  skip_if_not_installed("gamlss")
  skip_if_not_installed("gamlss.dist")

  grid <- gamlss.longitudinal:::.coverage_make_case_grid(
    families = "NO",
    copulas = "N",
    methods = "rs_separate",
    designs = "intercept"
  )

  results <- do.call(rbind, lapply(c("near_independent", "strong"), function(dependence) {
    gamlss.longitudinal:::.coverage_run_grid(
      grid,
      n = 24L,
      times = 1:3,
      seed = 3030L,
      max_outer_iter = 3L,
      max_inner_iter = 3L,
      max_elapsed_sec = 15,
      dependence = dependence
    )
  }))

  expect_equal(sort(unique(results$dependence)), c("near_independent", "strong"))
  expect_true(all(results$success))
  expect_true(all(is.finite(results$fitted_copula_tau)))
})
