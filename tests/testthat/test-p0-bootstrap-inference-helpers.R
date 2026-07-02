test_that("bootstrap inference argument helper validates object R and fit args", {
  object <- structure(list(), class = "gamlss.longitudinal")

  out <- .gl_validate_bootstrap_args(object, R = 2.9, fit_args = list(max_outer_iter = 1))

  expect_identical(out$R, 2L)
  expect_equal(out$fit_args, list(max_outer_iter = 1))
  expect_error(.gl_validate_bootstrap_args(list(), R = 1, fit_args = list()), "fitted 'gamlss.longitudinal'")
  expect_error(.gl_validate_bootstrap_args(object, R = 0, fit_args = list()), "positive integer")
  expect_error(.gl_validate_bootstrap_args(object, R = 1, fit_args = list(1)), "named list")
})

test_that("bootstrap simulation argument helper accepts only fitted copula simulation type", {
  expect_equal(
    .gl_normalize_bootstrap_simulation_args(list(simulation_type = "copula", nsim_seed = 1)),
    list(nsim_seed = 1)
  )
  expect_equal(.gl_normalize_bootstrap_simulation_args(list(foo = "bar")), list(foo = "bar"))
  expect_error(
    .gl_normalize_bootstrap_simulation_args(list(simulation_type = "margin")),
    "no longer supported"
  )
})

test_that("bootstrap seed helper preserves caller RNG state", {
  old_seed <- if (exists(".Random.seed", envir = .GlobalEnv, inherits = FALSE)) {
    get(".Random.seed", envir = .GlobalEnv, inherits = FALSE)
  } else {
    NULL
  }
  on.exit({
    if (is.null(old_seed)) {
      if (exists(".Random.seed", envir = .GlobalEnv, inherits = FALSE)) {
        rm(".Random.seed", envir = .GlobalEnv)
      }
    } else {
      assign(".Random.seed", old_seed, envir = .GlobalEnv)
    }
  }, add = TRUE)

  set.seed(1001)
  before <- get(".Random.seed", envir = .GlobalEnv, inherits = FALSE)
  expect_null(.gl_prepare_bootstrap_seed(NULL))
  expect_identical(get(".Random.seed", envir = .GlobalEnv, inherits = FALSE), before)

  restore <- .gl_prepare_bootstrap_seed(42)
  expect_type(restore, "closure")
  expect_false(identical(get(".Random.seed", envir = .GlobalEnv, inherits = FALSE), before))
  stats::runif(1)
  restore()
  expect_identical(get(".Random.seed", envir = .GlobalEnv, inherits = FALSE), before)
})

test_that("cluster bootstrap dataset resamples whole subjects", {
  dat <- data.frame(
    subject = rep(1:3, each = 2),
    time = rep(1:2, times = 3),
    response = rnorm(6)
  )
  object <- structure(list(dataset = dat), class = "gamlss.longitudinal")

  set.seed(1)
  out <- .gl_cluster_bootstrap_dataset(object)

  expect_equal(nrow(out), nrow(dat))
  expect_equal(length(unique(out$subject)), 3L)
  expect_true(all(table(out$subject) == 2L))
})

test_that("bootstrap summary helper computes finite replicate summaries", {
  estimates <- c(mu.age = 1, sigma.age = 2)
  boot_coef <- matrix(
    c(
      1.1, 2.1,
      0.9, NA,
      NA, NA
    ),
    nrow = 3,
    byrow = TRUE
  )
  colnames(boot_coef) <- names(estimates)

  summary <- .gl_bootstrap_summary(estimates, names(estimates), boot_coef, level = 0.5)

  expect_equal(summary$term, names(estimates))
  expect_equal(summary$estimate, c(1, 2))
  expect_equal(summary$bootstrap_mean, c(1, 2.1))
  expect_equal(summary$reps, c(2, 1))
  expect_true(is.finite(summary$conf.low[[1]]))
  expect_true(is.na(.gl_bootstrap_quantile(c(NA, Inf), 0.5)))
})

test_that("bootstrap result helper assembles classed result metadata", {
  boot_coef <- matrix(c(1, NA, 2, 3), nrow = 2, byrow = TRUE)
  colnames(boot_coef) <- c("a", "b")
  summary <- data.frame(term = c("a", "b"))
  errors <- c(NA_character_, "failed")

  out <- .gl_bootstrap_result(
    summary = summary,
    boot_coef = boot_coef,
    errors = errors,
    R = 2L,
    level = 0.95,
    fits = list(NULL, NULL)
  )

  expect_s3_class(out, "gamlss_longitudinal_bootstrap")
  expect_equal(out$successful_replicates, 1L)
  expect_equal(out$failed_replicates, 1L)
  expect_equal(out$simulation_type, "copula")
  expect_s3_class(out$replicates, "data.frame")

  out_cluster <- .gl_bootstrap_result(
    summary = summary,
    boot_coef = boot_coef,
    errors = errors,
    R = 2L,
    level = 0.95,
    simulation_type = "cluster"
  )
  expect_equal(out_cluster$simulation_type, "cluster")
})
