test_that("Wald contrast helper expands named contrast matrices", {
  estimates <- c(mu.age = 1, sigma.age = 2, theta.int = 3)
  L <- matrix(2, nrow = 1, dimnames = list("age effect", "mu.age"))

  out <- gamlss.longitudinal:::.gl_wald_contrast_matrix(
    L = L,
    terms = NULL,
    estimates = estimates,
    object = structure(list(), class = "gamlss.longitudinal"),
    joint = FALSE
  )

  expect_true(out$joint)
  expect_equal(colnames(out$L), names(estimates))
  expect_equal(rownames(out$L), "age effect")
  expect_equal(out$L[1, ], c(mu.age = 2, sigma.age = 0, theta.int = 0))
})

test_that("Wald contrast helper accepts full unnamed contrast matrices", {
  estimates <- c(mu.age = 1, sigma.age = 2, theta.int = 3)
  L <- matrix(c(1, 0, 0, 0, 1, -1), nrow = 2, byrow = TRUE)

  out <- gamlss.longitudinal:::.gl_wald_contrast_matrix(
    L = L,
    terms = NULL,
    estimates = estimates,
    object = structure(list(), class = "gamlss.longitudinal"),
    joint = FALSE
  )

  expect_true(out$joint)
  expect_equal(colnames(out$L), names(estimates))
  expect_equal(rownames(out$L), c("H1", "H2"))
  expect_equal(unname(out$L), L)
})

test_that("Wald contrast helper builds selector contrasts from terms", {
  estimates <- c(mu.age = 1, sigma.age = 2, theta.int = 3)

  out <- gamlss.longitudinal:::.gl_wald_contrast_matrix(
    L = NULL,
    terms = c("mu.age", "theta.int"),
    estimates = estimates,
    object = structure(list(), class = "gamlss.longitudinal"),
    joint = FALSE
  )

  expect_false(out$joint)
  expect_equal(rownames(out$L), c("mu.age", "theta.int"))
  expect_equal(colnames(out$L), names(estimates))
  expect_equal(unname(out$L[, "sigma.age"]), c(0, 0))
  expect_equal(unname(out$L[, c("mu.age", "theta.int")]), diag(2))
})

test_that("Wald contrast helper validates malformed contrasts", {
  estimates <- c(mu.age = 1, sigma.age = 2)
  object <- structure(list(), class = "gamlss.longitudinal")

  expect_error(
    gamlss.longitudinal:::.gl_wald_contrast_matrix(
      L = matrix(numeric(), nrow = 0, ncol = 2),
      terms = NULL,
      estimates = estimates,
      object = object,
      joint = FALSE
    ),
    "at least one contrast row"
  )
  expect_error(
    gamlss.longitudinal:::.gl_wald_contrast_matrix(
      L = matrix(1, nrow = 1, dimnames = list(NULL, "missing")),
      terms = NULL,
      estimates = estimates,
      object = object,
      joint = FALSE
    ),
    "Column names in 'L'"
  )
  expect_error(
    gamlss.longitudinal:::.gl_wald_contrast_matrix(
      L = matrix(1, nrow = 1, ncol = 1),
      terms = NULL,
      estimates = estimates,
      object = object,
      joint = FALSE
    ),
    "one column per fixed coefficient"
  )
})
