test_that("the BCPE Hessian boundary adapts smooth columns without changing inputs", {
  family <- list(family = c("BCPE", "Box-Cox power exponential"), mu.link = "log")
  eta <- list(mu = c(1, NA_real_, 1), sigma = c(1, -1, 2),
    nu = c(0, Inf, 0), tau = c(2, NA_real_, 2), theta = matrix(0, 2, 1))
  columns <- eta
  for (pn in c("mu", "sigma", "nu", "tau")) columns[[pn]] <- matrix(eta[[pn]], ncol = 1L)
  before <- columns
  rows <- .hessian_bcpe_curvature_eta_rows(columns, family, c(.5, NA_real_, 3))
  expect_identical(rows, eta)
  expect_identical(columns, before)
  expect_equal(.hessian_bcpe_log_mu_eta_curvature(rows, family, c(.5, NA_real_, 3)),
    c(-1, 0, -.25), tolerance = 1e-14)
  expect_error(.hessian_bcpe_log_mu_eta_curvature(columns, family, c(.5, NA_real_, 3)),
    "plain scalar or row vector")
  scalar <- eta; scalar$mu <- matrix(1, 1L, 1L)
  expect_identical(.hessian_bcpe_curvature_eta_rows(scalar, family, 1:3)$mu, 1)
  family$mu.link <- "identity"
  expect_identical(.hessian_bcpe_curvature_eta_rows(columns, family, 1:3), columns)
})

test_that("the BCPE Hessian boundary refuses incompatible matrix shapes and classes", {
  family <- list(family = "BCPE", mu.link = "log")
  bad <- list(matrix(1, 3, 2), matrix(1, 1, 3), matrix(1, 2, 1),
    array(1, c(3, 1, 1)), matrix("1", 3, 1),
    structure(matrix(1, 3, 1), class = c("custom", "matrix", "array")))
  for (pn in c("mu", "sigma", "nu", "tau")) for (value in bad) {
    eta <- list(mu = 1, sigma = 1, nu = 0, tau = 2)
    eta[[pn]] <- value
    expect_error(.hessian_bcpe_curvature_eta_rows(eta, family, 1:3),
      "response-sized column matrix")
  }
})

test_that("standalone Hessian loading includes the column boundary dependency", {
  run <- .gl_source_analytical_hessian_helpers
  fake <- new.env(parent = environment(run)); loaded <- character()
  fake$exists <- function(x, ...) x != ".hessian_bcpe_curvature_eta_rows"
  fake$file.exists <- function(paths) rep(TRUE, length(paths))
  fake$source <- function(path, ...) { loaded <<- c(loaded, basename(path)); invisible(NULL) }
  fake$vcov.gamlss.longitudinal <- function(...) NULL
  environment(run) <- fake
  expect_true(run())
  expect_true("hessian-bcpe-mu-curvature.R" %in% loaded)
  expect_lt(match("hessian-bcpe-mu-curvature.R", loaded), match("hessian-analytical.R", loaded))
  loaded <- character(); fake$exists <- function(...) TRUE
  expect_true(run())
  expect_length(loaded, 0L)
})
