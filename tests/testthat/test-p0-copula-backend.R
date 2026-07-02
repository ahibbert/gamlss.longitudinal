skip_if_no_vinecopula <- function() {
  if (!requireNamespace("VineCopula", quietly = TRUE)) {
    warning("Skipping VineCopula parity test because VineCopula is not installed.", call. = FALSE)
    testthat::skip("VineCopula is not installed.")
  }
}

test_that("copula backend wrappers preserve VineCopula outputs", {
  skip_if_no_vinecopula()

  old_backend <- getOption("gamlss.longitudinal.copula_backend")
  on.exit(options(gamlss.longitudinal.copula_backend = old_backend), add = TRUE)
  options(gamlss.longitudinal.copula_backend = "vinecopula")

  set.seed(101)
  u1 <- stats::runif(12, 0.05, 0.95)
  u2 <- stats::runif(12, 0.05, 0.95)

  cases <- list(
    N = list(par = 0.35, par2 = 0),
    C = list(par = 1.2, par2 = 0),
    F = list(par = 3, par2 = 0),
    G = list(par = 1.5, par2 = 0),
    J = list(par = 1.6, par2 = 0),
    t = list(par = 0.4, par2 = 5)
  )

  for (family in names(cases)) {
    par <- cases[[family]]$par
    par2 <- cases[[family]]$par2
    family_num <- .copula_family_number(family)

    expect_equal(
      .copula_pdf(u1, u2, family, par, par2),
      VineCopula::BiCopPDF(u1, u2, family_num, par, par2)
    )
    expect_equal(
      .copula_deriv(u1, u2, family, par, par2, deriv = "u1"),
      VineCopula::BiCopDeriv(u1, u2, family_num, par, par2, deriv = "u1")
    )
    expect_equal(
      .copula_deriv(u1, u2, family, par, par2, deriv = "par", log = TRUE),
      VineCopula::BiCopDeriv(u1, u2, family_num, par, par2, deriv = "par", log = TRUE)
    )
    expect_equal(
      .copula_deriv2(u1, u2, family, par, par2, deriv = "par"),
      VineCopula::BiCopDeriv2(u1, u2, family_num, par, par2, deriv = "par")
    )
  }
})

test_that("copula backend accepts package family codes only", {
  expect_equal(.copula_family_number("N"), 1)
  expect_equal(.copula_family_number("C"), 3)
  expect_error(.copula_family_number("Normal"), "Unsupported copula family code")
  expect_error(.copula_family_number(1), "single character code")
})

test_that("native Gaussian backend matches VineCopula for implemented operations", {
  skip_if_no_vinecopula()

  old_backend <- getOption("gamlss.longitudinal.copula_backend")
  on.exit(options(gamlss.longitudinal.copula_backend = old_backend), add = TRUE)

  u1 <- seq(0.1, 0.9, length.out = 7)
  u2 <- rev(u1)
  rho <- 0.42

  options(gamlss.longitudinal.copula_backend = "vinecopula")
  pdf_vine <- .copula_pdf(u1, u2, "N", rho)
  tau_vine <- .copula_par_to_tau("N", rho)
  par_vine <- .copula_tau_to_par("N", tau_vine)
  h_vine <- .copula_hfunc1(u1, u2, "N", rho)
  dcdu1_vine <- .copula_deriv(u1, u2, "N", rho, deriv = "u1")
  dcdu2_vine <- .copula_deriv(u1, u2, "N", rho, deriv = "u2")
  dcdrho_vine <- .copula_deriv(u1, u2, "N", rho, deriv = "par", log = FALSE)
  dlogdrho_vine <- .copula_deriv(u1, u2, "N", rho, deriv = "par", log = TRUE)
  d2cdu1_vine <- .copula_deriv2(u1, u2, "N", rho, deriv = "u1")
  d2cdu2_vine <- .copula_deriv2(u1, u2, "N", rho, deriv = "u2")
  d2cdrho_vine <- .copula_deriv2(u1, u2, "N", rho, deriv = "par")

  options(gamlss.longitudinal.copula_backend = "native")
  expect_equal(.copula_pdf(u1, u2, "N", rho), pdf_vine, tolerance = 1e-12)
  expect_equal(.copula_par_to_tau("N", rho), tau_vine, tolerance = 1e-12)
  expect_equal(.copula_tau_to_par("N", tau_vine), par_vine, tolerance = 1e-12)
  expect_equal(.copula_hfunc1(u1, u2, "N", rho), h_vine, tolerance = 1e-12)
  expect_equal(.copula_deriv(u1, u2, "N", rho, deriv = "u1"), dcdu1_vine, tolerance = 1e-10)
  expect_equal(.copula_deriv(u1, u2, "N", rho, deriv = "u2"), dcdu2_vine, tolerance = 1e-10)
  expect_equal(.copula_deriv(u1, u2, "N", rho, deriv = "par", log = FALSE), dcdrho_vine, tolerance = 1e-10)
  expect_equal(.copula_deriv(u1, u2, "N", rho, deriv = "par", log = TRUE), dlogdrho_vine, tolerance = 1e-10)
  expect_equal(.copula_deriv2(u1, u2, "N", rho, deriv = "u1"), d2cdu1_vine, tolerance = 1e-8)
  expect_equal(.copula_deriv2(u1, u2, "N", rho, deriv = "u2"), d2cdu2_vine, tolerance = 1e-8)
  expect_equal(.copula_deriv2(u1, u2, "N", rho, deriv = "par"), d2cdrho_vine, tolerance = 1e-8)
})

test_that("Gaussian rectangle probabilities match VineCopula corner differences", {
  skip_if_no_vinecopula()

  u1 <- c(0.20, 0.40, 0.40, 0.75, 1.00)
  u2 <- c(0.30, 0.60, 0.60, 0.85, 0.50)
  l1 <- c(0.00, 0.10, 0.10, 0.50, 0.80)
  l2 <- c(0.00, 0.20, 0.20, 0.40, 0.00)
  rho <- c(0.35, 0.35, 0.35, 0.20, 0.35)

  expected <- VineCopula::BiCopCDF(u1, u2, family = 1, par = rho) -
    VineCopula::BiCopCDF(l1, u2, family = 1, par = rho) -
    VineCopula::BiCopCDF(u1, l2, family = 1, par = rho) +
    VineCopula::BiCopCDF(l1, l2, family = 1, par = rho)

  expect_equal(.copula_rectangle_prob(u1, u2, l1, l2, "N", rho), pmax(expected, 1e-300), tolerance = 1e-10)
})

test_that("fast Poisson family shortcuts match gamlss.dist calls", {
  skip_if_not_installed("gamlss.dist")

  args <- list(x = 0:8, q = 0:8, mu = seq(0.8, 4.0, length.out = 9))

  expect_equal(
    gamlss.longitudinal:::.call_fast_count_family("d", "PO", args),
    gamlss.dist::dPO(args$x, mu = args$mu),
    tolerance = 1e-14
  )
  expect_equal(
    gamlss.longitudinal:::.call_fast_count_family("p", "PO", args),
    gamlss.dist::pPO(args$q, mu = args$mu),
    tolerance = 1e-14
  )
})

test_that("native Clayton backend matches VineCopula for implemented operations", {
  skip_if_no_vinecopula()

  old_backend <- getOption("gamlss.longitudinal.copula_backend")
  on.exit(options(gamlss.longitudinal.copula_backend = old_backend), add = TRUE)

  u1 <- seq(0.1, 0.9, length.out = 7)
  u2 <- rev(u1)
  theta <- 1.35

  options(gamlss.longitudinal.copula_backend = "vinecopula")
  pdf_vine <- .copula_pdf(u1, u2, "C", theta)
  cdf_vine <- .copula_cdf(u1, u2, "C", theta)
  tau_vine <- .copula_par_to_tau("C", theta)
  par_vine <- .copula_tau_to_par("C", tau_vine)
  h_vine <- .copula_hfunc1(u1, u2, "C", theta)
  dcdu1_vine <- .copula_deriv(u1, u2, "C", theta, deriv = "u1")
  dcdu2_vine <- .copula_deriv(u1, u2, "C", theta, deriv = "u2")
  dcdtheta_vine <- .copula_deriv(u1, u2, "C", theta, deriv = "par", log = FALSE)
  dlogdtheta_vine <- .copula_deriv(u1, u2, "C", theta, deriv = "par", log = TRUE)
  d2cdu1_vine <- .copula_deriv2(u1, u2, "C", theta, deriv = "u1")
  d2cdu2_vine <- .copula_deriv2(u1, u2, "C", theta, deriv = "u2")
  d2cdtheta_vine <- .copula_deriv2(u1, u2, "C", theta, deriv = "par")

  options(gamlss.longitudinal.copula_backend = "native")
  expect_equal(.copula_pdf(u1, u2, "C", theta), pdf_vine, tolerance = 1e-12)
  expect_equal(.copula_cdf(u1, u2, "C", theta), cdf_vine, tolerance = 1e-12)
  expect_equal(.copula_par_to_tau("C", theta), tau_vine, tolerance = 1e-12)
  expect_equal(.copula_tau_to_par("C", tau_vine), par_vine, tolerance = 1e-12)
  expect_equal(.copula_hfunc1(u1, u2, "C", theta), h_vine, tolerance = 1e-12)
  expect_equal(.copula_deriv(u1, u2, "C", theta, deriv = "u1"), dcdu1_vine, tolerance = 1e-10)
  expect_equal(.copula_deriv(u1, u2, "C", theta, deriv = "u2"), dcdu2_vine, tolerance = 1e-10)
  expect_equal(.copula_deriv(u1, u2, "C", theta, deriv = "par", log = FALSE), dcdtheta_vine, tolerance = 1e-10)
  expect_equal(.copula_deriv(u1, u2, "C", theta, deriv = "par", log = TRUE), dlogdtheta_vine, tolerance = 1e-10)
  expect_equal(.copula_deriv2(u1, u2, "C", theta, deriv = "u1"), d2cdu1_vine, tolerance = 1e-8)
  expect_equal(.copula_deriv2(u1, u2, "C", theta, deriv = "u2"), d2cdu2_vine, tolerance = 1e-8)
  expect_equal(.copula_deriv2(u1, u2, "C", theta, deriv = "par"), d2cdtheta_vine, tolerance = 1e-8)
})

test_that("native non-derivative backend matches VineCopula for remaining families", {
  skip_if_no_vinecopula()

  old_backend <- getOption("gamlss.longitudinal.copula_backend")
  on.exit(options(gamlss.longitudinal.copula_backend = old_backend), add = TRUE)

  u1 <- seq(0.12, 0.88, length.out = 6)
  u2 <- rev(seq(0.18, 0.82, length.out = 6))
  cases <- list(
    F = list(par = 3, par2 = 0, tolerance = 1e-7, tau_tolerance = 3e-3, par_tolerance = 1e-2),
    G = list(par = 1.7, par2 = 0, tolerance = 1e-7, tau_tolerance = 1e-7, par_tolerance = 1e-7),
    J = list(par = 1.8, par2 = 0, tolerance = 1e-7, tau_tolerance = 1e-6, par_tolerance = 1e-5),
    t = list(par = 0.45, par2 = 5, tolerance = 1e-5, tau_tolerance = 1e-7, par_tolerance = 1e-7)
  )

  for (family in names(cases)) {
    par <- cases[[family]]$par
    par2 <- cases[[family]]$par2
    family_num <- .copula_family_number(family)

    options(gamlss.longitudinal.copula_backend = "vinecopula")
    pdf_vine <- .copula_pdf(u1, u2, family, par, par2)
    cdf_vine <- .copula_cdf(u1, u2, family, par, par2)
    h_vine <- .copula_hfunc1(u1, u2, family, par, par2)
    tau_vine <- VineCopula::BiCopPar2Tau(family_num, par, par2)
    par_vine <- VineCopula::BiCopTau2Par(family_num, tau_vine)

    options(gamlss.longitudinal.copula_backend = "native")
    expect_equal(.copula_pdf(u1, u2, family, par, par2), pdf_vine, tolerance = cases[[family]]$tolerance)
    expect_equal(.copula_cdf(u1, u2, family, par, par2), cdf_vine, tolerance = cases[[family]]$tolerance)
    expect_equal(.copula_hfunc1(u1, u2, family, par, par2), h_vine, tolerance = cases[[family]]$tolerance)
    expect_equal(.copula_par_to_tau(family, par, par2), tau_vine, tolerance = cases[[family]]$tau_tolerance)
    expect_equal(.copula_tau_to_par(family, tau_vine), par_vine, tolerance = cases[[family]]$par_tolerance)
  }
})

test_that("native Gumbel first derivatives match VineCopula", {
  skip_if_no_vinecopula()

  old_backend <- getOption("gamlss.longitudinal.copula_backend")
  on.exit(options(gamlss.longitudinal.copula_backend = old_backend), add = TRUE)

  u1 <- seq(0.12, 0.88, length.out = 6)
  u2 <- rev(seq(0.18, 0.82, length.out = 6))
  theta <- 1.7

  options(gamlss.longitudinal.copula_backend = "vinecopula")
  dcdu1_vine <- .copula_deriv(u1, u2, "G", theta, deriv = "u1")
  dcdu2_vine <- .copula_deriv(u1, u2, "G", theta, deriv = "u2")
  dcdtheta_vine <- .copula_deriv(u1, u2, "G", theta, deriv = "par", log = FALSE)
  dlogdtheta_vine <- .copula_deriv(u1, u2, "G", theta, deriv = "par", log = TRUE)

  options(gamlss.longitudinal.copula_backend = "native")
  expect_equal(.copula_deriv(u1, u2, "G", theta, deriv = "u1"), dcdu1_vine, tolerance = 1e-8)
  expect_equal(.copula_deriv(u1, u2, "G", theta, deriv = "u2"), dcdu2_vine, tolerance = 1e-8)
  expect_equal(.copula_deriv(u1, u2, "G", theta, deriv = "par", log = FALSE), dcdtheta_vine, tolerance = 1e-8)
  expect_equal(.copula_deriv(u1, u2, "G", theta, deriv = "par", log = TRUE), dlogdtheta_vine, tolerance = 1e-8)
})

test_that("native Gumbel second derivatives match VineCopula", {
  skip_if_no_vinecopula()

  old_backend <- getOption("gamlss.longitudinal.copula_backend")
  on.exit(options(gamlss.longitudinal.copula_backend = old_backend), add = TRUE)

  u1 <- seq(0.12, 0.88, length.out = 6)
  u2 <- rev(seq(0.18, 0.82, length.out = 6))
  theta <- 1.7

  options(gamlss.longitudinal.copula_backend = "vinecopula")
  d2cdu1_vine <- .copula_deriv2(u1, u2, "G", theta, deriv = "u1")
  d2cdu2_vine <- .copula_deriv2(u1, u2, "G", theta, deriv = "u2")
  d2cdtheta_vine <- .copula_deriv2(u1, u2, "G", theta, deriv = "par")

  options(gamlss.longitudinal.copula_backend = "native")
  expect_equal(.copula_deriv2(u1, u2, "G", theta, deriv = "u1"), d2cdu1_vine, tolerance = 1e-5)
  expect_equal(.copula_deriv2(u1, u2, "G", theta, deriv = "u2"), d2cdu2_vine, tolerance = 1e-5)
  expect_equal(.copula_deriv2(u1, u2, "G", theta, deriv = "par"), d2cdtheta_vine, tolerance = 1e-5)
})

test_that("native remaining first derivatives match VineCopula", {
  skip_if_no_vinecopula()

  old_backend <- getOption("gamlss.longitudinal.copula_backend")
  on.exit(options(gamlss.longitudinal.copula_backend = old_backend), add = TRUE)

  u1 <- seq(0.12, 0.88, length.out = 6)
  u2 <- rev(seq(0.18, 0.82, length.out = 6))
  cases <- list(
    F = list(par = 3, par2 = 0, tolerance = 1e-8),
    J = list(par = 1.8, par2 = 0, tolerance = 1e-8),
    t = list(par = 0.45, par2 = 5, tolerance = 1e-8)
  )

  for (family in names(cases)) {
    par <- cases[[family]]$par
    par2 <- cases[[family]]$par2
    tolerance <- cases[[family]]$tolerance

    options(gamlss.longitudinal.copula_backend = "vinecopula")
    dcdu1_vine <- .copula_deriv(u1, u2, family, par, par2, deriv = "u1")
    dcdu2_vine <- .copula_deriv(u1, u2, family, par, par2, deriv = "u2")
    dcdpar_vine <- .copula_deriv(u1, u2, family, par, par2, deriv = "par", log = FALSE)
    dlogdpar_vine <- .copula_deriv(u1, u2, family, par, par2, deriv = "par", log = TRUE)

    options(gamlss.longitudinal.copula_backend = "native")
    expect_equal(.copula_deriv(u1, u2, family, par, par2, deriv = "u1"), dcdu1_vine, tolerance = tolerance)
    expect_equal(.copula_deriv(u1, u2, family, par, par2, deriv = "u2"), dcdu2_vine, tolerance = tolerance)
    expect_equal(.copula_deriv(u1, u2, family, par, par2, deriv = "par", log = FALSE), dcdpar_vine, tolerance = tolerance)
    expect_equal(.copula_deriv(u1, u2, family, par, par2, deriv = "par", log = TRUE), dlogdpar_vine, tolerance = tolerance)
  }
})

test_that("native t df first derivative matches VineCopula", {
  skip_if_no_vinecopula()

  old_backend <- getOption("gamlss.longitudinal.copula_backend")
  on.exit(options(gamlss.longitudinal.copula_backend = old_backend), add = TRUE)

  u1 <- seq(0.12, 0.88, length.out = 6)
  u2 <- rev(seq(0.18, 0.82, length.out = 6))
  rho <- 0.45
  df <- 5

  options(gamlss.longitudinal.copula_backend = "vinecopula")
  dcddf_vine <- .copula_deriv(u1, u2, "t", rho, df, deriv = "par2", log = FALSE)
  dlogddf_vine <- .copula_deriv(u1, u2, "t", rho, df, deriv = "par2", log = TRUE)

  options(gamlss.longitudinal.copula_backend = "native")
  expect_equal(.copula_deriv(u1, u2, "t", rho, df, deriv = "par2", log = FALSE), dcddf_vine, tolerance = 1e-5)
  expect_equal(.copula_deriv(u1, u2, "t", rho, df, deriv = "par2", log = TRUE), dlogddf_vine, tolerance = 1e-5)
})

test_that("native t multi-derivative helper matches scalar derivative calls", {
  u1 <- seq(0.12, 0.88, length.out = 8)
  u2 <- rev(seq(0.18, 0.82, length.out = 8))
  rho <- seq(0.25, 0.55, length.out = 8)
  df <- seq(3.5, 7, length.out = 8)

  out <- .copula_t_deriv_many(
    u1,
    u2,
    rho,
    df,
    derivs = c("u1", "u2", "par", "par2"),
    log = c(u1 = FALSE, u2 = FALSE, par = TRUE, par2 = TRUE)
  )

  expect_equal(out$u1, .copula_t_deriv(u1, u2, rho, df, deriv = "u1", log = FALSE))
  expect_equal(out$u2, .copula_t_deriv(u1, u2, rho, df, deriv = "u2", log = FALSE))
  expect_equal(out$par, .copula_t_deriv(u1, u2, rho, df, deriv = "par", log = TRUE))
  expect_equal(out$par2, .copula_t_deriv(u1, u2, rho, df, deriv = "par2", log = TRUE))
  expect_equal(attr(out, "density"), .copula_t_pdf(u1, u2, rho, df))
})

test_that("native t copula fit matches stored regression values", {
  suppressPackageStartupMessages({
    library(gamlss)
    library(gamlss.dist)
  })

  dat <- make_fixture_factor_time_interaction(n_subject = 12L)
  old_backend <- getOption("gamlss.longitudinal.copula_backend")
  on.exit(options(gamlss.longitudinal.copula_backend = old_backend), add = TRUE)

  options(gamlss.longitudinal.copula_backend = "native")
  invisible(utils::capture.output(
    fit <- suppressWarnings(
      gamlss.longitudinal::gamlss_longitudinal(
        dataset = dat,
        margin_dist = gamlss.dist::NO(),
        copula_dist = "t",
        time_var = "time_raw",
        subject_var = "id",
        mu.formula = "y ~ time_raw + gender + age",
        sigma.formula = "~ 1",
        nu.formula = "~ 1",
        tau.formula = "~ 1",
        theta.formula = "~ 1",
        zeta.formula = "~ 1",
        include_dlcopdpar = TRUE,
        method = "RS",
        use_backtracking = TRUE,
        warm_start_joint = FALSE,
        compute_vcov = FALSE,
        verbose = 0,
        max_outer_iter = 2L,
        max_inner_iter = 2L,
        outer_stop_crit = 1,
        inner_stop_crit = 1
      )
    )
  ))

  expected_par <- c(
    zeta.intercept = 22.3883186214636,
    theta.intercept = 1.38438531171454,
    sigma.intercept = -0.389898704959138,
    mu.intercept = 2.85252346840033,
    mu.time_covariatet2 = 0.133086932999785,
    mu.time_covariatet3 = 0.15661826117056,
    mu.genderM = -0.00860902868362784,
    mu.age = 0.000691628158450122
  )
  expected_log_lik <- c(
    marginal = -34.2812776488741,
    copula = 8.94746681767027,
    joint = -25.3338108312039
  )

  expect_equal(fit$par, expected_par, tolerance = 1e-8)
  expect_equal(fit$calc_lik_out_end$log_lik, expected_log_lik, tolerance = 5e-4)
})

test_that("native remaining second derivatives match VineCopula", {
  skip_if_no_vinecopula()

  old_backend <- getOption("gamlss.longitudinal.copula_backend")
  on.exit(options(gamlss.longitudinal.copula_backend = old_backend), add = TRUE)

  u1 <- seq(0.12, 0.88, length.out = 6)
  u2 <- rev(seq(0.18, 0.82, length.out = 6))
  cases <- list(
    F = list(par = 3, par2 = 0, tolerance = 1e-5),
    J = list(par = 1.8, par2 = 0, tolerance = 1e-5),
    t = list(par = 0.45, par2 = 5, tolerance = 1e-5)
  )

  for (family in names(cases)) {
    par <- cases[[family]]$par
    par2 <- cases[[family]]$par2
    tolerance <- cases[[family]]$tolerance

    options(gamlss.longitudinal.copula_backend = "vinecopula")
    d2cdu1_vine <- .copula_deriv2(u1, u2, family, par, par2, deriv = "u1")
    d2cdu2_vine <- .copula_deriv2(u1, u2, family, par, par2, deriv = "u2")
    d2cdpar_vine <- .copula_deriv2(u1, u2, family, par, par2, deriv = "par")

    options(gamlss.longitudinal.copula_backend = "native")
    expect_equal(.copula_deriv2(u1, u2, family, par, par2, deriv = "u1"), d2cdu1_vine, tolerance = tolerance)
    expect_equal(.copula_deriv2(u1, u2, family, par, par2, deriv = "u2"), d2cdu2_vine, tolerance = tolerance)
    expect_equal(.copula_deriv2(u1, u2, family, par, par2, deriv = "par"), d2cdpar_vine, tolerance = tolerance)
  }
})

test_that("native t df second derivatives match VineCopula", {
  skip_if_no_vinecopula()

  old_backend <- getOption("gamlss.longitudinal.copula_backend")
  on.exit(options(gamlss.longitudinal.copula_backend = old_backend), add = TRUE)

  u1 <- seq(0.12, 0.88, length.out = 6)
  u2 <- rev(seq(0.18, 0.82, length.out = 6))
  rho <- 0.45
  df <- 5

  options(gamlss.longitudinal.copula_backend = "vinecopula")
  d2cddf_vine <- .copula_deriv2(u1, u2, "t", rho, df, deriv = "par2")
  d2cdrhoddf_vine <- .copula_deriv2(u1, u2, "t", rho, df, deriv = "par1par2")

  options(gamlss.longitudinal.copula_backend = "native")
  expect_equal(.copula_deriv2(u1, u2, "t", rho, df, deriv = "par2"), d2cddf_vine, tolerance = 1e-4)
  expect_equal(.copula_deriv2(u1, u2, "t", rho, df, deriv = "par1par2"), d2cdrhoddf_vine, tolerance = 1e-4)
})

test_that("native boundary derivatives are finite without delegated fallback", {
  old_backend <- getOption("gamlss.longitudinal.copula_backend")
  on.exit(options(gamlss.longitudinal.copula_backend = old_backend), add = TRUE)
  options(gamlss.longitudinal.copula_backend = "native")

  u1 <- c(0.2, 0.5, 0.8)
  u2 <- c(0.7, 0.4, 0.3)
  boundary <- list(C = 0, F = 0, G = 1, J = 1)

  for (family in names(boundary)) {
    par <- boundary[[family]]
    expect_equal(.copula_deriv(u1, u2, family, par, deriv = "u1"), rep(0, length(u1)))
    expect_equal(.copula_deriv(u1, u2, family, par, deriv = "u2"), rep(0, length(u1)))
    expect_true(all(is.finite(.copula_deriv(u1, u2, family, par, deriv = "par"))))
    expect_equal(.copula_deriv2(u1, u2, family, par, deriv = "u1"), rep(0, length(u1)))
    expect_equal(.copula_deriv2(u1, u2, family, par, deriv = "u2"), rep(0, length(u1)))
    expect_true(all(is.finite(.copula_deriv2(u1, u2, family, par, deriv = "par"))))
  }
})

test_that("native independence limits are handled for all supported copulas", {
  old_backend <- getOption("gamlss.longitudinal.copula_backend")
  on.exit(options(gamlss.longitudinal.copula_backend = old_backend), add = TRUE)
  options(gamlss.longitudinal.copula_backend = "native")

  u1 <- c(0.15, 0.35, 0.65, 0.9)
  u2 <- c(0.8, 0.55, 0.25, 0.1)
  independence_cases <- list(
    N = list(par = 0, par2 = 0),
    C = list(par = 0, par2 = 0),
    F = list(par = 0, par2 = 0),
    G = list(par = 1, par2 = 0),
    J = list(par = 1, par2 = 0)
  )

  for (family in names(independence_cases)) {
    par <- independence_cases[[family]]$par
    par2 <- independence_cases[[family]]$par2

    expect_equal(.copula_pdf(u1, u2, family, par, par2), rep(1, length(u1)), tolerance = 1e-10)
    expect_equal(.copula_cdf(u1, u2, family, par, par2), u1 * u2, tolerance = 1e-8)
    expect_equal(.copula_hfunc1(u1, u2, family, par, par2), u2, tolerance = 1e-8)
    expect_equal(.copula_deriv(u1, u2, family, par, par2, deriv = "u1"), rep(0, length(u1)), tolerance = 1e-8)
    expect_equal(.copula_deriv(u1, u2, family, par, par2, deriv = "u2"), rep(0, length(u1)), tolerance = 1e-8)
    expect_true(all(is.finite(.copula_deriv(u1, u2, family, par, par2, deriv = "par"))))
    expect_equal(.copula_deriv2(u1, u2, family, par, par2, deriv = "u1"), rep(0, length(u1)), tolerance = 1e-6)
    expect_equal(.copula_deriv2(u1, u2, family, par, par2, deriv = "u2"), rep(0, length(u1)), tolerance = 1e-6)
    expect_true(all(is.finite(.copula_deriv2(u1, u2, family, par, par2, deriv = "par"))))
  }

  # The t-copula is not independent at rho = 0 for finite df because the
  # shared scale still induces tail dependence. Independence is recovered in
  # the Gaussian limit as df grows.
  t_df <- 1e6
  expect_equal(.copula_pdf(u1, u2, "t", 0, t_df), rep(1, length(u1)), tolerance = 1e-5)
  expect_equal(.copula_cdf(u1, u2, "t", 0, t_df), u1 * u2, tolerance = 1e-6)
  expect_equal(.copula_hfunc1(u1, u2, "t", 0, t_df), u2, tolerance = 1e-6)
  expect_equal(.copula_deriv(u1, u2, "t", 0, t_df, deriv = "u1"), rep(0, length(u1)), tolerance = 1e-4)
  expect_equal(.copula_deriv(u1, u2, "t", 0, t_df, deriv = "u2"), rep(0, length(u1)), tolerance = 1e-4)
  expect_true(all(is.finite(.copula_deriv(u1, u2, "t", 0, t_df, deriv = "par"))))
  expect_equal(.copula_deriv2(u1, u2, "t", 0, t_df, deriv = "u1"), rep(0, length(u1)), tolerance = 1e-4)
  expect_equal(.copula_deriv2(u1, u2, "t", 0, t_df, deriv = "u2"), rep(0, length(u1)), tolerance = 1e-4)
  expect_true(all(is.finite(.copula_deriv2(u1, u2, "t", 0, t_df, deriv = "par"))))
})
