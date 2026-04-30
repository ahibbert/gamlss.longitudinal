test_that("T009 analytical Hessian tracks numerical Hessian", {
  skip_on_cran()

  suppressPackageStartupMessages({
    library(gamlss)
    library(gamlss.dist)
    library(VineCopula)
  })

  dat <- make_fixture_factor_time_interaction(n_subject = 40L)
  fit <- gamlss.longitudinal::gamlss.longitudinal(
    dataset = dat,
    margin_dist = gamlss.dist::GG(),
    copula_dist = "N",
    time_var = "time_raw",
    subject_var = "id",
    mu.formula = "y ~ time_raw + gender + age",
    sigma.formula = "~ time_raw + gender",
    nu.formula = "~ 1",
    tau.formula = "~ 1",
    theta.formula = "~ time_raw",
    zeta.formula = "~ 1",
    include_dlcopdpar = TRUE,
    use_backtracking = TRUE,
    verbose = 0,
    max_outer_iter = 4L,
    max_inner_iter = 4L,
    outer_stop_crit = 0.5,
    inner_stop_crit = 0.5
  )

  vc_num <- suppressWarnings(vcov(fit, method = "numderiv", progress = FALSE))
  vc_ana <- suppressWarnings(vcov(fit, method = "analytical", progress = FALSE))

  se_num <- suppressWarnings(sqrt(diag(vc_num$vcov$overall)))
  se_ana <- suppressWarnings(sqrt(diag(vc_ana$vcov$overall)))

  valid_se <- is.finite(se_num) & is.finite(se_ana) & se_num > 0 & se_ana > 0
  expect_gte(sum(valid_se), ceiling(0.8 * length(se_num)))

  rel_se <- abs(se_ana[valid_se] - se_num[valid_se]) / pmax(abs(se_num[valid_se]), 1e-6)

  # Keep this threshold tight enough to catch Hessian regressions while
  # allowing small optimizer and finite-difference noise.
  expect_lte(max(rel_se, na.rm = TRUE), 0.20)

  H_num <- -solve(vc_num$vcov$overall)
  H_ana <- -solve(vc_ana$vcov$overall)

  block_rel_frob <- function(A, B) {
    num <- sqrt(sum((A - B)^2, na.rm = TRUE))
    den <- sqrt(sum(B^2, na.rm = TRUE))
    num / (den + 1e-12)
  }

  theta_rows <- grepl("theta", rownames(H_num))
  mu_rows <- grepl("^mu\\.", rownames(H_num))
  sigma_rows <- grepl("^sigma\\.", rownames(H_num))
  nu_rows <- grepl("^nu\\.", rownames(H_num))

  block_err <- c(
    theta_theta = block_rel_frob(H_ana[theta_rows, theta_rows, drop = FALSE], H_num[theta_rows, theta_rows, drop = FALSE]),
    mu_mu = block_rel_frob(H_ana[mu_rows, mu_rows, drop = FALSE], H_num[mu_rows, mu_rows, drop = FALSE]),
    sigma_sigma = block_rel_frob(H_ana[sigma_rows, sigma_rows, drop = FALSE], H_num[sigma_rows, sigma_rows, drop = FALSE]),
    nu_nu = block_rel_frob(H_ana[nu_rows, nu_rows, drop = FALSE], H_num[nu_rows, nu_rows, drop = FALSE]),
    mu_sigma = block_rel_frob(H_ana[mu_rows, sigma_rows, drop = FALSE], H_num[mu_rows, sigma_rows, drop = FALSE])
  )

  expect_true(all(is.finite(block_err)))
  expect_lte(max(block_err, na.rm = TRUE), 0.15)
})
