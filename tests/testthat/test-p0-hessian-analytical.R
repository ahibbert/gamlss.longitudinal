block_rel_frob <- function(A, B) {
  num <- sqrt(sum((A - B)^2, na.rm = TRUE))
  den <- sqrt(sum(B^2, na.rm = TRUE))
  num / (den + 1e-12)
}

make_del_count_fixture <- function(n_subject = 36L) {
  set.seed(2048)

  subject_tbl <- data.frame(
    id = seq_len(n_subject),
    gender = factor(sample(c("F", "M"), n_subject, replace = TRUE)),
    age = round(runif(n_subject, min = 20, max = 70), 1),
    stringsAsFactors = FALSE
  )

  time_levels <- c("t1", "t2", "t3")
  grid <- expand.grid(
    id = subject_tbl$id,
    time_raw = factor(time_levels, levels = time_levels, ordered = TRUE),
    KEEP.OUT.ATTRS = FALSE,
    stringsAsFactors = FALSE
  )

  dat <- merge(grid, subject_tbl, by = "id", sort = FALSE)
  dat <- dat[order(dat$id, dat$time_raw), ]
  rownames(dat) <- NULL

  t_num <- as.integer(dat$time_raw)
  g_num <- ifelse(dat$gender == "M", 1, 0)
  mu <- exp(1.2 + 0.12 * t_num + 0.18 * g_num + 0.003 * dat$age)
  sigma <- exp(-0.45 + 0.05 * t_num)
  nu <- plogis(-0.7 + 0.15 * g_num)
  dat$y <- gamlss.dist::rDEL(nrow(dat), mu = mu, sigma = sigma, nu = nu)

  dat
}

make_zip_count_fixture <- function(n_subject = 18L) {
  set.seed(4097)

  subject_tbl <- data.frame(
    id = seq_len(n_subject),
    gender = factor(sample(c("F", "M"), n_subject, replace = TRUE)),
    age = round(runif(n_subject, min = 20, max = 70), 1),
    stringsAsFactors = FALSE
  )

  time_levels <- c("t1", "t2", "t3")
  grid <- expand.grid(
    id = subject_tbl$id,
    time_raw = factor(time_levels, levels = time_levels, ordered = TRUE),
    KEEP.OUT.ATTRS = FALSE,
    stringsAsFactors = FALSE
  )

  dat <- merge(grid, subject_tbl, by = "id", sort = FALSE)
  dat <- dat[order(dat$id, dat$time_raw), ]
  rownames(dat) <- NULL

  t_num <- as.integer(dat$time_raw)
  g_num <- ifelse(dat$gender == "M", 1, 0)
  mu <- exp(0.9 + 0.08 * t_num + 0.10 * g_num + 0.002 * dat$age)
  sigma <- plogis(-0.5 + 0.2 * g_num)
  dat$y <- gamlss.dist::rZIP(nrow(dat), mu = mu, sigma = sigma)

  dat
}

test_that("T009 analytical Hessian warns for near-boundary GG and tracks numerical blocks", {
  skip_on_cran()

  suppressPackageStartupMessages({
    library(gamlss)
    library(gamlss.dist)
  })

  dat <- make_fixture_factor_time_interaction(n_subject = 40L)
  fit <- suppressWarnings(gamlss.longitudinal::gamlss_longitudinal(
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
  ))

  vc_num <- suppressWarnings(vcov(fit, method = "numderiv", progress = FALSE))
  vc_ana <- NULL
  expect_warning(
    vc_ana <- vcov(fit, method = "analytical", progress = FALSE),
    "GG may be numerically unstable because fitted nu is close to 0"
  )

  se_num <- suppressWarnings(sqrt(diag(vc_num$vcov$overall)))
  se_ana <- suppressWarnings(sqrt(diag(vc_ana$vcov$overall)))

  valid_se <- is.finite(se_num) & is.finite(se_ana) & se_num > 0 & se_ana > 0
  expect_gte(sum(valid_se), ceiling(0.8 * length(se_num)))

  rel_se <- abs(se_ana[valid_se] - se_num[valid_se]) / pmax(abs(se_num[valid_se]), 1e-6)

  # This fixture deliberately sits close to the GG nu = 0 limiting case, where
  # small Hessian differences can be amplified after matrix inversion.
  expect_lte(max(rel_se, na.rm = TRUE), 0.25)

  H_num <- -solve(vc_num$vcov$overall)
  H_ana <- -solve(vc_ana$vcov$overall)

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

test_that("T010 DEL analytical Hessian supports exact rectangle likelihood", {
  skip_on_cran()

  suppressPackageStartupMessages({
    library(gamlss)
    library(gamlss.dist)
  })

  fit <- gamlss.longitudinal::gamlss_longitudinal(
    dataset = make_del_count_fixture(n_subject = 36L),
    margin_dist = gamlss.dist::DEL(),
    copula_dist = "N",
    time_var = "time_raw",
    subject_var = "id",
    mu.formula = "y ~ time_raw + gender + age",
    sigma.formula = "~ time_raw",
    nu.formula = "~ gender",
    tau.formula = "~ 1",
    theta.formula = "~ time_raw",
    zeta.formula = "~ 1",
    include_dlcopdpar = TRUE,
    compute_vcov = FALSE,
    use_backtracking = TRUE,
    warm_start_joint = FALSE,
    verbose = 0,
    max_outer_iter = 5L,
    max_inner_iter = 5L,
    outer_stop_crit = 0.5,
    inner_stop_crit = 0.5
  )

  vc_ana <- suppressWarnings(vcov(fit, method = "analytical", progress = FALSE))

  expect_equal(vc_ana$method, "analytical")
  expect_true(all(is.finite(diag(vc_ana$vcov$overall))))
})

test_that("T011 analytical Hessian warns for zero-heavy discrete margins", {
  skip_on_cran()

  suppressPackageStartupMessages({
    library(gamlss)
    library(gamlss.dist)
  })

  fit <- suppressWarnings(gamlss.longitudinal::gamlss_longitudinal(
    dataset = make_zip_count_fixture(n_subject = 18L),
    margin_dist = gamlss.dist::ZIP(),
    copula_dist = "N",
    time_var = "time_raw",
    subject_var = "id",
    mu.formula = "y ~ time_raw + gender + age",
    sigma.formula = "~ gender",
    tau.formula = "~ 1",
    theta.formula = "~ time_raw",
    zeta.formula = "~ 1",
    include_dlcopdpar = TRUE,
    compute_vcov = FALSE,
    use_backtracking = TRUE,
    warm_start_joint = FALSE,
    verbose = 0,
    max_outer_iter = 3L,
    max_inner_iter = 3L,
    outer_stop_crit = 0.5,
    inner_stop_crit = 0.5
  ))

  vc_ana <- NULL
  expect_warning(
    vc_ana <- vcov(fit, method = "analytical", progress = FALSE),
    "zero-heavy discrete margins may be numerically delicate"
  )
  expect_equal(vc_ana$method, "analytical")
})
