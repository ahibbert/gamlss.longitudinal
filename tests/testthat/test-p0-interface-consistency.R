test_that("primary fitting interface is compact and uses clear names", {
  expect_identical(
    names(formals(gamlss_longitudinal)),
    c(
      "data", "margin_dist", "copula_dist", "time_var", "subject_var",
      "missingness", "mu.formula", "sigma.formula", "nu.formula",
      "tau.formula", "theta.formula", "zeta.formula", "method",
      "optimization", "optimizer_control", "start", "compute_vcov",
      "vcov_method", "margin_args", "verbose", "..."
    )
  )
})

test_that("joint and separate optimization names map to the established engine", {
  dat <- make_fixture_factor_time_interaction(n_subject = 12L)
  control <- gamlss_longitudinal_control(
    outer_tol = 1,
    max_outer_iter = 2L,
    rs = list(inner_tol = 1, max_inner_iter = 2L)
  )
  fit_one <- function(strategy) suppressWarnings(gamlss_longitudinal(
    data = dat,
    margin_dist = gamlss.dist::NO(),
    copula_dist = "N",
    time_var = "time_raw",
    subject_var = "id",
    mu.formula = y ~ time_raw,
    sigma.formula = ~1,
    theta.formula = ~1,
    optimization = strategy,
    optimizer_control = control,
    compute_vcov = FALSE
  ))

  joint <- fit_one("joint")
  separate <- fit_one("separate")
  expect_true(joint$include_dlcopdpar)
  expect_false(separate$include_dlcopdpar)
  expect_identical(joint$optimization, "joint")
  expect_identical(separate$convergence$optimization, "separate")
})

test_that("legacy fit aliases warn clearly and conflicting names error", {
  old <- options(gamlss.longitudinal.lifecycle_verbosity = "warning")
  on.exit(options(old), add = TRUE)
  .gl_deprecation_state$dataset_argument <- FALSE
  .gl_deprecation_state$include_dlcopdpar_argument <- FALSE

  dat <- make_fixture_factor_time_interaction(n_subject = 8L)
  args <- list(
    margin_dist = gamlss.dist::NO(), copula_dist = "N",
    time_var = "time_raw", subject_var = "id", mu.formula = y ~ 1,
    optimizer_control = gamlss_longitudinal_control(
      outer_tol = 1, max_outer_iter = 1L,
      rs = list(inner_tol = 1, max_inner_iter = 1L)
    ),
    compute_vcov = FALSE
  )

  expect_warning(
    suppressWarnings(do.call(gamlss_longitudinal, c(list(dataset = dat), args))),
    NA
  )
  .gl_deprecation_state$dataset_argument <- FALSE
  expect_warning(
    try(do.call(gamlss_longitudinal, c(list(dataset = dat), args)), silent = TRUE),
    "`dataset` is deprecated"
  )
  expect_error(
    do.call(gamlss_longitudinal, c(list(data = dat, dataset = dat), args)),
    "Supply only `data`"
  )

  start_legacy <- gamlss.longitudinal:::.gl_resolve_legacy_fit_dots(
    list(start_from = c(`mu.(Intercept)` = 1))
  )
  .gl_deprecation_state$start_from_argument <- FALSE
  expect_warning(
    resolved_start <- gamlss.longitudinal:::.gl_resolve_fit_start_alias(
      NULL, FALSE, start_legacy
    ),
    "`start_from` is deprecated"
  )
  expect_identical(resolved_start, c(`mu.(Intercept)` = 1))
  expect_error(
    gamlss.longitudinal:::.gl_resolve_fit_start_alias(
      c(`mu.(Intercept)` = 1), TRUE, start_legacy
    ),
    "Supply only `start`"
  )

  optimization_legacy <- gamlss.longitudinal:::.gl_resolve_legacy_fit_dots(
    list(include_dlcopdpar = FALSE)
  )
  .gl_deprecation_state$include_dlcopdpar_argument <- FALSE
  expect_warning(
    resolved_optimization <- gamlss.longitudinal:::.gl_resolve_optimization(
      "RS", c("joint", "separate"), FALSE, optimization_legacy
    ),
    "`include_dlcopdpar` is deprecated"
  )
  expect_identical(resolved_optimization, "separate")
  expect_error(
    gamlss.longitudinal:::.gl_resolve_optimization(
      "RS", "joint", TRUE, optimization_legacy
    ),
    "Supply only `optimization`"
  )
  for (bad in list(NA, "yes", c(TRUE, FALSE))) {
    invalid_legacy <- gamlss.longitudinal:::.gl_resolve_legacy_fit_dots(
      list(include_dlcopdpar = bad)
    )
    expect_error(
      gamlss.longitudinal:::.gl_resolve_optimization(
        "RS", c("joint", "separate"), FALSE, invalid_legacy
      ),
      "must be either TRUE or FALSE"
    )
  }
  expect_error(
    gamlss_longitudinal(
      data = dat, margin_dist = gamlss.dist::NO(), copula_dist = "N",
      time_var = "time_raw", subject_var = "id", method = "CG",
      optimization = "separate", compute_vcov = FALSE
    ),
    "CG always optimizes the joint objective"
  )
})

test_that("default fitting is quiet and fitted-model print is compact", {
  dat <- make_fixture_factor_time_interaction(n_subject = 10L)
  control <- gamlss_longitudinal_control(
    outer_tol = 1, max_outer_iter = 2L,
    rs = list(inner_tol = 1, max_inner_iter = 2L)
  )
  fit_output <- capture.output(
    fit <- suppressWarnings(gamlss_longitudinal(
      data = dat,
      margin_dist = gamlss.dist::NO(),
      copula_dist = "N",
      time_var = "time_raw",
      subject_var = "id",
      mu.formula = y ~ 1,
      optimizer_control = control,
      compute_vcov = FALSE
    ))
  )
  expect_identical(fit_output, character())
  verbose_output <- capture.output(suppressWarnings(gamlss_longitudinal(
    data = dat,
    margin_dist = gamlss.dist::NO(),
    copula_dist = "N",
    time_var = "time_raw",
    subject_var = "id",
    mu.formula = y ~ 1,
    optimizer_control = control,
    compute_vcov = FALSE,
    verbose = TRUE
  )))
  expect_gt(length(verbose_output), 0L)
  printed <- capture.output(print(fit))
  expect_lte(length(printed), 8L)
  expect_true(any(grepl("<gamlss_longitudinal>", printed, fixed = TRUE)))
  expect_true(any(grepl("Optimizer: RS, joint", printed, fixed = TRUE)))
  expect_false(any(grepl("par_history|dataset|calc_lik", printed)))
})

test_that("optimizer controls print as a narrow, scannable hierarchy", {
  control <- gamlss_longitudinal_control(
    outer_tol = 1e-4,
    max_outer_iter = 200L,
    max_elapsed_sec = 120,
    stop_on_convergence = FALSE,
    rs = list(
      inner_tol = 1e-5, max_inner_iter = 150L, max_steps = 20L,
      warm_start_joint = FALSE, backtracking_max_halves = 80L,
      update_lambda = FALSE, smooth_trust_radius = 2,
      discrete_score_method = "finite"
    )
  )
  printed <- capture.output(print(control))

  expect_lte(max(nchar(printed, type = "width")), 80L)
  expect_true(any(grepl("Shared:", printed, fixed = TRUE)))
  expect_true(any(grepl("RS overrides:", printed, fixed = TRUE)))
  expect_true(any(grepl("  inner_tol: 1e-05", printed, fixed = TRUE)))
})

test_that("selection displays are compact and best rows do not retain diagnostics", {
  set.seed(310)
  u1 <- stats::runif(80)
  u2 <- pmin(pmax(0.7 * u1 + 0.3 * stats::runif(80), 1e-6), 1 - 1e-6)
  result <- select_copula(u1 = u1, u2 = u2, families = c("N", "C"), min_pairs = 10)
  printed <- capture.output(print(result))
  best <- best_fit(result)

  expect_lte(max(nchar(printed, type = "width")), 100L)
  expect_true(any(grepl("Gaussian.*NA", printed)))
  expect_identical(attributes(best), list(names = names(best)))
  expect_true(all(c("rank", "family", "delta") %in% names(result)))
})

test_that("package plots share a theme and copula dashboard has focused views", {
  expect_s3_class(theme_gamlss_longitudinal(), "theme")
  fixture <- get_valid_plot_fixture()
  overview <- suppressWarnings(plot_copula_diagnostics(fixture$fit, plot = FALSE))
  all_panels <- suppressWarnings(plot_copula_diagnostics(fixture$fit, type = "all", plot = FALSE))
  positional_lag <- suppressWarnings(plot_copula_diagnostics(fixture$fit, 2, plot = FALSE))

  expect_s3_class(overview, "gamlss_longitudinal_plot")
  expect_identical(overview$type, "overview")
  expect_length(overview$displayed_plots, 4L)
  expect_length(all_panels$displayed_plots, 9L)
  expect_s3_class(positional_lag, "gamlss_longitudinal_plot")
  expect_true(all(vapply(overview$displayed_plots, inherits, logical(1), "ggplot")))
})

test_that("plot_dist temporarily accepts dataset without weakening argument checks", {
  old <- options(gamlss.longitudinal.lifecycle_verbosity = "warning")
  on.exit(options(old), add = TRUE)
  .gl_deprecation_state$plot_dist_dataset_argument <- FALSE
  dat <- data.frame(
    subject = rep(1:4, each = 2),
    time = rep(1:2, times = 4),
    response = seq_len(8)
  )

  expect_warning(
    out <- plot_dist(
      dataset = dat,
      subject_var = "subject",
      time_var = "time",
      response_var = "response"
    ),
    "`dataset` is deprecated"
  )
  expect_false(is.null(out))
  expect_error(
    plot_dist(
      data = dat,
      dataset = dat,
      subject_var = "subject",
      time_var = "time",
      response_var = "response"
    ),
    "Supply only `data`"
  )
  expect_error(plot_dist(data = dat, mystery = TRUE), "Unused argument: `mystery`")
  expect_error(plot_dist(data = dat, dist = gamlss.dist::NO()), "margin_dist")
  duplicate_call <- c(
    list(subject_var = "subject", time_var = "time", response_var = "response"),
    structure(list(dat, dat), names = c("dataset", "dataset"))
  )
  expect_error(do.call(plot_dist, duplicate_call), "must not be duplicated")
  expect_error(plot_dist(), "Supply `data`")
  expect_error(
    plot_dist(
      data = dat,
      subject_var = "absent",
      time_var = "time",
      response_var = "response"
    ),
    "not found in `data`"
  )
})
