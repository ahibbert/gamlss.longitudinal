test_that("capability registry is versioned, complete, and documentation-ready", {
  margins <- longitudinal_capabilities("margins")
  copulas <- longitudinal_capabilities("copulas")
  routes <- longitudinal_capabilities("routes")

  expect_s3_class(margins, "data.frame")
  expect_s3_class(copulas, "data.frame")
  expect_s3_class(routes, "data.frame")
  expect_true(all(nzchar(margins$registry_version)))
  expect_equal(unique(margins$registry_version), unique(routes$registry_version))
  expect_equal(unique(copulas$registry_version), unique(routes$registry_version))
  expect_named(margins, c(
    "registry_version", "family", "status", "family_type", "response_domain",
    "response_requirement", "parameters", "default_links", "inverse_links",
    "likelihood_route", "compatible_copulas", "hessian_path", "hessian_validity_guaranteed",
    "curvature_validation", "hessian_evidence", "randomized_pit",
    "diagnostics", "paper_route", "limitations"
  ))
  expect_named(copulas, c(
    "registry_version", "copula", "status", "parameters", "links",
    "inverse_links", "parameter_domain", "limitations"
  ))
  expect_true(all(c(
    "NO", "GA", "GG", "BCPE", "LOGNO", "PO", "NBI", "DEL", "ZIP", "ZAP", "ZINBI",
    "SICHEL", "ZISICHEL", "NBII", "PIG", "JSU", "SHASH"
  ) %in% margins$family[margins$status == "supported"]))
  expect_equal(
    sort(unique(routes$margin_family)),
    sort(margins$family[margins$status == "supported"])
  )
  expect_equal(copulas$copula, c("N", "C", "F", "G", "J", "t"))
  expect_true(all(nzchar(margins$response_requirement)))
  expect_true(all(nzchar(margins$limitations)))
  expect_true(all(nzchar(copulas$parameter_domain)))
  expect_true(all(routes$curvature_validation == "required_per_fit_jss002"))
  expect_false(any(routes$hessian_validity_guaranteed))
  expect_true(all(nzchar(routes$route_evidence)))
  expect_true(all(nzchar(routes$hessian_evidence)))
  expect_true(all(nzchar(routes$diagnostics)))
  expect_true(all(nzchar(routes$diagnostic_evidence)))
  expect_true(all(grepl("not guaranteed valid inference", routes$limitations, fixed = TRUE)))
})

test_that("registry identifies both mandatory public paper routes", {
  routes <- longitudinal_capabilities("routes")
  paper <- routes[routes$paper_route, c("margin_family", "copula", "likelihood_route")]

  expect_equal(nrow(paper), 2L)
  expect_true(any(paper$margin_family == "BCPE" & paper$copula == "t" & paper$likelihood_route == "continuous_density"))
  expect_true(any(paper$margin_family == "NBI" & paper$copula == "C" & paper$likelihood_route == "exact_discrete_rectangle"))
})

test_that("capability tables are direct projections of the registry specs", {
  margins <- longitudinal_capabilities("margins")
  routes <- longitudinal_capabilities("routes")
  supported_specs <- Filter(
    function(spec) identical(spec$status, "supported"),
    gamlss.longitudinal:::.gl_capability_all_margin_specs()
  )

  expect_equal(sort(margins$family[margins$status == "supported"]), sort(names(supported_specs)))
  for (family in names(supported_specs)) {
    expected <- supported_specs[[family]]$compatible_copulas
    observed <- routes$copula[routes$margin_family == family]
    expect_equal(observed, expected, info = family)
  }

  rd_path <- test_path("..", "..", "man", "longitudinal_capabilities.Rd")
  skip_if_not(file.exists(rd_path), "Raw source documentation is unavailable in installed-package tests.")
  rd <- readLines(rd_path, warn = FALSE)
  expect_true(any(grepl("BCPE", rd, fixed = TRUE)))
  expect_true(any(grepl("NBI", rd, fixed = TRUE)))
  expect_true(any(grepl("one homogeneous marginal family", tolower(rd), fixed = TRUE)))
})

test_that("registered routes and likelihood routing agree", {
  skip_if_not_installed("gamlss.dist")

  expect_invisible(gamlss.longitudinal:::.gl_validate_capability_route(
    gamlss.dist::BCPE(), "t", response = c(1, 2)
  ))
  expect_invisible(gamlss.longitudinal:::.gl_validate_capability_route(
    gamlss.dist::NBI(), "C", response = c(0, 2)
  ))
  expect_invisible(gamlss.longitudinal:::.gl_validate_capability_route(
    gamlss.dist::NBI(), "G", response = c(0, 2)
  ))
  expect_invisible(gamlss.longitudinal:::.gl_validate_capability_route(
    gamlss.dist::GG(), "F", response = c(1, 2)
  ))
  expect_identical(
    gamlss.longitudinal:::.gl_capability_likelihood_route(gamlss.dist::BCPE()),
    "continuous_density"
  )
  expect_identical(
    gamlss.longitudinal:::.gl_capability_likelihood_route(gamlss.dist::NBI()),
    "exact_discrete_rectangle"
  )
})

test_that("compatible pure gamlss.dist families use the generic Gaussian route", {
  skip_if_not_installed("gamlss.dist")

  for (family in c("SICHEL", "ZISICHEL", "NBII", "PIG", "JSU", "SHASH")) {
    object <- do.call(get(family, envir = asNamespace("gamlss.dist")), list())
    spec <- gamlss.longitudinal:::.gl_capability_margin_spec(family)
    expect_identical(spec$status, "supported", info = family)
    expect_identical(spec$compatible_copulas, "N", info = family)
    expect_invisible(gamlss.longitudinal:::.gl_validate_capability_route(
      object, "N"
    ))
  }

  expect_invisible(gamlss.longitudinal:::.gl_validate_capability_route(
    gamlss.dist::LNO(), "N"
  ))
  expect_error(
    gamlss.longitudinal:::.gl_validate_capability_route(
      gamlss.dist::SICHEL(), "N", response = c(0, -1)
    ),
    class = "gamlss_longitudinal_response_domain_error"
  )
})

test_that("mixed-measure and logarithmic-series families remain excluded", {
  skip_if_not_installed("gamlss.dist")

  expect_error(
    gamlss.longitudinal:::.gl_validate_capability_route(
      gamlss.dist::BEINF(), "N"
    ),
    "mixed continuous/discrete",
    class = "gamlss_longitudinal_mixed_margin_error"
  )
  expect_error(
    gamlss.longitudinal:::.gl_validate_capability_route(
      gamlss.dist::ZAGA(), "N"
    ),
    "mixed continuous/discrete",
    class = "gamlss_longitudinal_mixed_margin_error"
  )
  expect_error(
    gamlss.longitudinal:::.gl_validate_capability_route(
      gamlss.dist::LG(), "N"
    ),
    "logarithmic-series",
    class = "gamlss_longitudinal_unsupported_margin_error"
  )
})

test_that("missing family arguments, mixed routes, and invalid domains fail clearly", {
  skip_if_not_installed("gamlss.dist")

  expect_error(
    gamlss.longitudinal:::.gl_validate_capability_route(gamlss.dist::BI(), "N"),
    NA
  )
  expect_error(
    gamlss.longitudinal:::.gl_validate_capability_route(gamlss.dist::BCPE(), "C"),
    "tested allowlist", class = "gamlss_longitudinal_unsupported_route_error"
  )
  expect_error(
    gamlss.longitudinal:::.gl_validate_capability_route(list(gamlss.dist::NO(), gamlss.dist::PO()), "N"),
    "one homogeneous", class = "gamlss_longitudinal_mixed_margin_error"
  )
  expect_error(
    gamlss.longitudinal:::.gl_validate_capability_route(gamlss.dist::GA(), "N", response = c(0, 1)),
    "strictly greater", class = "gamlss_longitudinal_response_domain_error"
  )
  expect_error(
    gamlss.longitudinal:::.gl_validate_capability_route(gamlss.dist::NBI(), "C", response = c(0, 1.5)),
    "non-negative integers", class = "gamlss_longitudinal_response_domain_error"
  )
  expect_error(
    gamlss.longitudinal:::.gl_validate_capability_route(gamlss.dist::NO(), "N", response = c(1, Inf)),
    "Inf", class = "gamlss_longitudinal_response_domain_error"
  )
  expect_error(
    gamlss.longitudinal:::.gl_validate_capability_route(gamlss.dist::NO(), "N", response = c(-Inf, NA)),
    "Inf", class = "gamlss_longitudinal_response_domain_error"
  )
  expect_invisible(
    gamlss.longitudinal:::.gl_validate_capability_route(gamlss.dist::NO(), "N", response = c(NA, NA))
  )
})

test_that("fit preflight enforces binomial trials before the workflow", {
  skip_if_not_installed("gamlss.dist")

  dat <- data.frame(id = c(1, 1), time = c(1, 2), y = c(0, 1))
  expect_error(
    gamlss.longitudinal:::.gl_preflight_fit_capabilities(
      dat, gamlss.dist::BI(), "N", y ~ 1
    ),
    "explicit binomial trials",
    class = "gamlss_longitudinal_missing_margin_argument_error"
  )
  binomial <- gamlss.longitudinal:::.gl_attach_margin_args_spec(
    gamlss.dist::BI(), list(bd = "trials")
  )
  expect_invisible(gamlss.longitudinal:::.gl_preflight_fit_capabilities(
    transform(dat, trials = 2L), binomial, "N", y ~ 1
  ))
  expect_error(
    gamlss.longitudinal:::.gl_preflight_fit_capabilities(
      transform(dat, trials = c(1L, 0L)), binomial, "N", y ~ 1
    ),
    "positive integers"
  )
  expect_error(
    gamlss.longitudinal:::.gl_preflight_fit_capabilities(
      transform(dat, y = c(1, Inf)), gamlss.dist::NO(), "N", y ~ 1
    ),
    "Inf", class = "gamlss_longitudinal_response_domain_error"
  )
})

test_that("coverage catalog is audit evidence constrained by the registry", {
  skip_if_not_installed("gamlss.dist")

  catalog <- gamlss.longitudinal:::.coverage_family_catalog()
  supported <- sort(catalog$family[catalog$supported])
  registered <- sort(longitudinal_capabilities("margins")$family[
    longitudinal_capabilities("margins")$status == "supported"
  ])
  expect_true(all(supported %in% registered))
  expect_true(catalog$supported[match("BI", catalog$family)])
})

test_that("binomial margins use explicit constant or row-specific trials", {
  skip_if_not_installed("gamlss.dist")

  set.seed(842L)
  dat <- expand.grid(time = seq_len(4L), subject = seq_len(20L))
  dat$trials <- 5L + dat$subject %% 4L
  dat$y <- stats::rbinom(
    nrow(dat), dat$trials, stats::plogis(-0.5 + 0.15 * dat$time)
  )
  control <- gamlss_longitudinal_control(
    max_outer_iter = 1L,
    stop_on_convergence = FALSE,
    rs = list(
      max_inner_iter = 1L,
      warm_start_joint = FALSE,
      update_lambda = FALSE
    )
  )

  fit <- suppressWarnings(gamlss_longitudinal(
    dataset = dat,
    margin_dist = gamlss.dist::BI(),
    copula_dist = "N",
    time_var = "time",
    subject_var = "subject",
    mu.formula = y ~ time,
    optimizer_control = control,
    compute_vcov = FALSE,
    verbose = 0,
    margin_args = list(bd = "trials")
  ))

  expect_s3_class(fit, "gamlss.longitudinal")
  expect_true(is.finite(fit$convergence$logLik))
  expect_equal(
    attr(fit$margin_dist, "fixed_family_args")$bd,
    fit$dataset$trials
  )
  expected_density <- gamlss.dist::dBI(
    fit$response,
    bd = fit$dataset$trials,
    mu = fitted(fit, parameter = "mu")
  )
  predicted_density <- suppressWarnings(predict(fit, type = "density"))
  expect_equal(predicted_density$density, expected_density, tolerance = 1e-10)
  expect_equal(
    suppressWarnings(predict(fit, type = "mean")),
    fit$dataset$trials * fitted(fit, parameter = "mu"),
    tolerance = 1e-10
  )
  copula_data <- gamlss.longitudinal:::.copula_v2_fit_data(fit)
  expect_equal(
    copula_data$u,
    gamlss.longitudinal:::.copula_v2_clamp01(gamlss.dist::pBI(
      fit$response,
      bd = fit$dataset$trials,
      mu = fitted(fit, parameter = "mu")
    )),
    tolerance = 1e-10
  )
  simulated <- simulate(fit, nsim = 2L, seed = 844L)
  expect_true(all(simulated$sim_1 >= 0 & simulated$sim_1 <= fit$dataset$trials))
  expect_true(all(simulated$sim_2 >= 0 & simulated$sim_2 <= fit$dataset$trials))
  simulated_new <- simulate(fit, nsim = 1L, seed = 845L, newdata = dat[1:8, ])
  expect_true(all(simulated_new$sim_1 >= 0 & simulated_new$sim_1 <= dat$trials[1:8]))

  constant <- gamlss.longitudinal:::.gl_attach_margin_args_spec(
    gamlss.dist::BI(), list(bd = 8L)
  )
  expect_invisible(gamlss.longitudinal:::.gl_preflight_fit_capabilities(
    transform(dat, y = pmin(y, 8L)), constant, "N", y ~ time
  ))
})

test_that("native simulation passes binomial denominators to the GAMLSS quantile", {
  skip_if_not_installed("gamlss.dist")

  dat <- simulate_longitudinal_dataset(
    n = 12L,
    times = 1:4,
    margin_dist = gamlss.dist::BI(),
    copula_dist = "N",
    margin_params = list(mu = 0.35, bd = 6L),
    copula_params = list(theta = 0.3),
    seed = 846L
  )
  expect_true(all(dat$response >= 0 & dat$response <= 6L))
  expect_equal(unique(dat$true_bd), 6)
})

test_that("binomial margin and copula selection retain the trials contract", {
  skip_if_not_installed("gamlss")
  skip_if_not_installed("gamlss.dist")

  set.seed(843L)
  dat <- expand.grid(time = seq_len(4L), subject = seq_len(15L))
  dat$trials <- 10L
  subject_effect <- stats::rnorm(15L, sd = 0.7)
  dat$y <- stats::rbinom(
    nrow(dat), dat$trials,
    stats::plogis(-0.8 + subject_effect[dat$subject])
  )
  margins <- suppressWarnings(suppressMessages(select_margin(
    dat,
    response_var = "y",
    families = c("BI", "BB"),
    margin_args = list(bd = "trials"),
    trace = FALSE
  )))
  selected <- best_fit_family(margins)
  expect_true(selected$family[1] %in% c("BI", "BB"))
  expect_identical(attr(selected, "margin_args_spec")$bd, "trials")

  copulas <- suppressWarnings(suppressMessages(select_copula(
    data = dat,
    response_var = "y",
    margin_dist = margins,
    subject_var = "subject",
    time_var = "time",
    families = "N",
    margin_n_cycles = 10L,
    min_pairs = 10L
  )))
  expect_identical(copulas$family, "N")
  expect_true(all(is.finite(copulas$logLik)))
  expect_true(all(copulas$screen_likelihood == "exact_discrete_rectangle"))
})
