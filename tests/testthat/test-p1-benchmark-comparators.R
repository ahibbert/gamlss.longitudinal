test_that("benchmark comparator status records optional backends", {
  status <- benchmark_comparator_status()

  expect_s3_class(status, "data.frame")
  expect_true(all(c("comparator", "estimator", "package", "available", "role") %in% names(status)))
  expect_equal(sort(status$comparator), sort(c("gee", "glmm", "gam", "glmmTMB")))
  expect_type(status$available, "logical")
})

test_that("adoption benchmark scenarios define executable benchmark claims", {
  scenarios <- adoption_benchmark_scenarios()

  expect_s3_class(scenarios, "data.frame")
  expect_true(all(c(
    "scenario", "family", "copula", "design", "n_subject", "n_time",
    "dependence", "missingness", "claim", "primary_metrics", "methods"
  ) %in% names(scenarios)))
  expect_true(all(lengths(scenarios$primary_metrics) > 0L))
  expect_true(all(lengths(scenarios$methods) > 0L))

  one <- adoption_benchmark_scenarios("gaussian_heteroskedastic")
  expect_equal(nrow(one), 1L)
  expect_equal(one$family, "NO")
  expect_equal(one$design, "scale")
  dependence <- adoption_benchmark_scenarios("time_varying_dependence")
  expect_equal(dependence$design, "time_dependence")
  expect_true("benchmark_theta_time_abs_error" %in% dependence$primary_metrics[[1L]])
  expect_error(adoption_benchmark_scenarios("not_a_scenario"), "Unknown adoption benchmark")
})

test_that("run_adoption_benchmarks executes a small opt-in scenario", {
  scenario <- adoption_benchmark_scenarios("gaussian_heteroskedastic")
  scenario$n_subject <- 8L
  scenario$n_time <- 3L

  bench <- run_adoption_benchmarks(
    scenarios = scenario,
    reps = 1L,
    methods = "gam",
    metrics = "benchmark_mean_rmse",
    max_elapsed_sec = 10,
    write_results = FALSE
  )

  expect_s3_class(bench, "gamlss_longitudinal_adoption_benchmark")
  expect_s3_class(bench$results, "data.frame")
  expect_s3_class(bench$summary, "gamlss_longitudinal_benchmark_summary")
  expect_equal(unique(bench$results$benchmark_scenario), "gaussian_heteroskedastic")
  expect_equal(unique(bench$results$method), "gam")
  expect_true("benchmark_mean_rmse" %in% bench$metrics)
  expect_true(nrow(bench$summary$summary) >= 1L)

  report <- format_benchmark_report(bench)
  expect_type(report, "character")
  expect_true(any(grepl("Headline Results", report, fixed = TRUE)))
  expect_true(any(grepl("Scenario-Level Headline Results", report, fixed = TRUE)))
  expect_true(any(grepl("Scenario-Level Win/Tie/Loss Summary", report, fixed = TRUE)))
  expect_true(any(grepl("Interpretation Notes", report, fixed = TRUE)))

  report_path <- tempfile(fileext = ".md")
  written <- write_benchmark_report(bench, path = report_path)
  expect_true(file.exists(report_path))
  expect_equal(normalizePath(report_path, winslash = "/", mustWork = FALSE), written)
  report_text <- readLines(report_path, warn = FALSE)
  expect_true(any(grepl("Comparator Availability", report_text, fixed = TRUE)))
  expect_true(any(grepl("benchmark_scenario", report_text, fixed = TRUE)))
  expect_true(any(grepl("declared primary metrics", report_text, fixed = TRUE)))
})

test_that("package ships an opt-in adoption benchmark runner", {
  script <- system.file(
    "benchmarks",
    "run-adoption-benchmarks.R",
    package = "gamlss.longitudinal"
  )
  expect_true(file.exists(script))

  text <- readLines(script, warn = FALSE)
  expect_true(any(grepl("run_adoption_benchmarks", text, fixed = TRUE)))
  expect_true(any(grepl("adoption_benchmark_scenarios", text, fixed = TRUE)))
  expect_true(any(grepl("write_benchmark_report", text, fixed = TRUE)))
  expect_true(any(grepl("GAMLSS_LONGITUDINAL_ADOPTION_REPS", text, fixed = TRUE)))
})

test_that("benchmark_standard_models fits the always-available GAM comparator", {
  dat <- make_fixture_factor_time_interaction(n_subject = 10L)
  dat$id <- factor(dat$id)

  bench <- benchmark_standard_models(
    data = dat,
    formula = y ~ time_raw + gender + age,
    subject_var = "id",
    family = "gaussian",
    comparators = "gam"
  )

  expect_s3_class(bench, "gamlss_longitudinal_benchmark")
  expect_equal(bench$results$comparator, "gam")
  expect_true(bench$results$available)
  expect_true(bench$results$success)
  expect_true(is.finite(bench$results$rmse))
  expect_true(inherits(bench$fits$gam, "gam"))
})

test_that("benchmark_standard_models fits GEE and GLMM comparators when installed", {
  skip_if_not_installed("geepack")
  skip_if_not_installed("lme4")

  dat <- make_fixture_factor_time_interaction(n_subject = 10L)
  dat$id <- factor(dat$id)

  bench <- benchmark_standard_models(
    data = dat,
    formula = y ~ time_raw + gender + age,
    subject_var = "id",
    family = "gaussian",
    comparators = c("gee", "glmm")
  )

  expect_equal(sort(bench$results$comparator), sort(c("gee", "glmm")))
  expect_true(all(bench$results$available))
  expect_true(all(bench$results$success))
  expect_true(all(is.finite(bench$results$mae)))
})

test_that("benchmark_standard_models records unavailable optional comparators", {
  dat <- make_fixture_factor_time_interaction(n_subject = 8L)
  dat$id <- factor(dat$id)

  bench <- benchmark_standard_models(
    data = dat,
    formula = y ~ time_raw + gender,
    subject_var = "id",
    family = "gaussian",
    comparators = "glmmTMB"
  )

  expect_equal(bench$results$comparator, "glmmTMB")
  expect_type(bench$results$available, "logical")
  expect_true(is.logical(bench$results$success))
  if (!bench$results$available) {
    expect_false(bench$results$success)
    expect_match(bench$results$error, "not installed")
  }
})

test_that("summarise_benchmark_results creates win/tie/loss summaries", {
  results <- data.frame(
    family = "NO",
    copula = "N",
    design = "covariate",
    n_subject = 20L,
    n_time = 3L,
    dependence = "moderate",
    missingness = "none",
    start_mode = "default",
    method = c("gee", "glmm", "gam"),
    benchmark_mean_rmse = c(0.30, 0.20, 0.21),
    benchmark_interval_coverage_95 = c(0.80, 0.94, 0.96),
    elapsed_sec = c(0.1, 0.3, 0.2),
    stringsAsFactors = FALSE
  )

  out <- summarise_benchmark_results(results)

  expect_s3_class(out, "gamlss_longitudinal_benchmark_summary")
  expect_s3_class(out$summary, "data.frame")
  expect_s3_class(out$case_results, "data.frame")
  expect_true(all(c("metric", "estimand", "method", "wins", "ties", "losses", "win_or_tie_rate") %in% names(out$summary)))
  mean_cases <- out$case_results[out$case_results$metric == "benchmark_mean_rmse", , drop = FALSE]
  expect_equal(mean_cases$result[match("glmm", mean_cases$method)], "win")
  expect_equal(mean_cases$result[match("gee", mean_cases$method)], "loss")
  coverage_cases <- out$case_results[out$case_results$metric == "benchmark_interval_coverage_95", , drop = FALSE]
  expect_true(coverage_cases$result[match("glmm", coverage_cases$method)] %in% c("win", "tie"))
  expect_true(coverage_cases$result[match("gam", coverage_cases$method)] %in% c("win", "tie"))
})
