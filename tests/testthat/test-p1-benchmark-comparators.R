test_that("benchmark comparator status records optional backends", {
  status <- benchmark_comparator_status()

  expect_s3_class(status, "data.frame")
  expect_true(all(c("comparator", "estimator", "package", "available", "role") %in% names(status)))
  expect_equal(sort(status$comparator), sort(c("glm", "gee", "glmm", "gam", "gamm", "glmmTMB")))
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
  expect_true(any(grepl("Benchmark Summary", report, fixed = TRUE)))
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

test_that("benchmark_standard_models default comparators cover independence and longitudinal baselines", {
  dat <- make_fixture_factor_time_interaction(n_subject = 8L)
  dat$id <- factor(dat$id)

  bench <- benchmark_standard_models(
    data = dat,
    formula = y ~ time_raw + gender + age,
    subject_var = "id",
    family = "gaussian"
  )

  expect_equal(bench$results$comparator, c("glm", "gee", "glmm", "gam", "gamm"))
  expect_false("glmmTMB" %in% bench$results$comparator)
})

test_that("benchmark_standard_models fits GLM, GAM, and GAMM linear formula comparators", {
  dat <- make_fixture_factor_time_interaction(n_subject = 10L)
  dat$id <- factor(dat$id)

  bench <- benchmark_standard_models(
    data = dat,
    formula = y ~ time_raw + gender + age,
    subject_var = "id",
    family = "gaussian",
    comparators = c("glm", "gam", "gamm")
  )

  expect_equal(bench$results$comparator, c("glm", "gam", "gamm"))
  expect_true(all(bench$results$available))
  expect_true(all(bench$results$success))
  expect_true(inherits(bench$fits$glm, "glm"))
  expect_true(inherits(bench$fits$gam, "gam"))
  expect_true(inherits(bench$fits$gamm, "gam"))
})

test_that("benchmark_standard_models supports smooths for GAM and GAMM only", {
  dat <- make_fixture_factor_time_interaction(n_subject = 10L)
  dat$id <- factor(dat$id)

  bench <- benchmark_standard_models(
    data = dat,
    formula = y ~ time_raw + gender + s(age, bs = "ps"),
    subject_var = "id",
    family = "gaussian",
    comparators = c("glm", "gam", "gamm")
  )

  expect_false(bench$results$success[bench$results$comparator == "glm"])
  expect_match(
    bench$results$error[bench$results$comparator == "glm"],
    "Smooth terms via s\\(\\.\\.\\.\\) are only supported",
    perl = TRUE
  )
  expect_true(bench$results$success[bench$results$comparator == "gam"])
  expect_true(bench$results$success[bench$results$comparator == "gamm"])
})

test_that("benchmark_standard_models records smooth failures for GEE and GLMM when installed", {
  skip_if_not_installed("geepack")
  skip_if_not_installed("lme4")

  dat <- make_fixture_factor_time_interaction(n_subject = 10L)
  dat$id <- factor(dat$id)

  bench <- benchmark_standard_models(
    data = dat,
    formula = y ~ time_raw + gender + s(age, bs = "ps"),
    subject_var = "id",
    family = "gaussian",
    comparators = c("gee", "glmm")
  )

  expect_equal(bench$results$comparator, c("gee", "glmm"))
  expect_true(all(bench$results$available))
  expect_false(any(bench$results$success))
  expect_true(all(grepl("Smooth terms via s\\(\\.\\.\\.\\) are only supported", bench$results$error)))
})

test_that("benchmark_standard_models can score the supplied primary fit", {
  dat <- make_fixture_factor_time_interaction(n_subject = 10L)
  dat$id <- factor(dat$id)
  primary_fit <- stats::lm(y ~ time_raw + gender + age, data = dat)

  bench <- benchmark_standard_models(
    data = dat,
    formula = y ~ time_raw + gender + age,
    subject_var = "id",
    family = "gaussian",
    comparators = "gam",
    fit = primary_fit,
    fit_name = "primary"
  )

  expect_equal(bench$results$method, c("primary", "gam"))
  expect_equal(bench$results$comparator, c("primary", "gam"))
  expect_true(all(c("benchmark_mae", "benchmark_rmse") %in% names(bench$results)))
  expect_true(is.finite(bench$results$benchmark_rmse[bench$results$method == "primary"]))
  expect_true("primary" %in% names(bench$fits))
  expect_s3_class(bench$interpretation, "gamlss_longitudinal_benchmark_interpretation")
  expect_s3_class(bench$interpretation$leaders, "data.frame")
  expect_true(any(bench$interpretation$leaders$domain == "mean / marginal response fit"))
  expect_true(any(bench$interpretation$leaders$primary_leads_or_ties))
  expect_equal(bench$coefficients$parameter, "mu")
  expect_equal(bench$coefficients$reference_method, "primary")
  expect_s3_class(bench$coefficients$long, "data.frame")
  expect_s3_class(bench$coefficients$estimates, "data.frame")
  expect_s3_class(bench$coefficients$uncertainty, "data.frame")
  expect_true(all(c("term", "primary", "gam") %in% names(bench$coefficients$estimates)))
  expect_true(all(c(
    "estimate_diff_vs_reference",
    "se_ratio_vs_reference",
    "ci_width_ratio_vs_reference",
    "ci_overlap_with_reference"
  ) %in% names(bench$coefficients$uncertainty)))
  expect_true(any(bench$coefficients$long$term == "genderM"))

  summary <- summarise_benchmark_results(
    bench$results,
    metrics = c("benchmark_rmse", "elapsed_sec"),
    group_cols = NULL
  )
  expect_s3_class(summary, "gamlss_longitudinal_benchmark_summary")
  expect_true(any(summary$case_results$method == "primary"))

  printed <- utils::capture.output(print(bench))
  expect_true(any(grepl("Model Status", printed, fixed = TRUE)))
  expect_true(any(grepl("Mu Coefficients", printed, fixed = TRUE)))
  expect_true(any(grepl("Coefficient estimates", printed, fixed = TRUE)))
  expect_true(any(grepl("Uncertainty vs reference", printed, fixed = TRUE)))
  expect_true(any(grepl("Benchmark Metrics", printed, fixed = TRUE)))
  expect_true(any(grepl("Benchmark Summary", printed, fixed = TRUE)))
  expect_true(any(grepl("lead/tie", printed, fixed = TRUE)))
  expect_true(any(grepl("Observed response RMSE", printed, fixed = TRUE)))
})

test_that("benchmark interpretation groups metrics into readable domains", {
  domains <- .benchmark_metric_domain(
    c(
      "benchmark_mean_rmse",
      "benchmark_rmse",
      "benchmark_neg_log_score",
      "benchmark_q90_mae",
      "benchmark_interval_coverage_95",
      "benchmark_interval_width_95",
      "benchmark_theta_time_abs_error",
      "elapsed_sec"
    ),
    c("mean", "observed_response", "density", "quantile", "interval", "interval_width", "dependence", "runtime")
  )

  expect_equal(domains[1], "mean / marginal response fit")
  expect_equal(domains[2], "mean / marginal response fit")
  expect_equal(domains[3], "distributional prediction / shape")
  expect_equal(domains[4], "distributional prediction / shape")
  expect_equal(domains[5], "distributional prediction / shape")
  expect_equal(domains[6], "distributional prediction / shape")
  expect_equal(domains[7], "dependence")
  expect_equal(domains[8], "runtime")
})

test_that("benchmark interpretation excludes failed and unavailable rows from leaders", {
  results <- data.frame(
    method = c("gamlss.longitudinal", "glmmTMB"),
    comparator = c("gamlss.longitudinal", "glmmTMB"),
    available = c(TRUE, FALSE),
    success = c(TRUE, FALSE),
    benchmark_rmse = c(0.2, NA_real_),
    elapsed_sec = c(1.0, 0.001),
    error = c(NA_character_, "glmmTMB package is not installed"),
    stringsAsFactors = FALSE
  )

  interpretation <- .benchmark_interpretation(results, primary_method = "gamlss.longitudinal")

  expect_true(any(interpretation$unsupported$method == "glmmTMB"))
  expect_false(any(grepl("glmmTMB", interpretation$leaders$leading_methods, fixed = TRUE)))
  runtime <- interpretation$leaders[interpretation$leaders$metric == "elapsed_sec", , drop = FALSE]
  expect_true(runtime$primary_leads_or_ties)
  expect_false(runtime$standard_leads_or_ties)
})

test_that("benchmark_standard_models reports observed distribution diagnostics", {
  dat <- make_fixture_factor_time_interaction(n_subject = 10L)
  dat$id <- factor(dat$id)

  bench <- benchmark_standard_models(
    data = dat,
    formula = y ~ time_raw + gender + age,
    subject_var = "id",
    family = "gaussian",
    comparators = "gam"
  )

  metric_cols <- c(
    "benchmark_theta_time_abs_error",
    "benchmark_interval_coverage_95",
    "benchmark_interval_width_95",
    "benchmark_pit_ks_p_value",
    "benchmark_pit_mean_abs_error",
    "benchmark_tail_error_lower_05",
    "benchmark_tail_error_upper_05"
  )
  expect_true(all(metric_cols %in% names(bench$results)))
  observed_cols <- setdiff(metric_cols, "benchmark_theta_time_abs_error")
  expect_true(is.na(bench$results$benchmark_theta_time_abs_error))
  expect_true(all(is.finite(unlist(bench$results[observed_cols]))))
  expect_true(bench$results$benchmark_interval_coverage_95 >= 0)
  expect_true(bench$results$benchmark_interval_coverage_95 <= 1)
  expect_true(is.finite(bench$results$benchmark_interval_width_95))
})

test_that("benchmark_standard_models reports truth-aware quantile and tail metrics", {
  dat <- make_fixture_factor_time_interaction(n_subject = 12L)
  dat$id <- factor(dat$id)
  t_num <- as.integer(dat$time_raw)
  g_num <- ifelse(dat$gender == "M", 1, 0)
  dat$true_mu <- 1.5 + 0.25 * t_num + 0.4 * g_num + 0.3 * t_num * g_num + 0.01 * dat$age
  dat$true_sigma <- 0.15

  bench <- benchmark_standard_models(
    data = dat,
    formula = y ~ time_raw + gender + age,
    subject_var = "id",
    family = "gaussian",
    comparators = "gam",
    truth_family = "NO"
  )

  expect_true(all(c(
    "benchmark_mean_rmse",
    "benchmark_neg_log_score",
    "benchmark_q90_mae",
    "benchmark_interval_width_95",
    "benchmark_upper_tail_error_90"
  ) %in% names(bench$results)))
  expect_true(is.finite(bench$results$benchmark_mean_rmse))
  expect_true(is.finite(bench$results$benchmark_neg_log_score))
  expect_true(is.finite(bench$results$benchmark_q90_mae))
  expect_true(is.finite(bench$results$benchmark_interval_width_95))
  expect_true(is.finite(bench$results$benchmark_upper_tail_error_90))
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

test_that("summarise_benchmark_results recognises smooth recovery metrics", {
  results <- data.frame(
    family = "NO",
    copula = "N",
    design = "smooth",
    n_subject = 20L,
    n_time = 3L,
    dependence = "moderate",
    missingness = "none",
    start_mode = "default",
    method = c("rs_separate", "gam"),
    smooth_eta_rmse = c(0.12, NA_real_),
    smooth_eta_max_abs_error = c(0.20, NA_real_),
    stringsAsFactors = FALSE
  )

  out <- summarise_benchmark_results(
    results,
    metrics = c("smooth_eta_rmse", "smooth_eta_max_abs_error")
  )

  expect_true(all(c("smooth_eta_rmse", "smooth_eta_max_abs_error") %in% out$metric_catalog$metric))
  expect_true(all(out$metric_catalog$estimand == "smooth"))
  expect_true(any(out$case_results$method == "rs_separate"))
  expect_true(any(out$case_results$result == "missing"))
})
