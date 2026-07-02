test_that("benchmark report run summary and optional scenario section are stable", {
  inputs <- list(
    results = data.frame(method = c("a", "b")),
    reps = 2L,
    metrics = c("rmse", "coverage")
  )
  scenarios <- data.frame(
    scenario = "s1",
    label = "Scenario 1",
    family = "NO",
    copula = "N",
    design = "linear",
    n_subject = 10,
    n_time = 3,
    dependence = "moderate",
    missingness = "none",
    claim = "smoke",
    primary_metrics = I(list(c("rmse"))),
    methods = I(list(c("gamlss.longitudinal")))
  )

  summary_lines <- .benchmark_report_run_summary_section("Title", inputs)
  scenario_lines <- .benchmark_report_scenario_plan_section(scenarios)
  no_scenario_lines <- .benchmark_report_scenario_plan_section(NULL)

  expect_equal(summary_lines[1], "# Title")
  expect_true(any(grepl("Result rows: 2", summary_lines, fixed = TRUE)))
  expect_true(any(grepl("Metrics scored: rmse, coverage", summary_lines, fixed = TRUE)))
  expect_true(any(grepl("Scenario Plan", scenario_lines, fixed = TRUE)))
  expect_true(any(grepl("primary_metrics", scenario_lines, fixed = TRUE)))
  expect_identical(no_scenario_lines, character(0))
})

test_that("benchmark report table sections include their headings and selected columns", {
  comparator_status <- data.frame(
    comparator = "gam",
    comparator_class = "GAM",
    estimator = "mgcv::gam",
    package = "mgcv",
    available = TRUE,
    role = "smooth"
  )
  headlines <- data.frame(
    metric = "rmse",
    estimand = "mean",
    method = "gam",
    n = 2,
    n_finite = 2,
    wins = 1,
    ties = 0,
    losses = 1,
    missing = 0,
    win_or_tie_rate = 0.5,
    median_value = 1.25,
    median_score = -1.25
  )
  catalog <- data.frame(metric = "rmse", estimand = "mean", score_rule = "lower", target = "truth")

  comparator_lines <- .benchmark_report_comparator_section(comparator_status)
  headline_lines <- .benchmark_report_headline_section(headlines)
  metric_lines <- .benchmark_report_metric_definitions_section(catalog)

  expect_true(any(grepl("Comparator Availability", comparator_lines, fixed = TRUE)))
  expect_true(any(grepl("| comparator | comparator_class | estimator | package | available | role |", comparator_lines, fixed = TRUE)))
  expect_true(any(grepl("Headline Results", headline_lines, fixed = TRUE)))
  expect_true(any(grepl("| metric | estimand | method |", headline_lines, fixed = TRUE)))
  expect_true(any(grepl("Metric Definitions", metric_lines, fixed = TRUE)))
  expect_true(any(grepl("| metric | estimand | score_rule | target |", metric_lines, fixed = TRUE)))
})

test_that("benchmark report scenario and case-result sections are optional and bounded", {
  scenario_summary <- data.frame(
    benchmark_scenario = c("s1", "s1"),
    scenario_label = c("Scenario 1", "Scenario 1"),
    metric = c("rmse", "coverage"),
    estimand = c("mean", "interval"),
    method = c("gam", "gamlss.longitudinal"),
    n = c(2, 2),
    n_finite = c(2, 2),
    finite_rate = c(1, 1),
    wins = c(1, 0),
    ties = c(0, 1),
    losses = c(1, 1),
    missing = c(0, 0),
    win_or_tie_rate = c(0.5, 0.5),
    median_value = c(1.2, 0.9),
    median_score = c(-1.2, 0.9)
  )
  case_results <- data.frame(case_id = 1:3, method = c("a", "b", "c"))

  headline_lines <- .benchmark_report_scenario_headline_section(scenario_summary)
  summary_lines <- .benchmark_report_scenario_summary_section(scenario_summary)
  omitted_case_lines <- .benchmark_report_case_results_section(case_results, include_case_results = FALSE)
  included_case_lines <- .benchmark_report_case_results_section(
    case_results,
    include_case_results = TRUE,
    max_case_rows = 2L
  )

  expect_true(any(grepl("Scenario-Level Headline Results", headline_lines, fixed = TRUE)))
  expect_true(any(grepl("declared primary metrics", headline_lines, fixed = TRUE)))
  expect_true(any(grepl("Scenario-Level Win/Tie/Loss Summary", summary_lines, fixed = TRUE)))
  expect_identical(omitted_case_lines, character(0))
  expect_true(any(grepl("Per-Case Results", included_case_lines, fixed = TRUE)))
  expect_false(any(grepl("| 3 | c |", included_case_lines, fixed = TRUE)))
})

test_that("benchmark report headline helpers select best methods deterministically", {
  summary <- data.frame(
    metric = c("rmse", "rmse", "coverage", "coverage"),
    estimand = c("mean", "mean", "interval", "interval"),
    method = c("b", "a", "b", "a"),
    win_or_tie_rate = c(0.8, 0.8, 0.5, 0.7),
    median_score = c(1.2, 1.1, 0.4, 0.5),
    n_finite = c(10L, 10L, 8L, 8L),
    stringsAsFactors = FALSE
  )

  headlines <- .benchmark_report_headlines(summary)

  expect_equal(headlines$method[headlines$metric == "rmse"], "a")
  expect_equal(headlines$method[headlines$metric == "coverage"], "a")
  expect_equal(nrow(.benchmark_report_headlines(data.frame(metric = "rmse"))), 0L)
})

test_that("benchmark report group helpers summarize and rank scenario results", {
  case_results <- data.frame(
    benchmark_scenario = c("s1", "s1", "s1", "s2"),
    metric = c("rmse", "rmse", "rmse", "rmse"),
    estimand = c("mean", "mean", "mean", "mean"),
    method = c("a", "a", "b", "a"),
    value = c(1, 3, Inf, 2),
    score = c(1, 3, NA, 2),
    result = c("win", "tie", "missing", "loss"),
    stringsAsFactors = FALSE
  )

  group_summary <- .benchmark_report_group_summary(case_results, "benchmark_scenario")
  group_headlines <- .benchmark_report_group_headlines(group_summary, "benchmark_scenario")

  s1_a <- group_summary[group_summary$benchmark_scenario == "s1" & group_summary$method == "a", ]

  expect_equal(s1_a$n, 2L)
  expect_equal(s1_a$n_finite, 2L)
  expect_equal(s1_a$win_or_tie_rate, 1)
  expect_equal(s1_a$median_value, 2)
  expect_true(all(c("s1", "s2") %in% group_headlines$benchmark_scenario))
})

test_that("benchmark report metadata helpers filter scenarios and infer primary method", {
  scenario_summary <- data.frame(
    benchmark_scenario = c("s1", "s1", "s2", "unknown"),
    metric = c("rmse", "coverage", "elapsed_sec", "anything"),
    stringsAsFactors = FALSE
  )
  scenarios <- data.frame(
    scenario = c("s1", "s2"),
    label = c("Scenario 1", "Scenario 2"),
    primary_metrics = I(list("rmse", character(0))),
    stringsAsFactors = FALSE
  )

  filtered <- .benchmark_report_filter_scenario_metrics(scenario_summary, scenarios)
  domains <- .benchmark_metric_domain(
    metric = c("benchmark_rmse", "benchmark_neg_log_score", "elapsed_sec", "smooth_mu_edf", "other"),
    estimand = c("mean", "density", "runtime", "smooth", "custom")
  )

  expect_equal(filtered$metric, c("rmse", "elapsed_sec", "anything"))
  expect_equal(filtered$scenario_label, c("Scenario 1", "Scenario 2", NA))
  expect_equal(domains, c(
    "mean / marginal response fit",
    "distributional prediction / shape",
    "runtime",
    "smooth recovery",
    "other"
  ))
  expect_identical(.benchmark_infer_primary_method(data.frame(method = c("gam", "gamlss.longitudinal"))), "gamlss.longitudinal")
  expect_identical(.benchmark_infer_primary_method(data.frame(method = c("gam", "rs_separate"))), "rs_separate")
  expect_null(.benchmark_infer_primary_method(data.frame(method = "gam")))
})

test_that("benchmark interpretation notes remain reviewer-facing", {
  notes <- .benchmark_report_interpretation_notes_section()

  expect_true(any(grepl("Interpretation Notes", notes, fixed = TRUE)))
  expect_true(any(grepl("Benchmark claims should be reported by estimand", notes, fixed = TRUE)))
})

test_that("benchmark interpretation formatter reports leaders and unsupported rows", {
  interpretation <- list(
    primary_method = "gamlss.longitudinal",
    leaders = data.frame(
      domain = c("mean / marginal response fit", "runtime"),
      metric = c("benchmark_rmse", "elapsed_sec"),
      estimand = c("mean", "runtime"),
      primary_leads_or_ties = c(TRUE, FALSE),
      standard_leads_or_ties = c(FALSE, TRUE),
      primary_methods = c("gamlss.longitudinal", ""),
      standard_methods = c("", "gam"),
      leading_methods = c("gamlss.longitudinal", "gam"),
      stringsAsFactors = FALSE
    ),
    unsupported = data.frame(
      method = "glmmTMB",
      comparator = "glmmTMB",
      reason = "unavailable",
      error = "package not installed",
      stringsAsFactors = FALSE
    )
  )
  class(interpretation) <- "gamlss_longitudinal_benchmark_interpretation"

  lines <- .benchmark_interpretation_lines(interpretation)
  table <- .benchmark_interpretation_table(interpretation)

  expect_true(any(grepl("gamlss.longitudinal leads or ties on: benchmark_rmse", lines, fixed = TRUE)))
  expect_true(any(grepl("Standard comparators lead or tie on: elapsed_sec", lines, fixed = TRUE)))
  expect_true(any(grepl("glmmTMB [unavailable]", lines, fixed = TRUE)))
  expect_equal(table$area, c("Mean fit", "Runtime"))
  expect_equal(table$metric, c("benchmark_rmse", "elapsed_sec"))
  expect_equal(table$best_model, c("gamlss.longitudinal", "gam"))
})
