test_that("coverage family catalog separates supported and default-unsupported families", {
  skip_if_not_installed("gamlss.dist")

  catalog <- gamlss.longitudinal:::.coverage_family_catalog()

  expect_true(all(c("family", "type", "parameters", "supported", "unsupported_reason") %in% names(catalog)))
  expect_true(all(c("NO", "PO", "ZIP") %in% catalog$family))
  expect_true(catalog$supported[match("NO", catalog$family)])
  expect_true(catalog$supported[match("PO", catalog$family)])
  expect_true(catalog$supported[match("ZIP", catalog$family)])
  expect_false(catalog$supported[match("BI", catalog$family)])
  expect_false(catalog$supported[match("BB", catalog$family)])
  expect_false(catalog$supported[match("DBI", catalog$family)])
  expect_false(catalog$supported[match("LG", catalog$family)])
  expect_false(catalog$supported[match("MN3", catalog$family)])
  expect_match(catalog$unsupported_reason[match("BI", catalog$family)], "denominator")
  expect_match(catalog$unsupported_reason[match("BB", catalog$family)], "denominator")
  expect_match(catalog$unsupported_reason[match("LG", catalog$family)], "starting")
  expect_match(catalog$unsupported_reason[match("MN3", catalog$family)], "ordinal")

  mixed <- catalog[grepl("Mixed", catalog$type), , drop = FALSE]
  if (nrow(mixed) > 0L) {
    expect_false(any(mixed$supported))
    expect_true(all(nzchar(mixed$unsupported_reason)))
  }
})

test_that("coverage safe defaults keep constrained families inside support", {
  skip_if_not_installed("gamlss.dist")

  bi_params <- gamlss.longitudinal:::.coverage_default_margin_params(gamlss.dist::BI())
  lg_params <- gamlss.longitudinal:::.coverage_default_margin_params(gamlss.dist::LG())

  expect_equal(bi_params$mu, 0.5)
  expect_equal(lg_params$mu, 0.5)
  expect_true(is.finite(gamlss.dist::dBI(0, mu = bi_params$mu)))
  expect_true(is.finite(gamlss.dist::dLG(1, mu = lg_params$mu)))
})

test_that("coverage harness fits representative families and copulas", {
  skip_if_not_installed("gamlss")
  skip_if_not_installed("gamlss.dist")

  grid <- data.frame(
    family = c("NO", "GA", "BCPE", "PO", "NBI", "DEL", "ZIP", "ZAP", "ZINBI", "NO", "NO", "NO", "NO", "NO", "PO"),
    copula = c(rep("N", 9), "C", "F", "G", "J", "t", "C"),
    method = c(rep("rs_separate", 9), rep("rs_joint", 6)),
    design = "intercept",
    case_id = seq_len(15L),
    stringsAsFactors = FALSE
  )

  results <- gamlss.longitudinal:::.coverage_run_grid(
    grid,
    n = 36L,
    times = 1:3,
    seed = 930,
    max_outer_iter = 3L,
    max_inner_iter = 3L,
    max_elapsed_sec = 15
  )

  expect_equal(nrow(results), nrow(grid))
  expect_true(all(results$success))
  expect_true(all(is.finite(results$marginal_loglik)))
  expect_true(all(is.finite(results$joint_loglik)))
  expect_true(all(is.finite(results$elapsed_sec)))
  expect_true(all(is.na(results$true_copula_tau) | abs(results$fitted_copula_tau - results$true_copula_tau) < 0.95))
})

test_that("coverage harness records all fit methods on a small intercept case", {
  skip_if_not_installed("gamlss")
  skip_if_not_installed("gamlss.dist")

  grid <- gamlss.longitudinal:::.coverage_make_case_grid(
    families = "NO",
    copulas = "N",
    methods = c("rs_separate", "rs_joint"),
    designs = "intercept"
  )

  results <- gamlss.longitudinal:::.coverage_run_grid(
    grid,
    n = 40L,
    times = 1:3,
    seed = 940,
    max_outer_iter = 3L,
    max_inner_iter = 3L,
    max_elapsed_sec = 15
  )

  expect_equal(sort(results$method), sort(c("rs_separate", "rs_joint")))
  expect_true(all(results$success))
  expect_true(all(is.finite(results$marginal_loglik)))
  expect_true(all(results$failure_type == "ok"))
  expect_true(all(c(
    "margin_gap_pct_vs_gamlss2",
    "margin_review_class",
    "joint_delta_pct_vs_rs_separate",
    "joint_review_class",
    "elapsed_sec"
  ) %in% names(results)))
  expect_true(all(is.finite(results$elapsed_sec)))
  expect_true(all(results$joint_review_class != "missing"))
})

test_that("coverage harness includes CRAN-only gamlss baseline", {
  skip_if_not_installed("gamlss")
  skip_if_not_installed("gamlss.dist")

  grid <- gamlss.longitudinal:::.coverage_make_case_grid(
    families = "NO",
    copulas = "N",
    methods = c("gamlss", "rs_separate"),
    designs = "intercept"
  )

  results <- gamlss.longitudinal:::.coverage_run_grid(
    grid,
    n = 30L,
    times = 1:3,
    seed = 942,
    max_outer_iter = 2L,
    max_inner_iter = 3L,
    max_elapsed_sec = 15
  )

  expect_equal(sort(results$method), sort(c("gamlss", "rs_separate")))
  expect_true(all(results$success))
  expect_true(all(is.finite(results$marginal_loglik)))
  expect_true("reference_marginal_method" %in% names(results))
  expect_true(any(results$reference_marginal_method == "gamlss", na.rm = TRUE))
})

test_that("coverage gamlss2 baseline records fallback method when needed", {
  skip_if_not_installed("gamlss")
  skip_if_not_installed("gamlss.dist")

  grid <- gamlss.longitudinal:::.coverage_make_case_grid(
    families = "PIG2",
    copulas = "N",
    methods = "gamlss2",
    designs = "intercept"
  )

  results <- gamlss.longitudinal:::.coverage_run_grid(
    grid,
    n = 24L,
    times = 1:3,
    seed = 943,
    max_outer_iter = 2L,
    max_inner_iter = 2L,
    max_elapsed_sec = 15
  )

  expect_true(results$success)
  expect_true(results$marginal_fit_method %in% c("gamlss2", "gamlss"))
  expect_equal(results$reference_marginal_method, results$marginal_fit_method)
  if (isTRUE(results$marginal_fallback_used)) {
    expect_equal(results$marginal_fit_method, "gamlss")
    expect_true(nzchar(results$marginal_fallback_error))
    expect_equal(results$fit_attempt_trace, "gamlss2 > gamlss")
  }
})

test_that("coverage starts include exactly each margin parameter plus copula theta", {
  skip_if_not_installed("gamlss")
  skip_if_not_installed("gamlss.dist")

  for (family in c("LNO", "NET")) {
    sim <- gamlss.longitudinal:::.coverage_simulate_case(
      family = family,
      copula = "N",
      design = "intercept",
      n = 18L,
      times = 1:3,
      seed = 944,
      dependence = "moderate"
    )
    margin_dist <- gamlss.longitudinal:::.normalise_margin_dist_links(
      do.call(get(family, envir = asNamespace("gamlss.dist")), list())
    )
    starts <- suppressWarnings(suppressMessages(
      gamlss.longitudinal:::get_starting_values("N", margin_dist, sim, eta_transform = TRUE)
    ))

    expect_named(starts, c(names(margin_dist$parameters), "theta"), ignore.order = FALSE)
    expect_true(all(is.finite(starts)))
  }
})

test_that("coverage domain-sensitive defaults generate finite q, p, and d values", {
  skip_if_not_installed("gamlss.dist")

  for (family in c("BNB", "RGE", "SI", "ZABNB", "ZIBNB", "ZINBF")) {
    margin_dist <- do.call(get(family, envir = asNamespace("gamlss.dist")), list())
    params <- gamlss.longitudinal:::.coverage_default_margin_params(margin_dist)
    qfun <- get(paste0("q", family), envir = asNamespace("gamlss.dist"))
    pfun <- get(paste0("p", family), envir = asNamespace("gamlss.dist"))
    dfun <- get(paste0("d", family), envir = asNamespace("gamlss.dist"))

    y <- do.call(qfun, c(list(p = 0.5), params))
    p <- do.call(pfun, c(list(q = y), params))
    d <- do.call(dfun, c(list(x = y), params))

    expect_true(is.finite(y), info = family)
    expect_true(is.finite(p), info = family)
    expect_true(is.finite(d), info = family)
  }
})

test_that("discrete lower-bound CDFs short-circuit below support", {
  skip_if_not_installed("gamlss.dist")

  margin_dist <- gamlss.dist::ZINBF()
  response <- c(-1, 0, 2)
  eta_inv <- list(
    mu = rep(4, 3),
    sigma = rep(0.5, 3),
    nu = rep(2, 3),
    tau = rep(0.2, 3)
  )
  mm <- eta_inv

  out <- gamlss.longitudinal:::calc_F_x(eta_inv, mm, margin_dist, response)
  expect_equal(out[[1]], 0)
  expect_true(all(is.finite(out[2:3])))
})

test_that("coverage CG records fallback attempt metadata", {
  skip_if_not_installed("gamlss")
  skip_if_not_installed("gamlss.dist")

  grid <- gamlss.longitudinal:::.coverage_make_case_grid(
    families = "NO",
    copulas = "N",
    methods = "cg",
    designs = "intercept"
  )

  results <- gamlss.longitudinal:::.coverage_run_grid(
    grid,
    n = 22L,
    times = 1:3,
    seed = 946,
    max_outer_iter = 2L,
    max_inner_iter = 2L,
    max_elapsed_sec = 15
  )

  expect_true(all(c("fit_attempt", "fit_attempt_count", "fit_attempt_trace") %in% names(results)))
  expect_true(results$fit_attempt_count >= 1L)
  expect_true(nzchar(results$fit_attempt_trace))
  expect_true(is.finite(results$elapsed_sec))
})

test_that("coverage harness includes a tiny fixed-covariate RS subset", {
  skip_if_not_installed("gamlss")
  skip_if_not_installed("gamlss.dist")

  grid <- gamlss.longitudinal:::.coverage_make_case_grid(
    families = "NO",
    copulas = "N",
    methods = c("rs_separate", "rs_joint"),
    designs = "covariate"
  )

  results <- gamlss.longitudinal:::.coverage_run_grid(
    grid,
    n = 28L,
    times = 1:3,
    seed = 945,
    max_outer_iter = 3L,
    max_inner_iter = 3L,
    max_elapsed_sec = 15
  )

  expect_equal(sort(results$method), sort(c("rs_separate", "rs_joint")))
  expect_true(all(results$success))
  expect_true(all(is.finite(results$marginal_loglik)))
})

test_that("coverage harness fits smooths on every active parameter for representative families", {
  skip_on_cran()
  skip_if_not_installed("gamlss")
  skip_if_not_installed("gamlss.dist")

  grid <- gamlss.longitudinal:::.coverage_make_case_grid(
    families = c("NO", "GA", "PO", "ZIP", "BCPE"),
    copulas = "N",
    methods = c("rs_joint", "cg"),
    designs = "smooth"
  )

  results <- gamlss.longitudinal:::.coverage_run_grid(
    grid,
    n = 24L,
    times = 1:3,
    seed = 20260902,
    max_outer_iter = 2L,
    max_inner_iter = 2L,
    max_elapsed_sec = 30,
    method_max_outer_iter = list(cg = 2L),
    method_max_inner_iter = list(cg = 2L)
  )

  expect_true(all(results$success))
  expect_true(all(is.finite(results$smooth_eta_rmse)))
  expect_true(all(results$smooth_eta_rmse < 0.9))
  expect_true(all(results$smooth_eta_n > 0))
})

test_that("coverage harness supports scale-varying benchmark designs", {
  skip_if_not_installed("gamlss")
  skip_if_not_installed("gamlss.dist")

  grid <- gamlss.longitudinal:::.coverage_make_case_grid(
    families = "NO",
    copulas = "N",
    methods = c("rs_separate", "gam"),
    designs = "scale"
  )

  results <- gamlss.longitudinal:::.coverage_run_grid(
    grid,
    n = 24L,
    times = 1:3,
    seed = 947,
    max_outer_iter = 2L,
    max_inner_iter = 2L,
    max_elapsed_sec = 15
  )

  expect_equal(unique(results$design), "scale")
  expect_true(all(c("benchmark_interval_coverage_95", "benchmark_interval_width_95", "benchmark_tail_error_upper_05") %in% names(results)))
  parameter_results <- attr(results, "parameter_results")
  expect_s3_class(parameter_results, "data.frame")
  scale_truth <- parameter_results[
    parameter_results$method == "rs_separate" &
      parameter_results$parameter == "sigma" &
      parameter_results$term == "x",
    ,
    drop = FALSE
  ]
  expect_equal(nrow(scale_truth), 1L)
  expect_true(is.finite(scale_truth$true_eta))
})

test_that("coverage harness supports time-varying dependence designs", {
  skip_if_not_installed("gamlss")
  skip_if_not_installed("gamlss.dist")

  grid <- gamlss.longitudinal:::.coverage_make_case_grid(
    families = "NO",
    copulas = "N",
    methods = c("rs_separate", "gam"),
    designs = "time_dependence"
  )

  results <- gamlss.longitudinal:::.coverage_run_grid(
    grid,
    n = 28L,
    times = 1:4,
    seed = 948,
    max_outer_iter = 4L,
    max_inner_iter = 3L,
    max_elapsed_sec = 20
  )

  expect_equal(unique(results$design), "time_dependence")
  expect_true("benchmark_theta_time_abs_error" %in% names(results))
  longitudinal_row <- results[results$method == "rs_separate", , drop = FALSE]
  expect_true(longitudinal_row$success)
  expect_true(is.finite(longitudinal_row$benchmark_theta_time_abs_error))

  parameter_results <- attr(results, "parameter_results")
  theta_time <- parameter_results[
      parameter_results$method == "rs_separate" &
      parameter_results$parameter == "theta" &
      parameter_results$term == "time_covariate",
    ,
    drop = FALSE
  ]
  expect_equal(nrow(theta_time), 1L)
  expect_true(is.finite(theta_time$true_eta))
})

test_that("coverage harness records missingness and stress scenarios", {
  skip_if_not_installed("gamlss")
  skip_if_not_installed("gamlss.dist")

  grid <- gamlss.longitudinal:::.coverage_make_case_grid(
    families = "NO",
    copulas = "N",
    methods = "rs_separate",
    designs = "intercept"
  )

  missing_results <- do.call(rbind, lapply(c("mcar", "drop_rows"), function(missingness) {
    gamlss.longitudinal:::.coverage_run_grid(
      grid,
      n = 30L,
      times = 1:3,
      seed = 960,
      max_outer_iter = 3L,
      max_inner_iter = 3L,
      max_elapsed_sec = 15,
      missingness = missingness
    )
  }))

  expect_equal(sort(unique(missing_results$missingness)), c("drop_rows", "mcar"))
  expect_true(all(missing_results$success))
  expect_true(all(is.finite(missing_results$joint_loglik)))

  stress_results <- do.call(rbind, lapply(c("near_independent", "strong"), function(dependence) {
    gamlss.longitudinal:::.coverage_run_grid(
      grid,
      n = 30L,
      times = 1:3,
      seed = 970,
      max_outer_iter = 3L,
      max_inner_iter = 3L,
      max_elapsed_sec = 15,
      dependence = dependence
    )
  }))

  expect_equal(sort(unique(stress_results$dependence)), c("near_independent", "strong"))
  expect_true(all(stress_results$success))
  expect_true(all(is.finite(stress_results$true_copula_tau)))
  expect_true(all(is.finite(stress_results$fitted_copula_tau)))
})

test_that("coverage harness supports truth-adjacent starting values", {
  skip_if_not_installed("gamlss")
  skip_if_not_installed("gamlss.dist")

  grid <- gamlss.longitudinal:::.coverage_make_case_grid(
    families = "NO",
    copulas = "N",
    methods = "rs_separate",
    designs = "intercept"
  )

  results <- gamlss.longitudinal:::.coverage_run_grid(
    grid,
    n = 30L,
    times = 1:3,
    seed = 980,
    max_outer_iter = 3L,
    max_inner_iter = 3L,
    max_elapsed_sec = 15,
    start_mode = "truth_adjacent"
  )

  expect_equal(results$start_mode, "truth_adjacent")
  expect_true(results$success)
  expect_true(is.finite(results$marginal_loglik))
})

test_that("coverage harness can run an opt-in non-writing result set", {
  skip_if_not_installed("gamlss")
  skip_if_not_installed("gamlss.dist")

  results <- run_coverage_simulations(
    families = "NO",
    copulas = "N",
    methods = c("rs_separate", "rs_joint"),
    designs = "intercept",
    n = 20L,
    times = 1:3,
    seed = 950,
    max_outer_iter = 2L,
    max_inner_iter = 2L,
    max_elapsed_sec = 15,
    write_results = FALSE
  )

  expect_equal(nrow(results), 2L)
  expect_true(all(results$success))
  parameter_results <- attr(results, "parameter_results")
  runtime_summary <- attr(results, "runtime_summary")
  expect_s3_class(parameter_results, "data.frame")
  expect_s3_class(runtime_summary, "data.frame")
  expect_true(all(c("estimate_eta", "true_eta", "abs_eta_error", "eta_error_class") %in% names(parameter_results)))
  expect_true(all(c("elapsed_sec", "family", "method") %in% names(runtime_summary)))
})

test_that("coverage summary report is generated from multi-copula result rows", {
  results <- data.frame(
    family = rep(c("NO", "PO"), each = 4),
    copula = rep(c("N", "C"), each = 2, times = 2),
    design = "intercept",
    method = rep(c("gamlss2", "rs_separate"), times = 4),
    success = TRUE,
    converged = TRUE,
    failure_type = "ok",
    elapsed_sec = c(0.1, 0.5, 0.2, 0.8, 0.1, 0.6, 0.2, 0.9),
    marginal_loglik = -10,
    joint_loglik = -9,
    margin_gap_pct_vs_reference = c(NA, 0.5, NA, 10, NA, 0.7, NA, 1.2),
    joint_delta_pct_vs_rs_separate = c(NA, 0, NA, -3, NA, 0, NA, -0.5),
    margin_review_class = c("reference", "excellent", "reference", "review", "reference", "excellent", "reference", "excellent"),
    joint_review_class = c("not_applicable", "acceptable", "not_applicable", "review", "not_applicable", "acceptable", "not_applicable", "acceptable"),
    stringsAsFactors = FALSE
  )
  parameter_results <- data.frame(
    family = c("NO", "PO"),
    copula = c("N", "C"),
    method = c("rs_separate", "rs_separate"),
    parameter = c("mu", "mu"),
    true_eta = c(0, 1),
    estimate_eta = c(0.1, 1.2),
    abs_eta_error = c(0.1, 0.2),
    eta_error_class = c("acceptable", "concern"),
    stringsAsFactors = FALSE
  )
  out_tex <- file.path(tempdir(), "coverage-summary-test.tex")

  report <- write_coverage_summary_report(
    results = results,
    parameter_results = parameter_results,
    output_tex = out_tex,
    compile_pdf = FALSE,
    run_label = "unit-test"
  )

  expect_true(file.exists(out_tex))
  expect_false(report$compiled)
  txt <- readLines(out_tex, warn = FALSE)
  expect_true(any(grepl("Fit success by copula and method", txt, fixed = TRUE)))
  expect_true(any(grepl("N", txt, fixed = TRUE)))
  expect_true(any(grepl("C", txt, fixed = TRUE)))
  expect_true(any(grepl("Eta-scale absolute error summary", txt, fixed = TRUE)))
  expect_true(any(grepl("Detailed Review Rows", txt, fixed = TRUE)))
  expect_true(any(grepl("All parameter recovery concern rows", txt, fixed = TRUE)))
  expect_true(any(grepl("All likelihood review rows", txt, fixed = TRUE)))
})

test_that("coverage harness can include standard GEE/GLMM/GAM comparator rows", {
  skip_if_not_installed("gamlss.dist")
  skip_if_not_installed("geepack")
  skip_if_not_installed("lme4")

  grid <- gamlss.longitudinal:::.coverage_make_case_grid(
    families = "NO",
    copulas = "N",
    methods = c("rs_separate", "gee", "glmm", "gam"),
    designs = "covariate"
  )

  results <- gamlss.longitudinal:::.coverage_run_grid(
    grid,
    n = 24L,
    times = 1:3,
    seed = 990,
    max_outer_iter = 2L,
    max_inner_iter = 2L,
    max_elapsed_sec = 15
  )

  expect_equal(sort(results$method), sort(c("rs_separate", "gee", "glmm", "gam")))
  longitudinal_row <- results[results$method == "rs_separate", , drop = FALSE]
  expect_true(longitudinal_row$success)
  expect_equal(longitudinal_row$benchmark_class, "gamlss_longitudinal")
  expect_true(is.finite(longitudinal_row$benchmark_rmse))
  expect_true(is.finite(longitudinal_row$benchmark_mean_rmse))
  expect_true(is.finite(longitudinal_row$benchmark_neg_log_score))
  expect_true(is.finite(longitudinal_row$benchmark_q90_mae))
  expect_true(is.finite(longitudinal_row$benchmark_upper_tail_error_90))
  expect_true(is.finite(longitudinal_row$benchmark_interval_coverage_95))
  expect_true(is.finite(longitudinal_row$benchmark_interval_width_95))
  comparator_rows <- results[results$method %in% c("gee", "glmm", "gam"), , drop = FALSE]
  expect_true(all(comparator_rows$success))
  expect_true(all(is.finite(comparator_rows$benchmark_rmse)))
  expect_true(all(c(
    "benchmark_comparator",
    "benchmark_estimator",
    "benchmark_mae",
    "benchmark_rmse",
    "benchmark_mean_bias",
    "benchmark_mean_mae",
    "benchmark_mean_rmse",
    "benchmark_neg_log_score",
    "benchmark_q90_mae",
    "benchmark_upper_tail_error_90",
    "benchmark_interval_coverage_95",
    "benchmark_interval_width_95",
    "benchmark_pit_ks_p_value",
    "benchmark_tail_error_lower_05",
    "benchmark_tail_error_upper_05"
  ) %in% names(results)))
  expect_true(all(is.finite(comparator_rows$benchmark_mean_rmse)))
  expect_true(all(is.finite(comparator_rows$benchmark_neg_log_score)))
  expect_true(all(is.finite(comparator_rows$benchmark_q90_mae)))
  expect_true(all(is.finite(comparator_rows$benchmark_upper_tail_error_90)))
  expect_true(all(is.finite(comparator_rows$benchmark_interval_coverage_95)))
  expect_true(all(is.finite(comparator_rows$benchmark_interval_width_95)))
  expect_true(all(comparator_rows$benchmark_interval_coverage_95 >= 0 & comparator_rows$benchmark_interval_coverage_95 <= 1))

  summary <- summarise_benchmark_results(
    results,
    metrics = c("benchmark_mean_rmse", "benchmark_neg_log_score", "benchmark_interval_coverage_95", "benchmark_interval_width_95", "elapsed_sec")
  )
  expect_s3_class(summary, "gamlss_longitudinal_benchmark_summary")
  expect_true(all(c("summary", "case_results", "metric_catalog") %in% names(summary)))
  expect_true(any(summary$summary$metric == "benchmark_mean_rmse"))
  expect_true(any(summary$case_results$result %in% c("win", "tie", "loss")))
})

test_that("coverage benchmark metrics include non-Gaussian q90 and tail checks", {
  skip_if_not_installed("gamlss.dist")

  grid <- gamlss.longitudinal:::.coverage_make_case_grid(
    families = c("GA", "PO"),
    copulas = "N",
    methods = "gam",
    designs = "covariate"
  )

  results <- gamlss.longitudinal:::.coverage_run_grid(
    grid,
    n = 24L,
    times = 1:3,
    seed = 992,
    max_elapsed_sec = 15
  )

  expect_equal(sort(unique(results$family)), c("GA", "PO"))
  expect_true(all(results$success))
  expect_true(all(is.finite(results$benchmark_mean_rmse)))
  expect_true(all(is.finite(results$benchmark_neg_log_score)))
  expect_true(all(is.finite(results$benchmark_q90_mae)))
  expect_true(all(is.finite(results$benchmark_upper_tail_error_90)))
  expect_true(all(is.na(results$benchmark_interval_coverage_95)))

  summary <- summarise_benchmark_results(
    results,
    metrics = c("benchmark_neg_log_score", "benchmark_q90_mae", "benchmark_upper_tail_error_90")
  )
  expect_true(all(c("benchmark_neg_log_score", "benchmark_q90_mae", "benchmark_upper_tail_error_90") %in% summary$metric_catalog$metric))
})

test_that("coverage harness records unsupported standard comparator families", {
  skip_if_not_installed("gamlss.dist")

  grid <- gamlss.longitudinal:::.coverage_make_case_grid(
    families = "BCPE",
    copulas = "N",
    methods = "gee",
    designs = "intercept"
  )

  results <- gamlss.longitudinal:::.coverage_run_grid(
    grid,
    n = 12L,
    times = 1:3,
    seed = 991,
    max_elapsed_sec = 15
  )

  expect_equal(results$method, "gee")
  expect_false(results$success)
  expect_equal(results$failure_type, "unsupported comparator family")
  expect_match(results$benchmark_error, "No standard comparator family mapping")
  expect_true("benchmark_mean_rmse" %in% names(results))
  expect_true(is.na(results$benchmark_mean_rmse))
  expect_true("benchmark_q90_mae" %in% names(results))
  expect_true(is.na(results$benchmark_q90_mae))
})

test_that("package ships CRAN-safe JSS replication entry point", {
  script <- system.file("jss-replication", "run-replication.R", package = "gamlss.longitudinal")
  expect_true(file.exists(script))
  txt <- readLines(script, warn = FALSE)
  expect_true(any(grepl("run_coverage_simulations", txt, fixed = TRUE)))
  expect_true(any(grepl("\"gamlss\"", txt, fixed = TRUE)))
})
