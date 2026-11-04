test_that("coverage report LaTeX helpers escape special values", {
  escaped <- gamlss.longitudinal:::.coverage_tex_escape(c("a_b&c", NA, "x%y"))

  expect_equal(escaped[[1]], "a\\_b\\&c")
  expect_equal(escaped[[2]], "")
  expect_equal(escaped[[3]], "x\\%y")
  expect_equal(
    gamlss.longitudinal:::.coverage_tex_texttt("rs_separate"),
    "\\texttt{rs\\_separate}"
  )
})

test_that("coverage report numeric helper formats finite and missing values", {
  expect_equal(
    gamlss.longitudinal:::.coverage_format_number(c(1, 1.234, Inf, NA_real_), digits = 2),
    c("1.00", "1.23", "--", "--")
  )
})

test_that("coverage report table helpers render escaped tabular output", {
  df <- data.frame(
    Method = c("rs_separate", "gamlss2"),
    `A&B` = c("ok", "needs_review"),
    check.names = FALSE
  )

  table_lines <- gamlss.longitudinal:::.coverage_tex_table(df, "A_B caption")
  longtable_lines <- gamlss.longitudinal:::.coverage_tex_longtable(df, "Long_B caption")
  empty_lines <- gamlss.longitudinal:::.coverage_tex_table(df[0, ], "Empty_B caption")

  expect_true(any(grepl("\\\\caption\\{A\\\\_B caption\\}", table_lines)))
  expect_true(any(grepl("A\\\\&B", table_lines)))
  expect_true(any(grepl("rs\\\\_separate", table_lines)))
  expect_true(any(grepl("\\\\begin\\{longtable\\}", longtable_lines)))
  expect_true(any(grepl("Long\\\\_B caption continued", longtable_lines)))
  expect_equal(empty_lines, "\\paragraph{Empty\\_B caption} No rows available.")
})

test_that("coverage report fit and runtime table builders summarize rows", {
  results <- data.frame(
    family = c("NO", "NO", "PO", "PO"),
    copula = c("N", "C", "N", "C"),
    method = c("rs", "rs", "cg", "cg"),
    success = c(TRUE, FALSE, TRUE, TRUE),
    converged = c(TRUE, FALSE, TRUE, TRUE),
    elapsed_sec = c(1, 4, 2, 8),
    stringsAsFactors = FALSE
  )

  fit_tables <- gamlss.longitudinal:::.coverage_report_fit_tables(results)
  runtime_tables <- gamlss.longitudinal:::.coverage_report_runtime_tables(results, top_n = 2)

  expect_equal(fit_tables$fit_df$Method, c("cg", "rs"))
  expect_equal(fit_tables$fit_df$Successful, c(2L, 1L))
  expect_equal(fit_tables$fit_df$Failed, c(0L, 1L))
  expect_equal(nrow(fit_tables$by_copula), 4L)
  expect_equal(runtime_tables$runtime_df$Method, c("cg", "rs"))
  expect_equal(runtime_tables$runtime_df$`Total seconds`, c("10.00", "5.00"))
  expect_equal(nrow(runtime_tables$family_runtime_df), 2L)
  expect_equal(runtime_tables$slow_case_df$Seconds[[1]], "8.00")
})

test_that("coverage report eta and likelihood table builders isolate review rows", {
  parameter_results <- data.frame(
    family = c("NO", "PO", "NO"),
    copula = c("N", "C", "N"),
    method = c("rs", "rs", "gamlss"),
    parameter = c("mu", "sigma", "mu"),
    term = c("mu", "sigma", "mu"),
    estimate_eta = c(0.1, 1.2, 0),
    true_eta = c(0, 1, 0),
    abs_eta_error = c(0.1, 0.2, 0),
    eta_error_class = c("acceptable", "concern", "excellent"),
    stringsAsFactors = FALSE
  )
  results <- data.frame(
    family = c("NO", "PO"),
    copula = c("N", "C"),
    method = c("rs", "rs"),
    margin_review_class = c("excellent", "review"),
    joint_review_class = c("acceptable", "review"),
    margin_gap_pct_vs_reference = c(0.2, 5),
    joint_delta_pct_vs_rs_separate = c(0, -2),
    marginal_loglik = c(-10, -20),
    joint_loglik = c(-9, -18),
    elapsed_sec = c(1, 2),
    stringsAsFactors = FALSE
  )

  eta_tables <- gamlss.longitudinal:::.coverage_report_eta_tables(parameter_results, top_n = 2)
  likelihood_tables <- gamlss.longitudinal:::.coverage_report_likelihood_tables(results)

  expect_equal(eta_tables$eta_summary_df$Method, c("gamlss", "rs"))
  expect_equal(eta_tables$eta_class_df$concern[eta_tables$eta_class_df$Method == "rs"], 1L)
  expect_equal(nrow(eta_tables$eta_worst_df), 2L)
  expect_equal(eta_tables$eta_concern_df$Family, "PO")
  expect_true("review" %in% names(likelihood_tables$margin_df))
  expect_true("review" %in% names(likelihood_tables$joint_df))
  expect_equal(nrow(likelihood_tables$likelihood_review_df), 1L)
  expect_equal(likelihood_tables$likelihood_review_df$Family, "PO")
})

test_that("coverage report smooth table builder handles optional smooth results", {
  expect_null(gamlss.longitudinal:::.coverage_report_smooth_table(NULL))

  smooth_results <- data.frame(
    method = c("rs", "rs", "cg"),
    success = c(TRUE, FALSE, TRUE),
    smooth_eta_rmse = c(0.1, 0.3, 0.2),
    smooth_eta_max_abs_error = c(0.2, 0.5, 0.4),
    stringsAsFactors = FALSE
  )

  out <- gamlss.longitudinal:::.coverage_report_smooth_table(smooth_results)

  expect_equal(out$Method, c("cg", "rs"))
  expect_equal(out$`Successful fits`, c(1L, 1L))
  expect_equal(out$`Median eta RMSE`, c("0.200", "0.200"))
  expect_equal(out$`Max abs. eta error`, c("0.400", "0.500"))
})

test_that("coverage report workflow helper resolves prefix files and loads CSV input", {
  tmp <- tempfile("coverage-prefix-")
  results_file <- paste0(tmp, "_final_results.csv")
  parameter_file <- paste0(tmp, "_final_parameter_results.csv")
  results <- data.frame(
    family = "NO",
    copula = "N",
    design = "smoke",
    method = "rs",
    success = TRUE,
    elapsed_sec = 1,
    stringsAsFactors = FALSE
  )
  parameter_results <- data.frame(
    family = "NO",
    copula = "N",
    method = "rs",
    parameter = "mu",
    estimate_eta = 0,
    true_eta = 0,
    abs_eta_error = 0,
    eta_error_class = "excellent",
    stringsAsFactors = FALSE
  )
  utils::write.csv(results, results_file, row.names = FALSE)
  utils::write.csv(parameter_results, parameter_file, row.names = FALSE)

  inputs <- .coverage_report_prepare_inputs(
    prefix = tmp,
    results = NULL,
    parameter_results = NULL,
    smooth_results = NULL,
    results_file = NULL,
    parameter_results_file = NULL,
    smooth_results_file = NULL,
    output_tex = NULL
  )

  expect_equal(inputs$results_file, results_file)
  expect_equal(inputs$parameter_results_file, parameter_file)
  expect_equal(inputs$output_tex, paste0(tmp, "_summary.tex"))
  expect_equal(inputs$results$family, "NO")
  expect_equal(inputs$parameter_results$parameter, "mu")
})

test_that("coverage report document helper and writer produce reviewer tex", {
  results <- data.frame(
    family = c("NO", "PO"),
    copula = c("N", "C"),
    design = c("smoke", "smoke"),
    method = c("rs", "cg"),
    success = c(TRUE, FALSE),
    elapsed_sec = c(1, 2),
    margin_review_class = c("excellent", "review"),
    joint_review_class = c("excellent", "review"),
    margin_gap_pct_vs_reference = c(0, 2),
    joint_delta_pct_vs_rs_separate = c(0, -1),
    marginal_loglik = c(-10, -20),
    joint_loglik = c(-9, -18),
    stringsAsFactors = FALSE
  )
  tables <- .coverage_report_table_bundle(
    results = results,
    parameter_results = NULL,
    smooth_results = NULL,
    top_n = 1
  )

  lines <- .coverage_report_document_lines(
    results = results,
    tables = tables,
    title = "Coverage & Review",
    run_label = "unit-test",
    prefix = NULL,
    results_file = "results_file.csv",
    parameter_results_file = NULL
  )

  expect_true(any(grepl("\\\\title\\{Coverage \\\\& Review\\}", lines)))
  expect_true(any(grepl("\\\\date\\{unit-test\\}", lines)))
  expect_true(any(grepl("Total fitted cases: 2", lines)))

  output_tex <- tempfile(fileext = ".tex")
  out <- write_coverage_summary_report(
    results = results,
    output_tex = output_tex,
    compile_pdf = FALSE,
    title = "Coverage & Review",
    run_label = "unit-test",
    top_n = 1
  )

  expect_true(file.exists(output_tex))
  written <- readLines(output_tex, warn = FALSE)
  expect_true(any(grepl("\\\\section\\*\\{Scope\\}", written)))
  expect_false(out$compiled)
})
