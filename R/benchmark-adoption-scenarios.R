#' Adoption benchmark scenario grid
#'
#' Defines a small opt-in scenario plan for reviewer/adoption benchmarks. These
#' scenarios are wrappers around the coverage simulation harness, not CRAN-time
#' tests.
#'
#' @param profile Scenario profile. `"smoke"` is intentionally small; `"review"`
#'   adds a few more representative continuous/count and missingness cases.
#' @param include_mixed Include mixed-family coverage cases when supported by
#'   the coverage harness.
#'
#' @return A data frame with one row per benchmark scenario.
#' @export
adoption_benchmark_scenarios <- function(profile = c("smoke", "review"), include_mixed = FALSE) {
  profile <- match.arg(profile)

  scenarios <- data.frame(
    scenario = c("continuous_gaussian", "count_poisson"),
    label = c("Continuous Gaussian baseline", "Count Poisson baseline"),
    family = c("NO", "PO"),
    copula = c("N", "C"),
    design = c("covariate", "covariate"),
    n_subject = c(40L, 40L),
    n_time = c(3L, 3L),
    dependence = c("moderate", "moderate"),
    missingness = c("none", "none"),
    start_mode = c("default", "default"),
    claim = c("continuous-response adoption smoke", "count-response adoption smoke"),
    stringsAsFactors = FALSE
  )

  if (identical(profile, "review")) {
    scenarios <- rbind(
      scenarios,
      data.frame(
        scenario = c("positive_gamma", "dropout_sensitivity"),
        label = c("Positive continuous Gamma", "Dropout sensitivity"),
        family = c("GA", "NO"),
        copula = c("N", "N"),
        design = c("scale", "covariate"),
        n_subject = c(60L, 60L),
        n_time = c(4L, 4L),
        dependence = c("moderate", "moderate"),
        missingness = c("none", "dropout"),
        start_mode = c("default", "truth_adjacent"),
        claim = c("positive continuous adoption case", "missingness/dropout stress case"),
        stringsAsFactors = FALSE
      )
    )
  }

  scenarios$methods <- I(rep(list(c("gamlss", "rs_separate", "rs_joint", "cg")), nrow(scenarios)))
  scenarios$primary_metrics <- I(rep(list(c("benchmark_rmse", "benchmark_mae", "elapsed_sec")), nrow(scenarios)))
  scenarios$include_mixed <- isTRUE(include_mixed)

  scenarios
}

#' Run adoption benchmark scenarios
#'
#' Runs opt-in adoption benchmark scenarios through [run_coverage_simulations()]
#' and returns an object that can be printed or passed to
#' [format_benchmark_report()].
#'
#' @param scenarios Scenario data frame from [adoption_benchmark_scenarios()].
#' @param output_dir Directory for optional benchmark outputs.
#' @param write_results Write CSV/RDS outputs when `TRUE`.
#' @param write_report Write a Markdown benchmark report when `TRUE`.
#' @param report_path Markdown report path.
#' @param report_title Report title.
#' @param ... Additional arguments passed to [run_coverage_simulations()], such
#'   as `seed`, `max_outer_iter`, `max_inner_iter`, or `max_elapsed_sec`.
#'
#' @return An object of class `gamlss_longitudinal_adoption_benchmark`.
#' @export
run_adoption_benchmarks <- function(
    scenarios = adoption_benchmark_scenarios(),
    output_dir = file.path("results", "adoption_benchmarks"),
    write_results = TRUE,
    write_report = write_results,
    report_path = file.path(output_dir, "adoption_benchmark_report.md"),
    report_title = "gamlss.longitudinal Adoption Benchmark Report",
    ...) {
  scenarios <- as.data.frame(scenarios, stringsAsFactors = FALSE)
  required <- c("scenario", "label", "family", "copula", "design", "n_subject", "n_time")
  missing <- setdiff(required, names(scenarios))
  if (length(missing) > 0L) {
    stop("'scenarios' is missing required column(s): ", paste(missing, collapse = ", "), call. = FALSE)
  }

  if (nrow(scenarios) == 0L) {
    stop("'scenarios' must contain at least one row.", call. = FALSE)
  }

  results <- list()
  for (i in seq_len(nrow(scenarios))) {
    scenario <- scenarios[i, , drop = FALSE]
    methods <- .adoption_benchmark_list_value(scenarios, "methods", i, .coverage_default_methods())
    include_mixed <- isTRUE(.adoption_benchmark_scalar_value(scenarios, "include_mixed", i, FALSE))

    case_output_dir <- file.path(output_dir, scenario$scenario[[1]])
    case_results <- run_coverage_simulations(
      families = scenario$family[[1]],
      copulas = scenario$copula[[1]],
      methods = methods,
      designs = scenario$design[[1]],
      include_mixed = include_mixed,
      output_dir = case_output_dir,
      write_results = write_results,
      write_summary = FALSE,
      n = scenario$n_subject[[1]],
      times = seq_len(scenario$n_time[[1]]),
      dependence = .adoption_benchmark_scalar_value(scenarios, "dependence", i, "moderate"),
      missingness = .adoption_benchmark_scalar_value(scenarios, "missingness", i, "none"),
      start_mode = .adoption_benchmark_scalar_value(scenarios, "start_mode", i, "default"),
      ...
    )

    case_results$benchmark_scenario <- scenario$scenario[[1]]
    case_results$scenario_label <- scenario$label[[1]]
    results[[length(results) + 1L]] <- case_results
  }

  result_table <- do.call(rbind, results)
  rownames(result_table) <- NULL

  out <- list(
    results = result_table,
    scenarios = scenarios,
    output_dir = output_dir,
    report_path = if (isTRUE(write_report)) report_path else NULL
  )
  class(out) <- "gamlss_longitudinal_adoption_benchmark"

  if (isTRUE(write_results)) {
    dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)
    utils::write.csv(result_table, file.path(output_dir, "adoption_benchmark_results.csv"), row.names = FALSE)
    saveRDS(out, file.path(output_dir, "adoption_benchmark_results.rds"))
  }

  if (isTRUE(write_report)) {
    write_benchmark_report(out, path = report_path, title = report_title)
  }

  out
}

#' @export
print.gamlss_longitudinal_adoption_benchmark <- function(x, ...) {
  cat("gamlss.longitudinal adoption benchmark\n")
  cat("  Scenarios:", nrow(x$scenarios), "\n")
  cat("  Result rows:", nrow(x$results), "\n")
  if (!is.null(x$report_path)) {
    cat("  Report:", x$report_path, "\n")
  }
  invisible(x)
}

.adoption_benchmark_list_value <- function(scenarios, column, row, default) {
  if (!column %in% names(scenarios)) return(default)
  value <- scenarios[[column]]
  if (is.list(value)) {
    value <- value[[row]]
  } else {
    value <- value[row]
  }
  if (is.null(value) || length(value) == 0L || all(is.na(value))) default else as.character(value)
}

.adoption_benchmark_scalar_value <- function(scenarios, column, row, default) {
  if (!column %in% names(scenarios)) return(default)
  value <- scenarios[[column]]
  if (is.list(value)) {
    value <- value[[row]]
  } else {
    value <- value[row]
  }
  if (is.null(value) || length(value) == 0L || is.na(value[[1]])) default else value[[1]]
}
