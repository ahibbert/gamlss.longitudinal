#' List available standard longitudinal benchmark comparators
#'
#' @return A data frame describing optional comparator backends.
#' @export
benchmark_comparator_status <- function() {
  packages <- c("geepack", "lme4", "mgcv", "glmmTMB")
  data.frame(
    comparator = c("gee", "glmm", "gam", "glmmTMB"),
    comparator_class = c("gee", "glmm", "gamm", "glmm"),
    estimator = c("geepack::geeglm", "lme4::lmer/glmer", "mgcv::gam", "glmmTMB::glmmTMB"),
    package = packages,
    available = vapply(packages, requireNamespace, logical(1), quietly = TRUE),
    role = c(
      "marginal mean baseline with working correlation",
      "random-intercept conditional mean baseline",
      "smooth mean baseline with optional subject random effect",
      "optional flexible GLMM baseline"
    ),
    stringsAsFactors = FALSE
  )
}

#' Return named adoption benchmark scenarios
#'
#' `adoption_benchmark_scenarios()` defines a small, opinionated benchmark plan
#' for comparing `gamlss.longitudinal` with GEE/GLMM/GAM defaults. The scenarios
#' are designed to be executable through [run_adoption_benchmarks()] while also
#' documenting the applied claim each scenario is meant to test.
#'
#' @param scenarios Optional character vector of scenario names to keep.
#'
#' @return A data frame with one row per benchmark scenario.
#' @export
adoption_benchmark_scenarios <- function(scenarios = NULL) {
  out <- data.frame(
    scenario = c(
      "gaussian_heteroskedastic",
      "gamma_positive",
      "poisson_count",
      "time_varying_dependence",
      "missing_visits"
    ),
    label = c(
      "Gaussian outcome with heteroskedasticity",
      "Positive skewed outcome",
      "Count outcome with longitudinal dependence",
      "Gaussian outcome with time-varying dependence",
      "Gaussian outcome with missing visits"
    ),
    family = c("NO", "GA", "PO", "NO", "NO"),
    copula = c("N", "N", "N", "N", "N"),
    design = c("scale", "covariate", "covariate", "time_dependence", "covariate"),
    n_subject = c(80L, 80L, 80L, 80L, 80L),
    n_time = c(3L, 3L, 3L, 4L, 4L),
    dependence = c("moderate", "moderate", "moderate", "moderate", "moderate"),
    missingness = c("none", "none", "none", "none", "drop_rows"),
    claim = c(
      "Tests whether a true scale-varying GAMLSS margin improves calibration over mean-only longitudinal baselines.",
      "Tests whether a positive GAMLSS margin improves mean, quantile, and tail behaviour versus standard mean baselines.",
      "Tests whether count margins and copula dependence improve dispersion and upper-tail behaviour.",
      "Tests whether theta formulas recover changing adjacent-time dependence that standard exchangeable baselines cannot represent.",
      "Tests whether the workflow remains stable when common follow-up visits are absent."
    ),
    primary_metrics = I(list(
      c("benchmark_mean_rmse", "benchmark_q90_mae", "benchmark_upper_tail_error_90", "benchmark_interval_coverage_95", "elapsed_sec"),
      c("benchmark_mean_rmse", "benchmark_q90_mae", "benchmark_upper_tail_error_90", "elapsed_sec"),
      c("benchmark_mean_rmse", "benchmark_q90_mae", "benchmark_upper_tail_error_90", "elapsed_sec"),
      c("benchmark_theta_time_abs_error", "elapsed_sec"),
      c("benchmark_mean_rmse", "benchmark_interval_coverage_95", "elapsed_sec")
    )),
    methods = I(list(
      c("rs_separate", "gee", "glmm", "gam"),
      c("rs_separate", "gee", "glmm", "gam"),
      c("rs_separate", "gee", "glmm", "gam"),
      c("rs_separate", "gee", "glmm", "gam"),
      c("rs_separate", "gee", "glmm", "gam")
    )),
    stringsAsFactors = FALSE
  )

  if (!is.null(scenarios)) {
    scenarios <- as.character(scenarios)
    missing <- setdiff(scenarios, out$scenario)
    if (length(missing) > 0L) {
      stop("Unknown adoption benchmark scenario(s): ", paste(missing, collapse = ", "), call. = FALSE)
    }
    out <- out[match(scenarios, out$scenario), , drop = FALSE]
  }
  rownames(out) <- NULL
  out
}

.adoption_benchmark_times <- function(n_time) {
  seq_len(as.integer(n_time)[1L])
}

.adoption_benchmark_methods <- function(row, methods) {
  if (!is.null(methods)) {
    return(unique(as.character(methods)))
  }
  unique(as.character(row$methods[[1L]]))
}

.adoption_benchmark_metrics <- function(scenarios, metrics) {
  if (!is.null(metrics)) {
    return(unique(as.character(metrics)))
  }
  unique(unlist(scenarios$primary_metrics, use.names = FALSE))
}

#' Run repeated adoption benchmark scenarios
#'
#' `run_adoption_benchmarks()` executes the named scenarios from
#' [adoption_benchmark_scenarios()] through [run_coverage_simulations()] and
#' attaches a [summarise_benchmark_results()] win/tie/loss summary. It is meant
#' for opt-in benchmark runs, vignettes, and simulation reports rather than CRAN
#' tests.
#'
#' @param scenarios Scenario data frame from [adoption_benchmark_scenarios()] or
#'   character scenario names.
#' @param reps Number of simulation replicates per scenario.
#' @param methods Optional method vector overriding scenario defaults.
#' @param metrics Optional metrics passed to [summarise_benchmark_results()].
#' @param seed Base random seed.
#' @param write_results Write CSV/RDS outputs from each simulation call.
#' @param output_dir Directory for optional outputs.
#' @param ... Additional arguments passed to [run_coverage_simulations()], such
#'   as iteration limits or elapsed-time limits.
#'
#' @return An object of class `gamlss_longitudinal_adoption_benchmark` with
#'   `results`, `summary`, and `scenarios` components.
#' @export
run_adoption_benchmarks <- function(
  scenarios = adoption_benchmark_scenarios(),
  reps = 10L,
  methods = NULL,
  metrics = NULL,
  seed = 1L,
  write_results = FALSE,
  output_dir = file.path("results", "adoption_benchmarks"),
  ...
) {
  if (is.character(scenarios)) {
    scenarios <- adoption_benchmark_scenarios(scenarios)
  } else {
    scenarios <- as.data.frame(scenarios, stringsAsFactors = FALSE)
  }
  required <- c("scenario", "family", "copula", "design", "n_subject", "n_time", "dependence", "missingness", "methods")
  missing_cols <- setdiff(required, names(scenarios))
  if (length(missing_cols) > 0L) {
    stop("'scenarios' is missing required column(s): ", paste(missing_cols, collapse = ", "), call. = FALSE)
  }
  reps <- as.integer(reps)[1L]
  if (!is.finite(reps) || reps < 1L) {
    stop("'reps' must be a positive integer.", call. = FALSE)
  }

  rows <- list()
  for (scenario_idx in seq_len(nrow(scenarios))) {
    scenario_row <- scenarios[scenario_idx, , drop = FALSE]
    scenario_methods <- .adoption_benchmark_methods(scenario_row, methods)
    for (rep_idx in seq_len(reps)) {
      result <- run_coverage_simulations(
        families = scenario_row$family,
        copulas = scenario_row$copula,
        methods = scenario_methods,
        designs = scenario_row$design,
        n = scenario_row$n_subject,
        times = .adoption_benchmark_times(scenario_row$n_time),
        seed = seed + scenario_idx * 1000L + rep_idx,
        dependence = scenario_row$dependence,
        missingness = scenario_row$missingness,
        write_results = write_results,
        output_dir = file.path(output_dir, scenario_row$scenario),
        ...
      )
      result$benchmark_scenario <- scenario_row$scenario
      result$benchmark_label <- if ("label" %in% names(scenario_row)) scenario_row$label else scenario_row$scenario
      result$benchmark_rep <- rep_idx
      rows[[length(rows) + 1L]] <- result
    }
  }
  results <- do.call(rbind, rows)
  summary_metrics <- .adoption_benchmark_metrics(scenarios, metrics)
  summary <- summarise_benchmark_results(
    results,
    metrics = summary_metrics,
    group_cols = c(
      "benchmark_scenario",
      "benchmark_rep",
      "family",
      "copula",
      "design",
      "n_subject",
      "n_time",
      "dependence",
      "missingness",
      "start_mode"
    )
  )
  out <- list(
    results = results,
    summary = summary,
    scenarios = scenarios,
    reps = reps,
    metrics = summary_metrics
  )
  class(out) <- "gamlss_longitudinal_adoption_benchmark"
  out
}

#' @export
print.gamlss_longitudinal_adoption_benchmark <- function(x, digits = max(3, getOption("digits") - 3), ...) {
  cat("\nAdoption Benchmark Run\n")
  cat("----------------------\n")
  cat("Scenarios:", paste(x$scenarios$scenario, collapse = ", "), "\n")
  cat("Replicates per scenario:", x$reps, "\n")
  cat("Rows:", nrow(x$results), "\n")
  cat("Metrics:", paste(x$metrics, collapse = ", "), "\n\n")
  print(x$summary, digits = digits)
  invisible(x)
}

.benchmark_markdown_path <- function(path) {
  if (!is.character(path) || length(path) != 1L || is.na(path) || !nzchar(path)) {
    stop("'path' must be a single non-empty file path.", call. = FALSE)
  }
  if (!grepl("\\.(md|markdown)$", path, ignore.case = TRUE)) {
    path <- paste0(path, ".md")
  }
  path
}

.benchmark_report_escape <- function(x) {
  x <- as.character(x)
  x[is.na(x)] <- ""
  gsub("|", "\\|", x, fixed = TRUE)
}

.benchmark_report_flatten <- function(x, digits = 3) {
  x <- as.data.frame(x, stringsAsFactors = FALSE)
  for (nm in names(x)) {
    if (is.list(x[[nm]]) && !is.data.frame(x[[nm]])) {
      x[[nm]] <- vapply(x[[nm]], function(value) paste(value, collapse = "; "), character(1))
    } else if (is.numeric(x[[nm]])) {
      x[[nm]] <- ifelse(is.na(x[[nm]]), "", format(signif(x[[nm]], digits), scientific = FALSE, trim = TRUE))
    } else if (is.logical(x[[nm]])) {
      x[[nm]] <- ifelse(is.na(x[[nm]]), "", as.character(x[[nm]]))
    } else {
      x[[nm]] <- .benchmark_report_escape(x[[nm]])
    }
  }
  x
}

.benchmark_report_table <- function(x, columns = NULL, max_rows = Inf, digits = 3) {
  x <- as.data.frame(x, stringsAsFactors = FALSE)
  if (!is.null(columns)) {
    columns <- intersect(columns, names(x))
    x <- x[, columns, drop = FALSE]
  }
  if (nrow(x) == 0L || ncol(x) == 0L) {
    return("_No rows._")
  }
  if (is.finite(max_rows) && nrow(x) > max_rows) {
    x <- x[seq_len(max_rows), , drop = FALSE]
  }
  x <- .benchmark_report_flatten(x, digits = digits)
  header <- paste0("| ", paste(names(x), collapse = " | "), " |")
  rule <- paste0("| ", paste(rep("---", ncol(x)), collapse = " | "), " |")
  rows <- apply(x, 1L, function(row) paste0("| ", paste(row, collapse = " | "), " |"))
  c(header, rule, rows)
}

.benchmark_report_headlines <- function(summary) {
  summary <- as.data.frame(summary, stringsAsFactors = FALSE)
  if (nrow(summary) == 0L) return(summary)
  required <- c("metric", "estimand", "method", "win_or_tie_rate", "median_score", "n_finite")
  if (!all(required %in% names(summary))) return(data.frame())
  split_key <- interaction(summary$metric, drop = TRUE)
  rows <- lapply(split(seq_len(nrow(summary)), split_key), function(idx) {
    x <- summary[idx, , drop = FALSE]
    order_idx <- order(-x$win_or_tie_rate, x$median_score, -x$n_finite, x$method, na.last = TRUE)
    x[order_idx[1L], , drop = FALSE]
  })
  out <- do.call(rbind, rows)
  rownames(out) <- NULL
  out
}

.benchmark_report_group_summary <- function(case_results, group_col) {
  case_results <- as.data.frame(case_results, stringsAsFactors = FALSE)
  required <- c(group_col, "metric", "estimand", "method", "value", "score", "result")
  if (nrow(case_results) == 0L || !all(required %in% names(case_results))) {
    return(data.frame())
  }
  split_key <- interaction(
    case_results[[group_col]],
    case_results$metric,
    case_results$method,
    drop = TRUE,
    lex.order = TRUE
  )
  rows <- lapply(split(seq_len(nrow(case_results)), split_key), function(idx) {
    x <- case_results[idx, , drop = FALSE]
    finite <- is.finite(x$value)
    data.frame(
      group = x[[group_col]][[1L]],
      metric = x$metric[[1L]],
      estimand = x$estimand[[1L]],
      method = x$method[[1L]],
      n = nrow(x),
      n_finite = sum(finite),
      finite_rate = mean(finite),
      wins = sum(x$result == "win", na.rm = TRUE),
      ties = sum(x$result == "tie", na.rm = TRUE),
      losses = sum(x$result == "loss", na.rm = TRUE),
      missing = sum(x$result == "missing", na.rm = TRUE),
      win_or_tie_rate = mean(x$result %in% c("win", "tie"), na.rm = TRUE),
      median_value = if (any(finite)) stats::median(x$value[finite], na.rm = TRUE) else NA_real_,
      median_score = if (any(is.finite(x$score))) stats::median(x$score[is.finite(x$score)], na.rm = TRUE) else NA_real_,
      stringsAsFactors = FALSE
    )
  })
  out <- do.call(rbind, rows)
  names(out)[names(out) == "group"] <- group_col
  out <- out[order(out[[group_col]], out$estimand, out$metric, -out$win_or_tie_rate, out$median_score), , drop = FALSE]
  rownames(out) <- NULL
  out
}

.benchmark_report_group_headlines <- function(group_summary, group_col) {
  group_summary <- as.data.frame(group_summary, stringsAsFactors = FALSE)
  required <- c(group_col, "metric", "method", "win_or_tie_rate", "median_score", "n_finite")
  if (nrow(group_summary) == 0L || !all(required %in% names(group_summary))) {
    return(data.frame())
  }
  split_key <- interaction(group_summary[[group_col]], group_summary$metric, drop = TRUE, lex.order = TRUE)
  rows <- lapply(split(seq_len(nrow(group_summary)), split_key), function(idx) {
    x <- group_summary[idx, , drop = FALSE]
    order_idx <- order(-x$win_or_tie_rate, x$median_score, -x$n_finite, x$method, na.last = TRUE)
    x[order_idx[1L], , drop = FALSE]
  })
  out <- do.call(rbind, rows)
  rownames(out) <- NULL
  out
}

.benchmark_report_filter_scenario_metrics <- function(scenario_summary, scenarios) {
  scenario_summary <- as.data.frame(scenario_summary, stringsAsFactors = FALSE)
  if (
    nrow(scenario_summary) == 0L ||
      is.null(scenarios) ||
      !"benchmark_scenario" %in% names(scenario_summary) ||
      !"scenario" %in% names(scenarios) ||
      !"primary_metrics" %in% names(scenarios)
  ) {
    return(scenario_summary)
  }
  scenarios <- as.data.frame(scenarios, stringsAsFactors = FALSE)
  keep <- vapply(seq_len(nrow(scenario_summary)), function(i) {
    scenario_idx <- match(scenario_summary$benchmark_scenario[[i]], scenarios$scenario)
    if (is.na(scenario_idx)) {
      return(TRUE)
    }
    metric_set <- scenarios$primary_metrics[[scenario_idx]]
    if (is.null(metric_set) || length(metric_set) == 0L) {
      return(TRUE)
    }
    scenario_summary$metric[[i]] %in% as.character(metric_set)
  }, logical(1))
  out <- scenario_summary[keep, , drop = FALSE]
  if ("label" %in% names(scenarios)) {
    scenario_idx <- match(out$benchmark_scenario, scenarios$scenario)
    out$scenario_label <- scenarios$label[scenario_idx]
    out <- out[, c("benchmark_scenario", "scenario_label", setdiff(names(out), c("benchmark_scenario", "scenario_label"))), drop = FALSE]
  }
  rownames(out) <- NULL
  out
}

.benchmark_report_inputs <- function(benchmark, metrics = NULL) {
  if (inherits(benchmark, "gamlss_longitudinal_adoption_benchmark")) {
    return(list(
      results = benchmark$results,
      summary = benchmark$summary,
      scenarios = benchmark$scenarios,
      reps = benchmark$reps,
      metrics = benchmark$metrics
    ))
  }
  results <- as.data.frame(benchmark, stringsAsFactors = FALSE)
  summary <- summarise_benchmark_results(results, metrics = metrics)
  list(
    results = results,
    summary = summary,
    scenarios = NULL,
    reps = length(unique(results$benchmark_rep %||% 1L)),
    metrics = summary$metric_catalog$metric
  )
}

`%||%` <- function(x, y) {
  if (is.null(x)) y else x
}

#' Format an adoption benchmark report
#'
#' `format_benchmark_report()` turns a benchmark object or result table into
#' Markdown sections that can be used in an investigation note, package site, or
#' manuscript supplement. The report deliberately includes comparator
#' availability and caveats so benchmark claims remain tied to the estimand and
#' local software environment.
#'
#' @param benchmark A `gamlss_longitudinal_adoption_benchmark` object from
#'   [run_adoption_benchmarks()] or a results data frame accepted by
#'   [summarise_benchmark_results()].
#' @param title Report title.
#' @param comparator_status Optional comparator availability table. Defaults to
#'   [benchmark_comparator_status()].
#' @param metrics Optional metric vector used when `benchmark` is a raw results
#'   data frame.
#' @param include_case_results Logical; include per-case win/tie/loss rows.
#' @param max_case_rows Maximum case rows to include when
#'   `include_case_results = TRUE`.
#' @param digits Significant digits for numeric report tables.
#'
#' @return A character vector containing Markdown lines.
#' @export
format_benchmark_report <- function(
  benchmark,
  title = "gamlss.longitudinal Adoption Benchmark Report",
  comparator_status = benchmark_comparator_status(),
  metrics = NULL,
  include_case_results = FALSE,
  max_case_rows = 100L,
  digits = 3
) {
  inputs <- .benchmark_report_inputs(benchmark, metrics = metrics)
  summary <- inputs$summary$summary
  case_results <- inputs$summary$case_results
  catalog <- inputs$summary$metric_catalog
  scenarios <- inputs$scenarios
  headlines <- .benchmark_report_headlines(summary)
  scenario_summary <- .benchmark_report_group_summary(case_results, "benchmark_scenario")
  scenario_summary <- .benchmark_report_filter_scenario_metrics(scenario_summary, scenarios)
  scenario_headlines <- .benchmark_report_group_headlines(scenario_summary, "benchmark_scenario")

  lines <- c(
    paste0("# ", title),
    "",
    "## Run Summary",
    "",
    paste0("- Result rows: ", nrow(inputs$results)),
    paste0("- Replicates per scenario: ", inputs$reps),
    paste0("- Metrics scored: ", paste(inputs$metrics, collapse = ", ")),
    ""
  )

  if (!is.null(scenarios)) {
    lines <- c(
      lines,
      "## Scenario Plan",
      "",
      .benchmark_report_table(
        scenarios,
        columns = c("scenario", "label", "family", "copula", "design", "n_subject", "n_time", "dependence", "missingness", "claim", "primary_metrics", "methods"),
        digits = digits
      ),
      ""
    )
  }

  lines <- c(
    lines,
    "## Comparator Availability",
    "",
    .benchmark_report_table(
      comparator_status,
      columns = c("comparator", "comparator_class", "estimator", "package", "available", "role"),
      digits = digits
    ),
    "",
    "## Headline Results",
    "",
    .benchmark_report_table(
      headlines,
      columns = c("metric", "estimand", "method", "n", "n_finite", "wins", "ties", "losses", "missing", "win_or_tie_rate", "median_value", "median_score"),
      digits = digits
    ),
    "",
    "## Scenario-Level Headline Results",
    "",
    "These rows use each scenario's declared primary metrics.",
    "",
    .benchmark_report_table(
      scenario_headlines,
      columns = c("benchmark_scenario", "scenario_label", "metric", "estimand", "method", "n", "n_finite", "finite_rate", "wins", "ties", "losses", "missing", "win_or_tie_rate", "median_value", "median_score"),
      digits = digits
    ),
    "",
    "## Scenario-Level Win/Tie/Loss Summary",
    "",
    "These rows use each scenario's declared primary metrics.",
    "",
    .benchmark_report_table(
      scenario_summary,
      columns = c("benchmark_scenario", "scenario_label", "metric", "estimand", "method", "n", "n_finite", "finite_rate", "wins", "ties", "losses", "missing", "win_or_tie_rate", "median_value", "median_score"),
      digits = digits
    ),
    "",
    "## Full Win/Tie/Loss Summary",
    "",
    .benchmark_report_table(
      summary,
      columns = c("metric", "estimand", "method", "n", "n_finite", "wins", "ties", "losses", "missing", "win_rate", "win_or_tie_rate", "median_value", "median_score"),
      digits = digits
    ),
    "",
    "## Metric Definitions",
    "",
    .benchmark_report_table(
      catalog,
      columns = c("metric", "estimand", "score_rule", "target"),
      digits = digits
    ),
    "",
    "## Interpretation Notes",
    "",
    "- GEE rows target marginal mean estimands; GLMM rows target conditional mean estimands unless explicitly transformed.",
    "- Unsupported comparator-family combinations are evidence about method scope, not numerical failures.",
    "- Runtime rankings are local to the hardware, installed packages, and convergence limits used for this run.",
    "- Short smoke runs are useful for wiring checks but should not be cited as performance evidence.",
    "- Benchmark claims should be reported by estimand: mean, interval, quantile, tail, dependence, calibration, or runtime.",
    ""
  )

  if (isTRUE(include_case_results)) {
    lines <- c(
      lines,
      "## Per-Case Results",
      "",
      .benchmark_report_table(
        case_results,
        max_rows = max_case_rows,
        digits = digits
      ),
      ""
    )
  }

  lines
}

#' Write an adoption benchmark report
#'
#' @inheritParams format_benchmark_report
#' @param path Markdown output path.
#'
#' @return The normalized report path, invisibly.
#' @export
write_benchmark_report <- function(
  benchmark,
  path = file.path("results", "adoption_benchmarks", "adoption_benchmark_report.md"),
  title = "gamlss.longitudinal Adoption Benchmark Report",
  comparator_status = benchmark_comparator_status(),
  metrics = NULL,
  include_case_results = FALSE,
  max_case_rows = 100L,
  digits = 3
) {
  path <- .benchmark_markdown_path(path)
  lines <- format_benchmark_report(
    benchmark = benchmark,
    title = title,
    comparator_status = comparator_status,
    metrics = metrics,
    include_case_results = include_case_results,
    max_case_rows = max_case_rows,
    digits = digits
  )
  dir.create(dirname(path), recursive = TRUE, showWarnings = FALSE)
  writeLines(lines, con = path, useBytes = TRUE)
  invisible(normalizePath(path, winslash = "/", mustWork = FALSE))
}

.benchmark_family <- function(family) {
  if (inherits(family, "family")) {
    return(family)
  }
  family <- match.arg(as.character(family)[1], c("gaussian", "poisson", "binomial", "gamma"))
  switch(
    family,
    gaussian = stats::gaussian(),
    poisson = stats::poisson(),
    binomial = stats::binomial(),
    gamma = stats::Gamma(link = "log")
  )
}

.benchmark_backtick <- function(x) {
  if (grepl("^[.A-Za-z][.A-Za-z0-9_]*$", x)) {
    return(x)
  }
  paste0("`", gsub("`", "\\\\`", x), "`")
}

.benchmark_formula_with_random_intercept <- function(formula, subject_var) {
  stats::as.formula(paste(deparse(formula), "+ (1 |", .benchmark_backtick(subject_var), ")"))
}

.benchmark_formula_with_subject_re_smooth <- function(formula, subject_var) {
  stats::as.formula(paste(deparse(formula), "+ s(", .benchmark_backtick(subject_var), ", bs = 're')"))
}

.benchmark_response <- function(formula, data) {
  response_var <- all.vars(stats::as.formula(formula))[1L]
  if (!response_var %in% names(data)) {
    stop("The response variable from 'formula' was not found in 'data'.", call. = FALSE)
  }
  data[[response_var]]
}

.benchmark_predict_response <- function(fit, data) {
  pred <- tryCatch(
    stats::predict(fit, newdata = data, type = "response", allow.new.levels = TRUE),
    error = function(e) {
      tryCatch(stats::predict(fit, newdata = data, type = "response"), error = function(e2) rep(NA_real_, nrow(data)))
    }
  )
  as.numeric(pred)
}

.benchmark_score_predictions <- function(y, fitted) {
  ok <- is.finite(y) & is.finite(fitted)
  if (!any(ok)) {
    return(c(mae = NA_real_, rmse = NA_real_))
  }
  err <- fitted[ok] - y[ok]
  c(
    mae = mean(abs(err)),
    rmse = sqrt(mean(err^2))
  )
}

.benchmark_fit_one <- function(data, formula, subject_var, family, comparator, correlation, add_subject_re_to_gam, ...) {
  status <- benchmark_comparator_status()
  status_row <- status[match(comparator, status$comparator), , drop = FALSE]
  if (nrow(status_row) != 1L || is.na(status_row$comparator)) {
    stop("Unknown comparator: ", comparator, call. = FALSE)
  }

  start <- Sys.time()
  warnings <- character(0)
  fit <- NULL
  error <- NULL

  if (!isTRUE(status_row$available)) {
    elapsed <- as.numeric(difftime(Sys.time(), start, units = "secs"))
    return(list(
      fit = NULL,
      row = data.frame(
        comparator = comparator,
        comparator_class = status_row$comparator_class,
        estimator = status_row$estimator,
        package = status_row$package,
        available = FALSE,
        success = FALSE,
        elapsed_sec = elapsed,
        nobs = nrow(data),
        logLik = NA_real_,
        AIC = NA_real_,
        mae = NA_real_,
        rmse = NA_real_,
        warning = NA_character_,
        error = paste0("Package '", status_row$package, "' is not installed."),
        stringsAsFactors = FALSE
      )
    ))
  }

  captured <- withCallingHandlers(
    tryCatch({
      if (identical(comparator, "gee")) {
        data$.benchmark_subject_id <- data[[subject_var]]
        do.call(
          geepack::geeglm,
          c(
            list(
              formula = formula,
              id = as.name(".benchmark_subject_id"),
              data = data,
              family = family,
              corstr = correlation
            ),
            list(...)
          )
        )
      } else if (identical(comparator, "glmm")) {
        f_re <- .benchmark_formula_with_random_intercept(formula, subject_var)
        if (identical(family$family, "gaussian")) {
          lme4::lmer(f_re, data = data, ...)
        } else {
          lme4::glmer(f_re, data = data, family = family, ...)
        }
      } else if (identical(comparator, "gam")) {
        f_gam <- if (isTRUE(add_subject_re_to_gam)) {
          .benchmark_formula_with_subject_re_smooth(formula, subject_var)
        } else {
          formula
        }
        mgcv::gam(f_gam, data = data, family = family, method = "REML", ...)
      } else {
        f_tmb <- .benchmark_formula_with_random_intercept(formula, subject_var)
        glmmTMB_fit <- getExportedValue("glmmTMB", "glmmTMB")
        glmmTMB_fit(f_tmb, data = data, family = family, ...)
      }
    }, error = function(e) {
      error <<- conditionMessage(e)
      NULL
    }),
    warning = function(w) {
      warnings <<- c(warnings, conditionMessage(w))
      invokeRestart("muffleWarning")
    }
  )
  fit <- captured
  elapsed <- as.numeric(difftime(Sys.time(), start, units = "secs"))
  success <- !is.null(fit) && is.null(error)
  y <- .benchmark_response(formula, data)
  pred <- if (success) .benchmark_predict_response(fit, data) else rep(NA_real_, length(y))
  scores <- .benchmark_score_predictions(as.numeric(y), pred)
  ll <- if (success && !identical(comparator, "gee")) {
    as.numeric(tryCatch(stats::logLik(fit), error = function(e) NA_real_))
  } else {
    NA_real_
  }
  aic <- if (success && !identical(comparator, "gee")) {
    as.numeric(tryCatch(stats::AIC(fit), error = function(e) NA_real_))
  } else {
    NA_real_
  }

  list(
    fit = fit,
    row = data.frame(
      comparator = comparator,
      comparator_class = status_row$comparator_class,
      estimator = status_row$estimator,
      package = status_row$package,
      available = TRUE,
      success = success,
      elapsed_sec = elapsed,
      nobs = nrow(data),
      logLik = ll,
      AIC = aic,
      mae = unname(scores[["mae"]]),
      rmse = unname(scores[["rmse"]]),
      warning = if (length(warnings)) paste(unique(warnings), collapse = " | ") else NA_character_,
      error = if (is.null(error)) NA_character_ else error,
      stringsAsFactors = FALSE
    )
  )
}

#' Fit standard GEE/GLMM/GAM comparators for adoption benchmarks
#'
#' `benchmark_standard_models()` is an opt-in scaffold for comparing
#' `gamlss.longitudinal` with models users already know. It is intentionally
#' narrow: fit common mean-model baselines, record whether they ran, and return
#' simple timing and response-scale prediction metrics. Simulation studies can
#' build coverage, calibration, tail, and trajectory metrics on top of this.
#'
#' @param data Long-format data frame.
#' @param formula Mean-model formula used by the comparator models.
#' @param subject_var Subject identifier column.
#' @param family A base R family object or one of `"gaussian"`, `"poisson"`,
#'   `"binomial"`, or `"gamma"`.
#' @param comparators Character vector containing any of `"gee"`, `"glmm"`,
#'   `"gam"`, and `"glmmTMB"`.
#' @param correlation Working correlation passed to `geepack::geeglm()`.
#' @param add_subject_re_to_gam Logical; add `s(subject, bs = "re")` to the GAM
#'   comparator.
#' @param ... Additional arguments passed to each comparator fit.
#'
#' @return An object of class `gamlss_longitudinal_benchmark` with `results`
#'   and `fits` components.
#' @export
benchmark_standard_models <- function(
  data,
  formula,
  subject_var,
  family = "gaussian",
  comparators = c("gee", "glmm", "gam"),
  correlation = "exchangeable",
  add_subject_re_to_gam = TRUE,
  ...
) {
  data <- as.data.frame(data, stringsAsFactors = FALSE)
  if (!is.character(subject_var) || length(subject_var) != 1L || !subject_var %in% names(data)) {
    stop("'subject_var' must be a single column name in 'data'.", call. = FALSE)
  }
  formula <- stats::as.formula(formula)
  family <- .benchmark_family(family)
  comparators <- unique(as.character(comparators))
  valid <- benchmark_comparator_status()$comparator
  bad <- setdiff(comparators, valid)
  if (length(bad) > 0L) {
    stop("Unknown comparator(s): ", paste(bad, collapse = ", "), call. = FALSE)
  }

  runs <- lapply(comparators, function(comparator) {
    .benchmark_fit_one(
      data = data,
      formula = formula,
      subject_var = subject_var,
      family = family,
      comparator = comparator,
      correlation = correlation,
      add_subject_re_to_gam = add_subject_re_to_gam,
      ...
    )
  })
  out <- list(
    results = do.call(rbind, lapply(runs, `[[`, "row")),
    fits = stats::setNames(lapply(runs, `[[`, "fit"), comparators),
    formula = formula,
    subject_var = subject_var,
    family = family$family
  )
  class(out) <- "gamlss_longitudinal_benchmark"
  out
}

#' @export
print.gamlss_longitudinal_benchmark <- function(x, digits = max(3, getOption("digits") - 3), ...) {
  cat("\nStandard Longitudinal Comparator Benchmark\n")
  cat("-----------------------------------------\n")
  cat("Formula:", deparse(x$formula), "\n")
  cat("Subject:", x$subject_var, "\n")
  cat("Family:", x$family, "\n\n")
  print(x$results, digits = digits, row.names = FALSE)
  invisible(x)
}

.benchmark_metric_catalog <- function() {
  data.frame(
    metric = c(
      "benchmark_mean_rmse",
      "benchmark_mean_mae",
      "benchmark_mean_bias",
      "benchmark_rmse",
      "benchmark_mae",
      "benchmark_q90_mae",
      "benchmark_upper_tail_error_90",
      "benchmark_theta_time_abs_error",
      "benchmark_interval_coverage_95",
      "benchmark_pit_mean_abs_error",
      "benchmark_pit_ks_p_value",
      "benchmark_tail_error_lower_05",
      "benchmark_tail_error_upper_05",
      "elapsed_sec"
    ),
    estimand = c(
      "mean",
      "mean",
      "mean",
      "observed_response",
      "observed_response",
      "quantile",
      "tail",
      "dependence",
      "interval",
      "calibration",
      "calibration",
      "tail",
      "tail",
      "runtime"
    ),
    score_rule = c(
      "lower",
      "lower",
      "absolute",
      "lower",
      "lower",
      "lower",
      "absolute",
      "lower",
      "target",
      "lower",
      "higher",
      "absolute",
      "absolute",
      "lower"
    ),
    target = c(NA, NA, 0, NA, NA, NA, 0, NA, 0.95, NA, 1, 0, 0, NA),
    stringsAsFactors = FALSE
  )
}

.benchmark_metric_score <- function(value, rule, target) {
  if (identical(rule, "lower")) {
    return(value)
  }
  if (identical(rule, "higher")) {
    return(1 - value)
  }
  if (identical(rule, "absolute")) {
    return(abs(value))
  }
  abs(value - target)
}

.benchmark_case_label <- function(score, best_score, tie_tolerance, absolute_tolerance) {
  if (!is.finite(score) || !is.finite(best_score)) {
    return("missing")
  }
  gap <- score - best_score
  if (gap <= absolute_tolerance) {
    return("win")
  }
  scale <- max(abs(best_score), absolute_tolerance)
  if ((gap / scale) <= tie_tolerance) {
    return("tie")
  }
  "loss"
}

#' Summarise benchmark simulation results into win/tie/loss tables
#'
#' `summarise_benchmark_results()` turns rows from [run_coverage_simulations()]
#' or `benchmark_standard_models()`-style outputs into compact evidence tables.
#' Each metric is scored according to its estimand: lower error is better,
#' coverage targets 0.95, calibration/tail errors target 0, and runtime is
#' lower-is-better.
#'
#' @param results Data frame returned by [run_coverage_simulations()].
#' @param metrics Optional metric columns to summarise. Defaults to all known
#'   benchmark metrics present in `results`.
#' @param group_cols Columns defining one simulation case before method-level
#'   comparison.
#' @param tie_tolerance Relative gap from the case-best score treated as a tie.
#' @param absolute_tolerance Absolute score gap treated as a tie/win.
#'
#' @return An object of class `gamlss_longitudinal_benchmark_summary` with
#'   `summary`, `case_results`, and `metric_catalog` components.
#' @export
summarise_benchmark_results <- function(
  results,
  metrics = NULL,
  group_cols = c("family", "copula", "design", "n_subject", "n_time", "dependence", "missingness", "start_mode"),
  tie_tolerance = 0.05,
  absolute_tolerance = 1e-8
) {
  results <- as.data.frame(results, stringsAsFactors = FALSE)
  if (!"method" %in% names(results)) {
    stop("'results' must contain a 'method' column.", call. = FALSE)
  }
  catalog <- .benchmark_metric_catalog()
  if (is.null(metrics)) {
    metrics <- intersect(catalog$metric, names(results))
  } else {
    metrics <- intersect(as.character(metrics), names(results))
  }
  metrics <- metrics[metrics %in% catalog$metric]
  if (length(metrics) == 0L) {
    stop("No supported benchmark metric columns were found in 'results'.", call. = FALSE)
  }
  group_cols <- intersect(group_cols, names(results))
  if (length(group_cols) == 0L) {
    results$.benchmark_case <- 1L
    group_cols <- ".benchmark_case"
  }
  group_key <- interaction(results[group_cols], drop = TRUE, lex.order = TRUE)

  case_rows <- list()
  for (metric in metrics) {
    meta <- catalog[match(metric, catalog$metric), , drop = FALSE]
    value <- suppressWarnings(as.numeric(results[[metric]]))
    score <- .benchmark_metric_score(value, meta$score_rule, meta$target)
    for (key in unique(group_key)) {
      idx <- which(group_key == key)
      finite_idx <- idx[is.finite(score[idx])]
      if (length(finite_idx) == 0L) next
      best <- min(score[finite_idx], na.rm = TRUE)
      labels <- vapply(score[idx], .benchmark_case_label, character(1),
        best_score = best,
        tie_tolerance = tie_tolerance,
        absolute_tolerance = absolute_tolerance
      )
      case_rows[[length(case_rows) + 1L]] <- data.frame(
        results[idx, group_cols, drop = FALSE],
        method = results$method[idx],
        metric = metric,
        estimand = meta$estimand,
        value = value[idx],
        score = score[idx],
        best_score = best,
        result = labels,
        stringsAsFactors = FALSE,
        row.names = NULL
      )
    }
  }
  case_results <- if (length(case_rows) == 0L) {
    data.frame()
  } else {
    do.call(rbind, case_rows)
  }

  if (nrow(case_results) == 0L) {
    summary <- data.frame()
  } else {
    split_key <- interaction(case_results$metric, case_results$method, drop = TRUE, lex.order = TRUE)
    summary_rows <- lapply(split(seq_len(nrow(case_results)), split_key), function(idx) {
      x <- case_results[idx, , drop = FALSE]
      finite <- is.finite(x$value)
      data.frame(
        metric = x$metric[[1L]],
        estimand = x$estimand[[1L]],
        method = x$method[[1L]],
        n = nrow(x),
        n_finite = sum(finite),
        wins = sum(x$result == "win", na.rm = TRUE),
        ties = sum(x$result == "tie", na.rm = TRUE),
        losses = sum(x$result == "loss", na.rm = TRUE),
        missing = sum(x$result == "missing", na.rm = TRUE),
        win_rate = mean(x$result == "win", na.rm = TRUE),
        win_or_tie_rate = mean(x$result %in% c("win", "tie"), na.rm = TRUE),
        median_value = if (any(finite)) stats::median(x$value[finite], na.rm = TRUE) else NA_real_,
        median_score = if (any(is.finite(x$score))) stats::median(x$score[is.finite(x$score)], na.rm = TRUE) else NA_real_,
        stringsAsFactors = FALSE
      )
    })
    summary <- do.call(rbind, summary_rows)
    summary <- summary[order(summary$estimand, summary$metric, -summary$win_or_tie_rate, summary$median_score), , drop = FALSE]
    rownames(summary) <- NULL
  }

  out <- list(
    summary = summary,
    case_results = case_results,
    metric_catalog = catalog[catalog$metric %in% metrics, , drop = FALSE]
  )
  class(out) <- "gamlss_longitudinal_benchmark_summary"
  out
}

#' @export
print.gamlss_longitudinal_benchmark_summary <- function(x, digits = max(3, getOption("digits") - 3), ...) {
  cat("\nBenchmark Win/Tie/Loss Summary\n")
  cat("------------------------------\n")
  if (nrow(x$summary) == 0L) {
    cat("No finite benchmark metrics to summarise.\n")
  } else {
    print(x$summary, digits = digits, row.names = FALSE)
  }
  invisible(x)
}
