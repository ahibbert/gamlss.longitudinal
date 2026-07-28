# Consolidated support module for reviewer navigation.
# Source blocks are grouped mechanically from the files named below.

# ---- benchmark-report.R ----

#' Format a benchmark report

#'

#' `format_benchmark_report()` turns a benchmark object or result table into

#' Markdown sections that can be used in an investigation note, package site, or

#' manuscript supplement. The report deliberately includes comparator

#' availability and caveats so benchmark claims remain tied to the estimand and

#' local software environment.

#'

#' @param benchmark A `gamlss_longitudinal_benchmark` object from

#'   [benchmark_standard_models()] or a results data frame accepted by

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
    title = "gamlss.longitudinal Benchmark Report",
    comparator_status = benchmark_comparator_status(),
    metrics = NULL,
    include_case_results = FALSE,
    max_case_rows = 100L,
    digits = 3) {
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
    .benchmark_report_run_summary_section(title, inputs),
    .benchmark_report_scenario_plan_section(scenarios, digits = digits),
    .benchmark_report_comparator_section(comparator_status, digits = digits),
    .benchmark_report_headline_section(headlines, digits = digits),
    .benchmark_report_benchmark_summary_section(inputs$interpretation),
    .benchmark_report_scenario_headline_section(scenario_headlines, digits = digits),
    .benchmark_report_scenario_summary_section(scenario_summary, digits = digits),
    .benchmark_report_full_summary_section(summary, digits = digits),
    .benchmark_report_metric_definitions_section(catalog, digits = digits),
    .benchmark_report_interpretation_notes_section(),
    .benchmark_report_case_results_section(
      case_results,
      include_case_results = include_case_results,
      max_case_rows = max_case_rows,
      digits = digits
    )
  )

  lines
}

#' Write a benchmark report

#'

#' @inheritParams format_benchmark_report

#' @param path Markdown output path.

#'

#' @return The normalized report path, invisibly.

#' @export

write_benchmark_report <- function(
    benchmark,
    path = file.path("results", "benchmark_report.md"),
    title = "gamlss.longitudinal Benchmark Report",
    comparator_status = benchmark_comparator_status(),
    metrics = NULL,
    include_case_results = FALSE,
    max_case_rows = 100L,
    digits = 3) {
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

# ---- benchmark-report-inputs.R ----

.benchmark_standard_metric_columns <- function(results) {
  catalog_metrics <- .benchmark_metric_catalog()$metric

  requested <- c("elapsed_sec", "nobs", "logLik", "AIC", "mae", "rmse", catalog_metrics)

  unique(intersect(requested, names(results)))
}

.benchmark_model_status_table <- function(results) {
  results <- as.data.frame(results, stringsAsFactors = FALSE)

  columns <- c(
    "method",
    "estimator",
    "success"
  )

  columns <- intersect(columns, names(results))

  out <- results[, columns, drop = FALSE]

  rownames(out) <- NULL

  out
}

.benchmark_method_display_name <- function(method, estimator = NA_character_) {
  method <- as.character(method)

  estimator <- as.character(estimator %||% rep(NA_character_, length(method)))

  out <- method

  out[method %in% "gamlss.longitudinal"] <- "gamlss.longitudinal"

  out[method %in% "gamlss_longitudinal"] <- "gamlss.longitudinal"

  out[method %in% "gee"] <- "geeglm"

  out[method %in% "glmm"] <- "glmer"

  out[method %in% "gamm"] <- "gamm"

  out[method %in% "gam"] <- "gam"

  out[method %in% "glm"] <- "glm"

  if (length(estimator) == length(out)) {
    out[grepl("geeglm", estimator, fixed = TRUE)] <- "geeglm"

    out[grepl("lmer/glmer", estimator, fixed = TRUE)] <- "glmer"
  }

  out
}

.benchmark_method_display_list <- function(x) {
  vapply(strsplit(as.character(x), ",\\s*"), function(parts) {
    paste(.benchmark_method_display_name(parts), collapse = ", ")
  }, character(1))
}

.benchmark_included_methods_line <- function(results) {
  results <- as.data.frame(results, stringsAsFactors = FALSE)

  if (nrow(results) == 0L || !"method" %in% names(results)) {
    return("Methods included: none")
  }

  success <- if ("success" %in% names(results)) !is.na(results$success) & results$success == TRUE else rep(TRUE, nrow(results))

  included <- results[success, , drop = FALSE]

  if (nrow(included) == 0L) {
    return("Methods included: none")
  }

  names <- .benchmark_method_display_name(included$method, included$estimator %||% NA_character_)

  paste0("Methods included: ", paste(unique(names), collapse = ", "))
}

.benchmark_metric_matrix <- function(results) {
  results <- as.data.frame(results, stringsAsFactors = FALSE)

  if (!"method" %in% names(results)) {
    return(data.frame())
  }

  metrics <- .benchmark_standard_metric_columns(results)

  if (length(metrics) == 0L) {
    return(data.frame())
  }

  has_value <- vapply(results[, metrics, drop = FALSE], function(x) any(!is.na(x)), logical(1))

  metrics <- metrics[has_value]

  if (length(metrics) == 0L) {
    return(data.frame())
  }

  method_names <- make.unique(as.character(results$method))

  values <- t(as.matrix(results[, metrics, drop = FALSE]))

  colnames(values) <- method_names

  out <- data.frame(metric = rownames(values), values, check.names = FALSE, stringsAsFactors = FALSE)

  rownames(out) <- NULL

  out
}

.benchmark_report_inputs <- function(benchmark, metrics = NULL) {
  if (inherits(benchmark, "gamlss_longitudinal_adoption_benchmark")) {
    results <- as.data.frame(benchmark$results, stringsAsFactors = FALSE)
    summary <- summarise_benchmark_results(results, metrics = metrics)
    primary_method <- .benchmark_infer_primary_method(results)

    return(list(
      results = results,
      summary = summary,
      scenarios = benchmark$scenarios,
      reps = length(unique(results$benchmark_rep %||% 1L)),
      metrics = summary$metric_catalog$metric,
      interpretation = .benchmark_interpretation(results, primary_method = primary_method, metrics = metrics),
      primary_method = primary_method
    ))
  }

  if (inherits(benchmark, "gamlss_longitudinal_benchmark")) {
    interpretation <- benchmark$interpretation %||% .benchmark_interpretation(benchmark$results, primary_method = benchmark$primary_method)

    summary <- interpretation$summary %||% summarise_benchmark_results(benchmark$results, group_cols = NULL)

    return(list(
      results = benchmark$results,
      summary = summary,
      scenarios = NULL,
      reps = 1L,
      metrics = summary$metric_catalog$metric,
      interpretation = interpretation,
      primary_method = benchmark$primary_method
    ))
  }

  results <- as.data.frame(benchmark, stringsAsFactors = FALSE)

  summary <- summarise_benchmark_results(results, metrics = metrics)

  primary_method <- .benchmark_infer_primary_method(results)

  list(
    results = results,
    summary = summary,
    scenarios = NULL,
    reps = length(unique(results$benchmark_rep %||% 1L)),
    metrics = summary$metric_catalog$metric,
    interpretation = .benchmark_interpretation(results, primary_method = primary_method, metrics = metrics),
    primary_method = primary_method
  )
}

`%||%` <- function(x, y) {
  if (is.null(x)) y else x
}

# ---- benchmark-report-metadata.R ----

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

.benchmark_metric_domain <- function(metric, estimand = NULL) {
  metric <- as.character(metric)

  estimand <- as.character(estimand %||% rep(NA_character_, length(metric)))

  out <- rep("other", length(metric))

  out[grepl("^benchmark_mean_", metric) | metric %in% c("benchmark_rmse", "benchmark_mae")] <- "mean / marginal response fit"

  out[metric %in% c(
    "benchmark_q90_mae",
    "benchmark_neg_log_score",
    "benchmark_upper_tail_error_90",
    "benchmark_interval_coverage_95",
    "benchmark_interval_width_95",
    "benchmark_pit_mean_abs_error",
    "benchmark_pit_ks_p_value",
    "benchmark_tail_error_lower_05",
    "benchmark_tail_error_upper_05"
  ) | estimand %in% c("density", "quantile", "tail", "interval", "calibration")] <- "distributional prediction / shape"

  out[metric %in% "benchmark_theta_time_abs_error" | estimand %in% "dependence"] <- "dependence"

  out[metric %in% "elapsed_sec" | estimand %in% "runtime"] <- "runtime"

  out[grepl("^smooth_", metric) | estimand %in% "smooth"] <- "smooth recovery"

  out
}

.benchmark_infer_primary_method <- function(results) {
  results <- as.data.frame(results, stringsAsFactors = FALSE)

  if (!"method" %in% names(results)) {
    return(NULL)
  }

  methods <- unique(as.character(results$method))

  if ("gamlss.longitudinal" %in% methods) {
    return("gamlss.longitudinal")
  }

  if ("rs_separate" %in% methods) {
    return("rs_separate")
  }

  NULL
}

# ---- benchmark-report-headlines.R ----

.benchmark_report_headlines <- function(summary) {
  summary <- as.data.frame(summary, stringsAsFactors = FALSE)

  if (nrow(summary) == 0L) {
    return(summary)
  }

  required <- c("metric", "estimand", "method", "win_or_tie_rate", "median_score", "n_finite")

  if (!all(required %in% names(summary))) {
    return(data.frame())
  }

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

# ---- benchmark-report-sections.R ----

.benchmark_report_run_summary_section <- function(title, inputs) {
  c(
    paste0("# ", title),
    "",
    "## Run Summary",
    "",
    paste0("- Result rows: ", nrow(inputs$results)),
    paste0("- Replicates per scenario: ", inputs$reps),
    paste0("- Metrics scored: ", paste(inputs$metrics, collapse = ", ")),
    ""
  )
}

.benchmark_report_scenario_plan_section <- function(scenarios, digits = 3) {
  if (is.null(scenarios)) {
    return(character(0))
  }

  c(
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

.benchmark_report_comparator_section <- function(comparator_status, digits = 3) {
  c(
    "## Comparator Availability",
    "",
    .benchmark_report_table(
      comparator_status,
      columns = c("comparator", "comparator_class", "estimator", "package", "available", "role"),
      digits = digits
    ),
    ""
  )
}

.benchmark_report_headline_section <- function(headlines, digits = 3) {
  c(
    "## Headline Results",
    "",
    .benchmark_report_table(
      headlines,
      columns = c("metric", "estimand", "method", "n", "n_finite", "wins", "ties", "losses", "missing", "win_or_tie_rate", "median_value", "median_score"),
      digits = digits
    ),
    ""
  )
}

.benchmark_report_benchmark_summary_section <- function(interpretation) {
  c(
    "## Benchmark Summary",
    "",
    .benchmark_interpretation_lines(interpretation),
    ""
  )
}

.benchmark_report_scenario_headline_section <- function(scenario_headlines, digits = 3) {
  c(
    "## Scenario-Level Headline Results",
    "",
    "These rows use each scenario's declared primary metrics.",
    "",
    .benchmark_report_table(
      scenario_headlines,
      columns = c("benchmark_scenario", "scenario_label", "metric", "estimand", "method", "n", "n_finite", "finite_rate", "wins", "ties", "losses", "missing", "win_or_tie_rate", "median_value", "median_score"),
      digits = digits
    ),
    ""
  )
}

.benchmark_report_scenario_summary_section <- function(scenario_summary, digits = 3) {
  c(
    "## Scenario-Level Win/Tie/Loss Summary",
    "",
    "These rows use each scenario's declared primary metrics.",
    "",
    .benchmark_report_table(
      scenario_summary,
      columns = c("benchmark_scenario", "scenario_label", "metric", "estimand", "method", "n", "n_finite", "finite_rate", "wins", "ties", "losses", "missing", "win_or_tie_rate", "median_value", "median_score"),
      digits = digits
    ),
    ""
  )
}

.benchmark_report_full_summary_section <- function(summary, digits = 3) {
  c(
    "## Full Win/Tie/Loss Summary",
    "",
    .benchmark_report_table(
      summary,
      columns = c("metric", "estimand", "method", "n", "n_finite", "wins", "ties", "losses", "missing", "win_rate", "win_or_tie_rate", "median_value", "median_score"),
      digits = digits
    ),
    ""
  )
}

.benchmark_report_metric_definitions_section <- function(catalog, digits = 3) {
  c(
    "## Metric Definitions",
    "",
    .benchmark_report_table(
      catalog,
      columns = c("metric", "estimand", "score_rule", "target"),
      digits = digits
    ),
    ""
  )
}

.benchmark_report_interpretation_notes_section <- function() {
  c(
    "## Interpretation Notes",
    "",
    "- GEE rows target marginal mean estimands; GLMM benchmark predictions use population-level fixed-effect predictions when scored against marginal truth.",
    "- Unsupported comparator-family combinations are evidence about method scope, not numerical failures.",
    "- Runtime rankings are local to the hardware, installed packages, and convergence limits used for this run.",
    "- Short smoke runs are useful for wiring checks but should not be cited as performance evidence.",
    "- Benchmark claims should be reported by estimand: mean, interval, quantile, tail, dependence, calibration, or runtime.",
    ""
  )
}

.benchmark_report_case_results_section <- function(case_results,
                                                   include_case_results = FALSE,
                                                   max_case_rows = 100L,
                                                   digits = 3) {
  if (!isTRUE(include_case_results)) {
    return(character(0))
  }

  c(
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

# ---- benchmark-report-interpretation.R ----

.benchmark_interpretation_unsupported <- function(results) {
  results <- as.data.frame(results, stringsAsFactors = FALSE)

  if (nrow(results) == 0L || !"method" %in% names(results)) {
    return(data.frame())
  }

  available <- if ("available" %in% names(results)) !is.na(results$available) & results$available == TRUE else rep(TRUE, nrow(results))

  success <- if ("success" %in% names(results)) !is.na(results$success) & results$success == TRUE else rep(TRUE, nrow(results))

  eligible <- available & success

  if (all(eligible)) {
    return(data.frame())
  }

  reason <- rep("unsupported or failed", nrow(results))

  if ("available" %in% names(results)) {
    reason[!available] <- "unavailable"
  }

  if ("success" %in% names(results)) {
    reason[available & !success] <- "failed or unsupported"
  }

  error <- if ("error" %in% names(results)) as.character(results$error) else NA_character_

  out <- data.frame(
    method = as.character(results$method[!eligible]),
    comparator = if ("comparator" %in% names(results)) as.character(results$comparator[!eligible]) else as.character(results$method[!eligible]),
    reason = reason[!eligible],
    error = error[!eligible],
    stringsAsFactors = FALSE
  )

  out <- unique(out)

  rownames(out) <- NULL

  out
}

.benchmark_interpretation <- function(
    results,
    primary_method = NULL,
    metrics = NULL,
    tie_tolerance = 0.05,
    absolute_tolerance = 1e-8) {
  results <- as.data.frame(results, stringsAsFactors = FALSE)

  if (!"method" %in% names(results)) {
    stop("'results' must contain a 'method' column.", call. = FALSE)
  }

  if (!is.null(primary_method)) {
    primary_method <- as.character(primary_method)

    primary_method <- primary_method[primary_method %in% as.character(results$method)]
  }

  unsupported <- .benchmark_interpretation_unsupported(results)

  available <- if ("available" %in% names(results)) !is.na(results$available) & results$available == TRUE else rep(TRUE, nrow(results))

  success <- if ("success" %in% names(results)) !is.na(results$success) & results$success == TRUE else rep(TRUE, nrow(results))

  eligible <- results[available & success, , drop = FALSE]

  empty_leaders <- data.frame(
    domain = character(),
    metric = character(),
    estimand = character(),
    primary_leads_or_ties = logical(),
    standard_leads_or_ties = logical(),
    primary_methods = character(),
    standard_methods = character(),
    leading_methods = character(),
    stringsAsFactors = FALSE
  )

  summary <- NULL

  leaders <- empty_leaders

  if (nrow(eligible) > 0L) {
    summary <- tryCatch(
      summarise_benchmark_results(

        eligible,
        metrics = metrics,
        group_cols = NULL,
        tie_tolerance = tie_tolerance,
        absolute_tolerance = absolute_tolerance
      ),
      error = function(e) NULL
    )
  }

  if (!is.null(summary) && nrow(summary$case_results) > 0L) {
    case_results <- summary$case_results

    finite <- is.finite(case_results$value)

    leading <- case_results[finite & case_results$result %in% c("win", "tie"), , drop = FALSE]

    if (nrow(leading) > 0L) {
      rows <- lapply(split(seq_len(nrow(leading)), leading$metric, drop = TRUE), function(idx) {
        x <- leading[idx, , drop = FALSE]

        primary <- x$method[x$method %in% primary_method]

        standard <- x$method[!(x$method %in% primary_method)]

        data.frame(
          domain = .benchmark_metric_domain(x$metric[[1L]], x$estimand[[1L]]),
          metric = x$metric[[1L]],
          estimand = x$estimand[[1L]],
          primary_leads_or_ties = length(primary) > 0L,
          standard_leads_or_ties = length(standard) > 0L,
          primary_methods = paste(unique(primary), collapse = ", "),
          standard_methods = paste(unique(standard), collapse = ", "),
          leading_methods = paste(unique(x$method), collapse = ", "),
          stringsAsFactors = FALSE
        )
      })

      leaders <- do.call(rbind, rows)

      leaders <- leaders[order(leaders$domain, leaders$metric), , drop = FALSE]

      rownames(leaders) <- NULL
    }
  }

  out <- list(
    primary_method = primary_method,
    leaders = leaders,
    unsupported = unsupported,
    summary = summary
  )

  class(out) <- "gamlss_longitudinal_benchmark_interpretation"

  out
}

# ---- benchmark-report-interpretation-format.R ----

.benchmark_interpretation_join_metrics <- function(x) {
  if (length(x) == 0L) {
    return("none")
  }

  paste(unique(x), collapse = ", ")
}

.benchmark_interpretation_area <- function(domain) {
  out <- as.character(domain)

  out[out == "distributional prediction / shape"] <- "Distributional"

  out[out == "mean / marginal response fit"] <- "Mean fit"

  out[out == "dependence"] <- "Dependence"

  out[out == "runtime"] <- "Runtime"

  out[out == "smooth recovery"] <- "Smooth recovery"

  out[out == "other"] <- "Other"

  out
}

.benchmark_interpretation_lines <- function(interpretation) {
  if (is.null(interpretation) || !inherits(interpretation, "gamlss_longitudinal_benchmark_interpretation")) {
    return(character())
  }

  leaders <- as.data.frame(interpretation$leaders, stringsAsFactors = FALSE)

  primary_method <- interpretation$primary_method

  if (nrow(leaders) == 0L) {
    lines <- "- No finite supported benchmark metrics were available for interpretation."
  } else {
    domains <- unique(leaders$domain)

    lines <- unlist(lapply(domains, function(domain) {
      x <- leaders[leaders$domain == domain, , drop = FALSE]

      if (length(primary_method) > 0L) {
        c(
          paste0("- ", domain, ": ", primary_method[[1L]], " leads or ties on: ", .benchmark_interpretation_join_metrics(x$metric[x$primary_leads_or_ties])),
          paste0("- ", domain, ": Standard comparators lead or tie on: ", .benchmark_interpretation_join_metrics(x$metric[x$standard_leads_or_ties]))
        )
      } else {
        paste0(
          "- ",
          domain,
          ": leading method(s): ",
          paste(paste0(x$metric, " (", x$leading_methods, ")"), collapse = "; ")
        )
      }
    }), use.names = FALSE)
  }

  unsupported <- as.data.frame(interpretation$unsupported, stringsAsFactors = FALSE)

  if (nrow(unsupported) > 0L) {
    unsupported_label <- paste(unique(paste0(unsupported$method, " [", unsupported$reason, "]")), collapse = ", ")

    lines <- c(lines, paste0("- Unavailable or unsupported rows were not used as metric leaders: ", unsupported_label, "."))
  }

  lines
}

.benchmark_interpretation_table <- function(interpretation) {
  if (is.null(interpretation) || !inherits(interpretation, "gamlss_longitudinal_benchmark_interpretation")) {
    return(data.frame())
  }

  leaders <- as.data.frame(interpretation$leaders, stringsAsFactors = FALSE)

  if (nrow(leaders) == 0L) {
    return(data.frame())
  }

  out <- data.frame(
    area = .benchmark_interpretation_area(leaders$domain),
    metric = leaders$metric,
    best_model = .benchmark_method_display_list(leaders$leading_methods),
    stringsAsFactors = FALSE
  )

  rownames(out) <- NULL

  out
}

# ---- benchmark-report-utils.R ----

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

# ---- benchmark-report-print.R ----

.benchmark_metric_label <- function(metric) {
  metric <- as.character(metric)

  labels <- c(
    benchmark_mean_rmse = "Mean RMSE",
    benchmark_mean_mae = "Mean MAE",
    benchmark_mean_bias = "Mean bias",
    benchmark_rmse = "Observed response RMSE",
    benchmark_mae = "Observed response MAE",
    benchmark_q90_mae = "90th percentile MAE",
    benchmark_neg_log_score = "Negative log score",
    benchmark_upper_tail_error_90 = "90th percentile upper-tail error",
    benchmark_theta_time_abs_error = "Time-varying dependence absolute error",
    benchmark_interval_coverage_95 = "95% interval coverage",
    benchmark_interval_width_95 = "95% interval width",
    benchmark_pit_mean_abs_error = "PIT mean absolute error",
    benchmark_pit_ks_p_value = "PIT KS p-value",
    benchmark_tail_error_lower_05 = "Lower 5% tail error",
    benchmark_tail_error_upper_05 = "Upper 5% tail error",
    smooth_eta_rmse = "Smooth eta RMSE",
    smooth_eta_max_abs_error = "Smooth eta max absolute error",
    elapsed_sec = "Elapsed time (sec)",
    nobs = "Number of observations",
    logLik = "Log likelihood",
    AIC = "AIC",
    mae = "MAE",
    rmse = "RMSE"
  )

  out <- labels[metric]

  out[is.na(out)] <- metric[is.na(out)]

  unname(out)
}

.benchmark_print_metric_labels <- function(x) {
  x <- as.data.frame(x, stringsAsFactors = FALSE)

  if ("metric" %in% names(x)) {
    x$metric <- .benchmark_metric_label(x$metric)

    names(x)[names(x) == "metric"] <- "Metric"
  }

  if ("best_model" %in% names(x)) {
    names(x)[names(x) == "best_model"] <- "Best model"
  }

  x
}

.benchmark_print_method_columns <- function(x) {
  x <- as.data.frame(x, stringsAsFactors = FALSE)

  protected <- c("Metric", "metric", "metric_name", "term", "area", "Best model", "best_model")

  method_cols <- setdiff(names(x), protected)

  names(x)[match(method_cols, names(x))] <- .benchmark_method_display_name(method_cols)

  x
}

.benchmark_print_by_area <- function(x, digits = 3) {
  x <- as.data.frame(x, stringsAsFactors = FALSE)

  if (nrow(x) == 0L) {
    return(invisible(NULL))
  }

  if (!"area" %in% names(x)) {
    print(.benchmark_print_metric_labels(x), digits = digits, row.names = FALSE)

    return(invisible(NULL))
  }

  areas <- unique(x$area)

  for (area in areas) {
    area_rows <- x[x$area == area, setdiff(names(x), "area"), drop = FALSE]

    cat("\n", area, "\n", sep = "")

    cat(paste(rep("-", nchar(area)), collapse = ""), "\n", sep = "")

    print(.benchmark_print_metric_labels(area_rows), digits = digits, row.names = FALSE)
  }

  invisible(NULL)
}

.benchmark_print_coefficient_comparison <- function(coefficients, digits = 3) {
  if (is.null(coefficients) || nrow(coefficients$long %||% data.frame()) == 0L) {
    return(invisible(NULL))
  }

  cat("\nMu Coefficients\n")

  cat("---------------\n")

  cat("Coefficient Estimates\n\n")

  print(.benchmark_print_method_columns(coefficients$estimates), digits = digits, row.names = FALSE)

  long <- as.data.frame(coefficients$long, stringsAsFactors = FALSE)

  se <- .benchmark_wide_values(long, "std_error")

  lower <- .benchmark_wide_values(long, "conf.low")

  upper <- .benchmark_wide_values(long, "conf.high")

  cat("\nStandard Errors\n\n")

  print(.benchmark_print_method_columns(se), digits = digits, row.names = FALSE)

  cat("\n95% Lower Estimates\n\n")

  print(.benchmark_print_method_columns(lower), digits = digits, row.names = FALSE)

  cat("\n95% Upper Estimates\n\n")

  print(.benchmark_print_method_columns(upper), digits = digits, row.names = FALSE)

  invisible(NULL)
}
