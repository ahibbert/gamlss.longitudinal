#' List available standard longitudinal benchmark comparators
#'
#' @return A data frame describing optional comparator backends.
#' @export
benchmark_comparator_status <- function() {
  packages <- c("stats", "geepack", "lme4", "mgcv", "mgcv", "glmmTMB")
  data.frame(
    comparator = c("glm", "gee", "glmm", "gam", "gamm", "glmmTMB"),
    comparator_class = c("glm", "gee", "glmm", "gam", "gamm", "glmm"),
    estimator = c("stats::glm", "geepack::geeglm", "lme4::lmer/glmer", "mgcv::gam", "mgcv::gam + s(subject, bs = 're')", "glmmTMB::glmmTMB"),
    package = packages,
    available = vapply(packages, requireNamespace, logical(1), quietly = TRUE),
    role = c(
      "independence mean baseline",
      "marginal mean baseline with working correlation",
      "random-intercept conditional mean baseline",
      "independence smooth mean baseline",
      "smooth mean baseline with subject random-effect smooth",
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

.benchmark_metric_domain <- function(metric, estimand = NULL) {
  metric <- as.character(metric)
  estimand <- as.character(estimand %||% rep(NA_character_, length(metric)))
  out <- rep("other", length(metric))
  out[grepl("^benchmark_mean_", metric) | metric %in% c("benchmark_rmse", "benchmark_mae")] <- "mean / marginal response fit"
  out[metric %in% c(
    "benchmark_q90_mae",
    "benchmark_upper_tail_error_90",
    "benchmark_interval_coverage_95",
    "benchmark_pit_mean_abs_error",
    "benchmark_pit_ks_p_value",
    "benchmark_tail_error_lower_05",
    "benchmark_tail_error_upper_05"
  ) | estimand %in% c("quantile", "tail", "interval", "calibration")] <- "distributional prediction / shape"
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
  absolute_tolerance = 1e-8
) {
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

.benchmark_interpretation_join_metrics <- function(x) {
  if (length(x) == 0L) {
    return("none")
  }
  paste(unique(x), collapse = ", ")
}

.benchmark_interpretation_area <- function(domain) {
  out <- as.character(domain)
  out[out == "distributional prediction / shape"] <- "distributional"
  out[out == "mean / marginal response fit"] <- "mean fit"
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
  primary_method <- interpretation$primary_method
  if (length(primary_method) > 0L) {
    out <- data.frame(
      area = .benchmark_interpretation_area(leaders$domain),
      metric = leaders$metric,
      primary = ifelse(leaders$primary_leads_or_ties, "lead/tie", "-"),
      standard = ifelse(leaders$standard_leads_or_ties, "lead/tie", "-"),
      leaders = leaders$leading_methods,
      stringsAsFactors = FALSE
    )
  } else {
    out <- data.frame(
      area = .benchmark_interpretation_area(leaders$domain),
      metric = leaders$metric,
      leaders = leaders$leading_methods,
      stringsAsFactors = FALSE
    )
  }
  rownames(out) <- NULL
  out
}

.benchmark_standard_metric_columns <- function(results) {
  catalog_metrics <- .benchmark_metric_catalog()$metric
  requested <- c("elapsed_sec", "nobs", "logLik", "AIC", "mae", "rmse", catalog_metrics)
  unique(intersect(requested, names(results)))
}

.benchmark_model_status_table <- function(results) {
  results <- as.data.frame(results, stringsAsFactors = FALSE)
  columns <- c(
    "method",
    "comparator_class",
    "estimator",
    "package",
    "available",
    "success"
  )
  optional_columns <- c("warning", "error")
  optional_columns <- intersect(optional_columns, names(results))
  optional_columns <- optional_columns[vapply(results[optional_columns], function(x) any(!is.na(x)), logical(1))]
  columns <- c(columns, optional_columns)
  columns <- intersect(columns, names(results))
  out <- results[, columns, drop = FALSE]
  rownames(out) <- NULL
  out
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
    primary_method <- .benchmark_infer_primary_method(benchmark$results)
    return(list(
      results = benchmark$results,
      summary = benchmark$summary,
      scenarios = benchmark$scenarios,
      reps = benchmark$reps,
      metrics = benchmark$metrics,
      interpretation = .benchmark_interpretation(benchmark$results, primary_method = primary_method, metrics = benchmark$metrics),
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
    "## Benchmark Interpretation",
    "",
    .benchmark_interpretation_lines(inputs$interpretation),
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

.benchmark_formula_has_smooth <- function(formula) {
  term_obj <- tryCatch(
    stats::terms(stats::as.formula(formula), specials = "s"),
    error = function(e) NULL
  )
  specials <- attr(term_obj, "specials")
  length(specials$s) > 0L
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

.benchmark_empty_result_row <- function(
  status_row,
  comparator,
  available,
  success,
  elapsed_sec,
  nobs,
  warning = NA_character_,
  error = NA_character_
) {
  data.frame(
    method = comparator,
    comparator = comparator,
    comparator_class = status_row$comparator_class,
    estimator = status_row$estimator,
    package = status_row$package,
    available = available,
    success = success,
    elapsed_sec = elapsed_sec,
    nobs = nobs,
    logLik = NA_real_,
    AIC = NA_real_,
    mae = NA_real_,
    rmse = NA_real_,
    benchmark_mae = NA_real_,
    benchmark_rmse = NA_real_,
    as.list(.benchmark_distribution_metric_empty()),
    warning = warning,
    error = error,
    stringsAsFactors = FALSE
  )
}

.benchmark_distribution_metric_empty <- function() {
  c(
    benchmark_mean_bias = NA_real_,
    benchmark_mean_mae = NA_real_,
    benchmark_mean_rmse = NA_real_,
    benchmark_q90_mae = NA_real_,
    benchmark_upper_tail_error_90 = NA_real_,
    benchmark_theta_time_abs_error = NA_real_,
    benchmark_interval_coverage_95 = NA_real_,
    benchmark_pit_ks_p_value = NA_real_,
    benchmark_pit_mean_abs_error = NA_real_,
    benchmark_tail_error_lower_05 = NA_real_,
    benchmark_tail_error_upper_05 = NA_real_
  )
}

.benchmark_truth_family <- function(family, truth_family = NULL) {
  if (!is.null(truth_family)) {
    return(as.character(truth_family)[1L])
  }
  if (identical(family$family, "gaussian")) return("NO")
  if (identical(family$family, "Gamma")) return("GA")
  if (identical(family$family, "poisson")) return("PO")
  if (identical(family$family, "binomial")) return("BI")
  NA_character_
}

.benchmark_response_data <- function(data, formula) {
  response_var <- all.vars(stats::as.formula(formula))[1L]
  dat <- as.data.frame(data, stringsAsFactors = FALSE)
  dat$response <- dat[[response_var]]
  dat
}

.benchmark_predictive_distribution <- function(y, fitted, family, p = 0.9, interval_level = 0.95) {
  n <- length(fitted)
  empty <- list(
    q_p = rep(NA_real_, n),
    lower = rep(NA_real_, n),
    upper = rep(NA_real_, n),
    pit = rep(NA_real_, n)
  )
  ok <- is.finite(y) & is.finite(fitted)
  if (!any(ok)) {
    return(empty)
  }
  alpha <- (1 - interval_level) / 2

  if (identical(family$family, "gaussian")) {
    sigma_hat <- .coverage_comparator_dispersion(y, fitted, family)
    if (!is.finite(sigma_hat) || sigma_hat <= 0) return(empty)
    q_p <- fitted + stats::qnorm(p) * sigma_hat
    lower <- fitted + stats::qnorm(alpha) * sigma_hat
    upper <- fitted + stats::qnorm(1 - alpha) * sigma_hat
    pit <- stats::pnorm(y, mean = fitted, sd = sigma_hat)
  } else if (identical(family$family, "poisson")) {
    lambda <- pmax(fitted, .Machine$double.eps)
    q_p <- stats::qpois(p, lambda = lambda)
    lower <- stats::qpois(alpha, lambda = lambda)
    upper <- stats::qpois(1 - alpha, lambda = lambda)
    pit <- stats::ppois(y, lambda = lambda)
  } else if (identical(family$family, "binomial")) {
    prob <- pmin(pmax(fitted, .Machine$double.eps), 1 - .Machine$double.eps)
    q_p <- stats::qbinom(p, size = 1, prob = prob)
    lower <- stats::qbinom(alpha, size = 1, prob = prob)
    upper <- stats::qbinom(1 - alpha, size = 1, prob = prob)
    pit <- stats::pbinom(y, size = 1, prob = prob)
  } else if (identical(family$family, "Gamma")) {
    dispersion <- .coverage_comparator_dispersion(y, fitted, family)
    if (!is.finite(dispersion) || dispersion <= 0) return(empty)
    mu <- pmax(fitted, .Machine$double.eps)
    shape <- 1 / dispersion
    scale <- mu * dispersion
    q_p <- stats::qgamma(p, shape = shape, scale = scale)
    lower <- stats::qgamma(alpha, shape = shape, scale = scale)
    upper <- stats::qgamma(1 - alpha, shape = shape, scale = scale)
    pit <- stats::pgamma(y, shape = shape, scale = scale)
  } else {
    return(empty)
  }

  list(
    q_p = as.numeric(q_p),
    lower = as.numeric(lower),
    upper = as.numeric(upper),
    pit = pmin(pmax(as.numeric(pit), 0), 1)
  )
}

.benchmark_distribution_summaries <- function(y, fitted, family, data, formula, truth_family = NULL, p = 0.9, interval_level = 0.95) {
  out <- .benchmark_distribution_metric_empty()
  dat <- .benchmark_response_data(data, formula)

  if ("true_mu" %in% names(dat)) {
    ok_truth <- is.finite(dat$true_mu) & is.finite(fitted)
    if (any(ok_truth)) {
      err_mu <- fitted[ok_truth] - dat$true_mu[ok_truth]
      out["benchmark_mean_bias"] <- mean(err_mu)
      out["benchmark_mean_mae"] <- mean(abs(err_mu))
      out["benchmark_mean_rmse"] <- sqrt(mean(err_mu^2))
    }
  }

  pred_dist <- .benchmark_predictive_distribution(y, fitted, family, p = p, interval_level = interval_level)
  truth_dist <- .coverage_true_margin_distribution(
    dat,
    family = .benchmark_truth_family(family, truth_family = truth_family),
    p = p
  )

  ok_q <- is.finite(pred_dist$q_p) & is.finite(truth_dist$q)
  if (any(ok_q)) {
    out["benchmark_q90_mae"] <- mean(abs(pred_dist$q_p[ok_q] - truth_dist$q[ok_q]), na.rm = TRUE)
  }

  comparator_at_truth <- .coverage_comparator_distribution(y, fitted, family, truth_dist$q, p = p)
  ok_tail <- is.finite(comparator_at_truth$cdf_at_q) & is.finite(truth_dist$cdf)
  if (any(ok_tail)) {
    out["benchmark_upper_tail_error_90"] <- mean(
      (1 - comparator_at_truth$cdf_at_q[ok_tail]) - (1 - truth_dist$cdf[ok_tail]),
      na.rm = TRUE
    )
  }

  ok_interval <- is.finite(y) & is.finite(pred_dist$lower) & is.finite(pred_dist$upper)
  if (any(ok_interval)) {
    out["benchmark_interval_coverage_95"] <- mean(y[ok_interval] >= pred_dist$lower[ok_interval] & y[ok_interval] <= pred_dist$upper[ok_interval])
  }

  ok_pit <- is.finite(pred_dist$pit)
  if (sum(ok_pit) >= 3L) {
    pit <- pred_dist$pit[ok_pit]
    out["benchmark_pit_ks_p_value"] <- tryCatch(
      suppressWarnings(stats::ks.test(pit, "punif")$p.value),
      error = function(e) NA_real_
    )
    out["benchmark_pit_mean_abs_error"] <- abs(mean(pit, na.rm = TRUE) - 0.5)
    out["benchmark_tail_error_lower_05"] <- mean(pit <= 0.05, na.rm = TRUE) - 0.05
    out["benchmark_tail_error_upper_05"] <- mean(pit >= 0.95, na.rm = TRUE) - 0.05
  }

  out
}

.benchmark_predict_gamlss_vector <- function(fit, data, type, value_name = NULL, ...) {
  pred <- tryCatch(stats::predict(fit, newdata = data, type = type, ...), error = function(e) NULL)
  if (is.null(pred)) {
    return(rep(NA_real_, nrow(data)))
  }
  if (is.data.frame(pred)) {
    if (!is.null(value_name) && value_name %in% names(pred)) {
      pred <- pred[[value_name]]
    } else {
      pred <- pred[[ncol(pred)]]
    }
  }
  pred <- as.numeric(pred)
  if (length(pred) != nrow(data)) {
    return(rep(NA_real_, nrow(data)))
  }
  pred
}

.benchmark_gamlss_distribution_summaries <- function(fit, data, formula, truth_family = NULL, p = 0.9, interval_level = 0.95) {
  out <- .benchmark_distribution_metric_empty()
  y <- as.numeric(.benchmark_response(formula, data))
  fitted <- .benchmark_predict_response(fit, data)
  dat <- .benchmark_response_data(data, formula)

  if ("true_mu" %in% names(dat)) {
    ok_truth <- is.finite(dat$true_mu) & is.finite(fitted)
    if (any(ok_truth)) {
      err_mu <- fitted[ok_truth] - dat$true_mu[ok_truth]
      out["benchmark_mean_bias"] <- mean(err_mu)
      out["benchmark_mean_mae"] <- mean(abs(err_mu))
      out["benchmark_mean_rmse"] <- sqrt(mean(err_mu^2))
    }
  }

  alpha <- (1 - interval_level) / 2
  q_pred <- tryCatch(stats::predict(fit, newdata = data, type = "quantile", probs = c(alpha, 1 - alpha, p)), error = function(e) NULL)
  q_cols <- if (is.data.frame(q_pred) && ncol(q_pred) >= 6L) q_pred[(ncol(q_pred) - 2L):ncol(q_pred)] else NULL
  lower <- if (!is.null(q_cols)) as.numeric(q_cols[[1L]]) else rep(NA_real_, nrow(data))
  upper <- if (!is.null(q_cols)) as.numeric(q_cols[[2L]]) else rep(NA_real_, nrow(data))
  q_p <- if (!is.null(q_cols)) as.numeric(q_cols[[3L]]) else rep(NA_real_, nrow(data))

  pit <- .benchmark_predict_gamlss_vector(fit, data, type = "cdf", value_name = "cdf", q = y)
  truth_dist <- .coverage_true_margin_distribution(
    dat,
    family = .benchmark_truth_family(stats::gaussian(), truth_family = truth_family %||% fit$margin_dist$family[1]),
    p = p
  )
  cdf_at_truth_q <- .benchmark_predict_gamlss_vector(fit, data, type = "cdf", value_name = "cdf", q = truth_dist$q)

  ok_q <- is.finite(q_p) & is.finite(truth_dist$q)
  if (any(ok_q)) {
    out["benchmark_q90_mae"] <- mean(abs(q_p[ok_q] - truth_dist$q[ok_q]), na.rm = TRUE)
  }

  ok_tail <- is.finite(cdf_at_truth_q) & is.finite(truth_dist$cdf)
  if (any(ok_tail)) {
    out["benchmark_upper_tail_error_90"] <- mean((1 - cdf_at_truth_q[ok_tail]) - (1 - truth_dist$cdf[ok_tail]), na.rm = TRUE)
  }

  ok_interval <- is.finite(y) & is.finite(lower) & is.finite(upper)
  if (any(ok_interval)) {
    out["benchmark_interval_coverage_95"] <- mean(y[ok_interval] >= lower[ok_interval] & y[ok_interval] <= upper[ok_interval])
  }

  ok_pit <- is.finite(pit)
  if (sum(ok_pit) >= 3L) {
    pit <- pmin(pmax(pit[ok_pit], 0), 1)
    out["benchmark_pit_ks_p_value"] <- tryCatch(
      suppressWarnings(stats::ks.test(pit, "punif")$p.value),
      error = function(e) NA_real_
    )
    out["benchmark_pit_mean_abs_error"] <- abs(mean(pit, na.rm = TRUE) - 0.5)
    out["benchmark_tail_error_lower_05"] <- mean(pit <= 0.05, na.rm = TRUE) - 0.05
    out["benchmark_tail_error_upper_05"] <- mean(pit >= 0.95, na.rm = TRUE) - 0.05
  }

  out
}

.benchmark_fit_one <- function(
  data,
  formula,
  subject_var,
  family,
  comparator,
  correlation,
  add_subject_re_to_gamm,
  distributional_metrics = TRUE,
  truth_family = NULL,
  quantile_prob = 0.9,
  interval_level = 0.95,
  ...
) {
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
      row = .benchmark_empty_result_row(
        status_row = status_row,
        comparator = comparator,
        available = FALSE,
        success = FALSE,
        elapsed_sec = elapsed,
        nobs = nrow(data),
        error = paste0("Package '", status_row$package, "' is not installed.")
      )
    ))
  }

  smooth_comparators <- c("gam", "gamm")
  if (!comparator %in% smooth_comparators && .benchmark_formula_has_smooth(formula)) {
    elapsed <- as.numeric(difftime(Sys.time(), start, units = "secs"))
    return(list(
      fit = NULL,
      row = .benchmark_empty_result_row(
        status_row = status_row,
        comparator = comparator,
        available = TRUE,
        success = FALSE,
        elapsed_sec = elapsed,
        nobs = nrow(data),
        error = paste0(
          "Smooth terms via s(...) are only supported for 'gam' and 'gamm' comparators; ",
          "comparator '", comparator, "' cannot fit them."
        )
      )
    ))
  }

  captured <- withCallingHandlers(
    tryCatch({
      if (identical(comparator, "glm")) {
        do.call(
          stats::glm,
          c(
            list(
              formula = formula,
              data = data,
              family = family
            ),
            list(...)
          )
        )
      } else if (identical(comparator, "gee")) {
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
        mgcv::gam(formula, data = data, family = family, method = "REML", ...)
      } else if (identical(comparator, "gamm")) {
        f_gamm <- if (isTRUE(add_subject_re_to_gamm)) {
          .benchmark_formula_with_subject_re_smooth(formula, subject_var)
        } else {
          formula
        }
        mgcv::gam(f_gamm, data = data, family = family, method = "REML", ...)
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
  distribution_scores <- if (success && isTRUE(distributional_metrics)) {
    .benchmark_distribution_summaries(
      y = as.numeric(y),
      fitted = pred,
      family = family,
      data = data,
      formula = formula,
      truth_family = truth_family,
      p = quantile_prob,
      interval_level = interval_level
    )
  } else {
    .benchmark_distribution_metric_empty()
  }
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
      method = comparator,
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
      benchmark_mae = unname(scores[["mae"]]),
      benchmark_rmse = unname(scores[["rmse"]]),
      as.list(distribution_scores),
      warning = if (length(warnings)) paste(unique(warnings), collapse = " | ") else NA_character_,
      error = if (is.null(error)) NA_character_ else error,
      stringsAsFactors = FALSE
    )
  )
}

.benchmark_scalar_fit_stat <- function(expr, preferred = "joint") {
  value <- tryCatch(expr, error = function(e) NA_real_)
  if (is.data.frame(value)) {
    if (preferred %in% names(value)) {
      value <- value[[preferred]]
    } else {
      value <- value[[1L]]
    }
  } else if (is.matrix(value)) {
    if (preferred %in% colnames(value)) {
      value <- value[, preferred]
    } else if (preferred %in% rownames(value)) {
      value <- value[preferred, ]
    } else {
      value <- value[[1L]]
    }
  }
  if (length(value) == 0L) {
    return(NA_real_)
  }
  if (!is.null(names(value)) && preferred %in% names(value)) {
    value <- value[[preferred]]
  } else {
    value <- value[[1L]]
  }
  value <- suppressWarnings(as.numeric(value))
  if (length(value) == 0L || !is.finite(value[[1L]])) {
    return(NA_real_)
  }
  value[[1L]]
}

.benchmark_supplied_fit_aic <- function(fit) {
  aic <- .benchmark_scalar_fit_stat(stats::AIC(fit))
  if (is.finite(aic)) {
    return(aic)
  }
  if (inherits(fit, "gamlss.longitudinal")) {
    summary_obj <- tryCatch(summary(fit, include_vcov = FALSE), error = function(e) NULL)
    value <- summary_obj$fit$AIC %||% NA_real_
    value <- suppressWarnings(as.numeric(value))
    if (length(value) > 0L && is.finite(value[[1L]])) {
      return(value[[1L]])
    }
  }
  NA_real_
}

.benchmark_supplied_fit_one <- function(
  fit,
  data,
  formula,
  family,
  fit_name,
  distributional_metrics = TRUE,
  truth_family = NULL,
  quantile_prob = 0.9,
  interval_level = 0.95
) {
  if (is.null(fit)) {
    return(NULL)
  }
  fit_name <- as.character(fit_name)[1L]
  class_name <- class(fit)[1L]
  comparator_class <- if (inherits(fit, "gamlss.longitudinal")) {
    "gamlss_longitudinal"
  } else {
    class_name
  }
  estimator <- if (inherits(fit, "gamlss.longitudinal")) {
    "gamlss.longitudinal::gamlss_longitudinal"
  } else {
    class_name
  }
  package <- if (inherits(fit, "gamlss.longitudinal")) {
    "gamlss.longitudinal"
  } else {
    NA_character_
  }

  y <- .benchmark_response(formula, data)
  pred <- .benchmark_predict_response(fit, data)
  scores <- .benchmark_score_predictions(as.numeric(y), pred)
  success <- any(is.finite(pred))
  distribution_scores <- if (success && isTRUE(distributional_metrics) && inherits(fit, "gamlss.longitudinal")) {
    .benchmark_gamlss_distribution_summaries(
      fit = fit,
      data = data,
      formula = formula,
      truth_family = truth_family,
      p = quantile_prob,
      interval_level = interval_level
    )
  } else if (success && isTRUE(distributional_metrics)) {
    .benchmark_distribution_summaries(
      y = as.numeric(y),
      fitted = pred,
      family = family,
      data = data,
      formula = formula,
      truth_family = truth_family,
      p = quantile_prob,
      interval_level = interval_level
    )
  } else {
    .benchmark_distribution_metric_empty()
  }
  ll <- if (success) .benchmark_scalar_fit_stat(stats::logLik(fit)) else NA_real_
  aic <- if (success) .benchmark_supplied_fit_aic(fit) else NA_real_

  list(
    fit = fit,
    row = data.frame(
      method = fit_name,
      comparator = fit_name,
      comparator_class = comparator_class,
      estimator = estimator,
      package = package,
      available = TRUE,
      success = success,
      elapsed_sec = NA_real_,
      nobs = nrow(data),
      logLik = ll,
      AIC = aic,
      mae = unname(scores[["mae"]]),
      rmse = unname(scores[["rmse"]]),
      benchmark_mae = unname(scores[["mae"]]),
      benchmark_rmse = unname(scores[["rmse"]]),
      as.list(distribution_scores),
      warning = NA_character_,
      error = if (success) NA_character_ else "Prediction failed or returned no finite response-scale values.",
      stringsAsFactors = FALSE
    )
  )
}

#' Fit standard longitudinal comparator models for adoption benchmarks
#'
#' `benchmark_standard_models()` is an opt-in scaffold for comparing
#' `gamlss.longitudinal` with models users already know. It is intentionally
#' narrow: optionally score an already-fitted primary model, fit common
#' mean-model baselines, record whether they ran, and return timing,
#' response-scale prediction metrics, and distributional diagnostics when the
#' fitted model family supports them.
#'
#' @param data Long-format data frame.
#' @param formula Mean-model formula used by the comparator models.
#' @param subject_var Subject identifier column.
#' @param family A base R family object or one of `"gaussian"`, `"poisson"`,
#'   `"binomial"`, or `"gamma"`.
#' @param comparators Character vector containing any of `"glm"`, `"gee"`,
#'   `"glmm"`, `"gam"`, `"gamm"`, and `"glmmTMB"`.
#' @param correlation Working correlation passed to `geepack::geeglm()`.
#' @param add_subject_re_to_gamm Logical; add `s(subject, bs = "re")` to the
#'   GAMM comparator.
#' @param add_subject_re_to_gam Deprecated alias for `add_subject_re_to_gamm`.
#' @param fit Optional already-fitted primary model to score beside the standard
#'   comparators. For a `gamlss.longitudinal` fit, response-scale predictions
#'   are obtained with `predict(fit, newdata = data, type = "response")`.
#' @param fit_name Label used for the supplied `fit` row.
#' @param distributional_metrics Logical; compute interval coverage, PIT
#'   calibration, tail-frequency diagnostics, and truth-aware quantile/tail
#'   metrics when enough information is available.
#' @param truth_family Optional `gamlss.dist` family name used to derive true
#'   quantiles from `true_*` columns. Defaults are inferred from `family`.
#' @param quantile_prob Probability used for the benchmark quantile and upper
#'   tail metrics.
#' @param interval_level Prediction interval level used for empirical coverage.
#' @param ... Additional arguments passed to each comparator fit.
#'
#' @return An object of class `gamlss_longitudinal_benchmark` with `results`,
#'   `fits`, and `interpretation` components.
#' @export
benchmark_standard_models <- function(
  data,
  formula,
  subject_var,
  family = "gaussian",
  comparators = c("glm", "gee", "glmm", "gam", "gamm"),
  correlation = "exchangeable",
  add_subject_re_to_gamm = TRUE,
  add_subject_re_to_gam = NULL,
  fit = NULL,
  fit_name = "gamlss.longitudinal",
  distributional_metrics = TRUE,
  truth_family = NULL,
  quantile_prob = 0.9,
  interval_level = 0.95,
  ...
) {
  data <- as.data.frame(data, stringsAsFactors = FALSE)
  if (!is.character(subject_var) || length(subject_var) != 1L || !subject_var %in% names(data)) {
    stop("'subject_var' must be a single column name in 'data'.", call. = FALSE)
  }
  if (!is.null(add_subject_re_to_gam)) {
    warning(
      "'add_subject_re_to_gam' is deprecated; use 'add_subject_re_to_gamm' instead.",
      call. = FALSE
    )
    add_subject_re_to_gamm <- add_subject_re_to_gam
  }
  formula <- stats::as.formula(formula)
  family <- .benchmark_family(family)
  quantile_prob <- as.numeric(quantile_prob)[1L]
  interval_level <- as.numeric(interval_level)[1L]
  if (!is.finite(quantile_prob) || quantile_prob <= 0 || quantile_prob >= 1) {
    stop("'quantile_prob' must be a probability strictly between 0 and 1.", call. = FALSE)
  }
  if (!is.finite(interval_level) || interval_level <= 0 || interval_level >= 1) {
    stop("'interval_level' must be a probability strictly between 0 and 1.", call. = FALSE)
  }
  comparators <- unique(as.character(comparators))
  valid <- benchmark_comparator_status()$comparator
  bad <- setdiff(comparators, valid)
  if (length(bad) > 0L) {
    stop("Unknown comparator(s): ", paste(bad, collapse = ", "), call. = FALSE)
  }
  if (!is.null(fit)) {
    if (!is.character(fit_name) || length(fit_name) != 1L || is.na(fit_name) || !nzchar(fit_name)) {
      stop("'fit_name' must be a single non-empty character value.", call. = FALSE)
    }
    if (fit_name %in% comparators) {
      stop("'fit_name' must not duplicate a standard comparator name.", call. = FALSE)
    }
  }

  comparator_runs <- lapply(comparators, function(comparator) {
    .benchmark_fit_one(
      data = data,
      formula = formula,
      subject_var = subject_var,
      family = family,
      comparator = comparator,
      correlation = correlation,
      add_subject_re_to_gamm = add_subject_re_to_gamm,
      distributional_metrics = distributional_metrics,
      truth_family = truth_family,
      quantile_prob = quantile_prob,
      interval_level = interval_level,
      ...
    )
  })
  supplied_run <- .benchmark_supplied_fit_one(
    fit = fit,
    data = data,
    formula = formula,
    family = family,
    fit_name = fit_name,
    distributional_metrics = distributional_metrics,
    truth_family = truth_family,
    quantile_prob = quantile_prob,
    interval_level = interval_level
  )
  runs <- c(if (!is.null(supplied_run)) list(supplied_run), comparator_runs)
  fit_names <- c(if (!is.null(supplied_run)) fit_name, comparators)
  results <- do.call(rbind, lapply(runs, `[[`, "row"))
  primary_method <- if (!is.null(supplied_run)) fit_name else NULL
  out <- list(
    results = results,
    fits = stats::setNames(lapply(runs, `[[`, "fit"), fit_names),
    formula = formula,
    subject_var = subject_var,
    family = family$family,
    primary_method = primary_method,
    interpretation = .benchmark_interpretation(results, primary_method = primary_method)
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

  cat("Model Status\n")
  cat("------------\n")
  print(.benchmark_model_status_table(x$results), digits = digits, row.names = FALSE)

  metric_table <- .benchmark_metric_matrix(x$results)
  if (nrow(metric_table) > 0L) {
    cat("\nBenchmark Metrics\n")
    cat("-----------------\n")
    print(metric_table, digits = digits, row.names = FALSE)
  }

  interpretation_table <- .benchmark_interpretation_table(x$interpretation)
  if (nrow(interpretation_table) > 0L) {
    cat("\nBenchmark Interpretation\n")
    cat("------------------------\n")
    if (length(x$interpretation$primary_method) > 0L) {
      cat("Primary method:", x$interpretation$primary_method[[1L]], "\n\n")
    }
    print(interpretation_table, digits = digits, row.names = FALSE)
  } else {
    lines <- .benchmark_interpretation_lines(x$interpretation)
    if (length(lines) > 0L) {
      cat("\nBenchmark Interpretation\n")
      cat("------------------------\n")
      cat(paste(lines, collapse = "\n"), "\n", sep = "")
    }
  }

  unsupported <- if (!is.null(x$interpretation)) {
    as.data.frame(x$interpretation$unsupported, stringsAsFactors = FALSE)
  } else {
    data.frame()
  }
  if (nrow(unsupported) > 0L) {
    cat("\nRows Excluded From Interpretation\n")
    cat("---------------------------------\n")
    print(unsupported, digits = digits, row.names = FALSE)
  }
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
      "smooth_eta_rmse",
      "smooth_eta_max_abs_error",
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
      "smooth",
      "smooth",
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
      "lower",
      "lower",
      "lower"
    ),
    target = c(NA, NA, 0, NA, NA, NA, 0, NA, 0.95, NA, 1, 0, 0, NA, NA, NA),
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
