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

      c("benchmark_mean_rmse", "benchmark_neg_log_score", "benchmark_q90_mae", "benchmark_upper_tail_error_90", "benchmark_interval_coverage_95", "benchmark_interval_width_95", "elapsed_sec"),

      c("benchmark_mean_rmse", "benchmark_neg_log_score", "benchmark_q90_mae", "benchmark_upper_tail_error_90", "elapsed_sec"),

      c("benchmark_mean_rmse", "benchmark_neg_log_score", "benchmark_q90_mae", "benchmark_upper_tail_error_90", "elapsed_sec"),

      c("benchmark_theta_time_abs_error", "elapsed_sec"),

      c("benchmark_mean_rmse", "benchmark_neg_log_score", "benchmark_interval_coverage_95", "benchmark_interval_width_95", "elapsed_sec")

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


