source(file.path("paper", "R", "08-simulation-sensitivity-correlation-misspecification", "standard-model-benchmarking", "00-benchmark-setup.R"))

source_used <- bmk_load_package()
bmk_require_namespaces(c("gamlss.dist", "mvtnorm"), strict = TRUE)

base_run_dir <- bmk_env(
  "GAMLSS_LONGITUDINAL_SANDWICH_BASE_RUN_DIR",
  file.path(bmk_output_root, "run_20260619_t20_t50_combined")
)
base_run_dir <- normalizePath(base_run_dir, winslash = "/", mustWork = TRUE)

out_dir <- bmk_env(
  "GAMLSS_LONGITUDINAL_SANDWICH_OUT_DIR",
  file.path(base_run_dir, "sandwich_t20_grid")
)
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

scenario_names <- bmk_env_vector(
  "GAMLSS_LONGITUDINAL_SANDWICH_SCENARIOS",
  c(
    "external_exchangeable_moderate_t20",
    "external_exchangeable_high_t20",
    "external_ar1_moderate_t20",
    "external_ar1_high_t20",
    "internal_time_varying_high_t20",
    "internal_covariate_dependent_high_t20"
  )
)
family_names <- bmk_env_vector(
  "GAMLSS_LONGITUDINAL_SANDWICH_FAMILIES",
  c("gaussian", "gamma", "binary")
)
rep_limit <- bmk_env_int("GAMLSS_LONGITUDINAL_SANDWICH_REP_LIMIT", 20L)
n <- bmk_env_int("GAMLSS_LONGITUDINAL_BENCHMARK_N", 120L)
interval_level <- bmk_env_num("GAMLSS_LONGITUDINAL_BENCHMARK_INTERVAL_LEVEL", 0.95)
max_elapsed_sec <- bmk_env_num("GAMLSS_LONGITUDINAL_BENCHMARK_MAX_ELAPSED_SEC", 180)
max_outer_iter <- bmk_env_int("GAMLSS_LONGITUDINAL_BENCHMARK_MAX_OUTER_ITER", 100L)
max_inner_iter <- bmk_env_int("GAMLSS_LONGITUDINAL_BENCHMARK_MAX_INNER_ITER", 100L)
sandwich_h <- bmk_env_num("GAMLSS_LONGITUDINAL_SANDWICH_H", 1e-4)
sandwich_bread_method <- bmk_env("GAMLSS_LONGITUDINAL_SANDWICH_BREAD_METHOD", "analytical")
rerun_failures <- bmk_env_flag("GAMLSS_LONGITUDINAL_SANDWICH_RERUN_FAILURES", FALSE)

required_paths <- c(
  status = file.path(base_run_dir, "primary_status_by_rep.csv"),
  coef = file.path(base_run_dir, "coefficient_results_by_rep.csv"),
  benchmark_summary = file.path(base_run_dir, "benchmark_summary.csv"),
  se_summary = file.path(base_run_dir, "se_calibration_summary.csv"),
  dependence_summary = file.path(base_run_dir, "dependence_recovery_summary.csv"),
  complexity_summary = file.path(base_run_dir, "fit_complexity_summary.csv"),
  review = file.path(base_run_dir, "full_summary_review_table.csv"),
  status_summary = file.path(base_run_dir, "primary_status_summary.csv")
)
missing_paths <- required_paths[!file.exists(required_paths)]
if (length(missing_paths) > 0L) {
  stop("Missing required benchmark output(s): ", paste(missing_paths, collapse = ", "), call. = FALSE)
}

status <- read.csv(required_paths[["status"]], stringsAsFactors = FALSE, check.names = FALSE)
existing_coef <- read.csv(required_paths[["coef"]], stringsAsFactors = FALSE, check.names = FALSE)
benchmark_summary <- read.csv(required_paths[["benchmark_summary"]], stringsAsFactors = FALSE, check.names = FALSE)
se_summary <- read.csv(required_paths[["se_summary"]], stringsAsFactors = FALSE, check.names = FALSE)
dependence_summary <- read.csv(required_paths[["dependence_summary"]], stringsAsFactors = FALSE, check.names = FALSE)
complexity_summary <- read.csv(required_paths[["complexity_summary"]], stringsAsFactors = FALSE, check.names = FALSE)
review <- read.csv(required_paths[["review"]], stringsAsFactors = FALSE, check.names = FALSE)
status_summary <- read.csv(required_paths[["status_summary"]], stringsAsFactors = FALSE, check.names = FALSE)

ok <- function(x) x %in% c(TRUE, "TRUE", "True", "true", "1")

scenario_base <- sub("_t[0-9]+$", "", scenario_names)
timepoints <- suppressWarnings(as.integer(sub("^.*_t([0-9]+)$", "\\1", scenario_names)))
if (any(!is.finite(timepoints)) || any(timepoints != 20L)) {
  stop("This runner is restricted to scenario names with a '_t20' suffix.", call. = FALSE)
}
base_scenarios <- bmk_scenario_specs()
scenarios <- stats::setNames(vector("list", length(scenario_names)), scenario_names)
for (i in seq_along(scenario_names)) {
  if (!scenario_base[[i]] %in% names(base_scenarios)) {
    stop("Unknown scenario: ", scenario_names[[i]], call. = FALSE)
  }
  scenarios[[i]] <- bmk_resize_scenario_time(base_scenarios[[scenario_base[[i]]]], timepoints[[i]])
}
families <- bmk_select_named(bmk_family_specs(include_binary = TRUE), family_names)
if ("poisson" %in% names(families)) {
  stop("Poisson is intentionally excluded from the paper-facing sandwich grid.", call. = FALSE)
}

extract_sandwich_coefs <- function(fit, spec, level = 0.95) {
  s <- summary(
    fit,
    include_vcov = TRUE,
    vcov_method = "sandwich",
    sandwich_h = sandwich_h,
    sandwich_bread_method = sandwich_bread_method
  )
  tbl <- as.data.frame(s$coefficients, stringsAsFactors = FALSE)
  tbl <- tbl[tbl$parameter == "mu", , drop = FALSE]
  if (nrow(tbl) == 0L) return(data.frame())
  z <- stats::qnorm((1 + level) / 2)
  out <- data.frame(
    method = "rs_joint_sandwich",
    parameter = tbl$parameter,
    term = sub("^mu\\.", "", tbl$term),
    estimate = tbl$estimate,
    std_error = tbl$std_error,
    conf.low = tbl$estimate - z * tbl$std_error,
    conf.high = tbl$estimate + z * tbl$std_error,
    stringsAsFactors = FALSE
  )
  out$term <- gamlss.longitudinal:::.benchmark_normalize_coef_term(out$term, time_var = fit$time_var)
  bmk_annotate_coefficients(out, spec, level = level)
}

normalise_scenario_label <- function(x) {
  if (!nrow(x) || !"scenario" %in% names(x) || !"n_time" %in% names(x)) return(x)
  has_suffix <- grepl("_t[0-9]+$", x$scenario)
  x$scenario[!has_suffix] <- paste0(x$scenario[!has_suffix], "_t", x$n_time[!has_suffix])
  x
}

read_checkpoint_csv <- function(path) {
  if (!file.exists(path) || file.info(path)$size <= 4L) return(data.frame())
  tryCatch(
    read.csv(path, stringsAsFactors = FALSE, check.names = FALSE),
    error = function(e) data.frame()
  )
}

write_checkpoint_csv <- function(x, path) {
  if (!is.data.frame(x) || nrow(x) == 0L) return(invisible(path))
  bmk_write_csv(x, path)
}

case_key <- function(x) paste(x$scenario, x$family, x$rep, sep = "::")

candidate_cases <- status[
  status$method == "rs_joint" &
    ok(status$success) &
    status$scenario %in% scenario_names &
    status$family %in% family_names,
  ,
  drop = FALSE
]
ordered_cases <- list()
for (scenario_name in scenario_names) {
  for (family_name in family_names) {
    cell <- candidate_cases[
      candidate_cases$scenario == scenario_name &
        candidate_cases$family == family_name,
      ,
      drop = FALSE
    ]
    cell <- cell[order(cell$rep), , drop = FALSE]
    ordered_cases[[length(ordered_cases) + 1L]] <- utils::head(cell, rep_limit)
  }
}
candidate_cases <- bmk_bind_rows_fill(ordered_cases)
rownames(candidate_cases) <- NULL
if (nrow(candidate_cases) == 0L) stop("No candidate rs_joint cases found.", call. = FALSE)

coef_out_path <- file.path(out_dir, "rs_joint_sandwich_coefficients_by_rep.csv")
status_out_path <- file.path(out_dir, "rs_joint_sandwich_status_by_rep.csv")
complexity_out_path <- file.path(out_dir, "rs_joint_sandwich_complexity_by_rep.csv")
dependence_out_path <- file.path(out_dir, "rs_joint_sandwich_dependence_by_rep.csv")

sandwich_coef <- read_checkpoint_csv(coef_out_path)
sandwich_status <- read_checkpoint_csv(status_out_path)
sandwich_complexity <- read_checkpoint_csv(complexity_out_path)
sandwich_dependence <- read_checkpoint_csv(dependence_out_path)

completed_keys <- character()
if (nrow(sandwich_status) > 0L) {
  if (isTRUE(rerun_failures)) {
    completed_keys <- case_key(sandwich_status[ok(sandwich_status$success), , drop = FALSE])
  } else {
    completed_keys <- case_key(sandwich_status)
  }
}

flush_outputs <- function() {
  write_checkpoint_csv(normalise_scenario_label(sandwich_coef), coef_out_path)
  write_checkpoint_csv(normalise_scenario_label(sandwich_status), status_out_path)
  write_checkpoint_csv(normalise_scenario_label(sandwich_complexity), complexity_out_path)
  write_checkpoint_csv(normalise_scenario_label(sandwich_dependence), dependence_out_path)
}

for (case_idx in seq_len(nrow(candidate_cases))) {
  row <- candidate_cases[case_idx, , drop = FALSE]
  key <- case_key(row)
  if (key %in% completed_keys) {
    message("[", case_idx, "/", nrow(candidate_cases), "] skipping completed ", row$scenario, " / ", row$family, " / rep ", row$rep)
    next
  }

  scenario <- scenarios[[row$scenario]]
  scenario$scenario <- row$scenario
  spec <- families[[row$family]]
  message("[", case_idx, "/", nrow(candidate_cases), "] ", row$scenario, " / ", row$family, " / rep ", row$rep)

  dat <- bmk_simulate_dataset(spec, scenario, n = n, seed = row$seed)
  fit_start <- Sys.time()
  fit_capture <- bmk_capture_fit(
    bmk_fit_rs_joint(
      dat = dat,
      spec = spec,
      scenario = scenario,
      max_elapsed_sec = max_elapsed_sec,
      max_outer_iter = max_outer_iter,
      max_inner_iter = max_inner_iter
    )
  )
  fit_elapsed <- as.numeric(difftime(Sys.time(), fit_start, units = "secs"))
  fit_success <- inherits(fit_capture$value, "gamlss.longitudinal")
  sandwich_success <- FALSE
  sandwich_elapsed <- NA_real_
  error <- if (fit_success) NA_character_ else conditionMessage(fit_capture$value)
  warning <- if (length(fit_capture$warnings)) paste(fit_capture$warnings, collapse = " | ") else NA_character_

  if (fit_success) {
    sandwich_start <- Sys.time()
    sandwich_capture <- bmk_capture_fit(extract_sandwich_coefs(fit_capture$value, spec, level = interval_level))
    sandwich_elapsed <- as.numeric(difftime(Sys.time(), sandwich_start, units = "secs"))
    sandwich_success <- !inherits(sandwich_capture$value, "error")
    if (length(sandwich_capture$warnings)) {
      warning <- paste(c(warning, sandwich_capture$warnings), collapse = " | ")
      warning <- sub("^NA \\| ", "", warning)
    }
    if (!sandwich_success) {
      error <- conditionMessage(sandwich_capture$value)
    } else {
      coefs <- bmk_prefix_coefficients(
        sandwich_capture$value,
        scenario = scenario,
        spec = spec,
        rep_id = row$rep,
        seed = row$seed
      )
      coefs$method <- "rs_joint_sandwich"
      sandwich_coef <- bmk_bind_rows_fill(list(sandwich_coef, coefs))

      complexity <- bmk_fit_complexity_row(fit_capture$value, "rs_joint_sandwich", scenario, spec, row$rep, row$seed)
      dependence <- bmk_dependence_recovery_row(dat, fit_capture$value, "rs_joint_sandwich", scenario, spec, row$rep, row$seed)
      sandwich_complexity <- bmk_bind_rows_fill(list(sandwich_complexity, complexity))
      sandwich_dependence <- bmk_bind_rows_fill(list(sandwich_dependence, dependence))
    }
  }

  status_row <- data.frame(
    generator = row$generator,
    scenario = row$scenario,
    correlation = row$correlation,
    correlation_level = row$correlation_level,
    n_time = row$n_time,
    family = row$family,
    gamlss_family = row$gamlss_family,
    rep = row$rep,
    seed = row$seed,
    method = "rs_joint_sandwich",
    success = fit_success && sandwich_success,
    fit_success = fit_success,
    sandwich_success = sandwich_success,
    elapsed_sec = fit_elapsed + ifelse(is.finite(sandwich_elapsed), sandwich_elapsed, 0),
    fit_elapsed_sec = fit_elapsed,
    sandwich_elapsed_sec = sandwich_elapsed,
    error = error,
    warning = warning,
    stringsAsFactors = FALSE
  )
  sandwich_status <- bmk_bind_rows_fill(list(sandwich_status, status_row))
  completed_keys <- unique(c(completed_keys, key))
  flush_outputs()
}

sandwich_coef <- normalise_scenario_label(sandwich_coef)
sandwich_status <- normalise_scenario_label(sandwich_status)
sandwich_complexity <- normalise_scenario_label(sandwich_complexity)
sandwich_dependence <- normalise_scenario_label(sandwich_dependence)
flush_outputs()

to_numeric <- function(x, cols) {
  for (nm in intersect(cols, names(x))) x[[nm]] <- suppressWarnings(as.numeric(x[[nm]]))
  x
}
to_logical_cols <- function(x, cols) {
  for (nm in intersect(cols, names(x))) x[[nm]] <- ok(x[[nm]])
  x
}

make_coef_summary <- function(coef) {
  coef <- to_numeric(coef, c("estimate", "std_error", "conf.low", "conf.high", "truth", "bias", "ci_width", "p_value"))
  coef <- to_logical_cols(coef, c("ci_covers_truth", "false_positive"))
  group_cols <- c("generator", "scenario", "correlation", "correlation_level", "n_time", "family", "method", "term")
  if (nrow(coef) == 0L) return(data.frame())
  groups <- unique(coef[group_cols])
  rows <- vector("list", nrow(groups))
  for (i in seq_len(nrow(groups))) {
    idx <- rep(TRUE, nrow(coef))
    for (col in group_cols) idx <- idx & coef[[col]] == groups[[col]][i]
    sub <- coef[idx, , drop = FALSE]
    est <- sub$estimate
    se <- sub$std_error
    truth <- sub$truth[is.finite(sub$truth)][1L]
    empirical_sd <- stats::sd(est, na.rm = TRUE)
    mean_reported_se <- mean(se, na.rm = TRUE)
    rows[[i]] <- cbind(
      groups[i, , drop = FALSE],
      data.frame(
        truth = truth,
        n = nrow(sub),
        n_estimate = sum(is.finite(est)),
        mean_estimate = mean(est, na.rm = TRUE),
        bias = mean(est - sub$truth, na.rm = TRUE),
        rmse = sqrt(mean((est - sub$truth)^2, na.rm = TRUE)),
        empirical_sd = empirical_sd,
        mean_reported_se = mean_reported_se,
        se_calibration_ratio = mean_reported_se / empirical_sd,
        ci_coverage = mean(sub$ci_covers_truth, na.rm = TRUE),
        median_ci_width = stats::median(sub$ci_width, na.rm = TRUE),
        false_positive_rate = mean(sub$false_positive, na.rm = TRUE),
        stringsAsFactors = FALSE
      )
    )
  }
  bmk_bind_rows_fill(rows)
}

make_status_summary <- function(status) {
  if (nrow(status) == 0L) return(data.frame())
  group_cols <- c("generator", "scenario", "correlation", "correlation_level", "n_time", "family", "method")
  groups <- unique(status[group_cols])
  rows <- vector("list", nrow(groups))
  for (i in seq_len(nrow(groups))) {
    idx <- rep(TRUE, nrow(status))
    for (col in group_cols) idx <- idx & status[[col]] == groups[[col]][i]
    sub <- status[idx, , drop = FALSE]
    rows[[i]] <- cbind(
      groups[i, , drop = FALSE],
      data.frame(
        n = nrow(sub),
        success_rate = mean(ok(sub$success)),
        convergence_rate = NA_real_,
        median_elapsed_sec = stats::median(suppressWarnings(as.numeric(sub$elapsed_sec)), na.rm = TRUE),
        stringsAsFactors = FALSE
      )
    )
  }
  bmk_bind_rows_fill(rows)
}

sandwich_coef_summary <- make_coef_summary(sandwich_coef)
sandwich_status_summary <- make_status_summary(sandwich_status)

selected <- function(x) {
  x$scenario %in% scenario_names & x$family %in% family_names & x$n_time == 20L
}
drop_sandwich <- function(x) {
  if (!nrow(x) || !"method" %in% names(x)) return(x)
  x[x$method != "rs_joint_sandwich", , drop = FALSE]
}

augmented_benchmark_summary <- drop_sandwich(benchmark_summary)
benchmark_copy <- benchmark_summary[selected(benchmark_summary) & benchmark_summary$method == "rs_joint", , drop = FALSE]
benchmark_copy$method <- "rs_joint_sandwich"
augmented_benchmark_summary <- bmk_bind_rows_fill(list(augmented_benchmark_summary, benchmark_copy))

augmented_se_summary <- bmk_bind_rows_fill(list(drop_sandwich(se_summary), sandwich_coef_summary))

augmented_dependence_summary <- drop_sandwich(dependence_summary)
dependence_copy <- dependence_summary[selected(dependence_summary) & dependence_summary$method == "rs_joint", , drop = FALSE]
dependence_copy$method <- "rs_joint_sandwich"
augmented_dependence_summary <- bmk_bind_rows_fill(list(augmented_dependence_summary, dependence_copy))

augmented_complexity_summary <- drop_sandwich(complexity_summary)
complexity_copy <- complexity_summary[selected(complexity_summary) & complexity_summary$method == "rs_joint", , drop = FALSE]
complexity_copy$method <- "rs_joint_sandwich"
augmented_complexity_summary <- bmk_bind_rows_fill(list(augmented_complexity_summary, complexity_copy))

augmented_status_summary <- bmk_bind_rows_fill(list(drop_sandwich(status_summary), sandwich_status_summary))

augmented_review <- drop_sandwich(review)
review_copy <- review[selected(review) & review$method == "rs_joint", , drop = FALSE]
review_copy$method <- "rs_joint_sandwich"
x_summary <- sandwich_coef_summary[sandwich_coef_summary$term == "x", , drop = FALSE]
null_summary <- sandwich_coef_summary[sandwich_coef_summary$term == "z_null", , drop = FALSE]
summary_key <- function(x) paste(x$n_time, x$scenario, x$family, x$method, sep = "::")
review_key <- summary_key(review_copy)
x_key <- summary_key(x_summary)
null_key <- summary_key(null_summary)
review_copy$x_se_ratio <- x_summary$se_calibration_ratio[match(review_key, x_key)]
review_copy$x_ci_coverage <- x_summary$ci_coverage[match(review_key, x_key)]
review_copy$null_se_ratio <- null_summary$se_calibration_ratio[match(review_key, null_key)]
review_copy$null_ci_coverage <- null_summary$ci_coverage[match(review_key, null_key)]
review_copy$null_false_positive <- null_summary$false_positive_rate[match(review_key, null_key)]
augmented_review <- bmk_bind_rows_fill(list(augmented_review, review_copy[names(augmented_review)]))

existing_compare <- existing_coef[
  existing_coef$method == "rs_joint" &
    existing_coef$scenario %in% scenario_names &
    existing_coef$family %in% family_names &
    existing_coef$n_time == 20L,
  ,
  drop = FALSE
]
coef_compare <- merge(
  sandwich_coef[, c("scenario", "family", "rep", "seed", "term", "estimate", "std_error"), drop = FALSE],
  existing_compare[, c("scenario", "family", "rep", "seed", "term", "estimate", "std_error"), drop = FALSE],
  by = c("scenario", "family", "rep", "seed", "term"),
  all.x = TRUE,
  suffixes = c("_sandwich", "_model")
)
coef_compare <- to_numeric(coef_compare, c("estimate_sandwich", "estimate_model", "std_error_sandwich", "std_error_model"))

validation <- data.frame(
  expected_cases = nrow(candidate_cases),
  completed_status_cases = nrow(sandwich_status),
  successful_cases = sum(ok(sandwich_status$success)),
  failed_cases = sum(!ok(sandwich_status$success)),
  coefficient_rows = nrow(sandwich_coef),
  nonfinite_sandwich_se = sum(!is.finite(suppressWarnings(as.numeric(sandwich_coef$std_error)))),
  max_abs_estimate_difference_vs_original = max(abs(coef_compare$estimate_sandwich - coef_compare$estimate_model), na.rm = TRUE),
  median_elapsed_sec = stats::median(suppressWarnings(as.numeric(sandwich_status$elapsed_sec)), na.rm = TRUE),
  median_fit_elapsed_sec = stats::median(suppressWarnings(as.numeric(sandwich_status$fit_elapsed_sec)), na.rm = TRUE),
  median_sandwich_elapsed_sec = stats::median(suppressWarnings(as.numeric(sandwich_status$sandwich_elapsed_sec)), na.rm = TRUE),
  stringsAsFactors = FALSE
)

bmk_write_csv(augmented_benchmark_summary, file.path(out_dir, "augmented_benchmark_summary.csv"))
bmk_write_csv(augmented_se_summary, file.path(out_dir, "augmented_se_calibration_summary.csv"))
bmk_write_csv(augmented_dependence_summary, file.path(out_dir, "augmented_dependence_recovery_summary.csv"))
bmk_write_csv(augmented_complexity_summary, file.path(out_dir, "augmented_fit_complexity_summary.csv"))
bmk_write_csv(augmented_status_summary, file.path(out_dir, "augmented_primary_status_summary.csv"))
bmk_write_csv(augmented_review, file.path(out_dir, "augmented_full_summary_review_table.csv"))
bmk_write_csv(sandwich_coef_summary, file.path(out_dir, "rs_joint_sandwich_se_calibration_summary.csv"))
bmk_write_csv(sandwich_status_summary, file.path(out_dir, "rs_joint_sandwich_status_summary.csv"))
bmk_write_csv(coef_compare, file.path(out_dir, "rs_joint_sandwich_estimate_validation.csv"))
bmk_write_csv(validation, file.path(out_dir, "rs_joint_sandwich_validation_summary.csv"))

message("Wrote T=20 sandwich grid outputs to: ", out_dir)
print(validation, row.names = FALSE)
