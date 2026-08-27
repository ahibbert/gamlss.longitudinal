source(file.path("paper", "R", "08-simulation-sensitivity-correlation-misspecification", "standard-model-benchmarking", "00-benchmark-setup.R"))

source_used <- bmk_load_package()
bmk_require_namespaces(c("gamlss.dist", "mvtnorm"), strict = TRUE)

dir.create(bmk_output_root, recursive = TRUE, showWarnings = FALSE)

include_binary <- bmk_env_flag("GAMLSS_LONGITUDINAL_BENCHMARK_INCLUDE_BINARY", TRUE)
families <- bmk_select_named(
  bmk_family_specs(include_binary = include_binary),
  bmk_env_vector("GAMLSS_LONGITUDINAL_BENCHMARK_FAMILIES", character())
)
scenarios <- bmk_select_named(
  bmk_scenario_specs(),
  bmk_env_vector("GAMLSS_LONGITUDINAL_BENCHMARK_SCENARIOS", character())
)
timepoints <- bmk_env_int_vector("GAMLSS_LONGITUDINAL_BENCHMARK_TIMEPOINTS", integer())
scenarios <- bmk_expand_scenarios_for_timepoints(scenarios, timepoints)

reps <- bmk_env_int("GAMLSS_LONGITUDINAL_BENCHMARK_REPS", 20L)
rep_start <- bmk_env_int("GAMLSS_LONGITUDINAL_BENCHMARK_REP_START", 1L)
rep_end <- bmk_env_int("GAMLSS_LONGITUDINAL_BENCHMARK_REP_END", reps)
if (rep_start < 1L || rep_end < rep_start || rep_end > reps) {
  stop(
    "Invalid replicate range: require 1 <= GAMLSS_LONGITUDINAL_BENCHMARK_REP_START <= ",
    "GAMLSS_LONGITUDINAL_BENCHMARK_REP_END <= GAMLSS_LONGITUDINAL_BENCHMARK_REPS.",
    call. = FALSE
  )
}
rep_ids <- seq.int(rep_start, rep_end)
n <- bmk_env_int("GAMLSS_LONGITUDINAL_BENCHMARK_N", 120L)
seed_base <- bmk_env_int("GAMLSS_LONGITUDINAL_BENCHMARK_SEED", 20260612L)
interval_level <- bmk_env_num("GAMLSS_LONGITUDINAL_BENCHMARK_INTERVAL_LEVEL", 0.95)
max_elapsed_sec <- bmk_env_num("GAMLSS_LONGITUDINAL_BENCHMARK_MAX_ELAPSED_SEC", 180)
max_outer_iter <- bmk_env_int("GAMLSS_LONGITUDINAL_BENCHMARK_MAX_OUTER_ITER", 100L)
max_inner_iter <- bmk_env_int("GAMLSS_LONGITUDINAL_BENCHMARK_MAX_INNER_ITER", 100L)
checkpoint_every <- bmk_env_int("GAMLSS_LONGITUDINAL_BENCHMARK_CHECKPOINT_EVERY", 10L)

progress_path <- function(stem, completed = FALSE) {
  suffix <- if (isTRUE(completed)) "_by_rep.csv" else "_checkpoint.csv"
  file.path(run_dir, paste0(stem, suffix))
}

read_progress <- function(stem) {
  candidates <- c(progress_path(stem, completed = TRUE), progress_path(stem))
  path <- candidates[file.exists(candidates) & file.info(candidates)$size > 4L][1L]
  if (length(path) == 0L || is.na(path)) return(data.frame())
  tryCatch(
    read.csv(path, stringsAsFactors = FALSE, check.names = FALSE),
    error = function(e) stop("Could not read resumable benchmark progress: ", path, call. = FALSE)
  )
}

active_run_file <- file.path(bmk_output_root, "active_run_dir.txt")
latest_run_file <- file.path(bmk_output_root, "latest_run_dir.txt")
requested_run_dir <- bmk_env("GAMLSS_LONGITUDINAL_BENCHMARK_RESUME_RUN_DIR", "")
if (!nzchar(requested_run_dir) && file.exists(active_run_file)) {
  requested_run_dir <- trimws(readLines(active_run_file, warn = FALSE)[1L])
}
if (!nzchar(requested_run_dir) && file.exists(latest_run_file)) {
  requested_run_dir <- trimws(readLines(latest_run_file, warn = FALSE)[1L])
}
resuming <- nzchar(requested_run_dir) && dir.exists(requested_run_dir)
if (resuming) {
  run_dir <- normalizePath(requested_run_dir, winslash = "/", mustWork = TRUE)
  message("Resuming benchmark run: ", run_dir)
} else {
  stamp <- bmk_timestamp()
  run_dir <- file.path(bmk_output_root, paste0("run_", stamp))
  dir.create(run_dir, recursive = TRUE, showWarnings = FALSE)
  writeLines(normalizePath(run_dir, winslash = "/", mustWork = TRUE), active_run_file, useBytes = TRUE)
}

bmk_write_csv(bmk_session_info(source_used), file.path(run_dir, "session_info.csv"))
bmk_write_csv(bmk_available_methods(), file.path(run_dir, "comparator_status.csv"))

scenario_table <- do.call(rbind, lapply(scenarios, function(x) {
  data.frame(
    scenario = x$scenario,
    generator = x$generator,
    correlation = x$correlation,
    correlation_level = x$correlation_level,
    rho = x$rho,
    tau_base = x$tau_base,
    n_time = x$n_time,
    theta_formula = paste(deparse(x$theta_formula), collapse = ""),
    gee_correlations = paste(x$gee_correlations, collapse = ";"),
    stringsAsFactors = FALSE
  )
}))
scenario_grid_path <- file.path(run_dir, "scenario_grid.csv")
if (resuming && file.exists(scenario_grid_path)) {
  previous_scenarios <- read.csv(scenario_grid_path, stringsAsFactors = FALSE, check.names = FALSE)
  if (!identical(sort(previous_scenarios$scenario), sort(scenario_table$scenario))) {
    stop("The requested scenarios do not match the resumable benchmark run.", call. = FALSE)
  }
}
bmk_write_csv(scenario_table, scenario_grid_path)

family_table <- do.call(rbind, lapply(families, function(x) {
  data.frame(
    family = x$family,
    gamlss_family = x$gamlss_family,
    label = x$label,
    intercept = x$intercept,
    slope = x$slope,
    sigma = x$sigma,
    stringsAsFactors = FALSE
  )
}))
bmk_family_grid_path <- file.path(run_dir, "family_grid.csv")
if (resuming && file.exists(bmk_family_grid_path)) {
  previous_families <- read.csv(bmk_family_grid_path, stringsAsFactors = FALSE, check.names = FALSE)
  if (!identical(sort(previous_families$family), sort(family_table$family))) {
    stop("The requested families do not match the resumable benchmark run.", call. = FALSE)
  }
}
bmk_write_csv(family_table, bmk_family_grid_path)

as_progress_list <- function(x) if (is.data.frame(x) && nrow(x) > 0L) list(x) else list()
result_rows <- as_progress_list(read_progress("benchmark_results"))
coef_rows <- as_progress_list(read_progress("coefficient_results"))
complexity_rows <- as_progress_list(read_progress("fit_complexity"))
dependence_rows <- as_progress_list(read_progress("dependence_recovery"))
primary_status_rows <- as_progress_list(read_progress("primary_status"))

case_key <- function(x) paste(x$scenario, x$family, x$rep, sep = "::")
previous_status <- bmk_bind_rows_fill(primary_status_rows)
previous_results <- bmk_bind_rows_fill(result_rows)
completed_keys <- character()
if (nrow(previous_status) > 0L) {
  status_success <- previous_status$success %in% c(TRUE, "TRUE", "True", "true", "1")
  failed_keys <- case_key(previous_status[!status_success, , drop = FALSE])
  result_keys <- if (nrow(previous_results) > 0L) unique(case_key(previous_results)) else character()
  successful_keys <- case_key(previous_status[status_success, , drop = FALSE])
  completed_keys <- unique(c(failed_keys, intersect(successful_keys, result_keys)))
  message("Loaded ", length(completed_keys), " completed case(s) from resumable progress.")
}
case_index <- 0L
total_cases <- length(scenarios) * length(families) * length(rep_ids)

for (scenario_name in names(scenarios)) {
  scenario <- scenarios[[scenario_name]]
  for (family_name in names(families)) {
    spec <- families[[family_name]]
    for (rep_id in rep_ids) {
      case_index <- case_index + 1L
      seed <- seed_base + case_index
      current_key <- paste(scenario$scenario, spec$family, rep_id, sep = "::")
      if (current_key %in% completed_keys) {
        message(
          "[", case_index, "/", total_cases, "] skipping completed ",
          scenario$scenario, " / ", spec$family, " / rep ", rep_id
        )
        next
      }
      message(
        "[", case_index, "/", total_cases, "] ",
        scenario$scenario, " / ", spec$family, " / rep ", rep_id
      )

      sim <- bmk_capture_fit(bmk_simulate_dataset(spec, scenario, n = n, seed = seed))
      if (inherits(sim$value, "error")) {
        primary_status_rows[[length(primary_status_rows) + 1L]] <- data.frame(
          generator = scenario$generator,
          scenario = scenario$scenario,
          correlation = scenario$correlation,
          correlation_level = scenario$correlation_level,
          n_time = scenario$n_time,
          family = spec$family,
          gamlss_family = spec$gamlss_family,
          rep = rep_id,
          seed = seed,
          method = "rs_joint",
          success = FALSE,
          elapsed_sec = NA_real_,
          error = paste("Simulation failed:", conditionMessage(sim$value)),
          warning = if (length(sim$warnings)) paste(sim$warnings, collapse = " | ") else NA_character_,
          stringsAsFactors = FALSE
        )
        next
      }

      dat <- sim$value
      primary <- bmk_run_primary_fit(
        dat = dat,
        spec = spec,
        scenario = scenario,
        rep_id = rep_id,
        seed = seed,
        max_elapsed_sec = max_elapsed_sec,
        max_outer_iter = max_outer_iter,
        max_inner_iter = max_inner_iter
      )
      primary_status_rows[[length(primary_status_rows) + 1L]] <- primary$status

      bench <- bmk_capture_fit(
        bmk_run_benchmark_standard_models(
          dat = dat,
          primary_fit = primary$fit,
          spec = spec,
          scenario = scenario,
          rep_id = rep_id,
          seed = seed,
          interval_level = interval_level,
          primary_elapsed_sec = primary$status$elapsed_sec
        )
      )
      if (inherits(bench$value, "error")) {
        result_rows[[length(result_rows) + 1L]] <- data.frame(
          generator = scenario$generator,
          scenario = scenario$scenario,
          correlation = scenario$correlation,
          correlation_level = scenario$correlation_level,
          family = spec$family,
          gamlss_family = spec$gamlss_family,
          rep = rep_id,
          seed = seed,
          method = "benchmark_standard_models",
          success = FALSE,
          error = conditionMessage(bench$value),
          warning = if (length(bench$warnings)) paste(bench$warnings, collapse = " | ") else NA_character_,
          stringsAsFactors = FALSE
        )
        next
      }

      result_rows[[length(result_rows) + 1L]] <- bench$value$results
      coef_rows[[length(coef_rows) + 1L]] <- bench$value$coefficients
      complexity_rows[[length(complexity_rows) + 1L]] <- bench$value$complexity
      dependence_rows[[length(dependence_rows) + 1L]] <- bench$value$dependence

      if (checkpoint_every > 0L && case_index %% checkpoint_every == 0L) {
        if (length(result_rows) > 0L) {
          bmk_write_csv(bmk_bind_rows_fill(result_rows), file.path(run_dir, "benchmark_results_checkpoint.csv"))
        }
        if (length(coef_rows) > 0L) {
          bmk_write_csv(bmk_bind_rows_fill(coef_rows), file.path(run_dir, "coefficient_results_checkpoint.csv"))
        }
        if (length(complexity_rows) > 0L) {
          bmk_write_csv(bmk_bind_rows_fill(complexity_rows), file.path(run_dir, "fit_complexity_checkpoint.csv"))
        }
        if (length(dependence_rows) > 0L) {
          bmk_write_csv(bmk_bind_rows_fill(dependence_rows), file.path(run_dir, "dependence_recovery_checkpoint.csv"))
        }
        bmk_write_csv(bmk_bind_rows_fill(primary_status_rows), file.path(run_dir, "primary_status_checkpoint.csv"))
      }
    }
  }
}

results <- bmk_bind_rows_fill(result_rows)
coefficients <- bmk_bind_rows_fill(coef_rows)
complexity <- bmk_bind_rows_fill(complexity_rows)
dependence <- bmk_bind_rows_fill(dependence_rows)
primary_status <- bmk_bind_rows_fill(primary_status_rows)

bmk_write_csv(results, file.path(run_dir, "benchmark_results_by_rep.csv"))
bmk_write_csv(coefficients, file.path(run_dir, "coefficient_results_by_rep.csv"))
bmk_write_csv(complexity, file.path(run_dir, "fit_complexity_by_rep.csv"))
bmk_write_csv(dependence, file.path(run_dir, "dependence_recovery_by_rep.csv"))
bmk_write_csv(primary_status, file.path(run_dir, "primary_status_by_rep.csv"))

writeLines(run_dir, latest_run_file, useBytes = TRUE)
writeLines(run_dir, active_run_file, useBytes = TRUE)

message("Run complete: ", run_dir)
