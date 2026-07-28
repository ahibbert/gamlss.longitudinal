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

stamp <- bmk_timestamp()
run_dir <- file.path(bmk_output_root, paste0("run_", stamp))
dir.create(run_dir, recursive = TRUE, showWarnings = FALSE)

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
bmk_write_csv(scenario_table, file.path(run_dir, "scenario_grid.csv"))

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
bmk_write_csv(family_table, file.path(run_dir, "family_grid.csv"))

result_rows <- list()
coef_rows <- list()
complexity_rows <- list()
dependence_rows <- list()
primary_status_rows <- list()
case_index <- 0L
total_cases <- length(scenarios) * length(families) * length(rep_ids)

for (scenario_name in names(scenarios)) {
  scenario <- scenarios[[scenario_name]]
  for (family_name in names(families)) {
    spec <- families[[family_name]]
    for (rep_id in rep_ids) {
      case_index <- case_index + 1L
      seed <- seed_base + case_index
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

latest_file <- file.path(bmk_output_root, "latest_run_dir.txt")
writeLines(run_dir, latest_file, useBytes = TRUE)

message("Run complete: ", run_dir)
