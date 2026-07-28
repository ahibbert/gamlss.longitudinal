source(file.path("paper", "R", "08-simulation-sensitivity-correlation-misspecification", "standard-model-benchmarking", "00-benchmark-setup.R"))

source_used <- bmk_load_package()
bmk_require_namespaces(c("gamlss.dist", "mvtnorm", "geepack", "callr"), strict = TRUE)

source_dir <- bmk_env(
  "GAMLSS_LONGITUDINAL_BENCHMARK_SOURCE_DIR",
  trimws(readLines(file.path(bmk_output_root, "latest_combined_review_dir.txt"), warn = FALSE)[1L])
)
source_dir <- normalizePath(source_dir, winslash = "/", mustWork = TRUE)

stamp <- bmk_timestamp()
output_dir <- bmk_env(
  "GAMLSS_LONGITUDINAL_BENCHMARK_OUTPUT_DIR",
  file.path(bmk_output_root, paste0("combined_geepack_unstructured_", stamp))
)
dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)

message("Source: ", source_dir)
message("Output: ", output_dir)

read_required <- function(name) {
  path <- file.path(source_dir, name)
  if (!file.exists(path)) stop("Missing ", path, call. = FALSE)
  utils::read.csv(path, stringsAsFactors = FALSE, check.names = FALSE)
}

bind_rows_fill <- function(...) {
  pieces <- list(...)
  pieces <- pieces[vapply(pieces, nrow, integer(1L)) > 0L]
  if (length(pieces) == 0L) return(data.frame())
  cols <- unique(unlist(lapply(pieces, names), use.names = FALSE))
  pieces <- lapply(pieces, function(x) {
    missing <- setdiff(cols, names(x))
    for (col in missing) x[[col]] <- NA
    x[cols]
  })
  do.call(rbind, pieces)
}

gee_unstructured_failure_row <- function(scenario, spec, rep_id, seed, error, warning = NA_character_) {
  data.frame(
    generator = scenario$generator,
    scenario = scenario$scenario,
    correlation = scenario$correlation,
    correlation_level = scenario$correlation_level,
    family = spec$family,
    gamlss_family = spec$gamlss_family,
    rep = rep_id,
    seed = seed,
    method = "gee_unstructured",
    comparator = "gee_unstructured",
    comparator_class = "gee",
    estimator = "geepack::geeglm",
    package = "geepack",
    available = TRUE,
    success = FALSE,
    elapsed_sec = NA_real_,
    nobs = NA_integer_,
    logLik = NA_real_,
    AIC = NA_real_,
    mae = NA_real_,
    rmse = NA_real_,
    benchmark_mae = NA_real_,
    benchmark_rmse = NA_real_,
    benchmark_mean_bias = NA_real_,
    benchmark_mean_mae = NA_real_,
    benchmark_mean_rmse = NA_real_,
    benchmark_q90_mae = NA_real_,
    benchmark_neg_log_score = NA_real_,
    benchmark_upper_tail_error_90 = NA_real_,
    benchmark_theta_time_abs_error = NA_real_,
    benchmark_interval_coverage_95 = NA_real_,
    benchmark_interval_width_95 = NA_real_,
    benchmark_pit_ks_p_value = NA_real_,
    benchmark_pit_mean_abs_error = NA_real_,
    benchmark_tail_error_lower_05 = NA_real_,
    benchmark_tail_error_upper_05 = NA_real_,
    warning = warning,
    error = error,
    stringsAsFactors = FALSE
  )
}

run_gee_unstructured_case <- function(repo_root, scenario_name, family_name, rep_id, seed, n, interval_level, source_mode) {
  setwd(repo_root)
  Sys.setenv(GAMLSS_LONGITUDINAL_BENCHMARK_SOURCE = source_mode)
  source(file.path("paper", "R", "08-simulation-sensitivity-correlation-misspecification", "standard-model-benchmarking", "00-benchmark-setup.R"))
  source_used <- bmk_load_package()
  bmk_require_namespaces(c("gamlss.dist", "mvtnorm", "geepack"), strict = TRUE)
  scenario <- bmk_scenario_specs()[[scenario_name]]
  spec <- bmk_family_specs(include_binary = TRUE)[[family_name]]

  sim <- bmk_capture_fit(bmk_simulate_dataset(spec, scenario, n = n, seed = seed))
  if (inherits(sim$value, "error")) {
    return(list(
      results = gee_unstructured_failure_row(
        scenario = scenario,
        spec = spec,
        rep_id = rep_id,
        seed = seed,
        error = paste("Simulation failed:", conditionMessage(sim$value)),
        warning = if (length(sim$warnings)) paste(sim$warnings, collapse = " | ") else NA_character_
      ),
      coefficients = data.frame(),
      source_used = source_used
    ))
  }

  dat <- sim$value
  bench <- bmk_capture_fit(
    gamlss.longitudinal::benchmark_standard_models(
      data = dat,
      formula = response ~ x + z_null,
      subject_var = "subject",
      family = spec$standard_family,
      comparators = "gee",
      correlation = "unstructured",
      truth_family = spec$gamlss_family,
      interval_level = interval_level,
      waves = dat$time
    )
  )

  if (inherits(bench$value, "error")) {
    return(list(
      results = gee_unstructured_failure_row(
        scenario = scenario,
        spec = spec,
        rep_id = rep_id,
        seed = seed,
        error = conditionMessage(bench$value),
        warning = if (length(bench$warnings)) paste(bench$warnings, collapse = " | ") else NA_character_
      ),
      coefficients = data.frame(),
      source_used = source_used
    ))
  }

  list(
    results = bmk_prefix_results(
      bench$value$results,
      scenario = scenario,
      spec = spec,
      rep_id = rep_id,
      seed = seed,
      suffix = "_unstructured"
    ),
    coefficients = bmk_annotate_coefficients(
      bmk_prefix_coefficients(
        bench$value$coefficients$long,
        scenario = scenario,
        spec = spec,
        rep_id = rep_id,
        seed = seed,
        suffix = "_unstructured"
      ),
      spec = spec,
      level = interval_level
    ),
    source_used = source_used
  )
}

results <- read_required("benchmark_results_by_rep.csv")
coef <- read_required("coefficient_results_by_rep.csv")
status <- read_required("primary_status_by_rep.csv")

keep_methods <- c("rs_joint", "glm", "gee_exchangeable", "gee_ar1", "gee_unstructured")
filtered_results <- results[results$method %in% keep_methods, , drop = FALSE]
filtered_coef <- coef[coef$method %in% keep_methods, , drop = FALSE]

checkpoint_results_path <- file.path(output_dir, "benchmark_results_checkpoint.csv")
checkpoint_coef_path <- file.path(output_dir, "coefficient_results_checkpoint.csv")
if (file.exists(checkpoint_results_path)) {
  message("Resuming from checkpoint: ", checkpoint_results_path)
  filtered_results <- utils::read.csv(checkpoint_results_path, stringsAsFactors = FALSE, check.names = FALSE)
  if (file.exists(checkpoint_coef_path)) {
    filtered_coef <- utils::read.csv(checkpoint_coef_path, stringsAsFactors = FALSE, check.names = FALSE)
  }
}

case_cols <- c("scenario", "family", "rep", "seed")
cases <- unique(results[results$method == "rs_joint", case_cols, drop = FALSE])
existing_unstructured <- unique(filtered_results[filtered_results$method == "gee_unstructured", case_cols, drop = FALSE])
if (nrow(existing_unstructured) > 0L) {
  key <- do.call(paste, c(cases[case_cols], sep = "\r"))
  existing_key <- do.call(paste, c(existing_unstructured[case_cols], sep = "\r"))
  cases <- cases[!key %in% existing_key, , drop = FALSE]
}

max_cases <- bmk_env_int("GAMLSS_LONGITUDINAL_BENCHMARK_MAX_CASES", nrow(cases))
cases <- cases[seq_len(min(nrow(cases), max_cases)), , drop = FALSE]

scenarios <- bmk_scenario_specs()
families <- bmk_family_specs(include_binary = TRUE)
n <- bmk_env_int("GAMLSS_LONGITUDINAL_BENCHMARK_N", 120L)
interval_level <- bmk_env_num("GAMLSS_LONGITUDINAL_BENCHMARK_INTERVAL_LEVEL", 0.95)
checkpoint_every <- bmk_env_int("GAMLSS_LONGITUDINAL_BENCHMARK_CHECKPOINT_EVERY", 10L)
case_timeout_sec <- bmk_env_num("GAMLSS_LONGITUDINAL_BENCHMARK_GEE_UNSTRUCTURED_TIMEOUT_SEC", 90)
source_mode <- bmk_env("GAMLSS_LONGITUDINAL_BENCHMARK_SOURCE", "installed")

new_results <- list()
new_coef <- list()

write_checkpoint <- function() {
  bmk_write_csv(bind_rows_fill(filtered_results, do.call(bind_rows_fill, new_results)), checkpoint_results_path)
  bmk_write_csv(bind_rows_fill(filtered_coef, do.call(bind_rows_fill, new_coef)), checkpoint_coef_path)
}

for (i in seq_len(nrow(cases))) {
  scenario <- scenarios[[cases$scenario[i]]]
  spec <- families[[cases$family[i]]]
  rep_id <- cases$rep[i]
  seed <- cases$seed[i]

  message("[", i, "/", nrow(cases), "] gee_unstructured / ", scenario$scenario, " / ", spec$family, " / rep ", rep_id)

  case_result <- tryCatch(
    callr::r(
      run_gee_unstructured_case,
      args = list(
        repo_root = bmk_repo_root,
        scenario_name = cases$scenario[i],
        family_name = cases$family[i],
        rep_id = rep_id,
        seed = seed,
        n = n,
        interval_level = interval_level,
        source_mode = source_mode
      ),
      timeout = case_timeout_sec,
      spinner = FALSE,
      show = FALSE
    ),
    error = function(e) e
  )

  if (inherits(case_result, "error")) {
    new_results[[length(new_results) + 1L]] <- gee_unstructured_failure_row(
      scenario = scenario,
      spec = spec,
      rep_id = rep_id,
      seed = seed,
      error = paste("Unstructured GEE failed or timed out:", conditionMessage(case_result))
    )
  } else {
    new_results[[length(new_results) + 1L]] <- case_result$results
    new_coef[[length(new_coef) + 1L]] <- case_result$coefficients
  }

  if (checkpoint_every > 0L && i %% checkpoint_every == 0L) {
    write_checkpoint()
  }
}

augmented_results <- bind_rows_fill(filtered_results, do.call(bind_rows_fill, new_results))
augmented_coef <- bind_rows_fill(filtered_coef, do.call(bind_rows_fill, new_coef))

bmk_write_csv(augmented_results, file.path(output_dir, "benchmark_results_by_rep.csv"))
bmk_write_csv(augmented_coef, file.path(output_dir, "coefficient_results_by_rep.csv"))
bmk_write_csv(augmented_results, checkpoint_results_path)
bmk_write_csv(augmented_coef, checkpoint_coef_path)
bmk_write_csv(status, file.path(output_dir, "primary_status_by_rep.csv"))

for (name in c("scenario_grid.csv", "family_grid.csv", "comparator_status.csv", "combined_sources.csv")) {
  from <- file.path(source_dir, name)
  if (file.exists(from)) file.copy(from, file.path(output_dir, name), overwrite = TRUE)
}

bmk_write_csv(
  data.frame(
    item = c("timestamp", "source_dir", "source_used", "methods", "case_timeout_sec"),
    value = c(stamp, source_dir, source_used, paste(keep_methods, collapse = ","), case_timeout_sec),
    stringsAsFactors = FALSE
  ),
  file.path(output_dir, "augmentation_info.csv")
)

writeLines(output_dir, file.path(bmk_output_root, "latest_combined_review_dir.txt"), useBytes = TRUE)

message("Augmented review written: ", output_dir)
