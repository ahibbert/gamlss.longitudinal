source(file.path("paper", "R", "08-simulation-sensitivity-correlation-misspecification", "standard-model-benchmarking", "00-benchmark-setup.R"))

source_used <- bmk_load_package()
bmk_require_namespaces(c("gamlss.dist", "mvtnorm", "geepack"), strict = TRUE)

run_dir <- bmk_env("GAMLSS_LONGITUDINAL_BENCHMARK_REPORT_DIR", "")
if (!nzchar(run_dir)) {
  run_dir <- file.path(
    bmk_output_root,
    "rs-joint-standard-model-grid",
    "run_20260619_t20_t50_combined"
  )
}
run_dir <- normalizePath(run_dir, winslash = "/", mustWork = TRUE)

read_required <- function(name) {
  path <- file.path(run_dir, name)
  if (!file.exists(path)) stop("Missing ", path, call. = FALSE)
  utils::read.csv(path, stringsAsFactors = FALSE, check.names = FALSE)
}

status <- read_required("primary_status_by_rep.csv")
status <- status[status$method == "rs_joint", , drop = FALSE]

family_filter <- bmk_env_vector("GAMLSS_LONGITUDINAL_BENCHMARK_FAMILIES", character())
if (length(family_filter) > 0L) {
  status <- status[status$family %in% family_filter, , drop = FALSE]
}

scenario_filter <- bmk_env_vector("GAMLSS_LONGITUDINAL_BENCHMARK_SCENARIOS", character())
if (length(scenario_filter) > 0L) {
  status <- status[status$scenario %in% scenario_filter, , drop = FALSE]
}

max_cases <- bmk_env_int("GAMLSS_LONGITUDINAL_BENCHMARK_MAX_CASES", nrow(status))
status <- status[seq_len(min(nrow(status), max_cases)), , drop = FALSE]
if (nrow(status) == 0L) stop("No cases selected for all-pair dependence recomputation.", call. = FALSE)

case_start <- bmk_env_int("GAMLSS_LONGITUDINAL_BENCHMARK_CASE_START", 1L)
case_end <- bmk_env_int("GAMLSS_LONGITUDINAL_BENCHMARK_CASE_END", nrow(status))
case_start <- max(1L, case_start)
case_end <- min(nrow(status), case_end)
if (case_end < case_start) stop("Selected case range is empty.", call. = FALSE)
status <- status[seq.int(case_start, case_end), , drop = FALSE]

base_scenario <- function(x) sub("_t[0-9]+$", "", x)
timepoints <- sort(unique(as.integer(status$n_time)))
scenarios <- bmk_expand_scenarios_for_timepoints(bmk_scenario_specs(), timepoints)
families <- bmk_family_specs(include_binary = TRUE)

n <- bmk_env_int("GAMLSS_LONGITUDINAL_BENCHMARK_N", 120L)
interval_level <- bmk_env_num("GAMLSS_LONGITUDINAL_BENCHMARK_INTERVAL_LEVEL", 0.95)
max_elapsed_sec <- bmk_env_num("GAMLSS_LONGITUDINAL_BENCHMARK_MAX_ELAPSED_SEC", 180)
max_outer_iter <- bmk_env_int("GAMLSS_LONGITUDINAL_BENCHMARK_MAX_OUTER_ITER", 100L)
max_inner_iter <- bmk_env_int("GAMLSS_LONGITUDINAL_BENCHMARK_MAX_INNER_ITER", 100L)
checkpoint_every <- bmk_env_int("GAMLSS_LONGITUDINAL_BENCHMARK_CHECKPOINT_EVERY", 10L)

output_path <- bmk_env(
  "GAMLSS_LONGITUDINAL_BENCHMARK_DEPENDENCE_OUTPUT",
  file.path(run_dir, "dependence_recovery_by_rep.csv")
)
checkpoint_path <- bmk_env(
  "GAMLSS_LONGITUDINAL_BENCHMARK_DEPENDENCE_CHECKPOINT",
  sub("\\.csv$", "_checkpoint.csv", output_path)
)
dependence_rows <- list()
if (file.exists(checkpoint_path)) {
  message("Resuming from checkpoint: ", checkpoint_path)
  checkpoint <- utils::read.csv(checkpoint_path, stringsAsFactors = FALSE, check.names = FALSE)
  dependence_rows <- list(checkpoint)
  done <- unique(checkpoint[c("scenario", "family", "rep", "seed")])
  done_key <- do.call(paste, c(done, sep = "\r"))
  case_key <- do.call(paste, c(status[c("scenario", "family", "rep", "seed")], sep = "\r"))
  status <- status[!case_key %in% done_key, , drop = FALSE]
}

write_checkpoint <- function() {
  if (length(dependence_rows) > 0L) {
    bmk_write_csv(bmk_bind_rows_fill(dependence_rows), checkpoint_path)
  }
}

total_cases <- nrow(status)
for (i in seq_len(total_cases)) {
  scenario_name <- status$scenario[i]
  family_name <- status$family[i]
  scenario <- scenarios[[scenario_name]]
  if (is.null(scenario)) {
    scenario <- bmk_resize_scenario_time(bmk_scenario_specs()[[base_scenario(scenario_name)]], status$n_time[i])
    scenario$scenario <- scenario_name
  }
  spec <- families[[family_name]]
  rep_id <- status$rep[i]
  seed <- status$seed[i]

  message("[", i, "/", total_cases, "] all-pair dependence / ", scenario$scenario, " / ", spec$family, " / rep ", rep_id)

  sim <- bmk_capture_fit(bmk_simulate_dataset(spec, scenario, n = n, seed = seed))
  if (inherits(sim$value, "error")) next
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
  if (!is.null(primary$fit)) {
    dependence_rows[[length(dependence_rows) + 1L]] <- bmk_dependence_recovery_row(
      dat, primary$fit, "rs_joint", scenario, spec, rep_id, seed
    )
  }

  for (corstr in scenario$gee_correlations) {
    timeout <- if (identical(corstr, "unstructured")) {
      bmk_env_num("GAMLSS_LONGITUDINAL_BENCHMARK_GEE_UNSTRUCTURED_TIMEOUT_SEC", Inf)
    } else {
      bmk_env_num("GAMLSS_LONGITUDINAL_BENCHMARK_GEE_TIMEOUT_SEC", Inf)
    }
    gee_capture <- bmk_capture_fit(
      bmk_run_gee_benchmark(dat, spec, corstr, interval_level, timeout = timeout)
    )
    if (inherits(gee_capture$value, "error")) next
    gee_fit <- gee_capture$value$fits$gee
    dependence_rows[[length(dependence_rows) + 1L]] <- bmk_dependence_recovery_row(
      dat, gee_fit, paste0("gee_", corstr), scenario, spec, rep_id, seed, corstr = corstr
    )
  }

  if (checkpoint_every > 0L && i %% checkpoint_every == 0L) write_checkpoint()
}

all_pair_dependence <- bmk_bind_rows_fill(dependence_rows)
old_path <- file.path(run_dir, "dependence_recovery_by_rep.csv")
if (normalizePath(output_path, winslash = "/", mustWork = FALSE) == normalizePath(old_path, winslash = "/", mustWork = FALSE) &&
    file.exists(old_path)) {
  backup_path <- file.path(run_dir, "dependence_recovery_adjacent_by_rep_backup.csv")
  if (!file.exists(backup_path)) file.copy(old_path, backup_path, overwrite = FALSE)
}
bmk_write_csv(all_pair_dependence, output_path)
write_checkpoint()

message("All-pair dependence rows written: ", output_path)
message("Source used: ", source_used)
