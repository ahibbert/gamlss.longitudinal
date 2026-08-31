source(file.path("paper", "R", "09-simulation-multivariate-longitudinal", "00-multivariate-setup.R"))

mvt_load_package()
run_dir <- mvt_read_run_dir()
grid_path <- file.path(run_dir, "scenario_grid.csv")
if (!file.exists(grid_path)) {
  stop("Missing scenario grid: ", grid_path, call. = FALSE)
}

methods <- mvt_env_vector("GAMLSS_LONGITUDINAL_MVT_RERUN_METHODS", character())
methods <- intersect(methods, mvt_allowed_comparators())
methods[methods == "gamCopula"] <- "gamCopula_markov"
methods <- unique(methods)
if (length(methods) == 0L) {
  stop("Set GAMLSS_LONGITUDINAL_MVT_RERUN_METHODS to one or more valid methods.", call. = FALSE)
}

grid <- utils::read.csv(grid_path, stringsAsFactors = FALSE, check.names = FALSE)
families <- mvt_env_vector("GAMLSS_LONGITUDINAL_MVT_RERUN_FAMILIES", character())
if (length(families) > 0L && "family_name" %in% names(grid)) {
  grid <- grid[grid$family_name %in% families, , drop = FALSE]
}
dependencies <- mvt_env_vector("GAMLSS_LONGITUDINAL_MVT_RERUN_DEPENDENCE", character())
if (length(dependencies) > 0L && "dependence_name" %in% names(grid)) {
  grid <- grid[grid$dependence_name %in% dependencies, , drop = FALSE]
}
case_ids <- mvt_env_vector("GAMLSS_LONGITUDINAL_MVT_RERUN_CASE_IDS", character())
if (length(case_ids) > 0L) {
  grid <- grid[grid$case_id %in% case_ids, , drop = FALSE]
}
progress_key <- paste(
  "selected_methods",
  paste(methods, collapse = "-"),
  paste(families, collapse = "-"),
  paste(dependencies, collapse = "-"),
  sep = "__"
)
progress_key <- gsub("[^A-Za-z0-9_-]+", "_", progress_key)
progress_path <- file.path(run_dir, paste0("rerun_progress_", progress_key, ".csv"))
completed_progress <- mvt_read_optional_csv(progress_path)
completed_case_ids <- if (nrow(completed_progress) > 0L && "case_id" %in% names(completed_progress)) {
  unique(completed_progress$case_id)
} else {
  character()
}
if (length(completed_case_ids) > 0L) {
  grid <- grid[!grid$case_id %in% completed_case_ids, , drop = FALSE]
}
max_cases <- mvt_env_int("GAMLSS_LONGITUDINAL_MVT_RERUN_MAX_CASES", 0L)
if (max_cases > 0L && nrow(grid) > max_cases) {
  grid <- grid[seq_len(max_cases), , drop = FALSE]
}
if (nrow(grid) == 0L) {
  stop("No scenario rows matched the selected rerun filters.", call. = FALSE)
}

metadata <- mvt_read_optional_csv(file.path(run_dir, "run_metadata.csv"))
seed_base <- mvt_env_int("GAMLSS_LONGITUDINAL_MVT_SEED", 20260818L)
if (nrow(metadata) > 0L && all(c("name", "value") %in% names(metadata))) {
  metadata_seed <- suppressWarnings(as.integer(metadata$value[metadata$name == "seed_base"][[1L]]))
  if (length(metadata_seed) > 0L && is.finite(metadata_seed)) seed_base <- metadata_seed
}

old_comparators <- Sys.getenv("GAMLSS_LONGITUDINAL_MVT_COMPARATORS", unset = NA_character_)
on.exit({
  if (is.na(old_comparators)) {
    Sys.unsetenv("GAMLSS_LONGITUDINAL_MVT_COMPARATORS")
  } else {
    Sys.setenv(GAMLSS_LONGITUDINAL_MVT_COMPARATORS = old_comparators)
  }
}, add = TRUE)
Sys.setenv(GAMLSS_LONGITUDINAL_MVT_COMPARATORS = paste(methods, collapse = ","))

message("Rerunning selected methods in: ", run_dir)
message("Methods: ", paste(methods, collapse = ", "))
message("Cases: ", nrow(grid))
if (length(completed_case_ids) > 0L) {
  message("Already completed for this selection: ", length(completed_case_ids))
}

rows <- stats::setNames(vector("list", length(mvt_result_names())), mvt_result_names())
for (name in names(rows)) rows[[name]] <- list()
current <- stats::setNames(vector("list", length(mvt_result_names())), mvt_result_names())
for (name in names(current)) {
  current[[name]] <- mvt_read_optional_csv(file.path(run_dir, paste0(name, "_by_rep.csv")))
}

write_current_outputs <- function() {
  for (name in names(current)) {
    if (is.data.frame(current[[name]]) && nrow(current[[name]]) > 0L) {
      mvt_write_csv(current[[name]], file.path(run_dir, paste0(name, "_by_rep.csv")))
    }
  }
}

append_progress <- function(done_case_ids) {
  if (length(done_case_ids) == 0L) return(invisible(NULL))
  existing <- mvt_read_optional_csv(progress_path)
  progress <- unique(rbind(
    existing[, intersect(c("case_id", "completed_at"), names(existing)), drop = FALSE],
    data.frame(
      case_id = done_case_ids,
      completed_at = format(Sys.time(), "%Y-%m-%d %H:%M:%S %Z"),
      stringsAsFactors = FALSE
    )
  ))
  mvt_write_csv(progress, progress_path)
}

checkpoint_every <- mvt_env_int("GAMLSS_LONGITUDINAL_MVT_RERUN_CHECKPOINT_EVERY", 5L)
checkpoint_case_ids <- character()
for (i in seq_len(nrow(grid))) {
  message("[", i, "/", nrow(grid), "] ", grid$case_id[[i]])
  out <- mvt_run_case(grid[i, , drop = FALSE], seed_base = seed_base, require_gamcopula = any(startsWith(methods, "gamCopula")))
  for (name in names(rows)) {
    dat <- out[[name]]
    if (is.data.frame(dat) && nrow(dat) > 0L && "method" %in% names(dat)) {
      new_dat <- dat[dat$method %in% methods, , drop = FALSE]
      rows[[name]][[length(rows[[name]]) + 1L]] <- new_dat
      existing <- current[[name]]
      if (is.data.frame(existing) && nrow(existing) > 0L && all(c("case_id", "method") %in% names(existing))) {
        keep <- !(existing$case_id == grid$case_id[[i]] & existing$method %in% methods)
        current[[name]] <- mvt_bind_rows_fill(existing[keep, , drop = FALSE], new_dat)
      } else {
        current[[name]] <- mvt_bind_rows_fill(existing, new_dat)
      }
    }
  }
  checkpoint_case_ids <- c(checkpoint_case_ids, grid$case_id[[i]])
  if (checkpoint_every > 0L && i %% checkpoint_every == 0L) {
    write_current_outputs()
    append_progress(checkpoint_case_ids)
    checkpoint_case_ids <- character()
    mvt_summarise_results(run_dir)
  }
}

write_current_outputs()
append_progress(checkpoint_case_ids)
mvt_summarise_results(run_dir)
message("Selected-method rerun complete: ", run_dir)
