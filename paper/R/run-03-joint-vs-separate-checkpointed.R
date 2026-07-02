Sys.setenv(
  GAMLSS_LONGITUDINAL_JSS_PROFILE = Sys.getenv(
    "GAMLSS_LONGITUDINAL_JSS_PROFILE",
    unset = "expanded"
  ),
  GAMLSS_LONGITUDINAL_JSS_JVS_REPS = Sys.getenv(
    "GAMLSS_LONGITUDINAL_JSS_JVS_REPS",
    unset = "100"
  ),
  GAMLSS_LONGITUDINAL_JSS_JVS_MAX_OUTER = Sys.getenv(
    "GAMLSS_LONGITUDINAL_JSS_JVS_MAX_OUTER",
    unset = "16"
  ),
  GAMLSS_LONGITUDINAL_JSS_JVS_MAX_INNER = Sys.getenv(
    "GAMLSS_LONGITUDINAL_JSS_JVS_MAX_INNER",
    unset = "8"
  ),
  GAMLSS_LONGITUDINAL_JSS_JVS_MAX_ELAPSED = Sys.getenv(
    "GAMLSS_LONGITUDINAL_JSS_JVS_MAX_ELAPSED",
    unset = "120"
  )
)

library(gamlss.longitudinal)

source("paper/R/replication-helpers.R")
source("paper/R/03-joint-vs-separate-optimization.R")

jss_checkpoint_message <- function(...) {
  cat(format(Sys.time(), "%Y-%m-%d %H:%M:%S %z"), "-", ..., "\n")
  flush.console()
}

jss_checkpoint_path <- function(checkpoint_dir, rep_idx) {
  file.path(checkpoint_dir, sprintf("rep-%03d.csv", rep_idx))
}

jss_checkpoint_complete <- function(path, rep_idx, expected_rows) {
  if (!file.exists(path) || file.info(path)$size == 0L) {
    return(FALSE)
  }
  x <- tryCatch(utils::read.csv(path, stringsAsFactors = FALSE), error = function(e) NULL)
  if (is.null(x) || nrow(x) < expected_rows || !"joint_review_rep" %in% names(x)) {
    return(FALSE)
  }
  all(x$joint_review_rep == rep_idx)
}

jss_write_checkpoint <- function(result, path) {
  tmp <- paste0(path, ".tmp")
  if (file.exists(tmp)) {
    unlink(tmp)
  }
  utils::write.csv(result, tmp, row.names = FALSE)
  if (file.exists(path)) {
    unlink(path)
  }
  if (!file.rename(tmp, path)) {
    stop("Could not move checkpoint into place: ", path)
  }
}

jss_run_checkpointed_03 <- function(settings) {
  paths <- jss_joint_output_paths(settings)
  candidates <- jss_joint_candidate_selection(settings)
  utils::write.csv(candidates, paths$candidate_selection, row.names = FALSE)

  cfg <- jss_joint_simulation_settings(settings, candidates)
  checkpoint_dir <- file.path(
    settings$out_dir,
    "checkpoints",
    paste0(
      jss_joint_module_id(),
      "-",
      cfg$reps,
      "rep-outer",
      cfg$max_outer_iter
    )
  )
  dir.create(checkpoint_dir, recursive = TRUE, showWarnings = FALSE)

  expected_rows <- length(cfg$families) * length(cfg$copulas) *
    length(cfg$designs) * 2L

  jss_checkpoint_message(
    "starting checkpointed module 03 run;",
    "reps=", cfg$reps,
    " expected_rows_per_rep=", expected_rows,
    " checkpoint_dir=", checkpoint_dir,
    sep = ""
  )

  for (rep_idx in seq_len(cfg$reps)) {
    checkpoint <- jss_checkpoint_path(checkpoint_dir, rep_idx)
    if (jss_checkpoint_complete(checkpoint, rep_idx, expected_rows)) {
      jss_checkpoint_message("skipping completed rep ", rep_idx, sep = "")
      next
    }

    jss_checkpoint_message("running rep ", rep_idx, " of ", cfg$reps, sep = "")
    result <- gamlss.longitudinal::run_coverage_simulations(
      families = cfg$families,
      copulas = cfg$copulas,
      methods = c("rs_separate", "rs_joint"),
      designs = cfg$designs,
      n = cfg$n,
      times = cfg$times,
      seed = settings$seed + 3000L + rep_idx * 100L,
      max_outer_iter = cfg$max_outer_iter,
      max_inner_iter = cfg$max_inner_iter,
      max_elapsed_sec = cfg$max_elapsed_sec,
      dependence = "strong",
      missingness = "none",
      output_dir = settings$data_dir,
      write_results = FALSE,
      write_summary = FALSE
    )
    result$joint_review_rep <- rep_idx
    result$profile <- settings$profile
    jss_write_checkpoint(result, checkpoint)
    jss_checkpoint_message(
      "wrote rep ",
      rep_idx,
      " rows=",
      nrow(result),
      " path=",
      checkpoint,
      sep = ""
    )
    rm(result)
    invisible(gc())
  }

  checkpoints <- jss_checkpoint_path(checkpoint_dir, seq_len(cfg$reps))
  complete <- vapply(
    seq_along(checkpoints),
    function(i) jss_checkpoint_complete(checkpoints[[i]], i, expected_rows),
    logical(1)
  )
  if (!all(complete)) {
    stop("Missing incomplete checkpoints: ", paste(which(!complete), collapse = ", "))
  }

  results <- jss_joint_bind_rows(lapply(checkpoints, utils::read.csv, stringsAsFactors = FALSE))
  results$profile <- settings$profile
  utils::write.csv(results, paths$results, row.names = FALSE)

  deltas <- jss_joint_delta_table(results)
  summary <- jss_joint_summary_table(results, deltas, candidates)
  metric_wins <- jss_joint_metric_wins(results)
  utils::write.csv(summary, paths$summary, row.names = FALSE)
  utils::write.csv(metric_wins, paths$metric_wins, row.names = FALSE)

  jss_joint_write_delta_figure(deltas, paths$deltas_figure)
  jss_joint_write_metric_dashboard(metric_wins, paths$metric_dashboard)

  status_path <- file.path(
    settings$logs_dir,
    paste0(jss_joint_module_id(), "-checkpointed-status.csv")
  )
  status <- data.frame(
    finished_at = format(Sys.time(), "%Y-%m-%d %H:%M:%S %z"),
    reps = cfg$reps,
    rows = nrow(results),
    checkpoint_dir = checkpoint_dir,
    results_path = paths$results,
    stringsAsFactors = FALSE
  )
  utils::write.csv(status, status_path, row.names = FALSE)
  jss_checkpoint_message("completed module 03 checkpointed run; rows=", nrow(results), sep = "")

  invisible(list(paths = paths, status = status_path, checkpoint_dir = checkpoint_dir))
}

settings <- jss_settings()
jss_run_checkpointed_03(settings)
