#!/usr/bin/env Rscript

# Checkpointed full ACTG backward selection for the four zero-inflated paper
# routes. Each route is staged across all distributional parameters and then
# cycled until a complete pass leaves every selected formula unchanged.

script_argument <- grep("^--file=", commandArgs(trailingOnly = FALSE), value = TRUE)
script_path <- if (length(script_argument)) {
  normalizePath(sub("^--file=", "", script_argument[[1L]]), mustWork = TRUE)
} else {
  normalizePath(file.path(
    getwd(), "paper", "reviews", "2026-09-05-medical-candidates",
    "run-zero-inflated-four-route-build.R"
  ), mustWork = TRUE)
}
review_dir <- dirname(script_path)
package_root <- normalizePath(file.path(review_dir, "..", "..", ".."), mustWork = TRUE)
screen_script <- file.path(review_dir, "run-screening-v2.R")

argument_value <- function(name, default) {
  prefix <- paste0("--", name, "=")
  hit <- grep(paste0("^", prefix), commandArgs(trailingOnly = TRUE), value = TRUE)
  if (!length(hit)) return(default)
  sub(prefix, "", hit[[1L]], fixed = TRUE)
}

workers <- suppressWarnings(as.integer(argument_value("workers", "4")))
max_cycles <- suppressWarnings(as.integer(argument_value("max-cycles", "10")))
initial_inner_iter <- suppressWarnings(as.integer(argument_value("initial-inner-iter", "2")))
max_elapsed_sec <- suppressWarnings(as.numeric(argument_value("max-elapsed-sec", "1190")))
results_prefix <- argument_value("results-prefix", "actg-zero-inflated-defaults-2026-10-06")
dry_run <- tolower(argument_value("dry-run", "false")) %in% c("true", "1", "yes")
requested_routes <- trimws(strsplit(
  tolower(argument_value("routes", "zisichel-g,zisichel-t,zibnb-g,zibnb-t")),
  ",", fixed = TRUE
)[[1L]])

if (!is.finite(workers) || workers < 1L) stop("--workers must be positive.")
if (!is.finite(max_cycles) || max_cycles < 1L) stop("--max-cycles must be positive.")
if (!is.finite(initial_inner_iter) || initial_inner_iter < 1L) {
  stop("--initial-inner-iter must be positive.")
}
if (!is.finite(max_elapsed_sec) || max_elapsed_sec <= 0) {
  stop("--max-elapsed-sec must be positive.")
}
if (!grepl("^[A-Za-z0-9._-]+$", results_prefix)) {
  stop("--results-prefix may contain only letters, numbers, periods, underscores, and hyphens.")
}

route_registry <- list(
  `zisichel-g` = list(
    margin = "ZISICHEL", copula = "G",
    seed = file.path(
      package_root, "paper", "reviews", "2026-10-05-actg-zero-inflated-performance",
      "g-maximal-integrated-bt2.rds"
    )
  ),
  `zisichel-t` = list(
    margin = "ZISICHEL", copula = "t",
    seed = file.path(
      package_root, "paper", "reviews", "2026-10-05-actg-zero-inflated-performance",
      "zisichel-t-maximal-from-full-candidate-state-bfgs30-lambda5-bt2.rds"
    )
  ),
  `zibnb-g` = list(
    margin = "ZIBNB", copula = "G",
    seed = file.path(
      package_root, "paper", "reviews", "2026-10-05-actg-zero-inflated-performance",
      "zibnb-g-maximal-integrated-bfgs170-after2-lambda5-bt2.rds"
    )
  ),
  `zibnb-t` = list(
    margin = "ZIBNB", copula = "t",
    seed = file.path(
      package_root, "paper", "reviews", "2026-10-05-actg-zero-inflated-performance",
      "zibnb-t-maximal-integrated-rescue50x2-refine2-lambda5-bt2.rds"
    )
  )
)
if (!length(requested_routes) || any(!requested_routes %in% names(route_registry))) {
  stop("--routes must contain one or more of: ", paste(names(route_registry), collapse = ","), ".")
}
requested_routes <- unique(requested_routes)

package_source_files <- c(
  file.path(package_root, "DESCRIPTION"),
  list.files(file.path(package_root, "R"), pattern = "\\.[Rr]$", full.names = TRUE)
)
package_source_fingerprint <- tools::md5sum(package_source_files)
rscript <- file.path(R.home("bin"), if (.Platform$OS.type == "windows") "Rscript.exe" else "Rscript")

atomic_save_rds <- function(value, path) {
  dir.create(dirname(path), recursive = TRUE, showWarnings = FALSE)
  temporary <- tempfile("four-route-", tmpdir = dirname(path), fileext = ".rds")
  on.exit(unlink(temporary), add = TRUE)
  saveRDS(value, temporary, compress = "gzip")
  if (!file.copy(temporary, path, overwrite = TRUE)) stop("Could not save ", path, ".")
}

atomic_write_csv <- function(value, path) {
  temporary <- tempfile("four-route-", tmpdir = dirname(path), fileext = ".csv")
  on.exit(unlink(temporary), add = TRUE)
  utils::write.csv(value, temporary, row.names = FALSE, na = "")
  if (!file.copy(temporary, path, overwrite = TRUE)) stop("Could not save ", path, ".")
}

extract_fit <- function(path) {
  if (!file.exists(path)) return(NULL)
  value <- tryCatch(readRDS(path), error = function(error) NULL)
  if (inherits(value, "gamlss.longitudinal")) return(value)
  if (is.list(value) && inherits(value$fit, "gamlss.longitudinal")) return(value$fit)
  if (is.list(value) && inherits(value$result$fit, "gamlss.longitudinal")) return(value$result$fit)
  NULL
}

formula_signature <- function(path) {
  formulas <- readRDS(path)
  vapply(formulas, function(formula) {
    if (is.null(formula)) return("<NULL>")
    paste(deparse(formula, width.cutoff = 500L), collapse = " ")
  }, character(1L))
}

copula_label <- function(code) {
  if (identical(code, "G")) "Gumbel (G)" else "Student-t (t)"
}

valid_completed_stage <- function(manifest_path, fit_path, route, stage, parameters = NULL) {
  if (!file.exists(manifest_path) || !file.exists(fit_path)) return(FALSE)
  manifest <- tryCatch(readRDS(manifest_path), error = function(error) NULL)
  fit_info <- file.info(fit_path)
  accepted_path <- if (is.list(manifest)) manifest$accepted_path else NULL
  final_path_row <- if (is.data.frame(accepted_path) && nrow(accepted_path)) {
    accepted_path[nrow(accepted_path), , drop = FALSE]
  } else {
    NULL
  }
  manifest_fit_path <- if (is.list(manifest)) manifest$final_fit_path else NULL
  same_fit_path <- is.character(manifest_fit_path) && length(manifest_fit_path) == 1L &&
    identical(
      normalizePath(manifest_fit_path, winslash = "/", mustWork = FALSE),
      normalizePath(fit_path, winslash = "/", mustWork = FALSE)
    )
  is.list(manifest) && isTRUE(is.finite(fit_info$size)) && fit_info$size > 0 &&
    same_fit_path && is.data.frame(final_path_row) &&
    isTRUE(final_path_row$converged[[1L]]) && isTRUE(final_path_row$accepted[[1L]]) &&
    identical(manifest$stage, stage) &&
    identical(manifest$margin_family, route$margin) &&
    identical(manifest$copula_family, copula_label(route$copula)) &&
    identical(manifest$package_source_fingerprint, package_source_fingerprint) &&
    (is.null(parameters) || identical(manifest$parameters, parameters))
}

all_stage_parameters <- function(route) {
  c("mu", "sigma", "nu", "tau", "theta", if (identical(route$copula, "t")) "zeta")
}

initial_stages <- function(route) {
  c("mu-smooth", "sigma", "nu", "tau", "theta", if (identical(route$copula, "t")) "zeta")
}

run_route <- function(route_name) {
  route <- route_registry[[route_name]]
  results_subdir <- paste0(results_prefix, "-", route_name)
  results_dir <- file.path(review_dir, "results", results_subdir)
  dir.create(results_dir, recursive = TRUE, showWarnings = FALSE)
  stage_log_path <- file.path(results_dir, "four-route-stage-log.csv")
  completion_path <- file.path(results_dir, "SELECTION-COMPLETE.rds")

  if (file.exists(completion_path)) {
    completion <- tryCatch(readRDS(completion_path), error = function(error) NULL)
    if (is.list(completion) && isTRUE(completion$stable) &&
        identical(completion$package_source_fingerprint, package_source_fingerprint)) {
      cat("Reusing completed route ", route_name, ".\n", sep = "")
      return(completion)
    }
  }

  stage_rows <- if (file.exists(stage_log_path)) {
    utils::read.csv(stage_log_path, stringsAsFactors = FALSE)
  } else {
    data.frame(
      route = character(), label = character(), started_at = character(),
      completed_at = character(), elapsed_sec = numeric(), status = integer(),
      formula_path = character(), fit_path = character(), reused = logical(),
      stringsAsFactors = FALSE
    )
  }

  family_stem <- paste0(
    "actg-", tolower(route$margin),
    if (identical(route$copula, "t")) "-t" else ""
  )

  run_screen <- function(label, stage, stem, previous = NULL, input_formulas = NULL,
                         parameter = NULL) {
    formula_path <- file.path(results_dir, paste0(stem, "-final-formulas.rds"))
    fit_path <- file.path(results_dir, paste0(stem, "-final-fit.rds"))
    manifest_path <- file.path(results_dir, paste0(stem, "-run-manifest.rds"))
    manifest_parameter <- if (identical(stage, "refine")) parameter else NULL
    if (valid_completed_stage(
      manifest_path, fit_path, route, stage,
      parameters = manifest_parameter
    )) {
      cat("Reusing ", route_name, " / ", label, ".\n", sep = "")
      return(list(formula_path = formula_path, fit_path = fit_path, reused = TRUE))
    }

    start_fit <- if (!is.null(previous)) previous$fit_path else route$seed
    if (!is.null(start_fit) && !file.exists(start_fit)) start_fit <- NULL
    args <- c(
      screen_script,
      paste0("--margin=", route$margin),
      paste0("--copula=", route$copula),
      paste0("--stage=", stage),
      "--criterion=AIC", "--hierarchy=strong", "--resume=true",
      paste0("--workers=", workers),
      paste0("--initial-inner-iter=", initial_inner_iter),
      paste0("--max-elapsed-sec=", format(max_elapsed_sec, scientific = FALSE)),
      paste0("--results-subdir=", results_subdir),
      paste0("--output-stem=", stem)
    )
    if (!is.null(parameter)) args <- c(args, paste0("--parameters=", parameter))
    if (!is.null(input_formulas)) {
      args <- c(args, paste0(
        "--input-formulas=", normalizePath(input_formulas, winslash = "/", mustWork = TRUE)
      ))
    }
    if (!is.null(start_fit)) {
      args <- c(args, paste0(
        "--initial-fit=", normalizePath(start_fit, winslash = "/", mustWork = TRUE)
      ))
    }

    stdout_path <- file.path(results_dir, paste0(stem, ".stdout.log"))
    stderr_path <- file.path(results_dir, paste0(stem, ".stderr.log"))
    started_at <- Sys.time()
    cat("Starting ", route_name, " / ", label, " at ", format(started_at), ".\n", sep = "")
    status <- system2(
      rscript, vapply(args, shQuote, character(1L)),
      stdout = stdout_path, stderr = stderr_path
    )
    completed_at <- Sys.time()
    stage_rows <<- rbind(stage_rows, data.frame(
      route = route_name, label = label,
      started_at = format(started_at, tz = "Australia/Sydney"),
      completed_at = format(completed_at, tz = "Australia/Sydney"),
      elapsed_sec = as.numeric(difftime(completed_at, started_at, units = "secs")),
      status = status, formula_path = formula_path, fit_path = fit_path,
      reused = FALSE, stringsAsFactors = FALSE
    ))
    atomic_write_csv(stage_rows, stage_log_path)
    if (!identical(status, 0L) || !valid_completed_stage(
      manifest_path, fit_path, route, stage,
      parameters = manifest_parameter
    )) {
      stderr_tail <- if (file.exists(stderr_path)) {
        paste(tail(readLines(stderr_path, warn = FALSE), 20L), collapse = "\n")
      } else {
        "No stderr log was produced."
      }
      stop("Stage ", route_name, " / ", label, " failed.\n", stderr_tail)
    }
    list(formula_path = formula_path, fit_path = fit_path, reused = FALSE)
  }

  previous <- NULL
  for (stage in initial_stages(route)) {
    stem <- paste0(family_stem, "-backward-selection-", stage)
    previous <- run_screen(stage, stage, stem, previous = previous)
  }

  parameters <- all_stage_parameters(route)
  cycle_rows <- list()
  stable <- FALSE
  for (cycle in seq_len(max_cycles)) {
    start_signature <- formula_signature(previous$formula_path)
    for (parameter in parameters) {
      stem <- sprintf("%s-refinement-cycle-%02d-%s", family_stem, cycle, parameter)
      previous <- run_screen(
        sprintf("refinement cycle %02d: %s", cycle, parameter),
        "refine", stem, previous = previous,
        input_formulas = previous$formula_path, parameter = parameter
      )
    }
    end_signature <- formula_signature(previous$formula_path)
    stable <- identical(start_signature, end_signature)
    cycle_rows[[cycle]] <- data.frame(
      cycle = cycle, stable = stable,
      started_signature = paste(start_signature, collapse = " | "),
      completed_signature = paste(end_signature, collapse = " | "),
      stringsAsFactors = FALSE
    )
    atomic_write_csv(
      do.call(rbind, cycle_rows),
      file.path(results_dir, "refinement-cycle-summary.csv")
    )
    if (stable) break
  }
  if (!stable) {
    stop("Route ", route_name, " did not stabilize after ", max_cycles, " complete cycles.")
  }

  final_formulas <- readRDS(previous$formula_path)
  final_fit <- extract_fit(previous$fit_path)
  selected_formulas_path <- file.path(results_dir, "selected-formulas.rds")
  selected_fit_path <- file.path(results_dir, "selected-fit.rds")
  atomic_save_rds(final_formulas, selected_formulas_path)
  atomic_save_rds(final_fit, selected_fit_path)
  completion <- list(
    completed_at = Sys.time(), route = route_name,
    margin_family = route$margin, copula_family = route$copula,
    stable = TRUE, cycles_completed = length(cycle_rows),
    selected_formulas = final_formulas,
    selected_formulas_path = selected_formulas_path,
    selected_fit_path = selected_fit_path,
    stage_log_path = stage_log_path,
    package_source_fingerprint = package_source_fingerprint
  )
  atomic_save_rds(completion, completion_path)
  cat("Completed stable route ", route_name, ".\n", sep = "")
  completion
}

started_at <- Sys.time()
if (dry_run) {
  cat("ACTG four-route build dry run\n")
  for (route_name in requested_routes) {
    route <- route_registry[[route_name]]
    seed_fit <- extract_fit(route$seed)
    seed_ok <- inherits(seed_fit, "gamlss.longitudinal") &&
      identical(seed_fit$margin_dist$family[[1L]], route$margin) &&
      identical(seed_fit$copula_dist, route$copula) &&
      isTRUE(seed_fit$convergence$converged)
    cat(
      "  ", route_name, ": stages=",
      paste(initial_stages(route), collapse = ","),
      "; refinement=", paste(all_stage_parameters(route), collapse = ","),
      "; warm seed=", if (seed_ok) normalizePath(route$seed, winslash = "/") else "none",
      "\n", sep = ""
    )
  }
  quit(save = "no", status = 0L)
}
route_results <- lapply(requested_routes, run_route)
names(route_results) <- requested_routes
summary_path <- file.path(
  review_dir, "results", paste0(results_prefix, "-ALL-ROUTES-COMPLETE.rds")
)
atomic_save_rds(list(
  started_at = started_at, completed_at = Sys.time(),
  routes = route_results,
  package_source_fingerprint = package_source_fingerprint
), summary_path)
cat("All requested ACTG zero-inflated routes are stable: ", summary_path, "\n", sep = "")
