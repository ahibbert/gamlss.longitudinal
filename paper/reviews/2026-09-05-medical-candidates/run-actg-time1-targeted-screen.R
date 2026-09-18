#!/usr/bin/env Rscript

# Targeted, checkpointed follow-up to the preserved ACTG ZISICHEL--Gumbel fit.
# The search addresses the reproducible visit-1 PIT concentration without
# changing or overwriting the accepted model or its review bundle.

script_argument <- grep("^--file=", commandArgs(trailingOnly = FALSE), value = TRUE)
script_path <- if (length(script_argument)) {
  normalizePath(sub("^--file=", "", script_argument[[1L]]), mustWork = TRUE)
} else {
  normalizePath(file.path(
    getwd(), "paper", "reviews", "2026-09-05-medical-candidates",
    "run-actg-time1-targeted-screen.R"
  ), mustWork = TRUE)
}
package_root <- normalizePath(file.path(dirname(script_path), "..", "..", ".."))

argument_value <- function(name, default) {
  prefix <- paste0("--", name, "=")
  hit <- grep(paste0("^", prefix), commandArgs(trailingOnly = TRUE), value = TRUE)
  if (!length(hit)) return(default)
  sub(prefix, "", hit[[1L]], fixed = TRUE)
}
source_results <- argument_value(
  "source-results-dir",
  file.path(dirname(script_path), "results", "actg-revised-model-screen")
)
results_root <- argument_value(
  "output-dir",
  file.path(dirname(script_path), "results", "actg-time1-investigation", "targeted-screen")
)
base_fit_path <- argument_value("base-fit", file.path(source_results, "revised-final-fit.rds"))
base_formulas_path <- argument_value(
  "base-formulas", file.path(source_results, "revised-final-formulas.rds")
)
margin_family <- toupper(argument_value("margin", "ZISICHEL"))
if (!margin_family %in% c("ZISICHEL", "ZIBNB")) {
  stop("--margin must be ZISICHEL or ZIBNB.")
}
checkpoint_dir <- file.path(results_root, "checkpoints")
dir.create(checkpoint_dir, recursive = TRUE, showWarnings = FALSE)

if (!requireNamespace("devtools", quietly = TRUE)) stop("Install devtools.")
devtools::load_all(package_root, quiet = TRUE)

workers <- as.integer(argument_value("workers", "4"))
max_elapsed_sec <- as.numeric(argument_value("max-elapsed-sec", "7200"))
if (!is.finite(workers) || workers < 1L) stop("--workers must be positive.")

load_package_data <- function(name) {
  environment <- new.env(parent = emptyenv())
  utils::data(list = name, package = "gamlss.longitudinal", envir = environment)
  get(name, envir = environment, inherits = FALSE)
}
actg_data <- load_package_data("actg_cd4")
actg_data$week_factor <- factor(actg_data$week)
actg_data$week_index <- as.numeric(actg_data$week_factor) - 1
actg_data$visit1 <- as.integer(actg_data$week == min(actg_data$week))

base_fit <- readRDS(base_fit_path)
base_formulas <- lapply(
  readRDS(base_formulas_path),
  stats::as.formula
)

optimizer_control <- gamlss_longitudinal_control(
  outer_tol = 1e-4,
  max_outer_iter = 50L,
  max_elapsed_sec = max_elapsed_sec,
  rs = list(
    inner_tol = 1e-5, max_inner_iter = 100L,
    warm_start_joint = TRUE, warm_start_joint_iter = 5L,
    discrete_score_method = "analytical"
  )
)

formula_text <- function(x) paste(deparse(x, width.cutoff = 500L), collapse = " ")
safe_name <- function(x) gsub("(^-+|-+$)", "", tolower(gsub("[^A-Za-z0-9]+", "-", x)))
object_md5 <- function(x) {
  path <- tempfile(fileext = ".rds")
  on.exit(unlink(path), add = TRUE)
  saveRDS(x, path, version = 2)
  unname(tools::md5sum(path))
}
atomic_save_rds <- function(x, path) {
  temporary <- tempfile(tmpdir = dirname(path), fileext = ".rds")
  saveRDS(x, temporary, compress = "gzip")
  ok <- if (file.exists(path)) file.copy(temporary, path, overwrite = TRUE) else file.rename(temporary, path)
  if (file.exists(temporary)) unlink(temporary)
  if (!ok) stop("Could not write ", path)
}
atomic_write_csv <- function(x, path) {
  temporary <- tempfile(tmpdir = dirname(path), fileext = ".csv")
  utils::write.csv(x, temporary, row.names = FALSE, na = "")
  if (!file.copy(temporary, path, overwrite = TRUE)) stop("Could not write ", path)
  unlink(temporary)
}

cvm_uniform <- function(x) {
  x <- sort(x[is.finite(x)])
  n <- length(x)
  if (!n) return(NA_real_)
  1 / (12 * n) + sum((x - (2 * seq_len(n) - 1) / (2 * n))^2)
}
pit_diagnostics <- function(fit, seeds = 20260910L + 0:19) {
  observed <- is.finite(fit$dataset$response)
  times <- fit$dataset$time[observed]
  rows <- lapply(seeds, function(seed) {
    pit <- pithist(fit, randomize = TRUE, seed = seed, plot = FALSE)$pit
    if (length(pit) != length(times)) stop("PIT alignment failure.")
    do.call(rbind, lapply(sort(unique(times)), function(time) {
      value <- pit[times == time]
      data.frame(
        seed = seed, time = time, n = length(value), mean = mean(value), sd = stats::sd(value),
        sd_error = abs(stats::sd(value) - sqrt(1 / 12)), cvm = cvm_uniform(value),
        lower_10 = mean(value <= 0.1), upper_10 = mean(value >= 0.9)
      )
    }))
  })
  raw <- do.call(rbind, rows)
  by_time <- do.call(rbind, lapply(split(raw, raw$time), function(value) data.frame(
    time = value$time[[1L]], n = value$n[[1L]], pit_mean = mean(value$mean),
    pit_sd = mean(value$sd), pit_sd_error = mean(value$sd_error),
    pit_cvm = mean(value$cvm), lower_10 = mean(value$lower_10), upper_10 = mean(value$upper_10)
  )))
  rownames(by_time) <- NULL
  list(raw = raw, by_time = by_time)
}

make_task <- function(stage, label, formulas) {
  metadata <- list(
    stage = stage, label = label,
    formulas = vapply(formulas, formula_text, character(1L)),
    margin = margin_family, copula = "G", n_rows = nrow(actg_data),
    optimizer_control = optimizer_control
  )
  signature <- substr(object_md5(metadata), 1L, 12L)
  list(
    stage = stage, label = label, formulas = formulas, metadata = metadata,
    checkpoint = file.path(checkpoint_dir, paste0(safe_name(stage), "-", safe_name(label), "-", signature, ".rds"))
  )
}

fit_candidate <- function(task) {
  if (file.exists(task$checkpoint)) {
    cached <- tryCatch(readRDS(task$checkpoint), error = identity)
    if (is.list(cached) && identical(cached$metadata, task$metadata) && isTRUE(cached$converged)) {
      cached$from_checkpoint <- TRUE
      return(cached)
    }
  }
  started <- proc.time()[["elapsed"]]
  warnings <- character()
  fit <- tryCatch(withCallingHandlers(
    do.call(gamlss_longitudinal, list(
      data = actg_data, time_var = "week_factor", subject_var = "subject",
      margin_dist = getExportedValue("gamlss.dist", margin_family)(), copula_dist = "G",
      mu.formula = task$formulas$mu, sigma.formula = task$formulas$sigma,
      nu.formula = task$formulas$nu, tau.formula = task$formulas$tau,
      theta.formula = task$formulas$theta,
      missingness = "segment", method = "RS", optimization = "joint",
      optimizer_control = optimizer_control, compute_vcov = FALSE, verbose = 0
    )), warning = function(warning) {
      warnings <<- c(warnings, conditionMessage(warning))
      invokeRestart("muffleWarning")
    }
  ), error = identity)
  elapsed <- proc.time()[["elapsed"]] - started
  returned_fit <- inherits(fit, "gamlss.longitudinal")
  converged <- returned_fit && isTRUE(fit$convergence$converged)
  diagnostics <- if (converged) pit_diagnostics(fit) else NULL
  result <- list(
    metadata = task$metadata, fit = if (returned_fit) fit else NULL,
    converged = converged, elapsed_sec = elapsed,
    logLik = if (returned_fit) as.numeric(logLik(fit)) else NA_real_,
    AIC = if (returned_fit) AIC(fit) else NA_real_, BIC = if (returned_fit) BIC(fit) else NA_real_,
    EDF = if (returned_fit) as.numeric(fit$model_selection["EDF", "joint"]) else NA_real_,
    diagnostics = diagnostics,
    outer_iterations = if (returned_fit) fit$convergence$outer_iterations else NA_integer_,
    stop_reason = if (returned_fit) fit$convergence$stop_reason else "error",
    warnings = unique(warnings), error = if (inherits(fit, "error")) conditionMessage(fit) else NA_character_,
    checkpoint = task$checkpoint, from_checkpoint = FALSE
  )
  atomic_save_rds(result, task$checkpoint)
  result
}

result_row <- function(task, result) {
  time1 <- if (!is.null(result$diagnostics)) result$diagnostics$by_time[result$diagnostics$by_time$time == 1, ] else NULL
  data.frame(
    stage = task$stage, label = task$label,
    mu_formula = formula_text(task$formulas$mu), sigma_formula = formula_text(task$formulas$sigma),
    nu_formula = formula_text(task$formulas$nu), tau_formula = formula_text(task$formulas$tau),
    theta_formula = formula_text(task$formulas$theta),
    logLik = result$logLik, AIC = result$AIC, BIC = result$BIC, EDF = result$EDF,
    time1_pit_mean = if (nrow(time1)) time1$pit_mean else NA_real_,
    time1_pit_sd = if (nrow(time1)) time1$pit_sd else NA_real_,
    time1_pit_sd_error = if (nrow(time1)) time1$pit_sd_error else NA_real_,
    time1_pit_cvm = if (nrow(time1)) time1$pit_cvm else NA_real_,
    converged = result$converged, elapsed_sec = result$elapsed_sec,
    outer_iterations = result$outer_iterations, stop_reason = result$stop_reason,
    warnings = paste(result$warnings, collapse = " | "), error = result$error,
    checkpoint = result$checkpoint, from_checkpoint = result$from_checkpoint,
    stringsAsFactors = FALSE
  )
}

cluster <- NULL
if (workers > 1L) {
  cluster <- parallel::makePSOCKcluster(workers, outfile = file.path(results_root, "workers.log"), setup_strategy = "parallel")
  on.exit(parallel::stopCluster(cluster), add = TRUE)
  parallel::clusterCall(cluster, function(root) {
    setwd(root)
    Sys.setenv(OMP_NUM_THREADS = "1", OPENBLAS_NUM_THREADS = "1", MKL_NUM_THREADS = "1")
    devtools::load_all(root, quiet = TRUE)
    TRUE
  }, package_root)
  parallel::clusterExport(cluster, c(
    "actg_data", "optimizer_control", "cvm_uniform", "pit_diagnostics",
    "atomic_save_rds", "fit_candidate", "margin_family"
  ), envir = environment())
}
run_tasks <- function(tasks) {
  if (is.null(cluster)) lapply(tasks, fit_candidate) else parallel::parLapplyLB(cluster, tasks, fit_candidate)
}

additions <- list(
  mu_baseline_visit1 = function(f) {
    f$mu <- update(f$mu, . ~ . + baseline_cd4:visit1)
    f
  },
  sigma_visit1 = function(f) {
    f$sigma <- update(f$sigma, ~ . + visit1)
    f
  },
  sigma_baseline_visit1 = function(f) {
    f$sigma <- update(f$sigma, ~ . + visit1 + baseline_cd4:visit1)
    f
  },
  nu_visit1 = function(f) {
    f$nu <- update(f$nu, ~ . + visit1)
    f
  },
  tau_visit1 = function(f) {
    f$tau <- update(f$tau, ~ . + visit1)
    f
  },
  tau_baseline_visit1 = function(f) {
    f$tau <- update(f$tau, ~ . + visit1 + baseline_cd4:visit1)
    f
  }
)

base_task <- make_task("single-addition", "preserved-model", base_formulas)
base_diagnostics <- pit_diagnostics(base_fit)
base_result <- list(
  metadata = base_task$metadata, fit = base_fit, converged = isTRUE(base_fit$convergence$converged),
  elapsed_sec = 0, logLik = as.numeric(logLik(base_fit)), AIC = AIC(base_fit), BIC = BIC(base_fit),
  EDF = as.numeric(base_fit$model_selection["EDF", "joint"]), diagnostics = base_diagnostics,
  outer_iterations = base_fit$convergence$outer_iterations, stop_reason = base_fit$convergence$stop_reason,
  warnings = character(), error = NA_character_, checkpoint = file.path(checkpoint_dir, "preserved-model.rds"),
  from_checkpoint = TRUE
)
if (!file.exists(base_result$checkpoint)) atomic_save_rds(base_result, base_result$checkpoint)

stage1_tasks <- lapply(names(additions), function(label) make_task("single-addition", label, additions[[label]](base_formulas)))
cat("Fitting ", length(stage1_tasks), " single-addition candidates...\n", sep = "")
stage1_results <- run_tasks(stage1_tasks)
stage1_rows <- do.call(rbind, c(list(result_row(base_task, base_result)), Map(result_row, stage1_tasks, stage1_results)))
stage1_rows$delta_AIC <- stage1_rows$AIC - min(stage1_rows$AIC[stage1_rows$converged], na.rm = TRUE)
stage1_rows <- stage1_rows[order(stage1_rows$AIC), ]
stage1_rows$selected <- seq_len(nrow(stage1_rows)) == 1L
atomic_write_csv(stage1_rows, file.path(results_root, "single-addition-results.csv"))

winner_label <- stage1_rows$label[[1L]]
winner_result <- if (winner_label == "preserved-model") base_result else stage1_results[[match(winner_label, names(additions))]]
winner_formulas <- if (winner_label == "preserved-model") base_formulas else additions[[winner_label]](base_formulas)

# One greedy second step: add each remaining targeted correction to the AIC winner.
remaining <- setdiff(names(additions), winner_label)
stage2_tasks <- lapply(remaining, function(label) make_task(
  "second-addition", paste(winner_label, label, sep = "+"), additions[[label]](winner_formulas)
))
cat("Fitting ", length(stage2_tasks), " second-addition candidates from ", winner_label, "...\n", sep = "")
stage2_results <- run_tasks(stage2_tasks)
winner_task <- make_task("second-addition", winner_label, winner_formulas)
winner_result$checkpoint <- winner_result$checkpoint
stage2_rows <- do.call(rbind, c(list(result_row(winner_task, winner_result)), Map(result_row, stage2_tasks, stage2_results)))
stage2_rows$delta_AIC <- stage2_rows$AIC - min(stage2_rows$AIC[stage2_rows$converged], na.rm = TRUE)
stage2_rows <- stage2_rows[order(stage2_rows$AIC), ]
stage2_rows$selected <- seq_len(nrow(stage2_rows)) == 1L
atomic_write_csv(stage2_rows, file.path(results_root, "second-addition-results.csv"))

final_label <- stage2_rows$label[[1L]]
if (final_label == winner_label) {
  final_result <- winner_result
  final_formulas <- winner_formulas
} else {
  index <- match(final_label, vapply(stage2_tasks, `[[`, character(1L), "label"))
  final_result <- stage2_results[[index]]
  final_formulas <- stage2_tasks[[index]]$formulas
}

# Check whether allowing each visit its own count-shape parameter avoids
# improving visit 1 at the expense of visits 2--4, and whether baseline CD4
# explains residual shape heterogeneity.
nu_refinement_formulas <- list(
  nu_week_factor = within(final_formulas, nu <- ~ week_factor),
  nu_visit1_baseline = within(
    final_formulas,
    nu <- ~ baseline_cd4 + visit1 + baseline_cd4:visit1
  ),
  nu_week_factor_baseline = within(
    final_formulas,
    nu <- ~ week_factor + baseline_cd4
  ),
  nu_week_factor_baseline_visit1 = within(
    final_formulas,
    nu <- ~ week_factor + baseline_cd4 + baseline_cd4:visit1
  )
)
stage3_tasks <- lapply(names(nu_refinement_formulas), function(label) make_task(
  "nu-refinement", label, nu_refinement_formulas[[label]]
))
cat("Fitting ", length(stage3_tasks), " ", margin_family, " nu-refinement candidates...\n", sep = "")
stage3_results <- run_tasks(stage3_tasks)
stage2_winner_task <- make_task("nu-refinement", final_label, final_formulas)
stage3_rows <- do.call(rbind, c(
  list(result_row(stage2_winner_task, final_result)),
  Map(result_row, stage3_tasks, stage3_results)
))
stage3_rows$delta_AIC <- stage3_rows$AIC - min(stage3_rows$AIC[stage3_rows$converged], na.rm = TRUE)
stage3_rows <- stage3_rows[order(stage3_rows$AIC), ]
stage3_rows$selected <- seq_len(nrow(stage3_rows)) == 1L
atomic_write_csv(stage3_rows, file.path(results_root, "nu-refinement-results.csv"))

nu_winner_label <- stage3_rows$label[[1L]]
if (nu_winner_label != final_label) {
  index <- match(nu_winner_label, vapply(stage3_tasks, `[[`, character(1L), "label"))
  final_result <- stage3_results[[index]]
  final_formulas <- stage3_tasks[[index]]$formulas
  final_label <- nu_winner_label
}

atomic_save_rds(final_result$fit, file.path(results_root, "time1-revised-candidate-fit.rds"))
atomic_save_rds(final_formulas, file.path(results_root, "time1-revised-candidate-formulas.rds"))
all_rows <- rbind(stage1_rows, stage2_rows, stage3_rows)
atomic_write_csv(all_rows, file.path(results_root, "all-targeted-results.csv"))
if (!is.null(final_result$diagnostics)) {
  atomic_write_csv(final_result$diagnostics$by_time, file.path(results_root, "selected-pit-by-time.csv"))
  atomic_write_csv(final_result$diagnostics$raw, file.path(results_root, "selected-pit-by-time-seeds.csv"))
}
atomic_save_rds(list(
  completed_at = Sys.time(), preserved_fit = normalizePath(base_fit_path),
  selected_label = final_label, selected_formulas = final_formulas,
  AIC = final_result$AIC, BIC = final_result$BIC, source_AIC = AIC(base_fit), source_BIC = BIC(base_fit)
), file.path(results_root, "RUN-COMPLETE.rds"))

cat("\nTargeted visit-1 screen complete.\n")
cat("Selected: ", final_label, "\n", sep = "")
cat("AIC: ", final_result$AIC, " (preserved model ", AIC(base_fit), ")\n", sep = "")
cat("BIC: ", final_result$BIC, " (preserved model ", BIC(base_fit), ")\n", sep = "")
print(final_result$diagnostics$by_time)
