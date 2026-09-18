#!/usr/bin/env Rscript

# Final narrow polishing screen around the shortlisted nu ~ visit1 model.
# The reference fit and its full report are never overwritten.

script_argument <- grep("^--file=", commandArgs(trailingOnly = FALSE), value = TRUE)
script_path <- if (length(script_argument)) {
  normalizePath(sub("^--file=", "", script_argument[[1L]]), mustWork = TRUE)
} else {
  normalizePath("paper/reviews/actg-shortlist/run-time1-polish-screen.R", mustWork = TRUE)
}
review_root <- dirname(script_path)
package_root <- normalizePath(file.path(review_root, "..", "..", ".."))
results_root <- file.path(review_root, "time1-polish-screen")
checkpoint_dir <- file.path(results_root, "checkpoints")
dir.create(checkpoint_dir, recursive = TRUE, showWarnings = FALSE)

if (!requireNamespace("devtools", quietly = TRUE)) stop("Install devtools.")
devtools::load_all(package_root, quiet = TRUE)

argument_value <- function(name, default) {
  prefix <- paste0("--", name, "=")
  hit <- grep(paste0("^", prefix), commandArgs(trailingOnly = TRUE), value = TRUE)
  if (!length(hit)) return(default)
  sub(prefix, "", hit[[1L]], fixed = TRUE)
}
workers <- as.integer(argument_value("workers", "4"))
max_elapsed_sec <- as.numeric(argument_value("max-elapsed-sec", "7200"))

load_package_data <- function(name) {
  environment <- new.env(parent = emptyenv())
  utils::data(list = name, package = "gamlss.longitudinal", envir = environment)
  get(name, envir = environment, inherits = FALSE)
}
actg_data <- load_package_data("actg_cd4")
actg_data$week_factor <- factor(actg_data$week)
actg_data$week_index <- as.numeric(actg_data$week_factor) - 1
actg_data$visit1 <- as.integer(actg_data$week == min(actg_data$week))

reference_fit <- readRDS(file.path(
  review_root, "03-mu-baseline-visit1-nu-visit1", "objects", "final-model.rds"
))
reference_formulas <- reference_fit$formulas

optimizer_control <- gamlss_longitudinal_control(
  outer_tol = 1e-4, max_outer_iter = 50L, max_elapsed_sec = max_elapsed_sec,
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
  1 / (12 * n) + sum((x - (2 * seq_len(n) - 1) / (2 * n))^2)
}
worm_metrics <- function(pit) {
  z <- sort(stats::qnorm(pmin(pmax(pit, .Machine$double.eps), 1 - .Machine$double.eps)))
  n <- length(z)
  theoretical <- stats::qnorm((seq_len(n) - 0.5) / n)
  detrended <- z - theoretical
  central <- abs(theoretical) <= 2.5
  model <- stats::lm(detrended[central] ~ theoretical[central] + I(theoretical[central]^2) + I(theoretical[central]^3))
  coefficients <- stats::coef(model)
  c(
    worm_rms = sqrt(mean(detrended[central]^2)),
    worm_intercept = coefficients[[1L]], worm_linear = coefficients[[2L]],
    worm_quadratic = coefficients[[3L]], worm_cubic = coefficients[[4L]]
  )
}
fit_diagnostics <- function(fit, seeds = 20260910L + 0:19) {
  observed <- is.finite(fit$dataset$response)
  times <- fit$dataset$time[observed]
  values <- lapply(seeds, function(seed) {
    pit <- pithist(fit, randomize = TRUE, seed = seed, plot = FALSE)$pit
    if (length(pit) != length(times)) stop("PIT alignment failure.")
    do.call(rbind, lapply(sort(unique(times)), function(time) {
      value <- pit[times == time]
      data.frame(
        seed = seed, time = time, n = length(value), pit_mean = mean(value),
        pit_sd = stats::sd(value), pit_cvm = cvm_uniform(value),
        lower_10 = mean(value <= 0.1), upper_10 = mean(value >= 0.9),
        as.list(worm_metrics(value)), check.names = FALSE
      )
    }))
  })
  raw <- do.call(rbind, values)
  by_time <- do.call(rbind, lapply(split(raw, raw$time), function(value) {
    numeric_columns <- setdiff(names(value), c("seed", "time", "n"))
    result <- data.frame(time = value$time[[1L]], n = value$n[[1L]])
    for (column in numeric_columns) result[[column]] <- mean(value[[column]])
    result
  }))
  rownames(by_time) <- NULL
  list(raw = raw, by_time = by_time)
}

make_task <- function(label, formulas) {
  metadata <- list(
    label = label, formulas = vapply(formulas, formula_text, character(1L)),
    margin = "ZISICHEL", copula = "G", optimizer_control = optimizer_control,
    n_rows = nrow(actg_data)
  )
  signature <- substr(object_md5(metadata), 1L, 12L)
  list(
    label = label, formulas = formulas, metadata = metadata,
    checkpoint = file.path(checkpoint_dir, paste0(safe_name(label), "-", signature, ".rds"))
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
    gamlss_longitudinal(
      data = actg_data, time_var = "week_factor", subject_var = "subject",
      margin_dist = gamlss.dist::ZISICHEL(), copula_dist = "G",
      mu.formula = task$formulas$mu, sigma.formula = task$formulas$sigma,
      nu.formula = task$formulas$nu, tau.formula = task$formulas$tau,
      theta.formula = task$formulas$theta, zeta.formula = task$formulas$zeta,
      missingness = "segment", method = "RS", optimization = "joint",
      optimizer_control = optimizer_control, compute_vcov = FALSE, verbose = 0
    ), warning = function(warning) {
      warnings <<- c(warnings, conditionMessage(warning))
      invokeRestart("muffleWarning")
    }
  ), error = identity)
  elapsed <- proc.time()[["elapsed"]] - started
  returned_fit <- inherits(fit, "gamlss.longitudinal")
  converged <- returned_fit && isTRUE(fit$convergence$converged)
  result <- list(
    metadata = task$metadata, fit = if (returned_fit) fit else NULL,
    converged = converged, elapsed_sec = elapsed,
    logLik = if (returned_fit) as.numeric(logLik(fit)) else NA_real_,
    AIC = if (returned_fit) AIC(fit) else NA_real_, BIC = if (returned_fit) BIC(fit) else NA_real_,
    EDF = if (returned_fit) as.numeric(fit$model_selection["EDF", "joint"]) else NA_real_,
    diagnostics = if (converged) fit_diagnostics(fit) else NULL,
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
    label = task$label, mu_formula = formula_text(task$formulas$mu),
    sigma_formula = formula_text(task$formulas$sigma), nu_formula = formula_text(task$formulas$nu),
    tau_formula = formula_text(task$formulas$tau), theta_formula = formula_text(task$formulas$theta),
    logLik = result$logLik, AIC = result$AIC, BIC = result$BIC, EDF = result$EDF,
    time1_pit_mean = if (nrow(time1)) time1$pit_mean else NA_real_,
    time1_pit_sd = if (nrow(time1)) time1$pit_sd else NA_real_,
    time1_pit_cvm = if (nrow(time1)) time1$pit_cvm else NA_real_,
    time1_worm_rms = if (nrow(time1)) time1$worm_rms else NA_real_,
    time1_worm_intercept = if (nrow(time1)) time1$worm_intercept else NA_real_,
    time1_worm_linear = if (nrow(time1)) time1$worm_linear else NA_real_,
    time1_worm_quadratic = if (nrow(time1)) time1$worm_quadratic else NA_real_,
    time1_worm_cubic = if (nrow(time1)) time1$worm_cubic else NA_real_,
    converged = result$converged, elapsed_sec = result$elapsed_sec,
    warnings = paste(result$warnings, collapse = " | "), error = result$error,
    checkpoint = result$checkpoint, from_checkpoint = result$from_checkpoint,
    stringsAsFactors = FALSE
  )
}

cluster <- parallel::makePSOCKcluster(
  workers, outfile = file.path(results_root, "workers.log"), setup_strategy = "parallel"
)
on.exit(parallel::stopCluster(cluster), add = TRUE)
parallel::clusterCall(cluster, function(root) {
  setwd(root)
  Sys.setenv(OMP_NUM_THREADS = "1", OPENBLAS_NUM_THREADS = "1", MKL_NUM_THREADS = "1")
  devtools::load_all(root, quiet = TRUE)
  TRUE
}, package_root)
parallel::clusterExport(cluster, c(
  "actg_data", "optimizer_control", "cvm_uniform", "worm_metrics",
  "fit_diagnostics", "atomic_save_rds", "fit_candidate"
), envir = environment())

modify_formula <- function(formulas, parameter, formula) {
  formulas[[parameter]] <- formula
  formulas
}
specifications <- list(
  sigma_visit1 = modify_formula(
    reference_formulas, "sigma", update(reference_formulas$sigma, ~ . + visit1)
  ),
  sigma_baseline_visit1 = modify_formula(
    reference_formulas, "sigma",
    update(reference_formulas$sigma, ~ . + visit1 + baseline_cd4:visit1)
  ),
  tau_visit1 = modify_formula(
    reference_formulas, "tau", update(reference_formulas$tau, ~ . + visit1)
  ),
  tau_baseline_visit1 = modify_formula(
    reference_formulas, "tau",
    update(reference_formulas$tau, ~ . + visit1 + baseline_cd4:visit1)
  ),
  mu_quadratic_baseline_visit1 = modify_formula(
    reference_formulas, "mu",
    update(reference_formulas$mu, . ~ . + I(baseline_cd4^2):visit1)
  )
)

reference_task <- make_task("reference_nu_visit1", reference_formulas)
reference_result <- list(
  metadata = reference_task$metadata, fit = reference_fit,
  converged = isTRUE(reference_fit$convergence$converged), elapsed_sec = 0,
  logLik = as.numeric(logLik(reference_fit)), AIC = AIC(reference_fit), BIC = BIC(reference_fit),
  EDF = as.numeric(reference_fit$model_selection["EDF", "joint"]), diagnostics = fit_diagnostics(reference_fit),
  outer_iterations = reference_fit$convergence$outer_iterations,
  stop_reason = reference_fit$convergence$stop_reason, warnings = character(), error = NA_character_,
  checkpoint = file.path(checkpoint_dir, "reference-nu-visit1.rds"), from_checkpoint = TRUE
)
if (!file.exists(reference_result$checkpoint)) atomic_save_rds(reference_result, reference_result$checkpoint)

tasks <- lapply(names(specifications), function(label) make_task(label, specifications[[label]]))
cat("Fitting ", length(tasks), " final polishing candidates...\n", sep = "")
results <- parallel::parLapplyLB(cluster, tasks, fit_candidate)
rows <- do.call(rbind, c(list(result_row(reference_task, reference_result)), Map(result_row, tasks, results)))
rows$delta_AIC_from_reference <- rows$AIC - reference_result$AIC
rows$delta_BIC_from_reference <- rows$BIC - reference_result$BIC
rows <- rows[order(rows$AIC), ]
rows$rank_AIC <- seq_len(nrow(rows))
atomic_write_csv(rows, file.path(results_root, "polish-results.csv"))

best_label <- rows$label[[1L]]
if (best_label == reference_task$label) {
  best_result <- reference_result
  best_formulas <- reference_formulas
} else {
  index <- match(best_label, vapply(tasks, `[[`, character(1L), "label"))
  best_result <- results[[index]]
  best_formulas <- tasks[[index]]$formulas
}
atomic_save_rds(best_result$fit, file.path(results_root, "aic-best-polish-fit.rds"))
atomic_save_rds(best_formulas, file.path(results_root, "aic-best-polish-formulas.rds"))
if (!is.null(best_result$diagnostics)) {
  atomic_write_csv(best_result$diagnostics$by_time, file.path(results_root, "aic-best-diagnostics-by-time.csv"))
}
atomic_save_rds(list(
  completed_at = Sys.time(), reference_AIC = reference_result$AIC,
  reference_BIC = reference_result$BIC, best_label = best_label,
  best_AIC = best_result$AIC, best_BIC = best_result$BIC
), file.path(results_root, "RUN-COMPLETE.rds"))

cat("Polishing screen complete. Best AIC: ", best_label, "\n", sep = "")
print(rows[, c(
  "label", "AIC", "BIC", "delta_AIC_from_reference", "delta_BIC_from_reference",
  "time1_pit_sd", "time1_pit_cvm", "time1_worm_rms", "converged"
)], row.names = FALSE)
