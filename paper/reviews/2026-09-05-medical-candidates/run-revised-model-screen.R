#!/usr/bin/env Rscript

# Compact follow-up screen for the ACTG model. The stages deliberately remain
# sequential, while independent candidates within each stage run in parallel:
#   1. baseline-CD4 functional forms in mu and theta;
#   2. targeted time effects in mu, sigma, and nu;
#   3. alternative count margins, including ZISICHEL;
#   4. alternative copula families.
# Every candidate is checkpointed and the AIC winner is carried forward.

script_argument <- grep("^--file=", commandArgs(trailingOnly = FALSE), value = TRUE)
script_path <- if (length(script_argument)) {
  normalizePath(sub("^--file=", "", script_argument[[1L]]), mustWork = TRUE)
} else {
  normalizePath(file.path(
    getwd(), "paper", "reviews", "2026-09-05-medical-candidates",
    "run-revised-model-screen.R"
  ), mustWork = TRUE)
}
package_root <- normalizePath(file.path(dirname(script_path), "..", "..", ".."))
results_root <- file.path(dirname(script_path), "results", "actg-revised-model-screen")
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

argument_flag <- function(name, default) {
  value <- tolower(argument_value(name, if (default) "true" else "false"))
  if (!value %in% c("true", "false", "1", "0", "yes", "no")) {
    stop("--", name, " must be true or false.")
  }
  value %in% c("true", "1", "yes")
}

workers <- as.integer(argument_value("workers", "4"))
subject_argument <- tolower(argument_value("subjects", "all"))
resume <- argument_flag("resume", TRUE)
dry_run <- argument_flag("dry-run", FALSE)
max_elapsed_sec <- as.numeric(argument_value("max-elapsed-sec", "7200"))
if (!is.finite(workers) || workers < 1L) stop("--workers must be positive.")
if (!is.finite(max_elapsed_sec) || max_elapsed_sec <= 0) {
  stop("--max-elapsed-sec must be positive.")
}

load_package_data <- function(name) {
  environment <- new.env(parent = emptyenv())
  utils::data(list = name, package = "gamlss.longitudinal", envir = environment)
  get(name, envir = environment, inherits = FALSE)
}

actg_data <- load_package_data("actg_cd4")
actg_data$week_factor <- factor(actg_data$week)
actg_data$week_index <- as.numeric(actg_data$week_factor) - 1
actg_data$baseline_log_z <- as.numeric(scale(log1p(actg_data$baseline_cd4)))
actg_data$age_z <- as.numeric(scale(actg_data$age))
if (!identical(subject_argument, "all")) {
  n_subjects <- as.integer(subject_argument)
  if (!is.finite(n_subjects) || n_subjects < 1L) stop("Invalid --subjects value.")
  set.seed(20260910L)
  ids <- sample(unique(actg_data$subject), min(n_subjects, length(unique(actg_data$subject))))
  actg_data <- actg_data[actg_data$subject %in% ids, , drop = FALSE]
}

optimizer_control <- gamlss_longitudinal_control(
  outer_tol = 1e-4,
  max_outer_iter = 50L,
  max_elapsed_sec = max_elapsed_sec,
  rs = list(
    inner_tol = 1e-5,
    max_inner_iter = 100L,
    warm_start_joint = TRUE,
    warm_start_joint_iter = 5L,
    discrete_score_method = "analytical"
  )
)

formula_text <- function(x) paste(deparse(x, width.cutoff = 500L), collapse = " ")
formula_signature <- function(x) vapply(x, formula_text, character(1L))
safe_name <- function(x) {
  x <- tolower(gsub("[^A-Za-z0-9]+", "-", x))
  gsub("(^-+|-+$)", "", x)
}
object_md5 <- function(x) {
  path <- tempfile("actg-revision-signature-", fileext = ".rds")
  on.exit(unlink(path), add = TRUE)
  saveRDS(x, path, version = 2)
  unname(tools::md5sum(path))
}
atomic_save_rds <- function(x, path) {
  temporary <- tempfile("candidate-", tmpdir = dirname(path), fileext = ".rds")
  saveRDS(x, temporary, compress = "gzip")
  copied <- if (file.exists(path)) {
    file.copy(temporary, path, overwrite = TRUE)
  } else {
    file.rename(temporary, path)
  }
  if (file.exists(temporary)) unlink(temporary)
  if (!copied) stop("Could not write ", path)
  invisible(path)
}
atomic_write_csv <- function(x, path) {
  temporary <- tempfile("table-", tmpdir = dirname(path), fileext = ".csv")
  utils::write.csv(x, temporary, row.names = FALSE, na = "")
  if (!file.copy(temporary, path, overwrite = TRUE)) stop("Could not write ", path)
  unlink(temporary)
  invisible(path)
}

source_files <- c(
  file.path(package_root, "DESCRIPTION"),
  list.files(file.path(package_root, "R"), pattern = "\\.[Rr]$", full.names = TRUE)
)
source_fingerprint <- tools::md5sum(source_files)

base_formula_path <- file.path(
  dirname(results_root), "actg-covariate-refinement-final-formulas.rds"
)
base_fit_path <- file.path(
  dirname(results_root), "actg-covariate-refinement-final-fit.rds"
)
if (!file.exists(base_formula_path) || !file.exists(base_fit_path)) {
  stop("The refined ACTG formulas and fit must exist before this screen.")
}
base_formulas <- readRDS(base_formula_path)
base_formulas$tau <- base_formulas$tau %||% ~ 1
base_formulas$zeta <- base_formulas$zeta %||% ~ 1

empty_metrics <- function() {
  c(logLik = NA_real_, AIC = NA_real_, BIC = NA_real_, EDF = NA_real_)
}

fit_diagnostics <- function(fit, seed = 20260910L) {
  pit_object <- tryCatch(
    gamlss.longitudinal:::.gl_pit(fit, seed = seed),
    error = function(error) NULL
  )
  pit <- if (is.null(pit_object)) numeric() else pit_object$pit
  pit <- pit[is.finite(pit)]
  lag_summary <- tryCatch(
    gamlss.longitudinal:::.gl_residual_dependence_summary(fit, 1:3),
    error = function(error) data.frame()
  )
  lag_value <- function(lag) {
    if (!nrow(lag_summary) || !lag %in% lag_summary$lag) return(NA_real_)
    lag_summary$normal_score_cor[match(lag, lag_summary$lag)]
  }
  c(
    pit_mean = if (length(pit)) mean(pit) else NA_real_,
    pit_sd = if (length(pit) > 1L) stats::sd(pit) else NA_real_,
    pit_ks_p = if (length(pit)) stats::ks.test(pit, "punif")$p.value else NA_real_,
    pit_lower_10 = if (length(pit)) mean(pit <= 0.1) else NA_real_,
    pit_upper_10 = if (length(pit)) mean(pit >= 0.9) else NA_real_,
    residual_cor_lag1 = lag_value(1L),
    residual_cor_lag2 = lag_value(2L),
    residual_cor_lag3 = lag_value(3L)
  )
}

make_task <- function(stage, label, formulas, margin, copula) {
  metadata <- list(
    stage = stage,
    label = label,
    formulas = formula_signature(formulas),
    margin = margin,
    copula = copula,
    missingness = "segment",
    optimizer_control = optimizer_control,
    source_fingerprint = source_fingerprint,
    r_version = R.version.string,
    n_rows = nrow(actg_data),
    n_subjects = length(unique(actg_data$subject)),
    n_observed = sum(is.finite(actg_data$cd4))
  )
  signature <- substr(object_md5(metadata), 1L, 12L)
  list(
    stage = stage,
    label = label,
    formulas = formulas,
    margin = margin,
    copula = copula,
    metadata = metadata,
    checkpoint = file.path(
      checkpoint_dir,
      paste0(safe_name(stage), "-", safe_name(label), "-", signature, ".rds")
    )
  )
}

fit_candidate <- function(task) {
  if (isTRUE(resume) && file.exists(task$checkpoint)) {
    cached <- tryCatch(readRDS(task$checkpoint), error = identity)
    if (is.list(cached) && identical(cached$metadata, task$metadata) &&
        isTRUE(cached$converged) && inherits(cached$fit, "gamlss.longitudinal")) {
      cached$from_checkpoint <- TRUE
      return(cached)
    }
  }

  started <- proc.time()[["elapsed"]]
  warnings <- character()
  margin_function <- get(task$margin, envir = asNamespace("gamlss.dist"), mode = "function")
  fit <- tryCatch(
    withCallingHandlers(
      gamlss_longitudinal(
        data = actg_data,
        time_var = "week_factor",
        subject_var = "subject",
        margin_dist = margin_function(),
        copula_dist = task$copula,
        mu.formula = task$formulas$mu,
        sigma.formula = task$formulas$sigma,
        nu.formula = task$formulas$nu,
        tau.formula = task$formulas$tau,
        theta.formula = task$formulas$theta,
        zeta.formula = task$formulas$zeta,
        missingness = "segment",
        method = "RS",
        optimization = "joint",
        optimizer_control = optimizer_control,
        compute_vcov = FALSE,
        verbose = 0
      ),
      warning = function(warning) {
        warnings <<- c(warnings, conditionMessage(warning))
        invokeRestart("muffleWarning")
      }
    ),
    error = identity
  )
  elapsed <- proc.time()[["elapsed"]] - started
  returned_fit <- inherits(fit, "gamlss.longitudinal")
  converged <- returned_fit && isTRUE(fit$convergence$converged)
  metrics <- empty_metrics()
  diagnostics <- rep(NA_real_, 8L)
  names(diagnostics) <- c(
    "pit_mean", "pit_sd", "pit_ks_p", "pit_lower_10", "pit_upper_10",
    "residual_cor_lag1", "residual_cor_lag2", "residual_cor_lag3"
  )
  if (returned_fit && !is.null(fit$model_selection)) {
    metrics <- c(
      logLik = as.numeric(fit$model_selection["LogLik", "joint"]),
      AIC = as.numeric(fit$model_selection["AIC", "joint"]),
      BIC = as.numeric(fit$model_selection["BIC", "joint"]),
      EDF = as.numeric(fit$model_selection["EDF", "joint"])
    )
    if (converged) diagnostics <- fit_diagnostics(fit)
  }
  result <- list(
    metadata = task$metadata,
    fit = if (returned_fit) fit else NULL,
    metrics = metrics,
    diagnostics = diagnostics,
    converged = converged,
    elapsed_sec = elapsed,
    outer_iterations = if (returned_fit) fit$convergence$outer_iterations else NA_integer_,
    stop_reason = if (returned_fit) fit$convergence$stop_reason else "error",
    warnings = unique(warnings),
    error = if (inherits(fit, "error")) conditionMessage(fit) else NA_character_,
    checkpoint = task$checkpoint,
    from_checkpoint = FALSE
  )
  atomic_save_rds(result, task$checkpoint)
  result
}

result_row <- function(task, result) {
  values <- c(result$metrics, result$diagnostics)
  data.frame(
    stage = task$stage,
    label = task$label,
    margin = task$margin,
    copula = task$copula,
    mu_formula = formula_text(task$formulas$mu),
    sigma_formula = formula_text(task$formulas$sigma),
    nu_formula = formula_text(task$formulas$nu),
    tau_formula = formula_text(task$formulas$tau),
    theta_formula = formula_text(task$formulas$theta),
    zeta_formula = formula_text(task$formulas$zeta),
    as.list(values),
    converged = result$converged,
    elapsed_sec = result$elapsed_sec,
    outer_iterations = result$outer_iterations,
    stop_reason = result$stop_reason,
    from_checkpoint = result$from_checkpoint,
    warnings = paste(result$warnings, collapse = " | "),
    error = result$error,
    checkpoint = task$checkpoint,
    stringsAsFactors = FALSE,
    check.names = FALSE
  )
}

cluster <- NULL
if (workers > 1L && !dry_run) {
  cluster <- parallel::makePSOCKcluster(
    workers,
    outfile = file.path(results_root, "workers.log"),
    setup_strategy = "parallel"
  )
  on.exit(parallel::stopCluster(cluster), add = TRUE)
  parallel::clusterCall(cluster, function(root) {
    setwd(root)
    Sys.setenv(OMP_NUM_THREADS = "1", OPENBLAS_NUM_THREADS = "1", MKL_NUM_THREADS = "1")
    devtools::load_all(root, quiet = TRUE)
    TRUE
  }, package_root)
  parallel::clusterExport(
    cluster,
    c(
      "actg_data", "optimizer_control", "resume", "empty_metrics",
      "fit_diagnostics", "atomic_save_rds", "fit_candidate"
    ),
    envir = environment()
  )
}

run_tasks <- function(tasks) {
  if (dry_run) return(vector("list", length(tasks)))
  if (is.null(cluster)) lapply(tasks, fit_candidate) else parallel::parLapplyLB(cluster, tasks, fit_candidate)
}

all_rows <- list()
stage_winners <- list()
run_stage <- function(stage, specifications, reuse = list()) {
  tasks <- lapply(specifications, function(specification) {
    make_task(
      stage = stage,
      label = specification$label,
      formulas = specification$formulas,
      margin = specification$margin,
      copula = specification$copula
    )
  })
  cat("\n=== ", stage, " (", length(tasks), " candidates) ===\n", sep = "")
  if (dry_run) {
    for (task in tasks) {
      cat(task$label, ": ", task$margin, "/", task$copula, "\n", sep = "")
      print(task$formulas)
    }
    return(NULL)
  }
  results <- vector("list", length(tasks))
  labels <- vapply(tasks, `[[`, character(1L), "label")
  reused <- labels %in% names(reuse)
  if (any(reused)) {
    for (ii in which(reused)) {
      results[[ii]] <- reuse[[labels[[ii]]]]
      tasks[[ii]]$checkpoint <- results[[ii]]$checkpoint
    }
  }
  if (any(!reused)) {
    fitted <- run_tasks(tasks[!reused])
    results[!reused] <- fitted
  }
  rows <- do.call(rbind, Map(result_row, tasks, results))
  rows$rank_AIC <- NA_integer_
  eligible <- which(rows$converged & is.finite(rows$AIC))
  if (!length(eligible)) stop("No candidate converged in stage ", stage, ".")
  rows$rank_AIC[eligible] <- rank(rows$AIC[eligible], ties.method = "first")
  rows <- rows[order(rows$rank_AIC, rows$AIC, na.last = TRUE), , drop = FALSE]
  winner_index <- match(1L, rows$rank_AIC)
  winner_checkpoint <- rows$checkpoint[[winner_index]]
  winner_result <- readRDS(winner_checkpoint)
  winner_task <- tasks[[match(rows$label[[winner_index]], vapply(tasks, `[[`, character(1L), "label"))]]
  rows$selected <- rows$rank_AIC == 1L
  atomic_write_csv(rows, file.path(results_root, paste0(safe_name(stage), "-results.csv")))
  all_rows[[stage]] <<- rows
  stage_winners[[stage]] <<- list(task = winner_task, result = winner_result)
  progress <- do.call(rbind, all_rows)
  atomic_write_csv(progress, file.path(results_root, "all-candidates.csv"))
  atomic_save_rds(
    list(updated_at = Sys.time(), completed_stages = names(stage_winners), winners = stage_winners),
    file.path(results_root, "progress.rds")
  )
  cat("Winner: ", rows$label[[winner_index]], " | AIC=", format(rows$AIC[[winner_index]], digits = 12),
      " | BIC=", format(rows$BIC[[winner_index]], digits = 12), "\n", sep = "")
  stage_winners[[stage]]
}

current_fit <- readRDS(base_fit_path)
base_checkpoint <- file.path(checkpoint_dir, "00-refined-current-fit.rds")
base_task <- make_task("smooth", "current-k8-k8", base_formulas, "ZIPIG", "G")
base_task$checkpoint <- base_checkpoint
base_result <- list(
  metadata = base_task$metadata,
  fit = current_fit,
  metrics = c(
    logLik = as.numeric(logLik(current_fit)), AIC = AIC(current_fit),
    BIC = BIC(current_fit), EDF = current_fit$model_selection["EDF", "joint"]
  ),
  diagnostics = fit_diagnostics(current_fit),
  converged = isTRUE(current_fit$convergence$converged),
  elapsed_sec = 0,
  outer_iterations = current_fit$convergence$outer_iterations,
  stop_reason = current_fit$convergence$stop_reason,
  warnings = character(), error = NA_character_, checkpoint = base_checkpoint,
  from_checkpoint = TRUE
)
if (!dry_run && !file.exists(base_checkpoint)) atomic_save_rds(base_result, base_checkpoint)

with_formula <- function(formulas, parameter, formula) {
  formulas[[parameter]] <- formula
  formulas
}
mu_k4 <- cd4 ~ treatment + age + s(baseline_cd4, bs = "ps", k = 4) + week_factor
theta_k4 <- ~ treatment + s(baseline_cd4, bs = "ps", k = 4) + week_factor
mu_log <- cd4 ~ treatment + age + baseline_log_z + week_factor
theta_log <- ~ treatment + baseline_log_z + week_factor

smooth_specs <- list(
  list(label = "current-k8-k8", formulas = base_formulas, margin = "ZIPIG", copula = "G"),
  list(label = "mu-k4", formulas = with_formula(base_formulas, "mu", mu_k4), margin = "ZIPIG", copula = "G"),
  list(label = "theta-k4", formulas = with_formula(base_formulas, "theta", theta_k4), margin = "ZIPIG", copula = "G"),
  list(label = "both-k4", formulas = with_formula(with_formula(base_formulas, "mu", mu_k4), "theta", theta_k4), margin = "ZIPIG", copula = "G"),
  list(label = "mu-log", formulas = with_formula(base_formulas, "mu", mu_log), margin = "ZIPIG", copula = "G"),
  list(label = "theta-log", formulas = with_formula(base_formulas, "theta", theta_log), margin = "ZIPIG", copula = "G"),
  list(label = "both-log", formulas = with_formula(with_formula(base_formulas, "mu", mu_log), "theta", theta_log), margin = "ZIPIG", copula = "G")
)

if (dry_run) {
  run_stage("smooth", smooth_specs)
  cat("\nDry run stops before data-dependent later stages.\n")
  quit(save = "no", status = 0L)
}

# Reuse the existing refined fit rather than fitting the unchanged candidate.
smooth_other <- smooth_specs[vapply(smooth_specs, `[[`, character(1L), "label") != "current-k8-k8"]
smooth_tasks <- lapply(smooth_other, function(specification) make_task(
  "smooth", specification$label, specification$formulas, specification$margin, specification$copula
))
cat("\n=== smooth (7 candidates; current fit reused) ===\n")
smooth_results <- c(list(base_result), run_tasks(smooth_tasks))
smooth_tasks <- c(list(base_task), smooth_tasks)
smooth_rows <- do.call(rbind, Map(result_row, smooth_tasks, smooth_results))
eligible <- which(smooth_rows$converged & is.finite(smooth_rows$AIC))
smooth_rows$rank_AIC <- NA_integer_
smooth_rows$rank_AIC[eligible] <- rank(smooth_rows$AIC[eligible], ties.method = "first")
smooth_rows <- smooth_rows[order(smooth_rows$rank_AIC, smooth_rows$AIC, na.last = TRUE), , drop = FALSE]
smooth_rows$selected <- smooth_rows$rank_AIC == 1L
atomic_write_csv(smooth_rows, file.path(results_root, "smooth-results.csv"))
smooth_winner_checkpoint <- smooth_rows$checkpoint[match(1L, smooth_rows$rank_AIC)]
smooth_winner_result <- readRDS(smooth_winner_checkpoint)
smooth_winner_task <- smooth_tasks[[match(smooth_rows$label[match(1L, smooth_rows$rank_AIC)], vapply(smooth_tasks, `[[`, character(1L), "label"))]]
all_rows$smooth <- smooth_rows
stage_winners$smooth <- list(task = smooth_winner_task, result = smooth_winner_result)
cat("Winner: ", smooth_winner_task$label, " | AIC=", format(AIC(smooth_winner_result$fit), digits = 12), "\n", sep = "")

f <- smooth_winner_task$formulas
covariate_specs <- list(
  list(label = "smooth-winner", formulas = f, margin = "ZIPIG", copula = "G"),
  list(label = "sigma-week", formulas = with_formula(f, "sigma", ~ treatment + baseline_cd4 + week_factor), margin = "ZIPIG", copula = "G"),
  list(label = "nu-week", formulas = with_formula(f, "nu", ~ treatment + sex + baseline_cd4 + week_factor), margin = "ZIPIG", copula = "G"),
  list(label = "sigma-nu-week", formulas = with_formula(with_formula(f, "sigma", ~ treatment + baseline_cd4 + week_factor), "nu", ~ treatment + sex + baseline_cd4 + week_factor), margin = "ZIPIG", copula = "G"),
  list(label = "mu-treatment-linear-time", formulas = with_formula(f, "mu", update(f$mu, . ~ . + treatment:week_index)), margin = "ZIPIG", copula = "G")
)
covariate_winner <- run_stage(
  "targeted-covariates", covariate_specs,
  reuse = list(`smooth-winner` = smooth_winner_result)
)

f <- covariate_winner$task$formulas
margin_names <- c("ZIPIG", "ZISICHEL", "ZINBI", "DEL", "NBI")
margin_specs <- lapply(margin_names, function(margin) {
  margin_formulas <- f
  if (identical(margin, "ZISICHEL")) {
    # ZIPIG's nu is the zero-inflation probability. In ZISICHEL, nu is the
    # Sichel shape and tau is the zero-inflation probability, so transfer the
    # selected ZIPIG nu predictor to tau and estimate the new shape initially
    # as an intercept. This makes the comparison structurally fair.
    margin_formulas$nu <- ~ 1
    margin_formulas$tau <- f$nu
  }
  list(
    label = paste0("margin-", tolower(margin)), formulas = margin_formulas,
    margin = margin, copula = "G"
  )
})
margin_winner <- run_stage(
  "margin", margin_specs,
  reuse = list(`margin-zipig` = covariate_winner$result)
)

f <- margin_winner$task$formulas
margin <- margin_winner$task$margin
copula_names <- c("G", "N", "F", "C", "J", "t")
copula_specs <- lapply(copula_names, function(copula) list(
  label = paste0("copula-", tolower(copula)), formulas = f, margin = margin, copula = copula
))
copula_winner <- run_stage(
  "copula", copula_specs,
  reuse = list(`copula-g` = margin_winner$result)
)

final_fit <- copula_winner$result$fit
final_formulas <- copula_winner$task$formulas
final_fit_path <- file.path(results_root, "revised-final-fit.rds")
final_formula_path <- file.path(results_root, "revised-final-formulas.rds")
atomic_save_rds(final_fit, final_fit_path)
atomic_save_rds(final_formulas, final_formula_path)
all_candidates <- do.call(rbind, all_rows)
atomic_write_csv(all_candidates, file.path(results_root, "all-candidates.csv"))
winner_rows <- all_candidates[all_candidates$selected %in% TRUE, , drop = FALSE]
atomic_write_csv(winner_rows, file.path(results_root, "stage-winners.csv"))
atomic_save_rds(list(
  completed_at = Sys.time(), final_fit = final_fit_path, final_formulas = final_formula_path,
  selected_margin = copula_winner$task$margin, selected_copula = copula_winner$task$copula,
  AIC = AIC(final_fit), BIC = BIC(final_fit), stage_winners = stage_winners,
  source_fingerprint = source_fingerprint
), file.path(results_root, "RUN-COMPLETE.rds"))

cat("\nRevised screen complete.\n")
cat("Margin/copula: ", copula_winner$task$margin, "/", copula_winner$task$copula, "\n", sep = "")
print(final_formulas)
cat("AIC: ", format(AIC(final_fit), digits = 12), "\nBIC: ", format(BIC(final_fit), digits = 12), "\n", sep = "")
