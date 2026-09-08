#!/usr/bin/env Rscript

# Bounded profiler for exact-discrete ACTG fits. This is intentionally outside
# the manuscript targets graph and never runs an unrestricted optimizer.

if (!requireNamespace("devtools", quietly = TRUE)) {
  stop("Install devtools to profile the package checkout.")
}
devtools::load_all(quiet = TRUE)

argument_value <- function(name, default) {
  prefix <- paste0("--", name, "=")
  hit <- grep(paste0("^", prefix), commandArgs(trailingOnly = TRUE), value = TRUE)
  if (!length(hit)) return(default)
  sub(prefix, "", hit[[1L]], fixed = TRUE)
}

n_subjects <- as.integer(argument_value("subjects", "150"))
model <- match.arg(argument_value("model", "full"), c("minimal", "full"))
run_mode <- match.arg(argument_value("run", "bounded"), c("bounded", "converge"))
seed <- as.integer(argument_value("seed", "20260908"))
if (!is.finite(n_subjects) || n_subjects < 20L) {
  stop("--subjects must be an integer of at least 20.")
}

environment <- new.env(parent = emptyenv())
utils::data("actg_cd4", package = "gamlss.longitudinal", envir = environment)
data <- get("actg_cd4", envir = environment)
data$week_factor <- factor(data$week)
data$baseline_log_z <- as.numeric(scale(log1p(data$baseline_cd4)))
data$age_z <- as.numeric(scale(data$age))

subject_data <- data[!duplicated(data$subject), c("subject", "treatment")]
subject_groups <- split(subject_data$subject, subject_data$treatment)
set.seed(seed)
if (n_subjects >= nrow(subject_data)) {
  selected <- subject_data$subject
} else {
  allocation <- rep(n_subjects %/% length(subject_groups), length(subject_groups))
  allocation[seq_len(n_subjects %% length(subject_groups))] <-
    allocation[seq_len(n_subjects %% length(subject_groups))] + 1L
  selected <- unlist(Map(
    function(ids, size) sample(ids, min(size, length(ids))),
    subject_groups,
    allocation
  ), use.names = FALSE)
}
data <- data[data$subject %in% selected, , drop = FALSE]

formulas <- if (identical(model, "minimal")) {
  list(
    mu = cd4 ~ week_factor,
    sigma = ~ 1,
    nu = ~ 1,
    theta = ~ 1
  )
} else {
  list(
    mu = cd4 ~ week_factor + treatment + baseline_log_z + age_z + sex +
      treatment:week_factor,
    sigma = ~ treatment + baseline_log_z,
    nu = ~ 1,
    theta = ~ treatment
  )
}

control <- if (identical(run_mode, "bounded")) {
  gamlss_longitudinal_control(
    outer_tol = 1e-8,
    max_outer_iter = 1L,
    max_elapsed_sec = 600,
    stop_on_convergence = FALSE,
    rs = list(
      max_inner_iter = 1L,
      max_steps = 1L,
      warm_start_joint = FALSE,
      use_backtracking = FALSE,
      update_lambda = FALSE,
      discrete_score_method = "analytical"
    )
  )
} else {
  gamlss_longitudinal_control(
    outer_tol = 1e-4,
    max_outer_iter = 50L,
    max_elapsed_sec = 600,
    rs = list(
      inner_tol = 1e-5,
      max_inner_iter = 100L,
      warm_start_joint = TRUE,
      warm_start_joint_iter = 5L,
      discrete_score_method = "analytical"
    )
  )
}

output_dir <- file.path(
  "paper", "reviews", "2026-09-08-discrete-runtime", "results"
)
dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)
stem <- sprintf("actg-%s-%s-n%d", model, run_mode, length(selected))
profile_path <- file.path(output_dir, paste0(stem, ".out"))

warnings <- character()
started <- proc.time()[["elapsed"]]
Rprof(profile_path, interval = 0.005, line.profiling = TRUE)
fit <- tryCatch(
  withCallingHandlers(
    gamlss_longitudinal(
      dataset = data,
      time_var = "week_factor",
      subject_var = "subject",
      margin_dist = gamlss.dist::DEL(),
      copula_dist = "N",
      mu.formula = formulas$mu,
      sigma.formula = formulas$sigma,
      nu.formula = formulas$nu,
      theta.formula = formulas$theta,
      missingness = "segment",
      optimizer_control = control,
      compute_vcov = FALSE,
      verbose = 0
    ),
    warning = function(w) {
      warnings <<- c(warnings, conditionMessage(w))
      invokeRestart("muffleWarning")
    }
  ),
  error = identity
)
Rprof(NULL)
elapsed <- proc.time()[["elapsed"]] - started

profile <- summaryRprof(profile_path, lines = "show")
utils::write.csv(
  utils::head(profile$by.total, 40L),
  file.path(output_dir, paste0(stem, "-by-total.csv"))
)
utils::write.csv(
  utils::head(profile$by.self, 40L),
  file.path(output_dir, paste0(stem, "-by-self.csv"))
)
summary <- data.frame(
  model = model,
  run_mode = run_mode,
  subjects = length(selected),
  rows = nrow(data),
  observed = sum(is.finite(data$cd4)),
  elapsed_sec = elapsed,
  returned_fit = inherits(fit, "gamlss.longitudinal"),
  error = if (inherits(fit, "error")) conditionMessage(fit) else "",
  warnings = paste(unique(warnings), collapse = " | "),
  stringsAsFactors = FALSE
)
utils::write.csv(summary, file.path(output_dir, paste0(stem, "-summary.csv")), row.names = FALSE)
print(summary)
print(utils::head(profile$by.self, 20L))
