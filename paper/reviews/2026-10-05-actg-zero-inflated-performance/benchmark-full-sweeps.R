#!/usr/bin/env Rscript

script_argument <- grep("^--file=", commandArgs(trailingOnly = FALSE), value = TRUE)
script_path <- normalizePath(sub("^--file=", "", script_argument[[1L]]), winslash = "/", mustWork = TRUE)
root <- normalizePath(file.path(dirname(script_path), "..", "..", ".."), winslash = "/", mustWork = TRUE)
setwd(root)

suppressPackageStartupMessages(devtools::load_all(root, quiet = TRUE))
load(file.path(root, "data", "actg_cd4.rda"))
data <- actg_cd4
data$week_factor <- factor(data$week, levels = sort(unique(data$week)), ordered = FALSE)
data$visit1 <- as.integer(data$week == min(data$week, na.rm = TRUE))
data <- data[order(data$subject, data$week), , drop = FALSE]

rhs <- ~ treatment + s(age, bs = "ps", k = 8) + sex +
  s(baseline_cd4, bs = "ps", k = 8) + week_factor +
  treatment:week_factor + baseline_cd4:visit1
mu <- update(rhs, cd4 ~ .)
args <- commandArgs(TRUE)
copula <- if (length(args) >= 1L) args[[1L]] else "G"
max_inner <- if (length(args) >= 2L) as.integer(args[[2L]]) else 5L
max_outer <- if (length(args) >= 3L) as.integer(args[[3L]]) else 1L
outer_tol <- if (length(args) >= 4L) as.numeric(args[[4L]]) else 1e-20
start_step <- if (length(args) >= 6L) as.numeric(args[[6L]]) else 0.5
verbose <- if (length(args) >= 7L) as.integer(args[[7L]]) else 0L
warm_start <- if (length(args) >= 8L) as.logical(args[[8L]]) else FALSE
backtracking_halves <- if (length(args) >= 9L) as.integer(args[[9L]]) else 50L

control <- gamlss_longitudinal_control(
  outer_tol = outer_tol,
  max_outer_iter = max_outer,
  max_elapsed_sec = 1200,
  stop_on_convergence = TRUE,
  rs = list(
    inner_tol = 1e-20,
    max_inner_iter = max_inner,
    start_step_size = start_step,
    max_steps = 5L,
    warm_start_joint = warm_start,
    warm_start_joint_iter = if (warm_start) 5L else 0L,
    use_backtracking = TRUE,
    backtracking_max_halves = backtracking_halves,
    update_lambda = TRUE,
    discrete_score_method = "analytical"
  )
)

timing <- system.time(fit <- gamlss_longitudinal(
  data = data,
  time_var = "week_factor",
  subject_var = "subject",
  margin_dist = gamlss.dist::ZISICHEL(),
  copula_dist = copula,
  mu.formula = mu,
  sigma.formula = rhs,
  nu.formula = rhs,
  tau.formula = rhs,
  theta.formula = rhs,
  zeta.formula = ~ 1,
  missingness = "segment",
  method = "RS",
  optimization = "joint",
  optimizer_control = control,
  compute_vcov = FALSE,
  verbose = verbose
))
cat(sprintf(
  paste0("copula=%s inner=%d step=%.3f warm=%s outer_max=%d outer_done=%d elapsed=%.3f ",
         "logLik=%.15f trace_rows=%d converged=%s stop=%s\n"),
  copula, max_inner, start_step, warm_start, max_outer, fit$convergence$outer_iterations,
  timing[["elapsed"]], as.numeric(logLik(fit)), nrow(fit$rs_block_trace),
  fit$convergence$converged, fit$convergence$stop_reason
))
if (length(args) >= 5L && nzchar(args[[5L]])) {
  saveRDS(fit, args[[5L]], version = 3L)
}
