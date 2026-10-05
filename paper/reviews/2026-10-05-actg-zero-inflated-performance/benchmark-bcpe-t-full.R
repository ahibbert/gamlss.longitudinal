#!/usr/bin/env Rscript

script_argument <- grep("^--file=", commandArgs(trailingOnly = FALSE), value = TRUE)
script_path <- normalizePath(sub("^--file=", "", script_argument[[1L]]), winslash = "/", mustWork = TRUE)
root <- normalizePath(file.path(dirname(script_path), "..", "..", ".."), winslash = "/", mustWork = TRUE)
setwd(root)
suppressPackageStartupMessages(devtools::load_all(root, quiet = TRUE))

args <- commandArgs(TRUE)
template_container <- readRDS(args[[1L]])
template <- if (inherits(template_container, "gamlss.longitudinal")) {
  template_container
} else if (inherits(template_container$result$fit, "gamlss.longitudinal")) {
  template_container$result$fit
} else {
  template_container$fit %||% template_container$value
}
output_path <- args[[2L]]
stopifnot(inherits(template, "gamlss.longitudinal"))

data <- template$dataset_original
response_name <- all.vars(template$formulas$mu)[[1L]]
data[[response_name]] <- ifelse(
  is.na(data[[response_name]]), NA_real_, data[[response_name]] + 1
)
control <- gamlss_longitudinal_control(
  outer_tol = 1e-4,
  max_outer_iter = 100L,
  max_elapsed_sec = 1200,
  rs = list(
    inner_tol = 1e-5,
    max_inner_iter = 20L,
    initial_inner_iter = 5L,
    start_step_size = 0.5,
    warm_start_joint = FALSE,
    warm_start_joint_iter = 5L,
    use_backtracking = TRUE,
    backtracking_max_halves = 10L,
    adaptive_step_size = FALSE,
    update_lambda = TRUE,
    lambda_update_every = 5L,
    bfgs_refine = FALSE,
    bfgs_refine_after = 2L,
    bfgs_method = "BFGS",
    bfgs_max_refinements = 1L,
    bfgs_max_iter = 50L,
    bfgs_reltol = 1e-12,
    discrete_score_method = "analytical"
  )
)

timing <- system.time(fit <- gamlss_longitudinal(
  data = data,
  time_var = template$time_var,
  subject_var = template$subject_var,
  margin_dist = gamlss.dist::BCPE(),
  copula_dist = "t",
  mu.formula = template$formulas$mu,
  sigma.formula = template$formulas$sigma,
  nu.formula = template$formulas$nu,
  tau.formula = template$formulas$tau,
  theta.formula = template$formulas$theta,
  zeta.formula = template$formulas$zeta,
  missingness = template$missingness$requested,
  method = "RS",
  optimization = "joint",
  optimizer_control = control,
  compute_vcov = FALSE,
  verbose = 1L
))
saveRDS(fit, output_path, version = 3L)
cat(sprintf(
  "BCPE/t full elapsed=%.3f outer=%d rows=%d logLik=%.15f converged=%s stop=%s\n",
  timing[["elapsed"]], fit$convergence$outer_iterations,
  nrow(fit$rs_block_trace), as.numeric(logLik(fit)), fit$convergence$converged,
  fit$convergence$stop_reason
))
