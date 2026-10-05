#!/usr/bin/env Rscript

script_argument <- grep("^--file=", commandArgs(trailingOnly = FALSE), value = TRUE)
script_path <- normalizePath(sub("^--file=", "", script_argument[[1L]]), winslash = "/", mustWork = TRUE)
root <- normalizePath(file.path(dirname(script_path), "..", "..", ".."), winslash = "/", mustWork = TRUE)
setwd(root)
suppressPackageStartupMessages(devtools::load_all(root, quiet = TRUE))

args <- commandArgs(TRUE)
input_path <- normalizePath(args[[1L]], winslash = "/", mustWork = TRUE)
output_path <- args[[2L]]
bfgs_maxit <- if (length(args) >= 3L) as.integer(args[[3L]]) else 30L
rs_max_outer <- if (length(args) >= 4L) as.integer(args[[4L]]) else 30L
backtracking_halves <- if (length(args) >= 5L) as.integer(args[[5L]]) else 2L
joint_method <- if (length(args) >= 6L) args[[6L]] else "BFGS"

input <- readRDS(input_path)
fit <- if (inherits(input, "gamlss.longitudinal")) input else input$value
stopifnot(inherits(fit, "gamlss.longitudinal"))

augmented <- gamlss.longitudinal:::.gl_build_cg_model(
  fit$model_matrix,
  fit$par,
  fit$par_s
)
beta_start <- augmented$beta
penalty <- gamlss.longitudinal:::.gl_build_cg_penalty(
  beta_names = names(beta_start),
  lambda_current = fit$lambda_s,
  par_s = fit$par_s,
  mm = fit$model_matrix
)
copula_link <- gamlss.longitudinal:::get_copula_dist(fit$copula_dist)$copula_link
dataset <- data.frame(
  response = fit$response,
  time = fit$response_margin,
  subject = fit$response_subject
)
pair_cache <- gamlss.longitudinal:::build_copula_pair_cache(
  fit$response,
  fit$response_margin,
  fit$response_subject
)
margin_eval_cache <- gamlss.longitudinal:::.build_margin_eval_cache(
  fit$margin_dist,
  calc_d2 = FALSE
)
evaluate <- function(beta) {
  gamlss.longitudinal:::.gl_evaluate_cg_beta(
    beta_vec = beta,
    mm_cg = augmented$mm,
    par_cov_template = fit$par,
    par_s_template = fit$par_s,
    margin_dist = fit$margin_dist,
    copula_link = copula_link,
    copula_dist = fit$copula_dist,
    dataset = dataset,
    pair_cache = pair_cache,
    margin_eval_cache = margin_eval_cache
  )
}

last_valid <- new.env(parent = emptyenv())
last_valid$beta <- beta_start
last_valid$evaluation <- evaluate(beta_start)
if (is.null(last_valid$evaluation)) stop("The BFGS starting state is not likelihood-valid.")

objective <- function(beta) {
  evaluation <- evaluate(beta)
  if (is.null(evaluation)) return(1e100)
  last_valid$beta <- beta
  last_valid$evaluation <- evaluation
  -as.numeric(
    evaluation$loglik - 0.5 * sum(beta * as.numeric(penalty %*% beta))
  )
}
gradient <- function(beta) {
  evaluation <- evaluate(beta)
  if (is.null(evaluation)) {
    stop("BFGS requested a gradient at an invalid likelihood state.", call. = FALSE)
  }
  score <- gamlss.longitudinal:::.cg_analytical_gradient(
    beta,
    augmented$mm,
    evaluation$eta_out,
    evaluation$calc_lik,
    fit$margin_dist,
    fit$copula_dist,
    isTRUE(fit$include_dlcopdpar),
    fit$response,
    fit$response_margin,
    fit$response_subject
  )
  -stats::setNames(as.numeric(score - penalty %*% beta), names(beta))
}

bfgs_timing <- system.time({
  joint_control <- if (identical(joint_method, "L-BFGS-B")) {
    list(maxit = bfgs_maxit, factr = 1e4, pgtol = 1e-5, lmm = 20L, trace = 1L, REPORT = 1L)
  } else {
    list(maxit = bfgs_maxit, reltol = 1e-12, trace = 1L, REPORT = 1L)
  }
  bfgs <- stats::optim(
    par = beta_start,
    fn = objective,
    gr = gradient,
    method = joint_method,
    control = joint_control
  )
})
if (!is.finite(bfgs$value) || bfgs$value >= 1e99) {
  stop("BFGS did not return a finite likelihood-valid state.")
}

unpacked <- gamlss.longitudinal:::.gl_unpack_cg_beta(
  bfgs$par,
  par_cov_template = fit$par,
  par_s_template = fit$par_s
)
bfgs_start <- fit
bfgs_start$par <- unpacked$par_cov
bfgs_start$par_s <- unpacked$par_s

rs_control <- gamlss_longitudinal_control(
  outer_tol = 1e-4,
  max_outer_iter = rs_max_outer,
  max_elapsed_sec = 1200,
  rs = list(
    inner_tol = 1e-5,
    max_inner_iter = 20L,
    initial_inner_iter = 5L,
    start_step_size = 1,
    warm_start_joint = TRUE,
    warm_start_joint_iter = 5L,
    use_backtracking = TRUE,
    backtracking_max_halves = backtracking_halves,
    adaptive_step_size = FALSE,
    update_lambda = TRUE,
    lambda_update_every = 10L,
    discrete_score_method = "analytical"
  )
)

fit_args <- list(
  data = fit$dataset_original,
  time_var = fit$time_var,
  subject_var = fit$subject_var,
  margin_dist = fit$margin_dist,
  copula_dist = fit$copula_dist,
  missingness = fit$missingness$requested,
  method = "RS",
  optimization = "joint",
  optimizer_control = rs_control,
  start = bfgs_start,
  compute_vcov = FALSE,
  verbose = 1L
)
for (parameter in intersect(
  c("mu", "sigma", "nu", "tau", "theta", "zeta"),
  names(fit$formulas)
)) {
  fit_args[[paste0(parameter, ".formula")]] <- fit$formulas[[parameter]]
}

rs_timing <- system.time(refined <- do.call(gamlss_longitudinal, fit_args))
result <- list(
  input_path = input_path,
  bfgs = bfgs,
  bfgs_timing = bfgs_timing,
  rs_timing = rs_timing,
  fit = refined
)
saveRDS(result, output_path, version = 3L)
cat(sprintf(
  paste0("joint_method=%s bfgs_elapsed=%.3f bfgs_convergence=%d bfgs_fn=%d bfgs_gr=%d ",
         "rs_elapsed=%.3f outer=%d logLik=%.15f converged=%s stop=%s\n"),
  joint_method, bfgs_timing[["elapsed"]], bfgs$convergence, bfgs$counts[["function"]],
  bfgs$counts[["gradient"]], rs_timing[["elapsed"]],
  refined$convergence$outer_iterations, as.numeric(logLik(refined)),
  refined$convergence$converged, refined$convergence$stop_reason
))
