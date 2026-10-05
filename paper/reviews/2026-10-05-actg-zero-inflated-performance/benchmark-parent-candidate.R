#!/usr/bin/env Rscript

script_argument <- grep("^--file=", commandArgs(trailingOnly = FALSE), value = TRUE)
script_path <- normalizePath(sub("^--file=", "", script_argument[[1L]]), winslash = "/", mustWork = TRUE)
root <- normalizePath(file.path(dirname(script_path), "..", "..", ".."), winslash = "/", mustWork = TRUE)
setwd(root)
suppressPackageStartupMessages(devtools::load_all(root, quiet = TRUE))

args <- commandArgs(TRUE)
parent_path <- normalizePath(args[[1L]], winslash = "/", mustWork = TRUE)
output_path <- args[[2L]]
update_lambda <- if (length(args) >= 3L) as.logical(args[[3L]]) else TRUE
max_inner <- if (length(args) >= 4L) as.integer(args[[4L]]) else 5L
lambda_update_every <- if (length(args) >= 5L) as.integer(args[[5L]]) else 1L
start_step <- if (length(args) >= 6L) as.numeric(args[[6L]]) else 0.5
method <- if (length(args) >= 7L) toupper(args[[7L]]) else "RS"
verbose <- if (length(args) >= 8L) as.integer(args[[8L]]) else 1L
backtracking_halves <- if (length(args) >= 9L) as.integer(args[[9L]]) else 50L
max_outer <- if (length(args) >= 10L) as.integer(args[[10L]]) else 100L
max_elapsed <- if (length(args) >= 11L) as.numeric(args[[11L]]) else 1200
initial_inner <- if (length(args) >= 12L) as.integer(args[[12L]]) else max_inner
adaptive_step <- if (length(args) >= 13L) as.logical(args[[13L]]) else FALSE
bfgs_refine <- if (length(args) >= 14L) as.logical(args[[14L]]) else FALSE
bfgs_after <- if (length(args) >= 15L) as.integer(args[[15L]]) else 10L
bfgs_max_iter <- if (length(args) >= 16L) as.integer(args[[16L]]) else 100L
model_mode <- if (length(args) >= 17L) args[[17L]] else "drop"
joint_method <- if (length(args) >= 21L) args[[21L]] else "BFGS"
joint_max_refinements <- if (length(args) >= 22L) as.integer(args[[22L]]) else 1L
outer_tol <- if (length(args) >= 23L) as.numeric(args[[23L]]) else 1e-4
coordinate_rescue <- if (length(args) >= 24L) as.logical(args[[24L]]) else FALSE
coordinate_repolish_iter <- if (length(args) >= 25L) as.integer(args[[25L]]) else 50L
start_path <- if (length(args) >= 20L) args[[20L]] else parent_path
cold_start <- identical(tolower(start_path), "cold")
parent_container <- readRDS(parent_path)
extract_fit <- function(container) {
  if (inherits(container, "gamlss.longitudinal")) return(container)
  if (is.list(container) && inherits(container$fit, "gamlss.longitudinal")) {
    return(container$fit)
  }
  if (is.list(container) && inherits(container$result$fit, "gamlss.longitudinal")) {
    return(container$result$fit)
  }
  if (is.list(container) && inherits(container$value, "gamlss.longitudinal")) {
    return(container$value)
  }
  NULL
}
parent <- extract_fit(parent_container)
stopifnot(inherits(parent, "gamlss.longitudinal"))
target_margin <- if (length(args) >= 18L) {
  gamlss.longitudinal:::.margin_family_object(args[[18L]])
} else {
  parent$margin_dist
}
target_copula <- if (length(args) >= 19L) args[[19L]] else parent$copula_dist
fit_start <- if (cold_start) {
  NULL
} else {
  start_path <- normalizePath(start_path, winslash = "/", mustWork = TRUE)
  extract_fit(readRDS(start_path))
}
if (!cold_start) stopifnot(inherits(fit_start, "gamlss.longitudinal"))
same_route <- identical(parent$margin_dist$family[[1L]], target_margin$family[[1L]]) &&
  identical(parent$copula_dist, target_copula)
same_margin <- identical(parent$margin_dist$family[[1L]], target_margin$family[[1L]])
copula_coefficient <- function(x) grepl("^(theta|zeta)\\.", names(x))
merge_nested_state <- function(template, supplied, zero_missing = FALSE) {
  output <- template
  for (parameter in names(output)) {
    if (isTRUE(zero_missing)) {
      output[[parameter]] <- lapply(output[[parameter]], function(value) value * 0)
    }
    if (!parameter %in% names(supplied)) next
    for (index in seq_along(output[[parameter]])) {
      if (index > length(supplied[[parameter]])) next
      candidate <- supplied[[parameter]][[index]]
      if (!is.null(candidate) && length(candidate) == length(output[[parameter]][[index]])) {
        output[[parameter]][[index]] <- candidate
      }
    }
  }
  output
}
if (cold_start) {
  merged_par <- NA
  merged_par_s <- NULL
  merged_lambda_s <- NULL
} else if (same_route) {
  merged_par <- parent$par * 0
  merged_par[names(fit_start$par)] <- fit_start$par
  merged_par_s <- merge_nested_state(parent$par_s, fit_start$par_s, zero_missing = TRUE)
  merged_lambda_s <- merge_nested_state(parent$lambda_s, fit_start$lambda_s)
} else if (same_margin) {
  merged_par <- parent$par
  merged_par[copula_coefficient(merged_par)] <- 0
  route_copula_par <- fit_start$par[copula_coefficient(fit_start$par)]
  merged_par[names(route_copula_par)] <- route_copula_par
  merged_par_s <- parent$par_s
  for (parameter in intersect(c("theta", "zeta"), names(merged_par_s))) {
    merged_par_s[[parameter]] <- lapply(
      merged_par_s[[parameter]],
      function(beta) beta * 0
    )
  }
  merged_lambda_s <- parent$lambda_s
} else {
  merged_par <- parent$par * 0
  merged_par[names(fit_start$par)] <- fit_start$par
  merged_par_s <- lapply(
    parent$par_s,
    function(block) lapply(block, function(beta) beta * 0)
  )
  merged_lambda_s <- parent$lambda_s
}
fit_start <- if (cold_start) NA else list(
  par = merged_par,
  par_s = merged_par_s,
  lambda_s = merged_lambda_s
)

formulas <- parent$formulas
drop_term <- "none"
if (!identical(model_mode, "maximal")) {
  terms_mu <- stats::terms(formulas$mu)
  labels <- attr(terms_mu, "term.labels")
  drop_term <- if ("sex" %in% labels) "sex" else labels[[length(labels)]]
  formulas$mu <- stats::reformulate(
    labels[labels != drop_term],
    response = deparse(formulas$mu[[2L]], width.cutoff = 500L),
    intercept = attr(terms_mu, "intercept") == 1L,
    env = environment(formulas$mu)
  )
}

control <- if (identical(method, "RS")) {
  gamlss_longitudinal_control(
    outer_tol = outer_tol,
    max_outer_iter = max_outer,
    max_elapsed_sec = max_elapsed,
    rs = list(
      inner_tol = 1e-5,
      max_inner_iter = max_inner,
      initial_inner_iter = initial_inner,
      start_step_size = start_step,
      warm_start_joint = TRUE,
      warm_start_joint_iter = 5L,
      use_backtracking = TRUE,
      backtracking_max_halves = backtracking_halves,
      adaptive_step_size = adaptive_step,
      update_lambda = update_lambda,
      lambda_update_every = lambda_update_every,
      bfgs_refine = bfgs_refine,
      bfgs_refine_after = bfgs_after,
      bfgs_method = joint_method,
      bfgs_max_refinements = joint_max_refinements,
      bfgs_max_iter = bfgs_max_iter,
      bfgs_reltol = 1e-12,
      bfgs_coordinate_rescue = coordinate_rescue,
      bfgs_repolish_iter = coordinate_repolish_iter,
      discrete_score_method = "analytical"
    )
  )
} else {
  gamlss_longitudinal_control(
    outer_tol = outer_tol,
    max_outer_iter = max_outer,
    max_elapsed_sec = max_elapsed,
    cg = list(
      grad_tol = 1e-3,
      step_tol = 1e-5,
      gradient_method = "analytical",
      hessian_method = "analytical",
      update_lambda = update_lambda,
      lambda_update_every = lambda_update_every,
      line_search = "best",
      use_backtracking = TRUE
    )
  )
}

timing <- system.time(fit <- gamlss_longitudinal(
  data = parent$dataset_original,
  time_var = parent$time_var,
  subject_var = parent$subject_var,
  margin_dist = target_margin,
  copula_dist = target_copula,
  mu.formula = formulas$mu,
  sigma.formula = formulas$sigma,
  nu.formula = formulas$nu,
  tau.formula = formulas$tau,
  theta.formula = formulas$theta,
  zeta.formula = formulas$zeta,
  missingness = parent$missingness$requested,
  method = method,
  optimization = "joint",
  optimizer_control = control,
  start = fit_start,
  compute_vcov = FALSE,
  verbose = verbose
))
saveRDS(fit, output_path, version = 3L)
cat(sprintf(
  paste0("drop=%s method=%s update_lambda=%s inner=%d initial_inner=%d lambda_every=%d step=%.3f ",
         "backtracks=%d adaptive_step=%s bfgs=%s elapsed=%.3f outer=%d rows=%d ",
         "logLik=%.15f converged=%s stop=%s\n"),
  drop_term, method, update_lambda, max_inner, initial_inner, lambda_update_every, start_step,
  backtracking_halves, adaptive_step, bfgs_refine,
  timing[["elapsed"]], fit$convergence$outer_iterations,
  nrow(fit$rs_block_trace), as.numeric(logLik(fit)), fit$convergence$converged,
  fit$convergence$stop_reason
))
