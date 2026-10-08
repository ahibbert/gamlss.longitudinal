#!/usr/bin/env Rscript

# Refit one exact backward-selection candidate from its fitted parent while
# varying only the bounded joint BFGS schedule. This is diagnostic tooling for
# route-default selection; it does not alter paper checkpoints.

script_argument <- grep("^--file=", commandArgs(trailingOnly = FALSE), value = TRUE)
script_path <- normalizePath(sub("^--file=", "", script_argument[[1L]]),
                             winslash = "/", mustWork = TRUE)
root <- normalizePath(file.path(dirname(script_path), "..", "..", ".."),
                      winslash = "/", mustWork = TRUE)
setwd(root)
suppressPackageStartupMessages(devtools::load_all(root, quiet = TRUE))

args <- commandArgs(trailingOnly = TRUE)
if (length(args) < 5L) {
  stop(
    "Usage: diagnose-specific-candidate.R PARENT OUTPUT PARAMETER TERM BFGS_MAXIT ",
    "[MAX_ELAPSED] [MAX_REFINEMENTS] [START_MODE]",
    call. = FALSE
  )
}
parent_path <- normalizePath(args[[1L]], winslash = "/", mustWork = TRUE)
output_path <- args[[2L]]
parameter <- args[[3L]]
term <- args[[4L]]
bfgs_max_iter <- as.integer(args[[5L]])
max_elapsed_sec <- if (length(args) >= 6L) as.numeric(args[[6L]]) else 1800
bfgs_max_refinements <- if (length(args) >= 7L) as.integer(args[[7L]]) else 1L
start_mode <- if (length(args) >= 8L) args[[8L]] else "parent"
if (!start_mode %in% c("parent", "zero-parameter")) {
  stop("START_MODE must be parent or zero-parameter.", call. = FALSE)
}

container <- readRDS(parent_path)
parent <- if (inherits(container, "gamlss.longitudinal")) {
  container
} else if (is.list(container) && inherits(container$fit, "gamlss.longitudinal")) {
  container$fit
} else {
  NULL
}
stopifnot(inherits(parent, "gamlss.longitudinal"))
if (!parameter %in% names(parent$formulas)) {
  stop("Unknown formula parameter: ", parameter, call. = FALSE)
}

remove_formula_term <- function(formula, term) {
  formula_terms <- stats::terms(formula)
  labels <- attr(formula_terms, "term.labels")
  if (!term %in% labels) stop("Term is not present in the parent formula: ", term)
  response <- if (attr(formula_terms, "response") == 1L) {
    deparse(formula[[2L]], width.cutoff = 500L)
  } else {
    NULL
  }
  stats::reformulate(
    labels[labels != term],
    response = response,
    intercept = attr(formula_terms, "intercept") == 1L,
    env = environment(formula)
  )
}

formulas <- parent$formulas
formulas[[parameter]] <- remove_formula_term(formulas[[parameter]], term)
fit_start <- parent
if (identical(start_mode, "zero-parameter")) {
  fixed_names <- grep(paste0("^", parameter, "\\."), names(fit_start$par), value = TRUE)
  fit_start$par[fixed_names] <- 0
  if (parameter %in% names(fit_start$par_s)) {
    fit_start$par_s[[parameter]] <- lapply(
      fit_start$par_s[[parameter]],
      function(value) value * 0
    )
  }
}
control <- gamlss_longitudinal_control(
  outer_tol = 1e-4,
  max_outer_iter = 50L,
  max_elapsed_sec = max_elapsed_sec,
  rs = list(
    inner_tol = 1e-5,
    max_inner_iter = 20L,
    initial_inner_iter = 2L,
    start_step_size = 2,
    warm_start_joint = TRUE,
    warm_start_joint_iter = 2L,
    use_backtracking = TRUE,
    backtracking_max_halves = 2L,
    adaptive_step_size = FALSE,
    update_lambda = TRUE,
    lambda_update_every = 5L,
    bfgs_refine = TRUE,
    bfgs_refine_after = 1L,
    bfgs_method = "BFGS",
    bfgs_max_refinements = bfgs_max_refinements,
    bfgs_max_iter = bfgs_max_iter,
    bfgs_reltol = 1e-12,
    discrete_score_method = "analytical"
  )
)

fit_args <- list(
  data = parent$dataset_original,
  time_var = parent$time_var,
  subject_var = parent$subject_var,
  margin_dist = parent$margin_dist,
  copula_dist = parent$copula_dist,
  missingness = parent$missingness$requested,
  method = "RS",
  optimization = "joint",
  optimizer_control = control,
  start = fit_start,
  compute_vcov = FALSE,
  verbose = 0
)
for (name in intersect(
  c("mu", "sigma", "nu", "tau", "theta", "zeta"),
  names(formulas)
)) {
  fit_args[[paste0(name, ".formula")]] <- formulas[[name]]
}

warnings <- character()
timing <- system.time({
  fit <- withCallingHandlers(
    do.call(gamlss_longitudinal, fit_args),
    warning = function(warning) {
      warnings <<- c(warnings, conditionMessage(warning))
      invokeRestart("muffleWarning")
    }
  )
})
result <- list(
  parent_path = parent_path,
  parameter = parameter,
  dropped_term = term,
  bfgs_max_iter = bfgs_max_iter,
  bfgs_max_refinements = bfgs_max_refinements,
  start_mode = start_mode,
  elapsed_sec = unname(timing[["elapsed"]]),
  warnings = unique(warnings),
  fit = fit
)
saveRDS(result, output_path, version = 3L)
d <- fit$rs_bfgs_diagnostics
cat(sprintf(
  paste0(
    "elapsed=%.2f min outer=%d converged=%s logLik=%.12f AIC=%.12f ",
    "bfgs_convergence=%d bfgs_fn=%d bfgs_gr=%d cache_hits=%d\n"
  ),
  result$elapsed_sec / 60,
  fit$convergence$outer_iterations,
  fit$convergence$converged,
  as.numeric(fit$model_selection["LogLik", "joint"]),
  as.numeric(fit$model_selection["AIC", "joint"]),
  d$convergence,
  d$counts[["function"]],
  d$counts[["gradient"]],
  d$evaluation_cache_hits
))
