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
    paste0(
      "[MAX_ELAPSED] [MAX_REFINEMENTS] [START_MODE] [MAX_OUTER] ",
      paste0(
        "[BFGS_METHOD] [BFGS_PARSCALE_MAX] [MAX_INNER] [BFGS_REFINE] ",
        paste0(
          "[ADAPTIVE_STEP] [BFGS_ACTIVE_PARAMETERS] ",
          paste0(
            "[NEWTON_AUDIT_RDS] [NEWTON_STEP] [REPLACEMENT_TERM] ",
            "[BFGS_REFINE_AFTER] [INITIAL_INNER_ITER]"
          )
        )
      )
    ),
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
max_outer_iter <- if (length(args) >= 9L) as.integer(args[[9L]]) else 50L
bfgs_method <- if (length(args) >= 10L) args[[10L]] else "BFGS"
bfgs_parscale_max <- if (length(args) >= 11L) as.numeric(args[[11L]]) else 1
max_inner_iter <- if (length(args) >= 12L) as.integer(args[[12L]]) else 20L
bfgs_refine <- if (length(args) >= 13L) {
  identical(tolower(args[[13L]]), "true")
} else {
  TRUE
}
adaptive_step_size <- if (length(args) >= 14L) {
  identical(tolower(args[[14L]]), "true")
} else {
  FALSE
}
bfgs_active_parameters <- if (length(args) >= 15L &&
    !identical(tolower(args[[15L]]), "all")) {
  trimws(strsplit(args[[15L]], ",", fixed = TRUE)[[1L]])
} else {
  NULL
}
newton_audit_path <- if (length(args) >= 16L && nzchar(args[[16L]]) &&
    !identical(tolower(args[[16L]]), "none")) {
  normalizePath(args[[16L]], winslash = "/", mustWork = TRUE)
} else {
  NULL
}
newton_step <- if (length(args) >= 17L) as.numeric(args[[17L]]) else 1
replacement_term <- if (length(args) >= 18L && nzchar(args[[18L]]) &&
    !identical(tolower(args[[18L]]), "none")) {
  args[[18L]]
} else {
  NULL
}
bfgs_refine_after <- if (length(args) >= 19L) as.integer(args[[19L]]) else 1L
if (length(bfgs_refine_after) != 1L || is.na(bfgs_refine_after) ||
    bfgs_refine_after < 1L) {
  stop("BFGS_REFINE_AFTER must be a positive integer.", call. = FALSE)
}
initial_inner_iter <- if (length(args) >= 20L) as.integer(args[[20L]]) else 2L
if (length(initial_inner_iter) != 1L || is.na(initial_inner_iter) ||
    initial_inner_iter < 1L || initial_inner_iter > max_inner_iter) {
  stop(
    "INITIAL_INNER_ITER must be a positive integer no larger than MAX_INNER.",
    call. = FALSE
  )
}
if (length(newton_step) != 1L || !is.finite(newton_step) ||
    newton_step <= 0 || newton_step > 1) {
  stop("NEWTON_STEP must be in (0, 1].", call. = FALSE)
}
if (!start_mode %in% c(
    "parent", "mapped-list", "zero-parameter", "project-parameter")) {
  stop(
    paste0(
      "START_MODE must be parent, mapped-list, zero-parameter, ",
      "or project-parameter."
    ),
    call. = FALSE
  )
}
if (!bfgs_method %in% c("BFGS", "L-BFGS-B")) {
  stop("BFGS_METHOD must be BFGS or L-BFGS-B.", call. = FALSE)
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

replace_formula_term <- function(formula, term, replacement) {
  formula_terms <- stats::terms(formula)
  labels <- attr(formula_terms, "term.labels")
  if (!term %in% labels) stop("Term is not present in the parent formula: ", term)
  labels[labels == term] <- replacement
  labels <- unique(labels)
  response <- if (attr(formula_terms, "response") == 1L) {
    deparse(formula[[2L]], width.cutoff = 500L)
  } else {
    NULL
  }
  stats::reformulate(
    labels,
    response = response,
    intercept = attr(formula_terms, "intercept") == 1L,
    env = environment(formula)
  )
}

formulas <- parent$formulas
if (!identical(term, "__none__")) {
  formulas[[parameter]] <- if (is.null(replacement_term)) {
    remove_formula_term(formulas[[parameter]], term)
  } else {
    replace_formula_term(formulas[[parameter]], term, replacement_term)
  }
}
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
} else if (identical(start_mode, "project-parameter")) {
  fixed_design <- as.matrix(fit_start$model_matrix$x[[parameter]])
  fixed_names <- paste(parameter, colnames(fixed_design), sep = ".")
  smooth_design <- fit_start$model_matrix$s[[parameter]]
  smooth_coefficients <- fit_start$par_s[[parameter]]
  eta <- as.numeric(fixed_design %*% as.numeric(fit_start$par[fixed_names]))
  smooth_names <- intersect(names(smooth_design), names(smooth_coefficients))
  for (name in smooth_names) {
    eta <- eta + as.numeric(
      as.matrix(smooth_design[[name]]) %*% as.numeric(smooth_coefficients[[name]])
    )
  }
  retained_smooths <- setdiff(smooth_names, term)
  retained_designs <- lapply(retained_smooths, function(name) {
    as.matrix(smooth_design[[name]])
  })
  candidate_design <- do.call(cbind, c(list(fixed_design), retained_designs))
  projected <- stats::lm.fit(candidate_design, eta)$coefficients
  if (any(!is.finite(projected))) {
    stop("Projected candidate start was rank deficient.", call. = FALSE)
  }
  offset <- ncol(fixed_design)
  fit_start$par[fixed_names] <- projected[seq_len(offset)]
  for (name in retained_smooths) {
    width <- ncol(as.matrix(smooth_design[[name]]))
    index <- offset + seq_len(width)
    value <- projected[index]
    names(value) <- names(fit_start$par_s[[parameter]][[name]])
    fit_start$par_s[[parameter]][[name]] <- value
    offset <- offset + width
  }
}
if (identical(start_mode, "mapped-list")) {
  # A plain parameter list exercises the package's name-based coefficient and
  # compatible-smooth mapping without fitted-predictor projection. This is the
  # legacy fitted-parent behavior needed for an exact start-strategy A/B test.
  fit_start <- list(
    par = fit_start$par,
    par_s = fit_start$par_s,
    lambda_s = fit_start$lambda_s
  )
}
newton_delta <- NULL
if (!is.null(newton_audit_path)) {
  audit <- readRDS(newton_audit_path)
  hessian <- audit$hessian[names(fit_start$par), names(fit_start$par), drop = FALSE]
  information <- -(hessian + t(hessian)) / 2
  gradient <- audit$analytical_gradient[names(fit_start$par)]
  newton_delta <- tryCatch(
    as.numeric(solve(information, gradient)),
    error = function(error) stop("Could not solve audited Newton step: ", conditionMessage(error))
  )
  names(newton_delta) <- names(fit_start$par)
  fit_start$par <- fit_start$par + newton_step * newton_delta
}
control <- gamlss_longitudinal_control(
  outer_tol = 1e-4,
  max_outer_iter = max_outer_iter,
  max_elapsed_sec = max_elapsed_sec,
  rs = list(
    inner_tol = 1e-5,
    max_inner_iter = max_inner_iter,
    initial_inner_iter = initial_inner_iter,
    start_step_size = 2,
    warm_start_joint = TRUE,
    warm_start_joint_iter = 2L,
    use_backtracking = TRUE,
    backtracking_max_halves = 2L,
    adaptive_step_size = adaptive_step_size,
    update_lambda = TRUE,
    lambda_update_every = 5L,
    bfgs_refine = bfgs_refine,
    bfgs_refine_after = bfgs_refine_after,
    bfgs_method = bfgs_method,
    bfgs_max_refinements = bfgs_max_refinements,
    bfgs_max_iter = bfgs_max_iter,
    bfgs_reltol = 1e-12,
    bfgs_parscale_max = bfgs_parscale_max,
    bfgs_active_parameters = bfgs_active_parameters,
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
  bfgs_method = bfgs_method,
  bfgs_parscale_max = bfgs_parscale_max,
  bfgs_active_parameters = bfgs_active_parameters,
  max_inner_iter = max_inner_iter,
  bfgs_refine = bfgs_refine,
  bfgs_refine_after = bfgs_refine_after,
  initial_inner_iter = initial_inner_iter,
  adaptive_step_size = adaptive_step_size,
  start_mode = start_mode,
  replacement_term = replacement_term,
  newton_audit_path = newton_audit_path,
  newton_step = if (is.null(newton_audit_path)) NULL else newton_step,
  newton_delta = newton_delta,
  elapsed_sec = unname(timing[["elapsed"]]),
  warnings = unique(warnings),
  fit = fit
)
saveRDS(result, output_path, version = 3L)
d <- fit$rs_bfgs_diagnostics
cat(sprintf(
  paste0(
    "method=%s elapsed=%.2f min outer=%d converged=%s logLik=%.12f AIC=%.12f ",
    "bfgs_convergence=%d bfgs_fn=%d bfgs_gr=%d cache_hits=%d\n"
  ),
  result$bfgs_method,
  result$elapsed_sec / 60,
  fit$convergence$outer_iterations,
  fit$convergence$converged,
  as.numeric(fit$model_selection["LogLik", "joint"]),
  as.numeric(fit$model_selection["AIC", "joint"]),
  d$convergence %||% NA_integer_,
  (d$counts %||% stats::setNames(NA_integer_, "function"))[["function"]],
  (d$counts %||% stats::setNames(NA_integer_, "gradient"))[["gradient"]],
  d$evaluation_cache_hits %||% NA_integer_
))
