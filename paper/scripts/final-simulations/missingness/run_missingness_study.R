#!/usr/bin/env Rscript

# Simulation comparison for gamlss.longitudinal vs standard gamlss.
# Defaults are deliberately small so the script can be smoke-tested quickly.

safe_source <- function(path) {
  txt <- readLines(path, warn = FALSE, encoding = "UTF-8")
  if (length(txt) > 0) {
    bom <- intToUtf8(65279)
    txt[1] <- sub(paste0("^", bom), "", txt[1], useBytes = FALSE)
  }
  eval(parse(text = txt), envir = .GlobalEnv)
}

project_root <- normalizePath(getwd(), winslash = "/", mustWork = TRUE)
if (!requireNamespace("pkgload", quietly = TRUE)) stop("pkgload is required for verified checkout execution.")
pkgload::load_all(project_root, quiet = TRUE, export_all = TRUE, helpers = FALSE)
suppressPackageStartupMessages(library(gamlss.longitudinal))
list2env(as.list(getNamespace("gamlss.longitudinal"), all.names = TRUE), envir = .GlobalEnv)

suppressPackageStartupMessages({
  library(parallel)
  library(VineCopula)
  library(gamlss)
  library(gamlss2)
  library(gamlss.dist)
  library(ggplot2)
})

set.seed(20260513)

out_dir <- Sys.getenv("OUT_DIR", unset = file.path("results", "bcpe_t_missingness_comparison"))
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)
safe_source(file.path(project_root, "paper", "R", "missingness-study-helpers.R"))
expected_package_identity <- jss_missing_checkout_identity(project_root)
parent_package_identity <- jss_missing_verify_checkout(project_root, expected_package_identity, load_checkout = FALSE)
assign(".jss_missing_verified_worker_identity", parent_package_identity, envir = .GlobalEnv)
producer_sources <- c(
  file.path(project_root, "paper", "R", "missingness-study-helpers.R"),
  file.path(project_root, "paper", "scripts", "final-simulations", "missingness", "run_missingness_study.R")
)
producer_sha256 <- jss_missing_producer_sha256(producer_sources)
task_result_dir <- file.path(out_dir, "rep_results")
worker_log_dir <- file.path(out_dir, "worker_logs")
dir.create(task_result_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(worker_log_dir, recursive = TRUE, showWarnings = FALSE)

n_fits <- as.integer(Sys.getenv("N_FITS", unset = "20"))
rep_ids_env <- Sys.getenv("REP_IDS", unset = "")
rep_ids <- if (nzchar(rep_ids_env)) {
  as.integer(strsplit(rep_ids_env, ",", fixed = TRUE)[[1]])
} else {
  seq_len(n_fits)
}
rep_ids <- rep_ids[is.finite(rep_ids) & rep_ids >= 1L]
if (length(rep_ids) == 0L) {
  stop("REP_IDS must contain at least one positive integer replicate id.")
}
n_cores <- as.integer(Sys.getenv("N_CORES", unset = as.character(max(1, parallel::detectCores() - 2))))
smooth_k <- as.integer(Sys.getenv("SMOOTH_K", unset = "10"))
verbose_level <- as.integer(Sys.getenv("VERBOSE_FITS", unset = "0"))
verbose <- verbose_level > 0
parallel_setup_only <- as.logical(as.integer(Sys.getenv("PARALLEL_SETUP_ONLY", unset = "0")))
resume_checkpoints <- as.logical(as.integer(Sys.getenv("RESUME_CHECKPOINTS", unset = "1")))

missing_levels_env <- Sys.getenv("MISSING_LEVELS", unset = "0,0.1,0.2,0.3,0.4,0.5")
missing_levels <- as.numeric(strsplit(missing_levels_env, ",", fixed = TRUE)[[1]])
missing_levels <- missing_levels[is.finite(missing_levels) & missing_levels >= 0 & missing_levels < 1]
if (length(missing_levels) == 0L) {
  stop("MISSING_LEVELS must contain at least one value in [0, 1).", call. = FALSE)
}

missing_mechanisms_env <- Sys.getenv(
  "MISSING_MECHANISMS",
  unset = "monotone_dropout,time_dependent_intermittent_mar"
)
missing_mechanisms <- trimws(strsplit(missing_mechanisms_env, ",", fixed = TRUE)[[1]])
missing_mechanisms <- missing_mechanisms[nzchar(missing_mechanisms)]
missing_mechanisms <- jss_missing_validate_mechanisms(missing_mechanisms)
missingness_registry <- jss_missing_mechanism_registry()
missingness_registry <- missingness_registry[
  missingness_registry$missing_mechanism %in% missing_mechanisms, , drop = FALSE
]
write.csv(missingness_registry, file.path(out_dir, "missingness_design_registry.csv"), row.names = FALSE)
write.csv(jss_missing_estimand_registry(), file.path(out_dir, "missingness_estimand_registry.csv"), row.names = FALSE)

base_scenarios <- data.frame(
  n = 500L,
  d = 4L,
  stringsAsFactors = FALSE
)
scenarios <- merge(
  base_scenarios,
  expand.grid(
    missing_mechanism = missing_mechanisms,
    missing_rate = missing_levels,
    stringsAsFactors = FALSE
  ),
  by = NULL,
  all = TRUE
)

params_margin <- c("mu", "sigma", "nu", "tau")
params_copula <- c("theta", "zeta")
params_all <- c(params_margin, params_copula)
fixed_terms <- c("intercept", "x1", "x2", "t")
smooth_params_longitudinal <- c("mu", "sigma", "theta")
smooth_params_gamlss <- c("mu", "sigma")
compute_se <- as.logical(as.integer(Sys.getenv("COMPUTE_SE", unset = "1")))
save_fits <- as.logical(as.integer(Sys.getenv("SAVE_FITS", unset = "0")))
vcov_method_longitudinal <- Sys.getenv("VCOV_METHOD_LONGITUDINAL", unset = "analytical")
include_dlcopdpar <- as.logical(as.integer(Sys.getenv("INCLUDE_DLCOPDPAR", unset = "1")))
optim_method <- Sys.getenv("OPT_METHOD", unset = "RS")
max_outer_iter <- as.integer(Sys.getenv("MAX_OUTER_ITER", unset = "1000"))
max_inner_iter <- as.integer(Sys.getenv("MAX_INNER_ITER", unset = "100"))
max_elapsed_sec <- as.numeric(Sys.getenv("MAX_ELAPSED_SEC", unset = "180"))
outer_stop_crit_env <- Sys.getenv("OUTER_STOP_CRIT", unset = NA_character_)
outer_stop_crit <- if (is.na(outer_stop_crit_env) || !nzchar(outer_stop_crit_env)) {
  NA_real_
} else {
  as.numeric(outer_stop_crit_env)
}
inner_stop_crit_env <- Sys.getenv("INNER_STOP_CRIT", unset = NA_character_)
inner_stop_crit <- if (is.na(inner_stop_crit_env) || !nzchar(inner_stop_crit_env)) {
  NA_real_
} else {
  as.numeric(inner_stop_crit_env)
}
use_backtracking <- as.logical(as.integer(Sys.getenv("USE_BACKTRACKING", unset = "1")))
backtracking_max_halves <- as.integer(Sys.getenv("BACKTRACKING_MAX_HALVES", unset = "50"))
start_step_size <- as.numeric(Sys.getenv("START_STEP_SIZE", unset = "0.5"))
step_adjustment_env <- Sys.getenv("STEP_ADJUSTMENT", unset = NA_character_)
step_adjustment <- if (is.na(step_adjustment_env) || !nzchar(step_adjustment_env)) {
  NA_real_
} else {
  as.numeric(step_adjustment_env)
}
max_steps <- as.integer(Sys.getenv("MAX_STEPS", unset = "5"))
cg_max_delta <- as.numeric(Sys.getenv("CG_MAX_DELTA", unset = "0.5"))
cg_armijo_c1 <- as.numeric(Sys.getenv("CG_ARMIJO_C1", unset = "1e-4"))
cg_max_stall <- as.integer(Sys.getenv("CG_MAX_STALL", unset = "5"))
cg_update_lambda <- as.logical(as.integer(Sys.getenv("CG_UPDATE_LAMBDA", unset = "1")))
cg_line_search <- Sys.getenv("CG_LINE_SEARCH", unset = "best")
cg_max_line_search_evals_env <- Sys.getenv("CG_MAX_LINE_SEARCH_EVALS", unset = "60")
cg_max_line_search_evals <- if (is.na(cg_max_line_search_evals_env) || !nzchar(cg_max_line_search_evals_env)) {
  NA_integer_
} else {
  as.integer(cg_max_line_search_evals_env)
}
cg_gradient_method <- Sys.getenv("CG_GRADIENT_METHOD", unset = "forward")
cg_zeta_hessian <- Sys.getenv("CG_ZETA_HESSIAN", unset = "analytical")
cg_lambda_update_every <- as.integer(Sys.getenv("CG_LAMBDA_UPDATE_EVERY", unset = "10"))
cg_max_lambda_updates_env <- Sys.getenv("CG_MAX_LAMBDA_UPDATES", unset = NA_character_)
cg_max_lambda_updates <- if (is.na(cg_max_lambda_updates_env) || !nzchar(cg_max_lambda_updates_env)) {
  NA_integer_
} else {
  as.integer(cg_max_lambda_updates_env)
}
cg_raw_loglik_drop_tol_env <- Sys.getenv("CG_RAW_LOGLIK_DROP_TOL", unset = "10")
cg_raw_loglik_drop_tol <- if (is.na(cg_raw_loglik_drop_tol_env) || !nzchar(cg_raw_loglik_drop_tol_env)) {
  NA_real_
} else {
  as.numeric(cg_raw_loglik_drop_tol_env)
}
warm_start_joint <- as.logical(as.integer(Sys.getenv("WARM_START_JOINT", unset = "1")))
compute_predictive_scores <- as.logical(as.integer(Sys.getenv("COMPUTE_PREDICTIVE_SCORES", unset = "1")))
predictive_nsim <- as.integer(Sys.getenv("PREDICTIVE_NSIM", unset = "200"))
variogram_p <- as.numeric(Sys.getenv("VARIOGRAM_P", unset = "0.5"))
lambda_start_env <- Sys.getenv("LAMBDA_START", unset = NA_character_)
lambda_start <- if (is.na(lambda_start_env) || !nzchar(lambda_start_env)) {
  NA_real_
} else {
  as.numeric(lambda_start_env)
}
margin_family <- gamlss.dist::BCPE(mu.link = "log")
copula_link <- get_copula_dist("t")$copula_link

true_intercepts_natural <- list(
  mu = 2.20,
  sigma = 0.30,
  nu = 1.00,
  tau = margin_family$tau.linkinv(0.50),
  theta = copula_link$theta.linkinv(0.75),
  zeta = copula_link$zeta.linkinv(0.50)
)

true_beta <- list(
  mu = c(
    intercept = margin_family$mu.linkfun(true_intercepts_natural$mu),
    x1 = 0.12, x2 = -0.10, t = 0.18
  ),
  sigma = c(
    intercept = margin_family$sigma.linkfun(true_intercepts_natural$sigma),
    x1 = 0.18, x2 = 0.12, t = 0.15
  ),
  nu = c(intercept = 0.00, x1 = 0.10, x2 = -0.08, t = 0.12),
  tau = c(
    intercept = margin_family$tau.linkfun(true_intercepts_natural$tau),
    x1 = -0.08, x2 = 0.10, t = 0.11
  ),
  theta = c(
    intercept = copula_link$theta.linkfun(true_intercepts_natural$theta),
    x1 = 0.25, x2 = -0.10, t = 0.20
  ),
  zeta = c(
    intercept = copula_link$zeta.linkfun(true_intercepts_natural$zeta),
    x1 = 0.12, x2 = -0.10, t = 0.10
  )
)

smooth_truth <- list(
  mu = function(s1) 0.55 * sin(2 * pi * s1),
  sigma = function(s1) 0.30 * cos(2 * pi * s1) - 0.10 * (s1 - 0.50)^2,
  theta = function(s1) 2 * (0.25 * sin(pi * s1 + 0.40) + 0.12 * (s1 - 0.50))
)

checkpoint_configuration <- list(
  design_version = "missingness-v5",
  package_identity = expected_package_identity,
  producer_sha256 = producer_sha256,
  producer_sources = normalizePath(producer_sources, winslash = "/", mustWork = TRUE),
  estimand_registry = jss_missing_estimand_registry(),
  family = "BCPE_log_mu",
  copula = "t",
  true_beta = true_beta,
  smooth_truth_expression = lapply(smooth_truth, function(fn) paste(deparse(body(fn)), collapse = " ")),
  smooth_k = smooth_k,
  compute_se = compute_se,
  vcov_method_longitudinal = vcov_method_longitudinal,
  include_dlcopdpar = include_dlcopdpar,
  optim_method = optim_method,
  max_outer_iter = max_outer_iter,
  max_inner_iter = max_inner_iter,
  max_elapsed_sec = max_elapsed_sec,
  outer_stop_crit = outer_stop_crit,
  inner_stop_crit = inner_stop_crit,
  use_backtracking = use_backtracking,
  backtracking_max_halves = backtracking_max_halves,
  start_step_size = start_step_size,
  step_adjustment = step_adjustment,
  max_steps = max_steps,
  cg = list(
    max_delta = cg_max_delta, armijo_c1 = cg_armijo_c1, max_stall = cg_max_stall,
    update_lambda = cg_update_lambda, line_search = cg_line_search,
    max_line_search_evals = cg_max_line_search_evals, gradient_method = cg_gradient_method,
    zeta_hessian = cg_zeta_hessian, lambda_update_every = cg_lambda_update_every,
    max_lambda_updates = cg_max_lambda_updates, raw_loglik_drop_tol = cg_raw_loglik_drop_tol
  ),
  warm_start_joint = warm_start_joint,
  lambda_start = lambda_start,
  predictive = list(compute = compute_predictive_scores, nsim = predictive_nsim, variogram_p = variogram_p),
  mechanism_registry = jss_missing_mechanism_registry(),
  seed_formula_version = "simulation=100000+rep; mechanism offset + rate offset",
  rng_kind = RNGkind()
)
checkpoint_configuration_key <- paste0("missing-v5-", jss_missing_short_hash(checkpoint_configuration))

linpred <- function(beta, x1, x2, t) {
  beta[["intercept"]] + beta[["x1"]] * x1 + beta[["x2"]] * x2 + beta[["t"]] * t
}

clip_u <- function(u) {
  pmax(pmin(u, 1 - 1e-8), 1e-8)
}

safe_quantile <- function(x, p) {
  x <- x[is.finite(x)]
  if (length(x) == 0) {
    return(NA_real_)
  }
  as.numeric(stats::quantile(x, probs = p, names = FALSE, type = 7))
}

mc_mean_interval <- function(x, level = 0.95) {
  x <- as.numeric(x)
  x <- x[is.finite(x)]
  n <- length(x)
  estimate <- if (n) mean(x) else NA_real_
  mcse <- if (n > 1L) stats::sd(x) / sqrt(n) else NA_real_
  critical <- if (n > 1L) stats::qt(1 - (1 - level) / 2, df = n - 1L) else NA_real_
  c(n = n, estimate = estimate, mcse = mcse,
    conf_low = estimate - critical * mcse, conf_high = estimate + critical * mcse)
}

extract_fit_metric <- function(fit_obj, metric = c("logLik", "df")) {
  metric <- match.arg(metric)
  ll <- tryCatch(stats::logLik(fit_obj), error = function(e) NULL)
  if (is.null(ll)) {
    return(NA_real_)
  }

  pick_model_value <- function(x) {
    if (is.null(x) || length(x) == 0) {
      return(NA_real_)
    }
    x_num <- as.numeric(x)
    x_names <- names(x)
    if (!is.null(x_names) && "joint" %in% x_names) {
      return(x_num[match("joint", x_names)])
    }
    x_num[length(x_num)]
  }

  if (metric == "logLik") {
    return(pick_model_value(ll))
  }
  df_val <- attr(ll, "df")
  if (!is.null(df_val)) {
    return(pick_model_value(df_val))
  }
  if (!is.null(fit_obj$df)) {
    return(pick_model_value(fit_obj$df))
  }
  if (!is.null(fit_obj$edf)) {
    return(pick_model_value(fit_obj$edf))
  }
  if (!is.null(fit_obj$selection_criteria) && "EDF" %in% rownames(fit_obj$selection_criteria)) {
    return(pick_model_value(fit_obj$selection_criteria["EDF", ]))
  }
  if (!is.null(fit_obj$model_selection) && "EDF" %in% rownames(fit_obj$model_selection)) {
    return(pick_model_value(fit_obj$model_selection["EDF", ]))
  }
  if (!is.null(fit_obj$criteria) && "EDF" %in% rownames(fit_obj$criteria)) {
    return(pick_model_value(fit_obj$criteria["EDF", ]))
  }
  if (!is.null(fit_obj$df_s)) {
    df_s_vals <- suppressWarnings(as.numeric(unlist(fit_obj$df_s, use.names = FALSE)))
    df_s_vals <- df_s_vals[is.finite(df_s_vals)]
    if (length(df_s_vals) > 0 && !is.null(fit_obj$par)) {
      return(length(fit_obj$par) + sum(df_s_vals))
    }
  }
  if (!is.null(fit_obj$par)) {
    return(length(fit_obj$par))
  }
  NA_real_
}

extract_convergence_info <- function(fit_obj) {
  conv <- fit_obj$convergence
  if (is.null(conv)) {
    return(list(
      converged = FALSE,
      hit_outer_limit = NA,
      hit_max_stall = NA,
      hit_raw_loglik_deterioration = NA,
      stop_reason = "fit_error",
      grad_inf = NA_real_,
      step_l2 = NA_real_,
      best_raw_loglik = NA_real_,
      best_raw_loglik_iteration = NA_integer_,
      raw_loglik_drop_from_best = NA_real_,
      raw_loglik_drop_tol = NA_real_,
      outer_iterations = NA_integer_,
      outer_log_lik_change = NA_real_,
      outer_stop_crit = NA_real_
    ))
  }
  list(
    converged = if (!is.null(conv$converged)) isTRUE(conv$converged) else NA,
    hit_outer_limit = if (!is.null(conv$hit_outer_limit)) isTRUE(conv$hit_outer_limit) else NA,
    hit_max_stall = if (!is.null(conv$hit_max_stall)) isTRUE(conv$hit_max_stall) else NA,
    hit_raw_loglik_deterioration = if (!is.null(conv$hit_raw_loglik_deterioration)) isTRUE(conv$hit_raw_loglik_deterioration) else NA,
    stop_reason = if (!is.null(conv$stop_reason)) as.character(conv$stop_reason) else NA_character_,
    grad_inf = jss_missing_finite_scalar_or_na(conv$grad_inf),
    step_l2 = jss_missing_finite_scalar_or_na(conv$step_l2),
    best_raw_loglik = jss_missing_finite_scalar_or_na(conv$best_raw_loglik),
    best_raw_loglik_iteration = jss_missing_finite_scalar_or_na(
      conv$best_raw_loglik_iteration, integer = TRUE
    ),
    raw_loglik_drop_from_best = jss_missing_finite_scalar_or_na(conv$raw_loglik_drop_from_best),
    raw_loglik_drop_tol = jss_missing_finite_scalar_or_na(conv$raw_loglik_drop_tol),
    outer_iterations = jss_missing_finite_scalar_or_na(conv$outer_iterations, integer = TRUE),
    outer_log_lik_change = jss_missing_finite_scalar_or_na(conv$outer_log_lik_change),
    outer_stop_crit = jss_missing_finite_scalar_or_na(conv$outer_stop_crit)
  )
}

extract_gamlss2_convergence_info <- function(fit_obj) {
  if (is.null(fit_obj)) {
    return(list(converged = FALSE, stop_reason = "fit_error", iterations = NA_integer_))
  }
  iterations <- suppressWarnings(as.integer(fit_obj$iterations[[1L]]))
  maxit <- fit_obj$control$maxit
  max_outer <- if (is.null(maxit) || !length(maxit)) 20L else as.integer(maxit[[1L]])
  converged <- is.finite(iterations) && iterations < max_outer &&
    is.finite(extract_fit_metric(fit_obj, "logLik"))
  list(
    converged = isTRUE(converged),
    stop_reason = if (isTRUE(converged)) "relative_deviance_tolerance" else "outer_iteration_limit_or_invalid_loglik",
    iterations = iterations
  )
}

bind_non_null <- function(x) {
  keep <- x[!vapply(x, is.null, logical(1))]
  if (length(keep) == 0) {
    return(NULL)
  }
  do.call(rbind, keep)
}

task_result_path <- function(task) {
  file.path(
    task_result_dir,
    sprintf(
      "%s_scenario%d_%s_miss%02d_n%d_d%d_rep%d.rds",
      task$checkpoint_configuration_key,
      task$scenario_id,
      task$missing_mechanism,
      round(100 * task$missing_rate),
      task$n,
      task$d,
      task$rep
    )
  )
}

reverify_missingness_checkpoint_or_quarantine <- function(path) {
  tryCatch(
    jss_missing_reverify_sources(project_root, expected_package_identity,
      producer_sources, producer_sha256),
    error = function(e) {
      quarantine <- paste0(path, ".source-identity-changed-", Sys.getpid(), "-",
        format(Sys.time(), "%Y%m%d%H%M%S"))
      if (file.exists(path) && !file.rename(path, quarantine)) file.remove(path)
      stop(conditionMessage(e), "; checkpoint quarantined or removed.", call. = FALSE)
    }
  )
}

run_one_rep_and_save <- function(task) {
  final_path <- task_result_path(task)
  jss_missing_reverify_sources(project_root, expected_package_identity,
    producer_sources, producer_sha256)
  cat(sprintf(
    "[pid %s] starting scenario %d %s missing=%.0f%% n=%d d=%d rep=%d\n",
    Sys.getpid(), task$scenario_id, task$missing_mechanism, 100 * task$missing_rate,
    task$n, task$d, task$rep
  ))
  flush.console()
  result <- run_one_rep(task)
  jss_missing_reverify_sources(project_root, expected_package_identity,
    producer_sources, producer_sha256)
  if (!jss_missing_checkpoint_valid(result, task, checkpoint_configuration)) {
    stop("Refusing invalid in-memory missingness checkpoint payload.")
  }
  temporary_path <- paste0(final_path, ".", Sys.getpid(), ".tmp")
  jss_missing_reverify_sources(project_root, expected_package_identity,
    producer_sources, producer_sha256)
  saveRDS(result, temporary_path)
  tryCatch(jss_missing_reverify_sources(project_root, expected_package_identity,
    producer_sources, producer_sha256), error = function(e) {
      file.remove(temporary_path)
      stop(conditionMessage(e), "; temporary checkpoint removed.", call. = FALSE)
    })
  if (file.exists(final_path)) file.remove(final_path)
  if (!file.rename(temporary_path, final_path)) stop("Could not atomically install checkpoint: ", final_path)
  reverify_missingness_checkpoint_or_quarantine(final_path)
  durable <- tryCatch(readRDS(final_path), error = function(e) NULL)
  reverify_missingness_checkpoint_or_quarantine(final_path)
  if (!isTRUE(jss_missing_checkpoint_valid(durable, task, checkpoint_configuration)) ||
      !identical(durable$checkpoint_content_sha256, result$checkpoint_content_sha256) ||
      !identical(durable, result)) {
    quarantine <- paste0(final_path, ".invalid-", Sys.getpid(), "-", format(Sys.time(), "%Y%m%d%H%M%S"))
    if (!file.rename(final_path, quarantine)) file.remove(final_path)
    stop("Missingness checkpoint failed worker-side durable validation and was quarantined.")
  }
  cat(sprintf(
    "[pid %s] finished scenario %d %s missing=%.0f%% n=%d d=%d rep=%d -> %s\n",
    Sys.getpid(), task$scenario_id, task$missing_mechanism, 100 * task$missing_rate,
    task$n, task$d, task$rep, task_result_path(task)
  ))
  flush.console()
  result
}

center_curve <- function(y) {
  y - mean(y, na.rm = TRUE)
}

calc_smooth_mean <- function(data_used, parameter) {
  # Registered population target; deliberately independent of realized missingness.
  jss_missing_population_smooth_mean(smooth_truth[[parameter]])
}

fitted_population_smooth_mean_longitudinal <- function(fit_obj, data_used, parameter) {
  B_list <- fit_obj$model_matrix$s[[parameter]]
  par_s_list <- fit_obj$par_s[[parameter]]
  if (is.null(B_list) || !length(B_list) || is.null(par_s_list) || !length(par_s_list)) return(0)
  B <- B_list[[1L]]; beta <- par_s_list[[1L]]
  smooth_obs <- as.numeric(B %*% beta)
  x_obs <- attr(B, "smooth_x")
  if (is.null(x_obs)) {
    data_sub <- data_used
    if (parameter %in% params_copula) data_sub <- data_sub[data_sub$time < max(data_sub$time), , drop = FALSE]
    x_obs <- data_sub$s1
  }
  smooth_by_x <- tapply(smooth_obs, x_obs, mean)
  x_unique <- as.numeric(names(smooth_by_x)); y_unique <- as.numeric(smooth_by_x)
  if (length(x_unique) < 2L) return(NA_real_)
  mean(stats::approx(x_unique[order(x_unique)], y_unique[order(x_unique)],
    xout = seq(0, 1, length.out = 10001L), rule = 2)$y)
}

fitted_population_level_gamlss2 <- function(fit_obj, data_used, parameter) {
  grid <- seq(0, 1, length.out = 10001L)
  newdata <- data.frame(s1 = grid, x1 = 0, x2 = 0, t = 0,
    id = data_used$id[[1L]], time = data_used$time[[1L]],
    response = stats::median(data_used$response, na.rm = TRUE))
  pred <- tryCatch(stats::predict(fit_obj, newdata = newdata, type = "link"), error = function(e) NULL)
  if (is.null(pred)) return(NA_real_)
  mean(as.numeric(pred[, parameter]), na.rm = TRUE)
}

simulate_dataset <- function(n, d, seed) {
  set.seed(seed)

  subjects <- data.frame(
    id = seq_len(n),
    x1 = rnorm(n),
    x2 = rbinom(n, size = 1, prob = 0.5),
    s1 = runif(n, min = 0, max = 1),
    stringsAsFactors = FALSE
  )

  dat <- merge(subjects, data.frame(time_index = seq_len(d)), by = NULL, all = TRUE)
  dat <- dat[order(dat$id, dat$time_index), ]
  dat$t <- if (d > 1) (dat$time_index - 1) / (d - 1) else 0
  dat$time <- dat$time_index

  simulate_response_given_covariates(dat, seed = seed)
}

calibrate_missing_intercept <- function(linear_part, target_rate) {
  if (target_rate <= 0) {
    return(-Inf)
  }
  if (target_rate >= 1) {
    return(Inf)
  }
  stats::uniroot(
    function(a) mean(stats::plogis(a + linear_part), na.rm = TRUE) - target_rate,
    interval = c(-30, 30)
  )$root
}

apply_missingness <- function(dat, rate, mechanism, seed) {
  jss_missing_apply(dat, rate = rate, mechanism = mechanism, seed = seed)
}

missingness_summary <- function(dat, mechanism, target_rate) {
  jss_missing_summarize(dat, mechanism = mechanism, target_rate = target_rate)
}

simulate_response_given_covariates <- function(dat, seed) {
  set.seed(seed)

  dat <- dat[order(dat$id, dat$time), , drop = FALSE]
  if (!"time_index" %in% names(dat)) {
    dat$time_index <- as.integer(factor(dat$time, levels = sort(unique(dat$time))))
  }
  d <- length(unique(dat$time_index))
  n <- length(unique(dat$id))

  eta_mu <- linpred(true_beta$mu, dat$x1, dat$x2, dat$t) + smooth_truth$mu(dat$s1)
  eta_sigma <- linpred(true_beta$sigma, dat$x1, dat$x2, dat$t) + smooth_truth$sigma(dat$s1)
  eta_nu <- linpred(true_beta$nu, dat$x1, dat$x2, dat$t)
  eta_tau <- linpred(true_beta$tau, dat$x1, dat$x2, dat$t)
  eta_theta <- linpred(true_beta$theta, dat$x1, dat$x2, dat$t) + smooth_truth$theta(dat$s1)
  eta_zeta <- linpred(true_beta$zeta, dat$x1, dat$x2, dat$t)

  mu <- margin_family$mu.linkinv(eta_mu)
  sigma <- margin_family$sigma.linkinv(eta_sigma)
  nu <- margin_family$nu.linkinv(eta_nu)
  tau <- margin_family$tau.linkinv(eta_tau)
  theta <- copula_link$theta.linkinv(eta_theta)
  zeta <- copula_link$zeta.linkinv(eta_zeta)

  fam_t <- as.numeric(VineCopula::BiCopName("t"))
  U <- matrix(NA_real_, nrow = n, ncol = d)

  for (i in seq_len(n)) {
    u_prev <- runif(1)
    U[i, 1] <- u_prev

    if (d > 1) {
      for (j in 2:d) {
        row_prev <- which(dat$id == i & dat$time_index == (j - 1))
        u_new <- tryCatch(
          VineCopula::BiCopCondSim(
            N = 1,
            cond.val = u_prev,
            cond.var = 1,
            family = fam_t,
            par = theta[row_prev],
            par2 = zeta[row_prev]
          ),
          error = function(e) runif(1)
        )
        u_prev <- clip_u(as.numeric(u_new))
        U[i, j] <- u_prev
      }
    }
  }

  dat$u <- as.vector(t(U))
  dat$response <- gamlss.dist::qBCPE(dat$u, mu = mu, sigma = sigma, nu = nu, tau = tau)

  dat[, c("id", "time", "response", "s1", "x1", "x2", "t")]
}

fit_longitudinal_model <- function(dat) {
  miss <- attr(dat, "missingness_summary")
  segmented <- !is.null(miss) && miss$n_subjects_with_interior_gaps[[1L]] > 0L
  gamlss.longitudinal(
    dataset = dat,
    margin_dist = gamlss.dist::BCPE(mu.link = "log"),
    copula_dist = "t",
    time_var = "time",
    subject_var = "id",
    missingness = if (segmented) "segment" else "error",
    mu.formula = sprintf("response ~ x1 + x2 + t + s(s1, k = %d)", smooth_k),
    sigma.formula = sprintf("~ x1 + x2 + t + s(s1, k = %d)", smooth_k),
    nu.formula = "~ x1 + x2 + t",
    tau.formula = "~ x1 + x2 + t",
    theta.formula = sprintf("~ x1 + x2 + t + s(s1, k = %d)", smooth_k),
    zeta.formula = "~ x1 + x2 + t",
    verbose = verbose_level,
    compute_vcov = isTRUE(compute_se) && !segmented,
    vcov_method = vcov_method_longitudinal,
    include_dlcopdpar = include_dlcopdpar,
    inner_stop_crit = inner_stop_crit,
    outer_stop_crit = outer_stop_crit,
    max_outer_iter = max_outer_iter,
    max_inner_iter = max_inner_iter,
    max_elapsed_sec = max_elapsed_sec,
    start_step_size = start_step_size,
    step_adjustment = step_adjustment,
    max_steps = max_steps,
    use_backtracking = use_backtracking,
    backtracking_max_halves = backtracking_max_halves,
    method = optim_method,
    cg_max_delta = cg_max_delta,
    cg_armijo_c1 = cg_armijo_c1,
    cg_max_stall = cg_max_stall,
    cg_update_lambda = cg_update_lambda,
    cg_lambda_update_every = cg_lambda_update_every,
    cg_max_lambda_updates = cg_max_lambda_updates,
    cg_raw_loglik_drop_tol = cg_raw_loglik_drop_tol,
    cg_line_search = cg_line_search,
    cg_max_line_search_evals = cg_max_line_search_evals,
    cg_gradient_method = cg_gradient_method,
    cg_zeta_hessian = cg_zeta_hessian,
    warm_start_joint = warm_start_joint,
    lambda_start = lambda_start
  )
}

fit_gamlss2_model <- function(dat) {
  gamlss2::gamlss2(
    formula = stats::as.formula(sprintf(
      "response ~ x1 + x2 + t + s(s1, k = %d) | x1 + x2 + t + s(s1, k = %d) | x1 + x2 + t | x1 + x2 + t",
      smooth_k,
      smooth_k
    )),
    family = gamlss.dist::BCPE(mu.link = "log"),
    data = dat,
    trace = isTRUE(verbose)
  )
}

extract_one_longitudinal_term <- function(par_vec, parameter, term_name) {
  nm <- names(par_vec)
  pref <- paste0("^", parameter, "\\.")
  idx <- grep(pref, nm)
  if (length(idx) == 0) {
    return(NA_real_)
  }

  sub_names <- sub(pref, "", nm[idx])
  values <- as.numeric(par_vec[idx])
  target_pattern <- switch(
    term_name,
    intercept = "^(intercept|\\(Intercept\\))$",
    x1 = "^x1$",
    x2 = "^x2$",
    t = "^t$"
  )
  hit <- grep(target_pattern, sub_names)
  if (length(hit) == 0) {
    return(NA_real_)
  }
  values[hit[1]]
}

extract_one_gamlss_term <- function(coef_vec, term_name) {
  if (is.null(coef_vec) || length(coef_vec) == 0) {
    return(NA_real_)
  }
  nm <- names(coef_vec)
  target_pattern <- switch(
    term_name,
    intercept = "^(intercept|\\(Intercept\\))$",
    x1 = "^x1$",
    x2 = "^x2$",
    t = "^t$"
  )
  hit <- grep(target_pattern, nm)
  if (length(hit) == 0) {
    return(NA_real_)
  }
  as.numeric(coef_vec[hit[1]])
}

extract_fixed_estimates_longitudinal <- function(fit_obj, data_used) {
  rows <- list()
  par_vec <- fit_obj$par
  se_vec <- rep(NA_real_, length(par_vec))
  names(se_vec) <- names(par_vec)

  if (isTRUE(compute_se)) {
    sum_obj <- tryCatch(
      summary(fit_obj, include_vcov = TRUE, vcov_method = vcov_method_longitudinal),
      error = function(e) NULL
    )
    if (!is.null(sum_obj) && !is.null(sum_obj$coefficients)) {
      coef_tab <- sum_obj$coefficients
      if (all(c("term", "std_error") %in% names(coef_tab))) {
        se_tmp <- coef_tab$std_error
        names(se_tmp) <- coef_tab$term
        se_vec[names(se_vec)] <- se_tmp[names(se_vec)]
      }
    }
  }

  for (p in params_all) {
    smooth_mean <- if (p %in% names(smooth_truth)) calc_smooth_mean(data_used, p) else 0
    fitted_smooth_mean <- 0
    if (
      p %in% names(fit_obj$model_matrix$s) &&
        p %in% names(fit_obj$par_s) &&
        length(fit_obj$model_matrix$s[[p]]) > 0 &&
        length(fit_obj$par_s[[p]]) > 0
    ) {
      fitted_smooth_mean <- fitted_population_smooth_mean_longitudinal(fit_obj, data_used, p)
    }
    for (tm in fixed_terms) {
      true_value <- true_beta[[p]][[tm]]
      estimate <- extract_one_longitudinal_term(par_vec, p, tm)
      if (tm == "intercept") {
        true_value <- true_value + smooth_mean
        estimate <- estimate + fitted_smooth_mean
      }
      population_intercept <- tm == "intercept" && p %in% names(smooth_truth)
      rows[[length(rows) + 1]] <- data.frame(
        model = "gamlss.longitudinal",
        parameter = p,
        term = tm,
        estimate = estimate,
        # The population intercept is a contrast of the fixed intercept and
        # the fitted population smooth mean.  A coefficient-only SE omits
        # their covariance and is therefore not a valid SE for this estimand.
        std_error = if (population_intercept) NA_real_ else
          extract_one_longitudinal_term(se_vec, p, tm),
        true_value = true_value,
        intercept_includes_fitted_smooth_mean = population_intercept,
        inference_status = if (population_intercept) {
          "not_available_without_joint_fixed_smooth_covariance"
        } else {
          "coefficient_covariance_available_when_finite"
        },
        stringsAsFactors = FALSE
      )
    }
  }

  do.call(rbind, rows)
}

extract_fixed_estimates_gamlss2 <- function(fit_obj, data_used) {
  rows <- list()
  coef_by_param <- fit_obj$coefficients
  se_by_param <- setNames(vector("list", length(params_margin)), params_margin)
  for (p in params_margin) {
    coef_names <- names(coef_by_param[[p]])
    se_by_param[[p]] <- rep(NA_real_, length(coef_names))
    names(se_by_param[[p]]) <- coef_names
  }

  if (isTRUE(compute_se)) {
    se_all <- tryCatch(
      stats::vcov(fit_obj, type = "se", full = TRUE),
      error = function(e) NULL
    )
    if (is.null(se_all)) {
      V <- tryCatch(stats::vcov(fit_obj, full = TRUE), error = function(e) NULL)
      if (!is.null(V)) {
        se_all <- sqrt(pmax(0, diag(V)))
      }
    }
    if (!is.null(se_all)) {
      for (p in params_margin) {
        coef_names <- names(coef_by_param[[p]])
        for (tm in coef_names) {
          candidates <- c(
            paste0(p, ".p.", tm),
            paste0(p, ".", tm),
            paste0(p, ":", tm),
            tm
          )
          hit <- candidates[candidates %in% names(se_all)]
          if (length(hit) > 0) {
            se_by_param[[p]][tm] <- se_all[hit[1]]
          }
        }
      }
    }
  }

  for (p in params_all) {
    smooth_mean <- if (p %in% smooth_params_gamlss) calc_smooth_mean(data_used, p) else 0
    for (tm in fixed_terms) {
      true_value <- true_beta[[p]][[tm]]
      if (tm == "intercept") {
        true_value <- true_value + smooth_mean
      }
      estimate <- if (p %in% params_margin) {
        extract_one_gamlss_term(coef_by_param[[p]], tm)
      } else {
        NA_real_
      }
      if (tm == "intercept" && p %in% smooth_params_gamlss) {
        estimate <- fitted_population_level_gamlss2(fit_obj, data_used, p)
      }
      population_intercept <- tm == "intercept" && p %in% smooth_params_gamlss
      std_error <- if (p %in% params_margin) {
        extract_one_gamlss_term(se_by_param[[p]], tm)
      } else {
        NA_real_
      }
      if (population_intercept) std_error <- NA_real_
      rows[[length(rows) + 1]] <- data.frame(
        model = "gamlss2",
        parameter = p,
        term = tm,
        estimate = estimate,
        std_error = std_error,
        true_value = if (p %in% params_margin) true_value else NA_real_,
        intercept_includes_fitted_smooth_mean = tm == "intercept" && p %in% smooth_params_gamlss,
        inference_status = if (population_intercept) {
          "not_available_without_joint_fixed_smooth_covariance"
        } else {
          "coefficient_covariance_available_when_finite"
        },
        stringsAsFactors = FALSE
      )
    }
  }

  do.call(rbind, rows)
}

uniform_ks_stat <- function(x) {
  x <- x[is.finite(x)]
  if (length(x) < 2) return(NA_real_)
  x <- sort(pmax(pmin(x, 1), 0))
  n <- length(x)
  max(max(abs(x - seq_len(n) / n)), max(abs(x - (seq_len(n) - 1) / n)))
}

uniform_cvm_stat <- function(x) {
  x <- x[is.finite(x)]
  if (length(x) < 2) return(NA_real_)
  x <- sort(pmax(pmin(x, 1), 0))
  n <- length(x)
  1 / (12 * n) + sum((x - (2 * seq_len(n) - 1) / (2 * n))^2)
}

adjacent_pair_rows <- function(dat) {
  dat <- dat[order(dat$id, dat$time), ]
  split_dat <- split(seq_len(nrow(dat)), dat$id)
  out <- lapply(split_dat, function(idx) {
    if (length(idx) < 2) return(NULL)
    data.frame(row1 = idx[-length(idx)], row2 = idx[-1])
  })
  out <- out[!vapply(out, is.null, logical(1))]
  if (length(out) == 0) return(data.frame(row1 = integer(0), row2 = integer(0)))
  do.call(rbind, out)
}

rosenblatt_lag1_cor <- function(r, dat, normal_scores = FALSE) {
  dat_tmp <- dat[order(dat$id, dat$time), , drop = FALSE]
  r_tmp <- r[order(dat$id, dat$time)]
  if (normal_scores) {
    r_tmp <- stats::qnorm(pmax(pmin(r_tmp, 1 - 1e-8), 1e-8))
  }
  pairs <- adjacent_pair_rows(dat_tmp)
  suppressWarnings(stats::cor(r_tmp[pairs$row1], r_tmp[pairs$row2], use = "complete.obs"))
}

rosenblatt_subject_mean <- function(r, dat, normal_scores = FALSE) {
  dat_tmp <- dat[order(dat$id, dat$time), , drop = FALSE]
  r_tmp <- r[order(dat$id, dat$time)]
  if (normal_scores) {
    r_tmp <- stats::qnorm(pmax(pmin(r_tmp, 1 - 1e-8), 1e-8))
    center <- 0
  } else {
    center <- 0.5
  }
  mean_by_time <- stats::aggregate(r_tmp, by = list(time = dat_tmp$time), FUN = mean, na.rm = TRUE)
  names(mean_by_time)[2] <- "mean_residual"
  mean(abs(mean_by_time$mean_residual - center), na.rm = TRUE)
}

joint_metrics_from_rosenblatt <- function(model_name, scenario_label, n_val, d_val, rep_id, r, dat, logLik, df) {
  lag1_cor <- rosenblatt_lag1_cor(r, dat, normal_scores = FALSE)
  lag1_z_cor <- rosenblatt_lag1_cor(r, dat, normal_scores = TRUE)
  data.frame(
    scenario = scenario_label,
    n = n_val,
    d = d_val,
    rep = rep_id,
    model = model_name,
    logLik = logLik,
    df = df,
    rosenblatt_ks = uniform_ks_stat(r),
    rosenblatt_cvm = uniform_cvm_stat(r),
    rosenblatt_lag1_cor = lag1_cor,
    abs_rosenblatt_lag1_cor = abs(lag1_cor),
    rosenblatt_normal_lag1_cor = lag1_z_cor,
    abs_rosenblatt_normal_lag1_cor = abs(lag1_z_cor),
    rosenblatt_mean_abs_time_mean = rosenblatt_subject_mean(r, dat, normal_scores = FALSE),
    rosenblatt_normal_mean_abs_time_mean = rosenblatt_subject_mean(r, dat, normal_scores = TRUE),
    stringsAsFactors = FALSE
  )
}

joint_metrics_longitudinal <- function(fit_obj, dat, scenario_label, n_val, d_val, rep_id) {
  copula_link_fit <- get_copula_dist(fit_obj$copula_dist)$copula_link
  eta_out <- calc_eta(fit_obj$par, fit_obj$model_matrix, fit_obj$margin_dist, copula_link_fit, fit_obj$par_s)
  eta_inv <- eta_out$eta_inv
  pfun <- get(paste0("p", fit_obj$margin_dist$family[1]), envir = asNamespace("gamlss.dist"))
  u <- pfun(
    fit_obj$response,
    mu = eta_inv$mu,
    sigma = eta_inv$sigma,
    nu = eta_inv$nu,
    tau = eta_inv$tau
  )
  u <- pmax(pmin(as.numeric(u), 1 - 1e-8), 1e-8)

  pairs <- adjacent_pair_rows(data.frame(id = fit_obj$response_subject, time = fit_obj$response_margin))
  pair_u1 <- u[pairs$row1]
  pair_u2 <- u[pairs$row2]

  theta_rows <- which(fit_obj$response_margin < max(fit_obj$response_margin, na.rm = TRUE))
  theta_lookup <- rep(NA_integer_, length(fit_obj$response))
  theta_lookup[theta_rows] <- seq_along(theta_rows)
  theta_idx <- theta_lookup[pairs$row1]
  theta <- as.numeric(eta_inv$theta[theta_idx])
  zeta <- if (!is.null(eta_inv$zeta)) as.numeric(eta_inv$zeta[theta_idx]) else NA_real_
  fam_num <- as.numeric(VineCopula::BiCopName(fit_obj$copula_dist))
  conditional_rosenblatt <- tryCatch(
    # VineCopula convention: Hfunc1(u1, u2) is F(U2 <= u2 | U1 = u1).
    # For adjacent longitudinal pairs this gives the second residual conditional
    # on the previous time point.
    VineCopula::BiCopHfunc1(pair_u1, pair_u2, family = fam_num, par = theta, par2 = zeta),
    error = function(e) rep(NA_real_, length(pair_u1))
  )
  r <- u
  r[pairs$row2] <- pmax(pmin(as.numeric(conditional_rosenblatt), 1 - 1e-8), 1e-8)

  joint_metrics_from_rosenblatt(
    "gamlss.longitudinal", scenario_label, n_val, d_val, rep_id,
    r,
    data.frame(id = fit_obj$response_subject, time = fit_obj$response_margin),
    extract_fit_metric(fit_obj, "logLik"),
    extract_fit_metric(fit_obj, "df")
  )
}

joint_metrics_gamlss2 <- function(fit_obj, dat, scenario_label, n_val, d_val, rep_id) {
  pred <- stats::predict(fit_obj, newdata = dat, type = "link")
  u <- gamlss.dist::pBCPE(
    dat$response,
    mu = margin_family$mu.linkinv(pred[, "mu"]),
    sigma = margin_family$sigma.linkinv(pred[, "sigma"]),
    nu = margin_family$nu.linkinv(pred[, "nu"]),
    tau = margin_family$tau.linkinv(pred[, "tau"])
  )
  u <- pmax(pmin(as.numeric(u), 1 - 1e-8), 1e-8)

  # GAMLSS has no fitted serial dependence model. Its fitted joint model is
  # therefore the independence copula, for which F(U_t | U_{t-1}) = U_t.
  # Writing this explicitly keeps the diagnostic interpretation aligned with
  # the copula models, where later time points use conditional Rosenblatt
  # residuals from the fitted copula.
  pairs <- adjacent_pair_rows(dat)
  r <- u
  r[pairs$row2] <- u[pairs$row2]

  joint_metrics_from_rosenblatt(
    "gamlss2", scenario_label, n_val, d_val, rep_id,
    r,
    dat,
    extract_fit_metric(fit_obj, "logLik"),
    extract_fit_metric(fit_obj, "df")
  )
}

response_matrix <- function(dat, value_col = "response") {
  dat <- dat[order(dat$id, dat$time), , drop = FALSE]
  ids <- sort(unique(dat$id))
  times <- sort(unique(dat$time))
  out <- matrix(NA_real_, nrow = length(ids), ncol = length(times))
  rownames(out) <- ids
  colnames(out) <- times
  row_idx <- match(dat$id, ids)
  col_idx <- match(dat$time, times)
  out[cbind(row_idx, col_idx)] <- dat[[value_col]]
  out
}

variogram_score <- function(observed_mat, simulated_array, p = 0.5) {
  if (length(dim(simulated_array)) != 3L) {
    return(NA_real_)
  }
  n <- nrow(observed_mat)
  d <- ncol(observed_mat)
  if (d < 2 || dim(simulated_array)[2] != n || dim(simulated_array)[3] != d) {
    return(NA_real_)
  }
  pairs <- utils::combn(seq_len(d), 2)
  scores <- numeric(n)
  for (i in seq_len(n)) {
    y <- observed_mat[i, ]
    if (any(!is.finite(y))) {
      scores[i] <- NA_real_
      next
    }
    score_i <- 0
    for (k in seq_len(ncol(pairs))) {
      a <- pairs[1, k]
      b <- pairs[2, k]
      obs_diff <- abs(y[a] - y[b])^p
      sim_diff <- abs(simulated_array[, i, a] - simulated_array[, i, b])^p
      score_i <- score_i + (obs_diff - mean(sim_diff, na.rm = TRUE))^2
    }
    scores[i] <- score_i
  }
  mean(scores, na.rm = TRUE)
}

simulate_predictive_longitudinal <- function(fit_obj, nsim) {
  copula_link_fit <- get_copula_dist(fit_obj$copula_dist)$copula_link
  eta_out <- calc_eta(fit_obj$par, fit_obj$model_matrix, fit_obj$margin_dist, copula_link_fit, fit_obj$par_s)
  eta_inv <- eta_out$eta_inv
  ids <- sort(unique(fit_obj$response_subject))
  times <- sort(unique(fit_obj$response_margin))
  n <- length(ids)
  d <- length(times)
  row_map <- matrix(NA_integer_, nrow = n, ncol = d)
  for (i in seq_len(n)) {
    for (j in seq_len(d)) {
      row_map[i, j] <- which(fit_obj$response_subject == ids[i] & fit_obj$response_margin == times[j])[1]
    }
  }

  theta_rows <- which(fit_obj$response_margin < max(fit_obj$response_margin, na.rm = TRUE))
  theta_lookup <- rep(NA_integer_, length(fit_obj$response))
  theta_lookup[theta_rows] <- seq_along(theta_rows)
  fam_num <- as.numeric(VineCopula::BiCopName(fit_obj$copula_dist))
  out <- array(NA_real_, dim = c(nsim, n, d))

  for (s in seq_len(nsim)) {
    U <- matrix(NA_real_, nrow = n, ncol = d)
    for (i in seq_len(n)) {
      u_prev <- stats::runif(1)
      U[i, 1] <- u_prev
      if (d > 1) {
        for (j in 2:d) {
          prev_row <- row_map[i, j - 1]
          theta_idx <- theta_lookup[prev_row]
          theta <- as.numeric(eta_inv$theta[theta_idx])
          zeta <- if (!is.null(eta_inv$zeta)) as.numeric(eta_inv$zeta[theta_idx]) else NA_real_
          u_new <- tryCatch(
            VineCopula::BiCopCondSim(
              N = 1,
              cond.val = u_prev,
              cond.var = 1,
              family = fam_num,
              par = theta,
              par2 = zeta
            ),
            error = function(e) stats::runif(1)
          )
          u_prev <- clip_u(as.numeric(u_new))
          U[i, j] <- u_prev
        }
      }
    }
    u_vec <- as.vector(t(U))
    out[s, , ] <- matrix(
      gamlss.dist::qBCPE(
        u_vec,
        mu = eta_inv$mu,
        sigma = eta_inv$sigma,
        nu = eta_inv$nu,
        tau = eta_inv$tau
      ),
      nrow = n,
      ncol = d,
      byrow = TRUE
    )
  }
  out
}

simulate_predictive_gamlss2 <- function(fit_obj, dat, nsim) {
  dat <- dat[order(dat$id, dat$time), , drop = FALSE]
  pred <- stats::predict(fit_obj, newdata = dat, type = "link")
  n <- length(unique(dat$id))
  d <- length(unique(dat$time))
  out <- array(NA_real_, dim = c(nsim, n, d))
  for (s in seq_len(nsim)) {
    y <- gamlss.dist::qBCPE(
      stats::runif(nrow(dat)),
      mu = margin_family$mu.linkinv(pred[, "mu"]),
      sigma = margin_family$sigma.linkinv(pred[, "sigma"]),
      nu = margin_family$nu.linkinv(pred[, "nu"]),
      tau = margin_family$tau.linkinv(pred[, "tau"])
    )
    out[s, , ] <- matrix(y, nrow = n, ncol = d, byrow = TRUE)
  }
  out
}

predictive_scores_longitudinal <- function(fit_obj, test_dat, scenario_label, n_val, d_val, rep_id, nsim = predictive_nsim) {
  copula_link_fit <- get_copula_dist(fit_obj$copula_dist)$copula_link
  eta_out <- calc_eta(fit_obj$par, fit_obj$model_matrix, fit_obj$margin_dist, copula_link_fit, fit_obj$par_s)
  lik <- calc_likelihood_minimal(
    eta_out$eta_inv,
    mm = fit_obj$model_matrix$x,
    margin_dist = fit_obj$margin_dist,
    copula_dist = fit_obj$copula_dist,
    calc_d2 = FALSE,
    response = test_dat$response,
    response_margin = test_dat$time,
    response_subject = test_dat$id
  )
  sim_array <- simulate_predictive_longitudinal(fit_obj, nsim = nsim)
  y_obs <- response_matrix(test_dat)
  data.frame(
    scenario = scenario_label,
    n = n_val,
    d = d_val,
    rep = rep_id,
    model = "gamlss.longitudinal",
    test_log_score_joint = as.numeric(lik$log_lik["joint"]),
    test_log_score_marginal = as.numeric(lik$log_lik["marginal"]),
    test_log_score_copula = as.numeric(lik$log_lik["copula"]),
    test_log_score_per_obs = as.numeric(lik$log_lik["joint"]) / nrow(test_dat),
    variogram_score = variogram_score(y_obs, sim_array, p = variogram_p),
    variogram_p = variogram_p,
    predictive_nsim = nsim,
    stringsAsFactors = FALSE
  )
}

predictive_scores_gamlss2 <- function(fit_obj, test_dat, scenario_label, n_val, d_val, rep_id, nsim = predictive_nsim) {
  test_dat <- test_dat[order(test_dat$id, test_dat$time), , drop = FALSE]
  pred <- stats::predict(fit_obj, newdata = test_dat, type = "link")
  log_d <- gamlss.dist::dBCPE(
    test_dat$response,
    mu = margin_family$mu.linkinv(pred[, "mu"]),
    sigma = margin_family$sigma.linkinv(pred[, "sigma"]),
    nu = margin_family$nu.linkinv(pred[, "nu"]),
    tau = margin_family$tau.linkinv(pred[, "tau"]),
    log = TRUE
  )
  joint_log_score <- sum(log_d, na.rm = TRUE)
  sim_array <- simulate_predictive_gamlss2(fit_obj, test_dat, nsim = nsim)
  y_obs <- response_matrix(test_dat)
  data.frame(
    scenario = scenario_label,
    n = n_val,
    d = d_val,
    rep = rep_id,
    model = "gamlss2",
    test_log_score_joint = joint_log_score,
    test_log_score_marginal = joint_log_score,
    test_log_score_copula = 0,
    test_log_score_per_obs = joint_log_score / nrow(test_dat),
    variogram_score = variogram_score(y_obs, sim_array, p = variogram_p),
    variogram_p = variogram_p,
    predictive_nsim = nsim,
    stringsAsFactors = FALSE
  )
}

smooth_curve_longitudinal <- function(fit_obj, data_used, parameter, s1_grid) {
  B_list <- fit_obj$model_matrix$s[[parameter]]
  par_s_list <- fit_obj$par_s[[parameter]]

  if (is.null(B_list) || length(B_list) == 0 || is.null(par_s_list) || length(par_s_list) == 0) {
    return(rep(NA_real_, length(s1_grid)))
  }

  s_name <- names(B_list)[1]
  B <- B_list[[s_name]]
  b_s <- par_s_list[[s_name]]
  smooth_obs <- as.numeric(B %*% b_s)
  x_obs <- attr(B, "smooth_x")

  if (is.null(x_obs)) {
    data_sub <- data_used
    if (parameter %in% params_copula) {
      data_sub <- data_sub[data_sub$time < max(data_sub$time), , drop = FALSE]
    }
    x_obs <- data_sub$s1
  }

  smooth_by_x <- tapply(smooth_obs, x_obs, mean)
  x_unique <- as.numeric(names(smooth_by_x))
  y_unique <- as.numeric(smooth_by_x)
  ord <- order(x_unique)

  if (length(x_unique) < 2) {
    return(rep(NA_real_, length(s1_grid)))
  }

  center_curve(stats::approx(x_unique[ord], y_unique[ord], xout = s1_grid, rule = 2)$y)
}

smooth_curve_gamlss2 <- function(fit_obj, data_used, parameter, s1_grid) {
  newdata <- data.frame(
    s1 = s1_grid,
    x1 = 0,
    x2 = 0,
    t = 0,
    id = data_used$id[1],
    time = data_used$time[1],
    response = stats::median(data_used$response, na.rm = TRUE)
  )

  eta_grid <- tryCatch(
    {
      pred <- stats::predict(fit_obj, newdata = newdata, type = "link")
      as.numeric(pred[, parameter])
    },
    error = function(e) rep(NA_real_, length(s1_grid))
  )
  center_curve(eta_grid)
}

run_one_rep <- function(task) {
  n_val <- as.integer(task$n)
  d_val <- as.integer(task$d)
  rep_id <- as.integer(task$rep)
  missing_mechanism <- task$missing_mechanism
  missing_rate <- task$missing_rate
  scenario_label <- sprintf(
    "n%d_d%d_%s_miss%02d",
    n_val,
    d_val,
    missing_mechanism,
    round(100 * missing_rate)
  )
  seed <- 100000L + rep_id
  missingness_seed <- seed +
    jss_missing_seed_offset(missing_mechanism) +
    as.integer(round(1000 * missing_rate))
  s1_grid <- seq(0, 1, length.out = 101)

  dat <- simulate_dataset(n = n_val, d = d_val, seed = seed)
  dat <- apply_missingness(
    dat,
    rate = missing_rate,
    mechanism = missing_mechanism,
    seed = missingness_seed
  )
  miss_summary <- attr(dat, "missingness_summary")
  miss_summary$scenario <- scenario_label
  miss_summary$n <- n_val
  miss_summary$d <- d_val
  miss_summary$rep <- rep_id
  missingness_pattern <- data.frame(
    scenario = scenario_label,
    n = n_val,
    d = d_val,
    rep = rep_id,
    missing_mechanism = missing_mechanism,
    id = as.integer(dat$id),
    time = as.numeric(dat$time),
    response_observed = as.logical(dat$response_observed),
    stringsAsFactors = FALSE
  )

  run_rows <- list()
  fixed_rows <- list()
  smooth_rows <- list()
  joint_rows <- list()
  predictive_rows <- list()

  t0 <- Sys.time()
  err_long <- NA_character_
  captured_long <- jss_missing_capture_warnings(
    tryCatch(
      fit_longitudinal_model(dat),
      error = function(e) {
        err_long <<- conditionMessage(e)
        NULL
      }
    )
  )
  fit_long <- captured_long$value
  elapsed_long <- as.numeric(difftime(Sys.time(), t0, units = "secs"))

  t0 <- Sys.time()
  err_gamlss2 <- NA_character_
  captured_gamlss2 <- jss_missing_capture_warnings(
    tryCatch(
      fit_gamlss2_model(dat),
      error = function(e) {
        err_gamlss2 <<- conditionMessage(e)
        NULL
      }
    )
  )
  fit_gamlss2 <- captured_gamlss2$value
  elapsed_gamlss2 <- as.numeric(difftime(Sys.time(), t0, units = "secs"))

  conv_long <- if (!is.null(fit_long)) extract_convergence_info(fit_long) else extract_convergence_info(NULL)
  conv_gamlss2 <- extract_gamlss2_convergence_info(fit_gamlss2)
  if (is.null(fit_long)) conv_long$stop_reason <- jss_missing_failure_stop_reason(err_long)
  if (is.null(fit_gamlss2)) conv_gamlss2$stop_reason <- jss_missing_failure_stop_reason(err_gamlss2)
  retained_long <- !is.null(fit_long) && isTRUE(conv_long$converged)
  retained_gamlss2 <- !is.null(fit_gamlss2) && isTRUE(conv_gamlss2$converged)
  warning_events <- rbind(
    jss_missing_warning_events(captured_long, scenario_label, n_val, d_val, rep_id,
      missing_mechanism, missing_rate, "gamlss.longitudinal"),
    jss_missing_warning_events(captured_gamlss2, scenario_label, n_val, d_val, rep_id,
      missing_mechanism, missing_rate, "gamlss2")
  )

  run_rows[[1]] <- data.frame(
    scenario = scenario_label, n = n_val, d = d_val, rep = rep_id,
    simulation_seed = seed, missingness_seed = missingness_seed,
    missing_mechanism = missing_mechanism,
    missingness_label = miss_summary$missingness_label,
    missingness_pattern = miss_summary$missingness_pattern,
    analysis_role = miss_summary$analysis_role,
    target_missing_rate = missing_rate,
    observed_missing_rate = miss_summary$observed_missing_rate,
    complete_adjacent_pair_rate = miss_summary$complete_adjacent_pair_rate,
    complete_adjacent_pairs = miss_summary$complete_adjacent_pairs,
    n_dropout_subjects = miss_summary$n_dropout_subjects,
    n_subjects_with_interior_gaps = miss_summary$n_subjects_with_interior_gaps,
    no_observations_after_dropout = miss_summary$no_observations_after_dropout,
    model = "gamlss.longitudinal", success = !is.null(fit_long), retained = retained_long,
    logLik = if (!is.null(fit_long)) extract_fit_metric(fit_long, "logLik") else NA_real_,
    df = if (!is.null(fit_long)) extract_fit_metric(fit_long, "df") else NA_real_,
    converged = conv_long$converged,
    hit_outer_limit = conv_long$hit_outer_limit,
    hit_max_stall = conv_long$hit_max_stall,
    hit_raw_loglik_deterioration = conv_long$hit_raw_loglik_deterioration,
    stop_reason = conv_long$stop_reason,
    grad_inf = conv_long$grad_inf,
    step_l2 = conv_long$step_l2,
    best_raw_loglik = conv_long$best_raw_loglik,
    best_raw_loglik_iteration = conv_long$best_raw_loglik_iteration,
    raw_loglik_drop_from_best = conv_long$raw_loglik_drop_from_best,
    raw_loglik_drop_tol = conv_long$raw_loglik_drop_tol,
    outer_iterations = conv_long$outer_iterations,
    outer_log_lik_change = conv_long$outer_log_lik_change,
    outer_stop_crit = conv_long$outer_stop_crit,
    elapsed_sec = elapsed_long, error = if (is.null(fit_long)) err_long else NA_character_,
    stringsAsFactors = FALSE
  )
  run_rows[[2]] <- data.frame(
    scenario = scenario_label, n = n_val, d = d_val, rep = rep_id,
    simulation_seed = seed, missingness_seed = missingness_seed,
    missing_mechanism = missing_mechanism,
    missingness_label = miss_summary$missingness_label,
    missingness_pattern = miss_summary$missingness_pattern,
    analysis_role = miss_summary$analysis_role,
    target_missing_rate = missing_rate,
    observed_missing_rate = miss_summary$observed_missing_rate,
    complete_adjacent_pair_rate = miss_summary$complete_adjacent_pair_rate,
    complete_adjacent_pairs = miss_summary$complete_adjacent_pairs,
    n_dropout_subjects = miss_summary$n_dropout_subjects,
    n_subjects_with_interior_gaps = miss_summary$n_subjects_with_interior_gaps,
    no_observations_after_dropout = miss_summary$no_observations_after_dropout,
    model = "gamlss2", success = !is.null(fit_gamlss2), retained = retained_gamlss2,
    logLik = if (!is.null(fit_gamlss2)) extract_fit_metric(fit_gamlss2, "logLik") else NA_real_,
    df = if (!is.null(fit_gamlss2)) extract_fit_metric(fit_gamlss2, "df") else NA_real_,
    converged = conv_gamlss2$converged,
    hit_outer_limit = NA,
    hit_max_stall = NA,
    hit_raw_loglik_deterioration = NA,
    stop_reason = conv_gamlss2$stop_reason,
    grad_inf = NA_real_,
    step_l2 = NA_real_,
    best_raw_loglik = NA_real_,
    best_raw_loglik_iteration = NA_integer_,
    raw_loglik_drop_from_best = NA_real_,
    raw_loglik_drop_tol = NA_real_,
    outer_iterations = conv_gamlss2$iterations,
    outer_log_lik_change = NA_real_,
    outer_stop_crit = NA_real_,
    elapsed_sec = elapsed_gamlss2, error = if (is.null(fit_gamlss2)) err_gamlss2 else NA_character_,
    stringsAsFactors = FALSE
  )

  if (retained_long) {
    fixed_long <- extract_fixed_estimates_longitudinal(fit_long, dat)
    if (identical(miss_summary$missingness_pattern[[1L]], "intermittent")) {
      fixed_long$std_error <- NA_real_
      fixed_long$inference_status <- "not_applicable_segmented_model_hessian"
    }
    fixed_long$scenario <- scenario_label
    fixed_long$n <- n_val
    fixed_long$d <- d_val
    fixed_long$rep <- rep_id
    fixed_long$missing_mechanism <- missing_mechanism
    fixed_long$target_missing_rate <- missing_rate
    fixed_long$observed_missing_rate <- miss_summary$observed_missing_rate
    fixed_long$complete_adjacent_pair_rate <- miss_summary$complete_adjacent_pair_rate
    fixed_rows[[length(fixed_rows) + 1]] <- fixed_long

    for (p in smooth_params_longitudinal) {
      true_curve <- smooth_truth[[p]](s1_grid) - jss_missing_population_smooth_mean(smooth_truth[[p]])
      smooth_rows[[length(smooth_rows) + 1]] <- data.frame(
        scenario = scenario_label, n = n_val, d = d_val, rep = rep_id,
        missing_mechanism = missing_mechanism,
        target_missing_rate = missing_rate,
        observed_missing_rate = miss_summary$observed_missing_rate,
        complete_adjacent_pair_rate = miss_summary$complete_adjacent_pair_rate,
        model = "gamlss.longitudinal", parameter = p, s1 = s1_grid,
        smooth_hat = smooth_curve_longitudinal(fit_long, dat, p, s1_grid),
        smooth_true = true_curve,
        stringsAsFactors = FALSE
      )
    }

    joint_rows[[length(joint_rows) + 1]] <- tryCatch(
      joint_metrics_longitudinal(fit_long, dat, scenario_label, n_val, d_val, rep_id),
      error = function(e) NULL
    )

    if (isTRUE(compute_predictive_scores)) {
      test_dat <- simulate_response_given_covariates(dat, seed = seed + 500000)
      predictive_rows[[length(predictive_rows) + 1]] <- tryCatch(
        predictive_scores_longitudinal(fit_long, test_dat, scenario_label, n_val, d_val, rep_id),
        error = function(e) NULL
      )
    }
  }

  if (retained_gamlss2) {
    fixed_gamlss2 <- extract_fixed_estimates_gamlss2(fit_gamlss2, dat)
    fixed_gamlss2$scenario <- scenario_label
    fixed_gamlss2$n <- n_val
    fixed_gamlss2$d <- d_val
    fixed_gamlss2$rep <- rep_id
    fixed_gamlss2$missing_mechanism <- missing_mechanism
    fixed_gamlss2$target_missing_rate <- missing_rate
    fixed_gamlss2$observed_missing_rate <- miss_summary$observed_missing_rate
    fixed_gamlss2$complete_adjacent_pair_rate <- miss_summary$complete_adjacent_pair_rate
    fixed_rows[[length(fixed_rows) + 1]] <- fixed_gamlss2

    for (p in smooth_params_gamlss) {
      true_curve <- smooth_truth[[p]](s1_grid) - jss_missing_population_smooth_mean(smooth_truth[[p]])
      smooth_rows[[length(smooth_rows) + 1]] <- data.frame(
        scenario = scenario_label, n = n_val, d = d_val, rep = rep_id,
        missing_mechanism = missing_mechanism,
        target_missing_rate = missing_rate,
        observed_missing_rate = miss_summary$observed_missing_rate,
        complete_adjacent_pair_rate = miss_summary$complete_adjacent_pair_rate,
        model = "gamlss2", parameter = p, s1 = s1_grid,
        smooth_hat = smooth_curve_gamlss2(fit_gamlss2, dat, p, s1_grid),
        smooth_true = true_curve,
        stringsAsFactors = FALSE
      )
    }

    joint_rows[[length(joint_rows) + 1]] <- tryCatch(
      joint_metrics_gamlss2(fit_gamlss2, dat, scenario_label, n_val, d_val, rep_id),
      error = function(e) NULL
    )

    if (isTRUE(compute_predictive_scores)) {
      if (!exists("test_dat", inherits = FALSE)) {
        test_dat <- simulate_response_given_covariates(dat, seed = seed + 500000)
      }
      predictive_rows[[length(predictive_rows) + 1]] <- tryCatch(
        predictive_scores_gamlss2(fit_gamlss2, test_dat, scenario_label, n_val, d_val, rep_id),
        error = function(e) NULL
      )
    }
  }

  add_missingness_columns <- function(x) {
    if (is.null(x) || nrow(x) == 0L) {
      return(x)
    }
    x$missing_mechanism <- missing_mechanism
    x$target_missing_rate <- missing_rate
    x$observed_missing_rate <- miss_summary$observed_missing_rate
    x$complete_adjacent_pair_rate <- miss_summary$complete_adjacent_pair_rate
    x
  }

  runs <- do.call(rbind, run_rows)
  runs$failure_type <- ifelse(runs$success %in% TRUE & runs$converged %in% TRUE, "none",
    ifelse(runs$success %in% TRUE, paste0("optimizer_nonconvergence:", runs$stop_reason),
      ifelse(runs$stop_reason == "fit_error", "fit_error", paste0("fit_error:", runs$stop_reason))))
  out <- list(
    checkpoint_schema_version = 5L,
    checkpoint_configuration_key = checkpoint_configuration_key,
    checkpoint_spec = jss_missing_checkpoint_spec(task, checkpoint_configuration),
    fixed = bind_non_null(fixed_rows),
    smooth = bind_non_null(smooth_rows),
    joint = add_missingness_columns(bind_non_null(joint_rows)),
    predictive = add_missingness_columns(bind_non_null(predictive_rows)),
    runs = runs,
    missingness = miss_summary,
    missingness_pattern = missingness_pattern,
    warning_events = warning_events,
    checkpoint_provenance = c(jss_missing_runtime_identity(), list(
      package_identity_verified = isTRUE(get(".jss_missing_verified_worker_identity", envir = .GlobalEnv,
        inherits = TRUE)$verified),
      package_source_sha256 = expected_package_identity$source_sha256,
      producer_sha256_verified = producer_sha256
    ))
  )
  if (isTRUE(save_fits)) {
    out$fit_longitudinal <- fit_long
    out$fit_gamlss2 <- fit_gamlss2
    out$data <- dat
  }
  out$checkpoint_content_sha256 <- jss_missing_content_sha256(jss_missing_checkpoint_content(out))
  out
}

tasks <- do.call(
  rbind,
  lapply(seq_len(nrow(scenarios)), function(i) {
    data.frame(
      scenario_id = as.integer(i),
      n = scenarios$n[i],
      d = scenarios$d[i],
      missing_mechanism = scenarios$missing_mechanism[i],
      missing_rate = scenarios$missing_rate[i],
      rep = rep_ids,
      checkpoint_configuration_key = checkpoint_configuration_key,
      stringsAsFactors = FALSE
    )
  })
)
tasks <- split(tasks, seq_len(nrow(tasks)))
checkpoint_results <- if (isTRUE(resume_checkpoints)) {
  lapply(tasks, function(task) {
    path <- task_result_path(task)
    if (!file.exists(path)) return(NULL)
    reverify_missingness_checkpoint_or_quarantine(path)
    result <- tryCatch(readRDS(path), error = function(e) NULL)
    reverify_missingness_checkpoint_or_quarantine(path)
    if (jss_missing_checkpoint_valid(result, task, checkpoint_configuration)) return(result)
    quarantine <- paste0(path, ".invalid-parent-", Sys.getpid(), "-", format(Sys.time(), "%Y%m%d%H%M%S"))
    if (!file.rename(path, quarantine)) file.remove(path)
    NULL
  })
} else {
  rep(list(NULL), length(tasks))
}
completed <- !vapply(checkpoint_results, is.null, logical(1))
results <- checkpoint_results[completed]
tasks <- tasks[!completed]

cat(sprintf(
  "Running %d fit replicate(s) across %d missingness scenario(s): %d checkpoint(s) complete, %d pending, using up to %d core(s).\n",
  length(rep_ids), nrow(scenarios), sum(completed), length(tasks), min(n_cores, max(1L, length(tasks)))
))

if (length(tasks) == 0L) {
  new_results <- list()
} else if (length(tasks) == 1 || n_cores <= 1) {
  new_results <- lapply(tasks, run_one_rep_and_save)
} else {
  rscript_bin <- Sys.getenv(
    "R_SCRIPT",
    unset = file.path(R.home("bin"), if (.Platform$OS.type == "windows") "Rscript.exe" else "Rscript")
  )
  n_workers <- min(n_cores, length(tasks))
  worker_log <- file.path(worker_log_dir, "psock-workers.log")
  cat(sprintf(
    "Starting PSOCK cluster with %d worker(s). Worker log: %s\n",
    n_workers, worker_log
  ))
  cl <- parallel::makePSOCKcluster(
    n_workers,
    rscript = rscript_bin,
    rscript_args = "--vanilla",
    outfile = worker_log,
    setup_strategy = "parallel"
  )

  on.exit(if (!is.null(cl)) parallel::stopCluster(cl), add = TRUE)
  parallel::clusterCall(cl, setwd, project_root)
  worker_status <- parallel::clusterCall(cl, function(project_root, expected_identity, expected_producer_sha256) {
    safe_source <- function(path) {
      txt <- readLines(path, warn = FALSE, encoding = "UTF-8")
      if (length(txt) > 0) {
        bom <- intToUtf8(65279)
        txt[1] <- sub(paste0("^", bom), "", txt[1], useBytes = FALSE)
      }
      eval(parse(text = txt), envir = .GlobalEnv)
    }
    if (!requireNamespace("pkgload", quietly = TRUE)) stop("pkgload is required on missingness workers.")
    pkgload::load_all(project_root, quiet = TRUE, export_all = TRUE, helpers = FALSE)
    suppressPackageStartupMessages(library(gamlss.longitudinal))
    list2env(as.list(getNamespace("gamlss.longitudinal"), all.names = TRUE), envir = .GlobalEnv)
    safe_source(file.path(project_root, "paper", "R", "missingness-study-helpers.R"))
    verified <- jss_missing_verify_checkout(project_root, expected_identity, load_checkout = FALSE)
    producer_sources <- c(
      file.path(project_root, "paper", "R", "missingness-study-helpers.R"),
      file.path(project_root, "paper", "scripts", "final-simulations", "missingness", "run_missingness_study.R")
    )
    actual_producer <- digest::digest(
      paste(basename(producer_sources), vapply(producer_sources, jss_missing_sha256_file, character(1)), collapse = "\n"),
      algo = "sha256", serialize = FALSE
    )
    if (!identical(actual_producer, expected_producer_sha256)) {
      stop("Missingness producer changed between parent setup and worker attestation.")
    }
    assign(".jss_missing_verified_worker_identity", verified, envir = .GlobalEnv)
    suppressPackageStartupMessages({
      library(VineCopula)
      library(gamlss)
      library(gamlss2)
      library(gamlss.dist)
      library(ggplot2)
    })
    list(pid = Sys.getpid(), wd = getwd(), package_identity = verified,
      producer_sha256 = actual_producer)
  }, project_root = project_root, expected_identity = expected_package_identity,
  expected_producer_sha256 = producer_sha256)
  cat("Worker setup complete:\n")
  print(worker_status)

  export_names <- setdiff(
    ls(envir = .GlobalEnv),
    c("cl", "results", "fixed_all", "smooth_all", "joint_all", "runs_all")
  )
  parallel::clusterExport(cl, varlist = export_names, envir = .GlobalEnv)

  if (isTRUE(parallel_setup_only)) {
    task_test <- parallel::parLapplyLB(cl, seq_len(n_workers), function(i) {
      list(i = i, pid = Sys.getpid(), wd = getwd())
    })
    cat("Parallel export and scheduling test complete:\n")
    print(task_test)
    cat("PARALLEL_SETUP_ONLY=1; stopping after successful worker setup.\n")
    parallel::stopCluster(cl)
    cl <- NULL
    quit(save = "no", status = 0)
  }

  new_results <- parallel::parLapplyLB(cl, tasks, run_one_rep_and_save)
  parallel::stopCluster(cl)
  cl <- NULL
}
results <- c(results, new_results)

checkpoint_status <- do.call(rbind, lapply(results, function(result) {
  provenance <- result$checkpoint_provenance
  spec <- result$checkpoint_spec
  checkpoint_name <- sprintf(
    "%s_scenario%d_%s_miss%02d_n%d_d%d_rep%d.rds",
    result$checkpoint_configuration_key, spec$scenario_id, spec$missing_mechanism,
    round(100 * spec$missing_rate), spec$n, spec$d, spec$replicate)
  data.frame(
    checkpoint_schema_version = result$checkpoint_schema_version,
    scenario = spec$scenario, rep = spec$replicate,
    checkpoint = gsub("\\\\", "/", file.path("rep_results", checkpoint_name)),
    checkpoint_content_sha256 = result$checkpoint_content_sha256,
    public_payload_sha256 = jss_missing_portable_task_sha256(
      result$runs, result$fixed, result$smooth, result$missingness,
      result$missingness_pattern, result$warning_events
    ),
    package_source_sha256 = result$checkpoint_spec$configuration$package_identity$source_sha256,
    package_version = result$checkpoint_spec$configuration$package_identity$version,
    package_fingerprint_scope = result$checkpoint_spec$configuration$package_identity$fingerprint_scope,
    package_source_file_count = result$checkpoint_spec$configuration$package_identity$source_file_count,
    producer_sha256 = result$checkpoint_spec$configuration$producer_sha256,
    package_identity_verified = isTRUE(provenance$package_identity_verified),
    timestamp_utc = provenance$timestamp_utc, worker_pid = provenance$pid,
    host = provenance$host, os = provenance$os, platform = provenance$platform,
    r_version = provenance$r_version, rng_kind = provenance$rng_kind,
    blas = provenance$blas, lapack = provenance$lapack,
    rlibs_user_sha256 = provenance$rlibs_user_sha256,
    libpaths_sha256 = provenance$libpaths_sha256,
    stringsAsFactors = FALSE
  )
}))
write.csv(checkpoint_status, file.path(out_dir, "missingness_checkpoint_status.csv"), row.names = FALSE)
saveRDS(results, file.path(out_dir, "missingness_checkpoint_payloads.rds"), version = 3L)
checkpoint_content_manifest <- checkpoint_status[c("scenario", "rep", "checkpoint",
  "checkpoint_content_sha256", "public_payload_sha256")]
write.csv(checkpoint_content_manifest,
  file.path(out_dir, "missingness_checkpoint_content_manifest.csv"), row.names = FALSE)

runs_all <- bind_non_null(lapply(results, `[[`, "runs"))
missingness_all <- bind_non_null(lapply(results, `[[`, "missingness"))
missingness_pattern_all <- bind_non_null(lapply(results, `[[`, "missingness_pattern"))
warning_events_all <- jss_missing_collect_warning_events(results)
jss_missing_assert_warning_policy(warning_events_all)
runs_all$retained <- runs_all$success %in% TRUE & runs_all$converged %in% TRUE
fixed_all <- jss_missing_filter_payload(bind_non_null(lapply(results, `[[`, "fixed")), runs_all)
smooth_all <- jss_missing_filter_payload(bind_non_null(lapply(results, `[[`, "smooth")), runs_all)
joint_all <- jss_missing_filter_payload(bind_non_null(lapply(results, `[[`, "joint")), runs_all)
predictive_all <- jss_missing_filter_payload(bind_non_null(lapply(results, `[[`, "predictive")), runs_all)
public_hash_map <- setNames(checkpoint_status$public_payload_sha256,
  paste(checkpoint_status$scenario, checkpoint_status$rep, sep = "\r"))
attach_public_hash <- function(x) {
  if (is.null(x)) return(NULL)
  x$public_payload_sha256 <- unname(public_hash_map[paste(x$scenario, x$rep, sep = "\r")])
  x
}
runs_all <- attach_public_hash(runs_all)
missingness_all <- attach_public_hash(missingness_all)
fixed_all <- attach_public_hash(fixed_all)
smooth_all <- attach_public_hash(smooth_all)
missingness_pattern_all <- attach_public_hash(missingness_pattern_all)
warning_events_all <- attach_public_hash(warning_events_all)

saveRDS(
  list(
    fixed = fixed_all,
    smooth = smooth_all,
    joint = joint_all,
    predictive = predictive_all,
    runs = runs_all,
    missingness = missingness_all,
    missingness_pattern = missingness_pattern_all,
    warning_events = warning_events_all,
    truth = true_beta
  ),
  file.path(out_dir, "all_results.rds")
)
failure_rmse_penalties <- as.numeric(strsplit(
  Sys.getenv("FAILURE_RMSE_PENALTIES", unset = "0.5,1,2"),
  ",",
  fixed = TRUE
)[[1L]])
failure_rmse_penalties <- sort(unique(failure_rmse_penalties[is.finite(failure_rmse_penalties) & failure_rmse_penalties > 0]))
failure_rmse_penalty <- as.numeric(Sys.getenv("FAILURE_RMSE_PENALTY", unset = "1"))
if (!is.finite(failure_rmse_penalty) || failure_rmse_penalty <= 0) {
  stop("FAILURE_RMSE_PENALTY must be a positive finite link-scale value.", call. = FALSE)
}
failure_rmse_penalties <- sort(unique(c(failure_rmse_penalties, failure_rmse_penalty)))
penalty_registry <- data.frame(
  penalty = failure_rmse_penalties,
  role = ifelse(failure_rmse_penalties == failure_rmse_penalty, "primary", "sensitivity"),
  scale = "link-scale RMSE/IRMSE",
  application = "assigned only when an attempted fit has no retained estimate",
  justification = paste(
    "Penalty values are not estimands; they show how conclusions change when failed or",
    "nonconverged attempts receive prespecified moderate, primary, and severe losses."
  ),
  stringsAsFactors = FALSE
)
write.csv(penalty_registry, file.path(out_dir, "missingness_sensitivity_registry.csv"), row.names = FALSE)
write.csv(runs_all, file.path(out_dir, "fit_run_log.csv"), row.names = FALSE)
write.csv(warning_events_all, file.path(out_dir, "missingness_warning_events.csv"), row.names = FALSE)
fit_attempt_summary <- do.call(rbind, lapply(
  split(runs_all, list(runs_all$scenario, runs_all$model), drop = TRUE),
  function(df) data.frame(
    scenario = df$scenario[[1L]],
    model = df$model[[1L]],
    missing_mechanism = df$missing_mechanism[[1L]],
    missingness_label = df$missingness_label[[1L]],
    missingness_pattern = df$missingness_pattern[[1L]],
    analysis_role = df$analysis_role[[1L]],
    target_missing_rate = df$target_missing_rate[[1L]],
    attempted = nrow(df),
    fit_successful = sum(df$success %in% TRUE),
    converged = sum(df$converged %in% TRUE),
    retained = sum(df$retained %in% TRUE),
    failed = sum(!(df$retained %in% TRUE)),
    error_failures = sum(!(df$success %in% TRUE)),
    nonconvergence_failures = sum(df$success %in% TRUE & df$converged %in% FALSE),
    failure_inclusive_retention_rate = mean(df$retained %in% TRUE),
    retention_rate_mcse = sqrt(mean(df$retained %in% TRUE) * (1 - mean(df$retained %in% TRUE)) / nrow(df)),
    retention_rate_conf_low = max(
      0,
      mean(df$retained %in% TRUE) - 1.96 * sqrt(mean(df$retained %in% TRUE) * (1 - mean(df$retained %in% TRUE)) / nrow(df))
    ),
    retention_rate_conf_high = min(
      1,
      mean(df$retained %in% TRUE) + 1.96 * sqrt(mean(df$retained %in% TRUE) * (1 - mean(df$retained %in% TRUE)) / nrow(df))
    ),
    failure_rmse_penalty = failure_rmse_penalty,
    stringsAsFactors = FALSE
  )
))
write.csv(fit_attempt_summary, file.path(out_dir, "attempt_failure_summary.csv"), row.names = FALSE)
failure_reason <- ifelse(
  runs_all$retained,
  "retained",
  ifelse(
    !runs_all$success,
    runs_all$stop_reason,
    ifelse(runs_all$converged %in% FALSE, paste0("nonconverged: ", runs_all$stop_reason), "not_retained")
  )
)
failure_reason_summary <- stats::aggregate(
  rep(1L, nrow(runs_all)),
  by = data.frame(
    scenario = runs_all$scenario,
    model = runs_all$model,
    failure_reason = failure_reason,
    stringsAsFactors = FALSE
  ),
  FUN = sum
)
names(failure_reason_summary)[names(failure_reason_summary) == "x"] <- "attempts"
write.csv(failure_reason_summary, file.path(out_dir, "failure_reason_summary.csv"), row.names = FALSE)
if (!is.null(missingness_all) && nrow(missingness_all) > 0) {
  write.csv(missingness_all, file.path(out_dir, "missingness_by_rep.csv"), row.names = FALSE)
}
if (!is.null(missingness_pattern_all) && nrow(missingness_pattern_all) > 0) {
  write.csv(missingness_pattern_all, file.path(out_dir, "missingness_pattern_by_subject_visit.csv"), row.names = FALSE)
}
if (!is.null(fixed_all) && nrow(fixed_all) > 0) {
  write.csv(fixed_all, file.path(out_dir, "fixed_effects_by_rep.csv"), row.names = FALSE)
}
if (!is.null(smooth_all) && nrow(smooth_all) > 0) {
  write.csv(smooth_all, file.path(out_dir, "smooth_estimates_by_rep.csv"), row.names = FALSE)
}
if (!is.null(joint_all) && nrow(joint_all) > 0) {
  write.csv(joint_all, file.path(out_dir, "joint_distribution_metrics_by_rep.csv"), row.names = FALSE)
}
if (!is.null(predictive_all) && nrow(predictive_all) > 0) {
  write.csv(predictive_all, file.path(out_dir, "predictive_scores_by_rep.csv"), row.names = FALSE)
}

if (is.null(fixed_all) || nrow(fixed_all) == 0) {
  warning("No successful fits were available for fixed-effect summaries.")
  quit(save = "no", status = 1)
}

fixed_summary <- do.call(rbind, lapply(
  split(fixed_all, list(fixed_all$scenario, fixed_all$model, fixed_all$parameter, fixed_all$term), drop = TRUE),
  function(df) {
    err <- df$estimate - df$true_value
    ci_low <- df$estimate - stats::qnorm(0.975) * df$std_error
    ci_high <- df$estimate + stats::qnorm(0.975) * df$std_error
    coverage_available <- is.finite(df$true_value) & is.finite(ci_low) & is.finite(ci_high)
    covered <- rep(NA, nrow(df))
    covered[coverage_available] <-
      df$true_value[coverage_available] >= ci_low[coverage_available] &
      df$true_value[coverage_available] <= ci_high[coverage_available]
    sd_est <- stats::sd(df$estimate, na.rm = TRUE)
    mean_se <- mean(df$std_error, na.rm = TRUE)
    data.frame(
      scenario = df$scenario[1],
      model = df$model[1],
      n = df$n[1],
      d = df$d[1],
      parameter = df$parameter[1],
      term = df$term[1],
      true_value = df$true_value[1],
      mean_estimate = mean(df$estimate, na.rm = TRUE),
      bias = mean(err, na.rm = TRUE),
      rmse = sqrt(mean(err^2, na.rm = TRUE)),
      sd_estimate = sd_est,
      mean_std_error = mean_se,
      se_to_empirical_sd = mean_se / sd_est,
      coverage_95 = if (any(coverage_available)) mean(covered, na.rm = TRUE) else NA_real_,
      q05 = safe_quantile(df$estimate, 0.05),
      q95 = safe_quantile(df$estimate, 0.95),
      n_successful_fits = sum(is.finite(df$estimate)),
      n_se = sum(is.finite(df$std_error)),
      inference_status_raw = if ("inference_status" %in% names(df) &&
        length(unique(df$inference_status)) == 1L) unique(df$inference_status) else "mixed_or_unregistered",
      stringsAsFactors = FALSE
    )
  }
))

fit_metric_summary <- do.call(rbind, lapply(
  split(runs_all, list(runs_all$scenario, runs_all$model), drop = TRUE),
  function(df) {
    retained_df <- df[df$retained %in% TRUE, , drop = FALSE]
    data.frame(
      scenario = df$scenario[1],
      model = df$model[1],
      mean_logLik = if (nrow(retained_df)) mean(retained_df$logLik, na.rm = TRUE) else NA_real_,
      mean_df = if (nrow(retained_df)) mean(retained_df$df, na.rm = TRUE) else NA_real_,
      stringsAsFactors = FALSE
    )
  }
))

fixed_summary <- merge(
  fixed_summary,
  fit_metric_summary,
  by = c("scenario", "model"),
  all.x = TRUE,
  sort = FALSE
)
fixed_summary <- merge(
  fixed_summary,
  fit_attempt_summary[, c("scenario", "model", "attempted", "fit_successful", "retained", "failed")],
  by = c("scenario", "model"),
  all.x = TRUE,
  sort = FALSE
)
fixed_summary$conditional_rmse <- fixed_summary$rmse
fixed_summary$conditional_coverage_95 <- fixed_summary$coverage_95
fixed_summary$metric_failures <- pmax(fixed_summary$attempted - fixed_summary$n_successful_fits, 0)
fixed_summary$failure_inclusive_rmse_penalty <- failure_rmse_penalty
fixed_summary$failure_inclusive_rmse <- sqrt(
  (fixed_summary$rmse^2 * fixed_summary$n_successful_fits +
    failure_rmse_penalty^2 * fixed_summary$metric_failures) /
    fixed_summary$attempted
)
fixed_summary$failure_inclusive_coverage_95 <- ifelse(
  fixed_summary$n_se > 0,
  fixed_summary$coverage_95 * fixed_summary$n_se / fixed_summary$attempted,
  0
)

scenario_meta <- do.call(rbind, lapply(
  split(runs_all, runs_all$scenario, drop = TRUE),
  function(df) {
    data.frame(
      scenario = df$scenario[1],
      missing_mechanism = df$missing_mechanism[1],
      missingness_label = df$missingness_label[1],
      missingness_pattern = df$missingness_pattern[1],
      analysis_role = df$analysis_role[1],
      target_missing_rate = df$target_missing_rate[1],
      observed_missing_rate = mean(df$observed_missing_rate, na.rm = TRUE),
      complete_adjacent_pair_rate = mean(df$complete_adjacent_pair_rate, na.rm = TRUE),
      stringsAsFactors = FALSE
    )
  }
))
fixed_summary <- merge(
  fixed_summary,
  scenario_meta,
  by = "scenario",
  all.x = TRUE,
  sort = FALSE
)
fixed_summary$coverage_status <- ifelse(
  fixed_summary$model == "gamlss.longitudinal" & fixed_summary$missingness_pattern == "intermittent",
  "not_applicable_segmented_model_hessian",
  ifelse(
    fixed_summary$inference_status_raw == "not_available_without_joint_fixed_smooth_covariance",
    "not_available_without_joint_fixed_smooth_covariance",
    ifelse(fixed_summary$n_se > 0, "conditional_on_retained_fits", "not_available")
  )
)
coverage_not_applicable <- fixed_summary$coverage_status %in% c(
  "not_applicable_segmented_model_hessian",
  "not_available_without_joint_fixed_smooth_covariance"
)
fixed_summary$coverage_95[coverage_not_applicable] <- NA_real_
fixed_summary$conditional_coverage_95[coverage_not_applicable] <- NA_real_
fixed_summary$failure_inclusive_coverage_95[coverage_not_applicable] <- NA_real_
fixed_summary$coverage_mcse <- ifelse(
  fixed_summary$n_se > 0 & !coverage_not_applicable,
  sqrt(fixed_summary$coverage_95 * (1 - fixed_summary$coverage_95) / fixed_summary$n_se),
  NA_real_
)
fixed_summary$coverage_conf_low <- pmax(0, fixed_summary$coverage_95 - 1.96 * fixed_summary$coverage_mcse)
fixed_summary$coverage_conf_high <- pmin(1, fixed_summary$coverage_95 + 1.96 * fixed_summary$coverage_mcse)

fixed_rep_rmse <- fixed_all[fixed_all$term != "intercept", , drop = FALSE]
fixed_rep_rmse$squared_error <- (fixed_rep_rmse$estimate - fixed_rep_rmse$true_value)^2
fixed_rep_rmse <- stats::aggregate(
  fixed_rep_rmse$squared_error,
  by = fixed_rep_rmse[c("scenario", "model", "parameter", "rep")],
  FUN = function(x) sqrt(mean(x, na.rm = TRUE))
)
names(fixed_rep_rmse)[names(fixed_rep_rmse) == "x"] <- "replicate_rmse"
fixed_rep_rmse <- fixed_rep_rmse[is.finite(fixed_rep_rmse$replicate_rmse), , drop = FALSE]
headline_rows <- list()
headline_key <- interaction(
  fixed_rep_rmse$scenario, fixed_rep_rmse$model, fixed_rep_rmse$parameter,
  drop = TRUE, lex.order = TRUE
)
for (idx in split(seq_len(nrow(fixed_rep_rmse)), headline_key)) {
  df <- fixed_rep_rmse[idx, , drop = FALSE]
  attempts <- fit_attempt_summary[
    fit_attempt_summary$scenario == df$scenario[[1L]] &
      fit_attempt_summary$model == df$model[[1L]],
    ,
    drop = FALSE
  ]
  meta <- scenario_meta[scenario_meta$scenario == df$scenario[[1L]], , drop = FALSE]
  conditional <- mc_mean_interval(df$replicate_rmse)
  for (penalty in failure_rmse_penalties) {
    failures <- max(attempts$attempted[[1L]] - conditional[["n"]], 0L)
    inclusive <- mc_mean_interval(c(df$replicate_rmse, rep(penalty, failures)))
    headline_rows[[length(headline_rows) + 1L]] <- data.frame(
      scenario = df$scenario[[1L]], model = df$model[[1L]], parameter = df$parameter[[1L]],
      missing_mechanism = meta$missing_mechanism[[1L]],
      missingness_label = meta$missingness_label[[1L]],
      missingness_pattern = meta$missingness_pattern[[1L]], analysis_role = meta$analysis_role[[1L]],
      target_missing_rate = meta$target_missing_rate[[1L]],
      attempted = attempts$attempted[[1L]], retained = conditional[["n"]], failed = failures,
      conditional_mean_rmse = conditional[["estimate"]], conditional_rmse_mcse = conditional[["mcse"]],
      conditional_rmse_conf_low = conditional[["conf_low"]], conditional_rmse_conf_high = conditional[["conf_high"]],
      failure_penalty = penalty,
      failure_inclusive_mean_rmse = inclusive[["estimate"]],
      failure_inclusive_rmse_mcse = inclusive[["mcse"]],
      failure_inclusive_rmse_conf_low = inclusive[["conf_low"]],
      failure_inclusive_rmse_conf_high = inclusive[["conf_high"]],
      inference_status = if (
        identical(df$model[[1L]], "gamlss.longitudinal") &&
          identical(meta$missingness_pattern[[1L]], "intermittent")
      ) "not_applicable_segmented_model_hessian" else "conditional_on_retained_fits",
      stringsAsFactors = FALSE
    )
  }
}
missingness_headline_summary <- bind_non_null(headline_rows)
write.csv(missingness_headline_summary, file.path(out_dir, "missingness_headline_summary.csv"), row.names = FALSE)

write.csv(fixed_summary, file.path(out_dir, "fixed_effects_bias_rmse_table.csv"), row.names = FALSE)

fixed_compact <- reshape(
  fixed_summary[, c("scenario", "missing_mechanism", "target_missing_rate", "observed_missing_rate", "complete_adjacent_pair_rate", "n", "d", "parameter", "term", "model", "bias", "rmse", "mean_std_error", "sd_estimate", "se_to_empirical_sd", "coverage_95", "mean_logLik", "mean_df")],
  idvar = c("scenario", "missing_mechanism", "target_missing_rate", "observed_missing_rate", "complete_adjacent_pair_rate", "n", "d", "parameter", "term"),
  timevar = "model",
  direction = "wide"
)
write.csv(fixed_compact, file.path(out_dir, "fixed_effects_bias_rmse_compact.csv"), row.names = FALSE)

se_calibration <- fixed_summary[fixed_summary$term != "intercept", c(
  "scenario", "missing_mechanism", "target_missing_rate", "observed_missing_rate",
  "complete_adjacent_pair_rate", "model", "n", "d", "parameter", "term",
  "mean_std_error", "sd_estimate", "se_to_empirical_sd", "coverage_95", "n_se"
)]
write.csv(se_calibration, file.path(out_dir, "fixed_effects_se_calibration.csv"), row.names = FALSE)

method_summary <- do.call(rbind, lapply(
  split(runs_all, list(runs_all$scenario, runs_all$model), drop = TRUE),
  function(df) {
    fixed_df <- fixed_summary[fixed_summary$scenario == df$scenario[1] & fixed_summary$model == df$model[1], , drop = FALSE]
    fixed_margin_df <- fixed_df[fixed_df$parameter %in% params_margin & fixed_df$term != "intercept", , drop = FALSE]
    copula_df <- fixed_df[fixed_df$parameter %in% params_copula & fixed_df$term != "intercept", , drop = FALSE]
    data.frame(
      scenario = df$scenario[1],
      missing_mechanism = df$missing_mechanism[1],
      missingness_label = df$missingness_label[1],
      missingness_pattern = df$missingness_pattern[1],
      analysis_role = df$analysis_role[1],
      target_missing_rate = df$target_missing_rate[1],
      observed_missing_rate = mean(df$observed_missing_rate, na.rm = TRUE),
      complete_adjacent_pair_rate = mean(df$complete_adjacent_pair_rate, na.rm = TRUE),
      model = df$model[1],
      attempted = nrow(df),
      successful = sum(df$success == TRUE, na.rm = TRUE),
      retained = sum(df$retained %in% TRUE),
      failed = sum(!(df$retained %in% TRUE)),
      error_failures = sum(!(df$success %in% TRUE)),
      nonconvergence_failures = sum(df$success %in% TRUE & df$converged %in% FALSE),
      failure_inclusive_retention_rate = mean(df$retained %in% TRUE),
      convergence_rate = mean(df$converged, na.rm = TRUE),
      mean_elapsed_sec = mean(df$elapsed_sec, na.rm = TRUE),
      mean_outer_iterations = mean(df$outer_iterations, na.rm = TRUE),
      mean_margin_fixed_rmse = mean(fixed_margin_df$rmse, na.rm = TRUE),
      mean_margin_fixed_coverage_95 = mean(fixed_margin_df$coverage_95, na.rm = TRUE),
      mean_copula_fixed_rmse = mean(copula_df$rmse, na.rm = TRUE),
      mean_copula_fixed_coverage_95 = mean(copula_df$coverage_95, na.rm = TRUE),
      stringsAsFactors = FALSE
    )
  }
))
write.csv(method_summary, file.path(out_dir, "missingness_method_summary.csv"), row.names = FALSE)

if (!is.null(joint_all) && nrow(joint_all) > 0) {
  joint_summary <- do.call(rbind, lapply(
    split(joint_all, list(joint_all$scenario, joint_all$model), drop = TRUE),
    function(df) {
      data.frame(
        scenario = df$scenario[1],
        model = df$model[1],
        n = df$n[1],
        d = df$d[1],
        mean_logLik = mean(df$logLik, na.rm = TRUE),
        sd_logLik = stats::sd(df$logLik, na.rm = TRUE),
        mean_abs_rosenblatt_lag1_cor = mean(df$abs_rosenblatt_lag1_cor, na.rm = TRUE),
        mean_abs_rosenblatt_normal_lag1_cor = mean(df$abs_rosenblatt_normal_lag1_cor, na.rm = TRUE),
        mean_rosenblatt_mean_abs_time_mean = mean(df$rosenblatt_mean_abs_time_mean, na.rm = TRUE),
        mean_rosenblatt_normal_mean_abs_time_mean = mean(df$rosenblatt_normal_mean_abs_time_mean, na.rm = TRUE),
        stringsAsFactors = FALSE
      )
    }
  ))
  write.csv(joint_summary, file.path(out_dir, "joint_distribution_metrics_summary.csv"), row.names = FALSE)
}

if (!is.null(predictive_all) && nrow(predictive_all) > 0) {
  predictive_summary <- do.call(rbind, lapply(
    split(predictive_all, list(predictive_all$scenario, predictive_all$model), drop = TRUE),
    function(df) {
      data.frame(
        scenario = df$scenario[1],
        model = df$model[1],
        n = df$n[1],
        d = df$d[1],
        mean_test_log_score_joint = mean(df$test_log_score_joint, na.rm = TRUE),
        sd_test_log_score_joint = stats::sd(df$test_log_score_joint, na.rm = TRUE),
        mean_test_log_score_per_obs = mean(df$test_log_score_per_obs, na.rm = TRUE),
        mean_test_log_score_marginal = mean(df$test_log_score_marginal, na.rm = TRUE),
        mean_test_log_score_copula = mean(df$test_log_score_copula, na.rm = TRUE),
        mean_variogram_score = mean(df$variogram_score, na.rm = TRUE),
        sd_variogram_score = stats::sd(df$variogram_score, na.rm = TRUE),
        variogram_p = df$variogram_p[1],
        predictive_nsim = df$predictive_nsim[1],
        stringsAsFactors = FALSE
      )
    }
  ))
  write.csv(predictive_summary, file.path(out_dir, "predictive_scores_summary.csv"), row.names = FALSE)
}

if (!is.null(smooth_all) && nrow(smooth_all) > 0) {
  smooth_summary <- do.call(rbind, lapply(
    split(smooth_all, list(smooth_all$scenario, smooth_all$model, smooth_all$parameter, smooth_all$s1), drop = TRUE),
    function(df) {
      err <- df$smooth_hat - df$smooth_true
      data.frame(
        scenario = df$scenario[1],
        model = df$model[1],
        n = df$n[1],
        d = df$d[1],
        parameter = df$parameter[1],
        s1 = df$s1[1],
        smooth_true = df$smooth_true[1],
        smooth_median = stats::median(df$smooth_hat, na.rm = TRUE),
        smooth_q05 = safe_quantile(df$smooth_hat, 0.05),
        smooth_q95 = safe_quantile(df$smooth_hat, 0.95),
        bias = mean(err, na.rm = TRUE),
        rmse = sqrt(mean(err^2, na.rm = TRUE)),
        stringsAsFactors = FALSE
      )
    }
  ))
  write.csv(smooth_summary, file.path(out_dir, "smooth_pointwise_summary.csv"), row.names = FALSE)

  smooth_integrated <- do.call(rbind, lapply(
    split(smooth_all, list(smooth_all$scenario, smooth_all$model, smooth_all$parameter, smooth_all$rep), drop = TRUE),
    function(df) {
      err <- df$smooth_hat - df$smooth_true
      data.frame(
        scenario = df$scenario[1],
        model = df$model[1],
        n = df$n[1],
        d = df$d[1],
        parameter = df$parameter[1],
        rep = df$rep[1],
        bias_abs_integrated = mean(abs(err), na.rm = TRUE),
        irmse = sqrt(mean(err^2, na.rm = TRUE)),
        stringsAsFactors = FALSE
      )
    }
  ))
  smooth_integrated <- merge(
    smooth_integrated,
    fit_metric_summary,
    by = c("scenario", "model"),
    all.x = TRUE,
    sort = FALSE
  )
  write.csv(smooth_integrated, file.path(out_dir, "smooth_integrated_metrics.csv"), row.names = FALSE)

  smooth_integrated_meta <- merge(
    smooth_integrated,
    scenario_meta,
    by = "scenario",
    all.x = TRUE,
    sort = FALSE
  )
  smooth_irmse_summary <- do.call(rbind, lapply(
    split(
      smooth_integrated_meta,
      list(
        smooth_integrated_meta$missing_mechanism,
        smooth_integrated_meta$target_missing_rate,
        smooth_integrated_meta$model,
        smooth_integrated_meta$parameter
      ),
      drop = TRUE
    ),
    function(df) {
      attempt_row <- fit_attempt_summary[
        fit_attempt_summary$scenario == df$scenario[[1L]] &
          fit_attempt_summary$model == df$model[[1L]],
        ,
        drop = FALSE
      ]
      attempted <- attempt_row$attempted[[1L]]
      n_metric <- sum(is.finite(df$irmse))
      metric_failures <- pmax(attempted - n_metric, 0L)
      conditional_irmse <- mc_mean_interval(df$irmse)
      do.call(rbind, lapply(failure_rmse_penalties, function(penalty) {
        inclusive_irmse <- mc_mean_interval(c(df$irmse, rep(penalty, metric_failures)))
        data.frame(
        missing_mechanism = df$missing_mechanism[1],
        target_missing_rate = df$target_missing_rate[1],
        observed_missing_rate = mean(df$observed_missing_rate, na.rm = TRUE),
        complete_adjacent_pair_rate = mean(df$complete_adjacent_pair_rate, na.rm = TRUE),
        model = df$model[1],
        parameter = df$parameter[1],
        mean_irmse = conditional_irmse[["estimate"]],
        irmse_mcse = conditional_irmse[["mcse"]],
        irmse_conf_low = conditional_irmse[["conf_low"]],
        irmse_conf_high = conditional_irmse[["conf_high"]],
        sd_irmse = stats::sd(df$irmse, na.rm = TRUE),
        attempted = attempted,
        n_metric_successful = n_metric,
        metric_failures = metric_failures,
        failure_inclusive_mean_irmse = inclusive_irmse[["estimate"]],
        failure_inclusive_irmse_mcse = inclusive_irmse[["mcse"]],
        failure_inclusive_irmse_conf_low = inclusive_irmse[["conf_low"]],
        failure_inclusive_irmse_conf_high = inclusive_irmse[["conf_high"]],
        failure_rmse_penalty = penalty,
        stringsAsFactors = FALSE
        )
      }))
    }
  ))
  write.csv(smooth_irmse_summary, file.path(out_dir, "smooth_irmse_summary.csv"), row.names = FALSE)

  smooth_irmse_plot_data <- smooth_irmse_summary[
    smooth_irmse_summary$failure_rmse_penalty == failure_rmse_penalty, , drop = FALSE]
  smooth_irmse_plot_data$target_missing_pct <- 100 * smooth_irmse_plot_data$target_missing_rate
  smooth_irmse_plot_data$missing_mechanism <- factor(
    smooth_irmse_plot_data$missing_mechanism,
    levels = missingness_registry$missing_mechanism
  )
  smooth_irmse_plot_data$model <- factor(
    smooth_irmse_plot_data$model,
    levels = c("gamlss2", "gamlss.longitudinal")
  )
  plot_smooth_irmse <- ggplot(
    smooth_irmse_plot_data,
    aes(x = target_missing_pct, y = mean_irmse, colour = model, shape = model)
  ) +
    geom_line(linewidth = 0.75) +
    geom_point(size = 2) +
    facet_grid(parameter ~ missing_mechanism, scales = "free_y") +
    scale_x_continuous(breaks = seq(0, 50, 10)) +
    labs(x = "Target missingness (%)", y = "Smooth IRMSE", colour = NULL, shape = NULL) +
    theme_bw(base_size = 10) +
    theme(
      legend.position = "top",
      panel.spacing = grid::unit(0.8, "lines")
    )
  ggsave(file.path(out_dir, "smooth_irmse_by_missingness.png"), plot_smooth_irmse, width = 8.5, height = 5.2, dpi = 180)

  selected_scenarios <- scenario_meta[
    scenario_meta$target_missing_rate %in% c(0, 0.3, 0.5),
    ,
    drop = FALSE
  ]
  selected_smooth <- merge(
    smooth_summary,
    selected_scenarios,
    by = "scenario",
    all = FALSE,
    sort = FALSE
  )
  selected_smooth <- merge(
    selected_smooth,
    fit_attempt_summary[, c(
      "scenario", "model", "attempted", "retained", "failed",
      "failure_inclusive_retention_rate", "retention_rate_mcse",
      "retention_rate_conf_low", "retention_rate_conf_high"
    )],
    by = c("scenario", "model"),
    all.x = TRUE,
    sort = FALSE
  )
  selected_smooth$missing_label <- paste0(round(100 * selected_smooth$target_missing_rate), "%")
  selected_smooth$missing_label <- factor(selected_smooth$missing_label, levels = c("0%", "30%", "50%"))
  selected_smooth$model <- factor(selected_smooth$model, levels = c("gamlss2", "gamlss.longitudinal"))
  write.csv(selected_smooth, file.path(out_dir, "smooth_selected_plot_data.csv"), row.names = FALSE)
  plot_smooth_selected <- ggplot(selected_smooth, aes(x = s1)) +
    geom_line(aes(y = smooth_true), colour = "black", linewidth = 0.75) +
    geom_line(aes(y = smooth_median, colour = model), linewidth = 0.75, linetype = "dashed") +
    facet_grid(parameter + missing_mechanism ~ missing_label, scales = "free_y") +
    labs(x = "s1", y = "Centered smooth", colour = NULL) +
    theme_bw(base_size = 10) +
    theme(
      legend.position = "top",
      panel.spacing = grid::unit(0.65, "lines")
    )
  ggsave(file.path(out_dir, "smooth_selected_recovery_curves.png"), plot_smooth_selected, width = 8.5, height = 7, dpi = 180)
}

fixed_plot_summary <- fixed_summary[fixed_summary$term != "intercept", , drop = FALSE]
fixed_plot_summary$missing_mechanism <- factor(
  fixed_plot_summary$missing_mechanism,
  levels = missingness_registry$missing_mechanism
)
fixed_plot_summary$target_missing_pct <- 100 * fixed_plot_summary$target_missing_rate
fixed_plot_summary$model <- factor(
  fixed_plot_summary$model,
  levels = c("gamlss2", "gamlss.longitudinal")
)
fixed_term_summary <- do.call(rbind, lapply(
  split(
    fixed_plot_summary,
    list(
      fixed_plot_summary$missing_mechanism,
      fixed_plot_summary$target_missing_rate,
      fixed_plot_summary$model,
      fixed_plot_summary$parameter
    ),
    drop = TRUE
  ),
  function(df) {
    data.frame(
      missing_mechanism = df$missing_mechanism[1],
      target_missing_rate = df$target_missing_rate[1],
      observed_missing_rate = mean(df$observed_missing_rate, na.rm = TRUE),
      complete_adjacent_pair_rate = mean(df$complete_adjacent_pair_rate, na.rm = TRUE),
      model = df$model[1],
      parameter = df$parameter[1],
      attempted = max(df$attempted, na.rm = TRUE),
      retained = min(df$n_successful_fits, na.rm = TRUE),
      failed = max(df$attempted - df$n_successful_fits, na.rm = TRUE),
      mean_rmse = mean(df$rmse, na.rm = TRUE),
      mean_abs_bias = mean(abs(df$bias), na.rm = TRUE),
      mean_coverage_95 = mean(df$coverage_95, na.rm = TRUE),
      inference_status = if (df$model[[1L]] == "gamlss.longitudinal" &&
        df$missingness_pattern[[1L]] == "intermittent")
        "not_applicable_segmented_model_hessian" else "conditional_on_retained_fits",
      stringsAsFactors = FALSE
    )
  }
))
fixed_term_summary <- fixed_term_summary[is.finite(fixed_term_summary$mean_rmse), , drop = FALSE]
fixed_term_summary$target_missing_pct <- 100 * fixed_term_summary$target_missing_rate
fixed_term_summary$model <- factor(
  fixed_term_summary$model,
  levels = c("gamlss2", "gamlss.longitudinal")
)
write.csv(fixed_term_summary, file.path(out_dir, "fixed_term_summary_by_missingness.csv"), row.names = FALSE)

fixed_rmse_plot <- ggplot(
  fixed_term_summary[fixed_term_summary$parameter %in% params_margin, , drop = FALSE],
  aes(x = target_missing_pct, y = mean_rmse, colour = model, shape = model, group = model)
) +
  geom_line(linewidth = 0.7) +
  geom_point(size = 1.9) +
  facet_grid(parameter ~ missing_mechanism, scales = "free_y") +
  scale_x_continuous(breaks = seq(0, 50, 10)) +
  labs(x = "Target missingness (%)", y = "Fixed-effect RMSE", colour = NULL, shape = NULL) +
  theme_bw(base_size = 10) +
  theme(
    legend.position = "top",
    panel.spacing = grid::unit(0.75, "lines")
  )
ggsave(file.path(out_dir, "fixed_margin_rmse_by_missingness.png"), fixed_rmse_plot, width = 8.5, height = 5.5, dpi = 180)

coverage_plot_payload <- jss_missing_coverage_plot_payload(fixed_term_summary, params_margin)
warning_audit <- jss_missing_warning_audit(warning_events_all, runs_all, coverage_plot_payload)
jss_missing_assert_warning_policy(warning_events_all)
write.csv(warning_audit, file.path(out_dir, "missingness_warning_audit.csv"), row.names = FALSE)

fixed_coverage_plot <- ggplot(
  coverage_plot_payload$data,
  aes(x = target_missing_pct, y = mean_coverage_95, colour = model, shape = model, group = model)
) +
  geom_hline(yintercept = 0.95, colour = "gray45", linetype = "dotted") +
  geom_line(linewidth = 0.7) +
  geom_point(size = 1.9) +
  facet_grid(parameter ~ missing_mechanism) +
  scale_x_continuous(breaks = seq(0, 50, 10)) +
  coord_cartesian(ylim = c(0.75, 1.02)) +
  labs(x = "Target missingness (%)", y = "95% CI coverage", colour = NULL, shape = NULL) +
  theme_bw(base_size = 10) +
  theme(
    legend.position = "top",
    panel.spacing = grid::unit(0.75, "lines")
  )
ggsave(file.path(out_dir, "fixed_margin_coverage_by_missingness.png"), fixed_coverage_plot, width = 8.5, height = 5.5, dpi = 180)

copula_rmse_plot <- ggplot(
  fixed_term_summary[
    fixed_term_summary$model == "gamlss.longitudinal" &
      fixed_term_summary$parameter %in% params_copula,
    ,
    drop = FALSE
  ],
  aes(x = target_missing_pct, y = mean_rmse, colour = parameter, shape = parameter, group = parameter)
) +
  geom_line(linewidth = 0.75) +
  geom_point(size = 2) +
  facet_wrap(~ missing_mechanism, nrow = 1) +
  scale_x_continuous(breaks = seq(0, 50, 10)) +
  labs(x = "Target missingness (%)", y = "Copula fixed-effect RMSE", colour = NULL, shape = NULL) +
  theme_bw(base_size = 10) +
  theme(legend.position = "top")
ggsave(file.path(out_dir, "fixed_copula_rmse_by_missingness.png"), copula_rmse_plot, width = 8.5, height = 3.8, dpi = 180)

cat("\nSimulation comparison complete. Outputs written to:", out_dir, "\n")
