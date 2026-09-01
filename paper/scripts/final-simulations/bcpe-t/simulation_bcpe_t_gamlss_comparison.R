#!/usr/bin/env Rscript

# Simulation comparison for gamlss.longitudinal vs standard gamlss.
# Defaults are deliberately small so the script can be smoke-tested quickly.

suppressPackageStartupMessages({
  library(parallel)
  library(VineCopula)
  library(gamlss)
  library(gamlss2)
  library(gamlss.dist)
  library(ggplot2)
})

suppressPackageStartupMessages(library(gamlss.longitudinal))
list2env(as.list(getNamespace("gamlss.longitudinal"), all.names = TRUE), envir = .GlobalEnv)

set.seed(20260513)

out_dir <- Sys.getenv("OUT_DIR", unset = file.path("results", "bcpe_t_gamlss_comparison"))
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)
project_root <- normalizePath(getwd(), winslash = "/", mustWork = TRUE)
task_result_dir <- file.path(out_dir, "rep_results")
worker_log_dir <- file.path(out_dir, "worker_logs")
dir.create(task_result_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(worker_log_dir, recursive = TRUE, showWarnings = FALSE)

n_fits <- as.integer(Sys.getenv("N_FITS", unset = "1"))
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

scenarios <- data.frame(
  n = c(500),
  d = c(4),
  stringsAsFactors = FALSE
)

params_margin <- c("mu", "sigma", "nu", "tau")
params_copula <- c("theta", "zeta")
params_all <- c(params_margin, params_copula)
fixed_terms <- c("intercept", "x1", "x2", "t")
theta_effect_mode <- tolower(Sys.getenv("THETA_EFFECT_MODE", unset = "smooth"))
if (!theta_effect_mode %in% c("smooth", "binary")) {
  stop("THETA_EFFECT_MODE must be either 'smooth' or 'binary'.", call. = FALSE)
}
theta_binary_effect <- as.numeric(Sys.getenv("THETA_BINARY_EFFECT", unset = "1.2"))
theta_binary_prob <- as.numeric(Sys.getenv("THETA_BINARY_PROB", unset = "0.5"))
smooth_params_longitudinal <- if (theta_effect_mode == "smooth") c("mu", "sigma", "theta") else c("mu", "sigma")
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
variogram_p_values_env <- Sys.getenv("VARIOGRAM_P_VALUES", unset = "")
variogram_p_values <- if (nzchar(variogram_p_values_env)) {
  as.numeric(trimws(strsplit(variogram_p_values_env, "[|,;[:space:]]+")[[1]]))
} else {
  variogram_p
}
variogram_p_values <- unique(variogram_p_values[is.finite(variogram_p_values)])
if (length(variogram_p_values) == 0L) {
  variogram_p_values <- variogram_p
}
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

linpred <- function(beta, x1, x2, t) {
  beta[["intercept"]] + beta[["x1"]] * x1 + beta[["x2"]] * x2 + beta[["t"]] * t
}

theta_effect_truth <- function(dat) {
  if (theta_effect_mode == "binary") {
    if (!"theta_group" %in% names(dat)) {
      stop("Binary theta mode requires a theta_group column.", call. = FALSE)
    }
    return(theta_binary_effect * as.numeric(dat$theta_group == "high"))
  }
  smooth_truth$theta(dat$s1)
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
      converged = NA,
      hit_outer_limit = NA,
      hit_max_stall = NA,
      hit_raw_loglik_deterioration = NA,
      stop_reason = NA_character_,
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
    grad_inf = if (!is.null(conv$grad_inf)) as.numeric(conv$grad_inf) else NA_real_,
    step_l2 = if (!is.null(conv$step_l2)) as.numeric(conv$step_l2) else NA_real_,
    best_raw_loglik = if (!is.null(conv$best_raw_loglik)) as.numeric(conv$best_raw_loglik) else NA_real_,
    best_raw_loglik_iteration = if (!is.null(conv$best_raw_loglik_iteration)) as.integer(conv$best_raw_loglik_iteration) else NA_integer_,
    raw_loglik_drop_from_best = if (!is.null(conv$raw_loglik_drop_from_best)) as.numeric(conv$raw_loglik_drop_from_best) else NA_real_,
    raw_loglik_drop_tol = if (!is.null(conv$raw_loglik_drop_tol)) as.numeric(conv$raw_loglik_drop_tol) else NA_real_,
    outer_iterations = if (!is.null(conv$outer_iterations)) as.integer(conv$outer_iterations) else NA_integer_,
    outer_log_lik_change = if (!is.null(conv$outer_log_lik_change)) as.numeric(conv$outer_log_lik_change) else NA_real_,
    outer_stop_crit = if (!is.null(conv$outer_stop_crit)) as.numeric(conv$outer_stop_crit) else NA_real_
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
    sprintf("scenario%d_n%d_d%d_rep%d.rds", task$scenario_id, task$n, task$d, task$rep)
  )
}

run_one_rep_and_save <- function(task) {
  cat(sprintf(
    "[pid %s] starting scenario %d n=%d d=%d rep=%d\n",
    Sys.getpid(), task$scenario_id, task$n, task$d, task$rep
  ))
  flush.console()
  result <- run_one_rep(task)
  final_path <- task_result_path(task)
  temporary_path <- paste0(final_path, ".", Sys.getpid(), ".tmp")
  saveRDS(result, temporary_path)
  if (file.exists(final_path)) file.remove(final_path)
  if (!file.rename(temporary_path, final_path)) stop("Could not atomically install checkpoint: ", final_path)
  cat(sprintf(
    "[pid %s] finished scenario %d n=%d d=%d rep=%d -> %s\n",
    Sys.getpid(), task$scenario_id, task$n, task$d, task$rep, task_result_path(task)
  ))
  flush.console()
  result
}

center_curve <- function(y) {
  y - mean(y, na.rm = TRUE)
}

calc_smooth_mean <- function(data_used, parameter) {
  data_sub <- data_used
  if (parameter %in% params_copula) {
    data_sub <- data_sub[data_sub$time < max(data_sub$time), , drop = FALSE]
  }
  mean(smooth_truth[[parameter]](data_sub$s1), na.rm = TRUE)
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
  if (theta_effect_mode == "binary") {
    subjects$theta_group <- factor(
      ifelse(runif(n) < theta_binary_prob, "high", "low"),
      levels = c("low", "high")
    )
  }

  dat <- merge(subjects, data.frame(time_index = seq_len(d)), by = NULL, all = TRUE)
  dat <- dat[order(dat$id, dat$time_index), ]
  dat$t <- if (d > 1) (dat$time_index - 1) / (d - 1) else 0
  dat$time <- dat$time_index

  simulate_response_given_covariates(dat, seed = seed)
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
  eta_theta <- linpred(true_beta$theta, dat$x1, dat$x2, dat$t) + theta_effect_truth(dat)
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

  keep_cols <- c("id", "time", "response", "s1", "x1", "x2", "t")
  if ("theta_group" %in% names(dat)) {
    keep_cols <- c(keep_cols, "theta_group")
  }
  dat[, keep_cols]
}

fit_longitudinal_model <- function(dat) {
  theta_formula <- if (theta_effect_mode == "binary") {
    "~ x1 + x2 + t + theta_group"
  } else {
    sprintf("~ x1 + x2 + t + s(s1, k = %d)", smooth_k)
  }

  gamlss.longitudinal(
    dataset = dat,
    margin_dist = gamlss.dist::BCPE(mu.link = "log"),
    copula_dist = "t",
    time_var = "time",
    subject_var = "id",
    mu.formula = sprintf("response ~ x1 + x2 + t + s(s1, k = %d)", smooth_k),
    sigma.formula = sprintf("~ x1 + x2 + t + s(s1, k = %d)", smooth_k),
    nu.formula = "~ x1 + x2 + t",
    tau.formula = "~ x1 + x2 + t",
    theta.formula = theta_formula,
    zeta.formula = "~ x1 + x2 + t",
    verbose = verbose_level,
    compute_vcov = isTRUE(compute_se),
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
      fitted_smooth_parts <- Map(
        function(B, beta_s) as.numeric(B %*% beta_s),
        fit_obj$model_matrix$s[[p]],
        fit_obj$par_s[[p]]
      )
      fitted_smooth_mean <- mean(rowSums(do.call(cbind, fitted_smooth_parts)), na.rm = TRUE)
    }
    for (tm in fixed_terms) {
      true_value <- true_beta[[p]][[tm]]
      estimate <- extract_one_longitudinal_term(par_vec, p, tm)
      if (tm == "intercept") {
        true_value <- true_value + smooth_mean
        estimate <- estimate + fitted_smooth_mean
      }
      rows[[length(rows) + 1]] <- data.frame(
        model = "gamlss.longitudinal",
        parameter = p,
        term = tm,
        estimate = estimate,
        std_error = extract_one_longitudinal_term(se_vec, p, tm),
        true_value = true_value,
        intercept_includes_fitted_smooth_mean = tm == "intercept",
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
      std_error <- if (p %in% params_margin) {
        extract_one_gamlss_term(se_by_param[[p]], tm)
      } else {
        NA_real_
      }
      rows[[length(rows) + 1]] <- data.frame(
        model = "gamlss2",
        parameter = p,
        term = tm,
        estimate = estimate,
        std_error = std_error,
        true_value = if (p %in% params_margin) true_value else NA_real_,
        intercept_includes_fitted_smooth_mean = FALSE,
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

simulate_array_from_columns <- function(sim_df, dat, subject_col = "id") {
  dat <- dat[order(dat[[subject_col]], dat$time), , drop = FALSE]
  sim_mat <- as.matrix(sim_df)
  n <- length(unique(dat[[subject_col]]))
  d <- length(unique(dat$time))
  out <- array(NA_real_, dim = c(ncol(sim_mat), n, d))
  for (s in seq_len(ncol(sim_mat))) {
    out[s, , ] <- matrix(sim_mat[, s], nrow = n, ncol = d, byrow = TRUE)
  }
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

simulate_predictive_longitudinal <- function(fit_obj, dat, nsim) {
  dat <- dat[order(dat$id, dat$time), , drop = FALSE]
  sim_df <- stats::simulate(fit_obj, nsim = nsim, newdata = dat)
  simulate_array_from_columns(sim_df, dat, subject_col = "id")
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
  sim_array <- simulate_predictive_longitudinal(fit_obj, test_dat, nsim = nsim)
  y_obs <- response_matrix(test_dat)
  base <- data.frame(
    scenario = scenario_label,
    n = n_val,
    d = d_val,
    rep = rep_id,
    model = "gamlss.longitudinal",
    test_log_score_joint = as.numeric(lik$log_lik["joint"]),
    test_log_score_marginal = as.numeric(lik$log_lik["marginal"]),
    test_log_score_copula = as.numeric(lik$log_lik["copula"]),
    test_log_score_per_obs = as.numeric(lik$log_lik["joint"]) / nrow(test_dat),
    predictive_nsim = nsim,
    stringsAsFactors = FALSE
  )
  do.call(rbind, lapply(variogram_p_values, function(p_value) {
    cbind(
      base,
      variogram_score = variogram_score(y_obs, sim_array, p = p_value),
      variogram_p = p_value
    )
  }))
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
  base <- data.frame(
    scenario = scenario_label,
    n = n_val,
    d = d_val,
    rep = rep_id,
    model = "gamlss2",
    test_log_score_joint = joint_log_score,
    test_log_score_marginal = joint_log_score,
    test_log_score_copula = 0,
    test_log_score_per_obs = joint_log_score / nrow(test_dat),
    predictive_nsim = nsim,
    stringsAsFactors = FALSE
  )
  do.call(rbind, lapply(variogram_p_values, function(p_value) {
    cbind(
      base,
      variogram_score = variogram_score(y_obs, sim_array, p = p_value),
      variogram_p = p_value
    )
  }))
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
  n_val <- task$n
  d_val <- task$d
  rep_id <- task$rep
  scenario_label <- sprintf("n%d_d%d", n_val, d_val)
  seed <- 100000 + 10000 * task$scenario_id + rep_id
  s1_grid <- seq(0, 1, length.out = 101)

  dat <- simulate_dataset(n = n_val, d = d_val, seed = seed)

  run_rows <- list()
  fixed_rows <- list()
  smooth_rows <- list()
  joint_rows <- list()
  predictive_rows <- list()

  t0 <- Sys.time()
  err_long <- NA_character_
  fit_long <- tryCatch(
    fit_longitudinal_model(dat),
    error = function(e) {
      err_long <<- conditionMessage(e)
      NULL
    }
  )
  elapsed_long <- as.numeric(difftime(Sys.time(), t0, units = "secs"))

  t0 <- Sys.time()
  err_gamlss2 <- NA_character_
  fit_gamlss2 <- tryCatch(
    fit_gamlss2_model(dat),
    error = function(e) {
      err_gamlss2 <<- conditionMessage(e)
      NULL
    }
  )
  elapsed_gamlss2 <- as.numeric(difftime(Sys.time(), t0, units = "secs"))

  conv_long <- if (!is.null(fit_long)) extract_convergence_info(fit_long) else extract_convergence_info(NULL)

  run_rows[[1]] <- data.frame(
    scenario = scenario_label, n = n_val, d = d_val, rep = rep_id,
    model = "gamlss.longitudinal", success = !is.null(fit_long),
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
    model = "gamlss2", success = !is.null(fit_gamlss2),
    logLik = if (!is.null(fit_gamlss2)) extract_fit_metric(fit_gamlss2, "logLik") else NA_real_,
    df = if (!is.null(fit_gamlss2)) extract_fit_metric(fit_gamlss2, "df") else NA_real_,
    converged = NA,
    hit_outer_limit = NA,
    hit_max_stall = NA,
    hit_raw_loglik_deterioration = NA,
    stop_reason = NA_character_,
    grad_inf = NA_real_,
    step_l2 = NA_real_,
    best_raw_loglik = NA_real_,
    best_raw_loglik_iteration = NA_integer_,
    raw_loglik_drop_from_best = NA_real_,
    raw_loglik_drop_tol = NA_real_,
    outer_iterations = NA_integer_,
    outer_log_lik_change = NA_real_,
    outer_stop_crit = NA_real_,
    elapsed_sec = elapsed_gamlss2, error = if (is.null(fit_gamlss2)) err_gamlss2 else NA_character_,
    stringsAsFactors = FALSE
  )

  if (!is.null(fit_long)) {
    fixed_long <- extract_fixed_estimates_longitudinal(fit_long, dat)
    fixed_long$scenario <- scenario_label
    fixed_long$n <- n_val
    fixed_long$d <- d_val
    fixed_long$rep <- rep_id
    fixed_rows[[length(fixed_rows) + 1]] <- fixed_long

    for (p in smooth_params_longitudinal) {
      true_curve <- center_curve(smooth_truth[[p]](s1_grid))
      smooth_rows[[length(smooth_rows) + 1]] <- data.frame(
        scenario = scenario_label, n = n_val, d = d_val, rep = rep_id,
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

  if (!is.null(fit_gamlss2)) {
    fixed_gamlss2 <- extract_fixed_estimates_gamlss2(fit_gamlss2, dat)
    fixed_gamlss2$scenario <- scenario_label
    fixed_gamlss2$n <- n_val
    fixed_gamlss2$d <- d_val
    fixed_gamlss2$rep <- rep_id
    fixed_rows[[length(fixed_rows) + 1]] <- fixed_gamlss2

    for (p in smooth_params_gamlss) {
      true_curve <- center_curve(smooth_truth[[p]](s1_grid))
      smooth_rows[[length(smooth_rows) + 1]] <- data.frame(
        scenario = scenario_label, n = n_val, d = d_val, rep = rep_id,
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

  out <- list(
    fixed = bind_non_null(fixed_rows),
    smooth = bind_non_null(smooth_rows),
    joint = bind_non_null(joint_rows),
    predictive = bind_non_null(predictive_rows),
    runs = do.call(rbind, run_rows)
  )
  if (isTRUE(save_fits)) {
    out$fit_longitudinal <- fit_long
    out$fit_gamlss2 <- fit_gamlss2
    out$data <- dat
  }
  out
}

tasks <- do.call(
  rbind,
  lapply(seq_len(nrow(scenarios)), function(i) {
    data.frame(
      scenario_id = i,
      n = scenarios$n[i],
      d = scenarios$d[i],
      rep = rep_ids,
      stringsAsFactors = FALSE
    )
  })
)
tasks <- split(tasks, seq_len(nrow(tasks)))
checkpoint_results <- lapply(tasks, function(task) {
  path <- task_result_path(task)
  if (!file.exists(path)) return(NULL)
  tryCatch(readRDS(path), error = function(e) NULL)
})
completed <- !vapply(checkpoint_results, is.null, logical(1))
results <- checkpoint_results[completed]
tasks <- tasks[!completed]

cat(sprintf(
  "Running %d fit replicate(s) across %d scenario(s): %d checkpoint(s) complete, %d pending, using up to %d core(s).\n",
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

  on.exit(parallel::stopCluster(cl), add = TRUE)
  parallel::clusterCall(cl, setwd, project_root)
  worker_status <- parallel::clusterEvalQ(cl, {
    suppressPackageStartupMessages({
      library(VineCopula)
      library(gamlss)
      library(gamlss2)
      library(gamlss.dist)
      library(ggplot2)
    })
    suppressPackageStartupMessages(library(gamlss.longitudinal))
    list2env(as.list(getNamespace("gamlss.longitudinal"), all.names = TRUE), envir = .GlobalEnv)
    list(pid = Sys.getpid(), wd = getwd())
  })
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
    quit(save = "no", status = 0)
  }

  new_results <- parallel::parLapplyLB(cl, tasks, run_one_rep_and_save)
}
results <- c(results, new_results)

fixed_all <- bind_non_null(lapply(results, `[[`, "fixed"))
smooth_all <- bind_non_null(lapply(results, `[[`, "smooth"))
joint_all <- bind_non_null(lapply(results, `[[`, "joint"))
predictive_all <- bind_non_null(lapply(results, `[[`, "predictive"))
runs_all <- bind_non_null(lapply(results, `[[`, "runs"))

saveRDS(
  list(
    fixed = fixed_all,
    smooth = smooth_all,
    joint = joint_all,
    predictive = predictive_all,
    runs = runs_all,
    truth = true_beta,
    theta_effect_mode = theta_effect_mode,
    theta_binary_effect = theta_binary_effect,
    theta_binary_prob = theta_binary_prob
  ),
  file.path(out_dir, "all_results.rds")
)
write.csv(runs_all, file.path(out_dir, "fit_run_log.csv"), row.names = FALSE)
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
    covered <- is.finite(df$true_value) & is.finite(ci_low) & is.finite(ci_high) & df$true_value >= ci_low & df$true_value <= ci_high
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
      coverage_95 = mean(covered, na.rm = TRUE),
      q05 = safe_quantile(df$estimate, 0.05),
      q95 = safe_quantile(df$estimate, 0.95),
      n_successful_fits = sum(is.finite(df$estimate)),
      n_se = sum(is.finite(df$std_error)),
      stringsAsFactors = FALSE
    )
  }
))

fit_metric_summary <- do.call(rbind, lapply(
  split(runs_all, list(runs_all$scenario, runs_all$model), drop = TRUE),
  function(df) {
    data.frame(
      scenario = df$scenario[1],
      model = df$model[1],
      mean_logLik = mean(df$logLik, na.rm = TRUE),
      mean_df = mean(df$df, na.rm = TRUE),
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

write.csv(fixed_summary, file.path(out_dir, "fixed_effects_bias_rmse_table.csv"), row.names = FALSE)

fixed_compact <- reshape(
  fixed_summary[, c("scenario", "n", "d", "parameter", "term", "model", "bias", "rmse", "mean_std_error", "sd_estimate", "se_to_empirical_sd", "coverage_95", "mean_logLik", "mean_df")],
  idvar = c("scenario", "n", "d", "parameter", "term"),
  timevar = "model",
  direction = "wide"
)
write.csv(fixed_compact, file.path(out_dir, "fixed_effects_bias_rmse_compact.csv"), row.names = FALSE)

se_calibration <- fixed_summary[fixed_summary$term != "intercept", c(
  "scenario", "model", "n", "d", "parameter", "term",
  "mean_std_error", "sd_estimate", "se_to_empirical_sd", "coverage_95", "n_se"
)]
write.csv(se_calibration, file.path(out_dir, "fixed_effects_se_calibration.csv"), row.names = FALSE)

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
    split(predictive_all, list(predictive_all$scenario, predictive_all$model, predictive_all$variogram_p), drop = TRUE),
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

  plot_smooth <- ggplot(smooth_summary, aes(x = s1)) +
    geom_ribbon(aes(ymin = smooth_q05, ymax = smooth_q95), fill = "gray82", alpha = 0.85) +
    geom_line(aes(y = smooth_median), color = "gray30", linewidth = 0.65, linetype = "dashed") +
    geom_line(aes(y = smooth_true), color = "#B3262E", linewidth = 0.85) +
    facet_grid(model ~ parameter + scenario, scales = "free_y") +
    labs(x = "s1", y = "Centered smooth contribution") +
    theme_bw(base_size = 10) +
    theme(panel.spacing = grid::unit(0.7, "lines"))

  ggsave(file.path(out_dir, "smooth_recovery_plot.png"), plot_smooth, width = 13, height = 6.5, dpi = 180)
}

plot_fixed <- ggplot(fixed_summary, aes(x = term, y = mean_estimate)) +
  geom_errorbar(aes(ymin = q05, ymax = q95), width = 0.18, color = "gray55") +
  geom_point(size = 1.6, color = "black") +
  geom_point(aes(y = true_value), shape = 4, size = 2.7, stroke = 1.0, color = "#B3262E") +
  facet_grid(model + scenario ~ parameter, scales = "free_y") +
  labs(x = NULL, y = "Estimate on eta scale") +
  theme_bw(base_size = 10) +
  theme(
    axis.text.x = element_text(angle = 30, hjust = 1),
    panel.spacing = grid::unit(0.7, "lines")
  )

ggsave(file.path(out_dir, "fixed_effect_recovery_plot.png"), plot_fixed, width = 14, height = 7, dpi = 180)

cat("\nSimulation comparison complete. Outputs written to:", out_dir, "\n")
