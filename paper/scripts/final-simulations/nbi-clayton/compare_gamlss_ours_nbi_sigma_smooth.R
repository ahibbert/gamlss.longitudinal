#!/usr/bin/env Rscript

suppressPackageStartupMessages({
library(VineCopula)
library(gamlss)
library(gamlss.dist)
library(ggplot2)
})

project_root <- normalizePath(Sys.getenv("GAMLSS_LONGITUDINAL_SOURCE_ROOT", unset = getwd()), winslash = "/", mustWork = TRUE)
package_source_path <- project_root
package_source_identity <- "."
package_source_files <- c(file.path(project_root, "DESCRIPTION"), file.path(project_root, "NAMESPACE"), sort(list.files(file.path(project_root, "R"), pattern = "[.]R$", full.names = TRUE)))
if (!all(file.exists(package_source_files)) || !requireNamespace("pkgload", quietly = TRUE) || !requireNamespace("digest", quietly = TRUE)) stop("The checked-out package source, pkgload, and digest are required for authoritative recovery runs.", call. = FALSE)
source_file_hashes <- unname(vapply(package_source_files, digest::digest, character(1), file = TRUE, algo = "sha256", serialize = FALSE))
source_relative <- substring(normalizePath(package_source_files, winslash = "/", mustWork = TRUE), nchar(project_root) + 2L)
package_source_sha256 <- digest::digest(paste(source_relative, source_file_hashes, sep = "\t", collapse = "\n"), algo = "sha256", serialize = FALSE)
if ("package:gamlss.longitudinal" %in% search()) try(detach("package:gamlss.longitudinal", unload = TRUE, character.only = TRUE), silent = TRUE)
if ("gamlss.longitudinal" %in% loadedNamespaces()) try(unloadNamespace("gamlss.longitudinal"), silent = TRUE)
suppressPackageStartupMessages(pkgload::load_all(package_source_path, export_all = TRUE, helpers = FALSE, quiet = TRUE))
list2env(as.list(getNamespace("gamlss.longitudinal"), all.names = TRUE), envir = .GlobalEnv)

`%||%` <- function(x, y) if (is.null(x)) y else x

n_subject <- as.integer(Sys.getenv("NBI_COMPARE_N", unset = "500"))
reps <- as.integer(Sys.getenv("NBI_COMPARE_REPS", unset = "1"))
base_seed <- as.integer(Sys.getenv("NBI_COMPARE_SEED", unset = "20260529"))
runner_contract_version <- "nbi-clayton-main-recovery-2026-09-02.5"
phase1_contract_version <- "likelihood-jss001-2026-09-01|inference-2026.1|capability-2026.2"
nbi_copula_code <- "C"
nbi_copula_label <- "Clayton"
times <- seq(0, 1, length.out = 4)
sigma_signal_multiplier <- as.numeric(Sys.getenv("NBI_COMPARE_SIGMA_SIGNAL_MULTIPLIER", unset = "2"))
max_elapsed_sec <- as.numeric(Sys.getenv("NBI_COMPARE_MAX_ELAPSED", unset = "1800"))
max_outer_iter <- as.integer(Sys.getenv("NBI_COMPARE_MAX_OUTER", unset = "80"))
max_inner_iter <- as.integer(Sys.getenv("NBI_COMPARE_MAX_INNER", unset = "40"))
start_step_size <- as.numeric(Sys.getenv("NBI_COMPARE_START_STEP_SIZE", unset = "0.5"))
step_adjustment_env <- Sys.getenv("NBI_COMPARE_STEP_ADJUSTMENT", unset = "NA")
step_adjustment <- suppressWarnings(as.numeric(step_adjustment_env))
lambda_start <- as.numeric(Sys.getenv("NBI_COMPARE_LAMBDA_START", unset = "10"))
rs_update_lambda <- as.logical(Sys.getenv("NBI_COMPARE_RS_UPDATE_LAMBDA", unset = "FALSE"))
warm_start_joint_iter <- as.integer(Sys.getenv("NBI_COMPARE_WARM_START_JOINT_ITER", unset = "3"))
engine_env <- Sys.getenv("NBI_COMPARE_ENGINES", unset = "gamlss|ours_rs_separate|ours_rs_joint")
engines <- trimws(strsplit(engine_env, "\\|")[[1]])
engines <- engines[nzchar(engines)]
resume <- as.logical(Sys.getenv("NBI_COMPARE_RESUME", unset = "TRUE"))
max_attempts_per_fit <- as.integer(Sys.getenv("NBI_COMPARE_MAX_ATTEMPTS_PER_FIT", unset = "1"))
if (!is.finite(max_attempts_per_fit) || max_attempts_per_fit < 1L) {
  stop("NBI_COMPARE_MAX_ATTEMPTS_PER_FIT must be a positive integer.", call. = FALSE)
}
compute_se <- as.logical(Sys.getenv("NBI_COMPARE_COMPUTE_SE", unset = "TRUE"))
vcov_method_longitudinal <- Sys.getenv("NBI_COMPARE_VCOV_METHOD", unset = "analytical")
compute_predictive_scores <- as.logical(Sys.getenv("NBI_COMPARE_COMPUTE_PREDICTIVE", unset = "TRUE"))
predictive_nsim <- as.integer(Sys.getenv("NBI_COMPARE_PREDICTIVE_NSIM", unset = "100"))
variogram_p <- as.numeric(Sys.getenv("NBI_COMPARE_VARIOGRAM_P", unset = "0.5"))
variogram_p_values_env <- Sys.getenv("NBI_COMPARE_VARIOGRAM_P_VALUES", unset = "")
variogram_p_values <- if (nzchar(variogram_p_values_env)) {
  as.numeric(trimws(strsplit(variogram_p_values_env, "[|,;[:space:]]+")[[1]]))
} else {
  variogram_p
}
variogram_p_values <- unique(variogram_p_values[is.finite(variogram_p_values)])
if (length(variogram_p_values) == 0L) {
  variogram_p_values <- variogram_p
}
save_fits <- as.logical(Sys.getenv("NBI_COMPARE_SAVE_FITS", unset = "FALSE"))
theta_intercept <- as.numeric(Sys.getenv("NBI_COMPARE_THETA_INTERCEPT", unset = "0.70"))
theta_time_coef <- as.numeric(Sys.getenv("NBI_COMPARE_THETA_TIME_COEF", unset = "0.20"))
runner_git_sha <- tryCatch(system2("git", c("rev-parse", "HEAD"), stdout = TRUE, stderr = FALSE)[[1]], error = function(e) NA_character_)
runner_git_state <- tryCatch(if (length(system2("git", c("status", "--porcelain"), stdout = TRUE, stderr = FALSE))) "dirty" else "clean", error = function(e) "unknown")
runner_package_version <- as.character(read.dcf(file.path(package_source_path, "DESCRIPTION"))[1, "Version"])
runner_script_path <- normalizePath(file.path(project_root, "paper", "scripts", "final-simulations", "nbi-clayton", "compare_gamlss_ours_nbi_sigma_smooth.R"), winslash = "/", mustWork = TRUE)
runner_sha256 <- unname(digest::digest(runner_script_path, file = TRUE, algo = "sha256", serialize = FALSE))
nbi_rscript_path <- normalizePath(file.path(R.home("bin"), if (.Platform$OS.type == "windows") "Rscript.exe" else "Rscript"), winslash = "/", mustWork = TRUE)
nbi_rscript_sha256 <- unname(digest::digest(nbi_rscript_path, file = TRUE, algo = "sha256", serialize = FALSE))
runner_settings_signature <- paste(
  runner_contract_version, phase1_contract_version, n_subject, paste(times, collapse = ","), reps, base_seed,
  sigma_signal_multiplier, paste(engines, collapse = ","), max_elapsed_sec, max_outer_iter, max_inner_iter,
  start_step_size, step_adjustment_env, lambda_start, rs_update_lambda, warm_start_joint_iter,
  compute_se, vcov_method_longitudinal, compute_predictive_scores, predictive_nsim,
  paste(variogram_p_values, collapse = ","), theta_intercept, theta_time_coef,
  max_attempts_per_fit, 1L, "sequential", nbi_rscript_path, nbi_rscript_sha256, R.version.string,
  sep = "|"
)
runner_settings_sha256 <- digest::digest(runner_settings_signature, algo = "sha256", serialize = FALSE)

out_dir <- Sys.getenv(
  "NBI_COMPARE_OUT_DIR",
  unset = file.path(
    "results", "del_clayton_rectangle_diagnostics",
    sprintf(
      "nbi_sigma_compare_n%d_rep%d_signal%s",
      n_subject,
      reps,
      gsub("[^0-9A-Za-z]+", "p", as.character(sigma_signal_multiplier))
    )
  )
)
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)
if (identical(Sys.getenv("NBI_SOURCE_IDENTITY_AUDIT_ONLY", unset = "0"), "1")) {
  write.csv(data.frame(package_source_path = package_source_identity, package_version = runner_package_version,
    package_source_sha256 = package_source_sha256, runner_sha256 = runner_sha256,
    runner_settings_signature = runner_settings_signature, runner_settings_sha256 = runner_settings_sha256,
    loaded_namespace_path = normalizePath(getNamespaceInfo(asNamespace("gamlss.longitudinal"), "path"), winslash = "/", mustWork = TRUE)),
    file.path(out_dir, "source_identity_audit.csv"), row.names = FALSE)
  quit(save = "no", status = 0)
}
fit_dir <- file.path(out_dir, "fits")
if (isTRUE(save_fits)) {
  dir.create(fit_dir, recursive = TRUE, showWarnings = FALSE)
}

save_fit_object <- function(fit, rep_id, engine) {
  if (!isTRUE(save_fits)) {
    return(invisible(NULL))
  }
  file <- file.path(fit_dir, sprintf("fit_rep%03d_%s.rds", rep_id, engine))
  temporary <- paste0(file, ".tmp-", Sys.getpid())
  saveRDS(fit, file = temporary, compress = "gzip")
  if (file.exists(file) && !file.remove(file)) stop("Could not replace fit checkpoint: ", file, call. = FALSE)
  if (!file.rename(temporary, file)) stop("Could not atomically install fit checkpoint: ", file, call. = FALSE)
  invisible(file)
}

write_checkpoint_csv <- function(x, path) {
  temporary <- paste0(path, ".tmp-", Sys.getpid())
  write.csv(x, temporary, row.names = FALSE)
  if (file.exists(path) && !file.remove(path)) stop("Could not replace checkpoint: ", path, call. = FALSE)
  if (!file.rename(temporary, path)) stop("Could not atomically install checkpoint: ", path, call. = FALSE)
  invisible(path)
}

smooth_mu <- function(x) 0.55 * sin(2 * pi * x)
smooth_sigma <- function(x) (0.70 * sigma_signal_multiplier) * (2 * (x - 0.5)^2 - 1 / 6)
smooth_theta <- function(x) 0.22 * sin(pi * x) - 0.10

nbi_seed_registry <- function(rep_id, engine = NULL) {
  engine_offset <- if (is.null(engine)) 0L else match(engine, engines, nomatch = length(engines) + 1L)
  origin <- base_seed + 1000L * as.integer(rep_id)
  list(
    training_covariate_seed = origin + 11L,
    training_response_seed = origin + 12L,
    test_response_seed = origin + 13L,
    diagnostic_seed = origin + 100L + engine_offset,
    predictive_seed = origin + 200L + engine_offset
  )
}

with_preserved_seed <- function(seed, code) {
  had_seed <- exists(".Random.seed", envir = .GlobalEnv, inherits = FALSE)
  old_seed <- if (had_seed) get(".Random.seed", envir = .GlobalEnv, inherits = FALSE) else NULL
  on.exit({
    if (had_seed) assign(".Random.seed", old_seed, envir = .GlobalEnv)
    else if (exists(".Random.seed", envir = .GlobalEnv, inherits = FALSE)) rm(".Random.seed", envir = .GlobalEnv)
  }, add = TRUE)
  set.seed(as.integer(seed))
  force(code)
}

make_covariates <- function(seed) {
  with_preserved_seed(seed, {
    subject_x1 <- stats::rnorm(n_subject)
    subject_x2 <- stats::rbinom(n_subject, 1, 0.5)
    subject_s1 <- stats::runif(n_subject)
    data.frame(
      x1 = rep(subject_x1, each = length(times)),
      x2 = rep(subject_x2, each = length(times)),
      s1 = rep(subject_s1, each = length(times)),
      time_scaled = rep(sim_rescale01(times), times = n_subject),
      stringsAsFactors = FALSE
    )
  })
}

make_dat <- function(response_seed, covariates) {
  simulate_longitudinal_dataset(
    n = n_subject,
    times = times,
    margin_dist = gamlss.dist::NBI(),
    copula_dist = nbi_copula_code,
    covariates = covariates[, c("x1", "x2", "s1", "time_scaled"), drop = FALSE],
    margin_params = list(
      mu = function(d) exp(1.10 + 0.18 * d$x1 - 0.12 * d$x2 + 0.10 * d$time_scaled + smooth_mu(d$s1)),
      sigma = function(d) exp(-0.45 + 0.10 * d$x1 + 0.08 * d$time_scaled + smooth_sigma(d$s1))
    ),
    copula_params = list(
      theta = function(e) {
        time_left_scaled <- sim_rescale01(e$time_left)
        exp(theta_intercept + theta_time_coef * time_left_scaled + smooth_theta(e$s1))
      }
    ),
    seed = response_seed
  )
}

curve_from_gamlss <- function(fit, dat, grid, parameter) {
  ref <- data.frame(
    x1 = mean(dat$x1),
    x2 = 0,
    time_scaled = mean(dat$time_scaled),
    s1 = grid
  )
  ref0 <- ref
  ref0$s1 <- mean(dat$s1)
  eta <- as.numeric(predict(fit, what = parameter, type = "link", newdata = ref))
  eta0 <- as.numeric(predict(fit, what = parameter, type = "link", newdata = ref0))
  fitted <- eta - eta0
  truth_fun <- switch(parameter, mu = smooth_mu, sigma = smooth_sigma)
  truth <- truth_fun(grid) - truth_fun(mean(dat$s1))
  fitted <- fitted - mean(fitted)
  truth <- truth - mean(truth)
  data.frame(parameter = parameter, x = grid, fitted = fitted, truth = truth, error = fitted - truth)
}

curve_from_ours <- function(fit, dat, parameter) {
  s_name <- names(fit$par_s[[parameter]])[1]
  B <- fit$model_matrix$s[[parameter]][[s_name]]
  beta <- fit$par_s[[parameter]][[s_name]]
  x <- attr(B, "smooth_x")
  fitted <- as.numeric(B %*% beta)
  truth_fun <- switch(parameter, mu = smooth_mu, sigma = smooth_sigma, theta = smooth_theta)
  truth <- truth_fun(x)
  fitted <- fitted - mean(fitted, na.rm = TRUE)
  truth <- truth - mean(truth, na.rm = TRUE)
  out <- data.frame(parameter = parameter, x = x, fitted = fitted, truth = truth, error = fitted - truth)
  stats::aggregate(cbind(fitted, truth, error) ~ parameter + x, out, mean)
}

fit_gamlss <- function(dat) {
  gamlss_data <- dat[, c("response", "x1", "x2", "time_scaled", "s1")]
  assign(".nbi_gamlss_vcov_data", gamlss_data, envir = .GlobalEnv)
  gamlss(
    response ~ x1 + x2 + time_scaled + pb(s1, df = 6),
    sigma.formula = ~ x1 + time_scaled + pb(s1, df = 6),
    family = NBI(),
    data = .nbi_gamlss_vcov_data,
    method = RS(),
    control = gamlss.control(n.cyc = 80, trace = FALSE),
    i.control = glim.control(cyc = 80, bf.cyc = 80, bf.tol = 1e-4)
  )
}

gamlss_prediction_data <- function(dat) {
  dat[, c("x1", "x2", "time_scaled", "s1"), drop = FALSE]
}

fit_ours <- function(dat, engine) {
  is_joint <- identical(engine, "ours_rs_joint")
  gamlss.longitudinal::gamlss_longitudinal(
    dataset = dat,
    margin_dist = gamlss.dist::NBI(),
    copula_dist = nbi_copula_code,
    time_var = "time",
    subject_var = "subject",
    mu.formula = "response ~ x1 + x2 + time_scaled + s(s1, bs = 'ps', k = 8)",
    sigma.formula = "~ x1 + time_scaled + s(s1, bs = 'ps', k = 8)",
    theta.formula = "~ time_scaled + s(s1, bs = 'ps', k = 8)",
    zeta.formula = "~ 1",
    method = "RS",
    include_dlcopdpar = is_joint,
    warm_start_joint = is_joint,
    warm_start_joint_iter = warm_start_joint_iter,
    max_outer_iter = max_outer_iter,
    max_inner_iter = max_inner_iter,
    max_elapsed_sec = max_elapsed_sec,
    outer_stop_crit = 0.001,
    inner_stop_crit = 0.001,
    start_step_size = start_step_size,
    step_adjustment = step_adjustment,
    compute_vcov = FALSE,
    vcov_method = vcov_method_longitudinal,
    lambda_start = lambda_start,
    rs_update_lambda = rs_update_lambda,
    verbose = 0,
    discrete_score_method = "analytical"
  )
}

rows <- list()
curves <- list()
fixed_rows <- list()
joint_rows <- list()
predictive_rows <- list()
grid <- seq(0.01, 0.99, length.out = 100)

log_path <- file.path(out_dir, "nbi_sigma_compare_logs.csv")
curve_path <- file.path(out_dir, "nbi_sigma_compare_curves.csv")
fixed_path <- file.path(out_dir, "fixed_effects_by_rep.csv")
joint_path <- file.path(out_dir, "joint_distribution_metrics_by_rep.csv")
predictive_path <- file.path(out_dir, "predictive_scores_by_rep.csv")
attempt_checkpoint_dir <- file.path(out_dir, "attempt_checkpoints")
dir.create(attempt_checkpoint_dir, recursive = TRUE, showWarnings = FALSE)

attempt_checkpoint_path <- function(rep_id, engine, retry_index) {
  file.path(attempt_checkpoint_dir, sprintf("rep%03d_%s_try%02d.rds", rep_id, engine, retry_index))
}

attempt_checkpoint_issues <- function(x) {
  issues <- character()
  if (!is.list(x) || !identical(x$runner_contract_version, runner_contract_version)) return("runner_contract")
  if (!identical(x$phase1_contract_version, phase1_contract_version)) return("phase1_contract")
  if (!identical(x$runner_settings_signature, runner_settings_signature)) return("runner_settings")
  if (!identical(x$runner_settings_sha256, runner_settings_sha256) || !identical(x$runner_sha256, runner_sha256) ||
      !identical(x$package_source_path, package_source_identity) || !identical(x$package_version, runner_package_version) ||
      !identical(x$package_source_sha256, package_source_sha256)) return("source_identity")
  if (!is.data.frame(x$log) || nrow(x$log) != 1L) return("log_row")
  log <- x$log
  if (!identical(as.character(log$margin_family), "NBI") || !identical(as.character(log$copula_family), nbi_copula_label) ||
      !identical(as.character(log$copula_code), nbi_copula_code)) issues <- c(issues, "family_metadata")
  seed_columns <- c("training_covariate_seed", "training_response_seed", "test_response_seed", "diagnostic_seed", "predictive_seed")
  if (!all(seed_columns %in% names(log)) || any(!is.finite(as.numeric(log[1, intersect(seed_columns, names(log)), drop = TRUE])))) issues <- c(issues, "seed_registry")
  if (!isTRUE(log$success)) {
    empty_outputs <- all(vapply(c("curves", "fixed", "joint", "predictive"), function(name) is.data.frame(x[[name]]) && !nrow(x[[name]]), logical(1)))
    if (!empty_outputs) issues <- c(issues, "failed_attempt_has_outputs")
    return(unique(issues))
  }
  expected_parameters <- if (identical(as.character(log$engine), "gamlss")) c("mu", "sigma") else c("mu", "sigma", "theta")
  expected_fixed <- if (identical(as.character(log$engine), "gamlss")) 7L else 9L
  expected_curve_points <- if (identical(as.character(log$engine), "gamlss")) length(grid) else min(n_subject, length(grid))
  curve_ok <- is.data.frame(x$curves) && "parameter" %in% names(x$curves) &&
    setequal(unique(x$curves$parameter), expected_parameters) && all(table(x$curves$parameter) >= expected_curve_points)
  completeness_issues <- character()
  if (!curve_ok) completeness_issues <- c(completeness_issues, "curves")
  if (!is.data.frame(x$fixed) || nrow(x$fixed) != expected_fixed) completeness_issues <- c(completeness_issues, "fixed")
  if (!is.data.frame(x$joint) || nrow(x$joint) != 1L) completeness_issues <- c(completeness_issues, "diagnostic")
  predictive_ok <- is.data.frame(x$predictive) &&
    ((!isTRUE(compute_predictive_scores) && nrow(x$predictive) == 0L) ||
      (isTRUE(compute_predictive_scores) && nrow(x$predictive) == length(variogram_p_values) && "variogram_p" %in% names(x$predictive) &&
        setequal(as.numeric(x$predictive$variogram_p), as.numeric(variogram_p_values))))
  if (!predictive_ok) completeness_issues <- c(completeness_issues, "predictive")
  if (is.data.frame(x$fixed) && nrow(x$fixed) && (!all(is.finite(x$fixed$estimate)) || !all(is.finite(x$fixed$true_value)) || any(is.finite(x$fixed$std_error) & x$fixed$std_error < 0))) completeness_issues <- c(completeness_issues, "fixed_substantive_values")
  if (is.data.frame(x$curves) && nrow(x$curves) && (!all(is.finite(x$curves$fitted)) || !all(is.finite(x$curves$truth)))) completeness_issues <- c(completeness_issues, "smooth_substantive_values")
  if (is.data.frame(x$joint) && nrow(x$joint) && !all(vapply(c("logLik", "rosenblatt_ks", "rosenblatt_cvm", "abs_rosenblatt_lag1_cor", "abs_rosenblatt_normal_lag1_cor", "rosenblatt_mean_abs_time_mean", "rosenblatt_normal_mean_abs_time_mean"), function(name) name %in% names(x$joint) && all(is.finite(x$joint[[name]])), logical(1)))) completeness_issues <- c(completeness_issues, "diagnostic_substantive_values")
  if (is.data.frame(x$predictive) && nrow(x$predictive) && !all(vapply(c("test_log_score_joint", "test_log_score_marginal", "test_log_score_copula", "test_log_score_per_obs", "variogram_score", "predictive_nsim"), function(name) name %in% names(x$predictive) && all(is.finite(x$predictive[[name]])), logical(1)))) completeness_issues <- c(completeness_issues, "predictive_substantive_values")
  if (isTRUE(log$descriptive_outputs_complete) && length(completeness_issues)) issues <- c(issues, "declared_complete_but_incomplete", completeness_issues)
  if (!isTRUE(log$descriptive_outputs_complete) && !length(completeness_issues)) issues <- c(issues, "declared_incomplete_but_complete")
  unique(issues)
}

valid_attempt_checkpoint <- function(x) {
  length(attempt_checkpoint_issues(x)) == 0L
}

read_attempt_checkpoint <- function(path) {
  tryCatch({
    x <- readRDS(path)
    if (valid_attempt_checkpoint(x)) x else NULL
  }, error = function(e) NULL)
}

write_attempt_checkpoint <- function(x, path) {
  if (!valid_attempt_checkpoint(x)) stop("Refusing to write an incomplete NBI/Clayton attempt checkpoint.", call. = FALSE)
  temporary <- paste0(path, ".tmp-", Sys.getpid())
  saveRDS(x, temporary, compress = "gzip")
  if (file.exists(path) && !file.remove(path)) stop("Could not replace attempt checkpoint: ", path, call. = FALSE)
  if (!file.rename(temporary, path)) stop("Could not atomically install attempt checkpoint: ", path, call. = FALSE)
  invisible(path)
}

checkpoints <- if (isTRUE(resume)) lapply(list.files(attempt_checkpoint_dir, pattern = "[.]rds$", full.names = TRUE), read_attempt_checkpoint) else list()
checkpoints <- checkpoints[!vapply(checkpoints, is.null, logical(1))]
if (length(checkpoints)) {
  rows <- lapply(checkpoints, `[[`, "log")
  curves <- lapply(checkpoints, `[[`, "curves"); curves <- curves[vapply(curves, nrow, integer(1)) > 0L]
  fixed_rows <- lapply(checkpoints, `[[`, "fixed"); fixed_rows <- fixed_rows[vapply(fixed_rows, nrow, integer(1)) > 0L]
  joint_rows <- lapply(checkpoints, `[[`, "joint"); joint_rows <- joint_rows[vapply(joint_rows, nrow, integer(1)) > 0L]
  predictive_rows <- lapply(checkpoints, `[[`, "predictive"); predictive_rows <- predictive_rows[vapply(predictive_rows, nrow, integer(1)) > 0L]
}

prior_attempt_count <- function(rep_id, engine) {
  if (!length(rows)) return(0L)
  log_df <- do.call(rbind, rows)
  sum(log_df$rep == rep_id & log_df$engine == engine)
}

already_done <- function(rep_id, engine) {
  if (!isTRUE(resume) || !length(rows)) return(FALSE)
  log_df <- do.call(rbind, rows)
  prior <- log_df$rep == rep_id & log_df$engine == engine
  any(prior & log_df$publication_candidate %in% TRUE) || sum(prior) >= max_attempts_per_fit
}

fixed_from_gamlss <- function(fit, rep_id) {
  out <- list()
  co_mu <- stats::coef(fit, what = "mu")
  co_sigma <- stats::coef(fit, what = "sigma")
  se_all <- if (isTRUE(compute_se)) {
    tryCatch(stats::vcov(fit, type = "se", full = TRUE), error = function(e) NULL)
  } else {
    NULL
  }
  se_by_parameter <- list(mu = NULL, sigma = NULL)
  if (!is.null(se_all)) {
    se_all <- as.numeric(se_all)
    if (length(se_all) >= length(co_mu) + length(co_sigma)) {
      se_by_parameter$mu <- se_all[seq_along(co_mu)]
      se_by_parameter$sigma <- se_all[length(co_mu) + seq_along(co_sigma)]
      names(se_by_parameter$mu) <- sub("^\\(Intercept\\)$", "intercept", names(co_mu))
      names(se_by_parameter$sigma) <- sub("^\\(Intercept\\)$", "intercept", names(co_sigma))
    }
  }
  get_se <- function(parameter, term) {
    se_vec <- se_by_parameter[[parameter]]
    if (is.null(se_vec) || !term %in% names(se_vec)) {
      return(NA_real_)
    }
    as.numeric(se_vec[[term]])
  }
  map <- list(
    mu = list(coef = co_mu, truth = c(intercept = 1.10, x1 = 0.18, x2 = -0.12, time_scaled = 0.10)),
    sigma = list(coef = co_sigma, truth = c(intercept = -0.45, x1 = 0.10, time_scaled = 0.08))
  )
  for (parameter in names(map)) {
    coef_vec <- map[[parameter]]$coef
    truth <- map[[parameter]]$truth
    names(coef_vec) <- sub("^\\(Intercept\\)$", "intercept", names(coef_vec))
    for (term in intersect(names(truth), names(coef_vec))) {
      out[[length(out) + 1L]] <- data.frame(
        model = "gamlss",
        parameter = parameter,
        term = term,
        rep = rep_id,
        estimate = unname(coef_vec[[term]]),
        std_error = get_se(parameter, term),
        true_value = unname(truth[[term]]),
        stringsAsFactors = FALSE
      )
    }
  }
  do.call(rbind, out)
}

fixed_from_ours <- function(fit, rep_id, engine) {
  truth <- c(
    mu.intercept = 1.10,
    mu.x1 = 0.18,
    mu.x2 = -0.12,
    mu.time_scaled = 0.10,
    sigma.intercept = -0.45,
    sigma.x1 = 0.10,
    sigma.time_scaled = 0.08,
    theta.intercept = theta_intercept,
    theta.time_scaled = theta_time_coef
  )
  common <- intersect(names(truth), names(fit$par))
  se_vec <- rep(NA_real_, length(fit$par))
  names(se_vec) <- names(fit$par)
  if (isTRUE(compute_se)) {
    sum_obj <- tryCatch(
      summary(fit, include_vcov = TRUE, vcov_method = vcov_method_longitudinal),
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
  pieces <- strsplit(common, "\\.")
  data.frame(
    model = engine,
    parameter = vapply(pieces, `[`, character(1), 1),
    term = vapply(pieces, function(x) paste(x[-1], collapse = "."), character(1)),
    rep = rep_id,
    estimate = as.numeric(fit$par[common]),
    std_error = as.numeric(se_vec[common]),
    true_value = as.numeric(truth[common]),
    stringsAsFactors = FALSE
  )
}

scalar_numeric <- function(x, preferred_name = NULL) {
  if (is.null(x) || length(x) == 0L) {
    return(NA_real_)
  }
  if (!is.null(preferred_name) && preferred_name %in% names(x)) {
    return(as.numeric(x[[preferred_name]]))
  }
  as.numeric(x[[1L]])
}

scalar_integer <- function(x) {
  if (is.null(x) || length(x) == 0L) {
    return(NA_integer_)
  }
  as.integer(x[[1L]])
}

scalar_logical <- function(x) {
  if (is.null(x) || length(x) == 0L) {
    return(NA)
  }
  isTRUE(x[[1L]])
}

raw_convergence_evidence <- function(fit, engine, ok) {
  api <- if (engine == "gamlss") "gamlss::gamlss_result" else "gamlss.longitudinal::fit_convergence"
  basis <- if (engine == "gamlss") "gamlss_explicit_converged_v1" else "gamlss.longitudinal_explicit_convergence_v1"
  indicator_name <- if (engine == "gamlss") "fit$converged" else "fit$convergence$converged"
  if (!ok) return(list(api = api, basis = basis, indicator_name = indicator_name, status = "fit_execution_failed",
    indicator = NA, loglik = NA_real_, deviance = NA_real_, coefficient_count = 0L, coefficient_bad = 0L,
    fitted_count = 0L, fitted_bad = 0L, iterations = NA_integer_, cap = if (engine == "gamlss") 20L else max_outer_iter,
    hit_outer = FALSE, hit_stall = FALSE, hit_drop = FALSE, converged = FALSE))
  flatten <- function(x) suppressWarnings(as.numeric(unlist(x, recursive = TRUE, use.names = FALSE)))
  coefficients <- tryCatch(flatten(stats::coef(fit)), error = function(e) numeric())
  if (!length(coefficients)) coefficients <- flatten(fit$par)
  fitted <- tryCatch(flatten(stats::fitted(fit)), error = function(e) numeric())
  if (!length(fitted)) fitted <- flatten(fit$fitted.values)
  loglik <- extract_fit_metric(fit, "logLik")
  deviance <- suppressWarnings(as.numeric(fit$deviance)[1L]); if (!is.finite(deviance) && is.finite(loglik)) deviance <- -2 * loglik
  indicator <- if (engine == "gamlss") isTRUE(fit$converged) else scalar_logical(fit$convergence$converged)
  iterations <- if (engine == "gamlss") scalar_integer(fit$iter) else scalar_integer(fit$convergence$outer_iterations)
  cap <- if (engine == "gamlss") suppressWarnings(as.integer(fit$control$n.cyc)[1L]) else max_outer_iter
  if (!is.finite(cap) || cap < 1L) cap <- if (engine == "gamlss") 20L else max_outer_iter
  hit_outer <- if (engine == "gamlss") iterations >= cap else isTRUE(fit$convergence$hit_outer_limit)
  hit_stall <- engine != "gamlss" && isTRUE(fit$convergence$hit_max_stall)
  hit_drop <- engine != "gamlss" && isTRUE(fit$convergence$hit_raw_loglik_deterioration)
  status <- if (!is.finite(loglik) || !is.finite(deviance)) "nonfinite_objective" else if (!length(coefficients) || any(!is.finite(coefficients))) "nonfinite_or_missing_coefficients" else if (!length(fitted) || any(!is.finite(fitted))) "nonfinite_or_missing_fitted_values" else if (!is.finite(iterations) || iterations < 1L || !is.finite(cap) || cap < 1L) "missing_or_invalid_iterations" else if (hit_outer || iterations >= cap) "outer_iteration_cap_reached" else if (hit_stall) "maximum_stall_reached" else if (hit_drop) "raw_loglik_deterioration" else if (is.na(indicator)) "missing_explicit_convergence_indicator" else if (!indicator) "explicit_optimizer_nonconvergence" else "explicit_optimizer_convergence"
  list(api = api, basis = basis, indicator_name = indicator_name, status = status, indicator = indicator,
    loglik = loglik, deviance = deviance, coefficient_count = length(coefficients), coefficient_bad = sum(!is.finite(coefficients)),
    fitted_count = length(fitted), fitted_bad = sum(!is.finite(fitted)), iterations = iterations, cap = cap,
    hit_outer = hit_outer, hit_stall = hit_stall, hit_drop = hit_drop,
    converged = status == "explicit_optimizer_convergence")
}

extract_fit_metric <- function(fit, metric) {
  log_lik <- tryCatch(stats::logLik(fit), error = function(e) NA_real_)
  out <- scalar_numeric(log_lik, "joint")
  if (identical(metric, "logLik")) return(out)
  if (identical(metric, "df")) {
    df <- tryCatch(attr(stats::logLik(fit), "df"), error = function(e) NA_real_)
    return(scalar_numeric(df, "joint"))
  }
  NA_real_
}

clip_u <- function(x) pmax(pmin(as.numeric(x), 1 - 1e-8), 1e-8)

adjacent_pair_rows <- function(dat) {
  dat <- dat[order(dat$subject, dat$time), ]
  split_dat <- split(seq_len(nrow(dat)), dat$subject)
  out <- lapply(split_dat, function(idx) {
    if (length(idx) < 2L) return(NULL)
    data.frame(row1 = idx[-length(idx)], row2 = idx[-1])
  })
  out <- out[!vapply(out, is.null, logical(1))]
  if (length(out) == 0L) return(data.frame(row1 = integer(0), row2 = integer(0)))
  do.call(rbind, out)
}

uniform_ks_stat <- function(x) {
  x <- x[is.finite(x)]
  if (length(x) < 2L) return(NA_real_)
  x <- sort(pmax(pmin(x, 1), 0))
  n <- length(x)
  max(max(abs(x - seq_len(n) / n)), max(abs(x - (seq_len(n) - 1) / n)))
}

uniform_cvm_stat <- function(x) {
  x <- x[is.finite(x)]
  if (length(x) < 2L) return(NA_real_)
  x <- sort(pmax(pmin(x, 1), 0))
  n <- length(x)
  1 / (12 * n) + sum((x - (2 * seq_len(n) - 1) / (2 * n))^2)
}

rosenblatt_lag1_cor <- function(r, dat, normal_scores = FALSE) {
  dat_tmp <- dat[order(dat$subject, dat$time), , drop = FALSE]
  r_tmp <- r[order(dat$subject, dat$time)]
  if (normal_scores) {
    r_tmp <- stats::qnorm(clip_u(r_tmp))
  }
  pairs <- adjacent_pair_rows(dat_tmp)
  suppressWarnings(stats::cor(r_tmp[pairs$row1], r_tmp[pairs$row2], use = "complete.obs"))
}

rosenblatt_subject_mean <- function(r, dat, normal_scores = FALSE) {
  dat_tmp <- dat[order(dat$subject, dat$time), , drop = FALSE]
  r_tmp <- r[order(dat$subject, dat$time)]
  if (normal_scores) {
    r_tmp <- stats::qnorm(clip_u(r_tmp))
    center <- 0
  } else {
    center <- 0.5
  }
  mean_by_time <- stats::aggregate(r_tmp, by = list(time = dat_tmp$time), FUN = mean, na.rm = TRUE)
  names(mean_by_time)[2] <- "mean_residual"
  mean(abs(mean_by_time$mean_residual - center), na.rm = TRUE)
}

joint_metrics_from_rosenblatt <- function(model_name, rep_id, r, dat, logLik, df) {
  lag1_cor <- rosenblatt_lag1_cor(r, dat, normal_scores = FALSE)
  lag1_z_cor <- rosenblatt_lag1_cor(r, dat, normal_scores = TRUE)
  data.frame(
    scenario = sprintf("n%d_d%d_nbi_signal%s", n_subject, length(times), sigma_signal_multiplier),
    n = n_subject,
    d = length(times),
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

randomized_nbi_pit <- function(response, mu, sigma, seed) {
  upper <- gamlss.dist::pNBI(response, mu = mu, sigma = sigma)
  lower <- gamlss.dist::pNBI(response - 1, mu = mu, sigma = sigma)
  lower <- pmin(pmax(as.numeric(lower), 0), 1)
  upper <- pmin(pmax(as.numeric(upper), 0), 1)
  with_preserved_seed(seed, lower + stats::runif(length(upper)) * pmax(upper - lower, 0))
}

joint_metrics_gamlss <- function(fit, dat, rep_id, diagnostic_seed) {
  pred_data <- gamlss_prediction_data(dat)
  pred_mu <- as.numeric(predict(fit, what = "mu", type = "link", newdata = pred_data))
  pred_sigma <- as.numeric(predict(fit, what = "sigma", type = "link", newdata = pred_data))
  u <- randomized_nbi_pit(
    dat$response, gamlss.dist::NBI()$mu.linkinv(pred_mu),
    gamlss.dist::NBI()$sigma.linkinv(pred_sigma), diagnostic_seed
  )
  joint_metrics_from_rosenblatt("gamlss", rep_id, clip_u(u), dat, extract_fit_metric(fit, "logLik"), extract_fit_metric(fit, "df"))
}

joint_metrics_ours <- function(fit, rep_id, engine, diagnostic_seed) {
  copula_link_fit <- gamlss.longitudinal:::get_copula_dist(fit$copula_dist)$copula_link
  eta_inv <- gamlss.longitudinal:::calc_eta(fit$par, fit$model_matrix, fit$margin_dist, copula_link_fit, fit$par_s)$eta_inv
  u <- randomized_nbi_pit(fit$response, eta_inv$mu, eta_inv$sigma, diagnostic_seed)
  u <- clip_u(u)
  dat <- data.frame(subject = fit$response_subject, time = fit$response_margin)
  pairs <- adjacent_pair_rows(dat)
  theta_rows <- which(fit$response_margin < max(fit$response_margin, na.rm = TRUE))
  theta_lookup <- rep(NA_integer_, length(fit$response))
  theta_lookup[theta_rows] <- seq_along(theta_rows)
  theta_idx <- theta_lookup[pairs$row1]
  theta <- as.numeric(eta_inv$theta[theta_idx])
  if (!identical(fit$copula_dist, nbi_copula_code)) stop("NBI recovery fit used the wrong copula family.", call. = FALSE)
  fam_num <- gamlss.longitudinal:::.copula_family_number(fit$copula_dist)
  conditional_rosenblatt <- tryCatch(
    VineCopula::BiCopHfunc1(u[pairs$row1], u[pairs$row2], family = fam_num, par = theta, par2 = 0),
    error = function(e) rep(NA_real_, nrow(pairs))
  )
  r <- u
  r[pairs$row2] <- clip_u(conditional_rosenblatt)
  joint_metrics_from_rosenblatt(engine, rep_id, r, dat, extract_fit_metric(fit, "logLik"), extract_fit_metric(fit, "df"))
}

response_matrix <- function(dat, value_col = "response") {
  dat <- dat[order(dat$subject, dat$time), , drop = FALSE]
  ids <- sort(unique(dat$subject))
  time_values <- sort(unique(dat$time))
  out <- matrix(NA_real_, nrow = length(ids), ncol = length(time_values))
  row_idx <- match(dat$subject, ids)
  col_idx <- match(dat$time, time_values)
  out[cbind(row_idx, col_idx)] <- dat[[value_col]]
  out
}

simulate_array_from_columns <- function(sim_df, dat, subject_col = "subject") {
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
  if (length(dim(simulated_array)) != 3L) return(NA_real_)
  n <- nrow(observed_mat)
  d <- ncol(observed_mat)
  if (d < 2L || dim(simulated_array)[2] != n || dim(simulated_array)[3] != d) return(NA_real_)
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

simulate_predictive_gamlss <- function(fit, dat, nsim, seed) {
  dat <- dat[order(dat$subject, dat$time), , drop = FALSE]
  pred_data <- gamlss_prediction_data(dat)
  pred_mu <- as.numeric(predict(fit, what = "mu", type = "link", newdata = pred_data))
  pred_sigma <- as.numeric(predict(fit, what = "sigma", type = "link", newdata = pred_data))
  n <- length(unique(dat$subject))
  d <- length(unique(dat$time))
  out <- array(NA_real_, dim = c(nsim, n, d))
  draws <- with_preserved_seed(seed, lapply(seq_len(nsim), function(s) {
    y <- gamlss.dist::qNBI(
      stats::runif(nrow(dat)),
      mu = gamlss.dist::NBI()$mu.linkinv(pred_mu),
      sigma = gamlss.dist::NBI()$sigma.linkinv(pred_sigma)
    )
    matrix(y, nrow = n, ncol = d, byrow = TRUE)
  }))
  for (s in seq_len(nsim)) out[s, , ] <- draws[[s]]
  out
}

simulate_predictive_ours <- function(fit, dat, nsim, seed) {
  dat <- dat[order(dat$subject, dat$time), , drop = FALSE]
  sim_df <- with_preserved_seed(seed, stats::simulate(fit, nsim = nsim, newdata = dat))
  simulate_array_from_columns(sim_df, dat, subject_col = "subject")
}

predictive_scores_gamlss <- function(fit, test_dat, rep_id, predictive_seed, nsim = predictive_nsim) {
  test_dat <- test_dat[order(test_dat$subject, test_dat$time), , drop = FALSE]
  pred_data <- gamlss_prediction_data(test_dat)
  pred_mu <- as.numeric(predict(fit, what = "mu", type = "link", newdata = pred_data))
  pred_sigma <- as.numeric(predict(fit, what = "sigma", type = "link", newdata = pred_data))
  log_d <- gamlss.dist::dNBI(
    test_dat$response,
    mu = gamlss.dist::NBI()$mu.linkinv(pred_mu),
    sigma = gamlss.dist::NBI()$sigma.linkinv(pred_sigma),
    log = TRUE
  )
  joint_log_score <- sum(log_d, na.rm = TRUE)
  sim_array <- simulate_predictive_gamlss(fit, test_dat, nsim, predictive_seed)
  y_obs <- response_matrix(test_dat)
  base <- data.frame(
    scenario = sprintf("n%d_d%d_nbi_signal%s", n_subject, length(times), sigma_signal_multiplier),
    n = n_subject,
    d = length(times),
    rep = rep_id,
    model = "gamlss",
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

predictive_scores_ours <- function(fit, test_dat, rep_id, engine, predictive_seed, nsim = predictive_nsim) {
  copula_link_fit <- gamlss.longitudinal:::get_copula_dist(fit$copula_dist)$copula_link
  eta_inv <- gamlss.longitudinal:::calc_eta(fit$par, fit$model_matrix, fit$margin_dist, copula_link_fit, fit$par_s)$eta_inv
  lik <- gamlss.longitudinal:::calc_likelihood_minimal(
    eta_inv,
    mm = fit$model_matrix$x,
    margin_dist = fit$margin_dist,
    copula_dist = fit$copula_dist,
    calc_d2 = FALSE,
    response = test_dat$response,
    response_margin = test_dat$time,
    response_subject = test_dat$subject
  )
  sim_array <- simulate_predictive_ours(fit, test_dat, nsim, predictive_seed)
  y_obs <- response_matrix(test_dat)
  base <- data.frame(
    scenario = sprintf("n%d_d%d_nbi_signal%s", n_subject, length(times), sigma_signal_multiplier),
    n = n_subject,
    d = length(times),
    rep = rep_id,
    model = engine,
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

for (rep_id in seq_len(reps)) {
  data_seeds <- nbi_seed_registry(rep_id)
  fixed_covariates <- make_covariates(data_seeds$training_covariate_seed)
  dat <- make_dat(data_seeds$training_response_seed, fixed_covariates)
  test_dat <- make_dat(data_seeds$test_response_seed, fixed_covariates)

  for (engine in engines) {
    if (already_done(rep_id, engine)) {
      cat(sprintf("rep %d/%d | %s already complete, skipping\n", rep_id, reps, engine))
      next
    }
    cat(sprintf("rep %d/%d | %s\n", rep_id, reps, engine))
    attempt_seeds <- nbi_seed_registry(rep_id, engine)
    start <- Sys.time()
    fit <- tryCatch(
      if (engine == "gamlss") fit_gamlss(dat) else fit_ours(dat, engine),
      error = function(e) e
    )
    elapsed <- as.numeric(difftime(Sys.time(), start, units = "secs"))
    ok <- !inherits(fit, "error")
    if (ok) {
      save_fit_object(fit, rep_id, engine)
    }
    curve <- NULL
    fixed_piece <- data.frame()
    joint_piece <- data.frame()
    predictive_piece <- data.frame()
    if (ok) {
      curve <- if (engine == "gamlss") {
        rbind(
          curve_from_gamlss(fit, dat, grid, "mu"),
          curve_from_gamlss(fit, dat, grid, "sigma")
        )
      } else {
        rbind(
          curve_from_ours(fit, dat, "mu"),
          curve_from_ours(fit, dat, "sigma"),
          curve_from_ours(fit, dat, "theta")
        )
      }
      curve$rep <- rep_id
      curve$engine <- engine
      fixed_piece <- if (engine == "gamlss") {
        fixed_from_gamlss(fit, rep_id)
      } else {
        fixed_from_ours(fit, rep_id, engine)
      }
      joint_piece <- tryCatch(
        if (engine == "gamlss") joint_metrics_gamlss(fit, dat, rep_id, attempt_seeds$diagnostic_seed) else joint_metrics_ours(fit, rep_id, engine, attempt_seeds$diagnostic_seed),
        error = function(e) data.frame()
      )
      if (isTRUE(compute_predictive_scores)) {
        predictive_piece <- tryCatch(
          if (engine == "gamlss") predictive_scores_gamlss(fit, test_dat, rep_id, attempt_seeds$predictive_seed) else predictive_scores_ours(fit, test_dat, rep_id, engine, attempt_seeds$predictive_seed),
          error = function(e) data.frame()
        )
      }
    }
    retry_index <- prior_attempt_count(rep_id, engine) + 1L
    if (is.null(curve)) curve <- data.frame()
    add_retry <- function(x) {
      if (is.data.frame(x) && nrow(x)) x$retry_index <- retry_index
      x
    }
    curve <- add_retry(curve); fixed_piece <- add_retry(fixed_piece)
    joint_piece <- add_retry(joint_piece); predictive_piece <- add_retry(predictive_piece)
    if (nrow(fixed_piece)) {
      fixed_piece$inference_status <- ifelse(is.finite(fixed_piece$std_error), "available",
        if (isTRUE(compute_se)) "unavailable_computation_failed" else "not_requested_registered")
      fixed_piece$inference_denominator <- as.integer(is.finite(fixed_piece$std_error))
    }
    expected_fixed <- if (engine == "gamlss") 7L else 9L
    expected_parameters <- if (engine == "gamlss") c("mu", "sigma") else c("mu", "sigma", "theta")
    expected_curve_points <- if (engine == "gamlss") length(grid) else min(n_subject, length(grid))
    output_issues <- c(
      if (!nrow(curve) || !"parameter" %in% names(curve) || !setequal(unique(curve$parameter), expected_parameters) || any(table(curve$parameter) < expected_curve_points)) "curves" else NULL,
      if (nrow(fixed_piece) != expected_fixed) "fixed" else NULL,
      if (nrow(joint_piece) != 1L) "diagnostic" else NULL,
      if (isTRUE(compute_predictive_scores) && nrow(predictive_piece) != length(variogram_p_values)) "predictive" else NULL,
      if (nrow(fixed_piece) && (!all(is.finite(fixed_piece$estimate)) || !all(is.finite(fixed_piece$true_value)) || any(is.finite(fixed_piece$std_error) & fixed_piece$std_error < 0))) "fixed_substantive_values" else NULL,
      if (nrow(curve) && (!all(is.finite(curve$fitted)) || !all(is.finite(curve$truth)))) "smooth_substantive_values" else NULL,
      if (nrow(joint_piece) && !all(vapply(c("logLik", "rosenblatt_ks", "rosenblatt_cvm", "abs_rosenblatt_lag1_cor", "abs_rosenblatt_normal_lag1_cor", "rosenblatt_mean_abs_time_mean", "rosenblatt_normal_mean_abs_time_mean"), function(name) name %in% names(joint_piece) && all(is.finite(joint_piece[[name]])), logical(1)))) "diagnostic_substantive_values" else NULL,
      if (nrow(predictive_piece) && !all(vapply(c("test_log_score_joint", "test_log_score_marginal", "test_log_score_copula", "test_log_score_per_obs", "variogram_score", "predictive_nsim"), function(name) name %in% names(predictive_piece) && all(is.finite(predictive_piece[[name]])), logical(1)))) "predictive_substantive_values" else NULL
    )
    descriptive_outputs_complete <- ok && !length(output_issues)
    raw_convergence <- raw_convergence_evidence(fit, engine, ok)
    converged_value <- raw_convergence$converged
    log_piece <- data.frame(
      evidence_status = "post_phase1_production",
      study_id = "nbi_clayton_main_recovery",
      margin_family = "NBI",
      copula_family = nbi_copula_label,
      copula_code = nbi_copula_code,
      runner_contract_version = runner_contract_version,
      phase1_contract_version = phase1_contract_version,
      runner_settings_signature = runner_settings_signature,
      runner_settings_sha256 = runner_settings_sha256,
      runner_sha256 = runner_sha256,
      package_source_path = package_source_identity,
      package_version = runner_package_version,
      package_source_sha256 = package_source_sha256,
      scenario = sprintf("n%d_d%d_nbi_signal%s", n_subject, length(times), sigma_signal_multiplier),
      n = n_subject,
      d = length(times),
      target_replicates = reps,
      rep = rep_id,
      seed = data_seeds$training_response_seed,
      seed_source = "runner_metadata",
      training_covariate_seed = data_seeds$training_covariate_seed,
      training_response_seed = data_seeds$training_response_seed,
      test_response_seed = data_seeds$test_response_seed,
      diagnostic_seed = attempt_seeds$diagnostic_seed,
      predictive_seed = attempt_seeds$predictive_seed,
      attempted = TRUE,
      retry_index = retry_index,
      engine = engine,
      success = ok,
      execution_success = ok,
      descriptive_outputs_complete = descriptive_outputs_complete,
      output_failure_reason = if (descriptive_outputs_complete) NA_character_ else if (!ok) "fit_execution_failed" else paste(output_issues, collapse = "|"),
      runtime_n_cores = 1L, runtime_backend = "sequential", runtime_rscript_sha256 = nbi_rscript_sha256,
      publication_candidate = isTRUE(converged_value) && descriptive_outputs_complete,
      logLik = if (ok) extract_fit_metric(fit, "logLik") else NA_real_,
      df = if (ok) extract_fit_metric(fit, "df") else NA_real_,
      converged = converged_value,
      raw_convergence_schema = "raw-convergence-2026-09-01.1", raw_convergence_api = raw_convergence$api,
      raw_convergence_basis = raw_convergence$basis, raw_convergence_status = raw_convergence$status,
      raw_convergence_indicator_name = raw_convergence$indicator_name, raw_convergence_indicator_value = raw_convergence$indicator,
      raw_convergence_loglik = raw_convergence$loglik, raw_convergence_deviance = raw_convergence$deviance,
      raw_convergence_coefficient_count = raw_convergence$coefficient_count,
      raw_convergence_coefficient_nonfinite_count = raw_convergence$coefficient_bad,
      raw_convergence_fitted_count = raw_convergence$fitted_count, raw_convergence_fitted_nonfinite_count = raw_convergence$fitted_bad,
      raw_convergence_iteration_count = raw_convergence$iterations, raw_convergence_iteration_cap = raw_convergence$cap,
      raw_convergence_hit_outer_limit = raw_convergence$hit_outer, raw_convergence_hit_max_stall = raw_convergence$hit_stall,
      raw_convergence_hit_raw_loglik_deterioration = raw_convergence$hit_drop,
      stop_reason = if (ok && engine != "gamlss") as.character(fit$convergence$stop_reason %||% NA_character_) else if (ok && !isTRUE(fit$converged)) "optimizer_nonconvergence" else NA_character_,
      iter = if (ok && engine == "gamlss") scalar_integer(fit$iter) else if (ok) scalar_integer(fit$convergence$outer_iterations) else NA_integer_,
      mu_irmse = if (ok && any(curve$parameter == "mu")) sqrt(mean(curve$error[curve$parameter == "mu"]^2, na.rm = TRUE)) else NA_real_,
      sigma_irmse = if (ok && any(curve$parameter == "sigma")) sqrt(mean(curve$error[curve$parameter == "sigma"]^2, na.rm = TRUE)) else NA_real_,
      theta_irmse = if (ok && any(curve$parameter == "theta")) sqrt(mean(curve$error[curve$parameter == "theta"]^2, na.rm = TRUE)) else NA_real_,
      elapsed_sec = elapsed,
      execution_completed_at_utc = format(Sys.time(), "%Y-%m-%dT%H:%M:%SZ", tz = "UTC"),
      error = if (ok) NA_character_ else conditionMessage(fit),
      stringsAsFactors = FALSE
    )
    checkpoint <- list(
      runner_contract_version = runner_contract_version,
      phase1_contract_version = phase1_contract_version,
      runner_settings_signature = runner_settings_signature,
      runner_settings_sha256 = runner_settings_sha256,
      runner_sha256 = runner_sha256,
      package_source_path = package_source_identity,
      package_version = runner_package_version,
      package_source_sha256 = package_source_sha256,
      log = log_piece,
      curves = curve,
      fixed = fixed_piece,
      joint = joint_piece,
      predictive = predictive_piece
    )
    if (ok && !valid_attempt_checkpoint(checkpoint)) {
      incomplete_issues <- attempt_checkpoint_issues(checkpoint)
      stop("Refusing structurally invalid NBI/Clayton checkpoint: ", paste(incomplete_issues, collapse = "|"), call. = FALSE)
    }
    checkpoint_path <- attempt_checkpoint_path(rep_id, engine, retry_index)
    write_attempt_checkpoint(checkpoint, checkpoint_path)
    rows[[length(rows) + 1L]] <- checkpoint$log
    if (nrow(checkpoint$curves)) curves[[length(curves) + 1L]] <- checkpoint$curves
    if (nrow(checkpoint$fixed)) fixed_rows[[length(fixed_rows) + 1L]] <- checkpoint$fixed
    if (nrow(checkpoint$joint)) joint_rows[[length(joint_rows) + 1L]] <- checkpoint$joint
    if (nrow(checkpoint$predictive)) predictive_rows[[length(predictive_rows) + 1L]] <- checkpoint$predictive
    write_checkpoint_csv(do.call(rbind, rows), file.path(out_dir, "nbi_sigma_compare_logs.csv"))
    if (length(curves) > 0L) {
      write_checkpoint_csv(do.call(rbind, curves), curve_path)
    }
    if (length(fixed_rows) > 0L) {
      write_checkpoint_csv(do.call(rbind, fixed_rows), fixed_path)
    }
    if (length(joint_rows) > 0L) {
      write_checkpoint_csv(do.call(rbind, joint_rows), joint_path)
    }
    if (length(predictive_rows) > 0L) {
      write_checkpoint_csv(do.call(rbind, predictive_rows), predictive_path)
    }
  }
}

log_df <- do.call(rbind, rows)
curve_df <- if (length(curves) > 0L) do.call(rbind, curves) else data.frame()
fixed_df <- if (length(fixed_rows) > 0L) do.call(rbind, fixed_rows) else data.frame()
joint_df <- if (length(joint_rows) > 0L) do.call(rbind, joint_rows) else data.frame()
predictive_df <- if (length(predictive_rows) > 0L) do.call(rbind, predictive_rows) else data.frame()
write_checkpoint_csv(log_df, log_path)
if (nrow(curve_df)) write_checkpoint_csv(curve_df, curve_path)
if (nrow(fixed_df)) write_checkpoint_csv(fixed_df, fixed_path)
if (nrow(joint_df)) write_checkpoint_csv(joint_df, joint_path)
if (nrow(predictive_df)) write_checkpoint_csv(predictive_df, predictive_path)
write_checkpoint_csv(
  data.frame(
    runner_contract_version = runner_contract_version,
    phase1_contract_version = phase1_contract_version,
    runner_settings_signature = runner_settings_signature,
    runner_settings_sha256 = runner_settings_sha256,
    runner_sha256 = runner_sha256,
    package_source_path = package_source_identity,
    package_source_sha256 = package_source_sha256,
    evidence_status = "post_phase1_production",
    margin_family = "NBI",
    copula_family = nbi_copula_label,
    copula_code = nbi_copula_code,
    n_subject = n_subject,
    times = paste(times, collapse = "|"),
    reps = reps,
    base_seed = base_seed,
    sigma_signal_multiplier = sigma_signal_multiplier,
    engines = paste(engines, collapse = "|"),
    max_elapsed_sec = max_elapsed_sec,
    max_outer_iter = max_outer_iter,
    max_inner_iter = max_inner_iter,
    start_step_size = start_step_size,
    step_adjustment = if (is.na(step_adjustment)) NA_real_ else step_adjustment,
    step_adjustment_env = step_adjustment_env,
    lambda_start = lambda_start,
    rs_update_lambda = rs_update_lambda,
    warm_start_joint_iter = warm_start_joint_iter,
    max_attempts_per_fit = max_attempts_per_fit,
    runtime_n_cores = 1L,
    runtime_backend = "sequential",
    rscript_path = nbi_rscript_path,
    rscript_sha256 = nbi_rscript_sha256,
    rscript_version = R.version.string,
    compute_se = compute_se,
    vcov_method_longitudinal = vcov_method_longitudinal,
    compute_predictive_scores = compute_predictive_scores,
    predictive_nsim = predictive_nsim,
    variogram_p = variogram_p,
    variogram_p_values = paste(variogram_p_values, collapse = "|"),
    save_fits = save_fits,
    theta_intercept = theta_intercept,
    theta_time_coef = theta_time_coef,
    theta_link = "log",
    theta_inverse_link = "exp",
    seed_rule = "base_seed + 1000*rep + offsets: covariate=11,response=12,test=13,diagnostic=100+engine,predictive=200+engine",
    package_version = runner_package_version,
    git_sha = runner_git_sha,
    git_state = runner_git_state,
    r_version = R.version.string,
    platform = R.version$platform,
    os = paste(Sys.info()[c("sysname", "release", "machine")], collapse = "|")
  ),
  file.path(out_dir, "nbi_sigma_compare_settings.csv")
)
summary_df <- aggregate(
  cbind(success = success, converged = converged, mu_irmse = mu_irmse, sigma_irmse = sigma_irmse, theta_irmse = theta_irmse, elapsed_sec = elapsed_sec) ~ engine,
  log_df,
  function(x) {
    if (is.logical(x)) sum(x, na.rm = TRUE) else mean(x, na.rm = TRUE)
  },
  na.action = NULL
)
write_checkpoint_csv(summary_df, file.path(out_dir, "nbi_sigma_compare_summary.csv"))

if (nrow(curve_df) > 0L) {
  smooth_integrated <- aggregate(
    error ~ engine + rep + retry_index + parameter,
    curve_df,
    function(x) sqrt(mean(x^2, na.rm = TRUE))
  )
  names(smooth_integrated) <- c("model", "rep", "retry_index", "parameter", "irmse")
  smooth_integrated$scenario <- sprintf("n%d_d%d_nbi_signal%s", n_subject, length(times), sigma_signal_multiplier)
  smooth_integrated$n <- n_subject
  smooth_integrated$d <- length(times)
  smooth_integrated$bias_abs_integrated <- aggregate(
    abs(error) ~ engine + rep + retry_index + parameter,
    curve_df,
    mean
  )$`abs(error)`
  smooth_integrated <- smooth_integrated[, c("scenario", "model", "n", "d", "parameter", "rep", "retry_index", "bias_abs_integrated", "irmse")]
  write_checkpoint_csv(smooth_integrated, file.path(out_dir, "smooth_integrated_metrics.csv"))

  smooth_estimates <- data.frame(
    scenario = sprintf("n%d_d%d_nbi_signal%s", n_subject, length(times), sigma_signal_multiplier),
    n = n_subject,
    d = length(times),
    rep = curve_df$rep,
    retry_index = curve_df$retry_index,
    model = curve_df$engine,
    parameter = curve_df$parameter,
    s1 = curve_df$x,
    smooth_hat = curve_df$fitted,
    smooth_true = curve_df$truth,
    stringsAsFactors = FALSE
  )
  write_checkpoint_csv(smooth_estimates, file.path(out_dir, "smooth_estimates_by_rep.csv"))
}

if (nrow(fixed_df) > 0L) {
  fixed_df$error <- fixed_df$estimate - fixed_df$true_value
  fixed_summary <- do.call(rbind, lapply(
    split(fixed_df, list(fixed_df$model, fixed_df$parameter, fixed_df$term), drop = TRUE),
    function(df) {
      ci_low <- df$estimate - stats::qnorm(0.975) * df$std_error
      ci_high <- df$estimate + stats::qnorm(0.975) * df$std_error
      covered <- is.finite(df$true_value) & is.finite(ci_low) & is.finite(ci_high) &
        df$true_value >= ci_low & df$true_value <= ci_high
      sd_est <- stats::sd(df$estimate, na.rm = TRUE)
      mean_se <- mean(df$std_error, na.rm = TRUE)
      data.frame(
        model = df$model[1],
        parameter = df$parameter[1],
        term = df$term[1],
        true_value = df$true_value[1],
        mean_estimate = mean(df$estimate, na.rm = TRUE),
        bias = mean(df$error, na.rm = TRUE),
        rmse = sqrt(mean(df$error^2, na.rm = TRUE)),
        sd_estimate = sd_est,
        mean_std_error = mean_se,
        se_to_empirical_sd = mean_se / sd_est,
        coverage_95 = mean(covered, na.rm = TRUE),
        n_successful_fits = sum(is.finite(df$estimate)),
        n_se = sum(is.finite(df$std_error)),
        stringsAsFactors = FALSE
      )
    }
  ))
  write_checkpoint_csv(fixed_summary, file.path(out_dir, "fixed_effects_bias_rmse_table.csv"))
  se_calibration <- fixed_summary[fixed_summary$term != "intercept", c(
    "model", "parameter", "term", "mean_std_error", "sd_estimate",
    "se_to_empirical_sd", "coverage_95", "n_se"
  )]
  write_checkpoint_csv(se_calibration, file.path(out_dir, "fixed_effects_se_calibration.csv"))
}

if (nrow(joint_df) > 0L) {
  joint_summary <- do.call(rbind, lapply(
    split(joint_df, joint_df$model, drop = TRUE),
    function(df) {
      data.frame(
        model = df$model[1],
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
  write_checkpoint_csv(joint_summary, file.path(out_dir, "joint_distribution_metrics_summary.csv"))
}

if (nrow(predictive_df) > 0L) {
  predictive_summary <- do.call(rbind, lapply(
    split(predictive_df, list(predictive_df$model, predictive_df$variogram_p), drop = TRUE),
    function(df) {
      data.frame(
        model = df$model[1],
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
  write_checkpoint_csv(predictive_summary, file.path(out_dir, "predictive_scores_summary.csv"))
}

if (length(curves) > 0L) {
  curve_df$bin_x <- round(curve_df$x, 2)
  curve_summary <- aggregate(
    cbind(fitted, truth) ~ engine + parameter + bin_x,
    curve_df,
    function(x) c(mean = mean(x, na.rm = TRUE), q25 = stats::quantile(x, 0.25, na.rm = TRUE), q75 = stats::quantile(x, 0.75, na.rm = TRUE))
  )
  curve_plot_df <- data.frame(
    engine = curve_summary$engine,
    parameter = curve_summary$parameter,
    x = curve_summary$bin_x,
    fitted = curve_summary$fitted[, "mean"],
    fitted_q25 = curve_summary$fitted[, 2],
    fitted_q75 = curve_summary$fitted[, 3],
    truth = curve_summary$truth[, "mean"],
    stringsAsFactors = FALSE
  )
  write_checkpoint_csv(curve_plot_df, file.path(out_dir, "smooth_pointwise_summary.csv"))

  p <- ggplot(curve_plot_df, aes(x = x)) +
    geom_ribbon(aes(ymin = fitted_q25, ymax = fitted_q75), fill = "#56B4E9", alpha = 0.28, colour = NA) +
    geom_line(aes(y = truth), linetype = 2, linewidth = 0.75) +
    geom_line(aes(y = fitted), colour = "#D55E00", linewidth = 0.85) +
    facet_grid(engine ~ parameter) +
    labs(
      title = "NBI smooth recovery",
      x = "s1",
      y = "centred linear-predictor contribution"
    ) +
    theme_bw(base_size = 11)
  ggsave(file.path(out_dir, "nbi_sigma_compare_recovery.png"), p, width = 11, height = 7, dpi = 160)
}

print(summary_df)
