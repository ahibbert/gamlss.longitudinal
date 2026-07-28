#!/usr/bin/env Rscript

suppressPackageStartupMessages({
library(devtools)
library(VineCopula)
library(gamlss)
library(gamlss.dist)
library(ggplot2)
})

devtools::load_all(".", quiet = TRUE)

`%||%` <- function(x, y) if (is.null(x)) y else x

n_subject <- as.integer(Sys.getenv("NBI_COMPARE_N", unset = "500"))
reps <- as.integer(Sys.getenv("NBI_COMPARE_REPS", unset = "1"))
base_seed <- as.integer(Sys.getenv("NBI_COMPARE_SEED", unset = "20260529"))
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
fit_dir <- file.path(out_dir, "fits")
if (isTRUE(save_fits)) {
  dir.create(fit_dir, recursive = TRUE, showWarnings = FALSE)
}

save_fit_object <- function(fit, rep_id, engine) {
  if (!isTRUE(save_fits)) {
    return(invisible(NULL))
  }
  file <- file.path(fit_dir, sprintf("fit_rep%03d_%s.rds", rep_id, engine))
  saveRDS(fit, file = file, compress = "gzip")
  invisible(file)
}

smooth_mu <- function(x) 0.55 * sin(2 * pi * x)
smooth_sigma <- function(x) (0.70 * sigma_signal_multiplier) * (2 * (x - 0.5)^2 - 1 / 6)
smooth_theta <- function(x) 0.22 * sin(pi * x) - 0.10

make_dat <- function(rep_id, response_seed_offset = 0L, covariates = NULL) {
  simulate_longitudinal_dataset(
    n = n_subject,
    times = times,
    margin_dist = gamlss.dist::NBI(),
    copula_dist = "N",
    covariates = covariates %||% function(base) {
      simulate_longitudinal_covariates(
        base,
        subject = list(
          x1 = function(d) stats::rnorm(nrow(d)),
          x2 = function(d) stats::rbinom(nrow(d), 1, 0.5),
          s1 = function(d) stats::runif(nrow(d))
        ),
        observation = list(
          time_scaled = function(d) sim_rescale01(d$time)
        )
      )
    },
    margin_params = list(
      mu = function(d) exp(1.10 + 0.18 * d$x1 - 0.12 * d$x2 + 0.10 * d$time_scaled + smooth_mu(d$s1)),
      sigma = function(d) exp(-0.45 + 0.10 * d$x1 + 0.08 * d$time_scaled + smooth_sigma(d$s1))
    ),
    copula_params = list(
      theta = function(e) {
        time_left_scaled <- sim_rescale01(e$time_left)
        tanh(theta_intercept + theta_time_coef * time_left_scaled + smooth_theta(e$s1))
      }
    ),
    seed = base_seed + rep_id + response_seed_offset
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
    copula_dist = "N",
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
if (isTRUE(resume) && file.exists(log_path)) {
  old_logs <- read.csv(log_path, stringsAsFactors = FALSE)
  rows <- split(old_logs, seq_len(nrow(old_logs)))
}
if (isTRUE(resume) && file.exists(curve_path)) {
  old_curves <- read.csv(curve_path, stringsAsFactors = FALSE)
  curves <- split(old_curves, seq_len(nrow(old_curves)))
}
if (isTRUE(resume) && file.exists(fixed_path)) {
  old_fixed <- read.csv(fixed_path, stringsAsFactors = FALSE)
  fixed_rows <- split(old_fixed, seq_len(nrow(old_fixed)))
}
if (isTRUE(resume) && file.exists(joint_path)) {
  old_joint <- read.csv(joint_path, stringsAsFactors = FALSE)
  joint_rows <- split(old_joint, seq_len(nrow(old_joint)))
}
if (isTRUE(resume) && file.exists(predictive_path)) {
  old_predictive <- read.csv(predictive_path, stringsAsFactors = FALSE)
  predictive_rows <- split(old_predictive, seq_len(nrow(old_predictive)))
}

already_done <- function(rep_id, engine) {
  if (!isTRUE(resume) || length(rows) == 0L) return(FALSE)
  log_df <- do.call(rbind, rows)
  any(log_df$rep == rep_id & log_df$engine == engine & log_df$success)
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

joint_metrics_gamlss <- function(fit, dat, rep_id) {
  pred_data <- gamlss_prediction_data(dat)
  pred_mu <- as.numeric(predict(fit, what = "mu", type = "link", newdata = pred_data))
  pred_sigma <- as.numeric(predict(fit, what = "sigma", type = "link", newdata = pred_data))
  u <- gamlss.dist::pNBI(
    dat$response,
    mu = gamlss.dist::NBI()$mu.linkinv(pred_mu),
    sigma = gamlss.dist::NBI()$sigma.linkinv(pred_sigma)
  )
  joint_metrics_from_rosenblatt("gamlss", rep_id, clip_u(u), dat, extract_fit_metric(fit, "logLik"), extract_fit_metric(fit, "df"))
}

joint_metrics_ours <- function(fit, rep_id, engine) {
  copula_link_fit <- gamlss.longitudinal:::get_copula_dist(fit$copula_dist)$copula_link
  eta_inv <- gamlss.longitudinal:::calc_eta(fit$par, fit$model_matrix, fit$margin_dist, copula_link_fit, fit$par_s)$eta_inv
  u <- gamlss.dist::pNBI(fit$response, mu = eta_inv$mu, sigma = eta_inv$sigma)
  u <- clip_u(u)
  dat <- data.frame(subject = fit$response_subject, time = fit$response_margin)
  pairs <- adjacent_pair_rows(dat)
  theta_rows <- which(fit$response_margin < max(fit$response_margin, na.rm = TRUE))
  theta_lookup <- rep(NA_integer_, length(fit$response))
  theta_lookup[theta_rows] <- seq_along(theta_rows)
  theta_idx <- theta_lookup[pairs$row1]
  theta <- as.numeric(eta_inv$theta[theta_idx])
  fam_num <- as.numeric(VineCopula::BiCopName(fit$copula_dist))
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

simulate_predictive_gamlss <- function(fit, dat, nsim) {
  dat <- dat[order(dat$subject, dat$time), , drop = FALSE]
  pred_data <- gamlss_prediction_data(dat)
  pred_mu <- as.numeric(predict(fit, what = "mu", type = "link", newdata = pred_data))
  pred_sigma <- as.numeric(predict(fit, what = "sigma", type = "link", newdata = pred_data))
  n <- length(unique(dat$subject))
  d <- length(unique(dat$time))
  out <- array(NA_real_, dim = c(nsim, n, d))
  for (s in seq_len(nsim)) {
    y <- gamlss.dist::qNBI(
      stats::runif(nrow(dat)),
      mu = gamlss.dist::NBI()$mu.linkinv(pred_mu),
      sigma = gamlss.dist::NBI()$sigma.linkinv(pred_sigma)
    )
    out[s, , ] <- matrix(y, nrow = n, ncol = d, byrow = TRUE)
  }
  out
}

simulate_predictive_ours <- function(fit, dat, nsim) {
  dat <- dat[order(dat$subject, dat$time), , drop = FALSE]
  sim_df <- stats::simulate(fit, nsim = nsim, newdata = dat)
  simulate_array_from_columns(sim_df, dat, subject_col = "subject")
}

predictive_scores_gamlss <- function(fit, test_dat, rep_id, nsim = predictive_nsim) {
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
  sim_array <- simulate_predictive_gamlss(fit, test_dat, nsim)
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

predictive_scores_ours <- function(fit, test_dat, rep_id, engine, nsim = predictive_nsim) {
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
  sim_array <- simulate_predictive_ours(fit, test_dat, nsim)
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
  dat <- make_dat(rep_id)
  fixed_covariates <- dat[, c("x1", "x2", "s1", "time_scaled"), drop = FALSE]
  test_dat <- make_dat(rep_id, response_seed_offset = 500000L, covariates = fixed_covariates)

  for (engine in engines) {
    if (already_done(rep_id, engine)) {
      cat(sprintf("rep %d/%d | %s already complete, skipping\n", rep_id, reps, engine))
      next
    }
    cat(sprintf("rep %d/%d | %s\n", rep_id, reps, engine))
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
      curves[[length(curves) + 1L]] <- curve
      fixed_rows[[length(fixed_rows) + 1L]] <- if (engine == "gamlss") {
        fixed_from_gamlss(fit, rep_id)
      } else {
        fixed_from_ours(fit, rep_id, engine)
      }
      joint_rows[[length(joint_rows) + 1L]] <- tryCatch(
        if (engine == "gamlss") joint_metrics_gamlss(fit, dat, rep_id) else joint_metrics_ours(fit, rep_id, engine),
        error = function(e) NULL
      )
      if (isTRUE(compute_predictive_scores)) {
        predictive_rows[[length(predictive_rows) + 1L]] <- tryCatch(
          if (engine == "gamlss") predictive_scores_gamlss(fit, test_dat, rep_id) else predictive_scores_ours(fit, test_dat, rep_id, engine),
          error = function(e) NULL
        )
      }
    }
    rows[[length(rows) + 1L]] <- data.frame(
      rep = rep_id,
      engine = engine,
      success = ok,
      logLik = if (ok) extract_fit_metric(fit, "logLik") else NA_real_,
      df = if (ok) extract_fit_metric(fit, "df") else NA_real_,
      converged = if (ok && engine == "gamlss") isTRUE(fit$converged) else if (ok) scalar_logical(fit$convergence$converged) else NA,
      iter = if (ok && engine == "gamlss") scalar_integer(fit$iter) else if (ok) scalar_integer(fit$convergence$outer_iterations) else NA_integer_,
      mu_irmse = if (ok && any(curve$parameter == "mu")) sqrt(mean(curve$error[curve$parameter == "mu"]^2, na.rm = TRUE)) else NA_real_,
      sigma_irmse = if (ok && any(curve$parameter == "sigma")) sqrt(mean(curve$error[curve$parameter == "sigma"]^2, na.rm = TRUE)) else NA_real_,
      theta_irmse = if (ok && any(curve$parameter == "theta")) sqrt(mean(curve$error[curve$parameter == "theta"]^2, na.rm = TRUE)) else NA_real_,
      elapsed_sec = elapsed,
      error = if (ok) NA_character_ else conditionMessage(fit),
      stringsAsFactors = FALSE
    )
    write.csv(do.call(rbind, rows), file.path(out_dir, "nbi_sigma_compare_logs.csv"), row.names = FALSE)
    if (length(curves) > 0L) {
      write.csv(do.call(rbind, curves), curve_path, row.names = FALSE)
    }
    if (length(fixed_rows) > 0L) {
      write.csv(do.call(rbind, fixed_rows), fixed_path, row.names = FALSE)
    }
    if (length(joint_rows) > 0L) {
      write.csv(do.call(rbind, joint_rows), joint_path, row.names = FALSE)
    }
    if (length(predictive_rows) > 0L) {
      write.csv(do.call(rbind, predictive_rows), predictive_path, row.names = FALSE)
    }
  }
}

log_df <- do.call(rbind, rows)
curve_df <- if (length(curves) > 0L) do.call(rbind, curves) else data.frame()
fixed_df <- if (length(fixed_rows) > 0L) do.call(rbind, fixed_rows) else data.frame()
joint_df <- if (length(joint_rows) > 0L) do.call(rbind, joint_rows) else data.frame()
predictive_df <- if (length(predictive_rows) > 0L) do.call(rbind, predictive_rows) else data.frame()
write.csv(
  data.frame(
    n_subject = n_subject,
    reps = reps,
    sigma_signal_multiplier = sigma_signal_multiplier,
    engines = paste(engines, collapse = "|"),
    max_elapsed_sec = max_elapsed_sec,
    max_outer_iter = max_outer_iter,
    max_inner_iter = max_inner_iter,
    start_step_size = start_step_size,
    step_adjustment = if (is.na(step_adjustment)) NA_real_ else step_adjustment,
    lambda_start = lambda_start,
    rs_update_lambda = rs_update_lambda,
    warm_start_joint_iter = warm_start_joint_iter,
    compute_se = compute_se,
    vcov_method_longitudinal = vcov_method_longitudinal,
    compute_predictive_scores = compute_predictive_scores,
    predictive_nsim = predictive_nsim,
    variogram_p = variogram_p,
    variogram_p_values = paste(variogram_p_values, collapse = "|"),
    save_fits = save_fits,
    theta_intercept = theta_intercept,
    theta_time_coef = theta_time_coef
  ),
  file.path(out_dir, "nbi_sigma_compare_settings.csv"),
  row.names = FALSE
)
summary_df <- aggregate(
  cbind(success = success, converged = converged, mu_irmse = mu_irmse, sigma_irmse = sigma_irmse, theta_irmse = theta_irmse, elapsed_sec = elapsed_sec) ~ engine,
  log_df,
  function(x) {
    if (is.logical(x)) sum(x, na.rm = TRUE) else mean(x, na.rm = TRUE)
  },
  na.action = NULL
)
write.csv(summary_df, file.path(out_dir, "nbi_sigma_compare_summary.csv"), row.names = FALSE)

if (nrow(curve_df) > 0L) {
  smooth_integrated <- aggregate(
    error ~ engine + rep + parameter,
    curve_df,
    function(x) sqrt(mean(x^2, na.rm = TRUE))
  )
  names(smooth_integrated) <- c("model", "rep", "parameter", "irmse")
  smooth_integrated$scenario <- sprintf("n%d_d%d_nbi_signal%s", n_subject, length(times), sigma_signal_multiplier)
  smooth_integrated$n <- n_subject
  smooth_integrated$d <- length(times)
  smooth_integrated$bias_abs_integrated <- aggregate(
    abs(error) ~ engine + rep + parameter,
    curve_df,
    mean
  )$`abs(error)`
  smooth_integrated <- smooth_integrated[, c("scenario", "model", "n", "d", "parameter", "rep", "bias_abs_integrated", "irmse")]
  write.csv(smooth_integrated, file.path(out_dir, "smooth_integrated_metrics.csv"), row.names = FALSE)

  smooth_estimates <- data.frame(
    scenario = sprintf("n%d_d%d_nbi_signal%s", n_subject, length(times), sigma_signal_multiplier),
    n = n_subject,
    d = length(times),
    rep = curve_df$rep,
    model = curve_df$engine,
    parameter = curve_df$parameter,
    s1 = curve_df$x,
    smooth_hat = curve_df$fitted,
    smooth_true = curve_df$truth,
    stringsAsFactors = FALSE
  )
  write.csv(smooth_estimates, file.path(out_dir, "smooth_estimates_by_rep.csv"), row.names = FALSE)
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
  write.csv(fixed_summary, file.path(out_dir, "fixed_effects_bias_rmse_table.csv"), row.names = FALSE)
  se_calibration <- fixed_summary[fixed_summary$term != "intercept", c(
    "model", "parameter", "term", "mean_std_error", "sd_estimate",
    "se_to_empirical_sd", "coverage_95", "n_se"
  )]
  write.csv(se_calibration, file.path(out_dir, "fixed_effects_se_calibration.csv"), row.names = FALSE)
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
  write.csv(joint_summary, file.path(out_dir, "joint_distribution_metrics_summary.csv"), row.names = FALSE)
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
  write.csv(predictive_summary, file.path(out_dir, "predictive_scores_summary.csv"), row.names = FALSE)
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
  write.csv(curve_plot_df, file.path(out_dir, "smooth_pointwise_summary.csv"), row.names = FALSE)

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
