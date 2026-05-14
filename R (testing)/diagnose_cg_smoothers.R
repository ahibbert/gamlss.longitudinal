#!/usr/bin/env Rscript

args <- commandArgs(trailingOnly = TRUE)
if (length(args) < 1L) {
  stop("Usage: diagnose_cg_smoothers.R <result_dir> [<result_dir> ...]")
}

smooth_truth <- list(
  mu = function(s1) 0.55 * sin(2 * pi * s1),
  sigma = function(s1) 0.30 * cos(2 * pi * s1) - 0.10 * (s1 - 0.50)^2,
  theta = function(s1) 2 * (0.25 * sin(pi * s1 + 0.40) + 0.12 * (s1 - 0.50))
)

center_curve <- function(y) {
  y - mean(y, na.rm = TRUE)
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
    if (parameter %in% c("theta", "zeta")) {
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

smooth_penalty <- function(B, beta_s, lambda) {
  S <- attr(B, "penalty")
  if (is.null(S) || !is.matrix(S) || nrow(S) != ncol(B) || ncol(S) != ncol(B)) {
    if (ncol(B) >= 3) {
      D <- diff(diag(ncol(B)), differences = 2)
      S <- t(D) %*% D
    } else {
      S <- diag(ncol(B))
    }
  }
  as.numeric(lambda) * as.numeric(t(beta_s) %*% S %*% beta_s)
}

extract_smooth_diag <- function(label, fit_obj, data_used) {
  rows <- list()
  s1_grid <- seq(0, 1, length.out = 101)
  for (parameter in names(smooth_truth)) {
    B_list <- fit_obj$model_matrix$s[[parameter]]
    par_s_list <- fit_obj$par_s[[parameter]]
    if (is.null(B_list) || length(B_list) == 0) next
    for (s_name in names(B_list)) {
      B <- B_list[[s_name]]
      beta_s <- par_s_list[[s_name]]
      lambda <- fit_obj$lambda_s[[parameter]][[s_name]]
      df <- fit_obj$df_s[[parameter]][[s_name]]
      smooth_hat <- smooth_curve_longitudinal(fit_obj, data_used, parameter, s1_grid)
      smooth_true <- center_curve(smooth_truth[[parameter]](s1_grid))
      err <- smooth_hat - smooth_true
      rows[[length(rows) + 1L]] <- data.frame(
        run = label,
        optim_method = fit_obj$optim_method,
        include_dlcopdpar = fit_obj$include_dlcopdpar,
        parameter = parameter,
        smooth = s_name,
        lambda = as.numeric(lambda),
        df_s = suppressWarnings(as.numeric(df)),
        beta_l2 = sqrt(sum(as.numeric(beta_s)^2)),
        beta_max_abs = max(abs(as.numeric(beta_s))),
        penalty = smooth_penalty(B, beta_s, lambda),
        irmse = sqrt(mean(err^2, na.rm = TRUE)),
        iabs_bias = mean(abs(err), na.rm = TRUE),
        stringsAsFactors = FALSE
      )
    }
  }
  do.call(rbind, rows)
}

extract_loglik_diag <- function(label, fit_obj) {
  hist <- as.data.frame(fit_obj$log_lik_history)
  if (nrow(hist) == 0L) return(NULL)
  names(hist) <- c("copula", "margin", "joint")[seq_len(ncol(hist))]
  hist$iteration <- seq_len(nrow(hist))
  hist$run <- label
  hist$optim_method <- fit_obj$optim_method
  hist$include_dlcopdpar <- fit_obj$include_dlcopdpar
  hist
}

all_smooth <- list()
all_loglik <- list()

for (dir in args) {
  rep_files <- list.files(file.path(dir, "rep_results"), pattern = "\\.rds$", full.names = TRUE)
  if (length(rep_files) == 0L) next
  obj <- readRDS(rep_files[1])
  if (is.null(obj$fit_longitudinal)) {
    warning("No saved fit found in ", dir)
    next
  }
  label <- basename(normalizePath(dir, winslash = "/", mustWork = FALSE))
  all_smooth[[length(all_smooth) + 1L]] <- extract_smooth_diag(label, obj$fit_longitudinal, obj$data)
  all_loglik[[length(all_loglik) + 1L]] <- extract_loglik_diag(label, obj$fit_longitudinal)
}

smooth_out <- do.call(rbind, all_smooth)
loglik_out <- do.call(rbind, all_loglik)

out_dir <- file.path("results", "cg_smoother_diagnostics")
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)
write.csv(smooth_out, file.path(out_dir, "smooth_diagnostics.csv"), row.names = FALSE)
write.csv(loglik_out, file.path(out_dir, "loglik_histories.csv"), row.names = FALSE)

print(smooth_out)
if (!is.null(loglik_out)) {
  loglik_summary <- do.call(rbind, lapply(split(loglik_out, loglik_out$run), function(x) {
    data.frame(
      run = x$run[1],
      n_iter = nrow(x),
      start_joint = x$joint[1],
      end_joint = x$joint[nrow(x)],
      joint_gain = x$joint[nrow(x)] - x$joint[1],
      stringsAsFactors = FALSE
    )
  }))
  write.csv(loglik_summary, file.path(out_dir, "loglik_summary.csv"), row.names = FALSE)
  print(loglik_summary)
}
