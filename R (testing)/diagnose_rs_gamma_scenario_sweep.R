#!/usr/bin/env Rscript

# Scenario sweep for RS separate vs joint optimisation.
# Focus: continuous Gamma margin + Gaussian copula, no smooth terms.

safe_source <- function(path) {
  txt <- readLines(path, warn = FALSE, encoding = "UTF-8")
  if (length(txt) > 0) {
    bom <- intToUtf8(65279)
    txt[1] <- sub(paste0("^", bom), "", txt[1], useBytes = FALSE)
  }
  eval(parse(text = txt), envir = .GlobalEnv)
}

safe_source("R/common_functions.R")
safe_source("R/link_functions.R")

suppressPackageStartupMessages({
  library(gamlss)
  library(gamlss.dist)
  library(VineCopula)
})

out_dir <- file.path("results", "rs_gamma_scenario_sweep")
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

max_outer_iter <- max(2L, as.integer(Sys.getenv("RS_SWEEP_MAX_OUTER", unset = "35")))
max_inner_iter <- max(3L, as.integer(Sys.getenv("RS_SWEEP_MAX_INNER", unset = "12")))
verbose_fit <- as.integer(Sys.getenv("RS_SWEEP_VERBOSE", unset = "0"))

clip_u <- function(u) pmin(pmax(u, 1e-9), 1 - 1e-9)

scenarios <- list(
  baseline = list(
    n = 80, d = 4, sigma_base = 0.45, sigma_x = 0.12,
    theta_intercept = -0.55, theta_time = 0.28, theta_x = 0,
    theta_has_x = FALSE, missing = "none"
  ),
  strong_dependence = list(
    n = 80, d = 4, sigma_base = 0.45, sigma_x = 0.12,
    theta_intercept = 0.55, theta_time = 0.28, theta_x = 0,
    theta_has_x = FALSE, missing = "none"
  ),
  small_sample = list(
    n = 30, d = 4, sigma_base = 0.45, sigma_x = 0.12,
    theta_intercept = -0.55, theta_time = 0.28, theta_x = 0,
    theta_has_x = FALSE, missing = "none"
  ),
  skewed_margin = list(
    n = 80, d = 4, sigma_base = 0.90, sigma_x = 0.25,
    theta_intercept = -0.55, theta_time = 0.28, theta_x = 0,
    theta_has_x = FALSE, missing = "none"
  ),
  shared_sigma_theta_x = list(
    n = 80, d = 4, sigma_base = 0.60, sigma_x = 0.35,
    theta_intercept = -0.45, theta_time = 0.22, theta_x = 0.45,
    theta_has_x = TRUE, missing = "none"
  ),
  sharp_time_theta = list(
    n = 80, d = 5, sigma_base = 0.45, sigma_x = 0.12,
    theta_intercept = -1.05, theta_time = 0.55, theta_x = 0,
    theta_has_x = FALSE, missing = "none"
  ),
  uneven_missingness = list(
    n = 100, d = 4, sigma_base = 0.45, sigma_x = 0.12,
    theta_intercept = -0.55, theta_time = 0.28, theta_x = 0,
    theta_has_x = FALSE, missing = "late_time"
  )
)

simulate_gamma_gaussian <- function(cfg, seed) {
  set.seed(seed)
  n <- cfg$n
  d <- cfg$d

  subject_df <- data.frame(
    id = seq_len(n),
    x = rnorm(n),
    stringsAsFactors = FALSE
  )
  dat <- merge(
    subject_df,
    data.frame(time_index = seq_len(d), stringsAsFactors = FALSE),
    by = NULL,
    all = TRUE
  )
  dat <- dat[order(dat$id, dat$time_index), ]
  rownames(dat) <- NULL
  dat$t_scaled <- if (d > 1) (dat$time_index - 1) / (d - 1) else 0

  fam <- gamlss.dist::GA(mu.link = "log", sigma.link = "log")
  eta_mu <- 1.20 + 0.35 * dat$x + 0.25 * dat$t_scaled
  eta_sigma <- log(cfg$sigma_base) + cfg$sigma_x * dat$x
  mu <- fam$mu.linkinv(eta_mu)
  sigma <- fam$sigma.linkinv(eta_sigma)

  z <- matrix(NA_real_, nrow = n, ncol = d)
  z[, 1] <- rnorm(n)
  for (j in 2:d) {
    t_pair <- (j - 2) / max(1, d - 2)
    eta_theta_i <- cfg$theta_intercept + cfg$theta_time * t_pair + cfg$theta_x * subject_df$x
    rho_i <- tanh(eta_theta_i)
    z[, j] <- rho_i * z[, j - 1] + sqrt(pmax(1 - rho_i^2, 1e-8)) * rnorm(n)
  }

  dat$u_true <- clip_u(as.vector(t(pnorm(z))))
  dat$response <- gamlss.dist::qGA(dat$u_true, mu = mu, sigma = sigma)

  if (identical(cfg$missing, "late_time")) {
    p_miss <- c(0.02, 0.08, 0.20, 0.35, rep(0.35, max(0, d - 4)))
    p_miss <- p_miss[dat$time_index]
    dat$response[runif(nrow(dat)) < p_miss] <- NA_real_
  }

  dat
}

make_start <- function(dat, theta_has_x) {
  mu_start <- mean(dat$response, na.rm = TRUE)
  sigma_start <- stats::sd(dat$response, na.rm = TRUE) / mu_start
  start <- c(
    "mu.intercept" = log(mu_start),
    "mu.x" = 0,
    "mu.t_scaled" = 0,
    "sigma.intercept" = log(sigma_start),
    "sigma.x" = 0,
    "theta.intercept" = atanh(0.20),
    "theta.time_covariate" = 0
  )
  if (isTRUE(theta_has_x)) {
    start <- c(start, "theta.x" = 0)
  }
  start
}

fit_rs <- function(dat, cfg, include_dlcopdpar, start_from = NA) {
  theta_formula <- if (isTRUE(cfg$theta_has_x)) {
    ~ time_index + x
  } else {
    ~ time_index
  }
  if (all(is.na(start_from))) {
    start_from <- make_start(dat, cfg$theta_has_x)
  }

  fit <- NULL
  captured <- utils::capture.output({
    fit <- gamlss.longitudinal(
      dataset = dat,
      margin_dist = gamlss.dist::GA(mu.link = "log", sigma.link = "log"),
      copula_dist = "N",
      time_var = "time_index",
      subject_var = "id",
      mu.formula = response ~ x + t_scaled,
      sigma.formula = ~ x,
      theta.formula = theta_formula,
      include_dlcopdpar = include_dlcopdpar,
      start_from = start_from,
      method = "RS",
      compute_vcov = FALSE,
      max_outer_iter = max_outer_iter,
      max_inner_iter = max_inner_iter,
      inner_stop_crit = 1e-5,
      outer_stop_crit = 1e-5,
      start_step_size = 0.5,
      step_adjustment = 0.5,
      max_steps = 5,
      use_backtracking = TRUE,
      backtracking_max_halves = 50,
      max_negative_outer_streak = 20,
      verbose = verbose_fit,
      plot_results = FALSE
    )
  }, type = "output")
  attr(fit, "captured_output") <- captured
  fit
}

endpoint_row <- function(scenario_name, cfg, fit_name, fit) {
  ll <- fit$calc_lik_out_end$log_lik
  data.frame(
    scenario = scenario_name,
    fit = fit_name,
    n = cfg$n,
    d = cfg$d,
    sigma_base = cfg$sigma_base,
    sigma_x = cfg$sigma_x,
    theta_intercept = cfg$theta_intercept,
    theta_time = cfg$theta_time,
    theta_x = cfg$theta_x,
    missing = cfg$missing,
    marginal = as.numeric(ll["marginal"]),
    copula = as.numeric(ll["copula"]),
    joint = as.numeric(ll["joint"]),
    stringsAsFactors = FALSE
  )
}

coef_rows <- function(scenario_name, fits) {
  all_names <- unique(unlist(lapply(fits, function(fit) names(fit$par)), use.names = FALSE))
  rows <- lapply(names(fits), function(fit_name) {
    data.frame(
      scenario = scenario_name,
      fit = fit_name,
      parameter = all_names,
      estimate = as.numeric(fits[[fit_name]]$par[all_names]),
      stringsAsFactors = FALSE
    )
  })
  do.call(rbind, rows)
}

all_endpoints <- list()
all_coefs <- list()
all_traces <- list()

scenario_i <- 1L
for (scenario_name in names(scenarios)) {
  cfg <- scenarios[[scenario_name]]
  cat("\n=== Scenario: ", scenario_name, " ===\n", sep = "")
  dat <- simulate_gamma_gaussian(cfg, seed = 20260514 + scenario_i * 100)
  write.csv(dat, file.path(out_dir, paste0("data_", scenario_name, ".csv")), row.names = FALSE)

  fit_sep <- tryCatch(fit_rs(dat, cfg, include_dlcopdpar = FALSE), error = function(e) e)
  fit_joint <- tryCatch(fit_rs(dat, cfg, include_dlcopdpar = TRUE), error = function(e) e)

  if (inherits(fit_sep, "error") || inherits(fit_joint, "error")) {
    msg_sep <- if (inherits(fit_sep, "error")) conditionMessage(fit_sep) else ""
    msg_joint <- if (inherits(fit_joint, "error")) conditionMessage(fit_joint) else ""
    all_endpoints[[length(all_endpoints) + 1L]] <- data.frame(
      scenario = scenario_name,
      fit = c("separate", "joint", "joint_from_separate"),
      n = cfg$n,
      d = cfg$d,
      sigma_base = cfg$sigma_base,
      sigma_x = cfg$sigma_x,
      theta_intercept = cfg$theta_intercept,
      theta_time = cfg$theta_time,
      theta_x = cfg$theta_x,
      missing = cfg$missing,
      marginal = NA_real_,
      copula = NA_real_,
      joint = NA_real_,
      error = c(msg_sep, msg_joint, "not run"),
      stringsAsFactors = FALSE
    )
    scenario_i <- scenario_i + 1L
    next
  }

  fit_joint_from_sep <- tryCatch(
    fit_rs(dat, cfg, include_dlcopdpar = TRUE, start_from = fit_sep$par),
    error = function(e) e
  )

  fits <- list(separate = fit_sep, joint = fit_joint)
  if (!inherits(fit_joint_from_sep, "error")) {
    fits$joint_from_separate <- fit_joint_from_sep
  }

  writeLines(attr(fit_sep, "captured_output"), file.path(out_dir, paste0("fit_", scenario_name, "_separate.log")))
  writeLines(attr(fit_joint, "captured_output"), file.path(out_dir, paste0("fit_", scenario_name, "_joint.log")))
  if (!inherits(fit_joint_from_sep, "error")) {
    writeLines(attr(fit_joint_from_sep, "captured_output"), file.path(out_dir, paste0("fit_", scenario_name, "_joint_from_separate.log")))
  }

  rows <- do.call(rbind, lapply(names(fits), function(nm) endpoint_row(scenario_name, cfg, nm, fits[[nm]])))
  all_endpoints[[length(all_endpoints) + 1L]] <- rows
  all_coefs[[length(all_coefs) + 1L]] <- coef_rows(scenario_name, fits)

  trace_rows <- lapply(names(fits), function(nm) {
    lh <- as.matrix(fits[[nm]]$log_lik_history)
    if (nrow(lh) == 0) return(NULL)
    data.frame(
      scenario = scenario_name,
      fit = nm,
      step = seq_len(nrow(lh)),
      marginal = as.numeric(lh[, "marginal"]),
      copula = as.numeric(lh[, "copula"]),
      joint = as.numeric(lh[, "joint"]),
      stringsAsFactors = FALSE
    )
  })
  all_traces[[length(all_traces) + 1L]] <- do.call(rbind, trace_rows)

  print(rows[, c("scenario", "fit", "marginal", "copula", "joint")], row.names = FALSE)
  scenario_i <- scenario_i + 1L
}

endpoints <- do.call(rbind, all_endpoints)
coefs <- if (length(all_coefs) > 0) do.call(rbind, all_coefs) else data.frame()
traces <- if (length(all_traces) > 0) do.call(rbind, all_traces) else data.frame()

write.csv(endpoints, file.path(out_dir, "scenario_endpoint_loglik.csv"), row.names = FALSE)
write.csv(coefs, file.path(out_dir, "scenario_coefficients.csv"), row.names = FALSE)
write.csv(traces, file.path(out_dir, "scenario_likelihood_traces.csv"), row.names = FALSE)

valid_endpoints <- endpoints[is.finite(endpoints$joint), , drop = FALSE]
summary_rows <- lapply(split(valid_endpoints, valid_endpoints$scenario), function(x) {
  sep <- x[x$fit == "separate", , drop = FALSE]
  joint <- x[x$fit == "joint", , drop = FALSE]
  jfs <- x[x$fit == "joint_from_separate", , drop = FALSE]
  if (nrow(sep) == 0 || nrow(joint) == 0) return(NULL)
  data.frame(
    scenario = x$scenario[1],
    separate_joint = sep$joint[1],
    joint_joint = joint$joint[1],
    joint_from_separate_joint = if (nrow(jfs) > 0) jfs$joint[1] else NA_real_,
    delta_joint_vs_separate = joint$joint[1] - sep$joint[1],
    delta_jfs_vs_separate = if (nrow(jfs) > 0) jfs$joint[1] - sep$joint[1] else NA_real_,
    delta_jfs_vs_joint = if (nrow(jfs) > 0) jfs$joint[1] - joint$joint[1] else NA_real_,
    delta_joint_pct_abs = 100 * (joint$joint[1] - sep$joint[1]) / abs(sep$joint[1]),
    delta_jfs_pct_abs = if (nrow(jfs) > 0) 100 * (jfs$joint[1] - sep$joint[1]) / abs(sep$joint[1]) else NA_real_,
    stringsAsFactors = FALSE
  )
})
summary_df <- do.call(rbind, summary_rows)
summary_df <- summary_df[order(abs(summary_df$delta_jfs_vs_separate), decreasing = TRUE), ]
write.csv(summary_df, file.path(out_dir, "scenario_delta_summary.csv"), row.names = FALSE)

cat("\n=== Delta summary ===\n")
print(summary_df, row.names = FALSE)
