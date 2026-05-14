#!/usr/bin/env Rscript

# Focused diagnostic for the scenario that amplified joint-vs-separate
# differences most in diagnose_rs_gamma_scenario_sweep.R:
# a Gamma margin where x drives both sigma and Gaussian-copula theta.

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

out_dir <- file.path("results", "rs_shared_sigma_theta_x_focus")
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

max_outer_iter <- max(2L, as.integer(Sys.getenv("RS_FOCUS_MAX_OUTER", unset = "45")))
max_inner_iter <- max(3L, as.integer(Sys.getenv("RS_FOCUS_MAX_INNER", unset = "14")))
verbose_fit <- as.integer(Sys.getenv("RS_FOCUS_VERBOSE", unset = "0"))

variants <- list(
  base = list(n = 80, d = 4, sigma_base = 0.60, sigma_x = 0.35, theta_intercept = -0.45, theta_time = 0.22, theta_x = 0.45),
  stronger_x = list(n = 80, d = 4, sigma_base = 0.70, sigma_x = 0.55, theta_intercept = -0.45, theta_time = 0.22, theta_x = 0.70),
  strong_dep_x = list(n = 80, d = 4, sigma_base = 0.60, sigma_x = 0.45, theta_intercept = 0.15, theta_time = 0.22, theta_x = 0.65),
  small_strong_x = list(n = 35, d = 4, sigma_base = 0.70, sigma_x = 0.55, theta_intercept = -0.45, theta_time = 0.22, theta_x = 0.70)
)

clip_u <- function(u) pmin(pmax(u, 1e-9), 1 - 1e-9)

simulate_dat <- function(cfg, seed) {
  set.seed(seed)
  subject_df <- data.frame(id = seq_len(cfg$n), x = rnorm(cfg$n))
  dat <- merge(subject_df, data.frame(time_index = seq_len(cfg$d)), by = NULL, all = TRUE)
  dat <- dat[order(dat$id, dat$time_index), ]
  rownames(dat) <- NULL
  dat$t_scaled <- if (cfg$d > 1) (dat$time_index - 1) / (cfg$d - 1) else 0

  fam <- GA(mu.link = "log", sigma.link = "log")
  mu <- fam$mu.linkinv(1.20 + 0.35 * dat$x + 0.25 * dat$t_scaled)
  sigma <- fam$sigma.linkinv(log(cfg$sigma_base) + cfg$sigma_x * dat$x)

  z <- matrix(NA_real_, nrow = cfg$n, ncol = cfg$d)
  z[, 1] <- rnorm(cfg$n)
  for (j in 2:cfg$d) {
    t_pair <- (j - 2) / max(1, cfg$d - 2)
    eta_theta_i <- cfg$theta_intercept + cfg$theta_time * t_pair + cfg$theta_x * subject_df$x
    rho_i <- tanh(eta_theta_i)
    z[, j] <- rho_i * z[, j - 1] + sqrt(pmax(1 - rho_i^2, 1e-8)) * rnorm(cfg$n)
  }

  dat$response <- qGA(clip_u(as.vector(t(pnorm(z)))), mu = mu, sigma = sigma)
  dat
}

make_start <- function(dat) {
  mu_start <- mean(dat$response, na.rm = TRUE)
  sigma_start <- sd(dat$response, na.rm = TRUE) / mu_start
  c(
    "mu.intercept" = log(mu_start),
    "mu.x" = 0,
    "mu.t_scaled" = 0,
    "sigma.intercept" = log(sigma_start),
    "sigma.x" = 0,
    "theta.intercept" = atanh(0.20),
    "theta.time_covariate" = 0,
    "theta.x" = 0
  )
}

fit_rs <- function(dat, include_dlcopdpar, start_from = NA) {
  if (all(is.na(start_from))) start_from <- make_start(dat)
  fit <- NULL
  captured <- capture.output({
    fit <- gamlss.longitudinal(
      dataset = dat,
      margin_dist = GA(mu.link = "log", sigma.link = "log"),
      copula_dist = "N",
      time_var = "time_index",
      subject_var = "id",
      mu.formula = response ~ x + t_scaled,
      sigma.formula = ~ x,
      theta.formula = ~ time_index + x,
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

evaluate_coefficients <- function(fit, par_vec) {
  copula_link <- get_copula_dist(fit$copula_dist)$copula_link
  eta_out <- calc_eta(par_vec, fit$model_matrix, fit$margin_dist, copula_link, fit$par_s)
  pair_cache <- build_copula_pair_cache(fit$response, fit$response_margin, fit$response_subject)
  calc_likelihood_minimal(
    eta_inv = eta_out$eta_inv,
    mm = fit$model_matrix$x,
    margin_dist = fit$margin_dist,
    copula_dist = fit$copula_dist,
    calc_d2 = FALSE,
    response = fit$response,
    response_margin = fit$response_margin,
    response_subject = fit$response_subject,
    pair_cache = pair_cache
  )$log_lik
}

swap_par <- function(base, donor, prefixes) {
  out <- base
  idx <- Reduce(`|`, lapply(prefixes, function(p) startsWith(names(out), paste0(p, "."))))
  out[idx] <- donor[names(out)[idx]]
  out
}

endpoint_row <- function(variant, label, fit_ref, par_vec) {
  ll <- evaluate_coefficients(fit_ref, par_vec)
  data.frame(
    variant = variant,
    endpoint = label,
    marginal = as.numeric(ll["marginal"]),
    copula = as.numeric(ll["copula"]),
    joint = as.numeric(ll["joint"]),
    stringsAsFactors = FALSE
  )
}

all_endpoints <- list()
all_coefs <- list()
all_traces <- list()

i <- 1L
for (variant_name in names(variants)) {
  cfg <- variants[[variant_name]]
  cat("\n=== Focus variant: ", variant_name, " ===\n", sep = "")
  dat <- simulate_dat(cfg, seed = 20260514 + i * 1000)
  write.csv(dat, file.path(out_dir, paste0("data_", variant_name, ".csv")), row.names = FALSE)

  fit_sep <- fit_rs(dat, include_dlcopdpar = FALSE)
  fit_joint <- fit_rs(dat, include_dlcopdpar = TRUE)
  fit_jfs <- fit_rs(dat, include_dlcopdpar = TRUE, start_from = fit_sep$par)

  writeLines(attr(fit_sep, "captured_output"), file.path(out_dir, paste0("fit_", variant_name, "_separate.log")))
  writeLines(attr(fit_joint, "captured_output"), file.path(out_dir, paste0("fit_", variant_name, "_joint.log")))
  writeLines(attr(fit_jfs, "captured_output"), file.path(out_dir, paste0("fit_", variant_name, "_joint_from_separate.log")))

  par_sep <- fit_sep$par
  par_joint <- fit_joint$par
  par_jfs <- fit_jfs$par

  endpoints <- rbind(
    endpoint_row(variant_name, "separate", fit_sep, par_sep),
    endpoint_row(variant_name, "joint", fit_sep, par_joint),
    endpoint_row(variant_name, "joint_from_separate", fit_sep, par_jfs),
    endpoint_row(variant_name, "sep_mu_jfs_sigma_theta", fit_sep, swap_par(par_sep, par_jfs, c("sigma", "theta"))),
    endpoint_row(variant_name, "sep_sigma_theta_jfs_mu", fit_sep, swap_par(par_sep, par_jfs, c("mu"))),
    endpoint_row(variant_name, "sep_theta_jfs_mu_sigma", fit_sep, swap_par(par_sep, par_jfs, c("mu", "sigma"))),
    endpoint_row(variant_name, "sep_mu_sigma_jfs_theta", fit_sep, swap_par(par_sep, par_jfs, c("theta")))
  )
  endpoints$delta_vs_separate <- endpoints$joint - endpoints$joint[endpoints$endpoint == "separate"]
  all_endpoints[[length(all_endpoints) + 1L]] <- endpoints

  coef_names <- unique(c(names(par_sep), names(par_joint), names(par_jfs)))
  all_coefs[[length(all_coefs) + 1L]] <- data.frame(
    variant = variant_name,
    parameter = coef_names,
    separate = as.numeric(par_sep[coef_names]),
    joint = as.numeric(par_joint[coef_names]),
    joint_from_separate = as.numeric(par_jfs[coef_names]),
    stringsAsFactors = FALSE
  )

  traces <- lapply(
    list(separate = fit_sep, joint = fit_joint, joint_from_separate = fit_jfs),
    function(fit) {
      lh <- as.matrix(fit$log_lik_history)
      data.frame(
        variant = variant_name,
        step = seq_len(nrow(lh)),
        marginal = as.numeric(lh[, "marginal"]),
        copula = as.numeric(lh[, "copula"]),
        joint = as.numeric(lh[, "joint"]),
        stringsAsFactors = FALSE
      )
    }
  )
  traces <- do.call(rbind, Map(function(nm, x) transform(x, fit = nm), names(traces), traces))
  all_traces[[length(all_traces) + 1L]] <- traces

  print(endpoints[, c("variant", "endpoint", "marginal", "copula", "joint", "delta_vs_separate")], row.names = FALSE)
  i <- i + 1L
}

endpoint_df <- do.call(rbind, all_endpoints)
coef_df <- do.call(rbind, all_coefs)
trace_df <- do.call(rbind, all_traces)

write.csv(endpoint_df, file.path(out_dir, "focused_endpoint_decomposition.csv"), row.names = FALSE)
write.csv(coef_df, file.path(out_dir, "focused_coefficients.csv"), row.names = FALSE)
write.csv(trace_df, file.path(out_dir, "focused_likelihood_traces.csv"), row.names = FALSE)

summary_df <- endpoint_df[endpoint_df$endpoint %in% c("separate", "joint", "joint_from_separate"), ]
summary_wide <- reshape(
  summary_df[, c("variant", "endpoint", "joint")],
  idvar = "variant",
  timevar = "endpoint",
  direction = "wide"
)
names(summary_wide) <- sub("^joint\\.", "", names(summary_wide))
summary_wide$delta_joint_vs_separate <- summary_wide$joint - summary_wide$separate
summary_wide$delta_jfs_vs_separate <- summary_wide$joint_from_separate - summary_wide$separate
summary_wide$delta_jfs_vs_joint <- summary_wide$joint_from_separate - summary_wide$joint
summary_wide <- summary_wide[order(abs(summary_wide$delta_jfs_vs_separate), decreasing = TRUE), ]
write.csv(summary_wide, file.path(out_dir, "focused_delta_summary.csv"), row.names = FALSE)

cat("\n=== Focus summary ===\n")
print(summary_wide, row.names = FALSE)
