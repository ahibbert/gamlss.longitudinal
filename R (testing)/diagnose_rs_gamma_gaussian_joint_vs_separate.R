#!/usr/bin/env Rscript

# Diagnostic comparison of RS optimisation with and without dlcopdpar.
# The design is deliberately simple:
#   - continuous Gamma margin
#   - Gaussian pair copula
#   - no smooth terms
#   - theta has an intercept and a time coefficient
#
# Outputs are written to:
#   results/rs_gamma_gaussian_joint_vs_separate/

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

set.seed(20260514)

out_dir <- file.path("results", "rs_gamma_gaussian_joint_vs_separate")
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)
progress_file <- file.path(out_dir, "progress.log")
writeLines(character(0), progress_file)

mark_progress <- function(...) {
  msg <- paste0(format(Sys.time(), "%Y-%m-%d %H:%M:%S"), " | ", paste(..., collapse = ""))
  cat(msg, "\n", sep = "")
  flush.console()
  cat(msg, "\n", file = progress_file, append = TRUE, sep = "")
}

n_subject <- as.integer(Sys.getenv("RS_DIAG_N", unset = "80"))
n_time <- as.integer(Sys.getenv("RS_DIAG_D", unset = "4"))
max_outer_iter <- max(2L, as.integer(Sys.getenv("RS_DIAG_MAX_OUTER", unset = "20")))
max_inner_iter <- max(3L, as.integer(Sys.getenv("RS_DIAG_MAX_INNER", unset = "10")))
verbose_fit <- as.integer(Sys.getenv("RS_DIAG_VERBOSE", unset = "0"))

clip_u <- function(u) pmin(pmax(u, 1e-9), 1 - 1e-9)

simulate_gamma_gaussian <- function(n, d, seed = 20260514) {
  set.seed(seed)

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
  eta_sigma <- log(0.45) + 0.12 * dat$x
  mu <- fam$mu.linkinv(eta_mu)
  sigma <- fam$sigma.linkinv(eta_sigma)

  theta_eta_by_pair <- -0.55 + 0.28 * seq_len(d - 1)
  rho_by_pair <- tanh(theta_eta_by_pair)

  z <- matrix(NA_real_, nrow = n, ncol = d)
  z[, 1] <- rnorm(n)
  if (d > 1) {
    for (j in 2:d) {
      rho <- rho_by_pair[j - 1]
      z[, j] <- rho * z[, j - 1] + sqrt(1 - rho^2) * rnorm(n)
    }
  }

  u <- as.vector(t(pnorm(z)))
  dat$u_true <- clip_u(u)
  dat$response <- gamlss.dist::qGA(dat$u_true, mu = mu, sigma = sigma)
  dat$mu_true <- mu
  dat$sigma_true <- sigma
  dat$rho_true <- c(rho_by_pair, NA_real_)[dat$time_index]
  dat$theta_eta_true <- c(theta_eta_by_pair, NA_real_)[dat$time_index]

  dat[, c(
    "id", "time_index", "t_scaled", "x", "response",
    "mu_true", "sigma_true", "rho_true", "theta_eta_true", "u_true"
  )]
}

fit_rs <- function(dat, include_dlcopdpar, start_from = NA) {
  if (all(is.na(start_from))) {
    mu_start <- mean(dat$response, na.rm = TRUE)
    sigma_start <- stats::sd(dat$response, na.rm = TRUE) / mu_start
    start_from <- c(
      "mu.intercept" = log(mu_start),
      "mu.x" = 0,
      "mu.t_scaled" = 0,
      "sigma.intercept" = log(sigma_start),
      "sigma.x" = 0,
      "theta.intercept" = atanh(0.20),
      "theta.time_covariate" = 0
    )
  }

  fit_expr <- quote(gamlss.longitudinal(
    dataset = dat,
    margin_dist = gamlss.dist::GA(mu.link = "log", sigma.link = "log"),
    copula_dist = "N",
    time_var = "time_index",
    subject_var = "id",
    mu.formula = response ~ x + t_scaled,
    sigma.formula = ~ x,
    theta.formula = ~ time_index,
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
  ))

  fit <- NULL
  fit_output <- utils::capture.output({
    fit <- eval(fit_expr)
  }, type = "output")

  attr(fit, "captured_output") <- fit_output
  fit
}

evaluate_coefficients <- function(fit, par_vec) {
  copula_link <- get_copula_dist(fit$copula_dist)$copula_link
  eta_out <- calc_eta(
    par_cov = par_vec,
    mm = fit$model_matrix,
    margin_dist = fit$margin_dist,
    copula_link = copula_link,
    par_s = fit$par_s
  )
  pair_cache <- build_copula_pair_cache(
    response = fit$response,
    response_margin = fit$response_margin,
    response_subject = fit$response_subject
  )
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
  )
}

changed_prefix <- function(before, after, tol = 1e-10) {
  if (length(before) != length(after)) return(NA_character_)
  dif <- abs(after - before)
  changed <- names(dif)[is.finite(dif) & dif > tol]
  if (length(changed) == 0) return("")
  prefixes <- unique(sub("\\..*$", "", changed))
  paste(prefixes, collapse = "+")
}

trace_from_fit_history <- function(fit, label) {
  lh <- as.matrix(fit$log_lik_history)
  if (nrow(lh) == 0) {
    return(data.frame())
  }

  trace <- data.frame(
    fit = label,
    step = seq_len(nrow(lh)),
    marginal = as.numeric(lh[, "marginal"]),
    copula = as.numeric(lh[, "copula"]),
    joint = as.numeric(lh[, "joint"]),
    stringsAsFactors = FALSE
  )
  trace$delta_joint <- c(NA_real_, diff(trace$joint))
  trace$delta_marginal <- c(NA_real_, diff(trace$marginal))
  trace$delta_copula <- c(NA_real_, diff(trace$copula))
  trace
}

replay_path_from_par_history <- function(fit, label) {
  ph <- as.matrix(fit$par_history)
  if (nrow(ph) == 0) {
    return(data.frame())
  }

  rows <- vector("list", nrow(ph))
  for (i in seq_len(nrow(ph))) {
    par_i <- as.numeric(ph[i, ])
    names(par_i) <- colnames(ph)
    lik_i <- evaluate_coefficients(fit, par_i)$log_lik

    rows[[i]] <- data.frame(
      fit = label,
      step = i,
      marginal = as.numeric(lik_i["marginal"]),
      copula = as.numeric(lik_i["copula"]),
      joint = as.numeric(lik_i["joint"]),
      changed_since_previous = if (i == 1) "" else changed_prefix(ph[i - 1, ], ph[i, ]),
      stringsAsFactors = FALSE
    )
  }

  out <- do.call(rbind, rows)
  out$delta_joint <- c(NA_real_, diff(out$joint))
  out$delta_marginal <- c(NA_real_, diff(out$marginal))
  out$delta_copula <- c(NA_real_, diff(out$copula))
  out
}

coefficient_table <- function(fits) {
  all_names <- unique(unlist(lapply(fits, function(fit) names(fit$par)), use.names = FALSE))
  out <- data.frame(
    parameter = all_names,
    stringsAsFactors = FALSE
  )
  for (fit_name in names(fits)) {
    out[[fit_name]] <- as.numeric(fits[[fit_name]]$par[all_names])
  }
  if (all(c("separate", "joint") %in% names(fits))) {
    out$diff_joint_minus_separate <- out$joint - out$separate
  }
  if (all(c("separate", "joint_from_separate") %in% names(fits))) {
    out$diff_joint_from_separate_minus_separate <- out$joint_from_separate - out$separate
  }
  out
}

endpoint_rows <- function(fit_ref, par_list) {
  rows <- lapply(names(par_list), function(label) {
    lik <- evaluate_coefficients(fit_ref, par_list[[label]])$log_lik
    data.frame(
      endpoint = label,
      marginal = as.numeric(lik["marginal"]),
      copula = as.numeric(lik["copula"]),
      joint = as.numeric(lik["joint"]),
      stringsAsFactors = FALSE
    )
  })
  out <- do.call(rbind, rows)
  out$delta_joint_vs_separate <- out$joint - out$joint[out$endpoint == "separate_final"]
  out$delta_joint_vs_joint <- out$joint - out$joint[out$endpoint == "joint_final"]
  if ("joint_from_separate_final" %in% out$endpoint) {
    out$delta_joint_vs_joint_from_separate <-
      out$joint - out$joint[out$endpoint == "joint_from_separate_final"]
  }
  out
}

selected_gradient_parameters <- function(par_vec) {
  candidates <- c(
    grep("^mu\\.", names(par_vec), value = TRUE)[1],
    grep("^sigma\\.", names(par_vec), value = TRUE)[1],
    grep("^theta\\.", names(par_vec), value = TRUE)
  )
  unique(stats::na.omit(candidates))
}

finite_difference_gradients <- function(fit, par_vec, label, h = 1e-5) {
  params <- selected_gradient_parameters(par_vec)
  if (length(params) == 0) {
    return(data.frame())
  }

  rows <- lapply(params, function(pn) {
    p_plus <- par_vec
    p_minus <- par_vec
    p_plus[pn] <- p_plus[pn] + h
    p_minus[pn] <- p_minus[pn] - h

    lik_plus <- evaluate_coefficients(fit, p_plus)$log_lik
    lik_minus <- evaluate_coefficients(fit, p_minus)$log_lik

    grad <- (lik_plus - lik_minus) / (2 * h)
    data.frame(
      endpoint = label,
      parameter = pn,
      grad_marginal = as.numeric(grad["marginal"]),
      grad_copula = as.numeric(grad["copula"]),
      grad_joint = as.numeric(grad["joint"]),
      stringsAsFactors = FALSE
    )
  })

  do.call(rbind, rows)
}

plot_likelihood_trace <- function(trace, file) {
  png(file, width = 1100, height = 760)
  oldpar <- par(no.readonly = TRUE)
  on.exit({
    par(oldpar)
    dev.off()
  }, add = TRUE)

  par(mfrow = c(3, 1), mar = c(4, 4, 2, 1))
  cols <- c(
    separate = "firebrick",
    joint = "steelblue",
    joint_from_separate = "darkgreen"
  )
  for (component in c("joint", "marginal", "copula")) {
    yy <- trace[[component]]
    plot(
      NA,
      xlim = range(trace$step, na.rm = TRUE),
      ylim = range(yy, na.rm = TRUE),
      xlab = "Stored coefficient step",
      ylab = paste(component, "log-likelihood"),
      main = paste("RS path:", component)
    )
    for (lab in unique(trace$fit)) {
      idx <- trace$fit == lab
      lines(trace$step[idx], trace[[component]][idx], col = cols[[lab]], lwd = 2)
    }
    used_cols <- cols[names(cols) %in% unique(trace$fit)]
    legend("bottomright", legend = names(used_cols), col = used_cols, lwd = 2, bty = "n")
  }
}

plot_endpoint_bars <- function(endpoints, file) {
  png(file, width = 1000, height = 650)
  oldpar <- par(no.readonly = TRUE)
  on.exit({
    par(oldpar)
    dev.off()
  }, add = TRUE)

  mat <- t(as.matrix(endpoints[, c("marginal", "copula", "joint")]))
  colnames(mat) <- endpoints$endpoint
  barplot(
    mat,
    beside = TRUE,
    las = 2,
    col = c("grey55", "darkseagreen3", "steelblue"),
    ylab = "Log-likelihood",
    main = "Endpoint cross-evaluation",
    legend.text = rownames(mat),
    args.legend = list(x = "bottomright", bty = "n")
  )
}

mark_progress("Simulating Gamma/Gaussian diagnostic dataset...")
dat <- simulate_gamma_gaussian(n_subject, n_time)
write.csv(dat, file.path(out_dir, "simulated_gamma_gaussian_data.csv"), row.names = FALSE)

mark_progress("Fitting RS separate optimisation: include_dlcopdpar = FALSE")
fit_sep <- fit_rs(dat, include_dlcopdpar = FALSE)

mark_progress("Fitting RS joint optimisation: include_dlcopdpar = TRUE")
fit_joint <- fit_rs(dat, include_dlcopdpar = TRUE)

mark_progress("Fitting RS joint optimisation from separate endpoint")
fit_joint_from_sep <- fit_rs(
  dat,
  include_dlcopdpar = TRUE,
  start_from = fit_sep$par
)

writeLines(attr(fit_sep, "captured_output"), file.path(out_dir, "fit_separate_output.log"))
writeLines(attr(fit_joint, "captured_output"), file.path(out_dir, "fit_joint_output.log"))
writeLines(attr(fit_joint_from_sep, "captured_output"), file.path(out_dir, "fit_joint_from_separate_output.log"))

mark_progress("Writing likelihood paths...")
trace_sep <- trace_from_fit_history(fit_sep, "separate")
trace_joint <- trace_from_fit_history(fit_joint, "joint")
trace_joint_from_sep <- trace_from_fit_history(fit_joint_from_sep, "joint_from_separate")
trace <- rbind(trace_sep, trace_joint, trace_joint_from_sep)
write.csv(trace, file.path(out_dir, "likelihood_trace.csv"), row.names = FALSE)

par_replay <- rbind(
  replay_path_from_par_history(fit_sep, "separate"),
  replay_path_from_par_history(fit_joint, "joint"),
  replay_path_from_par_history(fit_joint_from_sep, "joint_from_separate")
)
write.csv(par_replay, file.path(out_dir, "par_history_replayed_likelihood_trace.csv"), row.names = FALSE)

mark_progress("Comparing coefficients and endpoints...")
fits <- list(
  separate = fit_sep,
  joint = fit_joint,
  joint_from_separate = fit_joint_from_sep
)
coef_cmp <- coefficient_table(fits)
write.csv(coef_cmp, file.path(out_dir, "coefficient_comparison.csv"), row.names = FALSE)

theta_names <- grep("^theta\\.", names(fit_sep$par), value = TRUE)
margin_names <- setdiff(names(fit_sep$par), theta_names)

par_sep <- fit_sep$par
par_joint <- fit_joint$par
par_joint_from_sep <- fit_joint_from_sep$par
par_sep_margin_joint_theta <- par_sep
par_sep_margin_joint_theta[theta_names] <- par_joint[theta_names]
par_joint_margin_sep_theta <- par_joint
par_joint_margin_sep_theta[theta_names] <- par_sep[theta_names]
par_sep_margin_joint_from_sep_theta <- par_sep
par_sep_margin_joint_from_sep_theta[theta_names] <- par_joint_from_sep[theta_names]
par_joint_from_sep_margin_sep_theta <- par_joint_from_sep
par_joint_from_sep_margin_sep_theta[theta_names] <- par_sep[theta_names]

endpoints <- endpoint_rows(
  fit_ref = fit_sep,
  par_list = list(
    separate_final = par_sep,
    joint_final = par_joint,
    joint_from_separate_final = par_joint_from_sep,
    separate_margin_joint_theta = par_sep_margin_joint_theta,
    joint_margin_separate_theta = par_joint_margin_sep_theta,
    separate_margin_joint_from_separate_theta = par_sep_margin_joint_from_sep_theta,
    joint_from_separate_margin_separate_theta = par_joint_from_sep_margin_sep_theta
  )
)
write.csv(endpoints, file.path(out_dir, "endpoint_cross_evaluation.csv"), row.names = FALSE)

mark_progress("Computing finite-difference endpoint gradients...")
grads <- rbind(
  finite_difference_gradients(fit_sep, par_sep, "separate_final"),
  finite_difference_gradients(fit_sep, par_joint, "joint_final"),
  finite_difference_gradients(fit_sep, par_joint_from_sep, "joint_from_separate_final")
)
write.csv(grads, file.path(out_dir, "endpoint_finite_difference_gradients.csv"), row.names = FALSE)

mark_progress("Writing plots...")
plot_likelihood_trace(trace, file.path(out_dir, "likelihood_trace.png"))
plot_endpoint_bars(endpoints, file.path(out_dir, "endpoint_cross_evaluation.png"))

summary_txt <- c(
  "RS Gamma/Gaussian joint vs separate diagnostic",
  paste("n_subject:", n_subject),
  paste("n_time:", n_time),
  "",
  "Final log-likelihoods:",
  capture.output(print(endpoints, row.names = FALSE)),
  "",
  "Largest absolute coefficient differences:",
  capture.output(print(
    coef_cmp[order(abs(coef_cmp$diff_joint_minus_separate), decreasing = TRUE), ][seq_len(min(10, nrow(coef_cmp))), ],
    row.names = FALSE
  )),
  "",
  "Finite-difference gradients at endpoints:",
  capture.output(print(grads, row.names = FALSE))
)
writeLines(summary_txt, file.path(out_dir, "summary.txt"))

mark_progress("Done. Outputs written to: ", out_dir)
cat("\nEndpoint comparison:\n")
print(endpoints, row.names = FALSE)
