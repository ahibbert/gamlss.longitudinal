#!/usr/bin/env Rscript

#region Setup and Libraries
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
library("VineCopula")
library("moments")
library("ggplot2")
library("latex2exp")
library("ggpubr")
library("Matrix")
library("MASS")
#library(gamlss.longitudinal); 
library(gamlss2); library(gamlss)
set.seed(100)


set.seed(20260512)

out_dir <- file.path("results", "recovery_suite")
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)
dist_plot_dir <- file.path(out_dir, "distribution_plots")
dir.create(dist_plot_dir, recursive = TRUE, showWarnings = FALSE)

n_reps <- 1
checkpoint_every <- 25
smooth_k <- 10
verbose <- TRUE
plot_distribution_reps <- 1

n_reps_env <- Sys.getenv("N_REPS", unset = "")
if (nzchar(n_reps_env)) {
  n_reps <- as.integer(n_reps_env)
}

checkpoint_env <- Sys.getenv("CHECKPOINT_EVERY", unset = "")
if (nzchar(checkpoint_env)) {
  checkpoint_every <- as.integer(checkpoint_env)
}

scenarios <- data.frame(
  n=c(100), #n = c(100, 100, 500, 500),
  d=c(4), #d = c(4, 10, 4, 10),
  stringsAsFactors = FALSE
)

fixed_terms <- c("intercept", "x1", "x2", "t")
params_all <- c("mu", "sigma", "nu", "tau", "theta", "zeta")
params_smooth <- c("mu", "sigma", "theta")
params_margin <- c("mu", "sigma", "nu", "tau")

margin_truth <- gamlss.dist::BCPE(mu.link = "log")
copula_truth <- get_copula_dist("t")$copula_link

# Keep natural-scale intercepts explicit, then map to eta-scale via link functions.
true_intercepts_natural <- list(
  mu = 2.20,
  sigma = margin_truth$sigma.linkinv(1),
  nu = 1,
  tau = margin_truth$tau.linkinv(1),
  theta = copula_truth$theta.linkinv(0.5),
  zeta = copula_truth$zeta.linkinv(2.30)
)

true_beta <- list(
  mu = c(intercept = margin_truth$mu.linkfun(true_intercepts_natural$mu), x1 = 0.12, x2 = -0.10, t = 0.18),
  sigma = c(intercept = margin_truth$sigma.linkfun(true_intercepts_natural$sigma), x1 = 0.18, x2 = 0.12, t = 0.15),
  nu = c(intercept = 0.00, x1 = 0.10, x2 = -0.08, t = 0.12),
  tau = c(intercept = margin_truth$tau.linkfun(true_intercepts_natural$tau), x1 = -0.08, x2 = 0.10, t = 0.11),
  theta = c(intercept = copula_truth$theta.linkfun(true_intercepts_natural$theta), x1 = 0.25, x2 = -0.10, t = 0.20),
  zeta = c(intercept = copula_truth$zeta.linkfun(true_intercepts_natural$zeta), x1 = 0.12, x2 = -0.10, t = 0.10)
)

f_mu <- function(s1) 0.30 * sin(2 * pi * s1)
f_sigma <- function(s1) 0.30 * cos(2 * pi * s1) - 0.10 * (s1 - 0.5)^2
f_theta <- function(s1) 0.25 * sin(pi * s1 + 0.4) + 0.12 * (s1 - 0.5)

smooth_truth <- list(mu = f_mu, sigma = f_sigma, theta = f_theta)

calc_smooth_mean <- function(data_used, parameter) {
  data_sub <- data_used
  if (parameter == "theta") {
    data_sub <- data_sub[data_sub$time < max(data_sub$time), , drop = FALSE]
  }
  mean(smooth_truth[[parameter]](data_sub$s1))
}

linpred <- function(beta, x1, x2, t) {
  beta[["intercept"]] + beta[["x1"]] * x1 + beta[["x2"]] * x2 + beta[["t"]] * t
}

clip_u <- function(u) {
  pmax(pmin(u, 1 - 1e-8), 1e-8)
}

plot_distribution_snapshot <- function(dat, scenario_label, rep_index, offdiag_scale = "pseudo") {
  dat_plot <- data.frame(
    subject = as.factor(dat$id),
    time = dat$time,
    response = dat$response,
    stringsAsFactors = FALSE
  )

  p <- tryCatch({
    plotDist(dat_plot, margin_truth, offdiag_scale = offdiag_scale)
  }, error = function(e) {
    warning(
      sprintf(
        "Distribution plot failed for %s rep %d: %s",
        scenario_label,
        rep_index,
        conditionMessage(e)
      )
    )
    NULL
  })

  if (!is.null(p)) {
    fname <- sprintf("dist_%s_rep%03d.png", scenario_label, rep_index)
    ggplot2::ggsave(
      filename = file.path(dist_plot_dir, fname),
      plot = p,
      width = 3.2 * length(unique(dat_plot$time)),
      height = 3.2 * length(unique(dat_plot$time)),
      dpi = 160
    )
  }
}

simulate_dataset <- function(n, d, seed) {
  set.seed(seed)

  ids <- seq_len(n)
  subj <- data.frame(
    id = ids,
    x1 = rnorm(n),
    x2 = rbinom(n, 1, 0.5),
    s1 = runif(n, 0, 1),
    s2 = runif(n, 0, 1),
    stringsAsFactors = FALSE
  )

  long <- merge(
    subj,
    data.frame(time_index = seq_len(d)),
    by = NULL,
    all = TRUE
  )
  long <- long[order(long$id, long$time_index), ]
  long$t <- if (d > 1) (long$time_index - 1) / (d - 1) else 0
  long$time <- long$time_index

  eta_mu <- linpred(true_beta$mu, long$x1, long$x2, long$t) + f_mu(long$s1)
  eta_sigma <- linpred(true_beta$sigma, long$x1, long$x2, long$t) + f_sigma(long$s1)
  eta_nu <- linpred(true_beta$nu, long$x1, long$x2, long$t)
  eta_tau <- linpred(true_beta$tau, long$x1, long$x2, long$t)
  eta_theta <- linpred(true_beta$theta, long$x1, long$x2, long$t) + f_theta(long$s1)
  eta_zeta <- linpred(true_beta$zeta, long$x1, long$x2, long$t)

  margin <- gamlss.dist::BCPE(mu.link = "log")
  cop_spec <- get_copula_dist("t")

  mu <- margin$mu.linkinv(eta_mu)
  sigma <- margin$sigma.linkinv(eta_sigma)
  nu <- margin$nu.linkinv(eta_nu)
  tau <- margin$tau.linkinv(eta_tau)

  theta <- cop_spec$copula_link$theta.linkinv(eta_theta)
  zeta <- cop_spec$copula_link$zeta.linkinv(eta_zeta)

  U <- matrix(NA_real_, nrow = n, ncol = d)
  fam_t <- as.numeric(VineCopula::BiCopName("t"))

  for (i in seq_len(n)) {
    u_prev <- runif(1)
    U[i, 1] <- u_prev

    if (d > 1) {
      for (j in 2:d) {
        row_prev <- which(long$id == i & long$time_index == (j - 1))
        th <- theta[row_prev]
        ze <- zeta[row_prev]

        u_new <- tryCatch({
          VineCopula::BiCopCondSim(
            N = 1,
            cond.val = u_prev,
            cond.var = 1,
            family = fam_t,
            par = th,
            par2 = ze
          )
        }, error = function(e) runif(1))

        u_prev <- clip_u(as.numeric(u_new))
        U[i, j] <- u_prev
      }
    }
  }

  long$u <- as.vector(t(U))
  long$response <- gamlss.dist::qBCPE(long$u, mu = mu, sigma = sigma, nu = nu, tau = tau)

  long[, c("id", "time", "response", "s1", "s2", "x1", "x2", "t")]
}

extract_one_term <- function(par_vec, parameter, term_name) {
  nm <- names(par_vec)
  if (is.null(nm) || length(nm) == 0) {
    return(NA_real_)
  }

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
    x2 = "^x2($|[^[:alnum:]_])|^x2[[:alnum:]_]+$",
    t = "^t$",
    `x1:t` = "^(x1:t|t:x1)$",
    `x2:t` = "^(x2:t|t:x2|x2[[:alnum:]_]+:t|t:x2[[:alnum:]_]+)$"
  )

  hit <- grep(target_pattern, sub_names)
  if (length(hit) == 0) {
    return(NA_real_)
  }

  values[hit[1]]
}

extract_fixed_estimates <- function(fit_obj, data_used) {
  par_vec <- fit_obj$par
  out <- list()
  k <- 1

  smooth_means <- lapply(params_smooth, function(p) calc_smooth_mean(data_used, p))
  names(smooth_means) <- params_smooth

  for (p in params_all) {
    for (tm in fixed_terms) {
      true_val <- true_beta[[p]][[tm]]
      if (tm == "intercept" && p %in% params_smooth) {
        # Fitted smooths are centered, so intercept absorbs average smooth level.
        true_val <- true_val + smooth_means[[p]]
      }

      out[[k]] <- data.frame(
        parameter = p,
        term = tm,
        estimate = extract_one_term(par_vec, p, tm),
        true_value = true_val,
        stringsAsFactors = FALSE
      )
      k <- k + 1
    }
  }

  do.call(rbind, out)
}

smooth_curve_from_fit <- function(fit_obj, data_used, parameter, s1_grid) {
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
    if (parameter == "theta") {
      data_sub <- data_sub[data_sub$time < max(data_sub$time), , drop = FALSE]
    }
    x_obs <- data_sub$s1
  }

  ord <- order(x_obs)
  x_obs <- x_obs[ord]
  smooth_obs <- smooth_obs[ord]

  # Collapse repeated x values so interpolation is stable.
  smooth_by_x <- tapply(smooth_obs, x_obs, mean)
  x_unique <- as.numeric(names(smooth_by_x))
  y_unique <- as.numeric(smooth_by_x)

  if (length(x_unique) < 2) {
    return(rep(NA_real_, length(s1_grid)))
  }

  smooth_grid_hat <- stats::approx(x = x_unique, y = y_unique, xout = s1_grid, rule = 2)$y

  # Compare centered smooth components to match identifiability constraints.
  smooth_grid_hat - mean(smooth_grid_hat, na.rm = TRUE)
}

fit_model <- function(dat) {
  fit <- NULL
  #invisible(capture.output({
    fit <- gamlss.longitudinal(
      dataset = dat,
      margin_dist = gamlss.dist::BCPE(mu.link = "log"),
      copula_dist = "t",
      time_var = "time",
      subject_var = "id",
      mu.formula = sprintf("response ~ x1 + x2 + t + s(s1, k = %d)", smooth_k),
      sigma.formula = sprintf("~ x1 + x2 + t + s(s1, k = %d)", smooth_k),
      nu.formula = "~ x1 + x2 + t",
      tau.formula = "~ x1 + x2 + t",
      theta.formula = sprintf("~ x1 + x2 + t + s(s1, k = %d)", smooth_k),
      zeta.formula = "~ x1 + x2 + t",
      verbose = if (verbose) 1 else 0,
      compute_vcov = FALSE,
      include_dlcopdpar = TRUE,
      method = "RS"
    )
  #}, type = "output"))
  fit
}

fit_model_gamlss <- function(dat) {
  gamlss::gamlss(
    formula = stats::as.formula(sprintf("response ~ x1 + x2 + t + pb(s1, df = %d)", smooth_k)),
    sigma.formula = stats::as.formula(sprintf("~ x1 + x2 + t + pb(s1, df = %d)", smooth_k)),
    nu.formula = ~ x1 + x2 + t,
    tau.formula = ~ x1 + x2 + t,
    family = gamlss.dist::BCPE(mu.link = "log"),
    data = dat,
    trace = isTRUE(verbose),
    control = gamlss::gamlss.control(n.cyc = 200)
  )
}

extract_fixed_estimates_gamlss <- function(fit_obj, data_used) {
  out <- list()
  k <- 1

  smooth_means <- lapply(params_smooth, function(p) calc_smooth_mean(data_used, p))
  names(smooth_means) <- params_smooth

  coef_by_param <- list(
    mu = stats::coef(fit_obj, what = "mu"),
    sigma = stats::coef(fit_obj, what = "sigma"),
    nu = stats::coef(fit_obj, what = "nu"),
    tau = stats::coef(fit_obj, what = "tau")
  )

  extract_one_from_named <- function(coef_vec, term_name) {
    nm <- names(coef_vec)
    if (is.null(nm) || length(nm) == 0) {
      return(NA_real_)
    }

    target_pattern <- switch(
      term_name,
      intercept = "^(intercept|\\(Intercept\\))$",
      x1 = "^x1$",
      x2 = "^x2($|[^[:alnum:]_])|^x2[[:alnum:]_]+$",
      t = "^t$"
    )
    hit <- grep(target_pattern, nm)
    if (length(hit) == 0) {
      return(NA_real_)
    }
    as.numeric(coef_vec[hit[1]])
  }

  for (p in params_margin) {
    for (tm in fixed_terms) {
      true_val <- true_beta[[p]][[tm]]
      if (tm == "intercept" && p %in% params_smooth) {
        true_val <- true_val + smooth_means[[p]]
      }

      out[[k]] <- data.frame(
        parameter = p,
        term = tm,
        estimate = extract_one_from_named(coef_by_param[[p]], tm),
        true_value = true_val,
        stringsAsFactors = FALSE
      )
      k <- k + 1
    }
  }

  do.call(rbind, out)
}

safe_quantile <- function(x, p) {
  x <- x[is.finite(x)]
  if (length(x) == 0) {
    return(NA_real_)
  }
  as.numeric(stats::quantile(x, probs = p, names = FALSE, type = 7))
}

all_fixed <- list()
all_smooth <- list()
fit_log <- list()

bind_non_null <- function(x) {
  keep <- x[!vapply(x, is.null, logical(1))]
  if (length(keep) == 0) {
    return(NULL)
  }
  do.call(rbind, keep)
}

for (s in seq_len(nrow(scenarios))) {
  n_val <- scenarios$n[s]
  d_val <- scenarios$d[s]

  cat(sprintf("\nScenario %d/%d: n=%d, d=%d\n", s, nrow(scenarios), n_val, d_val))

  fixed_rows <- list()
  smooth_rows <- list()
  run_rows <- list()

  s1_grid <- seq(0, 1, length.out = 101)

  for (r in seq_len(n_reps)) {
    seed <- 100000 + 10000 * s + r
    t0_gl <- Sys.time()

    dat <- simulate_dataset(n = n_val, d = d_val, seed = seed)
    scenario_label <- sprintf("n%d_d%d", n_val, d_val)

    if (r <= plot_distribution_reps) {
      plot_distribution_snapshot(dat, scenario_label, r, offdiag_scale = "pseudo")
    }

    fit_obj <- NULL
    err_gl <- NA_character_
    fit_gamlss <- NULL
    err_gamlss <- NA_character_

    fit_obj <- tryCatch({
      fit_model(dat)
    }, error = function(e) {
      err_gl <<- conditionMessage(e)
      NULL
    })

    elapsed_gl <- as.numeric(difftime(Sys.time(), t0_gl, units = "secs"))

    t0_g <- Sys.time()
    fit_gamlss <- tryCatch({
      fit_model_gamlss(dat)
    }, error = function(e) {
      err_gamlss <<- conditionMessage(e)
      NULL
    })

    elapsed_g <- as.numeric(difftime(Sys.time(), t0_g, units = "secs"))

    run_rows[[length(run_rows) + 1]] <- data.frame(
      scenario = scenario_label,
      model = "gamlss.longitudinal",
      rep = r,
      success = !is.null(fit_obj),
      elapsed_sec = elapsed_gl,
      error = if (is.null(fit_obj)) err_gl else NA_character_,
      stringsAsFactors = FALSE
    )

    run_rows[[length(run_rows) + 1]] <- data.frame(
      scenario = scenario_label,
      model = "gamlss",
      rep = r,
      success = !is.null(fit_gamlss),
      elapsed_sec = elapsed_g,
      error = if (is.null(fit_gamlss)) err_gamlss else NA_character_,
      stringsAsFactors = FALSE
    )

    if (!is.null(fit_obj)) {
      ftab <- extract_fixed_estimates(fit_obj, dat)
      ftab$model <- "gamlss.longitudinal"
      ftab$scenario <- scenario_label
      ftab$n <- n_val
      ftab$d <- d_val
      ftab$rep <- r
      fixed_rows[[length(fixed_rows) + 1]] <- ftab

      for (p in params_smooth) {
        curve_hat <- smooth_curve_from_fit(fit_obj, dat, p, s1_grid)
        smooth_mean <- calc_smooth_mean(dat, p)
        curve_true <- smooth_truth[[p]](s1_grid) - smooth_mean

        smooth_rows[[length(smooth_rows) + 1]] <- data.frame(
          scenario = scenario_label,
          n = n_val,
          d = d_val,
          rep = r,
          parameter = p,
          s1 = s1_grid,
          smooth_hat = curve_hat,
          smooth_true = curve_true,
          stringsAsFactors = FALSE
        )
      }
    }

    if (!is.null(fit_gamlss)) {
      ftab_g <- extract_fixed_estimates_gamlss(fit_gamlss, dat)
      ftab_g$model <- "gamlss"
      ftab_g$scenario <- scenario_label
      ftab_g$n <- n_val
      ftab_g$d <- d_val
      ftab_g$rep <- r
      fixed_rows[[length(fixed_rows) + 1]] <- ftab_g
    }

    if (r %% checkpoint_every == 0) {
      cat(sprintf("  Rep %d/%d complete\n", r, n_reps))

      saveRDS(
        list(
          fixed = if (length(fixed_rows) > 0) do.call(rbind, fixed_rows) else NULL,
          smooth = if (length(smooth_rows) > 0) do.call(rbind, smooth_rows) else NULL,
          runs = do.call(rbind, run_rows)
        ),
        file.path(out_dir, sprintf("checkpoint_n%d_d%d_rep%d.rds", n_val, d_val, r))
      )
    }
  }

  fixed_df <- if (length(fixed_rows) > 0) do.call(rbind, fixed_rows) else NULL
  smooth_df <- if (length(smooth_rows) > 0) do.call(rbind, smooth_rows) else NULL
  run_df <- do.call(rbind, run_rows)

  saveRDS(
    list(fixed = fixed_df, smooth = smooth_df, runs = run_df),
    file.path(out_dir, sprintf("scenario_n%d_d%d_results.rds", n_val, d_val))
  )

  all_fixed[[s]] <- fixed_df
  all_smooth[[s]] <- smooth_df
  fit_log[[s]] <- run_df
}

fixed_all <- bind_non_null(all_fixed)
smooth_all <- bind_non_null(all_smooth)
runs_all <- bind_non_null(fit_log)

saveRDS(list(fixed = fixed_all, smooth = smooth_all, runs = runs_all), file.path(out_dir, "all_results.rds"))
write.csv(runs_all, file.path(out_dir, "fit_run_log.csv"), row.names = FALSE)

if (is.null(fixed_all) || nrow(fixed_all) == 0) {
  warning("No successful fits were available for fixed-effect summaries. Check fit_run_log.csv for errors.")
  if (!is.null(runs_all) && nrow(runs_all) > 0) {
    failed <- runs_all[!runs_all$success, c("scenario", "rep", "error"), drop = FALSE]
    failed <- failed[!is.na(failed$error) & nzchar(failed$error), , drop = FALSE]
    if (nrow(failed) > 0) {
      cat("\nFirst fit errors:\n")
      print(utils::head(failed, 10))
    }
  }
  quit(save = "no", status = 1)
}

fixed_summary <- do.call(rbind, lapply(split(fixed_all, list(fixed_all$scenario, fixed_all$parameter, fixed_all$term), drop = TRUE), function(df) {
  err <- df$estimate - df$true_value
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
    sd_estimate = stats::sd(df$estimate, na.rm = TRUE),
    q05 = safe_quantile(df$estimate, 0.05),
    q95 = safe_quantile(df$estimate, 0.95),
    n_successful_fits = sum(is.finite(df$estimate)),
    stringsAsFactors = FALSE
  )
}))

write.csv(fixed_summary, file.path(out_dir, "fixed_effects_bias_rmse_table.csv"), row.names = FALSE)

if (is.null(smooth_all) || nrow(smooth_all) == 0) {
  warning("No successful fits were available for smooth summaries. Fixed-effect summary was still saved.")
  if (!is.null(runs_all) && nrow(runs_all) > 0) {
    failed <- runs_all[!runs_all$success, c("scenario", "rep", "error"), drop = FALSE]
    failed <- failed[!is.na(failed$error) & nzchar(failed$error), , drop = FALSE]
    if (nrow(failed) > 0) {
      cat("\nFirst fit errors:\n")
      print(utils::head(failed, 10))
    }
  }
  quit(save = "no", status = 1)
}

smooth_summary <- do.call(rbind, lapply(split(smooth_all, list(smooth_all$scenario, smooth_all$parameter, smooth_all$s1), drop = TRUE), function(df) {
  err <- df$smooth_hat - df$smooth_true
  data.frame(
    scenario = df$scenario[1],
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
}))

write.csv(smooth_summary, file.path(out_dir, "smooth_pointwise_summary.csv"), row.names = FALSE)

smooth_integrated <- do.call(rbind, lapply(split(smooth_all, list(smooth_all$scenario, smooth_all$parameter, smooth_all$rep), drop = TRUE), function(df) {
  err <- df$smooth_hat - df$smooth_true
  data.frame(
    scenario = df$scenario[1],
    n = df$n[1],
    d = df$d[1],
    parameter = df$parameter[1],
    rep = df$rep[1],
    irmse = sqrt(mean(err^2, na.rm = TRUE)),
    max_abs_error = max(abs(err), na.rm = TRUE),
    stringsAsFactors = FALSE
  )
}))

smooth_integrated_summary <- do.call(rbind, lapply(split(smooth_integrated, list(smooth_integrated$scenario, smooth_integrated$parameter), drop = TRUE), function(df) {
  data.frame(
    scenario = df$scenario[1],
    n = df$n[1],
    d = df$d[1],
    parameter = df$parameter[1],
    mean_irmse = mean(df$irmse, na.rm = TRUE),
    sd_irmse = stats::sd(df$irmse, na.rm = TRUE),
    mean_max_abs_error = mean(df$max_abs_error, na.rm = TRUE),
    stringsAsFactors = FALSE
  )
}))

write.csv(smooth_integrated_summary, file.path(out_dir, "smooth_integrated_metrics.csv"), row.names = FALSE)

plot_fixed <- ggplot(fixed_summary, aes(x = term, y = mean_estimate)) +
  geom_errorbar(aes(ymin = q05, ymax = q95), width = 0.2, color = "gray50") +
  geom_point(size = 1.7, color = "black") +
  geom_point(aes(y = true_value), shape = 21, size = 3.2, stroke = 1.1, fill = "white", color = "#D1495B") +
  geom_point(aes(y = true_value), size = 1.4, color = "#D1495B") +
  facet_wrap(vars(model, parameter, scenario), ncol = 3, scales = "free_y") +
  labs(
    title = "Fixed-effect recovery by scenario and model",
    subtitle = "Black points: mean estimate across replications, bars: 5%-95% quantiles, red bullseyes: true values",
    x = "Term",
    y = "Estimate"
  ) +
  theme_bw(base_size = 11) +
  theme(axis.text.x = element_text(angle = 30, hjust = 1))

ggsave(file.path(out_dir, "fixed_effect_recovery_plot.png"), plot_fixed, width = 15, height = 12, dpi = 160)

plot_smooth <- ggplot(smooth_summary, aes(x = s1)) +
  geom_ribbon(aes(ymin = smooth_q05, ymax = smooth_q95), fill = "gray80", alpha = 0.9) +
  geom_line(aes(y = smooth_median), color = "gray25", linewidth = 0.7, linetype = "dashed") +
  geom_line(aes(y = smooth_true), color = "black", linewidth = 0.9) +
  facet_wrap(vars(parameter, scenario), ncol = 3, scales = "free_y") +
  labs(
    title = "Smooth recovery for mu, sigma, theta",
    subtitle = "Black: true smooth, dashed: median fitted smooth, shaded: 5%-95% envelope",
    x = "s1",
    y = "Smooth contribution on eta scale"
  ) +
  theme_bw(base_size = 11)

ggsave(file.path(out_dir, "smooth_recovery_plot.png"), plot_smooth, width = 15, height = 10, dpi = 160)

cat("\nSimulation suite complete. Outputs written to:", out_dir, "\n")
