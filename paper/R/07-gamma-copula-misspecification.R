jss_misspec_copulas <- function() {
  c("N", "C", "F", "G", "J", "t")
}

jss_misspec_tau_levels <- function() {
  data.frame(
    tau_label = c("moderate", "high"),
    target_tau = c(0.25, 0.55),
    stringsAsFactors = FALSE
  )
}

jss_misspec_stage <- function(profile) {
  requested <- Sys.getenv("GAMLSS_LONGITUDINAL_JSS_MISSPEC_STAGE", unset = "")
  if (nzchar(requested)) {
    return(match.arg(requested, c("smoke", "pilot", "full")))
  }
  if (identical(profile, "expanded")) "pilot" else "smoke"
}

jss_misspec_config <- function(settings, stage = jss_misspec_stage(settings$profile)) {
  stage <- match.arg(stage, c("smoke", "pilot", "full"))
  reps_override <- Sys.getenv("GAMLSS_LONGITUDINAL_JSS_MISSPEC_REPS", unset = "")
  reps <- switch(stage, smoke = 1L, pilot = 10L, full = 100L)
  if (nzchar(reps_override)) {
    reps <- as.integer(reps_override)
    if (is.na(reps) || reps < 1L) {
      stop("GAMLSS_LONGITUDINAL_JSS_MISSPEC_REPS must be a positive integer.", call. = FALSE)
    }
  }

  list(
    stage = stage,
    copulas = jss_misspec_copulas(),
    tau_levels = if (identical(stage, "smoke")) {
      jss_misspec_tau_levels()[1L, , drop = FALSE]
    } else {
      jss_misspec_tau_levels()
    },
    sample_sizes = if (identical(stage, "smoke")) 20L else c(50L, 150L, 500L),
    times = 1:4,
    reps = reps,
    margin_mu = 2,
    margin_sigma = 0.5,
    t_zeta = 5,
    max_outer_iter = if (identical(stage, "smoke")) 3L else 100L,
    max_inner_iter = if (identical(stage, "smoke")) 3L else 100L,
    max_elapsed_sec = if (identical(stage, "smoke")) 20 else Inf,
    compute_vcov = !identical(stage, "smoke"),
    seed = settings$seed + 700000L
  )
}

jss_misspec_paths <- function(settings) {
  module_id <- "07-gamma-copula-misspecification"
  list(
    results = file.path(settings$data_dir, paste0(module_id, "-results.csv")),
    grid = file.path(settings$data_dir, paste0(module_id, "-grid.csv")),
    summary = file.path(settings$tables_dir, paste0(module_id, "-summary.csv")),
    selection = file.path(settings$tables_dir, paste0(module_id, "-selection.csv")),
    review = file.path(settings$tables_dir, paste0(module_id, "-review-gate.csv")),
    heatmap = file.path(settings$figures_dir, paste0(module_id, "-delta-heatmap.png")),
    margin_rmse_heatmap = file.path(settings$figures_dir, paste0(module_id, "-margin-rmse-heatmap.png")),
    tau_error_heatmap = file.path(settings$figures_dir, paste0(module_id, "-tau-error-heatmap.png")),
    paper_summary_heatmap = file.path(settings$figures_dir, paste0(module_id, "-paper-summary-heatmap.png")),
    convergence = file.path(settings$figures_dir, paste0(module_id, "-convergence.png")),
    checkpoints = file.path(settings$data_dir, paste0(module_id, "-checkpoints")),
    smoke_checkpoints = file.path(settings$data_dir, paste0(module_id, "-smoke-checkpoints"))
  )
}

jss_misspec_grid <- function(config) {
  base <- expand.grid(
    generating_copula = config$copulas,
    fitted_copula = config$copulas,
    tau_label = config$tau_levels$tau_label,
    n_subject = config$sample_sizes,
    rep = seq_len(config$reps),
    stringsAsFactors = FALSE
  )
  tau_map <- stats::setNames(config$tau_levels$target_tau, config$tau_levels$tau_label)
  base$target_tau <- unname(tau_map[base$tau_label])
  base$n_time <- length(config$times)
  base$seed <- mapply(
    jss_misspec_seed,
    base$generating_copula,
    base$fitted_copula,
    base$tau_label,
    base$n_subject,
    base$rep,
    MoreArgs = list(base_seed = config$seed),
    SIMPLIFY = TRUE,
    USE.NAMES = FALSE
  )
  base$fit_id <- mapply(
    jss_misspec_fit_id,
    base$generating_copula,
    base$fitted_copula,
    base$tau_label,
    base$n_subject,
    base$rep,
    SIMPLIFY = TRUE,
    USE.NAMES = FALSE
  )
  base
}

jss_misspec_seed <- function(generating_copula, fitted_copula, tau_label, n_subject, rep, base_seed) {
  copulas <- jss_misspec_copulas()
  gen_i <- match(generating_copula, copulas)
  fit_i <- match(fitted_copula, copulas)
  tau_i <- match(tau_label, jss_misspec_tau_levels()$tau_label)
  as.integer(base_seed + rep * 10000L + n_subject * 10L + gen_i * 100L + fit_i * 7L + tau_i)
}

jss_misspec_dataset_seed <- function(generating_copula, tau_label, n_subject, rep, base_seed) {
  jss_misspec_seed(generating_copula, generating_copula, tau_label, n_subject, rep, base_seed)
}

jss_misspec_fit_id <- function(generating_copula, fitted_copula, tau_label, n_subject, rep) {
  paste(
    paste0("gen-", generating_copula),
    paste0("fit-", fitted_copula),
    paste0("tau-", tau_label),
    paste0("n-", n_subject),
    sprintf("rep-%03d", as.integer(rep)),
    sep = "__"
  )
}

jss_misspec_checkpoint_dir <- function(paths, stage) {
  if (identical(stage, "smoke")) paths$smoke_checkpoints else paths$checkpoints
}

jss_misspec_checkpoint_path <- function(checkpoint_dir, fit_id) {
  file.path(checkpoint_dir, paste0(gsub("[^A-Za-z0-9_-]+", "_", fit_id), ".csv"))
}

jss_misspec_pending_grid <- function(grid, checkpoint_dir) {
  paths <- jss_misspec_checkpoint_path(checkpoint_dir, grid$fit_id)
  grid$checkpoint_path <- paths
  grid$checkpoint_exists <- file.exists(paths)
  grid[!grid$checkpoint_exists, , drop = FALSE]
}

jss_misspec_row_bind <- function(rows) {
  rows <- rows[!vapply(rows, is.null, logical(1))]
  if (!length(rows)) {
    return(data.frame())
  }
  all_names <- unique(unlist(lapply(rows, names), use.names = FALSE))
  rows <- lapply(rows, function(row) {
    missing <- setdiff(all_names, names(row))
    for (name in missing) row[[name]] <- NA
    row[all_names]
  })
  out <- do.call(rbind, rows)
  rownames(out) <- NULL
  out
}

jss_misspec_simulate <- function(generating_copula, target_tau, n_subject, config, seed) {
  gamlss.longitudinal::simulate_longitudinal_dataset(
    n = n_subject,
    times = config$times,
    margin_dist = gamlss.dist::GA(mu.link = "log", sigma.link = "log"),
    copula_dist = generating_copula,
    margin_params = list(mu = config$margin_mu, sigma = config$margin_sigma),
    copula_params = if (identical(generating_copula, "t")) {
      list(tau = target_tau, zeta = config$t_zeta)
    } else {
      list(tau = target_tau)
    },
    seed = seed,
    include_truth = TRUE
  )
}

jss_misspec_fit_one <- function(dat, row, config) {
  warnings <- character(0)
  start <- proc.time()[["elapsed"]]
  fit <- withCallingHandlers(
    tryCatch(
      gamlss.longitudinal::gamlss_longitudinal(
        dataset = dat,
        margin_dist = gamlss.dist::GA(mu.link = "log", sigma.link = "log"),
        copula_dist = row$fitted_copula,
        time_var = "time",
        subject_var = "subject",
        mu.formula = response ~ 1,
        sigma.formula = ~1,
        nu.formula = ~1,
        tau.formula = ~1,
        theta.formula = ~1,
        zeta.formula = ~1,
        method = "RS",
        max_outer_iter = config$max_outer_iter,
        max_inner_iter = config$max_inner_iter,
        max_elapsed_sec = config$max_elapsed_sec,
        compute_vcov = config$compute_vcov,
        verbose = 0
      ),
      error = function(e) e
    ),
    warning = function(w) {
      warnings <<- c(warnings, conditionMessage(w))
      invokeRestart("muffleWarning")
    }
  )
  elapsed <- proc.time()[["elapsed"]] - start

  jss_misspec_result_row(
    fit = fit,
    dat = dat,
    row = row,
    config = config,
    elapsed = elapsed,
    warnings = warnings
  )
}

jss_misspec_result_row <- function(fit, dat, row, config, elapsed, warnings) {
  is_fit <- inherits(fit, "gamlss.longitudinal")
  loglik <- if (is_fit) fit$calc_lik_out_end$log_lik else c(marginal = NA_real_, copula = NA_real_, joint = NA_real_)
  estimates <- if (is_fit) gamlss.longitudinal:::.coverage_natural_estimates(fit) else stats::setNames(numeric(0), character(0))
  fitted_tau <- if (is_fit && "theta" %in% names(estimates)) {
    zeta <- if ("zeta" %in% names(estimates)) estimates[["zeta"]] else 0
    suppressWarnings(gamlss.longitudinal:::.copula_par_to_tau(row$fitted_copula, estimates[["theta"]], zeta))
  } else {
    NA_real_
  }
  k <- if (is_fit) length(fit$par) else NA_integer_
  n_obs <- sum(is.finite(dat$response))
  joint_loglik <- as.numeric(loglik["joint"])
  aic_joint <- if (is.finite(joint_loglik) && is.finite(k)) -2 * joint_loglik + 2 * k else NA_real_
  bic_joint <- if (is.finite(joint_loglik) && is.finite(k)) -2 * joint_loglik + log(n_obs) * k else NA_real_
  truth_metrics <- if (is_fit) {
    gamlss.longitudinal:::.coverage_benchmark_gamlss_metrics(dat, fit, "GA")
  } else {
    gamlss.longitudinal:::.coverage_benchmark_gamlss_metrics(dat, NULL, "GA")
  }
  ci <- jss_misspec_margin_ci(fit)

  success <- is_fit && all(is.finite(as.numeric(loglik)))
  data.frame(
    fit_id = row$fit_id,
    stage = config$stage,
    generating_copula = row$generating_copula,
    fitted_copula = row$fitted_copula,
    correctly_specified = identical(row$generating_copula, row$fitted_copula),
    tau_label = row$tau_label,
    target_tau = row$target_tau,
    n_subject = row$n_subject,
    n_time = row$n_time,
    rep = row$rep,
    seed = row$seed,
    success = success,
    converged = if (is_fit && !is.null(fit$convergence$converged)) isTRUE(fit$convergence$converged) else FALSE,
    failure_type = if (success) "none" else if (inherits(fit, "error")) "error" else "nonfinite_or_no_fit",
    error = if (inherits(fit, "error")) conditionMessage(fit) else NA_character_,
    warnings = if (length(warnings)) paste(unique(warnings), collapse = " | ") else NA_character_,
    marginal_loglik = as.numeric(loglik["marginal"]),
    copula_loglik = as.numeric(loglik["copula"]),
    joint_loglik = joint_loglik,
    aic_joint = aic_joint,
    bic_joint = bic_joint,
    elapsed_sec = elapsed,
    n_parameters = k,
    fitted_mu = if ("mu" %in% names(estimates)) unname(estimates[["mu"]]) else NA_real_,
    fitted_sigma = if ("sigma" %in% names(estimates)) unname(estimates[["sigma"]]) else NA_real_,
    mu_bias = if ("mu" %in% names(estimates)) unname(estimates[["mu"]] - config$margin_mu) else NA_real_,
    sigma_bias = if ("sigma" %in% names(estimates)) unname(estimates[["sigma"]] - config$margin_sigma) else NA_real_,
    margin_param_rmse = jss_misspec_margin_rmse(estimates, config),
    mu_covered_95 = ci$mu_covered_95,
    sigma_covered_95 = ci$sigma_covered_95,
    fitted_copula_tau = as.numeric(fitted_tau)[1],
    tau_abs_error = abs(as.numeric(fitted_tau)[1] - row$target_tau),
    benchmark_mae = unname(truth_metrics[["benchmark_mae"]]),
    benchmark_rmse = unname(truth_metrics[["benchmark_rmse"]]),
    benchmark_mean_bias = unname(truth_metrics[["benchmark_mean_bias"]]),
    benchmark_mean_mae = unname(truth_metrics[["benchmark_mean_mae"]]),
    benchmark_mean_rmse = unname(truth_metrics[["benchmark_mean_rmse"]]),
    benchmark_q90_mae = unname(truth_metrics[["benchmark_q90_mae"]]),
    benchmark_neg_log_score = unname(truth_metrics[["benchmark_neg_log_score"]]),
    benchmark_upper_tail_error_90 = unname(truth_metrics[["benchmark_upper_tail_error_90"]]),
    benchmark_interval_coverage_95 = unname(truth_metrics[["benchmark_interval_coverage_95"]]),
    benchmark_interval_width_95 = unname(truth_metrics[["benchmark_interval_width_95"]]),
    benchmark_pit_ks_p_value = unname(truth_metrics[["benchmark_pit_ks_p_value"]]),
    benchmark_pit_mean_abs_error = unname(truth_metrics[["benchmark_pit_mean_abs_error"]]),
    benchmark_tail_error_lower_05 = unname(truth_metrics[["benchmark_tail_error_lower_05"]]),
    benchmark_tail_error_upper_05 = unname(truth_metrics[["benchmark_tail_error_upper_05"]]),
    stringsAsFactors = FALSE
  )
}

jss_misspec_margin_rmse <- function(estimates, config) {
  err <- c(
    if ("mu" %in% names(estimates)) estimates[["mu"]] - config$margin_mu else NA_real_,
    if ("sigma" %in% names(estimates)) estimates[["sigma"]] - config$margin_sigma else NA_real_
  )
  if (!any(is.finite(err))) NA_real_ else sqrt(mean(err[is.finite(err)]^2))
}

jss_misspec_margin_ci <- function(fit) {
  empty <- list(mu_covered_95 = NA, sigma_covered_95 = NA)
  if (!inherits(fit, "gamlss.longitudinal")) {
    return(empty)
  }
  ci <- tryCatch(
    stats::confint(fit, parm = c("mu.intercept", "sigma.intercept"), level = 0.95),
    error = function(e) NULL
  )
  if (is.null(ci)) {
    return(empty)
  }
  empty$mu_covered_95 <- isTRUE(ci["mu.intercept", 1] <= log(2) && ci["mu.intercept", 2] >= log(2))
  empty$sigma_covered_95 <- isTRUE(ci["sigma.intercept", 1] <= log(0.5) && ci["sigma.intercept", 2] >= log(0.5))
  empty
}

jss_misspec_run_checkpoints <- function(grid, config, checkpoint_dir) {
  dir.create(checkpoint_dir, recursive = TRUE, showWarnings = FALSE)
  pending <- jss_misspec_pending_grid(grid, checkpoint_dir)
  if (!nrow(pending)) {
    return(invisible(0L))
  }

  scenario_cols <- c("generating_copula", "tau_label", "target_tau", "n_subject", "n_time", "rep")
  scenario_key <- interaction(pending[scenario_cols], drop = TRUE, lex.order = TRUE)
  done <- 0L
  for (idx in split(seq_len(nrow(pending)), scenario_key)) {
    scenario <- pending[idx[1L], , drop = FALSE]
    dat_seed <- jss_misspec_dataset_seed(
      scenario$generating_copula,
      scenario$tau_label,
      scenario$n_subject,
      scenario$rep,
      config$seed
    )
    dat <- jss_misspec_simulate(
      generating_copula = scenario$generating_copula,
      target_tau = scenario$target_tau,
      n_subject = scenario$n_subject,
      config = config,
      seed = dat_seed
    )
    for (j in idx) {
      fit_row <- pending[j, , drop = FALSE]
      out <- jss_misspec_fit_one(dat, fit_row, config)
      utils::write.csv(out, fit_row$checkpoint_path, row.names = FALSE)
      done <- done + 1L
    }
  }
  invisible(done)
}

jss_misspec_read_checkpoints <- function(checkpoint_dir) {
  if (!dir.exists(checkpoint_dir)) {
    return(data.frame())
  }
  files <- list.files(checkpoint_dir, pattern = "[.]csv$", full.names = TRUE)
  if (!length(files)) {
    return(data.frame())
  }
  character_cols <- c(
    "fit_id", "stage", "generating_copula", "fitted_copula", "tau_label",
    "failure_type", "error", "warnings"
  )
  rows <- lapply(files, function(path) {
    header <- names(utils::read.csv(path, nrows = 0, stringsAsFactors = FALSE))
    classes <- stats::setNames(rep(NA_character_, length(header)), header)
    classes[intersect(character_cols, header)] <- "character"
    utils::read.csv(path, stringsAsFactors = FALSE, colClasses = classes)
  })
  jss_misspec_row_bind(rows)
}

jss_misspec_add_deltas <- function(results) {
  if (!nrow(results)) {
    return(results)
  }
  keys <- c("generating_copula", "tau_label", "n_subject", "rep")
  correct <- results[results$correctly_specified %in% TRUE, , drop = FALSE]
  keep <- c(keys, "joint_loglik", "copula_loglik", "tau_abs_error", "margin_param_rmse", "elapsed_sec")
  correct <- correct[keep]
  names(correct)[-(seq_along(keys))] <- paste0("correct_", names(correct)[-(seq_along(keys))])
  merged <- merge(results, correct, by = keys, all.x = TRUE, sort = FALSE)
  merged$delta_joint_loglik_vs_correct <- merged$joint_loglik - merged$correct_joint_loglik
  merged$delta_copula_loglik_vs_correct <- merged$copula_loglik - merged$correct_copula_loglik
  merged$delta_tau_abs_error_vs_correct <- merged$tau_abs_error - merged$correct_tau_abs_error
  merged$delta_margin_rmse_vs_correct <- merged$margin_param_rmse - merged$correct_margin_param_rmse
  merged$delta_elapsed_sec_vs_correct <- merged$elapsed_sec - merged$correct_elapsed_sec
  merged
}

jss_misspec_summary <- function(results) {
  if (!nrow(results)) {
    return(data.frame())
  }
  group_cols <- c("generating_copula", "fitted_copula", "tau_label", "target_tau", "n_subject")
  metric_cols <- c(
    "success", "converged", "joint_loglik", "copula_loglik", "aic_joint", "bic_joint",
    "elapsed_sec", "mu_bias", "sigma_bias", "margin_param_rmse", "tau_abs_error",
    "benchmark_neg_log_score", "benchmark_pit_mean_abs_error",
    "delta_joint_loglik_vs_correct", "delta_copula_loglik_vs_correct",
    "delta_tau_abs_error_vs_correct", "delta_margin_rmse_vs_correct",
    "delta_elapsed_sec_vs_correct"
  )
  split_key <- interaction(results[group_cols], drop = TRUE, lex.order = TRUE)
  rows <- lapply(split(seq_len(nrow(results)), split_key), function(idx) {
    row <- results[idx[1L], group_cols, drop = FALSE]
    row$n_attempted <- length(idx)
    row$success_rate <- mean(results$success[idx] %in% TRUE)
    row$convergence_rate <- mean(results$converged[idx] %in% TRUE)
    for (metric in setdiff(metric_cols, c("success", "converged"))) {
      row[[paste0("mean_", metric)]] <- mean(results[[metric]][idx], na.rm = TRUE)
    }
    row
  })
  jss_misspec_row_bind(rows)
}

jss_misspec_selection <- function(results) {
  if (!nrow(results)) {
    return(data.frame())
  }
  keys <- c("generating_copula", "tau_label", "target_tau", "n_subject", "rep")
  split_key <- interaction(results[keys], drop = TRUE, lex.order = TRUE)
  rows <- lapply(split(seq_len(nrow(results)), split_key), function(idx) {
    dat <- results[idx, , drop = FALSE]
    usable <- dat[dat$success %in% TRUE, , drop = FALSE]
    base <- dat[1L, keys, drop = FALSE]
    if (!nrow(usable)) {
      base$best_aic_copula <- NA_character_
      base$best_bic_copula <- NA_character_
      base$best_loglik_copula <- NA_character_
      return(base)
    }
    base$best_aic_copula <- jss_misspec_best_copula(usable, "aic_joint", decreasing = FALSE)
    base$best_bic_copula <- jss_misspec_best_copula(usable, "bic_joint", decreasing = FALSE)
    base$best_loglik_copula <- jss_misspec_best_copula(usable, "joint_loglik", decreasing = TRUE)
    base
  })
  selected <- jss_misspec_row_bind(rows)
  group_cols <- c("generating_copula", "tau_label", "target_tau", "n_subject")
  split_key <- interaction(selected[group_cols], drop = TRUE, lex.order = TRUE)
  out <- lapply(split(seq_len(nrow(selected)), split_key), function(idx) {
    row <- selected[idx[1L], group_cols, drop = FALSE]
    row$n_reps <- length(idx)
    row$aic_correct_selection_rate <- mean(selected$best_aic_copula[idx] == selected$generating_copula[idx], na.rm = TRUE)
    row$bic_correct_selection_rate <- mean(selected$best_bic_copula[idx] == selected$generating_copula[idx], na.rm = TRUE)
    row$loglik_correct_selection_rate <- mean(selected$best_loglik_copula[idx] == selected$generating_copula[idx], na.rm = TRUE)
    row
  })
  jss_misspec_row_bind(out)
}

jss_misspec_best_copula <- function(results, metric, decreasing = FALSE) {
  value <- results[[metric]]
  ok <- is.finite(value)
  if (!any(ok)) {
    return(NA_character_)
  }
  usable <- results[ok, , drop = FALSE]
  value <- value[ok]
  if (isTRUE(decreasing)) {
    usable$fitted_copula[which.max(value)]
  } else {
    usable$fitted_copula[which.min(value)]
  }
}

jss_misspec_review_gate <- function(results, grid, paths) {
  expected <- nrow(grid)
  actual <- nrow(results)
  combo_cols <- c("generating_copula", "fitted_copula", "tau_label", "n_subject")
  expected_combos <- unique(grid[combo_cols])
  actual_combos <- if (nrow(results)) unique(results[combo_cols]) else results[combo_cols]

  high_hard <- results[
    results$tau_label == "high" &
      results$generating_copula %in% c("C", "G", "J", "t") &
      results$fitted_copula %in% c("C", "G", "J", "t"),
    ,
    drop = FALSE
  ]
  correct <- results[results$correctly_specified %in% TRUE, , drop = FALSE]
  wrong <- results[!(results$correctly_specified %in% TRUE), , drop = FALSE]
  correct_success <- if (nrow(correct)) mean(correct$success %in% TRUE) else NA_real_
  wrong_success <- if (nrow(wrong)) mean(wrong$success %in% TRUE) else NA_real_

  data.frame(
    check = c(
      "row_count",
      "full_combination_coverage",
      "hard_high_tau_convergence_recorded",
      "finite_information_criteria",
      "correct_copula_not_systematically_worse",
      "summary_artifacts_rendered"
    ),
    status = c(
      if (actual == expected) "pass" else "review",
      if (nrow(actual_combos) == nrow(expected_combos)) "pass" else "review",
      if (!nrow(high_hard) || any(!is.na(high_hard$converged))) "pass" else "review",
      if (!nrow(results) || mean(is.finite(results$aic_joint) & is.finite(results$bic_joint), na.rm = TRUE) > 0.5) "pass" else "review",
      if (is.na(correct_success) || is.na(wrong_success) || correct_success + 0.05 >= wrong_success) "pass" else "review",
      if (all(file.exists(c(paths$summary, paths$selection, paths$heatmap, paths$convergence)))) "pass" else "review"
    ),
    detail = c(
      paste(actual, "of", expected, "fit rows available"),
      paste(nrow(actual_combos), "of", nrow(expected_combos), "scenario combinations represented"),
      paste(nrow(high_hard), "high-tau hard-family fit rows available"),
      paste("finite AIC/BIC rate:", round(mean(is.finite(results$aic_joint) & is.finite(results$bic_joint), na.rm = TRUE), 3)),
      paste("correct success:", round(correct_success, 3), "wrong success:", round(wrong_success, 3)),
      "summary, selection, heatmap, and convergence artifacts checked"
    ),
    stringsAsFactors = FALSE
  )
}

jss_misspec_write_heatmap <- function(summary, path) {
  if (!nrow(summary)) {
    return(jss_misspec_write_empty_plot(path, "Gamma copula mis-specification heatmap"))
  }
  p <- ggplot2::ggplot(
    summary,
    ggplot2::aes(x = fitted_copula, y = generating_copula, fill = mean_delta_joint_loglik_vs_correct)
  ) +
    ggplot2::geom_tile(color = "white", linewidth = 0.2) +
    ggplot2::facet_grid(tau_label ~ n_subject, labeller = ggplot2::label_both) +
    ggplot2::scale_fill_gradient2(low = "#b2182b", mid = "white", high = "#2166ac", midpoint = 0, na.value = "grey85") +
    ggplot2::labs(x = "Fitted copula", y = "Generating copula", fill = "Delta joint logLik") +
    ggplot2::theme_minimal(base_size = 10)
  ggplot2::ggsave(path, p, width = 10, height = 6, dpi = 320, bg = "white")
  path
}

jss_misspec_write_metric_heatmap <- function(summary, path, metric, title, fill_label) {
  if (!nrow(summary)) {
    return(jss_misspec_write_empty_plot(path, title))
  }
  p <- ggplot2::ggplot(
    summary,
    ggplot2::aes_string(x = "fitted_copula", y = "generating_copula", fill = metric)
  ) +
    ggplot2::geom_tile(color = "white", linewidth = 0.2) +
    ggplot2::facet_grid(tau_label ~ n_subject, labeller = ggplot2::label_both) +
    ggplot2::scale_fill_gradient(low = "#f7fbff", high = "#08519c", na.value = "grey85") +
    ggplot2::labs(title = title, x = "Fitted copula", y = "Generating copula", fill = fill_label) +
    ggplot2::theme_minimal(base_size = 10)
  ggplot2::ggsave(path, p, width = 10, height = 6, dpi = 320, bg = "white")
  path
}

jss_misspec_summary_panel <- function(plot_data, metric, title, fill_label, scale) {
  ggplot2::ggplot(
    plot_data,
    ggplot2::aes_string(x = "fitted_copula", y = "generating_copula", fill = metric)
  ) +
    ggplot2::geom_tile(color = "white", linewidth = 0.25) +
    scale +
    ggplot2::coord_equal() +
    ggplot2::labs(title = title, x = "Fitted copula", y = "Generating copula", fill = fill_label) +
    ggplot2::theme_minimal(base_size = 9) +
    ggplot2::theme(
      plot.title = ggplot2::element_text(face = "bold", size = 10),
      legend.position = "bottom",
      legend.key.width = grid::unit(1.25, "cm"),
      panel.grid = ggplot2::element_blank()
    )
}

jss_misspec_write_paper_summary_heatmap <- function(summary, path, n_subject = 500L, tau_label = "high") {
  title <- paste0("Gamma copula mis-specification summary: n = ", n_subject, ", tau = ", tau_label)
  if (!nrow(summary)) {
    return(jss_misspec_write_empty_plot(path, title))
  }

  plot_data <- summary[
    summary$n_subject == n_subject & summary$tau_label == tau_label,
    ,
    drop = FALSE
  ]
  if (!nrow(plot_data)) {
    return(jss_misspec_write_empty_plot(path, title))
  }

  copulas <- jss_misspec_copulas()
  plot_data$generating_copula <- factor(plot_data$generating_copula, levels = rev(copulas))
  plot_data$fitted_copula <- factor(plot_data$fitted_copula, levels = copulas)

  likelihood_panel <- jss_misspec_summary_panel(
    plot_data,
    metric = "mean_delta_joint_loglik_vs_correct",
    title = "Delta joint log-likelihood",
    fill_label = "Delta logLik",
    scale = ggplot2::scale_fill_gradient2(
      low = "#b2182b",
      mid = "white",
      high = "#2166ac",
      midpoint = 0,
      na.value = "grey85"
    )
  )
  margin_panel <- jss_misspec_summary_panel(
    plot_data,
    metric = "mean_margin_param_rmse",
    title = "Margin parameter RMSE",
    fill_label = "RMSE",
    scale = ggplot2::scale_fill_gradient(low = "#f7fbff", high = "#08519c", na.value = "grey85")
  )
  tau_panel <- jss_misspec_summary_panel(
    plot_data,
    metric = "mean_tau_abs_error",
    title = "Kendall tau absolute error",
    fill_label = "Abs. error",
    scale = ggplot2::scale_fill_gradient(low = "#f7fbff", high = "#08519c", na.value = "grey85")
  )

  grDevices::png(path, width = 11, height = 4.2, units = "in", res = 320, bg = "white")
  on.exit(grDevices::dev.off(), add = TRUE)
  grid::grid.newpage()
  grid::pushViewport(grid::viewport(layout = grid::grid.layout(1, 3)))
  print(likelihood_panel, vp = grid::viewport(layout.pos.row = 1, layout.pos.col = 1))
  print(margin_panel, vp = grid::viewport(layout.pos.row = 1, layout.pos.col = 2))
  print(tau_panel, vp = grid::viewport(layout.pos.row = 1, layout.pos.col = 3))
  grid::popViewport()
  path
}

jss_misspec_write_convergence_plot <- function(summary, path) {
  if (!nrow(summary)) {
    return(jss_misspec_write_empty_plot(path, "Gamma copula mis-specification convergence"))
  }
  p <- ggplot2::ggplot(
    summary,
    ggplot2::aes(x = fitted_copula, y = convergence_rate, fill = fitted_copula)
  ) +
    ggplot2::geom_col(width = 0.7, show.legend = FALSE) +
    ggplot2::facet_grid(generating_copula + tau_label ~ n_subject, labeller = ggplot2::label_both) +
    ggplot2::coord_cartesian(ylim = c(0, 1)) +
    ggplot2::scale_y_continuous(labels = scales::percent_format(accuracy = 1)) +
    ggplot2::labs(x = "Fitted copula", y = "Convergence rate") +
    ggplot2::theme_minimal(base_size = 9)
  ggplot2::ggsave(path, p, width = 11, height = 8, dpi = 320, bg = "white")
  path
}

jss_misspec_write_empty_plot <- function(path, title) {
  p <- ggplot2::ggplot(data.frame(x = 1, y = 1), ggplot2::aes(x, y)) +
    ggplot2::geom_blank() +
    ggplot2::labs(title = title, subtitle = "No result rows available yet.") +
    ggplot2::theme_void(base_size = 11)
  ggplot2::ggsave(path, p, width = 7, height = 4, dpi = 220, bg = "white")
  path
}

jss_run_07_gamma_copula_misspecification <- function(settings, stage = jss_misspec_stage(settings$profile)) {
  paths <- jss_misspec_paths(settings)
  config <- jss_misspec_config(settings, stage = stage)
  checkpoint_dir <- jss_misspec_checkpoint_dir(paths, config$stage)
  grid <- jss_misspec_grid(config)

  utils::write.csv(grid, paths$grid, row.names = FALSE)
  jss_misspec_run_checkpoints(grid, config, checkpoint_dir)

  results <- jss_misspec_read_checkpoints(checkpoint_dir)
  results <- results[results$fit_id %in% grid$fit_id, , drop = FALSE]
  results <- jss_misspec_add_deltas(results)
  summary <- jss_misspec_summary(results)
  selection <- jss_misspec_selection(results)

  utils::write.csv(results, paths$results, row.names = FALSE)
  utils::write.csv(summary, paths$summary, row.names = FALSE)
  utils::write.csv(selection, paths$selection, row.names = FALSE)
  jss_misspec_write_heatmap(summary, paths$heatmap)
  jss_misspec_write_metric_heatmap(
    summary,
    paths$margin_rmse_heatmap,
    metric = "mean_margin_param_rmse",
    title = "Gamma margin parameter RMSE by fitted copula",
    fill_label = "Margin RMSE"
  )
  jss_misspec_write_metric_heatmap(
    summary,
    paths$tau_error_heatmap,
    metric = "mean_tau_abs_error",
    title = "Kendall tau absolute error by fitted copula",
    fill_label = "Tau abs. error"
  )
  jss_misspec_write_paper_summary_heatmap(summary, paths$paper_summary_heatmap)
  jss_misspec_write_convergence_plot(summary, paths$convergence)
  review <- jss_misspec_review_gate(results, grid, paths)
  utils::write.csv(review, paths$review, row.names = FALSE)

  list(
    module_id = "07-gamma-copula-misspecification",
    title = "Gamma copula mis-specification simulation",
    status = "current",
    stage = config$stage,
    data = c(paths$grid, paths$results),
    tables = c(paths$summary, paths$selection, paths$review),
    figures = c(
      paths$heatmap,
      paths$margin_rmse_heatmap,
      paths$tau_error_heatmap,
      paths$paper_summary_heatmap,
      paths$convergence
    ),
    notes = paste("Checkpointed", config$stage, "run with", nrow(grid), "planned fits.")
  )
}
