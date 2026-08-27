jss_joint_module_id <- function() {
  "03-joint-vs-separate-optimization"
}

jss_joint_or <- function(x, y) {
  if (is.null(x) || length(x) == 0L || is.na(x[[1L]])) y else x[[1L]]
}

`%||%` <- function(a, b) if (!is.null(a)) a else b

jss_joint_env_int <- function(name, default, minimum = 1L) {
  value <- suppressWarnings(as.integer(Sys.getenv(name, unset = as.character(default))))
  if (!is.finite(value) || value < minimum) {
    return(as.integer(default))
  }
  as.integer(value)
}

jss_joint_env_num <- function(name, default, minimum = 0) {
  raw <- Sys.getenv(name, unset = as.character(default))
  if (tolower(raw) %in% c("inf", "infinity")) {
    return(Inf)
  }
  value <- suppressWarnings(as.numeric(raw))
  if (!is.finite(value) || value < minimum) {
    return(default)
  }
  value
}

jss_joint_output_paths <- function(settings) {
  module_id <- jss_joint_module_id()
  list(
    candidate_selection = file.path(
      settings$data_dir,
      paste0(module_id, "-candidate-selection.csv")
    ),
    results = file.path(
      settings$data_dir,
      paste0(module_id, "-results.csv")
    ),
    summary = file.path(
      settings$tables_dir,
      paste0(module_id, "-summary.csv")
    ),
    metric_wins = file.path(
      settings$tables_dir,
      paste0(module_id, "-metric-wins.csv")
    ),
    hypothesis_summary = file.path(
      settings$tables_dir,
      paste0(module_id, "-hypothesis-summary.csv")
    ),
    deltas_figure = file.path(
      settings$figures_dir,
      paste0(module_id, "-deltas.png")
    ),
    metric_dashboard = file.path(
      settings$figures_dir,
      paste0(module_id, "-metric-dashboard.png")
    )
  )
}

jss_joint_candidate_sources <- function(settings) {
  base <- file.path(
    settings$root,
    "results",
    "jss-exploratory",
    "03-joint-vs-separate-optimization"
  )
  discrete_base <- file.path(
    settings$root,
    "results",
    "jss-exploratory",
    "02-discrete-delaporte-clayton"
  )
  legacy <- file.path(settings$root, "results")

  c(
    file.path(base, "rs_joint_balanced_tau_grid", "balanced_tau_balanced_candidates.csv"),
    file.path(base, "rs_joint_balanced_tau_grid", "balanced_tau_joint_vs_separate_summary.csv"),
    file.path(base, "rs_joint_stress_effect_search", "stress_balanced_candidates.csv"),
    file.path(base, "rs_joint_stress_effect_search", "stress_max_loglik_candidates.csv"),
    file.path(base, "rs_joint_stress_effect_search", "stress_joint_vs_separate_summary.csv"),
    file.path(base, "rs_joint_controlled_effect_sweep", "controlled_top_candidates.csv"),
    file.path(base, "rs_joint_controlled_effect_sweep", "controlled_joint_vs_separate_summary.csv"),
    file.path(base, "rs_joint_broad_copula_screen", "broad_top_candidates.csv"),
    file.path(base, "rs_joint_broad_copula_screen", "top_cases_loglik_logscore_all_rmse_improved.csv"),
    file.path(base, "rs_joint_broad_copula_screen", "broad_joint_vs_separate_summary.csv"),
    file.path(base, "rs_joint_tau_calibrated_pilot", "tau_calibrated_balanced_candidates.csv"),
    file.path(base, "rs_joint_tau_calibrated_pilot", "tau_calibrated_joint_vs_separate_summary.csv"),
    file.path(base, "rs_joint_selected_copula_report", "selected_joint_vs_separate_summary.csv"),
    file.path(discrete_base, "rs_clayton_showcase_fit_quality", "joint_vs_separate_delta_summary.csv"),
    file.path(legacy, "rs_joint_balanced_tau_grid", "balanced_tau_balanced_candidates.csv"),
    file.path(legacy, "rs_joint_stress_effect_search", "stress_balanced_candidates.csv"),
    file.path(legacy, "rs_joint_controlled_effect_sweep", "controlled_top_candidates.csv"),
    file.path(legacy, "rs_joint_broad_copula_screen", "broad_top_candidates.csv"),
    file.path(legacy, "rs_joint_selected_copula_report", "selected_joint_vs_separate_summary.csv"),
    file.path(legacy, "rs_clayton_showcase_fit_quality", "joint_vs_separate_delta_summary.csv")
  )
}

jss_joint_read_candidate_source <- function(path) {
  if (!file.exists(path)) {
    return(NULL)
  }
  x <- tryCatch(utils::read.csv(path, stringsAsFactors = FALSE), error = function(e) NULL)
  if (is.null(x) || nrow(x) == 0L) {
    return(NULL)
  }
  x$source_file <- normalizePath(path, winslash = "/", mustWork = FALSE)
  x$source_name <- tools::file_path_sans_ext(basename(path))
  x$source_row <- seq_len(nrow(x))
  x
}

jss_joint_bind_rows <- function(rows) {
  rows <- rows[!vapply(rows, is.null, logical(1))]
  if (length(rows) == 0L) {
    return(data.frame())
  }
  all_names <- unique(unlist(lapply(rows, names), use.names = FALSE))
  rows <- lapply(rows, function(x) {
    missing <- setdiff(all_names, names(x))
    for (nm in missing) {
      x[[nm]] <- NA
    }
    x[all_names]
  })
  out <- do.call(rbind, rows)
  rownames(out) <- NULL
  out
}

jss_joint_numeric <- function(x) {
  suppressWarnings(as.numeric(x))
}

jss_joint_scaled <- function(x, higher_is_better = TRUE) {
  x <- jss_joint_numeric(x)
  if (!higher_is_better) {
    x <- -x
  }
  ok <- is.finite(x)
  out <- rep(0, length(x))
  if (sum(ok) < 2L) {
    out[ok] <- x[ok]
    return(out)
  }
  spread <- stats::sd(x[ok])
  if (!is.finite(spread) || spread == 0) {
    spread <- max(abs(x[ok]), na.rm = TRUE)
  }
  if (!is.finite(spread) || spread == 0) {
    spread <- 1
  }
  out[ok] <- (x[ok] - stats::median(x[ok], na.rm = TRUE)) / spread
  out
}

jss_joint_first_present <- function(x, names, default = NA_real_) {
  idx <- match(names, names(x))
  idx <- idx[!is.na(idx)]
  if (length(idx) == 0L) {
    return(rep(default, nrow(x)))
  }
  x[[idx[[1L]]]]
}

jss_joint_classify_candidate <- function(x) {
  test_gain <- jss_joint_numeric(jss_joint_first_present(x, "delta_test_log_score_per_obs_mean"))
  train_gain <- jss_joint_numeric(jss_joint_first_present(x, "delta_train_joint_loglik_mean"))
  theta_gain <- -jss_joint_numeric(jss_joint_first_present(x, "delta_train_rmse_theta_mean"))
  tau_gain <- -jss_joint_numeric(jss_joint_first_present(x, "delta_train_rmse_tau_mean"))
  mu_gain <- -jss_joint_numeric(jss_joint_first_present(x, "delta_train_rmse_mu_mean"))
  sigma_gain <- -jss_joint_numeric(jss_joint_first_present(x, "delta_train_rmse_sigma_mean"))

  dependence_gain <- pmax(theta_gain, tau_gain, na.rm = TRUE)
  dependence_gain[!is.finite(dependence_gain)] <- NA_real_
  margin_gain <- pmax(mu_gain, sigma_gain, na.rm = TRUE)
  margin_gain[!is.finite(margin_gain)] <- NA_real_

  type <- rep("other", nrow(x))
  type[
    is.finite(test_gain) & test_gain > 0 &
      is.finite(train_gain) & train_gain > 0 &
      is.finite(dependence_gain) & dependence_gain > 0
  ] <- "joint_win"
  type[
    abs(test_gain) <= 0.005 &
      abs(train_gain) <= 1 &
      (is.na(dependence_gain) | abs(dependence_gain) <= 0.02)
  ] <- "tie_control"
  type[
    is.finite(train_gain) & train_gain > 0 &
      (is.finite(test_gain) & test_gain <= 0 |
        is.finite(margin_gain) & margin_gain < 0 |
        is.finite(dependence_gain) & dependence_gain < 0)
  ] <- "cautionary"
  type
}

jss_joint_candidate_signature <- function(x) {
  parts <- lapply(
    c(
      "copula", "scenario", "time_shape", "profile_id", "profile_type",
      "tau_path", "theta_shape", "tau_max", "sigma_strength", "mu_strength",
      "theta_strength", "n", "time_points"
    ),
    function(nm) {
      if (nm %in% names(x)) as.character(x[[nm]]) else rep("", nrow(x))
    }
  )
  parts <- as.data.frame(parts, stringsAsFactors = FALSE)
  parts[is.na(parts)] <- ""
  apply(parts, 1L, paste, collapse = "|")
}

jss_joint_candidate_selection <- function(settings) {
  sources <- unique(jss_joint_candidate_sources(settings))
  candidates <- jss_joint_bind_rows(lapply(sources, jss_joint_read_candidate_source))

  required <- c(
    "candidate_id", "review_role", "candidate_score", "source_name",
    "source_file", "source_row", "candidate_signature", "copula", "scenario", "time_shape",
    "profile_id", "profile_type", "tau_path", "theta_shape", "tau_max",
    "sigma_strength", "mu_strength", "theta_strength", "n", "time_points",
    "delta_train_joint_loglik_mean", "delta_test_log_score_per_obs_mean",
    "delta_train_rmse_mu_mean", "delta_train_rmse_sigma_mean",
    "delta_train_rmse_theta_mean", "delta_train_rmse_tau_mean", "n_rows"
  )

  if (nrow(candidates) == 0L) {
    out <- as.data.frame(stats::setNames(rep(list(character()), length(required)), required))
    return(out)
  }

  if ("source_name" %in% names(candidates) && "copula" %in% names(candidates)) {
    missing_copula <- is.na(candidates$copula) | !nzchar(as.character(candidates$copula))
    candidates$copula[missing_copula & candidates$source_name == "joint_vs_separate_delta_summary"] <- "C"
  }

  candidates$review_role <- jss_joint_classify_candidate(candidates)
  candidates$candidate_score <-
    jss_joint_scaled(jss_joint_first_present(candidates, "delta_test_log_score_per_obs_mean")) +
    jss_joint_scaled(jss_joint_first_present(candidates, "delta_train_joint_loglik_mean")) +
    jss_joint_scaled(jss_joint_first_present(candidates, "delta_train_rmse_theta_mean"), higher_is_better = FALSE) +
    jss_joint_scaled(jss_joint_first_present(candidates, "delta_train_rmse_tau_mean"), higher_is_better = FALSE) +
    0.5 * jss_joint_scaled(jss_joint_first_present(candidates, "delta_train_rmse_mu_mean"), higher_is_better = FALSE) +
    0.5 * jss_joint_scaled(jss_joint_first_present(candidates, "delta_train_rmse_sigma_mean"), higher_is_better = FALSE)

  n_rows <- jss_joint_numeric(jss_joint_first_present(candidates, "n_rows"))
  candidates$candidate_score[is.finite(n_rows) & n_rows < 3L] <-
    candidates$candidate_score[is.finite(n_rows) & n_rows < 3L] - 2
  candidates$candidate_signature <- jss_joint_candidate_signature(candidates)
  candidates <- candidates[
    order(-candidates$candidate_score, candidates$source_name, candidates$source_row),
    ,
    drop = FALSE
  ]
  candidates <- candidates[!duplicated(candidates$candidate_signature), , drop = FALSE]
  rownames(candidates) <- NULL

  chosen <- list()
  take_role <- function(role, n) {
    idx <- which(candidates$review_role == role)
    if (length(idx) == 0L || n <= 0L) {
      return(NULL)
    }
    idx <- idx[order(-candidates$candidate_score[idx], candidates$source_name[idx], candidates$source_row[idx])]
    candidates[idx[seq_len(min(n, length(idx)))], , drop = FALSE]
  }
  chosen[[length(chosen) + 1L]] <- take_role("joint_win", 5L)
  chosen[[length(chosen) + 1L]] <- take_role("tie_control", 1L)
  chosen[[length(chosen) + 1L]] <- take_role("cautionary", 1L)

  selected <- jss_joint_bind_rows(chosen)
  if (nrow(selected) < 7L) {
    remaining <- candidates[
      !(seq_len(nrow(candidates)) %in% match(
        paste(selected$source_file, selected$source_row),
        paste(candidates$source_file, candidates$source_row)
      )),
      ,
      drop = FALSE
    ]
    remaining <- remaining[
      order(-remaining$candidate_score, remaining$source_name, remaining$source_row),
      ,
      drop = FALSE
    ]
    selected <- jss_joint_bind_rows(list(
      selected,
      remaining[seq_len(min(7L - nrow(selected), nrow(remaining))), , drop = FALSE]
    ))
  }

  selected <- selected[seq_len(min(7L, nrow(selected))), , drop = FALSE]
  selected$candidate_id <- sprintf("jvs-%02d", seq_len(nrow(selected)))

  missing <- setdiff(required, names(selected))
  for (nm in missing) {
    selected[[nm]] <- NA
  }
  selected <- selected[, required, drop = FALSE]
  rownames(selected) <- NULL
  selected
}

jss_joint_case_definitions <- function() {
  data.frame(
    candidate_id = sprintf("jvs-%02d", 1:7),
    case_id = sprintf("JVS%02d", 1:7),
    hypothesis_role = c(
      "theta_off_control",
      "weak_dependence",
      "central_joint_win",
      "fewer_time_points_same_rows",
      "weaker_shared_shape",
      "larger_sample_same_T",
      "discrete_analogue"
    ),
    family = c("NO", "NO", "NO", "NO", "NO", "NO", "NBI"),
    copula = "N",
    n = c(40L, 40L, 40L, 150L, 40L, 100L, 40L),
    time_points = c(75L, 75L, 75L, 20L, 75L, 75L, 75L),
    total_observations = c(40L * 75L, 40L * 75L, 40L * 75L, 150L * 20L, 40L * 75L, 100L * 75L, 40L * 75L),
    mu_strength = 0,
    sigma_strength = c(1.5, 1.5, 1.5, 1.5, 0.5, 1.5, 1.5),
    theta_strength = c(0, 1.5, 3.5, 3.5, 3.5, 3.5, 3.5),
    time_shape = "sigmoid",
    purpose = c(
      "Tests whether theta covariates are needed.",
      "Correlation-strength contrast.",
      "Main high-signal continuous case.",
      "Same total rows as JVS03, fewer time points.",
      "Shared shape/theta covariate contrast.",
      "Sample-size contrast against JVS03.",
      "Discrete-vs-continuous contrast."
    ),
    review_role = c("tie_control", "joint_win", "joint_win", "contrast", "contrast", "contrast", "contrast"),
    candidate_score = NA_real_,
    source_name = "fixed_seven_case_design",
    source_file = NA_character_,
    source_row = seq_len(7L),
    candidate_signature = NA_character_,
    scenario = NA_character_,
    profile_id = NA_character_,
    profile_type = NA_character_,
    tau_path = NA_character_,
    theta_shape = "sigmoid",
    tau_max = NA_real_,
    delta_train_joint_loglik_mean = NA_real_,
    delta_test_log_score_per_obs_mean = NA_real_,
    delta_train_rmse_mu_mean = NA_real_,
    delta_train_rmse_sigma_mean = NA_real_,
    delta_train_rmse_theta_mean = NA_real_,
    delta_train_rmse_tau_mean = NA_real_,
    n_rows = NA_integer_,
    stringsAsFactors = FALSE
  )
}

jss_joint_simulation_settings <- function(settings) {
  if (settings$profile %in% c("expanded", "full")) {
    return(list(
      reps = jss_joint_env_int("GAMLSS_LONGITUDINAL_JSS_JVS_REPS", 20L),
      max_outer_iter = jss_joint_env_int("GAMLSS_LONGITUDINAL_JSS_JVS_MAX_OUTER", 250L),
      max_inner_iter = jss_joint_env_int("GAMLSS_LONGITUDINAL_JSS_JVS_MAX_INNER", 100L),
      max_elapsed_sec = jss_joint_env_num("GAMLSS_LONGITUDINAL_JSS_JVS_MAX_ELAPSED", Inf),
      method = toupper(Sys.getenv("GAMLSS_LONGITUDINAL_JSS_JVS_METHOD", unset = "RS")),
      outer_stop_crit = jss_joint_env_num("GAMLSS_LONGITUDINAL_JSS_JVS_OUTER_STOP", 0.01),
      inner_stop_crit = jss_joint_env_num("GAMLSS_LONGITUDINAL_JSS_JVS_INNER_STOP", 0.002),
      warm_start_joint = as.logical(jss_joint_env_int("GAMLSS_LONGITUDINAL_JSS_JVS_WARM_START_JOINT", 1L, minimum = 0L)),
      warm_start_joint_iter = jss_joint_env_int("GAMLSS_LONGITUDINAL_JSS_JVS_WARM_START_JOINT_ITER", 5L, minimum = 0L),
      discrete_score_method = Sys.getenv("GAMLSS_LONGITUDINAL_JSS_JVS_DISCRETE_SCORE_METHOD", unset = "analytical"),
      cg_max_delta = jss_joint_env_num("GAMLSS_LONGITUDINAL_JSS_JVS_CG_MAX_DELTA", 0.5, minimum = 1e-8),
      cg_line_search = Sys.getenv("GAMLSS_LONGITUDINAL_JSS_JVS_CG_LINE_SEARCH", unset = "best"),
      cg_gradient_method = Sys.getenv("GAMLSS_LONGITUDINAL_JSS_JVS_CG_GRADIENT_METHOD", unset = "forward"),
      cg_hessian_method = Sys.getenv("GAMLSS_LONGITUDINAL_JSS_JVS_CG_HESSIAN_METHOD", unset = "auto"),
      cg_raw_loglik_drop_tol = jss_joint_env_num("GAMLSS_LONGITUDINAL_JSS_JVS_CG_RAW_LOGLIK_DROP_TOL", 10, minimum = 0),
      cg_max_line_search_evals = jss_joint_env_int("GAMLSS_LONGITUDINAL_JSS_JVS_CG_MAX_LINE_SEARCH_EVALS", 60L, minimum = 0L),
      variogram_nsim = jss_joint_env_int("GAMLSS_LONGITUDINAL_JSS_JVS_VARIOGRAM_NSIM", 100L)
    ))
  }

  list(
    reps = jss_joint_env_int("GAMLSS_LONGITUDINAL_JSS_JVS_REPS", 1L),
    max_outer_iter = jss_joint_env_int("GAMLSS_LONGITUDINAL_JSS_JVS_MAX_OUTER", 3L),
    max_inner_iter = jss_joint_env_int("GAMLSS_LONGITUDINAL_JSS_JVS_MAX_INNER", 3L),
    max_elapsed_sec = jss_joint_env_num("GAMLSS_LONGITUDINAL_JSS_JVS_MAX_ELAPSED", 20),
    method = toupper(Sys.getenv("GAMLSS_LONGITUDINAL_JSS_JVS_METHOD", unset = "RS")),
    outer_stop_crit = jss_joint_env_num("GAMLSS_LONGITUDINAL_JSS_JVS_OUTER_STOP", 0.01),
    inner_stop_crit = jss_joint_env_num("GAMLSS_LONGITUDINAL_JSS_JVS_INNER_STOP", 0.002),
    warm_start_joint = as.logical(jss_joint_env_int("GAMLSS_LONGITUDINAL_JSS_JVS_WARM_START_JOINT", 1L, minimum = 0L)),
    warm_start_joint_iter = jss_joint_env_int("GAMLSS_LONGITUDINAL_JSS_JVS_WARM_START_JOINT_ITER", 5L, minimum = 0L),
    discrete_score_method = Sys.getenv("GAMLSS_LONGITUDINAL_JSS_JVS_DISCRETE_SCORE_METHOD", unset = "analytical"),
    cg_max_delta = jss_joint_env_num("GAMLSS_LONGITUDINAL_JSS_JVS_CG_MAX_DELTA", 0.5, minimum = 1e-8),
    cg_line_search = Sys.getenv("GAMLSS_LONGITUDINAL_JSS_JVS_CG_LINE_SEARCH", unset = "best"),
    cg_gradient_method = Sys.getenv("GAMLSS_LONGITUDINAL_JSS_JVS_CG_GRADIENT_METHOD", unset = "forward"),
    cg_hessian_method = Sys.getenv("GAMLSS_LONGITUDINAL_JSS_JVS_CG_HESSIAN_METHOD", unset = "auto"),
    cg_raw_loglik_drop_tol = jss_joint_env_num("GAMLSS_LONGITUDINAL_JSS_JVS_CG_RAW_LOGLIK_DROP_TOL", 10, minimum = 0),
    cg_max_line_search_evals = jss_joint_env_int("GAMLSS_LONGITUDINAL_JSS_JVS_CG_MAX_LINE_SEARCH_EVALS", 60L, minimum = 0L),
    variogram_nsim = jss_joint_env_int("GAMLSS_LONGITUDINAL_JSS_JVS_VARIOGRAM_NSIM", 20L)
  )
}

jss_joint_clip_u <- function(u) {
  pmin(pmax(as.numeric(u), 1e-9), 1 - 1e-9)
}

jss_joint_rmse <- function(x) {
  sqrt(mean(as.numeric(x)^2, na.rm = TRUE))
}

jss_joint_time_shape <- function(t, shape) {
  if (identical(shape, "sigmoid")) {
    raw <- stats::plogis(8 * (t - 0.55))
    lo <- stats::plogis(8 * (0 - 0.55))
    hi <- stats::plogis(8 * (1 - 0.55))
    return((raw - lo) / max(1e-8, hi - lo))
  }
  t
}

jss_joint_margin_dist <- function(family) {
  do.call(get(family, envir = asNamespace("gamlss.dist"), mode = "function"), list())
}

jss_joint_case_value <- function(case, name, default) {
  if (!name %in% names(case)) {
    return(default)
  }
  value <- case[[name]][[1L]]
  if (length(value) == 0L || is.na(value)) {
    return(default)
  }
  value
}

jss_joint_linkinv <- function(margin_dist, parameter, eta) {
  linkinv <- margin_dist[[paste0(parameter, ".linkinv")]]
  if (is.function(linkinv)) linkinv(eta) else eta
}

jss_joint_truth_coefficients <- function(case) {
  if (identical(case$family[[1L]], "NBI")) {
    mu_intercept <- log(4)
    sigma_intercept <- log(0.5)
  } else {
    mu_intercept <- 1.2
    sigma_intercept <- log(0.75)
  }
  mu_intercept <- as.numeric(jss_joint_case_value(case, "mu_intercept", mu_intercept))
  sigma_intercept <- as.numeric(jss_joint_case_value(case, "sigma_intercept", sigma_intercept))
  theta_base_tau <- as.numeric(jss_joint_case_value(case, "theta_base_tau", 0.12))

  theta_intercept <- gamlss.longitudinal:::get_copula_dist(case$copula[[1L]])$copula_link$theta.linkfun(
    gamlss.longitudinal:::.copula_tau_to_par(case$copula[[1L]], theta_base_tau)
  )

  list(
    mu = c(
      intercept = mu_intercept,
      x_cont = 0.25 * case$mu_strength[[1L]],
      x_bin = 0.20 * case$mu_strength[[1L]],
      time_effect = 0.25 * case$mu_strength[[1L]]
    ),
    sigma = c(
      intercept = sigma_intercept,
      x_cont = 0.22 * case$sigma_strength[[1L]],
      x_bin = -0.18 * case$sigma_strength[[1L]],
      time_effect = 0.35 * case$sigma_strength[[1L]]
    ),
    theta = c(
      intercept = theta_intercept,
      x_cont = 0.22 * case$theta_strength[[1L]],
      x_bin = 0.18 * case$theta_strength[[1L]],
      t_pair_effect = 0.32 * case$theta_strength[[1L]]
    )
  )
}

jss_joint_make_covariates <- function(case, seed) {
  set.seed(seed)
  n <- case$n[[1L]]
  time_points <- case$time_points[[1L]]
  subject_df <- data.frame(
    x_cont = stats::rnorm(n),
    x_bin = stats::rbinom(n, size = 1, prob = 0.5),
    stringsAsFactors = FALSE
  )
  long <- subject_df[rep(seq_len(n), each = time_points), , drop = FALSE]
  time_index <- rep(seq_len(time_points), times = n)
  long$t_scaled <- (time_index - 1) / max(1, time_points - 1)
  long$t_pair_scaled <- pmax(0, (time_index - 1) / max(1, time_points - 2))
  long$time_effect <- jss_joint_time_shape(long$t_scaled, case$time_shape[[1L]])
  long$t_pair_effect <- jss_joint_time_shape(long$t_pair_scaled, case$time_shape[[1L]])
  rownames(long) <- NULL
  long
}

jss_joint_linear_predictor <- function(data, coefs) {
  eta <- rep(unname(coefs[["intercept"]]), nrow(data))
  for (term in setdiff(names(coefs), "intercept")) {
    eta <- eta + unname(coefs[[term]]) * data[[term]]
  }
  eta
}

jss_joint_simulate_case_data <- function(case, seed, response_seed_offset = 0L, covariates = NULL) {
  margin_dist <- jss_joint_margin_dist(case$family[[1L]])
  truth <- jss_joint_truth_coefficients(case)
  if (is.null(covariates)) {
    covariates <- jss_joint_make_covariates(case, seed + 17L)
  }

  margin_params <- list(
    mu = function(data) {
      jss_joint_linkinv(margin_dist, "mu", jss_joint_linear_predictor(data, truth$mu))
    },
    sigma = function(data) {
      pmax(jss_joint_linkinv(margin_dist, "sigma", jss_joint_linear_predictor(data, truth$sigma)), .Machine$double.eps)
    }
  )

  copula_link <- gamlss.longitudinal:::get_copula_dist(case$copula[[1L]])$copula_link
  copula_params <- list(
    theta = function(edge_data) {
      copula_link$theta.linkinv(jss_joint_linear_predictor(edge_data, truth$theta))
    }
  )

  dat <- gamlss.longitudinal::simulate_longitudinal_dataset(
    n = case$n[[1L]],
    times = seq_len(case$time_points[[1L]]),
    margin_dist = margin_dist,
    copula_dist = case$copula[[1L]],
    margin_params = margin_params,
    copula_params = copula_params,
    covariates = covariates,
    seed = seed + response_seed_offset,
    subject_var = "id",
    time_var = "time_index",
    response_var = "response",
    include_truth = TRUE,
    u_bounds = NULL
  )

  edge_left <- dat$time_index < max(dat$time_index)
  edge_right <- dat$time_index > min(dat$time_index)
  true_theta <- rep(NA_real_, nrow(dat))
  true_zeta <- rep(0, nrow(dat))
  true_theta[edge_left] <- dat$true_theta[edge_right]
  true_zeta[edge_left] <- dat$true_zeta[edge_right]
  dat$true_theta <- true_theta
  dat$true_zeta <- true_zeta
  dat$true_tau <- NA_real_
  dat$true_tau[edge_left] <- gamlss.longitudinal:::.copula_par_to_tau(
    case$copula[[1L]],
    dat$true_theta[edge_left],
    dat$true_zeta[edge_left]
  )
  dat
}

jss_joint_fit_model <- function(dat, case, include_dlcopdpar, cfg) {
  margin_dist <- jss_joint_margin_dist(case$family[[1L]])
  fit <- NULL
  warn <- character()
  err <- NA_character_
  elapsed <- system.time({
    capture.output({
      fit <- withCallingHandlers(
        tryCatch(
          gamlss.longitudinal::gamlss_longitudinal(
            dataset = dat,
            margin_dist = margin_dist,
            copula_dist = case$copula[[1L]],
            time_var = "time_index",
            subject_var = "id",
            mu.formula = response ~ x_cont + x_bin + time_effect,
            sigma.formula = ~ x_cont + x_bin + time_effect,
            theta.formula = ~ x_cont + x_bin + t_pair_effect,
            include_dlcopdpar = include_dlcopdpar,
            method = cfg$method %||% "RS",
            compute_vcov = FALSE,
            warm_start_joint = isTRUE(cfg$warm_start_joint),
            warm_start_joint_iter = cfg$warm_start_joint_iter %||% 5L,
            discrete_score_method = cfg$discrete_score_method %||% "analytical",
            cg_max_delta = cfg$cg_max_delta %||% 0.5,
            cg_line_search = cfg$cg_line_search %||% "best",
            cg_gradient_method = cfg$cg_gradient_method %||% "forward",
            cg_hessian_method = cfg$cg_hessian_method %||% "auto",
            cg_raw_loglik_drop_tol = cfg$cg_raw_loglik_drop_tol %||% 10,
            cg_max_line_search_evals = cfg$cg_max_line_search_evals %||% 60L,
            outer_stop_crit = cfg$outer_stop_crit,
            inner_stop_crit = cfg$inner_stop_crit,
            max_outer_iter = cfg$max_outer_iter,
            max_inner_iter = cfg$max_inner_iter,
            max_elapsed_sec = cfg$max_elapsed_sec,
            verbose = 0,
            plot_results = FALSE
          ),
          error = function(e) {
            err <<- conditionMessage(e)
            NULL
          }
        ),
        warning = function(w) {
          warn <<- c(warn, conditionMessage(w))
          invokeRestart("muffleWarning")
        }
      )
    })
  })
  list(fit = fit, elapsed = as.numeric(elapsed[["elapsed"]]), warnings = unique(warn), error = err)
}

jss_joint_extract_eta <- function(fit) {
  copula_link <- gamlss.longitudinal:::get_copula_dist(fit$copula_dist)$copula_link
  gamlss.longitudinal:::calc_eta(fit$par, fit$model_matrix, fit$margin_dist, copula_link, fit$par_s)$eta_inv
}

jss_joint_named_value <- function(x, candidates) {
  if (is.null(x) || length(x) == 0L) return(NA_real_)
  hit <- candidates[candidates %in% names(x)]
  if (length(hit) == 0L) return(NA_real_)
  as.numeric(x[[hit[[1L]]]])
}

jss_joint_fixed_effect_errors <- function(fit, case) {
  truth <- jss_joint_truth_coefficients(case)
  errors <- numeric()
  for (parameter in names(truth)) {
    for (term in names(truth[[parameter]])) {
      estimate <- jss_joint_named_value(
        fit$par,
        c(paste0(parameter, ".", term), paste0(parameter, ".(Intercept)"), paste0(parameter, ".intercept"))
      )
      errors <- c(errors, estimate - unname(truth[[parameter]][[term]]))
    }
  }
  errors
}

jss_joint_variogram_metric_names <- function(power) {
  suffix <- paste0("p", gsub("\\.", "", format(power, trim = TRUE, scientific = FALSE)))
  c(
    score = paste0("score_", suffix),
    mean_subject_score = paste0("mean_subject_score_", suffix),
    n_subjects_scored = paste0("n_subjects_scored_", suffix)
  )
}

jss_joint_variogram_empty <- function(nsim, powers) {
  out <- list(nsim = as.integer(nsim))
  for (power in powers) {
    nms <- jss_joint_variogram_metric_names(power)
    for (nm in nms) {
      out[[nm]] <- NA_real_
    }
  }
  out
}

jss_joint_variogram_score <- function(fit, data, nsim = 100L, powers = c(0.5, 2), seed = NULL) {
  powers <- sort(unique(as.numeric(powers)))
  powers <- powers[is.finite(powers) & powers > 0]
  if (length(powers) == 0L) {
    powers <- 2
  }
  empty <- list(
    nsim = as.integer(nsim)
  )
  empty <- jss_joint_variogram_empty(nsim, powers)
  if (!inherits(fit, "gamlss.longitudinal") || nrow(data) == 0L) {
    return(empty)
  }
  nsim <- as.integer(nsim)
  if (!is.finite(nsim) || nsim < 1L) {
    return(empty)
  }
  if (!requireNamespace("scoringRules", quietly = TRUE)) {
    return(empty)
  }
  sim <- tryCatch(
    stats::simulate(fit, nsim = nsim, seed = seed, newdata = data),
    error = function(e) NULL
  )
  if (is.null(sim) || nrow(sim) != nrow(data)) {
    return(empty)
  }
  out <- list(nsim = nsim)
  for (power in powers) {
    nms <- jss_joint_variogram_metric_names(power)
    subject_scores <- numeric()
    for (id in unique(data$id)) {
      idx <- which(data$id == id)
      idx <- idx[order(data$time_index[idx])]
      y <- as.numeric(data$response[idx])
      dat <- as.matrix(sim[idx, , drop = FALSE])
      ok_rows <- is.finite(y) & rowSums(is.finite(dat)) == ncol(dat)
      if (sum(ok_rows) < 2L) {
        next
      }
      y <- y[ok_rows]
      dat <- dat[ok_rows, , drop = FALSE]
      score <- tryCatch(
        scoringRules::vs_sample(y = y, dat = dat, p = power),
        error = function(e) NA_real_
      )
      if (is.finite(score)) {
        subject_scores <- c(subject_scores, score)
      }
    }
    out[[nms[["score"]]]] <- if (length(subject_scores)) mean(subject_scores, na.rm = TRUE) else NA_real_
    out[[nms[["mean_subject_score"]]]] <- out[[nms[["score"]]]]
    out[[nms[["n_subjects_scored"]]]] <- length(subject_scores)
  }
  out
}

jss_joint_base_result_row <- function(case, rep_idx, method, attempt) {
  data.frame(
    case_id = case$case_id,
    hypothesis_role = case$hypothesis_role,
    joint_review_rep = rep_idx,
    family = case$family,
    copula = case$copula,
    design = "seven_case",
    n = case$n,
    time_points = case$time_points,
    n_subject = case$n,
    n_time = case$time_points,
    total_observations = case$total_observations,
    mu_strength = case$mu_strength,
    sigma_strength = case$sigma_strength,
    theta_strength = case$theta_strength,
    time_shape = case$time_shape,
    dependence = "case_specific",
    missingness = "none",
    start_mode = "default",
    method = method,
    elapsed_sec = attempt$elapsed,
    success = FALSE,
    converged = FALSE,
    outer_iterations = NA_integer_,
    outer_log_lik_change = NA_real_,
    outer_stop_crit = NA_real_,
    hit_outer_limit = NA,
    train_marginal_loglik = NA_real_,
    train_copula_loglik = NA_real_,
    train_joint_loglik = NA_real_,
    marginal_loglik = NA_real_,
    copula_loglik = NA_real_,
    joint_loglik = NA_real_,
    train_rmse_mu = NA_real_,
    train_rmse_sigma = NA_real_,
    train_rmse_theta = NA_real_,
    train_rmse_tau = NA_real_,
    test_marginal_log_score = NA_real_,
    test_copula_log_score = NA_real_,
    test_joint_log_score = NA_real_,
    test_log_score_per_obs = NA_real_,
    heldout_variogram_score_p05 = NA_real_,
    heldout_variogram_score_p2 = NA_real_,
    heldout_variogram_nsim = NA_integer_,
    benchmark_mean_rmse = NA_real_,
    benchmark_neg_log_score = NA_real_,
    benchmark_variogram_score_p05 = NA_real_,
    benchmark_variogram_score_p2 = NA_real_,
    benchmark_theta_time_abs_error = NA_real_,
    max_abs_param_error = NA_real_,
    max_rel_param_error = NA_real_,
    warnings = paste(attempt$warnings, collapse = " | "),
    error = attempt$error,
    failure_type = if (is.na(attempt$error)) "fit_unavailable" else "error",
    stringsAsFactors = FALSE
  )
}

jss_joint_score_fit <- function(fit, train_dat, test_dat, case, rep_idx, method, attempt, cfg) {
  row <- jss_joint_base_result_row(case, rep_idx, method, attempt)
  if (!inherits(fit, "gamlss.longitudinal")) {
    return(row)
  }

  eta <- jss_joint_extract_eta(fit)
  lik_test <- gamlss.longitudinal:::calc_likelihood_minimal(
    eta,
    mm = fit$model_matrix$x,
    margin_dist = fit$margin_dist,
    copula_dist = fit$copula_dist,
    calc_d2 = FALSE,
    response = test_dat$response,
    response_margin = test_dat$time_index,
    response_subject = test_dat$id
  )
  edge_rows <- train_dat$time_index < max(train_dat$time_index)
  theta_hat <- eta$theta
  if (length(theta_hat) == nrow(train_dat)) {
    theta_hat <- theta_hat[edge_rows]
  }
  zeta_hat <- if (!is.null(eta$zeta)) eta$zeta else rep(0, length(theta_hat))
  if (length(zeta_hat) == nrow(train_dat)) {
    zeta_hat <- zeta_hat[edge_rows]
  }
  tau_hat <- gamlss.longitudinal:::.copula_par_to_tau(case$copula[[1L]], theta_hat, zeta_hat)
  fixed_errors <- jss_joint_fixed_effect_errors(fit, case)

  row$success <- all(is.finite(as.numeric(fit$calc_lik_out_end$log_lik)))
  row$converged <- isTRUE(fit$convergence$converged)
  row$outer_iterations <- fit$convergence$outer_iterations %||% NA_integer_
  row$outer_log_lik_change <- fit$convergence$outer_log_lik_change %||% NA_real_
  row$outer_stop_crit <- fit$convergence$outer_stop_crit %||% NA_real_
  row$hit_outer_limit <- isTRUE(fit$convergence$hit_outer_limit)
  row$train_marginal_loglik <- as.numeric(fit$calc_lik_out_end$log_lik["marginal"])
  row$train_copula_loglik <- as.numeric(fit$calc_lik_out_end$log_lik["copula"])
  row$train_joint_loglik <- as.numeric(fit$calc_lik_out_end$log_lik["joint"])
  row$marginal_loglik <- row$train_marginal_loglik
  row$copula_loglik <- row$train_copula_loglik
  row$joint_loglik <- row$train_joint_loglik
  row$train_rmse_mu <- jss_joint_rmse(eta$mu - train_dat$true_mu)
  row$train_rmse_sigma <- jss_joint_rmse(eta$sigma - train_dat$true_sigma)
  row$train_rmse_theta <- jss_joint_rmse(theta_hat - train_dat$true_theta[edge_rows])
  row$train_rmse_tau <- jss_joint_rmse(tau_hat - train_dat$true_tau[edge_rows])
  row$test_marginal_log_score <- as.numeric(lik_test$log_lik["marginal"])
  row$test_copula_log_score <- as.numeric(lik_test$log_lik["copula"])
  row$test_joint_log_score <- as.numeric(lik_test$log_lik["joint"])
  row$test_log_score_per_obs <- row$test_joint_log_score / nrow(test_dat)
  variogram <- jss_joint_variogram_score(
    fit,
    test_dat,
    nsim = cfg$variogram_nsim %||% 100L,
    powers = c(0.5, 2),
    seed = 900000L + as.integer(sub("^JVS", "", case$case_id)) * 1000L + rep_idx * 10L
  )
  row$heldout_variogram_score_p05 <- variogram$score_p05
  row$heldout_variogram_score_p2 <- variogram$score_p2
  row$heldout_variogram_nsim <- variogram$nsim
  row$benchmark_mean_rmse <- row$train_rmse_mu
  row$benchmark_neg_log_score <- -row$test_log_score_per_obs
  row$benchmark_variogram_score_p05 <- row$heldout_variogram_score_p05
  row$benchmark_variogram_score_p2 <- row$heldout_variogram_score_p2
  row$benchmark_theta_time_abs_error <- row$train_rmse_tau
  row$max_abs_param_error <- if (length(fixed_errors)) max(abs(fixed_errors), na.rm = TRUE) else NA_real_
  row$max_rel_param_error <- if (length(fixed_errors)) {
    truth_vals <- unlist(jss_joint_truth_coefficients(case), use.names = FALSE)
    max(abs(fixed_errors) / pmax(abs(truth_vals), 1e-8), na.rm = TRUE)
  } else {
    NA_real_
  }
  row$failure_type <- if (isTRUE(row$success)) "none" else "invalid_loglik"
  row
}

jss_joint_run_case_rep <- function(case, rep_idx, settings, cfg) {
  seed <- settings$seed + 3000L + as.integer(sub("^JVS", "", case$case_id)) * 100000L + rep_idx * 100L
  covariates <- jss_joint_make_covariates(case, seed + 17L)
  train_dat <- jss_joint_simulate_case_data(case, seed, response_seed_offset = 0L, covariates = covariates)
  test_dat <- jss_joint_simulate_case_data(case, seed, response_seed_offset = 500000L, covariates = covariates)

  sep <- jss_joint_fit_model(train_dat, case, include_dlcopdpar = FALSE, cfg = cfg)
  joint <- jss_joint_fit_model(train_dat, case, include_dlcopdpar = TRUE, cfg = cfg)

  jss_joint_bind_rows(list(
    jss_joint_score_fit(sep$fit, train_dat, test_dat, case, rep_idx, "rs_separate", sep, cfg),
    jss_joint_score_fit(joint$fit, train_dat, test_dat, case, rep_idx, "rs_joint", joint, cfg)
  ))
}

jss_joint_run_simulation <- function(settings, cases) {
  cfg <- jss_joint_simulation_settings(settings)
  rows <- list()
  for (case_idx in seq_len(nrow(cases))) {
    case <- cases[case_idx, , drop = FALSE]
    for (rep_idx in seq_len(cfg$reps)) {
      rows[[length(rows) + 1L]] <- jss_joint_run_case_rep(case, rep_idx, settings, cfg)
    }
  }
  jss_joint_bind_rows(rows)
}

jss_joint_delta_table <- function(results) {
  if (nrow(results) == 0L) {
    return(data.frame())
  }
  group_cols <- intersect(
    c(
      "case_id", "hypothesis_role", "joint_review_rep", "family", "copula", "design",
      "n", "time_points", "n_subject", "n_time", "total_observations",
      "mu_strength", "sigma_strength", "theta_strength", "time_shape",
      "dependence", "missingness", "start_mode"
    ),
    names(results)
  )
  metrics <- intersect(
    c(
      "joint_loglik", "copula_loglik", "marginal_loglik", "elapsed_sec",
      "train_joint_loglik", "train_copula_loglik", "train_marginal_loglik",
      "train_rmse_mu", "train_rmse_sigma", "train_rmse_theta", "train_rmse_tau",
      "test_joint_log_score", "test_log_score_per_obs",
      "test_marginal_log_score", "test_copula_log_score",
      "heldout_variogram_score_p05", "heldout_variogram_score_p2",
      "benchmark_mean_rmse", "benchmark_neg_log_score",
      "benchmark_variogram_score_p05", "benchmark_variogram_score_p2",
      "benchmark_q90_mae", "benchmark_upper_tail_error_90",
      "benchmark_interval_coverage_95", "benchmark_interval_width_95",
      "benchmark_pit_mean_abs_error", "benchmark_tail_error_lower_05",
      "benchmark_tail_error_upper_05", "benchmark_theta_time_abs_error",
      "smooth_eta_rmse", "smooth_eta_max_abs_error"
    ),
    names(results)
  )
  key <- interaction(results[group_cols], drop = TRUE, lex.order = TRUE)
  rows <- lapply(unique(key), function(k) {
    x <- results[key == k, , drop = FALSE]
    joint <- x[x$method == "rs_joint", , drop = FALSE][1L, , drop = FALSE]
    separate <- x[x$method == "rs_separate", , drop = FALSE][1L, , drop = FALSE]
    if (nrow(joint) == 0L || nrow(separate) == 0L) {
      return(NULL)
    }
    row <- joint[, group_cols, drop = FALSE]
    row$rs_joint_success <- joint$success
    row$rs_separate_success <- separate$success
    for (metric in metrics) {
      row[[paste0("delta_", metric)]] <-
        jss_joint_numeric(joint[[metric]]) - jss_joint_numeric(separate[[metric]])
    }
    row
  })
  jss_joint_bind_rows(rows)
}

jss_joint_summary_table <- function(results, deltas, candidates) {
  if (nrow(results) == 0L) {
    return(data.frame())
  }
  split_key <- interaction(
    results$case_id,
    results$hypothesis_role,
    results$method,
    drop = TRUE,
    lex.order = TRUE
  )
  rows <- lapply(split(seq_len(nrow(results)), split_key), function(idx) {
    x <- results[idx, , drop = FALSE]
    data.frame(
      profile = jss_joint_or(x$profile, NA_character_),
      case_id = x$case_id[[1L]],
      hypothesis_role = x$hypothesis_role[[1L]],
      method = x$method[[1L]],
      family = x$family[[1L]],
      copula = x$copula[[1L]],
      design = x$design[[1L]],
      n_subject = x$n_subject[[1L]],
      n_time = x$n_time[[1L]],
      total_observations = x$total_observations[[1L]],
      mu_strength = x$mu_strength[[1L]],
      sigma_strength = x$sigma_strength[[1L]],
      theta_strength = x$theta_strength[[1L]],
      time_shape = x$time_shape[[1L]],
      n = nrow(x),
      success_rate = mean(x$success %in% TRUE, na.rm = TRUE),
      convergence_rate = mean(x$converged %in% TRUE, na.rm = TRUE),
      median_elapsed_sec = stats::median(jss_joint_numeric(x$elapsed_sec), na.rm = TRUE),
      mean_train_joint_loglik = mean(jss_joint_numeric(x$train_joint_loglik), na.rm = TRUE),
      mean_test_log_score_per_obs = mean(jss_joint_numeric(x$test_log_score_per_obs), na.rm = TRUE),
      mean_heldout_variogram_score_p05 = mean(jss_joint_numeric(x$heldout_variogram_score_p05), na.rm = TRUE),
      mean_heldout_variogram_score_p2 = mean(jss_joint_numeric(x$heldout_variogram_score_p2), na.rm = TRUE),
      mean_train_rmse_theta = mean(jss_joint_numeric(x$train_rmse_theta), na.rm = TRUE),
      mean_train_rmse_tau = mean(jss_joint_numeric(x$train_rmse_tau), na.rm = TRUE),
      mean_benchmark_neg_log_score = mean(jss_joint_numeric(x$benchmark_neg_log_score), na.rm = TRUE),
      mean_benchmark_variogram_score_p05 = mean(jss_joint_numeric(x$benchmark_variogram_score_p05), na.rm = TRUE),
      mean_benchmark_variogram_score_p2 = mean(jss_joint_numeric(x$benchmark_variogram_score_p2), na.rm = TRUE),
      mean_benchmark_mean_rmse = mean(jss_joint_numeric(x$benchmark_mean_rmse), na.rm = TRUE),
      selected_review_cases = nrow(candidates),
      stringsAsFactors = FALSE
    )
  })
  summary <- jss_joint_bind_rows(rows)

  if (nrow(deltas) > 0L) {
    delta_key <- interaction(deltas$case_id, deltas$hypothesis_role, drop = TRUE, lex.order = TRUE)
    delta_summary <- jss_joint_bind_rows(lapply(split(seq_len(nrow(deltas)), delta_key), function(idx) {
      x <- deltas[idx, , drop = FALSE]
      data.frame(
        profile = jss_joint_or(unique(results$profile), NA_character_),
        case_id = x$case_id[[1L]],
        hypothesis_role = x$hypothesis_role[[1L]],
        method = "rs_joint_minus_rs_separate",
        family = x$family[[1L]],
        copula = x$copula[[1L]],
        design = x$design[[1L]],
        n_subject = x$n_subject[[1L]],
        n_time = x$n_time[[1L]],
        total_observations = x$total_observations[[1L]],
        mu_strength = x$mu_strength[[1L]],
        sigma_strength = x$sigma_strength[[1L]],
        theta_strength = x$theta_strength[[1L]],
        time_shape = x$time_shape[[1L]],
        n = nrow(x),
        success_rate = mean(x$rs_joint_success %in% TRUE & x$rs_separate_success %in% TRUE, na.rm = TRUE),
        convergence_rate = NA_real_,
        median_elapsed_sec = stats::median(jss_joint_numeric(x$delta_elapsed_sec), na.rm = TRUE),
        mean_train_joint_loglik = mean(jss_joint_numeric(x$delta_train_joint_loglik), na.rm = TRUE),
        mean_test_log_score_per_obs = mean(jss_joint_numeric(x$delta_test_log_score_per_obs), na.rm = TRUE),
        mean_heldout_variogram_score_p05 = mean(jss_joint_numeric(x$delta_heldout_variogram_score_p05), na.rm = TRUE),
        mean_heldout_variogram_score_p2 = mean(jss_joint_numeric(x$delta_heldout_variogram_score_p2), na.rm = TRUE),
        mean_train_rmse_theta = mean(jss_joint_numeric(x$delta_train_rmse_theta), na.rm = TRUE),
        mean_train_rmse_tau = mean(jss_joint_numeric(x$delta_train_rmse_tau), na.rm = TRUE),
        mean_benchmark_neg_log_score = mean(jss_joint_numeric(x$delta_benchmark_neg_log_score), na.rm = TRUE),
        mean_benchmark_variogram_score_p05 = mean(jss_joint_numeric(x$delta_benchmark_variogram_score_p05), na.rm = TRUE),
        mean_benchmark_variogram_score_p2 = mean(jss_joint_numeric(x$delta_benchmark_variogram_score_p2), na.rm = TRUE),
        mean_benchmark_mean_rmse = mean(jss_joint_numeric(x$delta_benchmark_mean_rmse), na.rm = TRUE),
        selected_review_cases = nrow(candidates),
        stringsAsFactors = FALSE
      )
    }))
    summary <- jss_joint_bind_rows(list(summary, delta_summary))
  }

  summary[order(summary$case_id, summary$method), , drop = FALSE]
}

jss_joint_metric_wins <- function(results) {
  metrics <- intersect(
    c(
      "benchmark_mean_rmse", "benchmark_mean_mae", "benchmark_mean_bias",
      "benchmark_rmse", "benchmark_mae", "benchmark_q90_mae",
      "benchmark_neg_log_score", "benchmark_upper_tail_error_90",
      "benchmark_variogram_score_p05", "benchmark_variogram_score_p2",
      "benchmark_theta_time_abs_error", "benchmark_interval_coverage_95",
      "benchmark_interval_width_95", "benchmark_pit_mean_abs_error",
      "benchmark_pit_ks_p_value", "benchmark_tail_error_lower_05",
      "benchmark_tail_error_upper_05", "smooth_eta_rmse",
      "smooth_eta_max_abs_error", "elapsed_sec"
    ),
    names(results)
  )
  if (nrow(results) == 0L || length(metrics) == 0L) {
    return(data.frame())
  }

  summary <- tryCatch(
    gamlss.longitudinal::summarise_benchmark_results(
      results,
      metrics = metrics,
      group_cols = intersect(
        c("case_id", "joint_review_rep", "family", "copula", "design", "n_subject", "n_time", "dependence", "missingness", "start_mode"),
        names(results)
      )
    ),
    error = function(e) NULL
  )
  if (is.null(summary)) {
    return(data.frame())
  }
  summary$summary
}

jss_joint_write_delta_figure <- function(deltas, path) {
  if (nrow(deltas) == 0L) {
    deltas <- data.frame(metric = "no finite deltas", value = 0)
  } else {
    metrics <- intersect(
      c(
        "delta_train_joint_loglik", "delta_test_log_score_per_obs",
        "delta_heldout_variogram_score_p05", "delta_heldout_variogram_score_p2",
        "delta_train_rmse_theta", "delta_train_rmse_tau",
        "delta_elapsed_sec"
      ),
      names(deltas)
    )
    plot_rows <- lapply(metrics, function(metric) {
      data.frame(metric = metric, value = jss_joint_numeric(deltas[[metric]]), stringsAsFactors = FALSE)
    })
    deltas <- jss_joint_bind_rows(plot_rows)
    deltas <- deltas[is.finite(deltas$value), , drop = FALSE]
    if (nrow(deltas) == 0L) {
      deltas <- data.frame(metric = "no finite deltas", value = 0)
    }
  }

  p <- ggplot2::ggplot(deltas, ggplot2::aes(x = metric, y = value, fill = metric)) +
    ggplot2::geom_boxplot(width = 0.65, outlier.alpha = 0.45) +
    ggplot2::coord_flip() +
    ggplot2::labs(x = NULL, y = "Joint - separate", title = "Joint versus separate deltas") +
    ggplot2::theme_minimal(base_size = 10) +
    ggplot2::theme(legend.position = "none")
  ggplot2::ggsave(path, p, width = 8, height = 5, dpi = 320, bg = "white")
  path
}

jss_joint_case_delta_summary <- function(deltas) {
  if (nrow(deltas) == 0L) {
    return(data.frame())
  }
  split_key <- interaction(deltas$case_id, drop = TRUE, lex.order = TRUE)
  jss_joint_bind_rows(lapply(split(seq_len(nrow(deltas)), split_key), function(idx) {
    x <- deltas[idx, , drop = FALSE]
    dep_delta <- pmin(
      jss_joint_numeric(x$delta_train_rmse_theta),
      jss_joint_numeric(x$delta_train_rmse_tau),
      na.rm = TRUE
    )
    dep_delta[!is.finite(dep_delta)] <- NA_real_
    data.frame(
      case_id = x$case_id[[1L]],
      hypothesis_role = x$hypothesis_role[[1L]],
      family = x$family[[1L]],
      n_subject = x$n_subject[[1L]],
      n_time = x$n_time[[1L]],
      total_observations = x$total_observations[[1L]],
      mu_strength = x$mu_strength[[1L]],
      sigma_strength = x$sigma_strength[[1L]],
      theta_strength = x$theta_strength[[1L]],
      reps = nrow(x),
      paired_success_rate = mean(x$rs_joint_success %in% TRUE & x$rs_separate_success %in% TRUE, na.rm = TRUE),
      delta_train_joint_loglik_mean = mean(jss_joint_numeric(x$delta_train_joint_loglik), na.rm = TRUE),
      delta_test_log_score_per_obs_mean = mean(jss_joint_numeric(x$delta_test_log_score_per_obs), na.rm = TRUE),
      delta_heldout_variogram_score_p05_mean = mean(jss_joint_numeric(x$delta_heldout_variogram_score_p05), na.rm = TRUE),
      delta_heldout_variogram_score_p2_mean = mean(jss_joint_numeric(x$delta_heldout_variogram_score_p2), na.rm = TRUE),
      delta_train_rmse_theta_mean = mean(jss_joint_numeric(x$delta_train_rmse_theta), na.rm = TRUE),
      delta_train_rmse_tau_mean = mean(jss_joint_numeric(x$delta_train_rmse_tau), na.rm = TRUE),
      delta_dependence_rmse_mean = mean(dep_delta, na.rm = TRUE),
      stringsAsFactors = FALSE
    )
  }))
}

jss_joint_gate_status <- function(row) {
  if (is.null(row) || nrow(row) == 0L) {
    return(FALSE)
  }
  train_ok <- is.finite(row$delta_train_joint_loglik_mean) && row$delta_train_joint_loglik_mean > 0
  variogram_p05_ok <- is.finite(row$delta_heldout_variogram_score_p05_mean) &&
    row$delta_heldout_variogram_score_p05_mean < 0
  variogram_p2_ok <- is.finite(row$delta_heldout_variogram_score_p2_mean) &&
    row$delta_heldout_variogram_score_p2_mean < 0
  train_ok && variogram_p05_ok && variogram_p2_ok
}

jss_joint_hypothesis_row <- function(case_summary, hypothesis, focal_case, comparator_case,
                                     decision, rationale) {
  focal <- case_summary[case_summary$case_id == focal_case, , drop = FALSE]
  comparator <- case_summary[case_summary$case_id == comparator_case, , drop = FALSE]
  data.frame(
    hypothesis = hypothesis,
    focal_case = focal_case,
    comparator_case = comparator_case,
    decision = decision,
    focal_gate_pass = jss_joint_gate_status(focal),
    comparator_gate_pass = jss_joint_gate_status(comparator),
    focal_delta_train_joint_loglik = jss_joint_or(focal$delta_train_joint_loglik_mean, NA_real_),
    comparator_delta_train_joint_loglik = jss_joint_or(comparator$delta_train_joint_loglik_mean, NA_real_),
    focal_delta_test_log_score_per_obs = jss_joint_or(focal$delta_test_log_score_per_obs_mean, NA_real_),
    comparator_delta_test_log_score_per_obs = jss_joint_or(comparator$delta_test_log_score_per_obs_mean, NA_real_),
    focal_delta_heldout_variogram_score_p05 = jss_joint_or(focal$delta_heldout_variogram_score_p05_mean, NA_real_),
    comparator_delta_heldout_variogram_score_p05 = jss_joint_or(comparator$delta_heldout_variogram_score_p05_mean, NA_real_),
    focal_delta_heldout_variogram_score_p2 = jss_joint_or(focal$delta_heldout_variogram_score_p2_mean, NA_real_),
    comparator_delta_heldout_variogram_score_p2 = jss_joint_or(comparator$delta_heldout_variogram_score_p2_mean, NA_real_),
    focal_delta_dependence_rmse = jss_joint_or(focal$delta_dependence_rmse_mean, NA_real_),
    comparator_delta_dependence_rmse = jss_joint_or(comparator$delta_dependence_rmse_mean, NA_real_),
    rationale = rationale,
    stringsAsFactors = FALSE
  )
}

jss_joint_compare_cases <- function(case_summary, focal_case, comparator_case,
                                    require_variogram_gain = TRUE) {
  focal <- case_summary[case_summary$case_id == focal_case, , drop = FALSE]
  comparator <- case_summary[case_summary$case_id == comparator_case, , drop = FALSE]
  if (nrow(focal) == 0L || nrow(comparator) == 0L) {
    return("inconclusive")
  }
  if (!jss_joint_gate_status(focal)) {
    return("disconfirmed")
  }
  variogram_ok <- !require_variogram_gain ||
    (
      focal$delta_heldout_variogram_score_p05_mean < comparator$delta_heldout_variogram_score_p05_mean &&
        focal$delta_heldout_variogram_score_p2_mean < comparator$delta_heldout_variogram_score_p2_mean
    )
  if (isTRUE(variogram_ok)) "confirmed" else "inconclusive"
}

jss_joint_hypothesis_summary <- function(deltas) {
  case_summary <- jss_joint_case_delta_summary(deltas)
  if (nrow(case_summary) == 0L) {
    return(data.frame())
  }

  rows <- list()
  rows[[length(rows) + 1L]] <- jss_joint_hypothesis_row(
    case_summary,
    "More time points improve joint-vs-separate performance at similar total rows",
    "JVS03",
    "JVS04",
    jss_joint_compare_cases(case_summary, "JVS03", "JVS04", require_variogram_gain = TRUE),
    "JVS03 has T=75 and n=40; JVS04 has T=20 and n=150, keeping total observations at 3000."
  )
  rows[[length(rows) + 1L]] <- jss_joint_hypothesis_row(
    case_summary,
    "Stronger dependence improves joint-vs-separate performance",
    "JVS03",
    "JVS02",
    jss_joint_compare_cases(case_summary, "JVS03", "JVS02", require_variogram_gain = TRUE),
    "JVS03 increases theta strength from 1.5 to 3.5 with all other design settings fixed."
  )
  rows[[length(rows) + 1L]] <- jss_joint_hypothesis_row(
    case_summary,
    "Theta covariates improve joint-vs-separate performance",
    "JVS03",
    "JVS01",
    jss_joint_compare_cases(case_summary, "JVS03", "JVS01", require_variogram_gain = TRUE),
    "JVS01 removes theta covariate signal; JVS03 restores high theta signal."
  )
  rows[[length(rows) + 1L]] <- jss_joint_hypothesis_row(
    case_summary,
    "Shared shape and theta covariates improve performance",
    "JVS03",
    "JVS05",
    jss_joint_compare_cases(case_summary, "JVS03", "JVS05", require_variogram_gain = TRUE),
    "JVS05 weakens sigma covariate signal while retaining high theta signal."
  )
  sample_decision <- {
    focal <- case_summary[case_summary$case_id == "JVS06", , drop = FALSE]
    if (nrow(focal) == 0L) {
      "inconclusive"
    } else if (jss_joint_gate_status(focal)) {
      "confirmed"
    } else {
      "disconfirmed"
    }
  }
  rows[[length(rows) + 1L]] <- jss_joint_hypothesis_row(
    case_summary,
    "Larger sample size stabilizes or improves joint fits",
    "JVS06",
    "JVS03",
    sample_decision,
    "JVS06 raises n from 40 to 100 at T=75; interpret effect-size shrinkage separately from gate success."
  )
  discrete_decision <- {
    focal <- case_summary[case_summary$case_id == "JVS07", , drop = FALSE]
    comparator <- case_summary[case_summary$case_id == "JVS03", , drop = FALSE]
    if (nrow(focal) == 0L || nrow(comparator) == 0L || focal$paired_success_rate < 0.8) {
      "inconclusive"
    } else {
      gate_diff <- !identical(jss_joint_gate_status(focal), jss_joint_gate_status(comparator))
      variogram_diff <- abs(focal$delta_heldout_variogram_score_p2_mean - comparator$delta_heldout_variogram_score_p2_mean) > 0.005 ||
        abs(focal$delta_heldout_variogram_score_p05_mean - comparator$delta_heldout_variogram_score_p05_mean) > 0.005
      if (gate_diff || variogram_diff) "confirmed" else "disconfirmed"
    }
  }
  rows[[length(rows) + 1L]] <- jss_joint_hypothesis_row(
    case_summary,
    "Discrete margins show different joint-vs-separate behaviour",
    "JVS07",
    "JVS03",
    discrete_decision,
    "JVS07 changes the margin from NO to NBI while retaining the high-signal case design."
  )

  jss_joint_bind_rows(rows)
}

jss_joint_write_metric_dashboard <- function(metric_wins, path) {
  if (nrow(metric_wins) == 0L) {
    plot_data <- data.frame(method = "no finite metrics", metric = "missing", win_or_tie_rate = 0)
  } else {
    plot_data <- metric_wins
    plot_data$metric <- as.character(plot_data$metric)
    plot_data$method <- as.character(plot_data$method)
    plot_data$win_or_tie_rate <- jss_joint_numeric(plot_data$win_or_tie_rate)
    plot_data <- plot_data[is.finite(plot_data$win_or_tie_rate), , drop = FALSE]
    if (nrow(plot_data) == 0L) {
      plot_data <- data.frame(method = "no finite metrics", metric = "missing", win_or_tie_rate = 0)
    }
  }

  p <- ggplot2::ggplot(plot_data, ggplot2::aes(x = metric, y = win_or_tie_rate, fill = method)) +
    ggplot2::geom_col(position = "dodge", width = 0.7) +
    ggplot2::coord_flip(ylim = c(0, 1)) +
    ggplot2::labs(x = NULL, y = "Win or tie rate", title = "Benchmark metric wins and ties") +
    ggplot2::theme_minimal(base_size = 10) +
    ggplot2::theme(legend.position = "bottom")
  ggplot2::ggsave(path, p, width = 9, height = 6, dpi = 320, bg = "white")
  path
}

jss_run_03_joint_vs_separate <- function(settings) {
  paths <- jss_joint_output_paths(settings)
  candidates <- jss_joint_case_definitions()
  utils::write.csv(candidates, paths$candidate_selection, row.names = FALSE)

  results <- jss_joint_run_simulation(settings, candidates)
  results$profile <- settings$profile
  utils::write.csv(results, paths$results, row.names = FALSE)

  deltas <- jss_joint_delta_table(results)
  summary <- jss_joint_summary_table(results, deltas, candidates)
  metric_wins <- jss_joint_metric_wins(results)
  hypothesis_summary <- jss_joint_hypothesis_summary(deltas)
  utils::write.csv(summary, paths$summary, row.names = FALSE)
  utils::write.csv(metric_wins, paths$metric_wins, row.names = FALSE)
  utils::write.csv(hypothesis_summary, paths$hypothesis_summary, row.names = FALSE)

  jss_joint_write_delta_figure(deltas, paths$deltas_figure)
  jss_joint_write_metric_dashboard(metric_wins, paths$metric_dashboard)

  list(
    module_id = jss_joint_module_id(),
    title = "Joint versus separate optimisation",
    status = "current",
    data = c(paths$candidate_selection, paths$results),
    tables = c(paths$summary, paths$metric_wins, paths$hypothesis_summary),
    figures = c(paths$deltas_figure, paths$metric_dashboard),
    notes = paste(
      "Seven fixed case studies compare rs_separate and rs_joint against",
      "pre-specified likelihood, held-out score, and dependence recovery gates."
    )
  )
}
