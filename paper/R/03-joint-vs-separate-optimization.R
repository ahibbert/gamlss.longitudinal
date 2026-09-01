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

jss_joint_env_flag <- function(name, default = TRUE) {
  raw <- tolower(trimws(Sys.getenv(name, unset = if (isTRUE(default)) "1" else "0")))
  raw %in% c("1", "true", "yes", "on")
}

jss_joint_mc_precision_registry <- function(profile = "full") {
  confidence <- 0.95
  precision_raw <- Sys.getenv("GAMLSS_LONGITUDINAL_JSS_JVS_MC_HALF_WIDTH", unset = "")
  target_half_width <- if (nzchar(precision_raw)) suppressWarnings(as.numeric(precision_raw)) else 0.10
  if (!is.finite(target_half_width) || target_half_width <= 0) {
    stop("Optimizer Monte Carlo half-width must be finite and strictly positive.", call. = FALSE)
  }
  worst_case_probability <- 0.50
  z <- stats::qnorm(1 - (1 - confidence) / 2)
  required_retained_pairs <- ceiling(
    z^2 * worst_case_probability * (1 - worst_case_probability) /
      target_half_width^2
  )
  registered_required <- if (identical(profile, "smoke")) 1L else as.integer(required_retained_pairs)
  data.frame(
    precision_estimand = "paired binary win/sign probability",
    confidence_level = confidence,
    target_half_width = target_half_width,
    planning_probability = worst_case_probability,
    required_retained_pairs = registered_required,
    initial_attempts = registered_required,
    top_up_batch = if (identical(profile, "smoke")) 1L else jss_joint_env_int(
      "GAMLSS_LONGITUDINAL_JSS_JVS_TOP_UP_BATCH", 10L
    ),
    max_attempts = if (identical(profile, "smoke")) 1L else jss_joint_env_int(
      "GAMLSS_LONGITUDINAL_JSS_JVS_MAX_ATTEMPTS", 2L * as.integer(required_retained_pairs),
      minimum = as.integer(required_retained_pairs)
    ),
    registered_scenario_cells = 4L,
    registered_method_cells = 8L,
    attempt_count_status = "dynamic: initial registered batch plus deficient-cell top-ups only",
    production_requires_precision = !identical(profile, "smoke"),
    eligibility_rule = "both RS fits returned and converged; begin at the retained-pair target and top up deficient cells in registered batches until achieved width or hard cap",
    formula = "ceiling(z_(1-alpha/2)^2 * p * (1-p) / half_width^2)",
    stringsAsFactors = FALSE
  )
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
    precision_registry = file.path(
      settings$data_dir,
      paste0(module_id, "-mc-precision.csv")
    ),
    checkpoint_status = file.path(
      settings$data_dir,
      paste0(module_id, "-checkpoint-status.csv")
    ),
    figure_registry = file.path(
      settings$data_dir,
      paste0(module_id, "-figure-registry.csv")
    ),
    checkpoint_payload_archive = file.path(
      settings$data_dir, paste0(module_id, "-checkpoint-payloads.rds")
    ),
    checkpoint_content_manifest = file.path(
      settings$data_dir, paste0(module_id, "-checkpoint-content-manifest.csv")
    ),
    failure_summary = file.path(
      settings$tables_dir,
      paste0(module_id, "-failures.csv")
    ),
    difference_uncertainty = file.path(
      settings$tables_dir,
      paste0(module_id, "-difference-uncertainty.csv")
    ),
    hypothesis_evidence = file.path(
      settings$tables_dir,
      paste0(module_id, "-hypothesis-evidence.csv")
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

jss_joint_mean_finite <- function(x) {
  x <- jss_joint_numeric(x)
  x <- x[is.finite(x)]
  if (length(x)) mean(x) else NA_real_
}

jss_joint_median_finite <- function(x) {
  x <- jss_joint_numeric(x)
  x <- x[is.finite(x)]
  if (length(x)) stats::median(x) else NA_real_
}

jss_joint_eligible <- function(x, success = "success", converged = "converged") {
  if (!is.data.frame(x) || !all(c(success, converged) %in% names(x))) {
    return(rep(FALSE, if (is.data.frame(x)) nrow(x) else 0L))
  }
  x[[success]] %in% TRUE & x[[converged]] %in% TRUE
}

jss_joint_pair_eligible <- function(x) {
  jss_joint_eligible(x, "rs_joint_success", "rs_joint_converged") &
    jss_joint_eligible(x, "rs_separate_success", "rs_separate_converged")
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
  cases <- data.frame(
    candidate_id = sprintf("jvs-%02d", 1:4),
    case_id = sprintf("JVS%02d", 1:4),
    base_case_id = "JVS01",
    contrast_factor = c("base", "theta_strength", "sigma_strength", "family"),
    contrast_label = c("base", "dependence strength", "scale signal", "margin family"),
    contrast_level = c("reference", "strong", "weak", "count"),
    hypothesis_role = c(
      "base_design",
      "stronger_dependence",
      "weaker_scale_signal",
      "count_margin"
    ),
    family = c("NO", "NO", "NO", "NBI"),
    copula = "N",
    n = 120L,
    time_points = 8L,
    total_observations = 120L * 8L,
    mu_strength = 1,
    sigma_strength = c(1.5, 1.5, 0.5, 1.5),
    theta_strength = c(1.5, 3.5, 1.5, 1.5),
    time_shape = "sigmoid",
    purpose = c(
      "Reference continuous design for every paired contrast.",
      "Changes only the dependence-strength factor from the reference.",
      "Changes only the scale-signal factor from the reference.",
      "Changes only the margin-family factor from the reference."
    ),
    review_role = c("reference", "contrast", "contrast", "contrast"),
    candidate_score = NA_real_,
    source_name = "registered_paired_one_factor_design",
    source_file = NA_character_,
    source_row = seq_len(4L),
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
  jss_joint_validate_case_definitions(cases)
  cases
}

jss_joint_validate_case_definitions <- function(cases) {
  design_fields <- c(
    "family", "copula", "n", "time_points", "mu_strength",
    "sigma_strength", "theta_strength", "time_shape"
  )
  required <- c("case_id", "base_case_id", "contrast_factor", design_fields)
  missing <- setdiff(required, names(cases))
  if (length(missing)) {
    stop("Optimizer design is missing field(s): ", paste(missing, collapse = ", "), call. = FALSE)
  }
  base <- cases[cases$case_id == unique(cases$base_case_id), , drop = FALSE]
  if (nrow(base) != 1L || !identical(base$contrast_factor[[1L]], "base")) {
    stop("Optimizer design must contain exactly one declared base case.", call. = FALSE)
  }
  variants <- cases[cases$case_id != base$case_id[[1L]], , drop = FALSE]
  for (i in seq_len(nrow(variants))) {
    changed <- design_fields[vapply(design_fields, function(nm) {
      !identical(as.character(variants[[nm]][[i]]), as.character(base[[nm]][[1L]]))
    }, logical(1))]
    if (length(changed) != 1L || !identical(changed, variants$contrast_factor[[i]])) {
      stop(
        "Case ", variants$case_id[[i]], " must change exactly its declared factor; changed: ",
        paste(changed, collapse = ", "),
        call. = FALSE
      )
    }
  }
  invisible(TRUE)
}

jss_joint_simulation_settings <- function(settings) {
  precision <- jss_joint_mc_precision_registry(settings$profile)
  if (settings$profile %in% c("expanded", "paper", "full")) {
    return(list(
      reps = precision$max_attempts[[1L]],
      initial_attempts = precision$initial_attempts[[1L]],
      top_up_batch = precision$top_up_batch[[1L]],
      precision = precision,
      resume = jss_joint_env_flag("GAMLSS_LONGITUDINAL_JSS_JVS_RESUME", TRUE),
      max_outer_iter = jss_joint_env_int("GAMLSS_LONGITUDINAL_JSS_JVS_MAX_OUTER", 250L),
      max_inner_iter = jss_joint_env_int("GAMLSS_LONGITUDINAL_JSS_JVS_MAX_INNER", 100L),
      max_elapsed_sec = jss_joint_env_num("GAMLSS_LONGITUDINAL_JSS_JVS_MAX_ELAPSED", Inf),
      method = "RS",
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
      variogram_nsim = jss_joint_env_int("GAMLSS_LONGITUDINAL_JSS_JVS_VARIOGRAM_NSIM", 100L),
      workers = jss_joint_env_int("GAMLSS_LONGITUDINAL_JSS_JVS_WORKERS", 1L)
    ))
  }

  list(
    reps = precision$max_attempts[[1L]],
    initial_attempts = precision$initial_attempts[[1L]],
    top_up_batch = precision$top_up_batch[[1L]],
    precision = precision,
    resume = jss_joint_env_flag("GAMLSS_LONGITUDINAL_JSS_JVS_RESUME", TRUE),
    max_outer_iter = jss_joint_env_int("GAMLSS_LONGITUDINAL_JSS_JVS_MAX_OUTER", 3L),
    max_inner_iter = jss_joint_env_int("GAMLSS_LONGITUDINAL_JSS_JVS_MAX_INNER", 3L),
    max_elapsed_sec = jss_joint_env_num("GAMLSS_LONGITUDINAL_JSS_JVS_MAX_ELAPSED", 20),
    method = "RS",
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
    variogram_nsim = jss_joint_env_int("GAMLSS_LONGITUDINAL_JSS_JVS_VARIOGRAM_NSIM", 20L),
    workers = jss_joint_env_int("GAMLSS_LONGITUDINAL_JSS_JVS_WORKERS", 1L)
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
            optimizer_control = gamlss.longitudinal::gamlss_longitudinal_control(
              outer_tol = cfg$outer_stop_crit,
              max_outer_iter = cfg$max_outer_iter,
              max_elapsed_sec = cfg$max_elapsed_sec,
              rs = list(
                inner_tol = cfg$inner_stop_crit,
                max_inner_iter = cfg$max_inner_iter,
                warm_start_joint = isTRUE(cfg$warm_start_joint),
                warm_start_joint_iter = cfg$warm_start_joint_iter %||% 5L,
                discrete_score_method = cfg$discrete_score_method %||% "analytical"
              )
            ),
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
    base_case_id = case$base_case_id,
    contrast_factor = case$contrast_factor,
    contrast_label = case$contrast_label,
    contrast_level = case$contrast_level,
    hypothesis_role = case$hypothesis_role,
    joint_review_rep = rep_idx,
    family = case$family,
    copula = case$copula,
    design = "paired_one_factor",
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
    retained = FALSE,
    stop_reason = "fit_error",
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

# Lightweight contract fixture used by checkpoint/PSOCK tests. It deliberately
# exercises the complete public payload schema without fitting a model.
jss_joint_contract_fixture <- function(case, rep_idx, settings, deterministic_offset = 0) {
  attempt <- list(elapsed = 1, warnings = character(), error = NA_character_)
  rows <- lapply(c("rs_separate", "rs_joint"), function(method) {
    row <- jss_joint_base_result_row(case, rep_idx, method, attempt)
    row$success <- TRUE
    row$converged <- TRUE
    row$retained <- TRUE
    row$stop_reason <- "converged"
    row$outer_iterations <- 1L
    row$outer_log_lik_change <- 0
    row$outer_stop_crit <- 0.01
    row$hit_outer_limit <- FALSE
    for (nm in jss_joint_required_retained_metrics()) row[[nm]] <- 1 + deterministic_offset
    row$benchmark_neg_log_score <- -row$test_log_score_per_obs
    row$heldout_variogram_nsim <- 1L
    row$failure_type <- "none"
    row
  })
  out <- jss_joint_bind_rows(rows)
  out$paired_seed <- settings$seed + 3000L + rep_idx * 100L
  out
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
  row$retained <- isTRUE(row$success) && isTRUE(row$converged)
  row$stop_reason <- as.character(fit$convergence$stop_reason %||% NA_character_)
  row$outer_iterations <- as.integer(fit$convergence$outer_iterations %||% NA_integer_)
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
    seed = 900000L + rep_idx * 10L
  )
  row$heldout_variogram_score_p05 <- variogram$score_p05
  row$heldout_variogram_score_p2 <- variogram$score_p2
  row$heldout_variogram_nsim <- as.integer(variogram$nsim)
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
  row$failure_type <- if (!isTRUE(row$success)) {
    row$stop_reason <- "invalid_likelihood"
    "invalid_loglik"
  } else if (!isTRUE(row$converged)) {
    paste0("optimizer_nonconvergence:", row$stop_reason %||% "unspecified")
  } else {
    "none"
  }
  row
}

jss_joint_run_case_rep <- function(case, rep_idx, settings, cfg) {
  # All cells in a replicate share the data-generation seed.  This provides the
  # paired random numbers needed for one-factor design contrasts.
  seed <- settings$seed + 3000L + rep_idx * 100L
  covariates <- jss_joint_make_covariates(case, seed + 17L)
  train_dat <- jss_joint_simulate_case_data(case, seed, response_seed_offset = 0L, covariates = covariates)
  test_dat <- jss_joint_simulate_case_data(case, seed, response_seed_offset = 500000L, covariates = covariates)

  sep <- jss_joint_fit_model(train_dat, case, include_dlcopdpar = FALSE, cfg = cfg)
  joint <- jss_joint_fit_model(train_dat, case, include_dlcopdpar = TRUE, cfg = cfg)

  out <- jss_joint_bind_rows(list(
    jss_joint_score_fit(sep$fit, train_dat, test_dat, case, rep_idx, "rs_separate", sep, cfg),
    jss_joint_score_fit(joint$fit, train_dat, test_dat, case, rep_idx, "rs_joint", joint, cfg)
  ))
  out$paired_seed <- as.integer(seed)
  out
}

jss_joint_checkpoint_dir <- function(settings) {
  file.path(settings$out_dir, "checkpoints", jss_joint_module_id(), "paired-one-factor-v6")
}

jss_joint_short_hash <- function(x) {
  text <- paste(capture.output(dput(x)), collapse = "")
  ints <- utf8ToInt(enc2utf8(text))
  checksum <- if (length(ints)) sum((seq_along(ints) %% 1009L + 1L) * ints) %% 2147483629 else 0
  paste0(nchar(text), "-", format(checksum, scientific = FALSE, trim = TRUE))
}

jss_joint_sha256_file <- function(path) {
  if (!requireNamespace("digest", quietly = TRUE)) {
    stop("The digest package is required for optimizer checkpoint identity verification.", call. = FALSE)
  }
  digest::digest(file = path, algo = "sha256", serialize = FALSE)
}

jss_joint_content_sha256 <- function(x) {
  if (!requireNamespace("digest", quietly = TRUE)) {
    stop("The digest package is required for optimizer result hashing.", call. = FALSE)
  }
  digest::digest(x, algo = "sha256", serialize = TRUE)
}

jss_joint_portable_result_sha256 <- function(x) {
  x <- x[, setdiff(names(x), c("result_content_sha256", "result_portable_sha256", "profile")), drop = FALSE]
  x <- x[order(as.character(x$method)), sort(names(x)), drop = FALSE]
  encode <- function(value) {
    if (is.logical(value)) return(ifelse(is.na(value), "NA", ifelse(value, "TRUE", "FALSE")))
    if (is.numeric(value)) return(ifelse(is.na(value), "NA", sprintf("%.15g", value)))
    value <- as.character(value)
    value[is.na(value) | !nzchar(value)] <- "NA"
    enc2utf8(value)
  }
  text <- paste(vapply(x, function(column) paste(encode(column), collapse = "\r"), character(1)), collapse = "\n")
  digest::digest(text, algo = "sha256", serialize = FALSE)
}

jss_joint_portable_frame_sha256 <- function(x) {
  x <- x[, sort(names(x)), drop = FALSE]
  encode <- function(value) {
    if (is.logical(value)) return(ifelse(is.na(value), "NA", ifelse(value, "TRUE", "FALSE")))
    if (is.numeric(value)) return(ifelse(is.na(value), "NA", sprintf("%.15g", value)))
    value <- as.character(value); value[is.na(value) | !nzchar(value)] <- "NA"; enc2utf8(value)
  }
  rows <- if (nrow(x)) apply(as.data.frame(lapply(x, encode), stringsAsFactors = FALSE),
    1L, paste, collapse = "\r") else character()
  digest::digest(paste(sort(rows), collapse = "\n"), algo = "sha256", serialize = FALSE)
}

jss_joint_normalize_results_csv <- function(x) {
  prototype <- jss_joint_contract_fixture(jss_joint_case_definitions()[1L, , drop = FALSE],
    1L, list(seed = 1L))
  for (nm in intersect(names(prototype), names(x))) {
    if (is.integer(x[[nm]]) && is.double(prototype[[nm]])) {
      x[[nm]] <- as.numeric(x[[nm]])
    } else if (is.logical(x[[nm]]) && all(is.na(x[[nm]])) && !is.logical(prototype[[nm]])) {
      x[[nm]] <- switch(typeof(prototype[[nm]]),
        integer = as.integer(x[[nm]]), double = as.numeric(x[[nm]]),
        character = as.character(x[[nm]]), x[[nm]])
    }
  }
  if (is.logical(x$warnings) && all(is.na(x$warnings))) x$warnings <- rep("", nrow(x))
  if (is.logical(x$error) && all(is.na(x$error))) x$error <- rep(NA_character_, nrow(x))
  if (is.character(x$warnings)) x$warnings[is.na(x$warnings)] <- ""
  if (!is.character(x$warnings) || !is.character(x$error)) {
    stop("Optimizer results CSV diagnostics did not normalize to character.", call. = FALSE)
  }
  x
}

jss_joint_read_canonical_results_csv <- function(path) {
  jss_joint_normalize_results_csv(utils::read.csv(path, stringsAsFactors = FALSE))
}

jss_joint_runtime_identity <- function() {
  ext <- extSoftVersion()
  portable_library <- function(value) {
    if (is.na(value) || !nzchar(value)) return("R-default-or-unreported")
    if (grepl("^([A-Za-z]:[/\\\\]|/)", value)) basename(value) else value
  }
  blas <- if ("BLAS" %in% names(ext)) unname(ext[["BLAS"]]) else NA_character_
  lapack <- if ("LAPACK" %in% names(ext)) unname(ext[["LAPACK"]]) else NA_character_
  list(
    timestamp_utc = format(Sys.time(), "%Y-%m-%dT%H:%M:%SZ", tz = "UTC"),
    pid = as.integer(Sys.getpid()),
    host = unname(Sys.info()[["nodename"]]),
    os = paste(unname(Sys.info()[c("sysname", "release", "machine")]), collapse = " "),
    platform = R.version$platform,
    r_version = R.version.string,
    rng_kind = paste(RNGkind(), collapse = "/"),
    blas = portable_library(blas),
    lapack = portable_library(lapack)
  )
}

jss_joint_checkout_package_identity <- function(settings) {
  if (!is.null(settings$checkpoint_package_identity)) return(settings$checkpoint_package_identity)
  root <- normalizePath(
    settings$package_source_root %||% settings$root %||% getwd(),
    winslash = "/", mustWork = TRUE
  )
  files <- c(
    list.files(file.path(root, "R"), pattern = "[.]R$", full.names = TRUE),
    file.path(root, c("DESCRIPTION", "NAMESPACE"))
  )
  files <- sort(normalizePath(files, winslash = "/", mustWork = TRUE))
  if (!length(files)) stop("No checkout package sources were available for fingerprinting.", call. = FALSE)
  relative <- substring(files, nchar(root) + 2L)
  hashes <- vapply(files, jss_joint_sha256_file, character(1))
  manifest <- paste(relative, hashes, sep = "\t", collapse = "\n")
  description <- read.dcf(file.path(root, "DESCRIPTION"))
  list(
    package = "gamlss.longitudinal",
    version = unname(description[1L, "Version"]),
    checkout_path = root,
    source_sha256 = digest::digest(manifest, algo = "sha256", serialize = FALSE),
    fingerprint_scope = "sorted R/*.R + DESCRIPTION + NAMESPACE",
    source_file_count = length(files),
    load_strategy = "pkgload_checkout"
  )
}

jss_joint_producer_fingerprint <- function(settings) {
  if (!is.null(settings$checkpoint_producer_fingerprint)) {
    return(settings$checkpoint_producer_fingerprint)
  }
  source_path <- settings$producer_source_path %||% file.path(
    settings$root %||% getwd(), "paper", "R", "03-joint-vs-separate-optimization.R"
  )
  if (!file.exists(source_path)) {
    stop("Cannot fingerprint optimizer producer source: ", source_path, call. = FALSE)
  }
  jss_joint_sha256_file(source_path)
}

jss_joint_verify_checkout_package <- function(settings, expected_identity, load_checkout = TRUE) {
  root <- expected_identity$checkout_path
  if (isTRUE(load_checkout)) {
    if (!requireNamespace("pkgload", quietly = TRUE)) {
      stop("pkgload is required to load the verified checkout package on optimizer workers.", call. = FALSE)
    }
    pkgload::load_all(root, quiet = TRUE, export_all = TRUE, helpers = FALSE)
  }
  actual <- jss_joint_checkout_package_identity(utils::modifyList(
    settings,
    list(package_source_root = root, checkpoint_package_identity = NULL)
  ))
  namespace_path <- normalizePath(
    getNamespaceInfo(asNamespace(expected_identity$package), "path"),
    winslash = "/", mustWork = TRUE
  )
  loaded_version <- as.character(utils::packageVersion(expected_identity$package))
  valid <- identical(actual$source_sha256, expected_identity$source_sha256) &&
    identical(actual$version, expected_identity$version) &&
    identical(loaded_version, expected_identity$version) &&
    identical(namespace_path, root)
  if (!valid) {
    stop(
      "Verified checkout package identity mismatch: expected ", expected_identity$version,
      " at ", root, " with SHA-256 ", expected_identity$source_sha256,
      "; loaded ", loaded_version, " at ", namespace_path,
      " with checkout SHA-256 ", actual$source_sha256, ".",
      call. = FALSE
    )
  }
  c(expected_identity, list(
    verified = TRUE,
    verified_package_path = namespace_path,
    verified_package_version = loaded_version,
    verified_source_sha256 = actual$source_sha256,
    inherited_rlibs_user = Sys.getenv("R_LIBS_USER", unset = ""),
    verified_libpaths = paste(normalizePath(.libPaths(), winslash = "/", mustWork = FALSE), collapse = ";")
  ))
}

jss_joint_reverify_task_sources <- function(settings) {
  expected_package <- settings$checkpoint_package_identity %||% jss_joint_checkout_package_identity(settings)
  actual_package <- jss_joint_checkout_package_identity(utils::modifyList(settings, list(
    checkpoint_package_identity = NULL,
    package_source_root = expected_package$checkout_path
  )))
  expected_producer <- settings$checkpoint_producer_fingerprint %||% jss_joint_producer_fingerprint(settings)
  actual_producer <- jss_joint_producer_fingerprint(utils::modifyList(settings, list(
    checkpoint_producer_fingerprint = NULL
  )))
  if (!identical(actual_package$source_sha256, expected_package$source_sha256) ||
      !identical(actual_package$version, expected_package$version) ||
      !identical(actual_producer, expected_producer)) {
    stop("Optimizer checkout or producer changed after attestation.", call. = FALSE)
  }
  invisible(list(package = actual_package, producer_sha256 = actual_producer))
}

jss_joint_reverify_checkpoint_or_quarantine <- function(settings, checkpoint) {
  tryCatch(jss_joint_reverify_task_sources(settings), error = function(e) {
    quarantine <- jss_joint_quarantine_checkpoint(checkpoint, reason = "source-identity-changed")
    stop(conditionMessage(e), "; checkpoint quarantined at ", quarantine %||% "removed", call. = FALSE)
  })
}

jss_joint_checkpoint_spec <- function(settings, case, rep_idx, cfg) {
  design_fields <- c(
    "case_id", "base_case_id", "contrast_factor", "family", "copula", "n",
    "time_points", "mu_strength", "sigma_strength", "theta_strength", "time_shape"
  )
  list(
    schema_version = 6L,
    producer_fingerprint = jss_joint_producer_fingerprint(settings),
    producer_fingerprint_algorithm = "SHA-256",
    package_identity = jss_joint_checkout_package_identity(settings),
    design = as.list(case[1L, design_fields, drop = FALSE]),
    replicate = as.integer(rep_idx),
    package_seed = as.integer(settings$seed),
    paired_seed = as.integer(settings$seed + 3000L + rep_idx * 100L),
    covariate_seed_offset = 17L,
    test_response_seed_offset = 500000L,
    variogram_seed = 900000L + rep_idx * 10L,
    optimizer_and_scoring_controls = cfg[setdiff(names(cfg), c("resume", "workers", "reps"))],
    methods = c("rs_separate", "rs_joint"),
    rng_kind = RNGkind()
  )
}

jss_joint_allowed_stop_reasons <- function() {
  c("converged", "max_iterations", "max_stall", "objective_deterioration",
    "invalid_likelihood", "numerical_failure", "time_limit", "fit_error")
}

jss_joint_required_retained_metrics <- function() {
  c(
    "elapsed_sec", "train_marginal_loglik", "train_copula_loglik",
    "train_joint_loglik", "marginal_loglik", "copula_loglik", "joint_loglik",
    "train_rmse_mu", "train_rmse_sigma", "train_rmse_theta", "train_rmse_tau",
    "test_marginal_log_score", "test_copula_log_score", "test_joint_log_score",
    "test_log_score_per_obs", "heldout_variogram_score_p05",
    "heldout_variogram_score_p2", "benchmark_mean_rmse",
    "benchmark_neg_log_score", "benchmark_variogram_score_p05",
    "benchmark_variogram_score_p2", "benchmark_theta_time_abs_error",
    "max_abs_param_error", "max_rel_param_error"
  )
}

jss_joint_validate_result_payload <- function(x, expected_spec, case, rep_idx) {
  reasons <- character()
  required <- c(
    "case_id", "joint_review_rep", "method", "paired_seed", "success",
    "converged", "retained", "stop_reason", "failure_type", "warnings", "error"
  )
  if (!is.data.frame(x) || nrow(x) != 2L || !all(required %in% names(x))) {
    return(list(valid = FALSE, reasons = "result is not the registered two-row payload"))
  }
  prototype <- jss_joint_base_result_row(case, as.integer(rep_idx), "rs_separate",
    list(elapsed = 0, warnings = character(), error = NA_character_))
  prototype$paired_seed <- as.integer(expected_spec$paired_seed)
  if (!identical(names(x), names(prototype))) {
    reasons <- c(reasons, "result has extra, missing, or reordered checkpoint fields")
  } else {
    wrong_type <- names(prototype)[vapply(names(prototype), function(nm)
      !identical(typeof(x[[nm]]), typeof(prototype[[nm]])), logical(1))]
    if (length(wrong_type)) {
      reasons <- c(reasons, paste0("checkpoint fields have invalid declared types: ",
        paste(wrong_type, collapse = ", ")))
    }
  }
  if (!all(vapply(x[c("success", "converged", "retained")], is.logical, logical(1)))) {
    reasons <- c(reasons, "status fields must be logical")
  }
  if (!all(vapply(x[c("case_id", "method", "stop_reason", "failure_type", "warnings", "error")],
      is.character, logical(1)))) reasons <- c(reasons, "label and diagnostic fields must be character")
  if (anyNA(x$warnings)) reasons <- c(reasons, "warnings must be nonmissing character diagnostics")
  metric_schema <- intersect(jss_joint_required_retained_metrics(), names(x))
  if (length(metric_schema) && any(!vapply(x[metric_schema], is.numeric, logical(1)))) {
    reasons <- c(reasons, "registered metric fields must be numeric regardless of fit status")
  }
  integral_fields <- intersect(c("joint_review_rep", "paired_seed", "outer_iterations",
    "heldout_variogram_nsim"), names(x))
  if (any(vapply(integral_fields, function(nm) {
    value <- x[[nm]]
    !is.integer(value) || any(!is.na(value) & (!is.finite(value) | value != floor(value)))
  }, logical(1)))) reasons <- c(reasons, "registered key/count fields must be integer typed and integral")
  eligible <- x$success %in% TRUE & x$converged %in% TRUE
  if (any(!(x$success %in% TRUE) & x$converged %in% TRUE)) {
    reasons <- c(reasons, "a converged fit cannot be marked unsuccessful")
  }
  if (!identical(as.logical(x$retained), eligible)) reasons <- c(reasons, "retained status is inconsistent")
  if (!all(x$case_id == case$case_id[[1L]]) || !all(x$joint_review_rep == rep_idx)) {
    reasons <- c(reasons, "scenario or replicate key is inconsistent")
  }
  if (!identical(sort(as.character(x$method)), c("rs_joint", "rs_separate"))) {
    reasons <- c(reasons, "method keys are not registered")
  }
  if (length(unique(x$paired_seed)) != 1L ||
      !identical(as.integer(unique(x$paired_seed)), expected_spec$paired_seed)) {
    reasons <- c(reasons, "paired seed is inconsistent")
  }
  numeric_fields <- names(x)[vapply(x, is.numeric, logical(1))]
  numeric_fields <- setdiff(numeric_fields, c("joint_review_rep", "paired_seed"))
  for (nm in numeric_fields) {
    values <- x[[nm]]
    if (any(!is.na(values) & (!is.finite(values) | abs(values) > 1e12))) {
      reasons <- c(reasons, paste0("numeric field out of contract: ", nm))
    }
  }
  nonnegative <- intersect(c(
    "elapsed_sec", "outer_iterations", "heldout_variogram_nsim",
    "train_rmse_mu", "train_rmse_sigma", "train_rmse_theta", "train_rmse_tau",
    "heldout_variogram_score_p05", "heldout_variogram_score_p2",
    "benchmark_mean_rmse",
    "benchmark_variogram_score_p05", "benchmark_variogram_score_p2",
    "benchmark_theta_time_abs_error", "max_abs_param_error", "max_rel_param_error"
  ), names(x))
  if (any(vapply(nonnegative, function(nm) any(!is.na(x[[nm]]) & x[[nm]] < 0), logical(1)))) {
    reasons <- c(reasons, "nonnegative metric is negative")
  }
  if ("outer_stop_crit" %in% names(x) && any(!is.na(x$outer_stop_crit) & x$outer_stop_crit <= 0)) {
    reasons <- c(reasons, "outer stopping tolerance must be positive")
  }
  flag_fields <- intersect(c("hit_outer_limit"), names(x))
  if (any(!vapply(x[flag_fields], is.logical, logical(1)))) reasons <- c(reasons, "optimizer flags must be logical")
  if (any(eligible)) {
    missing_metrics <- setdiff(jss_joint_required_retained_metrics(), names(x))
    if (length(missing_metrics)) reasons <- c(reasons, "retained payload omits registered metrics")
    present <- intersect(jss_joint_required_retained_metrics(), names(x))
    if (length(present) && any(!vapply(x[present], is.numeric, logical(1)))) {
      reasons <- c(reasons, "retained payload metrics must be numeric")
    } else if (length(present) && any(vapply(present, function(nm) any(!is.finite(x[[nm]][eligible])), logical(1)))) {
      reasons <- c(reasons, "retained payload has missing or nonfinite metrics")
    }
  }
  allowed_stops <- jss_joint_allowed_stop_reasons()
  for (i in seq_len(nrow(x))) {
    stop_reason <- as.character(x$stop_reason[[i]])
    failure_type <- as.character(x$failure_type[[i]])
    if (eligible[[i]]) {
      if (!identical(failure_type, "none") || !identical(stop_reason, "converged")) {
        reasons <- c(reasons, "retained row has an unregistered status label")
      }
      if (!is.na(x$error[[i]]) && nzchar(x$error[[i]])) {
        reasons <- c(reasons, "retained row cannot carry an error diagnostic")
      }
    } else if (x$success[[i]] %in% TRUE) {
      expected_failure <- paste0("optimizer_nonconvergence:", stop_reason)
      if (is.na(stop_reason) || !stop_reason %in% setdiff(allowed_stops, "converged") ||
          !identical(failure_type, expected_failure)) {
        reasons <- c(reasons, "nonconverged row has an unregistered failure label")
      }
      if (!is.na(x$error[[i]]) && nzchar(x$error[[i]])) {
        reasons <- c(reasons, "nonconverged fit cannot carry a fit-error diagnostic")
      }
    } else if (x$converged[[i]] %in% TRUE) {
      reasons <- c(reasons, "failed row cannot claim convergence")
    } else if (is.na(failure_type) || !failure_type %in% c("invalid_loglik", "fit_unavailable", "error")) {
      reasons <- c(reasons, "failed row has an unregistered failure label")
    } else if (identical(failure_type, "error") && (is.na(x$error[[i]]) || !nzchar(x$error[[i]]))) {
      reasons <- c(reasons, "error failure lacks an error diagnostic")
    } else if (identical(failure_type, "invalid_loglik") && !identical(stop_reason, "invalid_likelihood")) {
      reasons <- c(reasons, "invalid-likelihood failure has an inconsistent stop reason")
    } else if (failure_type %in% c("fit_unavailable", "error") && !identical(stop_reason, "fit_error")) {
      reasons <- c(reasons, "fit failure has an inconsistent stop reason")
    }
  }
  if (any(eligible)) {
    close <- function(a, b) isTRUE(all.equal(as.numeric(a[eligible]), as.numeric(b[eligible]), tolerance = 1e-10))
    aliases <- list(c("marginal_loglik", "train_marginal_loglik"),
      c("copula_loglik", "train_copula_loglik"), c("joint_loglik", "train_joint_loglik"),
      c("benchmark_mean_rmse", "train_rmse_mu"),
      c("benchmark_variogram_score_p05", "heldout_variogram_score_p05"),
      c("benchmark_variogram_score_p2", "heldout_variogram_score_p2"),
      c("benchmark_theta_time_abs_error", "train_rmse_tau"))
    if (any(vapply(aliases, function(pair) !close(x[[pair[[1L]]]], x[[pair[[2L]]]]), logical(1)))) {
      reasons <- c(reasons, "registered metric aliases are inconsistent")
    }
    if (!close(x$benchmark_neg_log_score, -x$test_log_score_per_obs)) {
      reasons <- c(reasons, "negative-log-score alias is inconsistent")
    }
  }
  list(valid = !length(reasons), reasons = unique(reasons))
}

jss_joint_checkpoint_path <- function(settings, case, rep_idx, cfg) {
  spec <- jss_joint_checkpoint_spec(settings, case, rep_idx, cfg)
  file.path(
    jss_joint_checkpoint_dir(settings),
    paste0("run-", jss_joint_short_hash(spec$optimizer_and_scoring_controls)),
    sprintf(
      "%s-design-%s-rep-%04d.rds",
      case$case_id[[1L]],
      jss_joint_short_hash(spec$design),
      as.integer(rep_idx)
    )
  )
}

jss_joint_read_checkpoint <- function(path, settings, case, rep_idx, cfg) {
  if (!file.exists(path) || file.info(path)$size <= 0L) return(NULL)
  checkpoint <- tryCatch(readRDS(path), error = function(e) NULL)
  expected_spec <- jss_joint_checkpoint_spec(settings, case, rep_idx, cfg)
  if (!is.list(checkpoint) || !identical(checkpoint$schema_version, 6L) ||
      !identical(checkpoint$spec, expected_spec) || !is.list(checkpoint$provenance) ||
      !is.character(checkpoint$result_sha256) || length(checkpoint$result_sha256) != 1L) return(NULL)
  provenance <- checkpoint$provenance
  observed <- jss_joint_runtime_identity()
  runtime_fields <- c("host", "os", "platform", "r_version", "rng_kind", "blas", "lapack")
  expected_provenance_names <- c(names(observed), "package_identity_verified",
    "package_checkout_path", "package_source_sha256", "producer_path_verified",
    "producer_sha256_verified")
  if (!setequal(names(provenance), expected_provenance_names) ||
      !is.integer(provenance$pid) || length(provenance$pid) != 1L ||
      is.na(provenance$pid) || provenance$pid < 1L ||
      any(vapply(runtime_fields, function(field)
        !identical(provenance[[field]], observed[[field]]), logical(1))) ||
      !identical(provenance$package_source_sha256, expected_spec$package_identity$source_sha256) ||
      !identical(provenance$package_checkout_path, expected_spec$package_identity$checkout_path) ||
      !identical(provenance$producer_sha256_verified, expected_spec$producer_fingerprint)) return(NULL)
  validation <- jss_joint_validate_result_payload(checkpoint$result, expected_spec, case, rep_idx)
  if (!isTRUE(validation$valid) ||
      !identical(jss_joint_content_sha256(checkpoint$result), checkpoint$result_sha256) ||
      !identical(jss_joint_portable_result_sha256(checkpoint$result), checkpoint$result_portable_sha256)) return(NULL)
  checkpoint
}

jss_joint_checkpoint_valid <- function(path, settings, case, rep_idx, cfg) {
  is.list(jss_joint_read_checkpoint(path, settings, case, rep_idx, cfg))
}

jss_joint_write_checkpoint <- function(result, path, settings, case, rep_idx, cfg) {
  spec <- jss_joint_checkpoint_spec(settings, case, rep_idx, cfg)
  validation <- jss_joint_validate_result_payload(result, spec, case, rep_idx)
  if (!isTRUE(validation$valid)) {
    stop("Refusing invalid optimizer checkpoint payload: ", paste(validation$reasons, collapse = "; "), call. = FALSE)
  }
  dir.create(dirname(path), recursive = TRUE, showWarnings = FALSE)
  temporary <- paste0(path, ".tmp-", Sys.getpid(), "-", sprintf("%09d", sample.int(1e9, 1L)))
  verified_identity <- if (exists(".jss_joint_verified_worker_identity", envir = .GlobalEnv, inherits = FALSE)) {
    get(".jss_joint_verified_worker_identity", envir = .GlobalEnv)
  } else settings$checkpoint_verified_identity %||% list(verified = NA)
  provenance <- c(jss_joint_runtime_identity(), list(
    package_identity_verified = isTRUE(verified_identity$verified),
    package_checkout_path = spec$package_identity$checkout_path,
    package_source_sha256 = spec$package_identity$source_sha256,
    producer_path_verified = normalizePath(
      settings$producer_source_path %||% file.path(settings$root %||% getwd(), "paper", "R", "03-joint-vs-separate-optimization.R"),
      winslash = "/", mustWork = FALSE
    ),
    producer_sha256_verified = spec$producer_fingerprint
  ))
  saveRDS(
    list(
      schema_version = 6L,
      spec = spec,
      result_sha256 = jss_joint_content_sha256(result),
      result_portable_sha256 = jss_joint_portable_result_sha256(result),
      provenance = provenance,
      result = result
    ),
    temporary
  )
  if (file.exists(path) && !file.remove(path)) {
    stop("Could not replace optimizer checkpoint: ", path, call. = FALSE)
  }
  if (!file.rename(temporary, path)) {
    stop("Could not atomically install optimizer checkpoint: ", path, call. = FALSE)
  }
  invisible(path)
}

jss_joint_acquire_run_lock <- function(settings) {
  lock <- paste0(jss_joint_checkpoint_dir(settings), ".run-lock")
  dir.create(dirname(lock), recursive = TRUE, showWarnings = FALSE)
  if (!dir.create(lock, recursive = FALSE, showWarnings = FALSE)) {
    stop("Another optimizer runner holds the output lock: ", lock, call. = FALSE)
  }
  saveRDS(jss_joint_runtime_identity(), file.path(lock, "owner.rds"))
  lock
}

jss_joint_release_run_lock <- function(lock) {
  if (is.character(lock) && length(lock) == 1L && dir.exists(lock) && grepl("[.]run-lock$", lock)) {
    unlink(lock, recursive = TRUE, force = TRUE)
  }
  invisible(NULL)
}

jss_joint_task_specs <- function(settings, cases, cfg) {
  tasks <- list()
  task_id <- 0L
  for (case_idx in seq_len(nrow(cases))) {
    for (rep_idx in seq_len(cfg$reps)) {
      task_id <- task_id + 1L
      case <- cases[case_idx, , drop = FALSE]
      tasks[[task_id]] <- list(
        task_id = task_id,
        case = case,
        rep_idx = rep_idx,
        checkpoint = jss_joint_checkpoint_path(settings, case, rep_idx, cfg)
      )
    }
  }
  checkpoints <- vapply(tasks, `[[`, character(1), "checkpoint")
  if (anyDuplicated(checkpoints)) stop("Optimizer task design produced duplicate checkpoint paths.", call. = FALSE)
  tasks
}

jss_joint_requested_workers <- function(cfg) {
  requested <- suppressWarnings(as.integer(
    cfg$workers %||% Sys.getenv("GAMLSS_LONGITUDINAL_JSS_JVS_WORKERS", unset = "1")
  ))
  if (!is.finite(requested) || requested < 1L) requested <- 1L
  requested
}

jss_joint_worker_count <- function(cfg, pending_tasks) {
  requested <- jss_joint_requested_workers(cfg)
  if (pending_tasks < 1L) return(0L)
  min(requested, as.integer(pending_tasks))
}

jss_joint_initialize_worker <- function(root, parent_libpaths, rlibs_user, expected_identity,
                                        expected_producer_fingerprint) {
  old_working_directory <- getwd()
  on.exit(setwd(old_working_directory), add = TRUE)
  if (nzchar(rlibs_user)) Sys.setenv(R_LIBS_USER = rlibs_user)
  .libPaths(unique(c(parent_libpaths, .libPaths())))
  setwd(root)
  if (!requireNamespace("pkgload", quietly = TRUE)) {
    stop("pkgload is required for verified checkout loading.", call. = FALSE)
  }
  pkgload::load_all(root, quiet = TRUE, export_all = TRUE, helpers = FALSE)
  source(file.path(root, "paper", "R", "03-joint-vs-separate-optimization.R"), local = .GlobalEnv)
  actual_producer_fingerprint <- jss_joint_producer_fingerprint(list(root = root))
  if (!identical(actual_producer_fingerprint, expected_producer_fingerprint)) {
    stop("Optimizer producer source changed between parent fingerprinting and worker initialization.", call. = FALSE)
  }
  identity <- jss_joint_verify_checkout_package(
    list(root = root, package_source_root = root), expected_identity, load_checkout = FALSE
  )
  identity$verified_producer_path <- normalizePath(
    file.path(root, "paper", "R", "03-joint-vs-separate-optimization.R"),
    winslash = "/", mustWork = TRUE
  )
  identity$verified_producer_sha256 <- actual_producer_fingerprint
  assign(".jss_joint_verified_worker_identity", identity, envir = .GlobalEnv)
  identity
}

jss_joint_run_parallel_task <- function(task, run_case_rep_fn, settings, cfg) {
  jss_joint_reverify_task_sources(settings)
  result <- run_case_rep_fn(task$case, task$rep_idx, settings, cfg)
  jss_joint_reverify_task_sources(settings)
  spec <- jss_joint_checkpoint_spec(settings, task$case, task$rep_idx, cfg)
  validation <- jss_joint_validate_result_payload(result, spec, task$case, task$rep_idx)
  if (!isTRUE(validation$valid)) {
    quarantine <- jss_joint_quarantine_checkpoint(task$checkpoint, reason = "worker-forged-payload")
    return(list(ok = FALSE, result = NULL, result_sha256 = NA_character_,
      worker_pid = Sys.getpid(), quarantine = quarantine,
      validation_reasons = validation$reasons,
      package_identity = get(".jss_joint_verified_worker_identity", envir = .GlobalEnv)))
  }
  result_sha256 <- jss_joint_content_sha256(result)
  jss_joint_reverify_task_sources(settings)
  jss_joint_write_checkpoint(
    result, task$checkpoint, settings, task$case, task$rep_idx, cfg
  )
  jss_joint_reverify_checkpoint_or_quarantine(settings, task$checkpoint)
  durable <- jss_joint_read_checkpoint(task$checkpoint, settings, task$case, task$rep_idx, cfg)
  valid <- is.list(durable) && identical(durable$result_sha256, result_sha256) &&
    identical(durable$result, result)
  quarantine <- NA_character_
  if (!isTRUE(valid)) {
    quarantine <- jss_joint_quarantine_checkpoint(task$checkpoint, reason = "worker-invalid")
  }
  jss_joint_reverify_checkpoint_or_quarantine(settings, task$checkpoint)
  list(
    ok = isTRUE(valid), result = result, result_sha256 = result_sha256,
    provenance = if (is.list(durable)) durable$provenance else NULL,
    worker_pid = Sys.getpid(),
    quarantine = quarantine,
    package_identity = get(".jss_joint_verified_worker_identity", envir = .GlobalEnv)
  )
}

jss_joint_quarantine_checkpoint <- function(path, reason = "invalid") {
  if (!file.exists(path)) return(NA_character_)
  quarantine <- paste0(
    path, ".", reason, "-pid-", Sys.getpid(), "-",
    format(Sys.time(), "%Y%m%d%H%M%S")
  )
  if (file.rename(path, quarantine)) {
    return(normalizePath(quarantine, winslash = "/", mustWork = FALSE))
  }
  if (!file.remove(path)) {
    stop("Could not quarantine or remove invalid optimizer checkpoint: ", path, call. = FALSE)
  }
  NA_character_
}

jss_joint_run_simulation_fixed <- function(settings, cases, cfg = NULL,
                                           run_case_rep_fn = jss_joint_run_case_rep,
                                           acquire_lock = TRUE, validate_design = TRUE) {
  if (is.null(cfg)) cfg <- jss_joint_simulation_settings(settings)
  if (isTRUE(validate_design)) jss_joint_validate_case_definitions(cases)
  if (isTRUE(acquire_lock)) {
    run_lock <- jss_joint_acquire_run_lock(settings)
    on.exit(jss_joint_release_run_lock(run_lock), add = TRUE)
  }
  expected_package_identity <- jss_joint_checkout_package_identity(settings)
  settings$checkpoint_package_identity <- expected_package_identity
  settings$checkpoint_producer_fingerprint <- jss_joint_producer_fingerprint(settings)
  uses_production_runner <- identical(run_case_rep_fn, jss_joint_run_case_rep)
  parent_package_identity <- if (uses_production_runner) {
    jss_joint_verify_checkout_package(settings, expected_package_identity, load_checkout = TRUE)
  } else {
    c(expected_package_identity, list(
      verified = NA,
      verified_package_path = expected_package_identity$checkout_path,
      verified_package_version = expected_package_identity$version,
      verified_source_sha256 = expected_package_identity$source_sha256,
      inherited_rlibs_user = Sys.getenv("R_LIBS_USER", unset = ""),
      verified_libpaths = paste(normalizePath(.libPaths(), winslash = "/", mustWork = FALSE), collapse = ";")
    ))
  }
  settings$checkpoint_verified_identity <- parent_package_identity
  tasks <- jss_joint_task_specs(settings, cases, cfg)
  rows <- vector("list", length(tasks))
  resumed <- rep(FALSE, length(tasks))
  worker_pid <- rep(NA_integer_, length(tasks))
  task_package_identity <- rep(list(parent_package_identity), length(tasks))
  task_provenance <- vector("list", length(tasks))
  result_sha256 <- rep(NA_character_, length(tasks))

  for (i in seq_along(tasks)) {
    task <- tasks[[i]]
    checkpoint <- if (isTRUE(cfg$resume)) {
      jss_joint_read_checkpoint(task$checkpoint, settings, task$case, task$rep_idx, cfg)
    } else NULL
    if (isTRUE(cfg$resume) && file.exists(task$checkpoint) && !is.list(checkpoint)) {
      jss_joint_quarantine_checkpoint(task$checkpoint, reason = "resume-invalid")
    }
    resumed[[i]] <- is.list(checkpoint)
    if (resumed[[i]]) {
      rows[[i]] <- checkpoint$result
      result_sha256[[i]] <- checkpoint$result_sha256
      task_provenance[[i]] <- checkpoint$provenance
    }
  }

  pending <- which(!resumed)
  workers_requested <- jss_joint_requested_workers(cfg)
  workers_used <- jss_joint_worker_count(cfg, length(pending))
  if (length(pending) && workers_used == 1L) {
    for (i in pending) {
      task <- tasks[[i]]
      jss_joint_reverify_task_sources(settings)
      rows[[i]] <- run_case_rep_fn(task$case, task$rep_idx, settings, cfg)
      jss_joint_reverify_task_sources(settings)
      payload_validation <- jss_joint_validate_result_payload(
        rows[[i]], jss_joint_checkpoint_spec(settings, task$case, task$rep_idx, cfg),
        task$case, task$rep_idx
      )
      if (!isTRUE(payload_validation$valid)) {
        jss_joint_quarantine_checkpoint(task$checkpoint, reason = "parent-forged-payload")
        stop("Serial optimizer task returned an invalid payload: ",
          paste(payload_validation$reasons, collapse = "; "), call. = FALSE)
      }
      result_sha256[[i]] <- jss_joint_content_sha256(rows[[i]])
      jss_joint_reverify_task_sources(settings)
      jss_joint_write_checkpoint(
        rows[[i]], task$checkpoint, settings, task$case, task$rep_idx, cfg
      )
      jss_joint_reverify_checkpoint_or_quarantine(settings, task$checkpoint)
      durable <- jss_joint_read_checkpoint(task$checkpoint, settings, task$case, task$rep_idx, cfg)
      if (!is.list(durable) || !identical(durable$result, rows[[i]]) ||
          !identical(durable$result_sha256, result_sha256[[i]])) {
        quarantine <- jss_joint_quarantine_checkpoint(task$checkpoint, reason = "parent-invalid")
        stop("Serial optimizer task produced an invalid checkpoint; quarantined at ", quarantine, call. = FALSE)
      }
      jss_joint_reverify_checkpoint_or_quarantine(settings, task$checkpoint)
      task_provenance[[i]] <- durable$provenance
      worker_pid[[i]] <- Sys.getpid()
    }
  } else if (length(pending)) {
    root <- normalizePath(settings$root %||% getwd(), winslash = "/", mustWork = TRUE)
    worker_log <- file.path(
      settings$logs_dir %||% jss_joint_checkpoint_dir(settings),
      "03-joint-vs-separate-optimization-workers.log"
    )
    dir.create(dirname(worker_log), recursive = TRUE, showWarnings = FALSE)
    cluster <- parallel::makePSOCKcluster(workers_used, outfile = worker_log)
    on.exit(parallel::stopCluster(cluster), add = TRUE)
    initialized <- parallel::clusterCall(
      cluster,
      jss_joint_initialize_worker,
      root = root,
      parent_libpaths = .libPaths(),
      rlibs_user = Sys.getenv("R_LIBS_USER", unset = ""),
      expected_identity = expected_package_identity
      , expected_producer_fingerprint = settings$checkpoint_producer_fingerprint
    )
    if (!all(vapply(initialized, function(x) is.list(x) && isTRUE(x$verified), logical(1)))) {
      stop("One or more optimizer PSOCK workers failed initialization.", call. = FALSE)
    }
    completed <- parallel::parLapplyLB(
      cluster,
      tasks[pending],
      jss_joint_run_parallel_task,
      run_case_rep_fn = run_case_rep_fn,
      settings = settings,
      cfg = cfg
    )
    for (j in seq_along(pending)) {
      i <- pending[[j]]
      if (!isTRUE(completed[[j]]$ok)) {
        stop(
          "Optimizer worker produced an invalid checkpoint; quarantine: ",
          completed[[j]]$quarantine %||% "unavailable", call. = FALSE
        )
      }
      rows[[i]] <- completed[[j]]$result
      if (!identical(jss_joint_content_sha256(rows[[i]]), completed[[j]]$result_sha256)) {
        stop("Optimizer worker returned content inconsistent with its attested hash.", call. = FALSE)
      }
      result_sha256[[i]] <- completed[[j]]$result_sha256
      task_provenance[[i]] <- completed[[j]]$provenance
      worker_pid[[i]] <- completed[[j]]$worker_pid
      task_package_identity[[i]] <- completed[[j]]$package_identity
    }
  }

  producer_fingerprint <- jss_joint_producer_fingerprint(settings)
  status <- lapply(seq_along(tasks), function(i) {
    task <- tasks[[i]]
    result <- rows[[i]]
    verified_identity <- task_package_identity[[i]]
    provenance <- task_provenance[[i]]
    data.frame(
      task_id = task$task_id,
      case_id = task$case$case_id[[1L]],
      joint_review_rep = task$rep_idx,
      paired_seed = unique(result$paired_seed)[[1L]],
      checkpoint = normalizePath(task$checkpoint, winslash = "/", mustWork = FALSE),
      checkpoint_schema_version = 6L,
      result_content_sha256 = result_sha256[[i]],
      result_portable_sha256 = jss_joint_portable_result_sha256(result),
      producer_fingerprint = producer_fingerprint,
      producer_fingerprint_algorithm = "SHA-256",
      package_version = expected_package_identity$version,
      package_fingerprint_scope = expected_package_identity$fingerprint_scope,
      package_source_file_count = expected_package_identity$source_file_count,
      package_checkout_path = expected_package_identity$checkout_path,
      package_source_sha256 = expected_package_identity$source_sha256,
      verified_package_path = verified_identity$verified_package_path,
      verified_package_version = verified_identity$verified_package_version,
      verified_source_sha256 = verified_identity$verified_source_sha256,
      verified_producer_path = verified_identity$verified_producer_path %||%
        normalizePath(file.path(settings$root, "paper", "R", "03-joint-vs-separate-optimization.R"), winslash = "/", mustWork = FALSE),
      verified_producer_sha256 = verified_identity$verified_producer_sha256 %||% producer_fingerprint,
      package_identity_verified = isTRUE(verified_identity$verified),
      verified_rlibs_user = verified_identity$inherited_rlibs_user,
      verified_libpaths = verified_identity$verified_libpaths,
      package_load_strategy = expected_package_identity$load_strategy,
      resumed = resumed[[i]],
      executed = !resumed[[i]],
      workers_requested = workers_requested,
      workers_used = if (resumed[[i]]) 0L else workers_used,
      execution_mode = if (resumed[[i]]) "resume_only" else if (workers_used > 1L) "psock" else "serial",
      worker_pid = if (resumed[[i]]) provenance$pid %||% NA_integer_ else worker_pid[[i]],
      checkpoint_timestamp_utc = provenance$timestamp_utc %||% NA_character_,
      execution_host = provenance$host %||% NA_character_,
      execution_os = provenance$os %||% NA_character_,
      execution_platform = provenance$platform %||% NA_character_,
      execution_r_version = provenance$r_version %||% NA_character_,
      execution_rng_kind = provenance$rng_kind %||% NA_character_,
      execution_blas = provenance$blas %||% NA_character_,
      execution_lapack = provenance$lapack %||% NA_character_,
      stringsAsFactors = FALSE
    )
  })
  out <- jss_joint_bind_rows(rows)
  portable_map <- setNames(vapply(status, function(x) x$result_portable_sha256[[1L]], character(1)),
    vapply(tasks, function(task) paste(task$case$case_id[[1L]], task$rep_idx, sep = "\r"), character(1)))
  out$result_portable_sha256 <- unname(portable_map[paste(out$case_id, out$joint_review_rep, sep = "\r")])
  attr(out, "checkpoint_status") <- jss_joint_bind_rows(status)
  out
}

jss_joint_run_simulation <- function(settings, cases, cfg = NULL,
                                     run_case_rep_fn = jss_joint_run_case_rep) {
  if (is.null(cfg)) cfg <- jss_joint_simulation_settings(settings)
  jss_joint_validate_case_definitions(cases)
  run_lock <- jss_joint_acquire_run_lock(settings)
  on.exit(jss_joint_release_run_lock(run_lock), add = TRUE)
  precision <- cfg$precision
  if (is.null(precision) || !is.data.frame(precision)) {
    return(jss_joint_run_simulation_fixed(
      settings, cases, cfg, run_case_rep_fn, acquire_lock = FALSE
    ))
  }
  required <- as.integer(precision$required_retained_pairs[[1L]])
  initial <- as.integer(precision$initial_attempts[[1L]])
  batch <- as.integer(precision$top_up_batch[[1L]])
  hard_cap <- as.integer(precision$max_attempts[[1L]])
  if (any(!is.finite(c(required, initial, batch, hard_cap))) ||
      required < 1L || initial < required || batch < 1L || hard_cap < initial) {
    stop("Optimizer adaptive precision registry is invalid.", call. = FALSE)
  }
  case_results <- setNames(vector("list", nrow(cases)), as.character(cases$case_id))
  case_status <- setNames(vector("list", nrow(cases)), as.character(cases$case_id))
  active <- as.character(cases$case_id)
  target_attempts <- setNames(rep(initial, nrow(cases)), as.character(cases$case_id))
  round <- 0L
  achievement <- NULL
  while (length(active)) {
    round <- round + 1L
    active_targets <- unique(target_attempts[active])
    if (length(active_targets) != 1L) stop("Adaptive optimizer cells lost registered batch alignment.", call. = FALSE)
    active_cases <- cases[match(active, cases$case_id), , drop = FALSE]
    batch_cfg <- cfg
    batch_cfg$reps <- as.integer(active_targets[[1L]])
    batch_result <- jss_joint_run_simulation_fixed(
      settings, active_cases, batch_cfg, run_case_rep_fn, acquire_lock = FALSE,
      validate_design = FALSE
    )
    batch_status <- attr(batch_result, "checkpoint_status")
    for (case_id in active) {
      case_results[[case_id]] <- batch_result[batch_result$case_id == case_id, , drop = FALSE]
      status_rows <- batch_status[batch_status$case_id == case_id, , drop = FALSE]
      existing <- case_status[[case_id]]
      previous_max <- if (is.null(existing) || !nrow(existing)) 0L else max(existing$joint_review_rep)
      new_rows <- status_rows[status_rows$joint_review_rep > previous_max, , drop = FALSE]
      new_rows$adaptive_round <- round
      case_status[[case_id]] <- if (is.null(existing) || !nrow(existing)) new_rows else
        jss_joint_bind_rows(list(existing, new_rows))
    }
    current <- jss_joint_bind_rows(case_results[!vapply(case_results, is.null, logical(1))])
    achievement <- jss_joint_precision_achievement(jss_joint_delta_table(current), precision)
    deficient <- achievement$case_id[!achievement$precision_met]
    at_cap <- achievement$case_id[achievement$attempted_pairs >= hard_cap]
    active <- setdiff(deficient, at_cap)
    if (length(active)) target_attempts[active] <- pmin(hard_cap, target_attempts[active] + batch)
  }
  out <- jss_joint_bind_rows(case_results[as.character(cases$case_id)])
  status <- jss_joint_bind_rows(case_status[as.character(cases$case_id)])
  status$task_id <- seq_len(nrow(status))
  status$registered_initial_attempts <- initial
  status$registered_top_up_batch <- batch
  status$registered_hard_cap <- hard_cap
  status$actual_attempts_in_cell <- ave(status$joint_review_rep, status$case_id, FUN = length)
  attr(out, "checkpoint_status") <- status
  attr(out, "precision_achievement") <- achievement
  out
}

jss_joint_delta_table <- function(results) {
  if (nrow(results) == 0L) {
    return(data.frame())
  }
  group_cols <- intersect(
    c(
      "case_id", "hypothesis_role", "joint_review_rep", "family", "copula", "design",
      "base_case_id", "contrast_factor", "contrast_label", "contrast_level", "paired_seed",
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
    row$rs_joint_converged <- joint$converged
    row$rs_separate_converged <- separate$converged
    row$rs_joint_retained <- jss_joint_eligible(joint)
    row$rs_separate_retained <- jss_joint_eligible(separate)
    paired_retained <- isTRUE(row$rs_joint_retained[[1L]]) && isTRUE(row$rs_separate_retained[[1L]])
    for (metric in metrics) {
      row[[paste0("rs_joint_", metric)]] <- jss_joint_numeric(joint[[metric]])
      row[[paste0("rs_separate_", metric)]] <- jss_joint_numeric(separate[[metric]])
      row[[paste0("delta_", metric)]] <- if (paired_retained) {
        jss_joint_numeric(joint[[metric]]) - jss_joint_numeric(separate[[metric]])
      } else {
        NA_real_
      }
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
    retained <- jss_joint_eligible(x)
    y <- x[retained, , drop = FALSE]
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
      n_attempted = nrow(x),
      n_fit_returned = sum(x$success %in% TRUE),
      n_successful = sum(retained),
      n_converged = sum(x$converged %in% TRUE),
      n_failed = sum(!retained),
      success_rate = mean(x$success %in% TRUE, na.rm = TRUE),
      convergence_rate = mean(x$converged %in% TRUE, na.rm = TRUE),
      retention_rate = mean(retained),
      median_elapsed_sec = jss_joint_median_finite(y$elapsed_sec),
      mean_train_joint_loglik = jss_joint_mean_finite(y$train_joint_loglik),
      mean_test_log_score_per_obs = jss_joint_mean_finite(y$test_log_score_per_obs),
      mean_heldout_variogram_score_p05 = jss_joint_mean_finite(y$heldout_variogram_score_p05),
      mean_heldout_variogram_score_p2 = jss_joint_mean_finite(y$heldout_variogram_score_p2),
      mean_train_rmse_theta = jss_joint_mean_finite(y$train_rmse_theta),
      mean_train_rmse_tau = jss_joint_mean_finite(y$train_rmse_tau),
      mean_benchmark_neg_log_score = jss_joint_mean_finite(y$benchmark_neg_log_score),
      mean_benchmark_variogram_score_p05 = jss_joint_mean_finite(y$benchmark_variogram_score_p05),
      mean_benchmark_variogram_score_p2 = jss_joint_mean_finite(y$benchmark_variogram_score_p2),
      mean_benchmark_mean_rmse = jss_joint_mean_finite(y$benchmark_mean_rmse),
      selected_review_cases = nrow(candidates),
      stringsAsFactors = FALSE
    )
  })
  summary <- jss_joint_bind_rows(rows)

  if (nrow(deltas) > 0L) {
    delta_key <- interaction(deltas$case_id, deltas$hypothesis_role, drop = TRUE, lex.order = TRUE)
    delta_summary <- jss_joint_bind_rows(lapply(split(seq_len(nrow(deltas)), delta_key), function(idx) {
      x <- deltas[idx, , drop = FALSE]
      paired_retained <- jss_joint_pair_eligible(x)
      delta_fields <- grep("^delta_", names(x), value = TRUE)
      x[!paired_retained, delta_fields] <- NA_real_
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
        n_attempted = nrow(x),
        n_fit_returned = sum(x$rs_joint_success %in% TRUE & x$rs_separate_success %in% TRUE),
        n_successful = sum(paired_retained),
        n_converged = NA_integer_,
        n_failed = sum(!paired_retained),
        success_rate = mean(x$rs_joint_success %in% TRUE & x$rs_separate_success %in% TRUE, na.rm = TRUE),
        convergence_rate = NA_real_,
        retention_rate = mean(paired_retained),
        median_elapsed_sec = jss_joint_median_finite(x$delta_elapsed_sec),
        mean_train_joint_loglik = jss_joint_mean_finite(x$delta_train_joint_loglik),
        mean_test_log_score_per_obs = jss_joint_mean_finite(x$delta_test_log_score_per_obs),
        mean_heldout_variogram_score_p05 = jss_joint_mean_finite(x$delta_heldout_variogram_score_p05),
        mean_heldout_variogram_score_p2 = jss_joint_mean_finite(x$delta_heldout_variogram_score_p2),
        mean_train_rmse_theta = jss_joint_mean_finite(x$delta_train_rmse_theta),
        mean_train_rmse_tau = jss_joint_mean_finite(x$delta_train_rmse_tau),
        mean_benchmark_neg_log_score = jss_joint_mean_finite(x$delta_benchmark_neg_log_score),
        mean_benchmark_variogram_score_p05 = jss_joint_mean_finite(x$delta_benchmark_variogram_score_p05),
        mean_benchmark_variogram_score_p2 = jss_joint_mean_finite(x$delta_benchmark_variogram_score_p2),
        mean_benchmark_mean_rmse = jss_joint_mean_finite(x$delta_benchmark_mean_rmse),
        selected_review_cases = nrow(candidates),
        stringsAsFactors = FALSE
      )
    }))
    summary <- jss_joint_bind_rows(list(summary, delta_summary))
  }

  summary[order(summary$case_id, summary$method), , drop = FALSE]
}

jss_joint_failure_summary <- function(results) {
  if (nrow(results) == 0L) return(data.frame())
  eligible <- jss_joint_eligible(results)
  failure_reason <- ifelse(
    eligible,
    "none",
    ifelse(
      !(results$success %in% TRUE) & !is.na(results$error) & nzchar(results$error),
      paste0("optimizer_fit_error: ", results$error),
      ifelse(
        results$success %in% TRUE & !(results$converged %in% TRUE),
        paste0("optimizer_nonconvergence: ", ifelse(is.na(results$stop_reason), "unspecified", results$stop_reason)),
        as.character(results$failure_type)
      )
    )
  )
  x <- data.frame(
    case_id = results$case_id,
    method = results$method,
    failure_reason = failure_reason,
    stringsAsFactors = FALSE
  )
  counts <- stats::aggregate(
    rep(1L, nrow(x)),
    by = x,
    FUN = sum
  )
  names(counts)[names(counts) == "x"] <- "attempts"
  totals <- stats::aggregate(rep(1L, nrow(x)), by = x[c("case_id", "method")], FUN = sum)
  names(totals)[names(totals) == "x"] <- "cell_attempts"
  counts <- merge(counts, totals, by = c("case_id", "method"), all.x = TRUE, sort = FALSE)
  counts$proportion_of_attempts <- counts$attempts / counts$cell_attempts
  counts[order(counts$case_id, counts$method, counts$failure_reason), , drop = FALSE]
}

jss_joint_metric_wins <- function(results) {
  results <- results[jss_joint_eligible(results), , drop = FALSE]
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

jss_joint_mc_stats <- function(x, beneficial = c("positive", "negative"), confidence = 0.95) {
  beneficial <- match.arg(beneficial)
  x <- as.numeric(x)
  x <- x[is.finite(x)]
  n <- length(x)
  estimate <- if (n) mean(x) else NA_real_
  mcse <- if (n > 1L) stats::sd(x) / sqrt(n) else NA_real_
  critical <- if (n > 1L) stats::qt(1 - (1 - confidence) / 2, df = n - 1L) else NA_real_
  wins <- if (!n) logical() else if (identical(beneficial, "positive")) x > 0 else x < 0
  sign_probability <- if (n) mean(wins) else NA_real_
  sign_mcse <- if (n) sqrt(sign_probability * (1 - sign_probability) / n) else NA_real_
  z <- stats::qnorm(1 - (1 - confidence) / 2)
  wilson_denominator <- 1 + z^2 / max(n, 1L)
  wilson_centre <- if (n) (sign_probability + z^2 / (2 * n)) / wilson_denominator else NA_real_
  wilson_half <- if (n) z * sqrt(sign_probability * (1 - sign_probability) / n + z^2 / (4 * n^2)) / wilson_denominator else NA_real_
  data.frame(
    n_finite = n,
    estimate = estimate,
    mcse = mcse,
    conf_low = estimate - critical * mcse,
    conf_high = estimate + critical * mcse,
    sign_probability = sign_probability,
    sign_probability_mcse = sign_mcse,
    sign_probability_conf_low = if (n) max(0, wilson_centre - wilson_half) else NA_real_,
    sign_probability_conf_high = if (n) min(1, wilson_centre + wilson_half) else NA_real_,
    stringsAsFactors = FALSE
  )
}

jss_joint_difference_metric_registry <- function() {
  data.frame(
    metric = c(
      "elapsed_sec", "train_joint_loglik", "train_copula_loglik", "train_marginal_loglik",
      "test_joint_log_score", "test_log_score_per_obs", "test_marginal_log_score",
      "test_copula_log_score", "heldout_variogram_score_p05", "heldout_variogram_score_p2",
      "train_rmse_mu", "train_rmse_sigma", "train_rmse_theta", "train_rmse_tau"
    ),
    beneficial_direction = c(
      "negative", rep("positive", 7L), rep("negative", 6L)
    ),
    evidence_family = c(
      "runtime", rep("likelihood_or_score", 7L), rep("predictive_variogram", 2L), rep("recovery", 4L)
    ),
    failure_penalty = c(60, 100, 100, 100, 100, 10, 100, 100, 10, 10, 10, 10, 10, 10),
    penalty_justification = c(
      "registered one-minute runtime loss", rep("registered 100 log-unit loss", 4L),
      "registered 10 log-score-per-observation loss", rep("registered 100 log-unit loss", 2L),
      rep("registered 10-unit score/error loss", 6L)
    ),
    stringsAsFactors = FALSE
  )
}

jss_joint_failure_inclusive_values <- function(x, metric, beneficial_direction, penalty) {
  joint_retained <- jss_joint_eligible(x, "rs_joint_success", "rs_joint_converged")
  separate_retained <- jss_joint_eligible(x, "rs_separate_success", "rs_separate_converged")
  paired <- joint_retained & separate_retained
  values <- rep(0, nrow(x))
  values[paired] <- x[[paste0("delta_", metric)]][paired]
  joint_only <- joint_retained & !separate_retained
  separate_only <- !joint_retained & separate_retained
  sign <- if (identical(beneficial_direction, "positive")) 1 else -1
  values[joint_only] <- sign * penalty
  values[separate_only] <- -sign * penalty
  values
}

jss_joint_precision_achievement <- function(deltas, precision) {
  if (!nrow(deltas)) return(data.frame())
  required <- as.integer(precision$required_retained_pairs[[1L]])
  confidence <- precision$confidence_level[[1L]]
  target <- precision$target_half_width[[1L]]
  z <- stats::qnorm(1 - (1 - confidence) / 2)
  rows <- lapply(split(deltas, deltas$case_id), function(x) {
    retained <- jss_joint_pair_eligible(x)
    n <- sum(retained)
    wilson_half <- if (n) z * sqrt(0.25 / n + z^2 / (4 * n^2)) / (1 + z^2 / n) else Inf
    data.frame(
      case_id = x$case_id[[1L]], attempted_pairs = nrow(x), retained_pairs = n,
      required_retained_pairs = required, achieved_worst_case_half_width = wilson_half,
      target_half_width = target,
      precision_met = n >= required && is.finite(wilson_half) && wilson_half <= target,
      hard_cap_reached = nrow(x) >= as.integer(precision$max_attempts[[1L]]),
      stringsAsFactors = FALSE
    )
  })
  jss_joint_bind_rows(rows)
}

jss_joint_difference_uncertainty <- function(deltas) {
  if (!nrow(deltas)) return(data.frame())
  registry <- jss_joint_difference_metric_registry()
  registry <- registry[paste0("delta_", registry$metric) %in% names(deltas), , drop = FALSE]
  rows <- list()
  for (case_id in unique(deltas$case_id)) {
    x <- deltas[deltas$case_id == case_id, , drop = FALSE]
    paired_retained <- jss_joint_pair_eligible(x)
    for (i in seq_len(nrow(registry))) {
      metric <- registry$metric[[i]]
      metric_values <- x[[paste0("delta_", metric)]]
      metric_values[!paired_retained] <- NA_real_
      stats <- jss_joint_mc_stats(metric_values, registry$beneficial_direction[[i]])
      joint_retained <- jss_joint_eligible(x, "rs_joint_success", "rs_joint_converged")
      separate_retained <- jss_joint_eligible(x, "rs_separate_success", "rs_separate_converged")
      penalty <- registry$failure_penalty[[i]]
      failure_inclusive <- jss_joint_failure_inclusive_values(
        x, metric, registry$beneficial_direction[[i]], penalty
      )
      sensitivity <- jss_joint_mc_stats(failure_inclusive, registry$beneficial_direction[[i]])
      retention_stats <- jss_joint_mc_stats(
        as.integer(joint_retained) - as.integer(separate_retained), "positive"
      )
      convergence_stats <- jss_joint_mc_stats(
        as.integer(x$rs_joint_converged %in% TRUE) - as.integer(x$rs_separate_converged %in% TRUE),
        "positive"
      )
      rows[[length(rows) + 1L]] <- data.frame(
        case_id = case_id,
        base_case_id = x$base_case_id[[1L]],
        contrast_factor = x$contrast_factor[[1L]],
        metric = metric,
        evidence_family = registry$evidence_family[[i]],
        beneficial_direction = registry$beneficial_direction[[i]],
        attempted_pairs = nrow(x),
        retained_pairs = sum(paired_retained),
        failed_pairs = sum(!paired_retained),
        metric_finite_pairs = stats$n_finite,
        metric_failed_pairs = nrow(x) - stats$n_finite,
        conditional_difference = stats$estimate,
        difference_mcse = stats$mcse,
        difference_conf_low = stats$conf_low,
        difference_conf_high = stats$conf_high,
        joint_better_sign_probability = stats$sign_probability,
        sign_probability_mcse = stats$sign_probability_mcse,
        sign_probability_conf_low = stats$sign_probability_conf_low,
        sign_probability_conf_high = stats$sign_probability_conf_high,
        retention_difference = retention_stats$estimate,
        retention_difference_mcse = retention_stats$mcse,
        retention_difference_conf_low = retention_stats$conf_low,
        retention_difference_conf_high = retention_stats$conf_high,
        convergence_difference = convergence_stats$estimate,
        convergence_difference_mcse = convergence_stats$mcse,
        convergence_difference_conf_low = convergence_stats$conf_low,
        convergence_difference_conf_high = convergence_stats$conf_high,
        failure_penalty = penalty,
        failure_penalty_justification = registry$penalty_justification[[i]],
        failure_inclusive_difference = sensitivity$estimate,
        failure_inclusive_mcse = sensitivity$mcse,
        failure_inclusive_conf_low = sensitivity$conf_low,
        failure_inclusive_conf_high = sensitivity$conf_high,
        failure_inclusive_sign_probability = sensitivity$sign_probability,
        failure_inclusive_sign_probability_mcse = sensitivity$sign_probability_mcse,
        failure_inclusive_sign_probability_conf_low = sensitivity$sign_probability_conf_low,
        failure_inclusive_sign_probability_conf_high = sensitivity$sign_probability_conf_high,
        conditioning = "both fits returned and converged",
        stringsAsFactors = FALSE
      )
    }
  }
  jss_joint_bind_rows(rows)
}

jss_joint_hypothesis_evidence <- function(hypotheses, uncertainty) {
  if (!nrow(hypotheses) || !nrow(uncertainty)) return(data.frame())
  uncertainty_fields <- c(
    "attempted_pairs", "retained_pairs", "failed_pairs", "metric_finite_pairs", "metric_failed_pairs",
    "conditional_difference",
    "difference_mcse", "difference_conf_low", "difference_conf_high",
    "joint_better_sign_probability", "sign_probability_mcse",
    "sign_probability_conf_low", "sign_probability_conf_high",
    "retention_difference", "retention_difference_mcse", "retention_difference_conf_low",
    "retention_difference_conf_high", "convergence_difference", "convergence_difference_mcse",
    "convergence_difference_conf_low", "convergence_difference_conf_high",
    "failure_penalty", "failure_penalty_justification", "failure_inclusive_difference",
    "failure_inclusive_mcse", "failure_inclusive_conf_low", "failure_inclusive_conf_high",
    "failure_inclusive_sign_probability", "failure_inclusive_sign_probability_mcse",
    "failure_inclusive_sign_probability_conf_low", "failure_inclusive_sign_probability_conf_high",
    "conditioning"
  )
  hypothesis_fields <- c(
    "hypothesis", "focal_case", "comparator_case", "decision", "rationale", "metric",
    "contrast_n_finite", "contrast_estimate", "contrast_mcse", "contrast_conf_low",
    "contrast_conf_high", "interval_excludes_zero", "failure_inclusive_contrast",
    "failure_inclusive_contrast_mcse", "failure_inclusive_contrast_conf_low",
    "failure_inclusive_contrast_conf_high", "retention_contrast", "retention_contrast_mcse",
    "retention_contrast_conf_low", "retention_contrast_conf_high", "convergence_contrast",
    "convergence_contrast_mcse", "convergence_contrast_conf_low", "convergence_contrast_conf_high",
    "publication_suitable", "decision_zero_tolerance"
  )
  focal <- merge(
    hypotheses[, intersect(hypothesis_fields, names(hypotheses)), drop = FALSE],
    uncertainty,
    by.x = c("focal_case", "metric"),
    by.y = c("case_id", "metric"),
    all.x = TRUE,
    sort = FALSE
  )
  focal_present <- intersect(uncertainty_fields, names(focal))
  names(focal)[match(focal_present, names(focal))] <- paste0("focal_", focal_present)
  comparator <- uncertainty[, c("case_id", "metric", intersect(uncertainty_fields, names(uncertainty))), drop = FALSE]
  comparator_present <- intersect(uncertainty_fields, names(comparator))
  names(comparator)[match(comparator_present, names(comparator))] <- paste0("comparator_", comparator_present)
  merge(
    focal,
    comparator,
    by.x = c("comparator_case", "metric"),
    by.y = c("case_id", "metric"),
    all.x = TRUE,
    sort = FALSE
  )
}

jss_joint_case_delta_summary <- function(deltas) {
  if (nrow(deltas) == 0L) {
    return(data.frame())
  }
  split_key <- interaction(deltas$case_id, drop = TRUE, lex.order = TRUE)
  jss_joint_bind_rows(lapply(split(seq_len(nrow(deltas)), split_key), function(idx) {
    x <- deltas[idx, , drop = FALSE]
    joint_retained <- jss_joint_eligible(x, "rs_joint_success", "rs_joint_converged")
    separate_retained <- jss_joint_eligible(x, "rs_separate_success", "rs_separate_converged")
    paired_retained <- joint_retained & separate_retained
    delta_fields <- grep("^delta_", names(x), value = TRUE)
    x[!paired_retained, delta_fields] <- NA_real_
    dep_delta <- pmin(
      jss_joint_numeric(x$delta_train_rmse_theta),
      jss_joint_numeric(x$delta_train_rmse_tau),
      na.rm = TRUE
    )
    dep_delta[!is.finite(dep_delta)] <- NA_real_
    log_score_mc <- jss_joint_mc_stats(x$delta_test_log_score_per_obs, "positive")
    dependence_mc <- jss_joint_mc_stats(dep_delta, "negative")
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
      attempts = nrow(x),
      paired_successes = sum(paired_retained),
      paired_failures = sum(!paired_retained),
      paired_success_rate = mean(paired_retained),
      failure_inclusive_success_advantage = mean(
        as.integer(joint_retained) - as.integer(separate_retained)
      ),
      delta_train_joint_loglik_mean = jss_joint_mean_finite(x$delta_train_joint_loglik),
      delta_test_log_score_per_obs_mean = jss_joint_mean_finite(x$delta_test_log_score_per_obs),
      delta_heldout_variogram_score_p05_mean = jss_joint_mean_finite(x$delta_heldout_variogram_score_p05),
      delta_heldout_variogram_score_p2_mean = jss_joint_mean_finite(x$delta_heldout_variogram_score_p2),
      delta_train_rmse_theta_mean = jss_joint_mean_finite(x$delta_train_rmse_theta),
      delta_train_rmse_tau_mean = jss_joint_mean_finite(x$delta_train_rmse_tau),
      delta_dependence_rmse_mean = jss_joint_mean_finite(dep_delta),
      test_log_score_n_finite = log_score_mc$n_finite,
      test_log_score_mcse = log_score_mc$mcse,
      test_log_score_conf_low = log_score_mc$conf_low,
      test_log_score_conf_high = log_score_mc$conf_high,
      test_log_score_joint_better_probability = log_score_mc$sign_probability,
      dependence_rmse_n_finite = dependence_mc$n_finite,
      dependence_rmse_mcse = dependence_mc$mcse,
      dependence_rmse_conf_low = dependence_mc$conf_low,
      dependence_rmse_conf_high = dependence_mc$conf_high,
      dependence_rmse_joint_better_probability = dependence_mc$sign_probability,
      stringsAsFactors = FALSE
    )
  }))
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

jss_joint_interval_direction <- function(low, high) {
  if (!is.finite(low) || !is.finite(high) || low <= 0 && high >= 0) return(0L)
  if (low > 0) 1L else -1L
}

jss_joint_hypothesis_summary <- function(deltas) {
  case_summary <- jss_joint_case_delta_summary(deltas)
  if (nrow(case_summary) == 0L) {
    return(data.frame())
  }
  contrasts <- data.frame(
    focal_case = c("JVS02", "JVS03", "JVS04"),
    hypothesis = c(
      "Dependence strength changes the RS joint-versus-separate tradeoff",
      "Scale-signal strength changes the RS joint-versus-separate tradeoff",
      "A count margin changes the RS joint-versus-separate tradeoff"
    ),
    rationale = c(
      "JVS02 changes only theta_strength from the JVS01 base design.",
      "JVS03 changes only sigma_strength from the JVS01 base design.",
      "JVS04 changes only family from the JVS01 base design."
    ),
    stringsAsFactors = FALSE
  )
  key_metrics <- c("train_joint_loglik", "heldout_variogram_score_p05", "heldout_variogram_score_p2")
  directions <- c("positive", "negative", "negative")
  zero_tolerance <- 0
  jss_joint_bind_rows(lapply(seq_len(nrow(contrasts)), function(i) {
    focal <- deltas[deltas$case_id == contrasts$focal_case[[i]], , drop = FALSE]
    base <- deltas[deltas$case_id == "JVS01", , drop = FALSE]
    paired <- merge(
      focal[, c("joint_review_rep", paste0("delta_", key_metrics)), drop = FALSE],
      base[, c("joint_review_rep", paste0("delta_", key_metrics)), drop = FALSE],
      by = "joint_review_rep", suffixes = c("_focal", "_base"), all = FALSE
    )
    metric_rows <- lapply(seq_along(key_metrics), function(j) {
      nm <- key_metrics[[j]]
      registry <- jss_joint_difference_metric_registry()
      meta <- registry[registry$metric == nm, , drop = FALSE]
      contrast <- paired[[paste0("delta_", nm, "_focal")]] -
        paired[[paste0("delta_", nm, "_base")]]
      stat <- jss_joint_mc_stats(contrast, directions[[j]])
      excludes_zero <- is.finite(stat$conf_low) && is.finite(stat$conf_high) &&
        (stat$conf_low > zero_tolerance || stat$conf_high < -zero_tolerance)
      matched_rep <- intersect(focal$joint_review_rep, base$joint_review_rep)
      focal_order <- match(matched_rep, focal$joint_review_rep)
      base_order <- match(matched_rep, base$joint_review_rep)
      focal_fi <- jss_joint_failure_inclusive_values(
        focal, nm, meta$beneficial_direction[[1L]], meta$failure_penalty[[1L]]
      )[focal_order]
      base_fi <- jss_joint_failure_inclusive_values(
        base, nm, meta$beneficial_direction[[1L]], meta$failure_penalty[[1L]]
      )[base_order]
      failure_stat <- jss_joint_mc_stats(focal_fi - base_fi, directions[[j]])
      focal_retention <- as.integer(jss_joint_eligible(focal, "rs_joint_success", "rs_joint_converged")) -
        as.integer(jss_joint_eligible(focal, "rs_separate_success", "rs_separate_converged"))
      base_retention <- as.integer(jss_joint_eligible(base, "rs_joint_success", "rs_joint_converged")) -
        as.integer(jss_joint_eligible(base, "rs_separate_success", "rs_separate_converged"))
      retention_stat <- jss_joint_mc_stats(focal_retention[focal_order] - base_retention[base_order], "positive")
      focal_convergence <- as.integer(focal$rs_joint_converged %in% TRUE) -
        as.integer(focal$rs_separate_converged %in% TRUE)
      base_convergence <- as.integer(base$rs_joint_converged %in% TRUE) -
        as.integer(base$rs_separate_converged %in% TRUE)
      convergence_stat <- jss_joint_mc_stats(
        focal_convergence[focal_order] - base_convergence[base_order], "positive"
      )
      conditional_sign <- jss_joint_interval_direction(stat$conf_low, stat$conf_high)
      failure_sign <- jss_joint_interval_direction(failure_stat$conf_low, failure_stat$conf_high)
      retention_neutral <- is.finite(retention_stat$conf_low) && is.finite(retention_stat$conf_high) &&
        retention_stat$conf_low <= zero_tolerance && retention_stat$conf_high >= -zero_tolerance
      convergence_neutral <- is.finite(convergence_stat$conf_low) && is.finite(convergence_stat$conf_high) &&
        convergence_stat$conf_low <= zero_tolerance && convergence_stat$conf_high >= -zero_tolerance
      publication_suitable <- conditional_sign != 0L && conditional_sign == failure_sign &&
        retention_neutral && convergence_neutral
      data.frame(metric = nm, contrast_n_finite = stat$n_finite,
        contrast_estimate = stat$estimate, contrast_mcse = stat$mcse,
        contrast_conf_low = stat$conf_low, contrast_conf_high = stat$conf_high,
        interval_excludes_zero = excludes_zero,
        failure_inclusive_contrast = failure_stat$estimate,
        failure_inclusive_contrast_mcse = failure_stat$mcse,
        failure_inclusive_contrast_conf_low = failure_stat$conf_low,
        failure_inclusive_contrast_conf_high = failure_stat$conf_high,
        retention_contrast = retention_stat$estimate,
        retention_contrast_mcse = retention_stat$mcse,
        retention_contrast_conf_low = retention_stat$conf_low,
        retention_contrast_conf_high = retention_stat$conf_high,
        convergence_contrast = convergence_stat$estimate,
        convergence_contrast_mcse = convergence_stat$mcse,
        convergence_contrast_conf_low = convergence_stat$conf_low,
        convergence_contrast_conf_high = convergence_stat$conf_high,
        publication_suitable = publication_suitable,
        decision_zero_tolerance = zero_tolerance,
        decision = if (publication_suitable) "supported_all_registered_intervals" else
          "inconclusive_or_failure_sensitive",
        stringsAsFactors = FALSE)
    })
    metric_rows <- jss_joint_bind_rows(metric_rows)
    base_row <- jss_joint_hypothesis_row(
      case_summary, contrasts$hypothesis[[i]], contrasts$focal_case[[i]], "JVS01",
      "metric_specific_interval_decision",
      paste0(contrasts$rationale[[i]],
        " Publication requires coherent conditional and failure-inclusive intervals and neutral retention/convergence contrasts; point signs cannot drive claims.")
    )
    out <- merge(base_row, metric_rows, by = NULL, suffixes = c("_hypothesis", ""))
    out$decision_hypothesis <- NULL
    out
  }))
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
  cfg <- jss_joint_simulation_settings(settings)
  precision <- cfg$precision
  for (nm in names(precision)) candidates[[nm]] <- precision[[nm]][[1L]]
  utils::write.csv(candidates, paths$candidate_selection, row.names = FALSE)
  utils::write.csv(precision, paths$precision_registry, row.names = FALSE)

  results <- jss_joint_run_simulation(settings, candidates, cfg = cfg)
  checkpoint_status <- attr(results, "checkpoint_status")
  results$profile <- settings$profile
  utils::write.csv(results, paths$results, row.names = FALSE)
  utils::write.csv(checkpoint_status, paths$checkpoint_status, row.names = FALSE)
  archive_groups <- split(results,
    interaction(results$case_id, results$joint_review_rep, drop = TRUE, lex.order = TRUE))
  archive_groups <- lapply(archive_groups, function(x)
    x[, setdiff(names(x), c("profile", "result_portable_sha256")), drop = FALSE])
  saveRDS(archive_groups, paths$checkpoint_payload_archive, version = 3L)
  archive_manifest <- data.frame(
    case_id = vapply(archive_groups, function(x) x$case_id[[1L]], character(1)),
    joint_review_rep = vapply(archive_groups, function(x) x$joint_review_rep[[1L]], integer(1)),
    result_content_sha256 = vapply(archive_groups, jss_joint_content_sha256, character(1)),
    result_portable_sha256 = vapply(archive_groups, jss_joint_portable_result_sha256, character(1)),
    stringsAsFactors = FALSE)
  utils::write.csv(archive_manifest, paths$checkpoint_content_manifest, row.names = FALSE)

  # Publication summaries are derived from the exact CSV representation that
  # downstream validators and readers consume.
  results <- jss_joint_read_canonical_results_csv(paths$results)

  deltas <- jss_joint_delta_table(results)
  precision_achievement <- jss_joint_precision_achievement(deltas, precision)
  precision$achieved_min_retained_pairs <- min(precision_achievement$retained_pairs)
  precision$achieved_max_half_width <- max(precision_achievement$achieved_worst_case_half_width)
  precision$all_cells_precision_met <- all(precision_achievement$precision_met)
  utils::write.csv(precision, paths$precision_registry, row.names = FALSE)
  if (isTRUE(precision$production_requires_precision[[1L]]) &&
      !isTRUE(precision$all_cells_precision_met[[1L]])) {
    stop(
      "Optimizer production evidence is ineligible: at least one case failed the registered retained-pair precision gate at the hard cap.",
      call. = FALSE
    )
  }
  summary <- jss_joint_summary_table(results, deltas, candidates)
  metric_wins <- jss_joint_metric_wins(results)
  failure_summary <- jss_joint_failure_summary(results)
  difference_uncertainty <- jss_joint_difference_uncertainty(deltas)
  hypothesis_summary <- jss_joint_hypothesis_summary(deltas)
  hypothesis_evidence <- jss_joint_hypothesis_evidence(hypothesis_summary, difference_uncertainty)
  utils::write.csv(summary, paths$summary, row.names = FALSE)
  utils::write.csv(metric_wins, paths$metric_wins, row.names = FALSE)
  utils::write.csv(hypothesis_summary, paths$hypothesis_summary, row.names = FALSE)
  utils::write.csv(failure_summary, paths$failure_summary, row.names = FALSE)
  utils::write.csv(difference_uncertainty, paths$difference_uncertainty, row.names = FALSE)
  utils::write.csv(hypothesis_evidence, paths$hypothesis_evidence, row.names = FALSE)

  jss_joint_write_delta_figure(deltas, paths$deltas_figure)
  jss_joint_write_metric_dashboard(metric_wins, paths$metric_dashboard)
  figure_registry <- data.frame(
    figure = basename(c(paths$deltas_figure, paths$metric_dashboard)),
    png_sha256 = vapply(c(paths$deltas_figure, paths$metric_dashboard), jss_joint_sha256_file, character(1)),
    plotted_data_sha256 = c(jss_joint_portable_frame_sha256(deltas),
      jss_joint_portable_frame_sha256(metric_wins)),
    plot_spec_sha256 = c(
      digest::digest("delta density+zero line+case facet; width=10 height=7 dpi=320", "sha256", serialize = FALSE),
      digest::digest("metric win-or-tie bars; width=9 height=6 dpi=320", "sha256", serialize = FALSE)
    ),
    stringsAsFactors = FALSE
  )
  utils::write.csv(figure_registry, paths$figure_registry, row.names = FALSE)

  list(
    module_id = jss_joint_module_id(),
    title = "Joint versus separate optimisation",
    status = "current",
    data = c(paths$candidate_selection, paths$precision_registry, paths$checkpoint_status,
      paths$results, paths$figure_registry, paths$checkpoint_payload_archive,
      paths$checkpoint_content_manifest),
    tables = c(
      paths$summary, paths$metric_wins, paths$hypothesis_summary,
      paths$failure_summary, paths$difference_uncertainty, paths$hypothesis_evidence
    ),
    figures = c(paths$deltas_figure, paths$metric_dashboard),
    notes = paste(
      "One registered base design and three paired one-factor contrasts compare",
      "RS-separate with RS-joint. CG is excluded because no retained empirical",
      "guidance currently requires it. Attempts are derived from the registered",
      "Monte Carlo sign-probability precision target."
    )
  )
}
