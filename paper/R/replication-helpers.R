jss_settings <- function(create = TRUE) {
  profile <- Sys.getenv("GAMLSS_LONGITUDINAL_JSS_PROFILE", unset = "smoke")
  if (identical(profile, "expanded")) profile <- "paper"
  profile <- match.arg(profile, c("smoke", "paper", "full"))
  root <- Sys.getenv("GAMLSS_LONGITUDINAL_JSS_ROOT", unset = getwd())
  out_dir <- file.path(root, "results", "jss-replication", profile)
  dirs <- file.path(out_dir, c("data", "tables", "figures", "logs"))
  if (isTRUE(create)) invisible(lapply(dirs, dir.create, recursive = TRUE, showWarnings = FALSE))
  list(
    profile = profile,
    root = root,
    out_dir = out_dir,
    data_dir = file.path(out_dir, "data"),
    tables_dir = file.path(out_dir, "tables"),
    figures_dir = file.path(out_dir, "figures"),
    logs_dir = file.path(out_dir, "logs"),
    seed = 20260528L,
    workers = as.integer(Sys.getenv("GAMLSS_LONGITUDINAL_JSS_WORKERS", "1")),
    public_data_dir = file.path(root, "paper", "data", "public-derived")
  )
}

jss_module_paths <- function(settings, module_id) {
  list(
    data = file.path(settings$data_dir, paste0(module_id, "-data.csv")),
    table = file.path(settings$tables_dir, paste0(module_id, "-status.csv")),
    figure = file.path(settings$figures_dir, paste0(module_id, "-stub.png"))
  )
}

jss_external_data_status <- function(envvar = NULL) {
  if (is.null(envvar) || is.na(envvar) || !nzchar(envvar)) {
    return(list(envvar = NA_character_, status = "not_required"))
  }
  value <- Sys.getenv(envvar, unset = "")
  if (!nzchar(value)) {
    return(list(envvar = envvar, status = "not_set"))
  }
  list(
    envvar = envvar,
    status = if (file.exists(value)) "available" else "path_not_found"
  )
}

jss_write_stub_module <- function(settings, module_id, title, family, copula,
                                  focus, planned_outputs,
                                  external_envvar = NULL,
                                  external_data_note = "not required") {
  paths <- jss_module_paths(settings, module_id)
  external <- jss_external_data_status(external_envvar)
  data <- data.frame(
    module_id = module_id,
    title = title,
    profile = settings$profile,
    analysis_state = "stub",
    family = family,
    copula = copula,
    focus = focus,
    external_data_envvar = external$envvar,
    external_data_status = external$status,
    external_data_note = external_data_note,
    stringsAsFactors = FALSE
  )
  utils::write.csv(data, paths$data, row.names = FALSE)

  output_rows <- data.frame(
    module_id = module_id,
    artifact = c("data", "table", "figure"),
    path = unlist(paths, use.names = FALSE),
    status = "stub",
    planned_outputs = planned_outputs,
    stringsAsFactors = FALSE
  )
  utils::write.csv(output_rows, paths$table, row.names = FALSE)

  plot_data <- data.frame(
    stage = c("pipeline", "analysis", "paper output"),
    complete = c(1, 0, 0)
  )
  p <- ggplot2::ggplot(plot_data, ggplot2::aes(x = stage, y = complete, fill = stage)) +
    ggplot2::geom_col(width = 0.7, show.legend = FALSE) +
    ggplot2::coord_cartesian(ylim = c(0, 1)) +
    ggplot2::scale_y_continuous(breaks = c(0, 1), labels = c("pending", "stub ready")) +
    ggplot2::labs(
      title = paste("STUB:", title),
      subtitle = "Pipeline artifact only; final JSS analysis is pending.",
      x = NULL,
      y = NULL
    ) +
    ggplot2::theme_minimal(base_size = 11)
  ggplot2::ggsave(paths$figure, p, width = 7, height = 4.5, dpi = 220, bg = "white")

  list(
    module_id = module_id,
    title = title,
    status = "stub",
    data = paths$data,
    tables = paths$table,
    figures = paths$figure,
    notes = planned_outputs
  )
}

jss_copy_final_artifacts <- function(settings, module_id, title, source_dir, artifacts, notes) {
  required <- c("source_file", "output_file", "artifact_type", "role")
  missing_cols <- setdiff(required, names(artifacts))
  if (length(missing_cols) > 0L) {
    stop(
      "Artifact specification for ", module_id, " is missing column(s): ",
      paste(missing_cols, collapse = ", "),
      call. = FALSE
    )
  }

  if (!dir.exists(source_dir)) {
    stop("Missing final artifact source directory: ", source_dir, call. = FALSE)
  }

  dest_base <- stats::setNames(
    c(settings$data_dir, settings$tables_dir, settings$figures_dir),
    c("data", "table", "figure")
  )
  bad_types <- setdiff(unique(artifacts$artifact_type), names(dest_base))
  if (length(bad_types) > 0L) {
    stop(
      "Unknown artifact type(s) for ", module_id, ": ",
      paste(bad_types, collapse = ", "),
      call. = FALSE
    )
  }

  artifacts$module_id <- module_id
  artifacts$title <- title
  artifacts$status <- "current"
  artifacts$source_path <- file.path(source_dir, artifacts$source_file)
  artifacts$output_path <- file.path(
    unname(dest_base[artifacts$artifact_type]),
    artifacts$output_file
  )

  missing_sources <- artifacts$source_path[!file.exists(artifacts$source_path)]
  if (length(missing_sources) > 0L) {
    stop(
      "Missing final artifact(s) for ", module_id, ": ",
      paste(missing_sources, collapse = ", "),
      call. = FALSE
    )
  }

  invisible(lapply(unique(dirname(artifacts$output_path)), dir.create, recursive = TRUE, showWarnings = FALSE))
  copied <- file.copy(artifacts$source_path, artifacts$output_path, overwrite = TRUE)
  if (any(!copied)) {
    stop(
      "Failed to copy final artifact(s) for ", module_id, ": ",
      paste(artifacts$source_path[!copied], collapse = ", "),
      call. = FALSE
    )
  }

  index_path <- file.path(settings$data_dir, paste0(module_id, "-source-artifacts.csv"))
  index <- artifacts[
    ,
    c("module_id", "title", "status", "role", "artifact_type", "source_path", "output_path"),
    drop = FALSE
  ]
  index$source_path <- normalizePath(index$source_path, winslash = "/", mustWork = TRUE)
  index$output_path <- normalizePath(index$output_path, winslash = "/", mustWork = TRUE)
  utils::write.csv(index, index_path, row.names = FALSE)

  data_files <- c(index_path, artifacts$output_path[artifacts$artifact_type == "data"])
  table_files <- artifacts$output_path[artifacts$artifact_type == "table"]
  figure_files <- artifacts$output_path[artifacts$artifact_type == "figure"]

  list(
    module_id = module_id,
    title = title,
    status = "current",
    data = unname(data_files),
    tables = unname(table_files),
    figures = unname(figure_files),
    notes = notes
  )
}

jss_collect_module_files <- function(...) {
  modules <- list(...)
  files <- unlist(lapply(modules, function(x) c(x$data, x$tables, x$figures)), use.names = FALSE)
  unique(files[file.exists(files)])
}

jss_methods <- function(profile) {
  methods <- c("gamlss", "rs_separate", "rs_joint", "cg")
  if (profile %in% c("paper", "full") && requireNamespace("gamlss2", quietly = TRUE)) {
    methods <- c("gamlss2", methods)
  }
  methods
}

jss_grid_settings <- function(profile) {
  if (profile %in% c("paper", "full")) {
    return(list(
      families = c("NO", "GA", "BCPE", "PO", "NBI", "DEL", "ZIP", "ZINBI"),
      copulas = c("N", "C", "F", "G", "J", "t"),
      designs = c("intercept", "covariate", "smooth"),
      n = 120L,
      times = 1:4,
      max_outer_iter = 12L,
      max_inner_iter = 12L,
      max_elapsed_sec = 90,
      dependence = "moderate",
      missingness = "none"
    ))
  }

  list(
    families = c("NO", "NBI"),
    copulas = c("N", "C"),
    designs = c("intercept", "covariate"),
    n = 36L,
    times = 1:3,
    max_outer_iter = 3L,
    max_inner_iter = 3L,
    max_elapsed_sec = 20,
    dependence = "moderate",
    missingness = "none"
  )
}

jss_run_coverage_suite <- function(settings) {
  set.seed(settings$seed)
  grid <- jss_grid_settings(settings$profile)
  results <- withCallingHandlers(gamlss.longitudinal::run_coverage_simulations(
    families = grid$families,
    copulas = grid$copulas,
    methods = jss_methods(settings$profile),
    designs = grid$designs,
    n = grid$n,
    times = grid$times,
    seed = settings$seed,
    max_outer_iter = grid$max_outer_iter,
    max_inner_iter = grid$max_inner_iter,
    max_elapsed_sec = grid$max_elapsed_sec,
    dependence = grid$dependence,
    missingness = grid$missingness,
    output_dir = settings$data_dir,
    write_results = TRUE
  ), warning = function(w) {
    if (grepl("max_outer_iter|not converged", conditionMessage(w), ignore.case = TRUE)) invokeRestart("muffleWarning")
  })
  nonconverged <- !is.na(results$success) & results$success &
    (is.na(results$converged) | !results$converged)
  results$failure_type[nonconverged & results$failure_type %in% c("ok", "none", "")] <-
    "optimizer_nonconvergence"
  saveRDS(results, file.path(settings$data_dir, "coverage_results.rds"))
  results
}

jss_write_coverage_summary <- function(results, settings) {
  group_cols <- c("family", "copula", "design", "method", "success", "converged", "failure_type")
  metric_cols <- c(
    "marginal_loglik",
    "joint_loglik",
    "elapsed_sec",
    "max_abs_param_error",
    "max_rel_param_error",
    "margin_gap_pct_vs_reference",
    "joint_delta_pct_vs_rs_separate"
  )
  split_key <- interaction(results[group_cols], drop = TRUE, lex.order = TRUE)
  out <- do.call(rbind, lapply(split(seq_len(nrow(results)), split_key), function(idx) {
    row <- results[idx[1], group_cols, drop = FALSE]
    for (metric in metric_cols) {
      row[[metric]] <- mean(results[[metric]][idx], na.rm = TRUE)
    }
    row
  }))
  rownames(out) <- NULL
  path <- file.path(settings$tables_dir, "coverage_summary.csv")
  utils::write.csv(out, path, row.names = FALSE)
  path
}

jss_write_fit_event_audit <- function(results, settings) {
  event <- is.na(results$success) | !results$success | is.na(results$converged) | !results$converged
  keep <- intersect(
    c("family", "copula", "design", "method", "seed", "success", "converged",
      "failure_type", "warnings", "error", "elapsed_sec"),
    names(results)
  )
  audit <- results[event, keep, drop = FALSE]
  audit$event_type <- ifelse(
    is.na(audit$success) | !audit$success,
    "execution_failure",
    "optimizer_nonconvergence"
  )
  audit <- audit[c("event_type", setdiff(names(audit), "event_type"))]
  path <- file.path(settings$logs_dir, "fit-events.csv")
  utils::write.csv(audit, path, row.names = FALSE)
  path
}

jss_write_convergence_summary <- function(results, settings) {
  out <- as.data.frame(
    xtabs(~ family + copula + design + method + success + converged + failure_type, data = results),
    stringsAsFactors = FALSE
  )
  names(out)[names(out) == "Freq"] <- "n"
  out <- out[out$n > 0, , drop = FALSE]
  path <- file.path(settings$tables_dir, "convergence_summary.csv")
  utils::write.csv(out, path, row.names = FALSE)
  path
}

jss_write_runtime_figure <- function(results, settings) {
  path <- file.path(settings$figures_dir, "runtime_by_method.png")
  p <- ggplot2::ggplot(results, ggplot2::aes(x = method, y = elapsed_sec, fill = method)) +
    ggplot2::geom_boxplot(width = 0.65, outlier.alpha = 0.45, na.rm = TRUE) +
    ggplot2::facet_grid(family ~ copula, scales = "free_y") +
    ggplot2::labs(x = NULL, y = "Runtime (s)") +
    ggplot2::theme_minimal(base_size = 10) +
    ggplot2::theme(legend.position = "none", axis.text.x = ggplot2::element_text(angle = 35, hjust = 1))
  ggplot2::ggsave(path, p, width = 10, height = 7, dpi = 320, bg = "white")
  path
}

jss_write_convergence_figure <- function(results, settings) {
  path <- file.path(settings$figures_dir, "convergence_by_method.png")
  results$converged_flag <- !is.na(results$converged) & results$converged
  conv <- stats::aggregate(converged_flag ~ family + copula + design + method, data = results, FUN = mean)
  names(conv)[names(conv) == "converged_flag"] <- "converged"
  p <- ggplot2::ggplot(conv, ggplot2::aes(x = method, y = converged, fill = method)) +
    ggplot2::geom_col(width = 0.7) +
    ggplot2::facet_grid(family ~ copula) +
    ggplot2::scale_y_continuous(labels = scales::percent_format(accuracy = 1)) +
    ggplot2::coord_cartesian(ylim = c(0, 1)) +
    ggplot2::labs(x = NULL, y = "Convergence rate") +
    ggplot2::theme_minimal(base_size = 10) +
    ggplot2::theme(legend.position = "none", axis.text.x = ggplot2::element_text(angle = 35, hjust = 1))
  ggplot2::ggsave(path, p, width = 10, height = 7, dpi = 320, bg = "white")
  path
}

jss_write_session_info <- function(settings) {
  path <- file.path(settings$logs_dir, "session_info.txt")
  lines <- c(
    paste("profile:", settings$profile),
    paste("seed:", settings$seed),
    paste("package_version:", as.character(utils::packageVersion("gamlss.longitudinal"))),
    paste("git_sha:", jss_git_sha(settings$root)),
    "",
    capture.output(utils::sessionInfo())
  )
  writeLines(lines, path, useBytes = TRUE)
  path
}

jss_git_sha <- function(root) {
  git_root <- if (.Platform$OS.type == "windows") shQuote(root, type = "cmd") else root
  out <- tryCatch(
    suppressWarnings(system2("git", c("-C", git_root, "rev-parse", "HEAD"), stdout = TRUE, stderr = TRUE)),
    error = function(e) NA_character_
  )
  status <- attr(out, "status")
  if (length(out) == 0 || isTRUE(!is.na(status) && status != 0) || is.na(out[1])) "unknown" else out[1]
}

jss_write_manifest <- function(settings) {
  template <- utils::read.csv(file.path(settings$root, "paper", "manifest.csv"), stringsAsFactors = FALSE)
  active <- vapply(strsplit(template$profiles, "|", fixed = TRUE), function(x) settings$profile %in% x, logical(1))
  template <- template[active & template$access == "public" & template$publication_status == "active", , drop = FALSE]
  template$output_path <- gsub("\\\\", "/", template$generated_path)
  if (identical(settings$profile, "full")) {
    template$verification[template$access == "public" & template$verification == "canonical_sha256"] <- "documented_numeric_tolerance"
  }
  template$generated_at <- format(Sys.time(), "%Y-%m-%dT%H:%M:%S%z")
  template$profile_run <- settings$profile
  path <- file.path(settings$out_dir, "manifest.csv")
  utils::write.csv(template, path, row.names = FALSE)
  path
}

jss_manifest_output_files <- function(manifest, manifest_path) {
  paths <- manifest$output_path
  absolute <- grepl("^([A-Za-z]:[/\\\\]|/)", paths)
  paths[!absolute & nzchar(paths)] <- file.path(dirname(manifest_path), paths[!absolute & nzchar(paths)])
  paths
}

jss_write_output_hashes <- function(settings, manifest_path, ...) {
  manifest <- utils::read.csv(manifest_path, stringsAsFactors = FALSE)
  manifest_files <- jss_manifest_output_files(manifest, manifest_path)
  files <- c(manifest_files[nzchar(manifest_files)], manifest_path,
    list.files(settings$data_dir, full.names = TRUE, recursive = TRUE),
    list.files(settings$tables_dir, full.names = TRUE, recursive = TRUE),
    list.files(settings$figures_dir, full.names = TRUE, recursive = TRUE))
  files <- unique(files)
  files <- files[file.exists(files)]
  root_prefix <- paste0(normalizePath(settings$out_dir, winslash = "/", mustWork = TRUE), "/")
  normalized <- normalizePath(files, winslash = "/", mustWork = TRUE)
  hashes <- data.frame(
    file = ifelse(startsWith(normalized, root_prefix), substring(normalized, nchar(root_prefix) + 1L), normalized),
    sha256 = vapply(files, jss_file_sha256, character(1)),
    bytes = unname(file.info(files)$size),
    stringsAsFactors = FALSE
  )
  path <- file.path(settings$logs_dir, "output_hashes.csv")
  utils::write.csv(hashes, path, row.names = FALSE)
  path
}

jss_validate_manifest <- function(manifest_path, output_hashes_path) {
  manifest <- utils::read.csv(manifest_path, stringsAsFactors = FALSE)
  public <- manifest$access == "public" & manifest$publication_status == "active"
  if (any(!nzchar(manifest$producer[public]))) stop("Unclassified or producer-free public artifact in manifest.", call. = FALSE)
  resolved <- jss_manifest_output_files(manifest, manifest_path)
  expected <- resolved[public]
  missing <- expected[!file.exists(expected)]
  if (length(missing) > 0) {
    stop("Missing JSS replication output(s): ", paste(missing, collapse = ", "), call. = FALSE)
  }
  if (!file.exists(output_hashes_path)) {
    stop("Missing output hash log: ", output_hashes_path, call. = FALSE)
  }
  too_small <- expected[file.info(expected)$size < 100L]
  if (length(too_small)) stop("Blank or truncated output(s): ", paste(too_small, collapse = ", "), call. = FALSE)
  exact <- public & manifest$verification == "canonical_sha256" & nzchar(manifest$approved_sha256)
  actual <- vapply(resolved[exact], jss_file_sha256, character(1))
  mismatch <- manifest$manuscript_path[exact][actual != manifest$approved_sha256[exact]]
  if (length(mismatch)) stop("Canonical table hash mismatch: ", paste(mismatch, collapse = ", "), call. = FALSE)
  figure_rows <- public & manifest$artifact_type == "figure"
  figures <- resolved[figure_rows]
  valid_figures <- vapply(figures, jss_png_nonblank, logical(1))
  actual_figure_hash <- vapply(figures, jss_file_sha256, character(1))
  reference_hash <- manifest$approved_sha256[figure_rows]
  figure_audit <- data.frame(
    artifact_id = manifest$artifact_id[figure_rows],
    output_path = manifest$output_path[figure_rows],
    policy = manifest$verification[figure_rows],
    dimensions_nonblank = valid_figures,
    actual_sha256 = actual_figure_hash,
    approved_reference_sha256 = reference_hash,
    matches_approved_reference = nzchar(reference_hash) & actual_figure_hash == reference_hash,
    stringsAsFactors = FALSE
  )
  utils::write.csv(figure_audit, file.path(dirname(manifest_path), "logs", "figure-verification.csv"), row.names = FALSE)
  invalid_figures <- figures[!valid_figures]
  if (length(invalid_figures)) stop("Invalid or blank PNG output(s): ", paste(invalid_figures, collapse = ", "), call. = FALSE)
  TRUE
}

jss_metric_numeric <- function(x) {
  if (is.logical(x)) return(as.numeric(x))
  if (is.numeric(x)) return(as.numeric(x))
  value <- trimws(as.character(x))
  value[value %in% c("", "--", "NA", "NaN")] <- NA_character_
  value[tolower(value) == "true"] <- "1"
  value[tolower(value) == "false"] <- "0"
  suppressWarnings(as.numeric(value))
}

jss_metric_key <- function(x, columns) {
  if (!length(columns)) return(rep("all", nrow(x)))
  missing <- setdiff(columns, names(x))
  if (length(missing)) stop("Tolerance input is missing key column(s): ", paste(missing, collapse = ", "), call. = FALSE)
  pieces <- lapply(columns, function(column) paste0(column, "=", as.character(x[[column]])))
  do.call(paste, c(pieces, sep = "|"))
}

jss_metric_snapshot <- function(x, study, keys, metric, tolerance_group,
                                statistic = c("mean", "rmse", "rate"),
                                source_column = metric) {
  statistic <- match.arg(statistic)
  if (!source_column %in% names(x)) {
    stop("Tolerance input for ", study, " is missing metric column: ", source_column, call. = FALSE)
  }
  key <- jss_metric_key(x, keys)
  value <- jss_metric_numeric(x[[source_column]])
  groups <- split(seq_along(value), key, drop = TRUE)
  rows <- lapply(names(groups), function(group_key) {
    v <- value[groups[[group_key]]]
    v <- v[is.finite(v)]
    if (!length(v)) return(NULL)
    estimate <- if (identical(statistic, "rmse")) {
      sqrt(mean(v^2))
    } else {
      mean(v)
    }
    data.frame(
      study = study, key = group_key, metric = metric,
      tolerance_group = tolerance_group, statistic = statistic,
      value = estimate, n = length(v), stringsAsFactors = FALSE
    )
  })
  rows <- rows[!vapply(rows, is.null, logical(1))]
  if (!length(rows)) {
    return(data.frame(
      study = character(), key = character(), metric = character(),
      tolerance_group = character(), statistic = character(), value = numeric(),
      n = integer(), stringsAsFactors = FALSE
    ))
  }
  do.call(rbind, rows)
}

jss_metric_snapshots_from_csv <- function(path, study, keys, specs, preprocess = identity) {
  if (!file.exists(path)) stop("Missing full-profile tolerance input: ", path, call. = FALSE)
  x <- preprocess(utils::read.csv(path, stringsAsFactors = FALSE, check.names = FALSE))
  rows <- lapply(seq_len(nrow(specs)), function(i) {
    jss_metric_snapshot(
      x, study, keys, specs$metric[[i]], specs$tolerance_group[[i]],
      specs$statistic[[i]], specs$source_column[[i]]
    )
  })
  do.call(rbind, rows)
}

jss_metric_spec <- function(source_column, metric = source_column, tolerance_group,
                            statistic = "mean") {
  data.frame(
    source_column = source_column, metric = metric,
    tolerance_group = tolerance_group, statistic = statistic,
    stringsAsFactors = FALSE
  )
}

jss_full_metric_sources <- function(settings, reference = FALSE) {
  public <- settings$public_data_dir
  if (reference) {
    paths <- list(
      bcpe = file.path(public, "bcpe-t"),
      nbi = file.path(public, "nbi-clayton"),
      missingness = file.path(public, "missingness"),
      copula = file.path(public, "copula-misspecification"),
      jvs = file.path(public, "joint-vs-separate"),
      corr = file.path(public, "correlation-misspecification")
    )
  } else {
    paths <- list(
      bcpe = file.path(settings$data_dir, "bcpe-t-full"),
      nbi = file.path(settings$data_dir, "nbi-clayton-full"),
      missingness = file.path(settings$out_dir, "missingness"),
      copula = settings$data_dir,
      jvs = file.path(settings$data_dir, "joint-vs-separate-full"),
      corr = settings$tables_dir
    )
  }
  paths
}

jss_full_metric_snapshots <- function(settings, reference = FALSE) {
  p <- jss_full_metric_sources(settings, reference)
  out <- list()
  add <- function(path, study, keys, specs, preprocess = identity) {
    out[[length(out) + 1L]] <<- jss_metric_snapshots_from_csv(path, study, keys, specs, preprocess)
  }
  fixed_specs <- rbind(
    jss_metric_spec("bias", tolerance_group = "fixed_effect_bias"),
    jss_metric_spec("rmse", tolerance_group = "fixed_effect_rmse")
  )
  add(file.path(p$bcpe, "fixed_effects_bias_rmse_table.csv"), "bcpe_t_fixed",
    c("scenario", "model", "n", "d", "parameter", "term"), fixed_specs)
  add(file.path(p$bcpe, "smooth_integrated_metrics.csv"), "bcpe_t_smooth",
    c("scenario", "model", "n", "d", "parameter"),
    jss_metric_spec("irmse", tolerance_group = "smooth_recovery"))
  add(file.path(p$bcpe, "predictive_scores_by_rep.csv"), "bcpe_t_predictive",
    c("scenario", "model", "variogram_p"), rbind(
      jss_metric_spec("test_log_score_per_obs", tolerance_group = "predictive_scores"),
      jss_metric_spec("variogram_score", tolerance_group = "predictive_scores")
    ))
  add(file.path(p$bcpe, "fit_run_log.csv"), "bcpe_t_convergence",
    c("scenario", "model"), rbind(
      jss_metric_spec("success", tolerance_group = "convergence", statistic = "rate"),
      jss_metric_spec("converged", tolerance_group = "convergence", statistic = "rate")
    ))

  nbi_preprocess <- function(x) {
    x$estimation_error <- jss_metric_numeric(x$estimate) - jss_metric_numeric(x$true_value)
    x
  }
  add(file.path(p$nbi, "fixed_effects_by_rep.csv"), "nbi_clayton_fixed",
    c("model", "parameter", "term"), rbind(
      jss_metric_spec("estimation_error", "bias", "fixed_effect_bias", "mean"),
      jss_metric_spec("estimation_error", "rmse", "fixed_effect_rmse", "rmse")
    ), nbi_preprocess)
  add(file.path(p$nbi, "smooth_integrated_metrics.csv"), "nbi_clayton_smooth",
    c("scenario", "model", "parameter"),
    jss_metric_spec("irmse", tolerance_group = "smooth_recovery"))
  add(file.path(p$nbi, "predictive_scores_by_rep.csv"), "nbi_clayton_predictive",
    c("scenario", "model", "variogram_p"), rbind(
      jss_metric_spec("test_log_score_per_obs", tolerance_group = "predictive_scores"),
      jss_metric_spec("variogram_score", tolerance_group = "predictive_scores")
    ))
  add(file.path(p$nbi, "nbi_sigma_compare_logs.csv"), "nbi_clayton_convergence",
    c("engine"), rbind(
      jss_metric_spec("success", tolerance_group = "convergence", statistic = "rate"),
      jss_metric_spec("converged", tolerance_group = "convergence", statistic = "rate")
    ))

  add(file.path(p$missingness, "fixed_term_summary_by_missingness.csv"), "missingness_fixed",
    c("missing_mechanism", "target_missing_rate", "model", "parameter"), rbind(
      jss_metric_spec("mean_abs_bias", tolerance_group = "fixed_effect_bias"),
      jss_metric_spec("mean_rmse", tolerance_group = "fixed_effect_rmse"),
      jss_metric_spec("mean_coverage_95", tolerance_group = "convergence")
    ))

  copula_path <- if (reference) file.path(p$copula, "results.csv") else file.path(p$copula, "07-gamma-copula-misspecification-results.csv")
  add(copula_path, "copula_misspecification",
    c("generating_copula", "fitted_copula", "tau_label", "n_subject", "n_time"), rbind(
      jss_metric_spec("mu_bias", tolerance_group = "fixed_effect_bias"),
      jss_metric_spec("sigma_bias", tolerance_group = "fixed_effect_bias"),
      jss_metric_spec("margin_param_rmse", tolerance_group = "fixed_effect_rmse"),
      jss_metric_spec("tau_abs_error", tolerance_group = "smooth_recovery"),
      jss_metric_spec("benchmark_neg_log_score", tolerance_group = "predictive_scores"),
      jss_metric_spec("success", tolerance_group = "convergence", statistic = "rate"),
      jss_metric_spec("converged", tolerance_group = "convergence", statistic = "rate")
    ))

  for (bundle in c("normal", "gamma", "nbi")) {
    path <- file.path(p$jvs, bundle, "data", "03-joint-vs-separate-optimization-deltas.csv")
    add(path, paste0("joint_vs_separate_", bundle), c("case_id", "family"), rbind(
      jss_metric_spec("delta_test_log_score_per_obs", tolerance_group = "predictive_scores"),
      jss_metric_spec("delta_heldout_variogram_score_p05", tolerance_group = "predictive_scores"),
      jss_metric_spec("delta_heldout_variogram_score_p2", tolerance_group = "predictive_scores"),
      jss_metric_spec("rs_joint_success", tolerance_group = "convergence", statistic = "rate"),
      jss_metric_spec("rs_separate_success", tolerance_group = "convergence", statistic = "rate")
    ))
  }

  corr_files <- c(
    "table_2_external_correlation.csv" = "08-simulation-sensitivity-correlation-misspecification-external-correlation.csv",
    "table_3_flexible_correlation.csv" = "08-simulation-sensitivity-correlation-misspecification-flexible-correlation.csv"
  )
  for (reference_name in names(corr_files)) {
    path <- if (reference) file.path(p$corr, reference_name) else file.path(p$corr, corr_files[[reference_name]])
    x <- utils::read.csv(path, stringsAsFactors = FALSE, check.names = FALSE)
    keys <- names(x)[1:2]
    for (column in names(x)[-(1:2)]) {
      metric_label <- paste(x[[keys[[2L]]]], column, sep = "|")
      values <- jss_metric_numeric(x[[column]])
      tmp <- data.frame(row_key = paste0(keys[[1L]], "=", x[[keys[[1L]]]]), metric_value = values)
      metric_type <- ifelse(grepl("coverage", metric_label, ignore.case = TRUE), "convergence",
        ifelse(grepl("RMSE|MAE|SE ratio", metric_label, ignore.case = TRUE), "fixed_effect_rmse", "predictive_scores"))
      for (i in seq_len(nrow(tmp))) {
        if (!is.finite(tmp$metric_value[[i]])) next
        out[[length(out) + 1L]] <- data.frame(
          study = paste0("correlation_misspecification_", tools::file_path_sans_ext(reference_name)),
          key = tmp$row_key[[i]], metric = metric_label[[i]], tolerance_group = metric_type[[i]],
          statistic = "reported", value = tmp$metric_value[[i]], n = 1L, stringsAsFactors = FALSE
        )
      }
    }
  }
  result <- do.call(rbind, out)
  result[order(result$study, result$key, result$metric), , drop = FALSE]
}

jss_compare_full_metrics <- function(actual, reference, tolerances) {
  keys <- c("study", "key", "metric", "tolerance_group", "statistic")
  comparison <- merge(
    reference[, c(keys, "value", "n")], actual[, c(keys, "value", "n")],
    by = keys, all = TRUE, suffixes = c("_reference", "_actual"), sort = TRUE
  )
  idx <- match(comparison$tolerance_group, tolerances$artifact_group)
  if (anyNA(idx)) stop("No registered tolerance for: ", paste(unique(comparison$tolerance_group[is.na(idx)]), collapse = ", "), call. = FALSE)
  comparison$absolute_tolerance <- tolerances$absolute_tolerance[idx]
  comparison$relative_tolerance <- tolerances$relative_tolerance[idx]
  comparison$absolute_difference <- abs(comparison$value_actual - comparison$value_reference)
  comparison$allowed_difference <- comparison$absolute_tolerance +
    comparison$relative_tolerance * abs(comparison$value_reference)
  comparison$status <- ifelse(
    !is.finite(comparison$value_reference), "missing_reference",
    ifelse(!is.finite(comparison$value_actual), "missing_actual",
      ifelse(comparison$absolute_difference <= comparison$allowed_difference, "pass", "outside_tolerance"))
  )
  comparison
}

jss_validate_full_tolerances <- function(settings) {
  if (!identical(settings$profile, "full")) return(character())
  tolerances <- utils::read.csv(file.path(settings$root, "paper", "tolerances.csv"), stringsAsFactors = FALSE)
  actual <- jss_full_metric_snapshots(settings, reference = FALSE)
  reference <- jss_full_metric_snapshots(settings, reference = TRUE)
  comparison <- jss_compare_full_metrics(actual, reference, tolerances)
  path <- file.path(settings$logs_dir, "full-tolerance-audit.csv")
  utils::write.csv(comparison, path, row.names = FALSE)
  failed <- comparison$status != "pass"
  if (any(failed)) {
    preview <- paste(utils::head(paste(comparison$study[failed], comparison$key[failed], comparison$metric[failed], comparison$status[failed], sep = ":"), 12L), collapse = ", ")
    stop("Full-profile metric tolerance validation failed (", sum(failed), " rows): ", preview, call. = FALSE)
  }
  path
}

jss_png_nonblank <- function(path) {
  if (!file.exists(path) || file.info(path)$size < 1000L) return(FALSE)
  con <- file(path, "rb"); on.exit(close(con), add = TRUE)
  header <- readBin(con, "raw", n = 24L)
  if (length(header) < 24L || !identical(as.integer(header[1:8]), c(137L, 80L, 78L, 71L, 13L, 10L, 26L, 10L))) return(FALSE)
  value <- function(x) sum(as.integer(x) * 256^(3:0))
  value(header[17:20]) > 10L && value(header[21:24]) > 10L
}

jss_file_sha256 <- function(path, canonical_text = grepl("[.](csv|tex|txt|md)$", path, ignore.case = TRUE)) {
  if (!requireNamespace("digest", quietly = TRUE)) stop("Package 'digest' is required for SHA-256 verification.", call. = FALSE)
  if (canonical_text) {
    text <- paste(readLines(path, warn = FALSE, encoding = "UTF-8"), collapse = "\n")
    text <- sub("^\ufeff", "", text)
    return(digest::digest(charToRaw(text), algo = "sha256", serialize = FALSE))
  }
  digest::digest(file = path, algo = "sha256", serialize = FALSE)
}

jss_git_dirty <- function(root) {
  out <- tryCatch(system2("git", c("-C", shQuote(root), "status", "--porcelain"), stdout = TRUE, stderr = TRUE), error = function(e) "unknown")
  if (identical(out, "unknown")) "unknown" else if (length(out)) "dirty" else "clean"
}

jss_write_run_metadata <- function(settings, started, bootstrap, store) {
  meta <- targets::tar_meta(store = store, fields = c("name", "seconds", "warnings", "error"))
  timing <- data.frame(name = meta$name, seconds = as.numeric(meta$seconds), warnings = as.character(meta$warnings), error = as.character(meta$error), stringsAsFactors = FALSE)
  timing <- timing[!is.na(timing$seconds) | (!is.na(timing$warnings) & nzchar(timing$warnings)) | (!is.na(timing$error) & nzchar(timing$error)), , drop = FALSE]
  utils::write.csv(timing, file.path(settings$logs_dir, "target-timings.csv"), row.names = FALSE)
  events <- timing[(!is.na(timing$warnings) & nzchar(timing$warnings)) | (!is.na(timing$error) & nzchar(timing$error)), , drop = FALSE]
  if (!nrow(events)) {
    events <- data.frame(name = NA_character_, seconds = NA_real_, warnings = NA_character_, error = NA_character_, status = "no target warnings or errors", stringsAsFactors = FALSE)
  } else {
    events$status <- "target warning or error"
  }
  utils::write.csv(events, file.path(settings$logs_dir, "target-events.csv"), row.names = FALSE)
  lock <- file.path(settings$root, "paper", "renv.lock")
  run <- data.frame(
    profile = settings$profile, started = format(started, "%Y-%m-%dT%H:%M:%S%z"), finished = format(Sys.time(), "%Y-%m-%dT%H:%M:%S%z"),
    git_sha = jss_git_sha(settings$root), git_state = jss_git_dirty(settings$root), r_version = R.version.string,
    platform = R.version$platform, workers = settings$workers, seed = settings$seed,
    library_id = paste(tail(strsplit(bootstrap$library, "/", fixed = TRUE)[[1L]], 2L), collapse = "/"),
    lockfile_sha256 = if (file.exists(lock)) jss_file_sha256(lock) else NA_character_,
    restored = bootstrap$restored, source_installed = bootstrap$source_installed, stringsAsFactors = FALSE
  )
  utils::write.csv(run, file.path(settings$logs_dir, "run-metadata.csv"), row.names = FALSE)
  control_inputs <- c(
    file.path(settings$root, "DESCRIPTION"),
    list.files(file.path(settings$root, "R"), pattern = "[.]R$", full.names = TRUE, recursive = TRUE),
    file.path(settings$root, "paper", c("manifest.csv", "seeds.csv", "tolerances.csv", "renv.lock", "_targets.R", "replicate.R", "bootstrap.R")),
    file.path(settings$root, "paper", "R", c(
      "replication-helpers.R", "public-paper-producers.R",
      "01-simulation-bcpe-t.R", "02-simulation-delaporte-clayton.R",
      "03-joint-vs-separate-optimization.R", "04-missingness-dropout-sensitivity.R",
      "07-gamma-copula-misspecification.R",
      "08-simulation-sensitivity-correlation-misspecification.R"
    ))
  )
  inputs <- unique(control_inputs[file.exists(control_inputs)])
  if (!identical(settings$profile, "smoke")) {
    inputs <- c(inputs, list.files(settings$public_data_dir, full.names = TRUE, recursive = TRUE))
  }
  inputs <- inputs[file.info(inputs)$isdir %in% FALSE & !grepl("[.]log$", inputs, ignore.case = TRUE)]
  input_log <- data.frame(
    file = substring(normalizePath(inputs, winslash = "/", mustWork = TRUE), nchar(normalizePath(settings$root, winslash = "/", mustWork = TRUE)) + 2L),
    sha256 = vapply(inputs, jss_file_sha256, character(1)),
    raw_sha256 = vapply(inputs, jss_file_sha256, character(1), canonical_text = FALSE),
    hash_mode = ifelse(grepl("[.](csv|tex|txt|md)$", inputs, ignore.case = TRUE), "canonical_text_sha256", "raw_file_sha256"),
    bytes = file.info(inputs)$size,
    stringsAsFactors = FALSE
  )
  utils::write.csv(input_log, file.path(settings$logs_dir, "input-hashes.csv"), row.names = FALSE)

  network <- suppressWarnings(targets::tar_network(
    targets_only = TRUE, outdated = FALSE, callr_function = NULL,
    script = file.path(settings$root, "paper", "_targets.R"), store = store
  ))
  utils::write.csv(as.data.frame(network$vertices), file.path(settings$logs_dir, "target-graph-vertices.csv"), row.names = FALSE)
  utils::write.csv(as.data.frame(network$edges), file.path(settings$logs_dir, "target-graph-edges.csv"), row.names = FALSE)

  fit_path <- file.path(settings$logs_dir, "fit-events.csv")
  fit_audit <- if (file.exists(fit_path)) utils::read.csv(fit_path, stringsAsFactors = FALSE) else data.frame()
  fit_events <- nrow(fit_audit)
  nonconverged <- if ("event_type" %in% names(fit_audit)) sum(fit_audit$event_type == "optimizer_nonconvergence") else 0L
  execution_failures <- if ("event_type" %in% names(fit_audit)) sum(fit_audit$event_type == "execution_failure") else 0L
  summary <- c(
    paste0("# Reviewer run summary: ", settings$profile), "",
    paste0("- Status: completed"),
    paste0("- Started: ", run$started),
    paste0("- Finished: ", run$finished),
    paste0("- Git: `", run$git_sha, "` (", run$git_state, ")"),
    paste0("- R/platform: ", run$r_version, " / ", run$platform),
    paste0("- Lock restored: ", run$restored, "; checked-out source installed: ", run$source_installed),
    paste0("- Targets recorded: ", nrow(timing)),
    paste0("- Target warnings/errors: ", if (all(is.na(events$name))) 0L else nrow(events)),
    paste0("- Structured fit events: ", fit_events),
    paste0("- Optimizer nonconvergence events: ", nonconverged),
    paste0("- Execution failures: ", execution_failures),
    if (nonconverged > 0L) paste0("- Interpretation: ", nonconverged, " fits completed without an execution error but did not converge under this profile's iteration/time limits; do not treat them as converged estimates.") else NULL,
    "",
    "See `target-timings.csv`, `target-events.csv`, `fit-events.csv`,",
    "`target-graph-vertices.csv`, `target-graph-edges.csv`, `input-hashes.csv`,",
    "`output_hashes.csv`, and `figure-verification.csv` for the auditable detail."
  )
  writeLines(summary, file.path(settings$logs_dir, "reviewer-summary.md"), useBytes = TRUE)
  invisible(run)
}

jss_write_provenance_hashes <- function(settings) {
  path <- file.path(settings$logs_dir, "provenance-hashes.csv")
  files <- list.files(settings$logs_dir, full.names = TRUE, recursive = FALSE)
  files <- files[normalizePath(files, winslash = "/", mustWork = FALSE) != normalizePath(path, winslash = "/", mustWork = FALSE)]
  out <- data.frame(
    file = basename(files),
    sha256 = vapply(files, jss_file_sha256, character(1)),
    bytes = file.info(files)$size,
    stringsAsFactors = FALSE
  )
  utils::write.csv(out, path, row.names = FALSE)
  path
}

jss_check_target_warnings <- function(store, strict = FALSE) {
  meta <- targets::tar_meta(store = store, fields = c("name", "warnings"))
  warned <- meta[!is.na(meta$warnings) & nzchar(meta$warnings), c("name", "warnings"), drop = FALSE]
  if (nrow(warned) && strict) stop("Unexpected target warning(s): ", paste(warned$name, collapse = ", "), call. = FALSE)
  invisible(warned)
}
