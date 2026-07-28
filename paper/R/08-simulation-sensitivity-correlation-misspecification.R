jss_corr_misspec_module_id <- function() {
  "08-simulation-sensitivity-correlation-misspecification"
}

jss_corr_misspec_source_run_dir <- function(settings) {
  override <- Sys.getenv("GAMLSS_LONGITUDINAL_JSS_CORR_MISSPEC_RUN_DIR", unset = "")
  if (nzchar(override)) {
    return(normalizePath(override, winslash = "/", mustWork = TRUE))
  }
  normalizePath(
    file.path(
      settings$root,
      "results",
      "jss-exploratory",
      "05-standard-model-benchmarking",
      "rs-joint-standard-model-grid",
      "run_20260619_t20_t50_combined"
    ),
    winslash = "/",
    mustWork = TRUE
  )
}

jss_corr_misspec_paths <- function(settings) {
  module_id <- jss_corr_misspec_module_id()
  list(
    provenance = file.path(settings$data_dir, paste0(module_id, "-provenance.csv")),
    validation = file.path(settings$tables_dir, paste0(module_id, "-validation.csv")),
    scenario_design = file.path(settings$tables_dir, paste0(module_id, "-scenario-design.csv")),
    external_correlation = file.path(settings$tables_dir, paste0(module_id, "-external-correlation.csv")),
    flexible_correlation = file.path(settings$tables_dir, paste0(module_id, "-flexible-correlation.csv")),
    scenario_design_tex = file.path(settings$tables_dir, paste0(module_id, "-scenario-design.tex")),
    external_correlation_tex = file.path(settings$tables_dir, paste0(module_id, "-external-correlation.tex")),
    flexible_correlation_tex = file.path(settings$tables_dir, paste0(module_id, "-flexible-correlation.tex")),
    appendix_external_correlation_t50_tex = file.path(settings$tables_dir, paste0(module_id, "-appendix-external-correlation-t50.tex")),
    appendix_flexible_correlation_t50_tex = file.path(settings$tables_dir, paste0(module_id, "-appendix-flexible-correlation-t50.tex")),
    story_tables_tex = file.path(settings$tables_dir, paste0(module_id, "-story-tables.tex")),
    appendix_t50_tex = file.path(settings$tables_dir, paste0(module_id, "-appendix-t50-tables.tex")),
    se_ratio_figure = file.path(settings$figures_dir, paste0(module_id, "-slope-se-ratio.png"))
  )
}

jss_corr_misspec_source_paths <- function(source_run_dir) {
  story_dir <- file.path(source_run_dir, "story_tables")
  sandwich_dir <- file.path(source_run_dir, "sandwich_t20_grid")
  list(
    validation = file.path(sandwich_dir, "rs_joint_sandwich_validation_summary.csv"),
    status = file.path(sandwich_dir, "rs_joint_sandwich_status_by_rep.csv"),
    scenario_design = file.path(story_dir, "table_1_scenario_design.csv"),
    external_correlation = file.path(story_dir, "table_2_external_correlation.csv"),
    flexible_correlation = file.path(story_dir, "table_3_flexible_correlation.csv"),
    scenario_design_tex = file.path(story_dir, "table_1_scenario_design.tex"),
    external_correlation_tex = file.path(story_dir, "table_2_external_correlation.tex"),
    flexible_correlation_tex = file.path(story_dir, "table_3_flexible_correlation.tex"),
    appendix_external_correlation_t50_tex = file.path(story_dir, "appendix_table_1_external_correlation_t50.tex"),
    appendix_flexible_correlation_t50_tex = file.path(story_dir, "appendix_table_2_flexible_correlation_t50.tex"),
    story_tables_tex = file.path(story_dir, "story_tables.tex"),
    appendix_t50_tex = file.path(story_dir, "appendix_t50_tables.tex")
  )
}

jss_corr_misspec_read_csv <- function(path) {
  if (!file.exists(path)) {
    stop("Missing correlation misspecification source artifact: ", path, call. = FALSE)
  }
  utils::read.csv(path, stringsAsFactors = FALSE, check.names = FALSE)
}

jss_corr_misspec_copy <- function(from, to) {
  if (!file.exists(from)) {
    stop("Missing correlation misspecification source artifact: ", from, call. = FALSE)
  }
  dir.create(dirname(to), recursive = TRUE, showWarnings = FALSE)
  ok <- file.copy(from, to, overwrite = TRUE)
  if (!isTRUE(ok)) {
    stop("Failed to copy correlation misspecification artifact from ", from, " to ", to, call. = FALSE)
  }
  normalizePath(to, winslash = "/", mustWork = TRUE)
}

jss_corr_misspec_status_summary <- function(status, validation) {
  success <- status$success %in% c(TRUE, "TRUE", "True", "true", "1")
  data.frame(
    module_id = jss_corr_misspec_module_id(),
    expected_cases = validation$expected_cases[[1L]],
    completed_status_cases = validation$completed_status_cases[[1L]],
    successful_cases = validation$successful_cases[[1L]],
    failed_cases = validation$failed_cases[[1L]],
    nonfinite_sandwich_se = validation$nonfinite_sandwich_se[[1L]],
    max_abs_estimate_difference_vs_original = validation$max_abs_estimate_difference_vs_original[[1L]],
    median_elapsed_sec = validation$median_elapsed_sec[[1L]],
    median_fit_elapsed_sec = validation$median_fit_elapsed_sec[[1L]],
    median_sandwich_elapsed_sec = validation$median_sandwich_elapsed_sec[[1L]],
    status_rows = nrow(status),
    status_success_rows = sum(success),
    status_failed_rows = sum(!success),
    stringsAsFactors = FALSE
  )
}

jss_corr_misspec_provenance <- function(settings, source_run_dir, source_paths, output_paths) {
  outputs <- unlist(output_paths, use.names = TRUE)
  sources <- unlist(source_paths, use.names = TRUE)
  matched_outputs <- outputs[match(names(sources), names(outputs))]
  matched_outputs[is.na(matched_outputs)] <- NA_character_
  data.frame(
    module_id = jss_corr_misspec_module_id(),
    profile = settings$profile,
    source_run_dir = normalizePath(source_run_dir, winslash = "/", mustWork = TRUE),
    source_artifact = unname(names(sources)),
    source_path = unname(normalizePath(sources, winslash = "/", mustWork = FALSE)),
    output_path = unname(ifelse(is.na(matched_outputs), NA_character_, normalizePath(matched_outputs, winslash = "/", mustWork = FALSE))),
    generated_at = format(Sys.time(), "%Y-%m-%dT%H:%M:%S%z"),
    stringsAsFactors = FALSE
  )
}

jss_corr_misspec_plot_se_ratio <- function(external_table, flexible_table, path) {
  row_is_se <- function(x) identical(x, "Slope SE ratio (target: 1.0)")
  rows <- list()

  ext_se <- external_table[vapply(external_table$Metric, row_is_se, logical(1)), , drop = FALSE]
  ext_value_cols <- setdiff(names(ext_se), c("Level", "Metric"))
  for (col in ext_value_cols) {
    parts <- strsplit(col, " ", fixed = TRUE)[[1L]]
    scenario <- paste(parts[seq_len(length(parts) - 1L)], collapse = " ")
    method <- parts[[length(parts)]]
    rows[[length(rows) + 1L]] <- data.frame(
      scenario = paste(scenario, ext_se$Level),
      method = method,
      se_ratio = suppressWarnings(as.numeric(ext_se[[col]])),
      stringsAsFactors = FALSE
    )
  }

  flex_se <- flexible_table[vapply(flexible_table$Metric, row_is_se, logical(1)), , drop = FALSE]
  flex_value_cols <- setdiff(names(flex_se), c("Scenario", "Metric"))
  for (col in flex_value_cols) {
    rows[[length(rows) + 1L]] <- data.frame(
      scenario = flex_se$Scenario,
      method = sub(" T=20$", "", col),
      se_ratio = suppressWarnings(as.numeric(flex_se[[col]])),
      stringsAsFactors = FALSE
    )
  }

  plot_data <- do.call(rbind, rows)
  plot_data <- plot_data[is.finite(plot_data$se_ratio), , drop = FALSE]
  plot_data$method <- factor(
    plot_data$method,
    levels = c("gamlss.long", "robust", "exch.", "AR(1)", "unstr.", "glm")
  )
  p <- ggplot2::ggplot(plot_data, ggplot2::aes(x = method, y = se_ratio, fill = method)) +
    ggplot2::geom_col(width = 0.7, show.legend = FALSE) +
    ggplot2::geom_hline(yintercept = 1, linewidth = 0.3, linetype = 2) +
    ggplot2::facet_wrap(~ scenario, ncol = 2) +
    ggplot2::coord_cartesian(ylim = c(0, max(2.1, plot_data$se_ratio, na.rm = TRUE))) +
    ggplot2::labs(x = NULL, y = "Slope SE ratio", title = "Correlation misspecification sensitivity") +
    ggplot2::theme_minimal(base_size = 10) +
    ggplot2::theme(axis.text.x = ggplot2::element_text(angle = 35, hjust = 1))
  ggplot2::ggsave(path, p, width = 8.5, height = 7, dpi = 320, bg = "white")
  normalizePath(path, winslash = "/", mustWork = TRUE)
}

jss_run_08_simulation_sensitivity_correlation_misspecification <- function(settings) {
  source_run_dir <- jss_corr_misspec_source_run_dir(settings)
  source_paths <- jss_corr_misspec_source_paths(source_run_dir)
  paths <- jss_corr_misspec_paths(settings)

  validation <- jss_corr_misspec_read_csv(source_paths$validation)
  status <- jss_corr_misspec_read_csv(source_paths$status)
  external_table <- jss_corr_misspec_read_csv(source_paths$external_correlation)
  flexible_table <- jss_corr_misspec_read_csv(source_paths$flexible_correlation)

  copied <- c(
    scenario_design = jss_corr_misspec_copy(source_paths$scenario_design, paths$scenario_design),
    external_correlation = jss_corr_misspec_copy(source_paths$external_correlation, paths$external_correlation),
    flexible_correlation = jss_corr_misspec_copy(source_paths$flexible_correlation, paths$flexible_correlation),
    scenario_design_tex = jss_corr_misspec_copy(source_paths$scenario_design_tex, paths$scenario_design_tex),
    external_correlation_tex = jss_corr_misspec_copy(source_paths$external_correlation_tex, paths$external_correlation_tex),
    flexible_correlation_tex = jss_corr_misspec_copy(source_paths$flexible_correlation_tex, paths$flexible_correlation_tex),
    appendix_external_correlation_t50_tex = jss_corr_misspec_copy(source_paths$appendix_external_correlation_t50_tex, paths$appendix_external_correlation_t50_tex),
    appendix_flexible_correlation_t50_tex = jss_corr_misspec_copy(source_paths$appendix_flexible_correlation_t50_tex, paths$appendix_flexible_correlation_t50_tex),
    story_tables_tex = jss_corr_misspec_copy(source_paths$story_tables_tex, paths$story_tables_tex),
    appendix_t50_tex = jss_corr_misspec_copy(source_paths$appendix_t50_tex, paths$appendix_t50_tex)
  )

  validation_out <- jss_corr_misspec_status_summary(status, validation)
  utils::write.csv(validation_out, paths$validation, row.names = FALSE)

  provenance <- jss_corr_misspec_provenance(settings, source_run_dir, source_paths, paths)
  utils::write.csv(provenance, paths$provenance, row.names = FALSE)

  figure <- jss_corr_misspec_plot_se_ratio(external_table, flexible_table, paths$se_ratio_figure)

  list(
    module_id = jss_corr_misspec_module_id(),
    title = "Simulation sensitivity to correlation misspecification",
    status = "current",
    data = paths$provenance,
    tables = c(paths$validation, unname(copied)),
    figures = figure,
    source_run_dir = source_run_dir
  )
}
