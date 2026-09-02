jss_run_public_workflow_figures <- function(settings) {
  set.seed(settings$seed + 1L)
  intro <- file.path(settings$figures_dir, "intro - copula examples.png")
  grDevices::png(intro, width = 2400, height = 650, res = 220, bg = "white")
  old <- graphics::par(mfrow = c(1, 4), mar = c(3, 3, 2, 1))
  on.exit(graphics::par(old), add = TRUE)
  examples <- list(
    list(name = "Gaussian", family = 1, par = sin(pi / 4), par2 = 0),
    list(name = "Student t (df = 4)", family = 2, par = sin(pi / 4), par2 = 4),
    list(name = "Clayton", family = 3, par = 2, par2 = 0),
    list(name = "Gumbel", family = 4, par = 2, par2 = 0)
  )
  for (example in examples) {
    uv <- VineCopula::BiCopSim(450, family = example$family, par = example$par, par2 = example$par2)
    graphics::plot(uv[, 1], uv[, 2], pch = 16, cex = 0.35,
      col = grDevices::adjustcolor("#0072B2", 0.45), xlab = "u", ylab = "v",
      main = example$name)
  }
  grDevices::dev.off()

  set.seed(settings$seed + 2L)
  dat <- gamlss.longitudinal::simulate_longitudinal_dataset(
    n = if (identical(settings$profile, "smoke")) 35L else 80L, times = 1:3,
    margin_dist = gamlss.dist::BCPE(), copula_dist = "t",
    margin_params = list(mu = 2, sigma = 0.3, nu = 1, tau = 2),
    copula_params = list(tau = 0.4, zeta = 5), seed = settings$seed + 2L
  )
  dist <- file.path(settings$figures_dir, "software - plot_dist.png")
  dist_dat <- dat[dat$time %in% sort(unique(dat$time))[1:2], , drop = FALSE]
  p <- gamlss.longitudinal::plot_dist(
    dataset = dist_dat, margin_dist = gamlss.dist::BCPE(),
    subject_var = "subject", time_var = "time", response_var = "response"
  )
  ggplot2::ggsave(dist, p, width = 8.5, height = 7.0, dpi = 320, bg = "white")

  diag <- file.path(settings$figures_dir, "software - plot_copula_diagnostics.png")
  fit <- suppressWarnings(gamlss.longitudinal::gamlss_longitudinal(
    dataset = dat, margin_dist = gamlss.dist::BCPE(), copula_dist = "t",
    time_var = "time", subject_var = "subject", mu.formula = response ~ 1,
    sigma.formula = ~1, nu.formula = ~1, tau.formula = ~1,
    theta.formula = ~1, zeta.formula = ~1, max_outer_iter = 4L,
    max_inner_iter = 4L, max_elapsed_sec = 30, compute_vcov = FALSE, verbose = 0
  ))
  diagnostics <- suppressWarnings(gamlss.longitudinal::plot_copula_diagnostics(
    fit, residual_lags = 1:2, plot = TRUE
  ))
  wrap_label <- function(x, width) {
    if (is.null(x) || !length(x) || is.na(x)) return(x)
    paste(strwrap(x, width = width), collapse = "\n")
  }
  paper_panels <- lapply(diagnostics$plots, function(panel) {
    panel +
      ggplot2::labs(
        title = wrap_label(panel$labels$title, 29L),
        subtitle = wrap_label(panel$labels$subtitle, 48L)
      ) +
      ggplot2::theme(
        plot.title = ggplot2::element_text(size = 7.5, lineheight = 0.9),
        plot.subtitle = ggplot2::element_text(size = 4.8, lineheight = 0.85),
        axis.title = ggplot2::element_text(size = 6.5),
        axis.text = ggplot2::element_text(size = 5.5),
        legend.title = ggplot2::element_text(size = 5.5),
        legend.text = ggplot2::element_text(size = 5.2)
      )
  })
  paper_dashboard <- do.call(
    ggpubr::ggarrange,
    c(paper_panels, list(ncol = 3L, nrow = 3L))
  )
  ggplot2::ggsave(diag, paper_dashboard, width = 10, height = 7, dpi = 320, bg = "white")
  diag_csv <- file.path(settings$data_dir, "software-copula-diagnostic-kendall.csv")
  utils::write.csv(diagnostics$kendall_summary, diag_csv, row.names = FALSE)
  list(module_id = "00-public-workflow", status = "regenerated", data = diag_csv, tables = character(), figures = c(intro, dist, diag))
}

jss_jvs_full_bundle_specs <- function() {
  data.frame(
    bundle = c("normal", "gamma", "nbi"),
    output_name = c(
      "normal-joint-vs-separate-six-case-median-iqr-table.tex",
      "gamma-joint-vs-separate-six-case-median-iqr-table.tex",
      "negative-binomial-joint-vs-separate-six-case-median-iqr-table.tex"
    ),
    label = c("tab:jvs-small-sample-no", "tab:jvs-small-sample-ga", "tab:jvs-small-sample-nbi"),
    caption = c(
      "Normal-margin joint versus separate optimisation results.",
      "Gamma-margin joint versus separate optimisation results.",
      "Negative-binomial-margin joint versus separate optimisation results."
    ),
    stringsAsFactors = FALSE
  )
}

jss_jvs_full_cases <- function(settings, bundle) {
  reference <- file.path(
    settings$public_data_dir, "joint-vs-separate", bundle, "data",
    "03-joint-vs-separate-optimization-deltas.csv"
  )
  x <- utils::read.csv(reference, stringsAsFactors = FALSE)
  fields <- c(
    "case_id", "hypothesis_role", "family", "copula", "n", "time_points",
    "total_observations", "mu_strength", "sigma_strength", "theta_strength", "time_shape"
  )
  if (!all(fields %in% names(x))) {
    stop("Incomplete joint-versus-separate public design bundle: ", bundle, call. = FALSE)
  }
  cases <- unique(x[, fields, drop = FALSE])
  cases <- cases[order(cases$case_id), , drop = FALSE]
  expected <- sprintf("JVS%02d", 1:6)
  if (!identical(as.character(cases$case_id), expected)) {
    stop("Joint-versus-separate bundle must define JVS01-JVS06 exactly: ", bundle, call. = FALSE)
  }
  cases
}

jss_jvs_full_checkpoint_path <- function(directory, case_id, rep) {
  file.path(directory, sprintf("%s-rep-%03d.csv", tolower(case_id), rep))
}

jss_jvs_full_checkpoint_complete <- function(path, case_id, rep) {
  if (!file.exists(path) || file.info(path)$size < 100L) return(FALSE)
  x <- tryCatch(utils::read.csv(path, stringsAsFactors = FALSE), error = function(e) NULL)
  !is.null(x) && nrow(x) == 2L &&
    all(x$case_id == case_id) && all(x$joint_review_rep == rep) &&
    setequal(x$method, c("rs_separate", "rs_joint"))
}

jss_jvs_full_write_checkpoint <- function(x, path) {
  tmp <- paste0(path, ".tmp")
  utils::write.csv(x, tmp, row.names = FALSE)
  if (file.exists(path)) unlink(path)
  if (!file.rename(tmp, path)) stop("Could not atomically install checkpoint: ", path, call. = FALSE)
  path
}

jss_run_03_full_bundle <- function(settings, bundle, reps = NULL) {
  spec <- jss_jvs_full_bundle_specs()
  spec <- spec[spec$bundle == bundle, , drop = FALSE]
  if (nrow(spec) != 1L) stop("Unknown joint-versus-separate bundle: ", bundle, call. = FALSE)
  cases <- jss_jvs_full_cases(settings, bundle)
  reference_deltas <- utils::read.csv(
    file.path(settings$public_data_dir, "joint-vs-separate", bundle, "data",
      "03-joint-vs-separate-optimization-deltas.csv"),
    stringsAsFactors = FALSE
  )
  published_reps <- max(as.integer(reference_deltas$joint_review_rep), na.rm = TRUE)
  if (is.null(reps)) reps <- published_reps
  if (!is.finite(reps) || reps < 1L) stop("Invalid replicate count for bundle: ", bundle, call. = FALSE)
  base <- file.path(settings$data_dir, "joint-vs-separate-full", bundle)
  checkpoint_dir <- file.path(base, "checkpoints")
  data_dir <- file.path(base, "data")
  table_dir <- file.path(base, "tables")
  dir.create(checkpoint_dir, recursive = TRUE, showWarnings = FALSE)
  dir.create(data_dir, recursive = TRUE, showWarnings = FALSE)
  dir.create(table_dir, recursive = TRUE, showWarnings = FALSE)
  cfg <- jss_joint_simulation_settings(settings)
  cfg$reps <- as.integer(reps)

  checkpoint_paths <- character()
  for (case_idx in seq_len(nrow(cases))) {
    case <- cases[case_idx, , drop = FALSE]
    for (rep_idx in seq_len(cfg$reps)) {
      path <- jss_jvs_full_checkpoint_path(checkpoint_dir, case$case_id[[1L]], rep_idx)
      checkpoint_paths <- c(checkpoint_paths, path)
      if (jss_jvs_full_checkpoint_complete(path, case$case_id[[1L]], rep_idx)) next
      result <- jss_joint_run_case_rep(case, rep_idx, settings, cfg)
      jss_jvs_full_write_checkpoint(result, path)
    }
  }
  results <- jss_joint_bind_rows(lapply(checkpoint_paths, utils::read.csv, stringsAsFactors = FALSE))
  results$profile <- settings$profile
  deltas <- jss_joint_delta_table(results)
  summary <- jss_joint_summary_table(results, deltas, cases)
  results_path <- file.path(data_dir, "03-joint-vs-separate-optimization-results.csv")
  deltas_path <- file.path(data_dir, "03-joint-vs-separate-optimization-deltas.csv")
  summary_path <- file.path(table_dir, "03-joint-vs-separate-optimization-summary.csv")
  utils::write.csv(results, results_path, row.names = FALSE)
  utils::write.csv(deltas, deltas_path, row.names = FALSE)
  utils::write.csv(summary, summary_path, row.names = FALSE)

  manuscript_table <- file.path(settings$tables_dir, spec$output_name[[1L]])
  callr::rscript(
    file.path(settings$root, "paper", "R", "write-jvs-small-sample-latex-table.R"),
    cmdargs = c(base, file.path("tables", basename(summary_path)), manuscript_table,
      spec$label[[1L]], spec$caption[[1L]]),
    wd = settings$root, show = FALSE
  )
  list(
    bundle = bundle, published_reps = published_reps, checkpoints = checkpoint_paths, data = c(results_path, deltas_path),
    tables = c(summary_path, manuscript_table)
  )
}

jss_legacy_run_03_joint_vs_separate_full <- function(settings) {
  specs <- jss_jvs_full_bundle_specs()
  runs <- lapply(specs$bundle, function(bundle) jss_run_03_full_bundle(settings, bundle))
  list(
    module_id = "03-joint-vs-separate", status = "regenerated",
    data = unlist(lapply(runs, `[[`, "data"), use.names = FALSE),
    tables = unlist(lapply(runs, `[[`, "tables"), use.names = FALSE), figures = character()
  )
}

jss_legacy_run_03_joint_vs_separate <- function(settings) {
  if (identical(settings$profile, "full")) return(jss_legacy_run_03_joint_vs_separate_full(settings))
  src <- file.path(settings$public_data_dir, "joint-vs-separate")
  names <- c("normal-joint-vs-separate-six-case-median-iqr-table.tex", "gamma-joint-vs-separate-six-case-median-iqr-table.tex", "negative-binomial-joint-vs-separate-six-case-median-iqr-table.tex")
  out <- file.path(settings$tables_dir, names)
  specs <- data.frame(
    bundle = c("normal", "gamma", "nbi"),
    label = c("tab:jvs-small-sample-no", "tab:jvs-small-sample-ga", "tab:jvs-small-sample-nbi"),
    caption = c("Normal-margin joint versus separate optimisation results.", "Gamma-margin joint versus separate optimisation results.", "Negative-binomial-margin joint versus separate optimisation results."),
    stringsAsFactors = FALSE
  )
  script <- file.path(settings$root, "paper", "R", "write-jvs-small-sample-latex-table.R")
  for (i in seq_len(nrow(specs))) {
    run_dir <- file.path(src, specs$bundle[[i]])
    callr::rscript(script, cmdargs = c(run_dir, file.path("tables", "03-joint-vs-separate-optimization-summary.csv"), out[[i]], specs$label[[i]], specs$caption[[i]]), wd = settings$root, show = FALSE)
  }
  list(module_id = "03-joint-vs-separate", status = "regenerated", data = file.path(src, specs$bundle, "data", "03-joint-vs-separate-optimization-deltas.csv"), tables = out, figures = character())
}

jss_run_07_from_public_results <- function(settings) {
  paths <- jss_misspec_paths(settings)
  config <- jss_misspec_config(settings, stage = "full")
  grid <- jss_misspec_grid(config)
  public_results <- file.path(settings$public_data_dir, "copula-misspecification", "results.csv")
  approval <- jss_misspec_validate_approved_public_bundle(public_results, grid, config)
  results <- approval$results
  approval$results <- NULL
  jss_misspec_validate_public_full_bundle(results, grid, config)
  results <- jss_misspec_add_deltas(results)
  summary <- jss_misspec_summary(results)
  selection_attempts <- jss_misspec_selection_attempts(results)
  selection <- jss_misspec_selection(results)
  selection_confusion <- jss_misspec_selection_confusion(selection_attempts, criterion = "aic")
  selection_failures <- jss_misspec_selection_failures(selection_attempts)
  paired_effects <- jss_misspec_paired_effects(results)
  warning_audit <- jss_misspec_warning_audit(results)
  jss_misspec_write_csv_atomic(grid, paths$grid)
  jss_misspec_write_csv_atomic(results, paths$results)
  jss_misspec_write_csv_atomic(summary, paths$summary)
  jss_misspec_write_csv_atomic(selection, paths$selection)
  jss_misspec_write_csv_atomic(selection_attempts, paths$selection_attempts)
  jss_misspec_write_csv_atomic(selection_confusion, paths$selection_confusion)
  jss_misspec_write_csv_atomic(selection_failures, paths$selection_failures)
  jss_misspec_write_csv_atomic(paired_effects, paths$paired_effects)
  jss_misspec_write_csv_atomic(warning_audit, paths$warning_audit)
  source_manifest <- jss_misspec_evidence_paths(public_results)[["execution_manifest"]]
  manifest <- utils::read.csv(source_manifest, stringsAsFactors = FALSE, check.names = FALSE)
  jss_misspec_write_csv_atomic(manifest, paths$execution_manifest)
  jss_misspec_write_paper_summary_heatmap(summary, paths$paper_summary_heatmap)
  jss_misspec_validate_evidence_bundle(paths$results, config)
  installed_results <- utils::read.csv(paths$results, stringsAsFactors = FALSE, check.names = FALSE)
  review <- jss_misspec_binding_review_gate(
    installed_results, grid = grid, paths = paths,
    context = "paper-public-derived", config = config
  )
  jss_misspec_write_csv_atomic(review[c("check", "status", "detail")], paths$review)
  jss_misspec_revalidate_approved_source(public_results, approval, config)
  list(module_id = "07-gamma-copula-misspecification", status = "regenerated",
    data = c(paths$grid, paths$results, paths$selection_attempts, paths$execution_manifest),
    tables = c(paths$summary, paths$selection, paths$selection_confusion, paths$selection_failures,
      paths$paired_effects, paths$warning_audit, paths$review),
    figures = paths$paper_summary_heatmap, approved_identity = approval)
}

jss_run_phase2_multivariate_benchmark <- function(settings) {
  configured <- Sys.getenv("GAMLSS_LONGITUDINAL_MVT_APPROVED_RUN_DIR", unset = "")
  source <- if (nzchar(configured)) configured else
    file.path(settings$public_data_dir, "multivariate-benchmark")
  if (!dir.exists(source)) {
    stop("Approved Module 09 snapshot bundle is missing: ", source, call. = FALSE)
  }
  attestation <- if (!is.null(settings$multivariate_benchmark_attestation)) {
    settings$multivariate_benchmark_attestation
  } else {
    mvt_phase2_snapshot_attestation_path()
  }
  signature <- if (!is.null(settings$multivariate_benchmark_signature)) {
    settings$multivariate_benchmark_signature
  } else mvt_phase2_snapshot_signature_path()
  source_identity <- mvt_validate_phase2_claim_evidence(source, attestation, signature)
  destination <- file.path(settings$data_dir, "multivariate-benchmark")
  if (dir.exists(destination)) {
    installed_identity <- mvt_validate_phase2_claim_evidence(destination, attestation, signature)
    if (!identical(source_identity, installed_identity)) {
      stop("Existing Module 09 public integration does not match the approved source snapshot.", call. = FALSE)
    }
    integration <- c(installed_identity, list(integration_dir = destination))
  } else {
    integration <- mvt_integrate_approved_phase2_snapshot(source, destination, attestation, signature)
  }
  allowlist <- mvt_phase2_public_output_allowlist()
  attempt <- file.path(destination, allowlist$attempt_artifacts)
  evidence <- file.path(destination, allowlist$evidence_artifacts)
  if (!all(file.exists(c(attempt, evidence)))) {
    stop("Approved Module 09 integration is missing an exact allowlisted artifact.", call. = FALSE)
  }
  list(
    module_id = "09-simulation-multivariate-longitudinal", status = "current",
    data = attempt, tables = evidence, figures = character(),
    approved_identity = integration,
    notes = "Exact approved t20 x four-family x four-dependence x R=100 benchmark snapshot."
  )
}

jss_legacy_run_08_public <- function(settings) {
  if (identical(settings$profile, "full")) return(jss_legacy_run_08_full(settings))
  old <- Sys.getenv("GAMLSS_LONGITUDINAL_JSS_CORR_MISSPEC_RUN_DIR", unset = NA_character_)
  on.exit(if (is.na(old)) Sys.unsetenv("GAMLSS_LONGITUDINAL_JSS_CORR_MISSPEC_RUN_DIR") else Sys.setenv(GAMLSS_LONGITUDINAL_JSS_CORR_MISSPEC_RUN_DIR = old), add = TRUE)
  Sys.setenv(GAMLSS_LONGITUDINAL_JSS_CORR_MISSPEC_RUN_DIR = file.path(settings$public_data_dir, "correlation-misspecification"))
  jss_run_08_simulation_sensitivity_correlation_misspecification(settings)
}

jss_legacy_run_08_full <- function(settings) {
  scripts <- file.path(settings$root, "paper", "R", "08-simulation-sensitivity-correlation-misspecification", "standard-model-benchmarking")
  output_root <- file.path(settings$data_dir, "correlation-misspecification-full")
  base_env <- c(
    GAMLSS_LONGITUDINAL_BENCHMARK_OUTPUT_ROOT = output_root,
    GAMLSS_LONGITUDINAL_BENCHMARK_SOURCE = "installed",
    GAMLSS_LONGITUDINAL_BENCHMARK_FAMILIES = "gaussian,gamma,binary",
    GAMLSS_LONGITUDINAL_BENCHMARK_TIMEPOINTS = "20,50",
    GAMLSS_LONGITUDINAL_BENCHMARK_REPS = "20",
    GAMLSS_LONGITUDINAL_BENCHMARK_N = "120",
    GAMLSS_LONGITUDINAL_BENCHMARK_PRIMARY_TIMEOUT_SEC = "300",
    GAMLSS_LONGITUDINAL_BENCHMARK_GEE_TIMEOUT_SEC = "30",
    GAMLSS_LONGITUDINAL_BENCHMARK_GEE_UNSTRUCTURED_TIMEOUT_SEC = "30"
  )
  jss_run_script(file.path(scripts, "02-run-rs-joint-standard-model-grid.R"), base_env, settings$root)
  latest <- trimws(readLines(file.path(output_root, "latest_run_dir.txt"), warn = FALSE)[1L])
  follow_env <- c(base_env, GAMLSS_LONGITUDINAL_BENCHMARK_REPORT_DIR = latest)
  for (name in c("03-summarise-rs-joint-standard-model-grid.R", "05-add-gee-unstructured-rows.R", "08-recompute-all-pair-dependence.R", "11-run-t20-sandwich-grid.R", "07-write-story-tables.R")) {
    jss_run_script(file.path(scripts, name), follow_env, settings$root)
  }
  old <- Sys.getenv("GAMLSS_LONGITUDINAL_JSS_CORR_MISSPEC_RUN_DIR", unset = NA_character_)
  on.exit(if (is.na(old)) Sys.unsetenv("GAMLSS_LONGITUDINAL_JSS_CORR_MISSPEC_RUN_DIR") else Sys.setenv(GAMLSS_LONGITUDINAL_JSS_CORR_MISSPEC_RUN_DIR = old), add = TRUE)
  Sys.setenv(GAMLSS_LONGITUDINAL_JSS_CORR_MISSPEC_RUN_DIR = latest)
  jss_run_08_simulation_sensitivity_correlation_misspecification(settings)
}
