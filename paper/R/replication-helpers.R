jss_settings <- function() {
  profile <- Sys.getenv("GAMLSS_LONGITUDINAL_JSS_PROFILE", unset = "smoke")
  profile <- match.arg(profile, c("smoke", "expanded"))
  root <- Sys.getenv("GAMLSS_LONGITUDINAL_JSS_ROOT", unset = getwd())
  out_dir <- file.path(root, "results", "jss-replication", profile)
  dirs <- file.path(out_dir, c("data", "tables", "figures", "logs"))
  invisible(lapply(dirs, dir.create, recursive = TRUE, showWarnings = FALSE))
  list(
    profile = profile,
    root = root,
    out_dir = out_dir,
    data_dir = file.path(out_dir, "data"),
    tables_dir = file.path(out_dir, "tables"),
    figures_dir = file.path(out_dir, "figures"),
    logs_dir = file.path(out_dir, "logs"),
    seed = 20260528L
  )
}

jss_methods <- function(profile) {
  methods <- c("gamlss", "rs_separate", "rs_joint", "cg")
  if (identical(profile, "expanded") && requireNamespace("gamlss2", quietly = TRUE)) {
    methods <- c("gamlss2", methods)
  }
  methods
}

jss_grid_settings <- function(profile) {
  if (identical(profile, "expanded")) {
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
  results <- gamlss.longitudinal::run_coverage_simulations(
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
  )
  saveRDS(results, file.path(settings$data_dir, "coverage_results.rds"))
  results
}

jss_write_coverage_summary <- function(results, settings) {
  group_cols <- c("family", "copula", "design", "method", "success", "failure_type")
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
  template$output_path <- gsub("<profile>", settings$profile, template$output_path, fixed = TRUE)
  template$generated_at <- format(Sys.time(), "%Y-%m-%dT%H:%M:%S%z")
  template$profile_run <- settings$profile
  path <- file.path(settings$out_dir, "manifest.csv")
  utils::write.csv(template, path, row.names = FALSE)
  path
}

jss_write_output_hashes <- function(settings, manifest_path, ...) {
  explicit_files <- unlist(list(...), use.names = FALSE)
  files <- c(explicit_files, manifest_path)
  if (!length(explicit_files)) {
    files <- c(
      list.files(settings$tables_dir, full.names = TRUE, recursive = TRUE),
      list.files(settings$figures_dir, full.names = TRUE, recursive = TRUE),
      file.path(settings$logs_dir, "session_info.txt"),
      manifest_path
    )
  }
  files <- files[file.exists(files)]
  hashes <- data.frame(
    file = normalizePath(files, winslash = "/", mustWork = TRUE),
    md5 = unname(tools::md5sum(files)),
    stringsAsFactors = FALSE
  )
  path <- file.path(settings$logs_dir, "output_hashes.csv")
  utils::write.csv(hashes, path, row.names = FALSE)
  path
}

jss_validate_manifest <- function(manifest_path, output_hashes_path) {
  manifest <- utils::read.csv(manifest_path, stringsAsFactors = FALSE)
  expected <- manifest$output_path
  missing <- expected[!file.exists(expected)]
  if (length(missing) > 0) {
    stop("Missing JSS replication output(s): ", paste(missing, collapse = ", "), call. = FALSE)
  }
  if (!file.exists(output_hashes_path)) {
    stop("Missing output hash log: ", output_hashes_path, call. = FALSE)
  }
  TRUE
}
