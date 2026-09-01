source(file.path("paper", "R", "09-simulation-multivariate-longitudinal", "00-multivariate-setup.R"))

mvt_load_package()
run_dir <- mvt_read_run_dir()
tables_dir <- file.path(run_dir, "paper_tables")
dir.create(tables_dir, recursive = TRUE, showWarnings = FALSE)

summary_paths <- file.path(
  run_dir,
  c(
    "benchmark_summary.csv",
    "coefficient_summary.csv",
    "dependence_recovery_summary.csv",
    "variogram_summary.csv"
  )
)
if (any(!file.exists(summary_paths))) {
  mvt_summarise_results(run_dir)
}

latex_escape <- function(x) {
  x <- as.character(x)
  x <- gsub("\\\\", "\\\\textbackslash{}", x)
  x <- gsub("([&_#%$])", "\\\\\\1", x, perl = TRUE)
  x <- gsub("~", "\\\\textasciitilde{}", x, fixed = TRUE)
  x <- gsub("\\^", "\\\\textasciicircum{}", x)
  x
}

format_table_value <- function(x) {
  if (is.numeric(x)) {
    out <- ifelse(is.finite(x), formatC(x, format = "fg", digits = 4), "")
    return(out)
  }
  x[is.na(x)] <- ""
  as.character(x)
}

write_latex_tabular <- function(dat, path, max_rows = 60L) {
  dat <- head(dat, max_rows)
  if (nrow(dat) == 0L || ncol(dat) == 0L) {
    writeLines("% Empty table", path, useBytes = TRUE)
    return(invisible(path))
  }
  display <- as.data.frame(lapply(dat, format_table_value), stringsAsFactors = FALSE)
  display[] <- lapply(display, latex_escape)
  align <- paste0("l", paste(rep("r", max(0L, ncol(display) - 1L)), collapse = ""))
  header <- paste(latex_escape(names(display)), collapse = " & ")
  body <- apply(display, 1L, paste, collapse = " & ")
  lines <- c(
    paste0("\\begin{tabular}{", align, "}"),
    "\\hline",
    paste0(header, " \\\\"),
    "\\hline",
    paste0(body, " \\\\"),
    "\\hline",
    "\\end{tabular}"
  )
  writeLines(lines, path, useBytes = TRUE)
  invisible(path)
}

copy_table <- function(stem) {
  src <- file.path(run_dir, paste0(stem, ".csv"))
  if (!file.exists(src)) return(invisible(NULL))
  dat <- mvt_read_optional_csv(src)
  out_csv <- file.path(tables_dir, paste0(stem, ".csv"))
  mvt_write_csv(dat, out_csv)
  out_tex <- file.path(tables_dir, paste0(stem, ".tex"))
  write_latex_tabular(dat, out_tex)
  invisible(out_csv)
}

for (stem in c("benchmark_summary", "coefficient_summary", "dependence_recovery_summary", "variogram_summary")) {
  copy_table(stem)
}

read_run_table <- function(stem) {
  mvt_read_optional_csv(file.path(run_dir, paste0(stem, ".csv")))
}

wide_metric_summary <- function(dat, id_cols) {
  if (nrow(dat) == 0L || !"metric" %in% names(dat) || !"mean" %in% names(dat)) return(data.frame())
  keep <- intersect(c(id_cols, "metric", "mean", "n_finite"), names(dat))
  dat <- dat[, keep, drop = FALSE]
  value_wide <- stats::reshape(
    dat[, intersect(c(id_cols, "metric", "mean"), names(dat)), drop = FALSE],
    idvar = id_cols,
    timevar = "metric",
    direction = "wide"
  )
  names(value_wide) <- sub("^mean[.]", "", names(value_wide))
  if (!"n_finite" %in% names(dat)) return(value_wide)
  n_wide <- stats::reshape(
    dat[, intersect(c(id_cols, "metric", "n_finite"), names(dat)), drop = FALSE],
    idvar = id_cols,
    timevar = "metric",
    direction = "wide"
  )
  names(n_wide) <- sub("^n_finite[.]", "n_", names(n_wide))
  merge(value_wide, n_wide, by = id_cols, all = TRUE)
}

write_primary_table <- function(dat, stem) {
  out_csv <- file.path(tables_dir, paste0(stem, ".csv"))
  out_tex <- file.path(tables_dir, paste0(stem, ".tex"))
  mvt_write_csv(dat, out_csv)
  write_latex_tabular(dat, out_tex, max_rows = 80L)
  invisible(out_csv)
}

bench_summary <- read_run_table("benchmark_summary")
vario_summary <- read_run_table("variogram_summary")
coef_rep <- read_run_table("coefficient_results_by_rep")
dep_summary <- read_run_table("dependence_recovery_summary")

distribution_fit <- wide_metric_summary(
  bench_summary[bench_summary$metric %in% c("benchmark_neg_log_score", "logLik", "logLik_df", "AIC", "BIC"), , drop = FALSE],
  id_cols = intersect(c("scenario", "n_time", "family", "method"), names(bench_summary))
)
vario_fit <- wide_metric_summary(
  vario_summary[vario_summary$metric %in% c("variogram_score_p05", "variogram_score_p2"), , drop = FALSE],
  id_cols = intersect(c("scenario", "n_time", "family", "method"), names(vario_summary))
)
if (nrow(distribution_fit) > 0L && nrow(vario_fit) > 0L) {
  distribution_fit <- merge(distribution_fit, vario_fit, by = intersect(c("scenario", "n_time", "family", "method"), names(distribution_fit)), all = TRUE)
}
write_primary_table(distribution_fit, "primary_distributional_fit")

if (nrow(coef_rep) > 0L && all(c("scenario", "n_time", "family", "method", "term", "bias", "estimate", "std_error", "ci_covers_truth") %in% names(coef_rep))) {
  coef_rep <- coef_rep[coef_rep$term %in% c("time", "x", "z"), , drop = FALSE]
  groups <- unique(coef_rep[c("scenario", "n_time", "family", "method", "term")])
  coef_rows <- lapply(seq_len(nrow(groups)), function(i) {
    idx <- rep(TRUE, nrow(coef_rep))
    for (col in names(groups)) idx <- idx & coef_rep[[col]] == groups[[col]][i]
    sub <- coef_rep[idx, , drop = FALSE]
    est <- suppressWarnings(as.numeric(sub$estimate))
    se <- suppressWarnings(as.numeric(sub$std_error))
    bias <- suppressWarnings(as.numeric(sub$bias))
    data.frame(
      groups[i, , drop = FALSE],
      marginal_bias = mean(bias, na.rm = TRUE),
      marginal_rmse = sqrt(mean(bias^2, na.rm = TRUE)),
      marginal_mae = mean(abs(bias), na.rm = TRUE),
      empirical_sd = stats::sd(est, na.rm = TRUE),
      mean_reported_se = mean(se, na.rm = TRUE),
      se_ratio = mean(se, na.rm = TRUE) / stats::sd(est, na.rm = TRUE),
      ci_coverage = mean(sub$ci_covers_truth %in% c(TRUE, "TRUE", "true", "1"), na.rm = TRUE),
      mean_ci_width = mean(suppressWarnings(as.numeric(sub$ci_width)), na.rm = TRUE),
      n_reps = nrow(sub),
      stringsAsFactors = FALSE
    )
  })
  coef_primary <- mvt_bind_rows_fill(coef_rows)
} else {
  coef_primary <- data.frame()
}
write_primary_table(coef_primary, "primary_marginal_covariate_recovery")
write_primary_table(coef_primary[, intersect(c("scenario", "n_time", "family", "method", "term", "empirical_sd", "mean_reported_se", "se_ratio", "ci_coverage", "mean_ci_width", "n_reps"), names(coef_primary)), drop = FALSE], "primary_marginal_uncertainty")

dependence_primary <- wide_metric_summary(
  dep_summary[dep_summary$metric %in% c("theta_rmse", "tau_rmse", "theta_mae", "tau_mae"), , drop = FALSE],
  id_cols = intersect(c("scenario", "n_time", "family", "method"), names(dep_summary))
)
write_primary_table(dependence_primary, "primary_dependence_recovery")

message("Paper table fragments written to: ", tables_dir)
