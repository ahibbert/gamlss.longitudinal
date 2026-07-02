args <- commandArgs(trailingOnly = TRUE)

default_run_dir <- file.path(
  "results",
  "jss-exploratory",
  "03-joint-vs-separate-optimization",
  "six_case_full_n_10rep_heldout_vs100"
)
run_dir <- if (length(args) >= 1L && nzchar(args[[1L]])) args[[1L]] else default_run_dir
summary_path <- if (length(args) >= 2L && nzchar(args[[2L]])) {
  args[[2L]]
} else {
  file.path("tables", "six-case-full-n-10rep-heldout-vs100-summary.csv")
}
summary_path <- if (grepl("^[A-Za-z]:|^/", summary_path)) {
  summary_path
} else {
  file.path(run_dir, summary_path)
}
out_path <- if (length(args) >= 3L && nzchar(args[[3L]])) {
  args[[3L]]
} else {
  file.path(run_dir, "tables", "six-case-joint-vs-separate-small-sample-table.tex")
}
table_label <- if (length(args) >= 4L && nzchar(args[[4L]])) {
  args[[4L]]
} else {
  "tab:jvs-small-sample"
}
caption_prefix <- if (length(args) >= 5L && nzchar(args[[5L]])) {
  args[[5L]]
} else {
  "Small-sample joint versus separate optimisation results."
}
deltas_path <- file.path(run_dir, "data", "03-joint-vs-separate-optimization-deltas.csv")

if (!file.exists(summary_path)) {
  fallback_summary_paths <- file.path(
    run_dir,
    "tables",
    c(
      "family-calibration-summary.csv",
      "03-joint-vs-separate-optimization-summary.csv"
    )
  )
  fallback_summary_paths <- fallback_summary_paths[file.exists(fallback_summary_paths)]
  if (length(args) < 2L && length(fallback_summary_paths)) {
    summary_path <- fallback_summary_paths[[1L]]
  } else {
    stop("Summary CSV not found: ", summary_path)
  }
}
if (!file.exists(deltas_path)) {
  stop("Delta CSV not found: ", deltas_path)
}

fmt <- function(x) {
  x <- suppressWarnings(as.numeric(x))
  ifelse(is.finite(x), sprintf("%.3f", x), "--")
}

fmt1 <- function(x) {
  x <- suppressWarnings(as.numeric(x))
  x[is.finite(x) & abs(x) < 0.05] <- 0
  ifelse(is.finite(x), sprintf("%.1f", x), "--")
}

fmt_median_iqr <- function(median, q1, q3) {
  paste0(fmt1(median), " [", fmt1(q1), ", ", fmt1(q3), "]")
}

escape_latex <- function(x) {
  x <- gsub("\\\\", "\\\\textbackslash{}", x)
  x <- gsub("([#$%&_{}])", "\\\\\\1", x, perl = TRUE)
  x
}

summary <- utils::read.csv(summary_path, stringsAsFactors = FALSE)
summary <- summary[summary$case_id %in% sprintf("JVS%02d", 1:6), , drop = FALSE]
if (!nrow(summary)) {
  stop("No JVS01-JVS06 rows found in: ", summary_path)
}
summary$case_order <- match(summary$case_id, sprintf("JVS%02d", 1:6))
summary <- summary[order(summary$case_order), , drop = FALSE]
summary <- unique(summary[, c("case_id", "case_order"), drop = FALSE])

deltas <- utils::read.csv(deltas_path, stringsAsFactors = FALSE)
deltas <- deltas[deltas$case_id %in% summary$case_id, , drop = FALSE]
if (!nrow(deltas)) {
  stop("No matching JVS01-JVS06 rows found in: ", deltas_path)
}

metrics <- c(
  "delta_train_joint_loglik",
  "delta_heldout_variogram_score_p05",
  "delta_heldout_variogram_score_p2"
)
metric_stats <- do.call(rbind, lapply(split(deltas, deltas$case_id), function(x) {
  out <- data.frame(case_id = x$case_id[[1L]], stringsAsFactors = FALSE)
  for (metric in metrics) {
    values <- suppressWarnings(as.numeric(x[[metric]]))
    out[[paste0(metric, "_mean")]] <- mean(values, na.rm = TRUE)
    out[[paste0(metric, "_sd")]] <- stats::sd(values, na.rm = TRUE)
    qs <- stats::quantile(values, probs = c(0.25, 0.5, 0.75), na.rm = TRUE, names = FALSE)
    out[[paste0(metric, "_q1")]] <- qs[[1L]]
    out[[paste0(metric, "_median")]] <- qs[[2L]]
    out[[paste0(metric, "_q3")]] <- qs[[3L]]
  }
  out
}))
summary <- merge(
  summary,
  metric_stats,
  by = "case_id",
  all.x = TRUE,
  sort = FALSE
)
summary <- summary[order(summary$case_order), , drop = FALSE]

caption <- paste(
  caption_prefix,
  "Values are median joint-minus-separate improvements with interquartile",
  "ranges [Q1, Q3]. Positive train log-likelihood values and negative",
  "variogram score values favour joint optimisation."
)

lines <- c(
  "\\begin{table}[ht]",
  "\\centering",
  paste0("\\caption{", caption, "}"),
  paste0("\\label{", table_label, "}"),
  "\\begin{tabular}{lrrr}",
  "\\toprule",
  "Case & Train LL & Held-out VS, \\(p = 0.5\\) & Held-out VS, \\(p = 2\\) \\\\",
  "\\midrule"
)

for (i in seq_len(nrow(summary))) {
  lines <- c(
    lines,
    paste(
      escape_latex(summary$case_id[[i]]),
      "&",
      fmt_median_iqr(
        summary$delta_train_joint_loglik_median[[i]],
        summary$delta_train_joint_loglik_q1[[i]],
        summary$delta_train_joint_loglik_q3[[i]]
      ),
      "&",
      fmt_median_iqr(
        summary$delta_heldout_variogram_score_p05_median[[i]],
        summary$delta_heldout_variogram_score_p05_q1[[i]],
        summary$delta_heldout_variogram_score_p05_q3[[i]]
      ),
      "&",
      fmt_median_iqr(
        summary$delta_heldout_variogram_score_p2_median[[i]],
        summary$delta_heldout_variogram_score_p2_q1[[i]],
        summary$delta_heldout_variogram_score_p2_q3[[i]]
      ),
      "\\\\"
    )
  )
}

lines <- c(
  lines,
  "\\bottomrule",
  "\\end{tabular}",
  "\\end{table}",
  ""
)

dir.create(dirname(out_path), recursive = TRUE, showWarnings = FALSE)
writeLines(lines, out_path)
cat("Wrote ", normalizePath(out_path, winslash = "/", mustWork = FALSE), "\n", sep = "")
