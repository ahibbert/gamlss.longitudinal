# Lean supplementary publication layer for the RS joint-versus-separate study.
#
# The author-facing checkpoint runner remains in 03-joint-vs-separate-optimization.R.
# This file deliberately exports only reviewer-relevant attempt-level evidence.
# It does not expose or validate signatures, hash chains, checkout paths, hosts,
# process identifiers, or other machine-specific execution metadata.

jss_optimizer_lean_version <- function() "optimizer-publication-lean-v1"

jss_optimizer_lean_cases <- function() {
  data.frame(
    case_id = c("JVS01", "JVS02", "JVS03"),
    reference_case_id = c("JVS01", "JVS01", "JVS01"),
    case_role = c("base", "dependence_strength", "scale_complexity"),
    family = "NO",
    copula = "N",
    n_subject = 120L,
    n_time = 8L,
    changed_factor = c("none", "theta_strength", "sigma_strength"),
    theta_strength = c(1.5, 3.5, 1.5),
    sigma_strength = c(1.5, 1.5, 0.5),
    publication_scope = "supplementary",
    stringsAsFactors = FALSE
  )
}

jss_optimizer_lean_precision <- function() {
  data.frame(
    schema_version = jss_optimizer_lean_version(),
    target_usable_pairs = 97L,
    confidence_level = 0.95,
    estimand = "paired RS-joint minus RS-separate mean difference",
    runtime_interpretation = "descriptive_order_confounded",
    publication_scope = "supplementary",
    stringsAsFactors = FALSE
  )
}

jss_optimizer_lean_metrics <- function() {
  data.frame(
    metric = c(
      "convergence_rate", "elapsed_sec", "train_joint_loglik", "test_log_score_per_obs",
      "heldout_variogram_score_p05", "train_rmse_mu", "train_rmse_sigma",
      "train_rmse_theta", "train_rmse_tau"
    ),
    direction = c("higher", "descriptive", "higher", "higher", "lower", "lower", "lower", "lower", "lower"),
    interpretation = c(
      "paired_attempt_convergence", "descriptive_order_confounded",
      rep("conditional_on_usable_matched_pairs", 7L)
    ),
    requires_retained_pair = c(FALSE, rep(TRUE, 8L)),
    stringsAsFactors = FALSE
  )
}

jss_optimizer_lean_files <- function() {
  list(
    data = c(
      "03-joint-vs-separate-optimization-candidate-selection.csv",
      "03-joint-vs-separate-optimization-mc-precision.csv",
      "03-joint-vs-separate-optimization-results.csv",
      "03-joint-vs-separate-optimization-paired-differences.csv"
    ),
    tables = c(
      "03-joint-vs-separate-optimization-summary.csv",
      "03-joint-vs-separate-optimization-failures.csv",
      "03-joint-vs-separate-optimization-difference-uncertainty.csv"
    ),
    figures = "03-joint-vs-separate-optimization-deltas.png"
  )
}

jss_optimizer_lean_attempt_names <- function() {
  c(
    "schema_version", "case_id", "case_role", "joint_review_rep", "paired_seed",
    "method", "success", "converged", "retained", "failure_reason", "stop_reason",
    "outer_iterations", "elapsed_sec", "train_joint_loglik", "test_log_score_per_obs",
    "heldout_variogram_score_p05", "train_rmse_mu", "train_rmse_sigma",
    "train_rmse_theta", "train_rmse_tau", "publication_scope"
  )
}

jss_optimizer_lean_pair_names <- function() {
  c(
    "schema_version", "case_id", "case_role", "joint_review_rep", "paired_seed",
    "rs_separate_success", "rs_separate_converged", "rs_joint_success",
    "rs_joint_converged", "pair_retained", "pair_failure_reason",
    paste0("delta_", jss_optimizer_lean_metrics()$metric),
    "runtime_interpretation", "publication_scope"
  )
}

jss_optimizer_lean_summary_names <- function() {
  c("schema_version", "case_id", "case_role", "method", "attempted", "converged",
    "retained", "failed", "publication_scope")
}

jss_optimizer_lean_failure_names <- function() {
  c("schema_version", "case_id", "case_role", "method", "failure_reason",
    "failure_count", "publication_scope")
}

jss_optimizer_lean_uncertainty_names <- function() {
  c(
    "schema_version", "comparison_id", "comparison_type", "focal_case_id",
    "reference_case_id", "metric", "interpretation", "attempted_pairs",
    "usable_pairs", "failed_pairs", "estimate", "mcse", "conf_low", "conf_high",
    "confidence_level", "publication_scope"
  )
}

jss_optimizer_lean_bind <- function(rows) {
  if (!length(rows)) return(data.frame())
  out <- do.call(rbind, rows)
  rownames(out) <- NULL
  out
}

jss_optimizer_lean_failure_reason <- function(x) {
  retained <- x$success %in% TRUE & x$converged %in% TRUE
  reason <- as.character(x$failure_type)
  fallback <- as.character(x$error)
  use <- is.na(reason) | !nzchar(trimws(reason)) | reason == "none"
  reason[use] <- fallback[use]
  use <- is.na(reason) | !nzchar(trimws(reason)) | reason == "none"
  reason[use] <- as.character(x$stop_reason[use])
  use <- is.na(reason) | !nzchar(trimws(reason)) | reason == "none"
  reason[use] <- "unclassified_failure"
  reason[retained] <- "none"
  reason
}

jss_optimizer_lean_attempts <- function(raw) {
  required <- c(
    "case_id", "joint_review_rep", "paired_seed", "method", "success", "converged",
    "stop_reason", "outer_iterations", "elapsed_sec", "train_joint_loglik",
    "test_log_score_per_obs", "heldout_variogram_score_p05", "train_rmse_mu",
    "train_rmse_sigma", "train_rmse_theta", "train_rmse_tau", "failure_type", "error"
  )
  if (!is.data.frame(raw) || !all(required %in% names(raw))) {
    stop("Optimizer checkpoint rows lack the lean attempt fields.", call. = FALSE)
  }
  cases <- jss_optimizer_lean_cases()
  raw <- raw[raw$case_id %in% cases$case_id, , drop = FALSE]
  if (!nrow(raw)) stop("No authoritative optimizer attempt rows were found.", call. = FALSE)
  case_index <- match(raw$case_id, cases$case_id)
  out <- data.frame(
    schema_version = jss_optimizer_lean_version(),
    case_id = as.character(raw$case_id),
    case_role = cases$case_role[case_index],
    joint_review_rep = as.integer(raw$joint_review_rep),
    paired_seed = as.integer(raw$paired_seed),
    method = as.character(raw$method),
    success = raw$success %in% TRUE,
    converged = raw$converged %in% TRUE,
    retained = raw$success %in% TRUE & raw$converged %in% TRUE,
    failure_reason = jss_optimizer_lean_failure_reason(raw),
    stop_reason = as.character(raw$stop_reason),
    outer_iterations = as.integer(raw$outer_iterations),
    elapsed_sec = as.numeric(raw$elapsed_sec),
    train_joint_loglik = as.numeric(raw$train_joint_loglik),
    test_log_score_per_obs = as.numeric(raw$test_log_score_per_obs),
    heldout_variogram_score_p05 = as.numeric(raw$heldout_variogram_score_p05),
    train_rmse_mu = as.numeric(raw$train_rmse_mu),
    train_rmse_sigma = as.numeric(raw$train_rmse_sigma),
    train_rmse_theta = as.numeric(raw$train_rmse_theta),
    train_rmse_tau = as.numeric(raw$train_rmse_tau),
    publication_scope = "supplementary",
    stringsAsFactors = FALSE
  )
  out <- out[order(match(out$case_id, cases$case_id), out$joint_review_rep,
    match(out$method, c("rs_separate", "rs_joint"))), , drop = FALSE]
  rownames(out) <- NULL
  out
}

jss_optimizer_lean_read_checkpoints <- function(checkpoint_dir) {
  files <- list.files(checkpoint_dir, pattern = "[.]rds$", recursive = TRUE, full.names = TRUE)
  if (!length(files)) stop("Optimizer checkpoint directory contains no RDS files.", call. = FALSE)
  payloads <- lapply(files, function(path) {
    value <- tryCatch(readRDS(path), error = function(e) NULL)
    if (!is.list(value) || !is.data.frame(value$result)) {
      stop("Optimizer checkpoint is unreadable or lacks result rows: ", basename(path), call. = FALSE)
    }
    value$result
  })
  jss_optimizer_lean_attempts(jss_optimizer_lean_bind(payloads))
}

jss_optimizer_lean_pairs <- function(attempts) {
  groups <- split(attempts, paste(attempts$case_id, attempts$joint_review_rep, sep = "\r"))
  metrics <- jss_optimizer_lean_metrics()$metric
  rows <- lapply(groups, function(x) {
    x <- x[match(c("rs_separate", "rs_joint"), x$method), , drop = FALSE]
    if (nrow(x) != 2L || anyNA(x$method) ||
        !identical(as.character(x$method), c("rs_separate", "rs_joint"))) {
      stop("Optimizer replicate is missing one member of its method pair.", call. = FALSE)
    }
    if (x$paired_seed[[1L]] != x$paired_seed[[2L]]) {
      stop("Optimizer method pair does not share one paired seed.", call. = FALSE)
    }
    pair_retained <- all(x$retained)
    metric_contract <- jss_optimizer_lean_metrics()
    deltas <- vapply(seq_len(nrow(metric_contract)), function(i) {
      metric <- metric_contract$metric[[i]]
      values <- if (metric == "convergence_rate") as.numeric(x$converged) else as.numeric(x[[metric]])
      eligible <- !metric_contract$requires_retained_pair[[i]] || pair_retained
      if (eligible && all(is.finite(values))) values[[2L]] - values[[1L]] else NA_real_
    }, numeric(1L))
    names(deltas) <- metric_contract$metric
    failures <- if (pair_retained) "none" else paste(
      paste0(x$method[!x$retained], ":", x$failure_reason[!x$retained]), collapse = ";"
    )
    row <- data.frame(
      schema_version = jss_optimizer_lean_version(),
      case_id = x$case_id[[1L]], case_role = x$case_role[[1L]],
      joint_review_rep = x$joint_review_rep[[1L]], paired_seed = x$paired_seed[[1L]],
      rs_separate_success = x$success[[1L]], rs_separate_converged = x$converged[[1L]],
      rs_joint_success = x$success[[2L]], rs_joint_converged = x$converged[[2L]],
      pair_retained = pair_retained, pair_failure_reason = failures,
      runtime_interpretation = "descriptive_order_confounded",
      publication_scope = "supplementary", stringsAsFactors = FALSE
    )
    for (metric in metrics) row[[paste0("delta_", metric)]] <- deltas[[metric]]
    row[, jss_optimizer_lean_pair_names(), drop = FALSE]
  })
  out <- jss_optimizer_lean_bind(rows)
  cases <- jss_optimizer_lean_cases()
  out <- out[order(match(out$case_id, cases$case_id), out$joint_review_rep), , drop = FALSE]
  rownames(out) <- NULL
  out
}

jss_optimizer_lean_summary <- function(attempts) {
  cases <- jss_optimizer_lean_cases()
  groups <- split(attempts, paste(attempts$case_id, attempts$method, sep = "\r"))
  rows <- lapply(groups, function(x) data.frame(
    schema_version = jss_optimizer_lean_version(), case_id = x$case_id[[1L]],
    case_role = x$case_role[[1L]], method = x$method[[1L]], attempted = nrow(x),
    converged = sum(x$converged), retained = sum(x$retained), failed = sum(!x$retained),
    publication_scope = "supplementary", stringsAsFactors = FALSE
  ))
  out <- jss_optimizer_lean_bind(rows)
  out <- out[order(match(out$case_id, cases$case_id),
    match(out$method, c("rs_separate", "rs_joint"))), , drop = FALSE]
  rownames(out) <- NULL
  out
}

jss_optimizer_lean_failures <- function(attempts) {
  cases <- jss_optimizer_lean_cases()
  groups <- split(attempts, paste(attempts$case_id, attempts$method, sep = "\r"))
  rows <- unlist(lapply(groups, function(x) {
    failed <- x$failure_reason[!x$retained]
    counts <- table(failed)
    if (!length(counts)) counts <- structure(0L, names = "none")
    lapply(seq_along(counts), function(i) data.frame(
      schema_version = jss_optimizer_lean_version(), case_id = x$case_id[[1L]],
      case_role = x$case_role[[1L]], method = x$method[[1L]],
      failure_reason = names(counts)[[i]], failure_count = as.integer(counts[[i]]),
      publication_scope = "supplementary", stringsAsFactors = FALSE
    ))
  }), recursive = FALSE)
  out <- jss_optimizer_lean_bind(rows)
  out <- out[order(match(out$case_id, cases$case_id),
    match(out$method, c("rs_separate", "rs_joint")), out$failure_reason), , drop = FALSE]
  rownames(out) <- NULL
  out
}

jss_optimizer_lean_interval <- function(values, confidence) {
  values <- values[is.finite(values)]
  n <- length(values)
  estimate <- if (n) mean(values) else NA_real_
  mcse <- if (n > 1L) stats::sd(values) / sqrt(n) else NA_real_
  half <- if (is.finite(mcse)) stats::qnorm((1 + confidence) / 2) * mcse else NA_real_
  c(estimate = estimate, mcse = mcse, conf_low = estimate - half, conf_high = estimate + half)
}

jss_optimizer_lean_uncertainty <- function(pairs) {
  precision <- jss_optimizer_lean_precision()
  metrics <- jss_optimizer_lean_metrics()
  cases <- jss_optimizer_lean_cases()
  rows <- list()
  add <- function(comparison_id, comparison_type, focal, reference, metric, values, attempted,
      interpretation) {
    interval <- jss_optimizer_lean_interval(values, precision$confidence_level[[1L]])
    data.frame(
      schema_version = jss_optimizer_lean_version(), comparison_id = comparison_id,
      comparison_type = comparison_type, focal_case_id = focal, reference_case_id = reference,
      metric = metric, interpretation = interpretation, attempted_pairs = as.integer(attempted),
      usable_pairs = as.integer(sum(is.finite(values))),
      failed_pairs = as.integer(attempted - sum(is.finite(values))),
      estimate = unname(interval[["estimate"]]), mcse = unname(interval[["mcse"]]),
      conf_low = unname(interval[["conf_low"]]), conf_high = unname(interval[["conf_high"]]),
      confidence_level = precision$confidence_level[[1L]], publication_scope = "supplementary",
      stringsAsFactors = FALSE
    )
  }
  for (case_id in cases$case_id) {
    x <- pairs[pairs$case_id == case_id, , drop = FALSE]
    for (i in seq_len(nrow(metrics))) {
      metric <- metrics$metric[[i]]
      rows[[length(rows) + 1L]] <- add(
        case_id, "paired_optimizer_difference", case_id, case_id, metric,
        x[[paste0("delta_", metric)]], nrow(x), metrics$interpretation[[i]]
      )
    }
  }
  base <- pairs[pairs$case_id == "JVS01", , drop = FALSE]
  for (focal in c("JVS02", "JVS03")) {
    x <- merge(
      pairs[pairs$case_id == focal, c("joint_review_rep", paste0("delta_", metrics$metric)), drop = FALSE],
      base[, c("joint_review_rep", paste0("delta_", metrics$metric)), drop = FALSE],
      by = "joint_review_rep", suffixes = c("_focal", "_base"), all = FALSE, sort = TRUE
    )
    for (i in seq_len(nrow(metrics))) {
      metric <- metrics$metric[[i]]
      values <- x[[paste0("delta_", metric, "_focal")]] - x[[paste0("delta_", metric, "_base")]]
      interpretation <- if (metric == "elapsed_sec") {
        "descriptive_order_confounded"
      } else if (metric == "convergence_rate") {
        "paired_attempt_convergence"
      } else {
        "difference_in_paired_differences"
      }
      rows[[length(rows) + 1L]] <- add(
        paste0(focal, "-JVS01"), "difference_in_differences", focal, "JVS01",
        metric, values, nrow(x), interpretation
      )
    }
  }
  out <- jss_optimizer_lean_bind(rows)
  out[, jss_optimizer_lean_uncertainty_names(), drop = FALSE]
}

jss_optimizer_lean_compare <- function(actual, expected, keys, label, tolerance = 1e-10) {
  if (!identical(names(actual), names(expected))) stop(label, " has a noncanonical schema.", call. = FALSE)
  order_frame <- function(x) {
    if (!nrow(x)) return(x)
    x[do.call(order, unname(lapply(x[keys], as.character))), , drop = FALSE]
  }
  actual <- order_frame(actual); expected <- order_frame(expected)
  rownames(actual) <- NULL; rownames(expected) <- NULL
  if (nrow(actual) != nrow(expected)) stop(label, " has a noncanonical row count.", call. = FALSE)
  for (name in names(expected)) {
    ok <- if (is.numeric(expected[[name]])) {
      isTRUE(all.equal(as.numeric(actual[[name]]), as.numeric(expected[[name]]),
        tolerance = tolerance, check.attributes = FALSE))
    } else identical(as.character(actual[[name]]), as.character(expected[[name]]))
    if (!ok) stop(label, " disagrees with attempts in column ", name, ".", call. = FALSE)
  }
  invisible(TRUE)
}

jss_optimizer_lean_validate_attempts <- function(attempts, design, precision) {
  if (!identical(names(attempts), jss_optimizer_lean_attempt_names())) {
    stop("Optimizer attempt CSV has a noncanonical schema.", call. = FALSE)
  }
  logical_fields <- c("success", "converged", "retained")
  integer_fields <- c("joint_review_rep", "paired_seed", "outer_iterations")
  numeric_fields <- setdiff(jss_optimizer_lean_metrics()$metric, "convergence_rate")
  if (!all(vapply(attempts[logical_fields], is.logical, logical(1L))) ||
      !all(vapply(attempts[integer_fields], is.integer, logical(1L)))) {
    stop("Optimizer attempt CSV violates logical or integer types.", call. = FALSE)
  }
  if (any(vapply(attempts[numeric_fields], function(x) any(is.infinite(x)), logical(1L))) ||
      any(attempts$elapsed_sec < 0, na.rm = TRUE) || any(attempts$outer_iterations < 0L, na.rm = TRUE)) {
    stop("Optimizer attempt CSV contains impossible numeric values.", call. = FALSE)
  }
  if (!setequal(unique(attempts$case_id), design$case_id) ||
      any(attempts$case_id == "JVS04") ||
      !setequal(unique(attempts$method), c("rs_separate", "rs_joint")) ||
      any(attempts$retained != (attempts$success & attempts$converged)) ||
      any(attempts$retained & attempts$failure_reason != "none") ||
      any(!attempts$retained & (is.na(attempts$failure_reason) | !nzchar(attempts$failure_reason) |
        attempts$failure_reason == "none"))) {
    stop("Optimizer attempt lifecycle, cases, methods, or failure reasons are invalid.", call. = FALSE)
  }
  key <- paste(attempts$case_id, attempts$joint_review_rep, attempts$method, sep = "\r")
  if (anyDuplicated(key)) stop("Optimizer attempt keys are duplicated.", call. = FALSE)
  groups <- split(attempts, attempts$case_id)
  target <- precision$target_usable_pairs[[1L]]
  complete <- vapply(groups, function(x) {
    reps <- sort(unique(x$joint_review_rep))
    identical(reps, seq_len(max(reps))) && all(table(x$joint_review_rep) == 2L) &&
      sum(tapply(x$retained, x$joint_review_rep, all)) >= target
  }, logical(1L))
  if (length(groups) != nrow(design) || any(!complete)) {
    stop("Optimizer cases lack complete paired attempts at the registered precision target.", call. = FALSE)
  }
  invisible(TRUE)
}

jss_optimizer_lean_png_ok <- function(path) {
  if (!file.exists(path) || file.info(path)$size < 1000L) return(FALSE)
  bytes <- readBin(path, "raw", n = 24L)
  signature <- as.raw(c(137, 80, 78, 71, 13, 10, 26, 10))
  identical(bytes[seq_along(signature)], signature) && length(unique(bytes)) > 8L
}

jss_optimizer_lean_validate <- function(data_dir, tables_dir, figures_dir) {
  files <- jss_optimizer_lean_files()
  required <- c(file.path(data_dir, files$data), file.path(tables_dir, files$tables),
    file.path(figures_dir, files$figures))
  if (any(!file.exists(required))) {
    stop("Lean optimizer publication bundle is incomplete: ",
      paste(basename(required[!file.exists(required)]), collapse = ", "), call. = FALSE)
  }
  design <- utils::read.csv(file.path(data_dir, files$data[[1L]]), stringsAsFactors = FALSE)
  precision <- utils::read.csv(file.path(data_dir, files$data[[2L]]), stringsAsFactors = FALSE)
  attempts <- utils::read.csv(file.path(data_dir, files$data[[3L]]), stringsAsFactors = FALSE)
  pairs <- utils::read.csv(file.path(data_dir, files$data[[4L]]), stringsAsFactors = FALSE)
  summary <- utils::read.csv(file.path(tables_dir, files$tables[[1L]]), stringsAsFactors = FALSE)
  failures <- utils::read.csv(file.path(tables_dir, files$tables[[2L]]), stringsAsFactors = FALSE)
  uncertainty <- utils::read.csv(file.path(tables_dir, files$tables[[3L]]), stringsAsFactors = FALSE)
  jss_optimizer_lean_compare(design, jss_optimizer_lean_cases(), "case_id", "Optimizer design")
  jss_optimizer_lean_compare(precision, jss_optimizer_lean_precision(), "schema_version", "Optimizer precision")
  jss_optimizer_lean_validate_attempts(attempts, design, precision)
  expected_pairs <- jss_optimizer_lean_pairs(attempts)
  expected_summary <- jss_optimizer_lean_summary(attempts)
  expected_failures <- jss_optimizer_lean_failures(attempts)
  expected_uncertainty <- jss_optimizer_lean_uncertainty(expected_pairs)
  jss_optimizer_lean_compare(pairs, expected_pairs,
    c("case_id", "joint_review_rep"), "Optimizer paired differences")
  jss_optimizer_lean_compare(summary, expected_summary,
    c("case_id", "method"), "Optimizer lifecycle summary")
  jss_optimizer_lean_compare(failures, expected_failures,
    c("case_id", "method", "failure_reason"), "Optimizer failure summary")
  jss_optimizer_lean_compare(uncertainty, expected_uncertainty,
    c("comparison_id", "metric"), "Optimizer uncertainty")
  factor_rows <- uncertainty$comparison_type == "difference_in_differences"
  nonruntime <- uncertainty$metric != "elapsed_sec"
  if (any(uncertainty$usable_pairs[factor_rows & nonruntime] < precision$target_usable_pairs[[1L]]) ||
      any(uncertainty$interpretation[uncertainty$metric == "elapsed_sec"] !=
        "descriptive_order_confounded") || any(uncertainty$publication_scope != "supplementary")) {
    stop("Optimizer uncertainty lacks usable factor contrasts or supplementary runtime labeling.", call. = FALSE)
  }
  if (!jss_optimizer_lean_png_ok(file.path(figures_dir, files$figures[[1L]]))) {
    stop("Optimizer supplementary figure is missing, blank, or not a PNG.", call. = FALSE)
  }
  invisible(required)
}

jss_optimizer_lean_write_figure <- function(uncertainty, path) {
  colours <- jss_plot_colours()
  selected <- uncertainty[
    uncertainty$comparison_type == "paired_optimizer_difference" &
      uncertainty$metric %in% c("train_joint_loglik", "test_log_score_per_obs",
        "heldout_variogram_score_p05", "train_rmse_theta"), , drop = FALSE]
  metric_labels <- c(
    train_joint_loglik = "Training joint log-likelihood",
    test_log_score_per_obs = "Test log score per observation",
    heldout_variogram_score_p05 = "Held-out variogram score (p = 0.5)",
    train_rmse_theta = "Training dependence-parameter RMSE"
  )
  selected$metric <- factor(selected$metric, levels = names(metric_labels), labels = metric_labels)
  selected$focal_case_id <- factor(
    selected$focal_case_id,
    levels = unique(selected$focal_case_id)
  )
  p <- ggplot2::ggplot(
    selected,
    ggplot2::aes(.data$focal_case_id, .data$estimate)
  ) +
    ggplot2::geom_hline(
      yintercept = 0,
      linetype = "dashed",
      colour = colours[["reference"]],
      linewidth = 0.45
    ) +
    ggplot2::geom_errorbar(
      ggplot2::aes(ymin = .data$conf_low, ymax = .data$conf_high),
      width = 0.12,
      colour = colours[["observed"]]
    ) +
    ggplot2::geom_point(colour = colours[["fitted"]], size = 2.3) +
    ggplot2::facet_wrap(~metric, scales = "free_y", ncol = 2) +
    ggplot2::labs(
      x = "Scenario",
      y = "RS joint - RS separate"
    ) +
    theme_jss_paper(base_size = 9, legend_position = "none")
  dir.create(dirname(path), recursive = TRUE, showWarnings = FALSE)
  jss_save_figure(path, p, size = "standard")
  invisible(path)
}

jss_optimizer_lean_publish <- function(checkpoint_dir, data_dir, tables_dir, figures_dir) {
  attempts <- jss_optimizer_lean_read_checkpoints(checkpoint_dir)
  pairs <- jss_optimizer_lean_pairs(attempts)
  summary <- jss_optimizer_lean_summary(attempts)
  failures <- jss_optimizer_lean_failures(attempts)
  uncertainty <- jss_optimizer_lean_uncertainty(pairs)
  files <- jss_optimizer_lean_files()
  invisible(lapply(c(data_dir, tables_dir, figures_dir), dir.create,
    recursive = TRUE, showWarnings = FALSE))
  utils::write.csv(jss_optimizer_lean_cases(), file.path(data_dir, files$data[[1L]]), row.names = FALSE)
  utils::write.csv(jss_optimizer_lean_precision(), file.path(data_dir, files$data[[2L]]), row.names = FALSE)
  utils::write.csv(attempts, file.path(data_dir, files$data[[3L]]), row.names = FALSE, na = "")
  utils::write.csv(pairs, file.path(data_dir, files$data[[4L]]), row.names = FALSE, na = "")
  utils::write.csv(summary, file.path(tables_dir, files$tables[[1L]]), row.names = FALSE)
  utils::write.csv(failures, file.path(tables_dir, files$tables[[2L]]), row.names = FALSE)
  utils::write.csv(uncertainty, file.path(tables_dir, files$tables[[3L]]), row.names = FALSE, na = "")
  jss_optimizer_lean_write_figure(uncertainty, file.path(figures_dir, files$figures[[1L]]))
  jss_optimizer_lean_validate(data_dir, tables_dir, figures_dir)
  list(
    module_id = "03-joint-vs-separate-optimization", status = "current",
    data = file.path(data_dir, files$data), tables = file.path(tables_dir, files$tables),
    figures = file.path(figures_dir, files$figures),
    notes = "Three-scenario paired RS benchmark; supplementary; runtime descriptive and order-confounded."
  )
}

# Overrides used by the Phase 2 target graph. Reviewer-facing paper and full
# profiles both stage the versioned analysis-ready bundle. Running or resuming
# checkpoints remains a separate author workflow.
jss_phase2_optimizer_files <- jss_optimizer_lean_files

jss_phase2_optimizer_validate_candidate <- function(data_dir, tables_dir, figures_dir,
    root = getwd(), ..., require_promotion = FALSE) {
  jss_optimizer_lean_validate(data_dir, tables_dir, figures_dir)
}

jss_phase2_optimizer_validate <- function(data_dir, tables_dir, figures_dir,
    root = getwd(), ...) {
  jss_optimizer_lean_validate(data_dir, tables_dir, figures_dir)
}

jss_run_phase2_optimizer <- function(settings) {
  if (!settings$profile %in% c("paper", "full")) {
    return(list(module_id = "03-joint-vs-separate-optimization",
      status = "not_run_in_smoke", data = character(), tables = character(), figures = character()))
  }
  files <- jss_optimizer_lean_files()
  tracked <- file.path(settings$public_data_dir, "optimizer-benchmark")
  source <- c(file.path(tracked, "data", files$data), file.path(tracked, "tables", files$tables))
  destination <- c(file.path(settings$data_dir, files$data), file.path(settings$tables_dir, files$tables))
  jss_phase2_copy_exact(source, destination)
  uncertainty <- utils::read.csv(
    file.path(settings$tables_dir, files$tables[[3L]]),
    stringsAsFactors = FALSE
  )
  figure <- file.path(settings$figures_dir, files$figures[[1L]])
  jss_optimizer_lean_write_figure(uncertainty, figure)
  jss_optimizer_lean_validate(settings$data_dir, settings$tables_dir, settings$figures_dir)
  list(
    module_id = "03-joint-vs-separate-optimization", status = "current",
    data = file.path(settings$data_dir, files$data), tables = file.path(settings$tables_dir, files$tables),
    figures = figure,
    notes = "Three-scenario paired RS benchmark; supplementary; runtime descriptive and order-confounded."
  )
}
