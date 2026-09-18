#!/usr/bin/env Rscript

# Build a self-contained review bundle after the staged ACTG search completes.
# This is deliberately outside the manuscript targets graph.

script_argument <- grep("^--file=", commandArgs(trailingOnly = FALSE), value = TRUE)
script_path <- if (length(script_argument)) {
  normalizePath(sub("^--file=", "", script_argument[[1L]]), mustWork = TRUE)
} else {
  normalizePath(file.path(
    getwd(), "paper", "reviews", "2026-09-05-medical-candidates",
    "build-final-review.R"
  ), mustWork = TRUE)
}
package_root <- normalizePath(file.path(dirname(script_path), "..", "..", ".."))
default_results_dir <- file.path(dirname(script_path), "results")

argument_value <- function(name, default) {
  prefix <- paste0("--", name, "=")
  hit <- grep(paste0("^", prefix), commandArgs(trailingOnly = TRUE), value = TRUE)
  if (!length(hit)) return(default)
  sub(prefix, "", hit[[1L]], fixed = TRUE)
}

results_dir <- argument_value("source-results-dir", default_results_dir)

fit_path <- argument_value(
  "fit",
  file.path(results_dir, "actg-backward-selection-theta-final-fit.rds")
)
output_dir <- argument_value(
  "output-dir",
  file.path(results_dir, "actg-final-model-review")
)
pit_seed <- as.integer(argument_value("pit-seed", "20260909"))
skip_vcov <- tolower(argument_value("skip-vcov", "false")) %in% c("true", "1", "yes")

if (!file.exists(fit_path)) stop("Final fitted model not found: ", fit_path)
if (!is.finite(pit_seed)) stop("--pit-seed must be a finite integer.")
if (!requireNamespace("devtools", quietly = TRUE)) {
  stop("Install devtools with install.packages('devtools').")
}
if (!requireNamespace("ggplot2", quietly = TRUE)) {
  stop("Install ggplot2 with install.packages('ggplot2').")
}

devtools::load_all(package_root, quiet = TRUE)

`%||%` <- function(x, y) if (is.null(x)) y else x

tables_dir <- file.path(output_dir, "tables")
figures_dir <- file.path(output_dir, "figures")
objects_dir <- file.path(output_dir, "objects")
selection_dir <- file.path(output_dir, "model-selection")
logs_dir <- file.path(output_dir, "logs")
for (path in c(output_dir, tables_dir, figures_dir, objects_dir, selection_dir, logs_dir)) {
  dir.create(path, recursive = TRUE, showWarnings = FALSE)
}

artifact_status <- data.frame(
  artifact = character(), status = character(), started_at = character(),
  elapsed_sec = numeric(), warnings = character(), error = character(),
  stringsAsFactors = FALSE
)
status_path <- file.path(logs_dir, "artifact-status.csv")

run_artifact <- function(name, expr, required = FALSE) {
  started <- Sys.time()
  warnings_seen <- character()
  error_seen <- NA_character_
  value <- tryCatch(
    withCallingHandlers(
      force(expr),
      warning = function(warning) {
        warnings_seen <<- c(warnings_seen, conditionMessage(warning))
        invokeRestart("muffleWarning")
      }
    ),
    error = function(error) {
      error_seen <<- conditionMessage(error)
      NULL
    }
  )
  row <- data.frame(
    artifact = name,
    status = if (is.na(error_seen)) "complete" else "failed",
    started_at = format(started, "%Y-%m-%d %H:%M:%S %Z"),
    elapsed_sec = round(as.numeric(difftime(Sys.time(), started, units = "secs")), 2),
    warnings = paste(unique(warnings_seen), collapse = " | "),
    error = ifelse(is.na(error_seen), "", error_seen),
    stringsAsFactors = FALSE
  )
  artifact_status <<- rbind(artifact_status, row)
  utils::write.csv(artifact_status, status_path, row.names = FALSE, na = "")
  cat(row$status, ": ", name, " (", row$elapsed_sec, " sec)\n", sep = "")
  if (isTRUE(required) && !is.na(error_seen)) stop(name, ": ", error_seen)
  value
}

write_csv <- function(x, path) {
  utils::write.csv(as.data.frame(x), path, row.names = FALSE, na = "")
  invisible(path)
}

deparse_one <- function(x) paste(deparse(x, width.cutoff = 500L), collapse = " ")

save_plot <- function(plot, stem, width = 7, height = 5) {
  if (is.null(plot)) return(invisible(NULL))
  png_path <- file.path(figures_dir, paste0(stem, ".png"))
  pdf_path <- file.path(figures_dir, paste0(stem, ".pdf"))
  ggplot2::ggsave(png_path, plot = plot, width = width, height = height, units = "in", dpi = 320, bg = "white")
  ggplot2::ggsave(pdf_path, plot = plot, width = width, height = height, units = "in", bg = "white")
  invisible(c(png_path, pdf_path))
}

safe_name <- function(x) {
  x <- tolower(gsub("[^A-Za-z0-9]+", "-", x))
  gsub("(^-+|-+$)", "", x)
}

save_nested_plots <- function(x, prefix, width = 8, height = 5) {
  if (inherits(x, c("ggplot", "grob", "gtable", "patchwork"))) {
    save_plot(x, safe_name(prefix), width = width, height = height)
    return(invisible(NULL))
  }
  if (!is.list(x) || is.data.frame(x)) return(invisible(NULL))
  child_names <- names(x)
  if (is.null(child_names)) child_names <- paste0("plot-", seq_along(x))
  for (index in seq_along(x)) {
    save_nested_plots(x[[index]], paste(prefix, child_names[[index]], sep = "-"), width, height)
  }
  invisible(NULL)
}

html_escape <- function(x) {
  x <- gsub("&", "&amp;", as.character(x), fixed = TRUE)
  x <- gsub("<", "&lt;", x, fixed = TRUE)
  x <- gsub(">", "&gt;", x, fixed = TRUE)
  x <- gsub('"', "&quot;", x, fixed = TRUE)
  x
}

html_table <- function(x) {
  x <- as.data.frame(x, stringsAsFactors = FALSE)
  header <- paste0("<tr>", paste0("<th>", html_escape(names(x)), "</th>", collapse = ""), "</tr>")
  body <- apply(x, 1L, function(row) {
    paste0("<tr>", paste0("<td>", html_escape(row), "</td>", collapse = ""), "</tr>")
  })
  paste0("<table><thead>", header, "</thead><tbody>", paste(body, collapse = "\n"), "</tbody></table>")
}

fit <- run_artifact("load final selected fit", readRDS(fit_path), required = TRUE)
if (!inherits(fit, "gamlss.longitudinal")) stop("The final RDS is not a gamlss.longitudinal fit.")

run_artifact("copy final model and staged audit files", {
  file.copy(fit_path, file.path(objects_dir, "final-model.rds"), overwrite = TRUE)
  patterns <- c(
    "actg-(zibnb-)?backward-selection-(mu-smooth|sigma|nu|tau|theta)-(model-selection-table|path|final-formulas|run-manifest)\\.(csv|rds)$",
    "actg-(zibnb-)?(covariate-|)refinement-cycle-[0-9]+-(mu|sigma|nu|tau|theta)-(model-selection-table|path|final-formulas|run-manifest)\\.(csv|rds)$",
    "actg-covariate-refinement-(cycle-summary|final-formulas|run-manifest)\\.(csv|rds)$",
    "(actg-backward-selection-all-stages-complete|ZIBNB-SELECTION-COMPLETE)\\.rds$"
  )
  candidates <- list.files(results_dir, full.names = TRUE)
  keep <- Reduce(`|`, lapply(patterns, grepl, x = basename(candidates), perl = TRUE))
  if (any(keep)) file.copy(candidates[keep], selection_dir, overwrite = TRUE)
})

point_summary <- run_artifact(
  "point-estimate model summary",
  summary(fit, include_vcov = FALSE),
  required = TRUE
)
saveRDS(point_summary, file.path(objects_dir, "model-summary-point-estimates.rds"))

margin_label <- as.character(point_summary$model$margin_dist)
copula_code <- as.character(point_summary$model$copula_dist)
copula_names <- c(
  N = "Gaussian", C = "Clayton", F = "Frank", G = "Gumbel",
  J = "Joe", t = "Student-t"
)
copula_label <- unname(copula_names[copula_code])
if (!length(copula_label) || is.na(copula_label)) copula_label <- copula_code
model_label <- paste(margin_label, copula_label, sep = "--")

active_formula_parameters <- names(fit$model_matrix$x)
formula_parameters <- intersect(
  c("mu", "sigma", "nu", "tau", "theta", "zeta"),
  active_formula_parameters
)
formula_table <- data.frame(
  parameter = formula_parameters,
  formula = vapply(fit$formulas[formula_parameters], deparse_one, character(1)),
  stringsAsFactors = FALSE
)
write_csv(formula_table, file.path(tables_dir, "final-formulas.csv"))

model_table <- data.frame(
  statistic = c(
    "Margin distribution", "Copula distribution", "Observed responses",
    "Submitted rows", "Subjects", "Time points", "Fixed coefficients",
    "Smooth terms", "Smooth EDF", "Log likelihood", "AIC", "BIC",
    "Converged", "Outer iterations", "Stop reason", "Missingness likelihood"
  ),
  value = c(
    point_summary$model$margin_dist,
    point_summary$model$copula_dist,
    stats::nobs(fit, type = "observed"),
    stats::nobs(fit, type = "submitted"),
    point_summary$model$n_subjects,
    point_summary$model$n_timepoints,
    point_summary$model$n_fixed,
    point_summary$model$n_smooth_terms,
    point_summary$model$edf_smooth,
    point_summary$fit$logLik,
    point_summary$fit$AIC,
    point_summary$fit$BIC,
    fit$convergence$converged,
    fit$convergence$outer_iterations,
    fit$convergence$stop_reason,
    fit$likelihood_contract$objective %||% fit$missingness
  ),
  stringsAsFactors = FALSE
)
write_csv(model_table, file.path(tables_dir, "model-fit-statistics.csv"))
write_csv(point_summary$smooth_terms, file.path(tables_dir, "smooth-terms-edf.csv"))
write_csv(point_summary$coefficients, file.path(tables_dir, "coefficient-point-estimates.csv"))

run_artifact("combine staged model-selection tables", {
  model_paths <- list.files(
    results_dir,
    pattern = paste0(
      "^actg-(zibnb-)?(backward-selection-(mu-smooth|sigma|nu|tau|theta)|",
      "(covariate-|)refinement-cycle-[0-9]+-(mu|sigma|nu|tau|theta))-",
      "model-selection-table\\.csv$"
    ),
    full.names = TRUE
  )
  path_paths <- sub("model-selection-table\\.csv$", "path.csv", model_paths)
  all_candidates <- list()
  accepted_paths <- list()
  for (index in seq_along(model_paths)) {
    model_path <- model_paths[[index]]
    path_path <- path_paths[[index]]
    stage <- sub("^actg-(zibnb-)?", "", basename(model_path))
    stage <- sub("-(model-selection-table|path)\\.csv$", "", stage)
    if (file.exists(model_path)) {
      value <- utils::read.csv(model_path, stringsAsFactors = FALSE, check.names = FALSE)
      value$stage <- stage
      all_candidates[[stage]] <- value
    }
    if (file.exists(path_path)) {
      value <- utils::read.csv(path_path, stringsAsFactors = FALSE, check.names = FALSE)
      value$stage <- stage
      accepted_paths[[stage]] <- value
    }
  }
  if (length(all_candidates)) {
    write_csv(do.call(rbind, all_candidates), file.path(selection_dir, "all-candidate-models.csv"))
  }
  if (length(accepted_paths)) {
    write_csv(do.call(rbind, accepted_paths), file.path(selection_dir, "accepted-selection-paths.csv"))
  }
})

inference_summary <- NULL
fit_with_vcov <- fit
if (!skip_vcov) {
  inference_summary <- run_artifact(
    "subject-clustered sandwich coefficient inference",
    summary(fit, include_vcov = TRUE, vcov_method = "sandwich")
  )
  if (!is.null(inference_summary)) {
    saveRDS(inference_summary, file.path(objects_dir, "model-summary-sandwich.rds"))
    fit_with_vcov$vcov <- inference_summary$vcov
    fit_with_vcov$vcov_meta$precomputed <- TRUE
    fit_with_vcov$vcov_meta$numderiv <- FALSE
    fit_with_vcov$vcov_meta$cache_version <- 1L
    fit_with_vcov$vcov_meta$method <- "sandwich"
    fit_with_vcov$vcov_meta$method_used <- inference_summary$vcov$method %||% "sandwich_cluster"
    saveRDS(fit_with_vcov, file.path(objects_dir, "final-model-with-sandwich-vcov.rds"))
  }
}

model_table <- rbind(
  model_table,
  data.frame(
    statistic = "Primary coefficient VCOV",
    value = if (is.null(inference_summary)) "unavailable" else "subject-clustered sandwich",
    stringsAsFactors = FALSE
  )
)
write_csv(model_table, file.path(tables_dir, "model-fit-statistics.csv"))

coef_source <- if (is.null(inference_summary)) point_summary$coefficients else inference_summary$coefficients
coef_numeric <- as.data.frame(coef_source, stringsAsFactors = FALSE)
coef_numeric$coefficient <- coef_numeric$term
coef_numeric$term <- sub("^[^.]+\\.", "", coef_numeric$term)
coef_numeric$term[coef_numeric$term == "intercept"] <- "(Intercept)"
coef_numeric$component <- ifelse(coef_numeric$parameter %in% c("theta", "zeta"), "Copula", "Margin")
coef_numeric$statistic <- coef_numeric$estimate / coef_numeric$std_error
coef_numeric$conf_low <- coef_numeric$estimate - stats::qnorm(0.975) * coef_numeric$std_error
coef_numeric$conf_high <- coef_numeric$estimate + stats::qnorm(0.975) * coef_numeric$std_error
coef_numeric <- coef_numeric[, c(
  "component", "parameter", "term", "coefficient", "estimate", "std_error",
  "model_std_error", "robust_std_error", "statistic", "p_value", "conf_low",
  "conf_high", "signif"
)]
write_csv(coef_numeric, file.path(tables_dir, "coefficient-table-numeric.csv"))

format_number <- function(x, digits = 3L) ifelse(is.finite(x), formatC(x, format = "f", digits = digits), "")
format_p <- function(x) ifelse(
  !is.finite(x), "",
  ifelse(x < 0.001, "<0.001", formatC(x, format = "f", digits = 3L))
)
coef_publication <- data.frame(
  Component = coef_numeric$component,
  Parameter = coef_numeric$parameter,
  Term = coef_numeric$term,
  Estimate = format_number(coef_numeric$estimate),
  `Robust SE` = format_number(coef_numeric$robust_std_error),
  `Model-based SE` = format_number(coef_numeric$model_std_error),
  `95% CI` = ifelse(
    is.finite(coef_numeric$conf_low) & is.finite(coef_numeric$conf_high),
    paste0("[", format_number(coef_numeric$conf_low), ", ", format_number(coef_numeric$conf_high), "]"),
    ""
  ),
  `p-value` = format_p(coef_numeric$p_value),
  check.names = FALSE,
  stringsAsFactors = FALSE
)
write_csv(coef_publication, file.path(tables_dir, "coefficient-table-publication.csv"))

run_artifact("publication coefficient LaTeX table", {
  if (!requireNamespace("knitr", quietly = TRUE)) stop("knitr is not installed")
  latex <- knitr::kable(
    coef_publication, format = "latex", booktabs = TRUE,
    caption = paste0("Final ACTG longitudinal ", model_label, " model coefficients. Subject-clustered sandwich standard errors are primary.")
  )
  writeLines(as.character(latex), file.path(tables_dir, "coefficient-table-publication.tex"))
})

if (!is.null(inference_summary)) {
  run_artifact("write inference contract", {
    contract <- inference_summary$fit$inference_contract %||% attr(inference_summary, "inference_contract")
    writeLines(capture.output(str(contract, max.level = 5)), file.path(tables_dir, "inference-contract.txt"))
    diagnostics <- inference_summary$fit$hessian_diagnostics
    writeLines(capture.output(str(diagnostics, max.level = 5)), file.path(tables_dir, "sandwich-bread-diagnostics.txt"))
  })
}

model_check <- run_artifact(
  "marginal and residual model checks",
  check_model(
    fit_with_vcov,
    include_vcov = !is.null(inference_summary),
    include_plots = TRUE,
    vcov_method = if (is.null(inference_summary)) NULL else "sandwich",
    pit_seed = pit_seed,
    residual_lags = 1:3
  )
)
if (!is.null(model_check)) {
  saveRDS(model_check, file.path(objects_dir, "check-model.rds"))
  for (name in c("scores", "pit", "tail", "residual_dependence", "basic_checks", "checks", "flags")) {
    value <- model_check[[name]]
    if (is.data.frame(value)) write_csv(value, file.path(tables_dir, paste0("diagnostic-", name, ".csv")))
  }
  run_artifact("standard diagnostic figures", {
    for (name in names(model_check$plots)) {
      save_plot(model_check$plots[[name]], paste0("diagnostic-", name))
    }
  })
}

run_artifact("time-stratified diagnostic figures", {
  save_plot(pithist(fit, randomize = TRUE, seed = pit_seed, by_time = TRUE), "diagnostic-pit-by-time", 9, 7)
  save_plot(qqrplot(fit, randomize = TRUE, seed = pit_seed, by_time = TRUE), "diagnostic-qq-by-time", 9, 7)
  save_plot(wormplot(fit, randomize = TRUE, seed = pit_seed, by_time = TRUE), "diagnostic-worm-by-time", 9, 7)
  save_plot(rootogram(fit, by_time = TRUE), "diagnostic-rootogram-by-time", 9, 7)
})

copula_diagnostics <- run_artifact(
  "copula diagnostics",
  plot_copula_diagnostics(
    fit, plot = FALSE, type = "all", lags = 1,
    residual_lags = 1:3, tail_thresholds = c(0.05, 0.10, 0.20)
  )
)
if (!is.null(copula_diagnostics)) {
  saveRDS(copula_diagnostics, file.path(objects_dir, "copula-diagnostics.rds"))
  for (name in names(copula_diagnostics)) {
    value <- copula_diagnostics[[name]]
    if (is.data.frame(value)) {
      write_csv(value, file.path(tables_dir, paste0("copula-", safe_name(name), ".csv")))
    }
  }
  run_artifact("copula diagnostic figures", {
    save_plot(copula_diagnostics$dashboard, "copula-diagnostics-dashboard", 12, 16)
    for (name in names(copula_diagnostics$plots)) {
      save_plot(copula_diagnostics$plots[[name]], paste0("copula-", safe_name(name)), 7, 5)
    }
  })
}

run_artifact("row-level residual and calibration summaries", {
  dataset <- as.data.frame(fit$dataset, stringsAsFactors = FALSE)
  dataset$fitted_mu <- stats::fitted(fit, parameter = "mu", finite = FALSE)
  dataset$response_residual <- dataset$response - dataset$fitted_mu
  pearson_all <- stats::residuals(fit, type = "pearson", finite = FALSE)
  dataset$pearson_residual <- pearson_all
  observed <- dataset[is.finite(dataset$response) & is.finite(dataset$fitted_mu), , drop = FALSE]
  pit <- pithist(fit, randomize = TRUE, seed = pit_seed, plot = FALSE)
  if (nrow(pit) == nrow(observed)) {
    observed$randomized_pit <- pit$pit
    observed$randomized_quantile_residual <- stats::qnorm(pmin(pmax(pit$pit, .Machine$double.eps), 1 - .Machine$double.eps))
  }
  write_csv(observed, file.path(tables_dir, "observation-level-diagnostics.csv"))

  group <- interaction(observed$time, observed$treatment, drop = TRUE, lex.order = TRUE)
  pieces <- split(observed, group)
  by_time_treatment <- do.call(rbind, lapply(pieces, function(value) {
    data.frame(
      time = value$time[[1L]], treatment = as.character(value$treatment[[1L]]),
      n = nrow(value), observed_mean = mean(value$response),
      fitted_mean = mean(value$fitted_mu), mean_residual = mean(value$response_residual),
      rmse = sqrt(mean(value$response_residual^2)),
      pit_mean = if ("randomized_pit" %in% names(value)) mean(value$randomized_pit) else NA_real_,
      pit_sd = if ("randomized_pit" %in% names(value)) stats::sd(value$randomized_pit) else NA_real_,
      stringsAsFactors = FALSE
    )
  }))
  rownames(by_time_treatment) <- NULL
  write_csv(by_time_treatment, file.path(tables_dir, "calibration-by-time-and-treatment.csv"))

  observed$time_factor <- factor(observed$time)
  p_fit <- ggplot2::ggplot(by_time_treatment, ggplot2::aes(x = time, group = treatment, colour = treatment)) +
    ggplot2::geom_line(ggplot2::aes(y = observed_mean, linetype = "Observed"), linewidth = 0.8) +
    ggplot2::geom_point(ggplot2::aes(y = observed_mean), size = 2) +
    ggplot2::geom_line(ggplot2::aes(y = fitted_mean, linetype = "Fitted"), linewidth = 0.8) +
    ggplot2::geom_point(ggplot2::aes(y = fitted_mean), shape = 1, size = 2) +
    ggplot2::labs(x = "Week", y = "Mean CD4 count", colour = "Treatment", linetype = NULL, title = "Observed and fitted means") +
    theme_gamlss_longitudinal()
  save_plot(p_fit, "observed-vs-fitted-by-time-treatment", 9, 6)

  p_resid <- ggplot2::ggplot(observed, ggplot2::aes(x = fitted_mu, y = response_residual)) +
    ggplot2::geom_point(alpha = 0.25, size = 0.8) +
    ggplot2::geom_hline(yintercept = 0, linetype = "dashed") +
    ggplot2::geom_smooth(method = "loess", se = TRUE, linewidth = 0.7) +
    ggplot2::labs(x = "Fitted mu", y = "Response residual", title = "Residuals versus fitted values") +
    theme_gamlss_longitudinal()
  save_plot(p_resid, "residuals-vs-fitted", 7, 5)

  if ("randomized_quantile_residual" %in% names(observed)) {
    p_time <- ggplot2::ggplot(observed, ggplot2::aes(x = time_factor, y = randomized_quantile_residual)) +
      ggplot2::geom_boxplot(outlier.alpha = 0.15) +
      ggplot2::geom_hline(yintercept = 0, linetype = "dashed") +
      ggplot2::labs(x = "Week", y = "Randomized quantile residual", title = "Quantile residuals by time") +
      theme_gamlss_longitudinal()
    save_plot(p_time, "quantile-residuals-by-time", 7, 5)
  }
})

if (!is.null(inference_summary)) {
  term_plots <- run_artifact(
    "final model term-effect figures",
    plot_terms(
      fit_with_vcov, data = fit$dataset_original,
      vcov_method = "sandwich", ncol = 2, include_intercept = FALSE
    )
  )
  if (!is.null(term_plots)) {
    saveRDS(term_plots, file.path(objects_dir, "term-plots.rds"))
    run_artifact("save term-effect figures", {
      save_nested_plots(term_plots, "term-effects", width = 9, height = 7)
    })
  }
}

run_artifact("session information", {
  writeLines(capture.output(utils::sessionInfo()), file.path(logs_dir, "session-info.txt"))
})

generated_at <- format(Sys.time(), "%Y-%m-%d %H:%M:%S %Z")
formula_lines <- paste0("- `", formula_table$parameter, "`: `", formula_table$formula, "`")
inference_note <- if (is.null(inference_summary)) {
  "Cluster-robust covariance could not be produced; inspect `logs/artifact-status.csv`. The coefficient table therefore contains point estimates only."
} else {
  "The primary coefficient uncertainty uses subject-clustered sandwich standard errors. Model-based standard errors are included as a sensitivity comparison."
}
readme <- c(
  paste0("# Final ACTG ", margin_label, "–", copula_label, " model review bundle"), "",
  paste0("Generated: ", generated_at), "",
  "## Final formulas", "", formula_lines, "",
  "## Start here", "",
  "- `review-index.html`: compact browsable summary.",
  "- `tables/coefficient-table-publication.csv`: formatted publication coefficient table.",
  "- `tables/coefficient-table-numeric.csv`: unrounded numeric coefficient results.",
  "- `tables/model-fit-statistics.csv`: likelihood, information criteria, sample sizes, and convergence.",
  "- `tables/smooth-terms-edf.csv`: retained smooth terms and effective degrees of freedom.",
  "- `model-selection/all-candidate-models.csv`: every candidate evaluated across the staged and refinement searches.",
  "- `model-selection/accepted-selection-paths.csv`: the accepted greedy AIC path.",
  "- `figures/`: marginal, residual, copula, calibration, and term-effect diagnostics in PNG and PDF.",
  "- `logs/artifact-status.csv`: success, warnings, runtime, and any failure for every output step.", "",
  "## Interpretation notes", "",
  paste0("- ", inference_note),
  "- Confidence intervals and p-values are conditional on the final selected model and do not account for the preceding model-selection search.",
  "- AIC/BIC use the stated segmented likelihood: observations in separate contiguous segments for a subject are treated as independent.",
  paste0("- Discrete-margin PIT diagnostics use seeded randomized PIT values (seed ", pit_seed, ")."),
  "- Penalized smooths are summarized by EDF and effect plots; basis coefficients are not suitable as a publication coefficient table.",
  "- Diagnostic thresholds are descriptive aids, not automatic goodness-of-fit verdicts."
)
writeLines(readme, file.path(output_dir, "README.md"))

figure_pngs <- sort(list.files(figures_dir, pattern = "\\.png$"))
figure_html <- paste0(
  "<section><h3>", html_escape(sub("\\.png$", "", figure_pngs)), "</h3>",
  "<a href='figures/", html_escape(figure_pngs), "'><img src='figures/", html_escape(figure_pngs), "'></a></section>",
  collapse = "\n"
)
html <- c(
  "<!doctype html><html><head><meta charset='utf-8'><title>ACTG final model review</title>",
  "<style>body{font-family:Arial,sans-serif;max-width:1200px;margin:2rem auto;padding:0 1rem;color:#222}table{border-collapse:collapse;width:100%;margin:1rem 0 2rem}th,td{border:1px solid #ccc;padding:.45rem;text-align:left;font-size:.9rem}th{background:#f3f5f7}code{background:#f3f5f7;padding:.15rem .3rem}img{max-width:100%;height:auto;border:1px solid #ddd}section{margin:2rem 0}.note{background:#fff7df;border-left:4px solid #d99b00;padding:1rem}</style></head><body>",
  paste0("<h1>Final ACTG ", html_escape(margin_label), "–", html_escape(copula_label), " model review</h1>"),
  paste0("<p>Generated ", html_escape(generated_at), ".</p>"),
  paste0("<div class='note'><p>", html_escape(inference_note), "</p><p>Inference is conditional on the AIC-selected model and does not account for model selection. The segmented-likelihood assumption applies.</p></div>"),
  "<h2>Final formulas</h2>", html_table(formula_table),
  "<h2>Model fit</h2>", html_table(model_table),
  "<h2>Publication coefficients</h2>", html_table(coef_publication),
  "<h2>Diagnostic figures</h2>", figure_html,
  "</body></html>"
)
writeLines(html, file.path(output_dir, "review-index.html"))

completion <- list(
  completed_at = Sys.time(), fit_path = normalizePath(fit_path),
  output_dir = normalizePath(output_dir), pit_seed = pit_seed,
  inference_available = !is.null(inference_summary), artifact_status = artifact_status
)
saveRDS(completion, file.path(output_dir, "REVIEW-BUNDLE-COMPLETE.rds"))
cat("\nFinal review bundle written to: ", normalizePath(output_dir), "\n", sep = "")
