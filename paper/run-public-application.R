#!/usr/bin/env Rscript

# Candidate application runner. These models are intentionally separate from
# the default paper targets graph until a dataset is selected for the vignette
# and manuscript.

arguments <- commandArgs(trailingOnly = TRUE)

option_value <- function(flag, default = NULL) {
  exact <- which(arguments == flag)
  if (length(exact) && exact[[1L]] < length(arguments)) return(arguments[[exact[[1L]] + 1L]])
  prefix <- paste0(flag, "=")
  inline <- arguments[startsWith(arguments, prefix)]
  if (length(inline)) return(sub(prefix, "", inline[[1L]], fixed = TRUE))
  default
}

has_flag <- function(flag) flag %in% arguments

candidate <- option_value("--dataset", "patents")
stage <- match.arg(option_value("--stage", "fit"), c("screen", "fit", "refine", "all"))
output_root <- option_value("--output-dir", file.path("results", "public-applications"))
compute_vcov <- !has_flag("--no-vcov")
run_comparators <- has_flag("--comparators")
generate_outputs <- has_flag("--outputs")
force <- has_flag("--force")

if (requireNamespace("pkgload", quietly = TRUE) && file.exists("DESCRIPTION")) {
  suppressPackageStartupMessages(pkgload::load_all(".", quiet = TRUE))
} else {
  suppressPackageStartupMessages(library(gamlss.longitudinal))
}
source(file.path("paper", "R", "10-public-application-candidates.R"))

if (identical(candidate, "all")) {
  candidates <- jss_public_application_names()
} else {
  candidates <- match.arg(candidate, jss_public_application_names())
}

write_capture <- function(object, path) {
  writeLines(utils::capture.output(print(object)), path, useBytes = TRUE)
}

write_workflow_outputs <- function(fit, directory) {
  outputs <- tryCatch(
    jss_public_application_outputs(fit),
    error = function(e) structure(
      list(error = conditionMessage(e)),
      class = "candidate_output_error"
    )
  )
  saveRDS(outputs, file.path(directory, "workflow-outputs.rds"))
  write_capture(outputs, file.path(directory, "workflow-outputs.txt"))
  if (inherits(outputs, "candidate_output_error")) return(invisible(outputs))

  if (is.data.frame(outputs$coefficient_table)) {
    utils::write.csv(
      outputs$coefficient_table,
      file.path(directory, "coefficient-table.csv"),
      row.names = FALSE
    )
  }
  if (inherits(outputs$margin_diagnostics, c("gg", "ggplot"))) {
    ggplot2::ggsave(
      file.path(directory, "margin-diagnostics-4-timepoints.png"),
      outputs$margin_diagnostics,
      width = 10,
      height = 8,
      dpi = 180
    )
  }
  if (inherits(outputs$copula_diagnostics$dashboard, c("gg", "ggplot"))) {
    ggplot2::ggsave(
      file.path(directory, "copula-diagnostics-dashboard.png"),
      outputs$copula_diagnostics$dashboard,
      width = 12,
      height = 9,
      dpi = 180
    )
  }
  invisible(outputs)
}

for (name in candidates) {
  message("Public application candidate: ", name)
  directory <- file.path(output_root, name)
  dir.create(directory, recursive = TRUE, showWarnings = FALSE)

  contract <- jss_public_application_contract(name)
  saveRDS(contract, file.path(directory, "data-contract.rds"))
  write_capture(contract, file.path(directory, "data-contract.txt"))

  if (stage %in% c("screen", "all")) {
    screen_path <- file.path(directory, "margin-screen.csv")
    if (force || !file.exists(screen_path)) {
      screen <- jss_screen_public_application_margins(name)
      utils::write.csv(screen, screen_path, row.names = FALSE)
    } else {
      message("Reusing ", screen_path)
    }
  }

  if (stage %in% c("fit", "all")) {
    fit_path <- file.path(directory, "full-fit.rds")
    if (force || !file.exists(fit_path)) {
      fit <- jss_fit_public_application(name, compute_vcov = compute_vcov)
      saveRDS(fit, fit_path)
    } else {
      message("Reusing ", fit_path)
      fit <- readRDS(fit_path)
    }

    coefficient_estimates <- data.frame(
      term = names(stats::coef(fit)),
      estimate = unname(stats::coef(fit)),
      stringsAsFactors = FALSE
    )
    utils::write.csv(
      coefficient_estimates,
      file.path(directory, "coefficient-estimates.csv"),
      row.names = FALSE
    )
    if (isTRUE(fit$vcov_meta$precomputed)) {
      write_capture(summary(fit), file.path(directory, "full-fit-summary.txt"))
    }
    write_capture(fit$convergence, file.path(directory, "convergence.txt"))

    if (isTRUE(fit$convergence$converged)) {
      metrics <- data.frame(
        candidate = name,
        observations = nobs(fit),
        df = attr(logLik(fit), "df"),
        logLik = as.numeric(logLik(fit)),
        AIC = AIC(fit),
        BIC = BIC(fit),
        stringsAsFactors = FALSE
      )
      utils::write.csv(metrics, file.path(directory, "fit-metrics.csv"), row.names = FALSE)

      if (generate_outputs) {
        write_workflow_outputs(fit, directory)
      }
    }

    if (run_comparators) {
      spec <- jss_public_application_spec(name)
      for (comparison in names(spec$comparisons)) {
        comparison_path <- file.path(directory, paste0(comparison, "-fit.rds"))
        if (force || !file.exists(comparison_path)) {
          comparison_fit <- jss_fit_public_application_comparator(
            fit,
            component = comparison,
            compute_vcov = FALSE
          )
          saveRDS(comparison_fit, comparison_path)
        }
      }
    }
  }

  if (stage %in% c("refine", "all") && identical(name, "patents")) {
    refinement_directory <- file.path(directory, "refinement")
    dir.create(refinement_directory, recursive = TRUE, showWarnings = FALSE)

    screen_path <- file.path(refinement_directory, "refined-margin-screen.csv")
    if (force || !file.exists(screen_path)) {
      refined_screen <- jss_patents_refined_margin_screen()
      utils::write.csv(refined_screen, screen_path, row.names = FALSE)
    }

    base_path <- file.path(directory, "full-fit.rds")
    if (!file.exists(base_path)) {
      stop("Run the patents fit stage before the refinement stage.", call. = FALSE)
    }
    refined_path <- file.path(refinement_directory, "NBI_rd_dispersion-fit.rds")
    if (force || !file.exists(refined_path)) {
      base_fit <- readRDS(base_path)
      refined_spec <- jss_patents_joint_refinement_specs()$NBI_rd_dispersion
      refined_fit <- jss_fit_public_application(
        "patents",
        formulas = refined_spec$formulas,
        margin = refined_spec$margin,
        compute_vcov = compute_vcov,
        outer_tol = 0.01,
        inner_tol = 1e-5,
        max_outer_iter = 60L,
        max_inner_iter = 150L,
        max_elapsed_sec = 1200,
        start_from = jss_patents_nbi_rd_start(base_fit)
      )
      saveRDS(refined_fit, refined_path)
    } else {
      message("Reusing ", refined_path)
      refined_fit <- readRDS(refined_path)
    }

    if (generate_outputs && isTRUE(refined_fit$convergence$converged)) {
      write_workflow_outputs(refined_fit, refinement_directory)
    }
  }
}

message("Candidate application run complete: ", normalizePath(output_root, winslash = "/", mustWork = FALSE))
