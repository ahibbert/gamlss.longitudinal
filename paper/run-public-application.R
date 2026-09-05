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

write_distribution_selection <- function(selection, directory) {
  saveRDS(selection, file.path(directory, "distribution-selection.rds"))
  write_capture(selection, file.path(directory, "distribution-selection.txt"))
  utils::write.csv(
    as.data.frame(selection$margin),
    file.path(directory, "margin-selection.csv"),
    row.names = FALSE
  )
  utils::write.csv(
    as.data.frame(selection$copula),
    file.path(directory, "copula-selection.csv"),
    row.names = FALSE
  )
  invisible(selection)
}

load_or_run_distribution_selection <- function(name, directory, force = FALSE) {
  path <- file.path(directory, "distribution-selection.rds")
  if (!force && file.exists(path)) {
    selection <- readRDS(path)
    registry_version <- unique(
      gamlss.longitudinal::longitudinal_capabilities("routes")$registry_version
    )
    if (identical(selection$contract$registry_version, registry_version)) {
      message("Reusing ", path)
      return(selection)
    }
    message("Refreshing stale distribution selection after registry change.")
  }
  selection <- jss_select_public_application_distributions(name)
  write_distribution_selection(selection, directory)
}

for (name in candidates) {
  message("Public application candidate: ", name)
  directory <- file.path(output_root, name)
  dir.create(directory, recursive = TRUE, showWarnings = FALSE)

  contract <- jss_public_application_contract(name)
  saveRDS(contract, file.path(directory, "data-contract.rds"))
  write_capture(contract, file.path(directory, "data-contract.txt"))

  selection <- NULL
  if (stage %in% c("screen", "all")) {
    selection <- load_or_run_distribution_selection(name, directory, force = force)
  }

  if (stage %in% c("fit", "all")) {
    if (is.null(selection)) {
      selection <- load_or_run_distribution_selection(name, directory, force = FALSE)
    }
    fit_path <- file.path(directory, "full-fit.rds")
    fit_is_selected <- FALSE
    if (!force && file.exists(fit_path)) {
      fit <- readRDS(fit_path)
      fit_is_selected <- identical(
        fit$margin_dist$family[[1L]], selection$selected_margin
      ) && identical(fit$copula_dist, selection$selected_copula)
    }
    if (force || !file.exists(fit_path) || !fit_is_selected) {
      if (file.exists(fit_path) && !force) {
        message("Replacing stale fit whose distributions do not match the saved selection.")
      }
      fit <- jss_fit_public_application(
        name,
        margin = selection$selected_margin,
        copula = selection$selected_copula,
        compute_vcov = compute_vcov
      )
      saveRDS(fit, fit_path)
    } else {
      message("Reusing ", fit_path)
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
      write_capture(
        summary(fit, vcov_method = fit$vcov_meta$method),
        file.path(directory, "full-fit-summary.txt")
      )
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
    if (is.null(selection)) {
      selection <- load_or_run_distribution_selection(name, directory, force = FALSE)
    }
    base_fit <- readRDS(base_path)
    if (!identical(base_fit$margin_dist$family[[1L]], selection$selected_margin) ||
        !identical(base_fit$copula_dist, selection$selected_copula)) {
      stop(
        "The patents base fit predates the saved distribution selection; rerun the fit stage.",
        call. = FALSE
      )
    }
    refinement_name <- paste0(selection$selected_margin, "_rd_dispersion")
    refined_spec <- jss_patents_joint_refinement_specs()[[refinement_name]]
    if (is.null(refined_spec)) {
      stop("No patents refinement is defined for ", refinement_name, ".", call. = FALSE)
    }
    refined_path <- file.path(
      refinement_directory,
      paste0(refinement_name, "-", selection$selected_copula, "-fit.rds")
    )
    if (force || !file.exists(refined_path)) {
      refined_fit <- jss_fit_public_application(
        "patents",
        formulas = refined_spec$formulas,
        margin = refined_spec$margin,
        copula = selection$selected_copula,
        compute_vcov = compute_vcov,
        outer_tol = 0.01,
        inner_tol = 1e-5,
        max_outer_iter = 60L,
        max_inner_iter = 150L,
        max_elapsed_sec = 1200,
        start_from = jss_patents_rd_dispersion_start(base_fit)
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
