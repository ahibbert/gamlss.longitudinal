#!/usr/bin/env Rscript

# Recompute smooth EDF and information criteria for the original ACTG screening
# checkpoints without altering the archived fitted objects.

script_argument <- grep("^--file=", commandArgs(trailingOnly = FALSE), value = TRUE)
script_path <- if (length(script_argument)) {
  normalizePath(sub("^--file=", "", script_argument[[1L]]), mustWork = TRUE)
} else {
  normalizePath(file.path(
    getwd(), "paper", "reviews", "2026-09-05-medical-candidates",
    "audit-selection-edf.R"
  ), mustWork = TRUE)
}
package_root <- normalizePath(file.path(dirname(script_path), "..", "..", ".."))
results_dir <- file.path(dirname(script_path), "results")
output_path <- file.path(results_dir, "actg-original-selection-edf-audit.csv")

if (!requireNamespace("devtools", quietly = TRUE)) stop("devtools is required.")
devtools::load_all(package_root, quiet = TRUE)

stages <- c("mu-smooth", "sigma", "nu", "theta")
audit_rows <- list()

for (stage in stages) {
  table_path <- file.path(
    results_dir,
    paste0("actg-backward-selection-", stage, "-model-selection-table.csv")
  )
  table <- utils::read.csv(table_path, stringsAsFactors = FALSE, check.names = FALSE)
  corrected <- lapply(seq_len(nrow(table)), function(index) {
    checkpoint <- readRDS(table$checkpoint[[index]])
    fit <- checkpoint$fit
    refreshed_df <- .gl_refresh_rs_smooth_edf(
      df_s = fit$df_s,
      mm_s = fit$model_matrix$s,
      lambda_s = fit$lambda_s,
      weights_final = fit$weights
    )
    criteria <- .gl_fit_information_criteria(
      par_cov = fit$par,
      par_s = fit$par_s,
      df_s = refreshed_df,
      calc_lik_out_end = fit$calc_lik_out_end,
      dataset = fit$dataset
    )
    c(
      corrected_AIC = criteria["AIC", "joint"],
      corrected_BIC = criteria["BIC", "joint"],
      corrected_EDF = criteria["EDF", "joint"]
    )
  })
  corrected <- as.data.frame(do.call(rbind, corrected))
  table <- cbind(table, corrected)
  table$stage <- stage
  table$original_accepted <- table$accepted
  table$corrected_rank_within_step <- NA_integer_
  table$corrected_selected <- FALSE
  table$corrected_improves_parent <- FALSE

  for (step in sort(unique(table$step))) {
    step_index <- which(table$step == step)
    values <- table$corrected_AIC[step_index]
    table$corrected_rank_within_step[step_index] <- rank(values, ties.method = "min")
    best_index <- step_index[[which.min(values)]]
    if (step == 0L) {
      table$corrected_selected[[best_index]] <- TRUE
      table$corrected_improves_parent[[best_index]] <- TRUE
      next
    }
    prior_accepted <- which(table$step < step & table$original_accepted)
    if (!length(prior_accepted)) stop("No parent model found for ", stage, " step ", step, ".")
    parent_index <- prior_accepted[[which.max(table$step[prior_accepted])]]
    improves <- table$corrected_AIC[[best_index]] < table$corrected_AIC[[parent_index]]
    table$corrected_improves_parent[[best_index]] <- improves
    table$corrected_selected[[best_index]] <- improves
  }

  keep <- c(
    "stage", "step", "model_id", "parameter", "change_type", "dropped_term",
    "candidate_formula", "AIC", "corrected_AIC", "BIC", "corrected_BIC",
    "EDF", "corrected_EDF", "original_accepted", "corrected_selected",
    "criterion_rank_within_step", "corrected_rank_within_step", "checkpoint"
  )
  audit_rows[[stage]] <- table[, keep]
}

audit <- do.call(rbind, audit_rows)
rownames(audit) <- NULL
utils::write.csv(audit, output_path, row.names = FALSE, na = "")

changed <- audit[audit$original_accepted != audit$corrected_selected, , drop = FALSE]
cat("Audited ", nrow(audit), " saved candidate fits.\n", sep = "")
cat("Rows whose selection status changes after EDF repair: ", nrow(changed), "\n", sep = "")
if (nrow(changed)) {
  print(changed[, c(
    "stage", "step", "dropped_term", "AIC", "corrected_AIC",
    "original_accepted", "corrected_selected"
  )], row.names = FALSE)
}
cat("Saved audit: ", output_path, "\n", sep = "")
