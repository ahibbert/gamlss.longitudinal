#!/usr/bin/env Rscript

# plot_terms() deliberately omits continuous-by-indicator interactions. Add a
# publication-quality contribution plot for the shortlisted visit-1 interaction.

if (!requireNamespace("ggplot2", quietly = TRUE)) stop("Install ggplot2.")
script_argument <- grep("^--file=", commandArgs(trailingOnly = FALSE), value = TRUE)
script_path <- if (length(script_argument)) {
  normalizePath(sub("^--file=", "", script_argument[[1L]]), mustWork = TRUE)
} else {
  normalizePath("paper/reviews/actg-shortlist/add-interaction-term-plots.R", mustWork = TRUE)
}
root <- dirname(script_path)

for (model in c("02-mu-baseline-visit1", "03-mu-baseline-visit1-nu-visit1")) {
  model_dir <- file.path(root, model)
  fit <- readRDS(file.path(model_dir, "objects", "final-model-with-sandwich-vcov.rds"))
  coefficients <- utils::read.csv(
    file.path(model_dir, "tables", "coefficient-table-numeric.csv"),
    stringsAsFactors = FALSE, check.names = FALSE
  )
  row <- coefficients[
    coefficients$parameter == "mu" & coefficients$term == "baseline_cd4:visit1",
    , drop = FALSE
  ]
  if (nrow(row) != 1L) stop("Could not resolve the visit-1 baseline-CD4 coefficient for ", model)
  observed_baseline <- fit$dataset_original$baseline_cd4
  grid <- data.frame(
    baseline_cd4 = seq(min(observed_baseline, na.rm = TRUE), max(observed_baseline, na.rm = TRUE), length.out = 250L)
  )
  grid$effect <- row$estimate * grid$baseline_cd4
  grid$std_error <- abs(grid$baseline_cd4) * row$robust_std_error
  grid$conf_low <- grid$effect - stats::qnorm(0.975) * grid$std_error
  grid$conf_high <- grid$effect + stats::qnorm(0.975) * grid$std_error

  p <- ggplot2::ggplot(grid, ggplot2::aes(x = baseline_cd4, y = effect)) +
    ggplot2::geom_ribbon(
      ggplot2::aes(ymin = conf_low, ymax = conf_high),
      fill = "#56B4E9", alpha = 0.22
    ) +
    ggplot2::geom_line(colour = "#0072B2", linewidth = 1) +
    ggplot2::geom_hline(yintercept = 0, linetype = "dashed", colour = "#555555") +
    ggplot2::labs(
      title = "Additional baseline-CD4 effect at visit 1",
      subtitle = "Contribution to the mu linear predictor; 95% subject-clustered sandwich CI",
      x = "Baseline CD4", y = "Additional visit-1 contribution"
    ) +
    ggplot2::theme_minimal(base_size = 12)

  stem <- file.path(
    model_dir, "figures",
    "term-effects-fixed-terms-mu-baseline-cd4-visit1-plot"
  )
  ggplot2::ggsave(paste0(stem, ".png"), p, width = 7, height = 5, dpi = 320, bg = "white")
  ggplot2::ggsave(paste0(stem, ".pdf"), p, width = 7, height = 5, bg = "white")
}

cat("Interaction term plots written under ", root, "\n", sep = "")
