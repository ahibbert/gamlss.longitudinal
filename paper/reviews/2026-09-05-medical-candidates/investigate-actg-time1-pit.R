#!/usr/bin/env Rscript

# Read-only diagnostic audit of the preserved ACTG ZISICHEL--Gumbel fit.
# Outputs are deliberately written to a separate investigation directory.

script_argument <- grep("^--file=", commandArgs(trailingOnly = FALSE), value = TRUE)
script_path <- if (length(script_argument)) {
  normalizePath(sub("^--file=", "", script_argument[[1L]]), mustWork = TRUE)
} else {
  normalizePath(file.path(
    getwd(), "paper", "reviews", "2026-09-05-medical-candidates",
    "investigate-actg-time1-pit.R"
  ), mustWork = TRUE)
}
package_root <- normalizePath(file.path(dirname(script_path), "..", "..", ".."))
review_root <- file.path(
  dirname(script_path), "results", "actg-revised-model-screen", "review-full-inference"
)
output_dir <- file.path(dirname(script_path), "results", "actg-time1-investigation")
dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)

if (!requireNamespace("devtools", quietly = TRUE)) stop("Install devtools.")
if (!requireNamespace("ggplot2", quietly = TRUE)) stop("Install ggplot2.")
devtools::load_all(package_root, quiet = TRUE)

fit <- readRDS(file.path(review_root, "objects", "final-model.rds"))
observed <- utils::read.csv(file.path(
  review_root, "tables", "observation-level-diagnostics.csv"
))
observed <- observed[is.finite(observed$response), , drop = FALSE]

cvm_uniform <- function(x) {
  x <- sort(x[is.finite(x)])
  n <- length(x)
  if (n < 1L) return(NA_real_)
  1 / (12 * n) + sum((x - (2 * seq_len(n) - 1) / (2 * n))^2)
}

summarise_pit <- function(x) {
  x <- x[is.finite(x)]
  c(
    n = length(x), mean = mean(x), sd = stats::sd(x),
    sd_error = abs(stats::sd(x) - sqrt(1 / 12)),
    lower_10 = mean(x <= 0.1), upper_10 = mean(x >= 0.9),
    cvm = cvm_uniform(x),
    ks_p = suppressWarnings(stats::ks.test(x, "punif")$p.value)
  )
}

# Repeat randomized PIT calculations to distinguish reproducible lack of fit
# from the arbitrary within-mass randomization used for a discrete margin.
seeds <- 20260910L + 0:99
seed_rows <- lapply(seeds, function(seed) {
  pit <- pithist(fit, randomize = TRUE, seed = seed, plot = FALSE)$pit
  if (length(pit) != nrow(observed)) stop("PIT rows did not align with observed rows.")
  pieces <- split(pit, observed$time)
  do.call(rbind, lapply(names(pieces), function(time) {
    data.frame(seed = seed, time = as.integer(time), t(summarise_pit(pieces[[time]])))
  }))
})
seed_results <- do.call(rbind, seed_rows)
utils::write.csv(seed_results, file.path(output_dir, "pit-by-time-100-seeds.csv"), row.names = FALSE)

seed_summary <- do.call(rbind, lapply(split(seed_results, seed_results$time), function(value) {
  data.frame(
    time = value$time[[1L]], n = value$n[[1L]],
    pit_mean_average = mean(value$mean), pit_mean_range_low = min(value$mean),
    pit_mean_range_high = max(value$mean), pit_sd_average = mean(value$sd),
    pit_sd_range_low = min(value$sd), pit_sd_range_high = max(value$sd),
    cvm_average = mean(value$cvm), ks_rejection_rate_005 = mean(value$ks_p < 0.05),
    lower_10_average = mean(value$lower_10), upper_10_average = mean(value$upper_10)
  )
}))
rownames(seed_summary) <- NULL
utils::write.csv(seed_summary, file.path(output_dir, "pit-by-time-seed-summary.csv"), row.names = FALSE)

# Use the fixed review seed for covariate localization.
pit <- pithist(fit, randomize = TRUE, seed = 20260910L, plot = FALSE)$pit
observed$pit <- pit
observed$qres <- stats::qnorm(pmin(pmax(pit, .Machine$double.eps), 1 - .Machine$double.eps))
time1 <- observed[observed$time == min(observed$time), , drop = FALSE]
time1$age_group <- cut(time1$age, unique(stats::quantile(time1$age, 0:4 / 4)), include.lowest = TRUE)
time1$baseline_group <- cut(
  time1$baseline_cd4,
  unique(stats::quantile(time1$baseline_cd4, 0:4 / 4)), include.lowest = TRUE
)

stratum_summary <- function(variable) {
  pieces <- split(time1$pit, time1[[variable]], drop = TRUE)
  do.call(rbind, lapply(names(pieces), function(level) {
    data.frame(variable = variable, level = level, t(summarise_pit(pieces[[level]])))
  }))
}
strata <- do.call(rbind, lapply(
  c("treatment", "sex", "age_group", "baseline_group"), stratum_summary
))
rownames(strata) <- NULL
utils::write.csv(strata, file.path(output_dir, "time1-pit-by-covariate.csv"), row.names = FALSE)

# Location and scale regressions are screening diagnostics, not final inference.
# qres tests conditional location; qres^2 tests conditional residual spread.
time1$age_z_local <- as.numeric(scale(time1$age))
time1$baseline_log_z_local <- as.numeric(scale(log1p(time1$baseline_cd4)))
location_fit <- stats::lm(
  qres ~ treatment + age_z_local + sex + baseline_log_z_local, data = time1
)
scale_fit <- stats::lm(
  I(qres^2) ~ treatment + age_z_local + sex + baseline_log_z_local, data = time1
)
tidy_lm <- function(model, diagnostic) {
  value <- as.data.frame(summary(model)$coefficients)
  value$term <- rownames(value)
  rownames(value) <- NULL
  names(value)[1:4] <- c("estimate", "std_error", "statistic", "p_value")
  value$diagnostic <- diagnostic
  value[, c("diagnostic", "term", "estimate", "std_error", "statistic", "p_value")]
}
regression_table <- rbind(tidy_lm(location_fit, "qres_location"), tidy_lm(scale_fit, "qres_squared_scale"))
utils::write.csv(regression_table, file.path(output_dir, "time1-residual-regressions.csv"), row.names = FALSE)

plot_data <- time1
plot_data$baseline_group <- factor(plot_data$baseline_group)
plot_data$age_group <- factor(plot_data$age_group)
long <- rbind(
  data.frame(variable = "Treatment", level = plot_data$treatment, pit = plot_data$pit),
  data.frame(variable = "Sex", level = plot_data$sex, pit = plot_data$pit),
  data.frame(variable = "Age quartile", level = plot_data$age_group, pit = plot_data$pit),
  data.frame(variable = "Baseline CD4 quartile", level = plot_data$baseline_group, pit = plot_data$pit)
)
p <- ggplot2::ggplot(long, ggplot2::aes(x = pit)) +
  ggplot2::geom_histogram(
    ggplot2::aes(y = ggplot2::after_stat(density)),
    bins = 10, boundary = 0, colour = "white", fill = "#2C7FB8"
  ) +
  ggplot2::facet_grid(variable ~ level, scales = "free_y") +
  ggplot2::geom_hline(yintercept = 1, linetype = "dashed", colour = "#D95F02") +
  ggplot2::labs(
    title = "Visit-1 randomized PIT by observed covariates",
    subtitle = "Dashed reference denotes equal mass across ten PIT bins",
    x = "Randomized PIT", y = "Count"
  ) +
  ggplot2::theme_minimal(base_size = 11)
ggplot2::ggsave(file.path(output_dir, "time1-pit-by-covariate.png"), p, width = 13, height = 10, dpi = 180)

saveRDS(list(
  generated_at = Sys.time(), preserved_fit = normalizePath(file.path(review_root, "objects", "final-model.rds")),
  seed_summary = seed_summary, strata = strata,
  location_model = location_fit, scale_model = scale_fit
), file.path(output_dir, "time1-audit.rds"))

cat("Time-1 audit written to: ", normalizePath(output_dir), "\n", sep = "")
print(seed_summary)
cat("\nResidual regression terms with p < 0.05:\n")
print(regression_table[regression_table$p_value < 0.05, , drop = FALSE])
