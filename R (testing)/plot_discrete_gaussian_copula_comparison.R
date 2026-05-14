#!/usr/bin/env Rscript

suppressPackageStartupMessages({
  library(ggplot2)
})

out_dir <- file.path("results", "discrete_gaussian_copula_comparison")
csv_final <- file.path(out_dir, "discrete_gaussian_copula_comparison_results.csv")
csv_partial <- file.path(out_dir, "discrete_gaussian_copula_comparison_results_partial.csv")
csv_in <- if (file.exists(csv_partial)) csv_partial else csv_final

if (!file.exists(csv_in)) {
  stop("No discrete comparison results CSV found in ", out_dir)
}

res <- read.csv(csv_in, stringsAsFactors = FALSE)
res$method_label <- ifelse(res$model == "gamlss", "gamlss", res$method)
res$method_label[res$model == "truth"] <- "truth"
res$method_label <- factor(
  res$method_label,
  levels = c("truth", "gamlss", "rectangle", "midpoint", "distributional_transform")
)

fit_res <- res[res$model != "truth" & isTRUE(TRUE), , drop = FALSE]
truth <- res[res$model == "truth", , drop = FALSE]

coef_cols <- c("mu_intercept", "mu_x", "mu_t", "sigma_intercept", "sigma_x", "theta")
coef_long <- do.call(rbind, lapply(coef_cols, function(col) {
  data.frame(
    family = fit_res$family,
    replicate = fit_res$replicate,
    method_label = fit_res$method_label,
    coefficient = col,
    estimate = fit_res[[col]],
    stringsAsFactors = FALSE
  )
}))
coef_long <- coef_long[is.finite(coef_long$estimate), , drop = FALSE]

truth_long <- do.call(rbind, lapply(coef_cols, function(col) {
  data.frame(
    family = truth$family,
    coefficient = col,
    truth = truth[[col]],
    stringsAsFactors = FALSE
  )
}))
truth_long <- truth_long[is.finite(truth_long$truth), , drop = FALSE]

p_coef <- ggplot(coef_long, aes(x = method_label, y = estimate, color = method_label)) +
  geom_hline(
    data = truth_long,
    aes(yintercept = truth),
    inherit.aes = FALSE,
    linetype = 2,
    color = "black"
  ) +
  geom_point(position = position_jitter(width = 0.12, height = 0), alpha = 0.75, size = 2) +
  stat_summary(fun = mean, geom = "crossbar", width = 0.55, color = "black", fatten = 0) +
  facet_grid(coefficient ~ family, scales = "free_y") +
  labs(x = NULL, y = "Estimate", color = NULL, title = "Discrete Gaussian Copula Simulation: Parameter Estimates") +
  theme_bw(base_size = 11) +
  theme(axis.text.x = element_text(angle = 35, hjust = 1), legend.position = "bottom")

theta_long <- coef_long[coef_long$coefficient == "theta", , drop = FALSE]
theta_truth <- truth_long[truth_long$coefficient == "theta", , drop = FALSE]

p_theta <- ggplot(theta_long, aes(x = method_label, y = estimate, color = method_label)) +
  geom_hline(
    data = theta_truth,
    aes(yintercept = truth),
    inherit.aes = FALSE,
    linetype = 2,
    color = "black"
  ) +
  geom_point(position = position_jitter(width = 0.12, height = 0), alpha = 0.8, size = 2.2) +
  stat_summary(fun = mean, geom = "crossbar", width = 0.55, color = "black", fatten = 0) +
  facet_wrap(~ family, nrow = 1) +
  labs(x = NULL, y = "Gaussian copula rho", color = NULL, title = "Copula Parameter Estimates") +
  theme_bw(base_size = 11) +
  theme(axis.text.x = element_text(angle = 35, hjust = 1), legend.position = "bottom")

loglik_res <- fit_res[is.finite(fit_res$logLik), , drop = FALSE]
p_loglik <- ggplot(loglik_res, aes(x = method_label, y = logLik, color = method_label)) +
  geom_point(position = position_jitter(width = 0.12, height = 0), alpha = 0.8, size = 2.2) +
  stat_summary(fun = mean, geom = "crossbar", width = 0.55, color = "black", fatten = 0) +
  facet_wrap(~ family, scales = "free_y") +
  labs(x = NULL, y = "Joint log-likelihood", color = NULL, title = "Fitted Joint Log-Likelihood") +
  theme_bw(base_size = 11) +
  theme(axis.text.x = element_text(angle = 35, hjust = 1), legend.position = "bottom")

coef_png <- file.path(out_dir, "discrete_parameter_estimates.png")
theta_png <- file.path(out_dir, "discrete_theta_estimates.png")
loglik_png <- file.path(out_dir, "discrete_loglik_comparison.png")

ggsave(coef_png, p_coef, width = 13, height = 10, dpi = 180)
ggsave(theta_png, p_theta, width = 10, height = 4.8, dpi = 180)
ggsave(loglik_png, p_loglik, width = 10, height = 4.8, dpi = 180)

cat("Read results from: ", csv_in, "\n", sep = "")
cat("Wrote plots:\n")
cat("  ", coef_png, "\n", sep = "")
cat("  ", theta_png, "\n", sep = "")
cat("  ", loglik_png, "\n", sep = "")
