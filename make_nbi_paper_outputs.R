#!/usr/bin/env Rscript

suppressPackageStartupMessages({
  library(dplyr)
  library(ggplot2)
  library(knitr)
  library(tidyr)
})

root <- normalizePath(".", winslash = "/", mustWork = TRUE)
comparison_dir <- file.path(root, "results", "del_clayton_rectangle_diagnostics", "nbi_highsignal_n500_rep100_comparison")
prefix <- "paper_nbi_highsignal_rs_comparison"

model_labels <- c(
  "gamlss" = "gamlss",
  "ours_rs_separate" = "RS separate",
  "ours_rs_joint" = "RS joint"
)

method_order <- unname(model_labels)
method_colours <- c(
  "gamlss" = "#D55E00",
  "RS separate" = "#0072B2",
  "RS joint" = "#009E73"
)
method_fills <- c(
  "gamlss" = "#E69F00",
  "RS separate" = "#56B4E9",
  "RS joint" = "#7CAE00"
)

parameter_order <- c("mu", "sigma", "theta")

read_result <- function(name) {
  path <- file.path(comparison_dir, name)
  if (!file.exists(path)) {
    stop("Missing required result file: ", path, call. = FALSE)
  }
  read.csv(path, stringsAsFactors = FALSE)
}

as_method <- function(x) {
  out <- unname(model_labels[x])
  ifelse(is.na(out), x, out)
}

format_num <- function(x, digits = 3) {
  ifelse(is.na(x) | !is.finite(x), "--", formatC(x, digits = digits, format = "f"))
}

format_signed <- function(x, digits = 3) {
  ifelse(is.na(x) | !is.finite(x), "--", sprintf(paste0("%+.", digits, "f"), x))
}

format_mean_sd <- function(mean, sd, digits = 3) {
  ifelse(
    is.na(mean) | !is.finite(mean),
    "--",
    paste0(format_num(mean, digits), " (", format_num(sd, digits), ")")
  )
}

latex_escape <- function(x) {
  x <- gsub("\\\\", "\\\\textbackslash{}", x)
  x <- gsub("_", "\\\\_", x)
  x
}

param_latex <- function(x) {
  recode(
    x,
    "mu" = "$\\mu$",
    "sigma" = "$\\sigma$",
    "theta" = "$\\theta$",
    .default = latex_escape(x)
  )
}

term_latex <- function(x) {
  recode(
    x,
    "intercept" = "$\\beta_0$",
    "time_scaled" = "$t$",
    "x1" = "$x_1$",
    "x2" = "$x_2$",
    "smooth" = "$s_1$",
    .default = latex_escape(x)
  )
}

write_latex_table <- function(data, path, caption, label, align = NULL) {
  tex <- kable(
    data,
    format = "latex",
    booktabs = TRUE,
    escape = FALSE,
    caption = caption,
    label = label,
    align = align,
    linesep = ""
  )
  writeLines(tex, path, useBytes = TRUE)
}

theme_paper <- function(base_size = 11) {
  theme_minimal(base_size = base_size) +
    theme(
      plot.background = element_rect(fill = "white", colour = NA),
      panel.background = element_rect(fill = "white", colour = NA),
      panel.grid.minor = element_blank(),
      panel.grid.major.x = element_blank(),
      panel.grid.major.y = element_blank(),
      axis.text.x = element_text(angle = 35, hjust = 1, vjust = 1),
      legend.position = "top",
      legend.title = element_blank(),
      strip.text = element_text(face = "bold"),
      strip.text.y = element_text(face = "bold", angle = 0),
      strip.text.y.right = element_text(face = "bold", angle = 0),
      plot.title.position = "plot"
    )
}

fixed_by_rep <- read_result("fixed_effects_by_rep.csv") %>%
  mutate(
    method = factor(as_method(model), levels = method_order),
    parameter = factor(parameter, levels = parameter_order),
    term = factor(term, levels = c("intercept", "time_scaled", "x1", "x2")),
    term_math = recode(
      as.character(term),
      "intercept" = "beta[0]",
      "time_scaled" = "t",
      "x1" = "x[1]",
      "x2" = "x[2]"
    ),
    error = estimate - true_value
  ) %>%
  filter(is.finite(error))

fixed_summary <- fixed_by_rep %>%
  group_by(method, parameter, term, term_math) %>%
  summarise(
    mean_error = mean(error, na.rm = TRUE),
    lower_error = quantile(error, 0.05, na.rm = TRUE),
    upper_error = quantile(error, 0.95, na.rm = TRUE),
    .groups = "drop"
  ) %>%
  mutate(term_math = factor(term_math, levels = c("beta[0]", "t", "x[1]", "x[2]")))

fixed_limits <- fixed_summary %>%
  group_by(parameter) %>%
  summarise(limit = max(abs(c(lower_error, upper_error)), na.rm = TRUE), .groups = "drop") %>%
  mutate(term_math = factor("beta[0]", levels = c("beta[0]", "t", "x[1]", "x[2]"))) %>%
  tidyr::pivot_longer(cols = limit, names_to = "bound", values_to = "bias") %>%
  bind_rows(
    fixed_summary %>%
      group_by(parameter) %>%
      summarise(limit = -max(abs(c(lower_error, upper_error)), na.rm = TRUE), .groups = "drop") %>%
      mutate(
        bound = "lower",
        bias = limit,
        term_math = factor("beta[0]", levels = c("beta[0]", "t", "x[1]", "x[2]"))
      ) %>%
      select(parameter, term_math, bound, bias)
  )

fixed_plot <- ggplot(fixed_summary, aes(x = term_math, y = mean_error, colour = method)) +
  geom_blank(data = fixed_limits, aes(x = term_math, y = bias), inherit.aes = FALSE) +
  geom_hline(yintercept = 0, linewidth = 0.35, colour = "grey30") +
  geom_pointrange(
    aes(ymin = lower_error, ymax = upper_error),
    position = position_dodge(width = 0.55),
    linewidth = 0.5,
    fatten = 2.4
  ) +
  facet_wrap(~ parameter, scales = "free_y", nrow = 1, labeller = label_parsed) +
  scale_colour_manual(values = method_colours, drop = FALSE) +
  scale_x_discrete(labels = function(x) parse(text = x)) +
  labs(x = NULL, y = "Bias") +
  theme_paper()

ggsave(
  file.path(comparison_dir, paste0(prefix, "_fixed_effect_recovery.png")),
  fixed_plot,
  width = 12,
  height = 3.8,
  dpi = 320,
  bg = "white"
)

smooth_by_rep <- read_result("smooth_estimates_by_rep.csv") %>%
  mutate(
    method = factor(as_method(model), levels = method_order),
    parameter = factor(parameter, levels = parameter_order),
    s1_plot = round(s1, 2)
  ) %>%
  group_by(method, parameter, rep, s1_plot) %>%
  summarise(
    smooth_hat = mean(smooth_hat, na.rm = TRUE),
    smooth_true = mean(smooth_true, na.rm = TRUE),
    .groups = "drop"
  ) %>%
  rename(s1 = s1_plot)

smooth_summary <- smooth_by_rep %>%
  group_by(method, parameter, s1) %>%
  summarise(
    mean_hat = mean(smooth_hat, na.rm = TRUE),
    q25_hat = quantile(smooth_hat, 0.25, na.rm = TRUE),
    q75_hat = quantile(smooth_hat, 0.75, na.rm = TRUE),
    true = mean(smooth_true, na.rm = TRUE),
    .groups = "drop"
  )

truth_curves <- smooth_summary %>%
  distinct(method, parameter, s1, true)

smooth_plot <- ggplot(smooth_summary, aes(x = s1)) +
  geom_ribbon(aes(ymin = q25_hat, ymax = q75_hat, fill = method), colour = NA, alpha = 0.18) +
  geom_line(aes(y = mean_hat, colour = method), linewidth = 0.8) +
  geom_line(
    data = truth_curves,
    aes(x = s1, y = true),
    inherit.aes = FALSE,
    colour = "black",
    linewidth = 0.9,
    linetype = "dashed"
  ) +
  facet_grid(
    method ~ parameter,
    scales = "free_y",
    drop = FALSE,
    labeller = labeller(parameter = label_parsed, method = label_value)
  ) +
  scale_colour_manual(values = method_colours, drop = FALSE) +
  scale_fill_manual(values = method_fills, drop = FALSE) +
  labs(x = NULL, y = "Smooth effect", caption = "Ribbon shows interquartile range across replicates; dashed line is truth.") +
  theme_paper() +
  theme(strip.text.y = element_blank(), strip.text.y.right = element_blank())

ggsave(
  file.path(comparison_dir, paste0(prefix, "_smooth_effect_recovery.png")),
  smooth_plot,
  width = 11,
  height = 6.5,
  dpi = 320,
  bg = "white"
)

fixed_table_long <- read_result("fixed_effects_bias_rmse_table.csv") %>%
  mutate(
    Method = as_method(model),
    Parameter = as.character(parameter),
    Term = term
  ) %>%
  select(Parameter, Term, Method, Bias = bias, RMSE = rmse)

fixed_wide <- fixed_table_long %>%
  pivot_wider(
    names_from = Method,
    values_from = c(Bias, RMSE),
    names_glue = "{Method}_{.value}"
  )

smooth_table_long <- read_result("smooth_integrated_metrics.csv") %>%
  mutate(
    Method = as_method(model),
    Parameter = as.character(parameter)
  ) %>%
  group_by(Parameter, Method) %>%
  summarise(`Smooth IRMSE` = mean(irmse, na.rm = TRUE), .groups = "drop")

smooth_wide <- smooth_table_long %>%
  pivot_wider(
    names_from = Method,
    values_from = `Smooth IRMSE`,
    names_glue = "{Method}_Smooth_IRMSE"
  ) %>%
  mutate(Term = "smooth")

recovery_table_csv <- bind_rows(
  fixed_wide %>% mutate(
    gamlss_Smooth_IRMSE = NA_real_,
    `RS separate_Smooth_IRMSE` = NA_real_,
    `RS joint_Smooth_IRMSE` = NA_real_
  ),
  smooth_wide %>% mutate(
    gamlss_Bias = NA_real_,
    gamlss_RMSE = NA_real_,
    `RS separate_Bias` = NA_real_,
    `RS separate_RMSE` = NA_real_,
    `RS joint_Bias` = NA_real_,
    `RS joint_RMSE` = NA_real_
  )
) %>%
  arrange(
    factor(Parameter, levels = parameter_order),
    factor(Term, levels = c("intercept", "time_scaled", "x1", "x2", "smooth"))
  )

write.csv(
  recovery_table_csv,
  file.path(comparison_dir, paste0(prefix, "_fixed_parameter_bias_rmse.csv")),
  row.names = FALSE
)

recovery_latex <- recovery_table_csv %>%
  mutate(
    `gamlss RMSE/IRMSE` = ifelse(Term == "smooth", gamlss_Smooth_IRMSE, gamlss_RMSE),
    `RS separate RMSE/IRMSE` = ifelse(Term == "smooth", `RS separate_Smooth_IRMSE`, `RS separate_RMSE`),
    `RS joint RMSE/IRMSE` = ifelse(Term == "smooth", `RS joint_Smooth_IRMSE`, `RS joint_RMSE`)
  ) %>%
  transmute(
    Parameter = param_latex(Parameter),
    Term = term_latex(Term),
    `gamlss Bias` = format_num(gamlss_Bias),
    `gamlss RMSE` = format_num(`gamlss RMSE/IRMSE`),
    `RS separate Bias` = format_num(`RS separate_Bias`),
    `RS separate RMSE` = format_num(`RS separate RMSE/IRMSE`),
    `RS joint Bias` = format_num(`RS joint_Bias`),
    `RS joint RMSE` = format_num(`RS joint RMSE/IRMSE`)
  )

write_latex_table(
  recovery_latex,
  file.path(comparison_dir, paste0(prefix, "_fixed_parameter_bias_rmse.tex")),
  caption = "NBI high-signal parameter recovery. Smooth rows report IRMSE in the RMSE column.",
  label = "tab:nbi-parameter-recovery"
)

summarise_metric <- function(data, metric, label, digits = 3) {
  data %>%
    group_by(Method) %>%
    summarise(
      mean = mean(.data[[metric]], na.rm = TRUE),
      sd = sd(.data[[metric]], na.rm = TRUE),
      .groups = "drop"
    ) %>%
    mutate(Metric = label, Value = format_mean_sd(mean, sd, digits)) %>%
    select(Metric, Method, Value)
}

run_logs <- read_result("nbi_sigma_compare_logs.csv") %>%
  mutate(Method = as_method(engine))

fixed_rmse_by_rep <- fixed_by_rep %>%
  group_by(Method = method, rep) %>%
  summarise(fixed_rmse = sqrt(mean(error^2, na.rm = TRUE)), .groups = "drop") %>%
  mutate(Method = as.character(Method))

smooth_irmse_by_rep <- read_result("smooth_integrated_metrics.csv") %>%
  mutate(Method = as_method(model)) %>%
  group_by(Method, rep) %>%
  summarise(smooth_irmse = mean(irmse, na.rm = TRUE), .groups = "drop")

fit_table_long <- bind_rows(
  summarise_metric(fixed_rmse_by_rep, "fixed_rmse", "Fixed-effect RMSE", 3),
  summarise_metric(smooth_irmse_by_rep, "smooth_irmse", "Smooth IRMSE", 3),
  summarise_metric(run_logs, "elapsed_sec", "Runtime (s)", 1),
  summarise_metric(run_logs, "iter", "Outer iterations", 1)
)

fit_table_csv <- fit_table_long %>%
  mutate(Method = factor(Method, levels = method_order)) %>%
  pivot_wider(names_from = Method, values_from = Value) %>%
  select(Metric, all_of(method_order))

write.csv(
  fit_table_csv,
  file.path(comparison_dir, paste0(prefix, "_fit_characteristics.csv")),
  row.names = FALSE
)

write_latex_table(
  fit_table_csv,
  file.path(comparison_dir, paste0(prefix, "_fit_characteristics.tex")),
  caption = "NBI high-signal fit characteristics. Values are mean (SD) across replicates.",
  label = "tab:nbi-fit-characteristics"
)

runtime_plot <- ggplot(run_logs, aes(x = Method, y = elapsed_sec, fill = Method)) +
  geom_boxplot(width = 0.6, outlier.alpha = 0.45) +
  scale_fill_manual(values = method_fills, drop = FALSE) +
  labs(x = NULL, y = "Runtime (s)") +
  theme_paper() +
  theme(legend.position = "none")

ggsave(
  file.path(comparison_dir, paste0(prefix, "_runtime_comparison.png")),
  runtime_plot,
  width = 7,
  height = 4.2,
  dpi = 320,
  bg = "white"
)

iteration_plot <- ggplot(run_logs, aes(x = Method, y = iter, fill = Method)) +
  geom_boxplot(width = 0.6, outlier.alpha = 0.45) +
  scale_fill_manual(values = method_fills, drop = FALSE) +
  labs(x = NULL, y = "Outer iterations") +
  theme_paper() +
  theme(legend.position = "none")

ggsave(
  file.path(comparison_dir, paste0(prefix, "_iteration_comparison.png")),
  iteration_plot,
  width = 7,
  height = 4.2,
  dpi = 320,
  bg = "white"
)

convergence_counts <- run_logs %>%
  group_by(Method, success, converged) %>%
  summarise(n = n(), .groups = "drop")
write.csv(
  convergence_counts,
  file.path(comparison_dir, paste0(prefix, "_convergence_counts.csv")),
  row.names = FALSE
)

message("Wrote NBI paper outputs to: ", comparison_dir)
