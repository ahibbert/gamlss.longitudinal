#!/usr/bin/env Rscript

suppressPackageStartupMessages({
  library(dplyr)
  library(ggplot2)
  library(knitr)
  library(tidyr)
})

root <- normalizePath(".", winslash = "/", mustWork = TRUE)
input_dir_override <- Sys.getenv("NBI_PAPER_INPUT_DIR", unset = "")
output_dir_override <- Sys.getenv("NBI_PAPER_OUTPUT_DIR", unset = "")
comparison_dir_override <- Sys.getenv("NBI_PAPER_COMPARISON_DIR", unset = "")
if (nzchar(input_dir_override)) {
  input_dir <- normalizePath(input_dir_override, winslash = "/", mustWork = TRUE)
  comparison_dir <- if (nzchar(output_dir_override)) output_dir_override else input_dir
  dir.create(comparison_dir, recursive = TRUE, showWarnings = FALSE)
  comparison_dir <- normalizePath(comparison_dir, winslash = "/", mustWork = TRUE)
} else if (nzchar(comparison_dir_override)) {
  input_dir <- normalizePath(comparison_dir_override, winslash = "/", mustWork = TRUE)
  comparison_dir <- input_dir
} else {
  comparison_dir <- file.path(
    root,
    "results",
    "jss-exploratory",
    "02-discrete-delaporte-clayton",
    "nbi_highsignal_n500_rep100_se_diagnostics"
  )
  if (!dir.exists(comparison_dir)) {
    comparison_dir <- file.path(
      root,
      "results",
      "jss-exploratory",
      "02-discrete-delaporte-clayton",
      "nbi_highsignal_n500_rep100_comparison"
    )
  }
  if (!dir.exists(comparison_dir)) {
    comparison_dir <- file.path(
      root,
      "results",
      "del_clayton_rectangle_diagnostics",
      "nbi_highsignal_n500_rep100_comparison"
    )
  }
  input_dir <- comparison_dir
}
prefix <- "paper_simulation_nbi_clayton_highsignal"

model_labels <- c(
  "gamlss" = "gamlss2",
  "ours_rs_joint" = "gamlss.longitudinal"
)

method_order <- unname(model_labels)
method_colours <- c(
  "gamlss2" = "#D55E00",
  "gamlss.longitudinal" = "#009E73"
)
method_fills <- c(
  "gamlss2" = "#E69F00",
  "gamlss.longitudinal" = "#7CAE00"
)

parameter_order <- c("mu", "sigma", "theta")
marginal_parameters <- c("mu", "sigma")

read_result <- function(name) {
  path <- file.path(input_dir, name)
  if (!file.exists(path)) {
    stop("Missing required result file: ", path, call. = FALSE)
  }
  read.csv(path, stringsAsFactors = FALSE)
}

read_optional_result <- function(name) {
  path <- file.path(input_dir, name)
  if (!file.exists(path)) {
    return(data.frame(model = character()))
  }
  read.csv(path, stringsAsFactors = FALSE)
}

as_method <- function(x) {
  out <- unname(model_labels[x])
  ifelse(is.na(out), x, out)
}

adjust_gamlss_mu_intercept <- function(data) {
  idx <- data$model == "gamlss" & data$parameter == "mu" & data$term == "intercept"
  # The retained gamlss coefficient stores the smooth centering shift in beta0.
  shift <- mean(data$estimate[idx] - data$true_value[idx], na.rm = TRUE)
  if (is.finite(shift)) {
    data$estimate[idx] <- data$estimate[idx] - shift
  }
  attr(data, "gamlss_mu_intercept_centering_shift") <- shift
  data
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

write_recovery_latex_table <- function(data, path) {
  lines <- c(
    "\\begin{table}",
    "\\centering",
    "\\caption{\\label{tab:nbi-parameter-recovery}NBI high-signal parameter recovery for the gamlss2 reference and gamlss.longitudinal model. Smooth rows report IRMSE in the RMSE column.}",
    "\\begin{tabular}[t]{llrrrrrr}",
    "\\toprule",
    "Parameter & Term & \\multicolumn{3}{c}{gamlss2} & \\multicolumn{3}{c}{gamlss.longitudinal}\\\\",
    "\\cmidrule(lr){3-5} \\cmidrule(lr){6-8}",
    " &  & Bias & RMSE & 95\\% CI & Bias & RMSE & 95\\% CI\\\\",
    "\\midrule"
  )

  body <- character()
  previous_parameter <- NULL
  for (i in seq_len(nrow(data))) {
    row <- data[i, , drop = FALSE]
    current_parameter <- row[["Parameter"]]
    if (!is.null(previous_parameter) && !identical(current_parameter, previous_parameter)) {
      body <- c(body, "\\addlinespace")
    }

    parameter_cell <- if (is.null(previous_parameter) || !identical(current_parameter, previous_parameter)) {
      param_latex(current_parameter)
    } else {
      ""
    }

    body <- c(
      body,
      paste(
        parameter_cell,
        term_latex(row[["Term"]]),
        row[["gamlss2 bias"]],
        row[["gamlss2 RMSE"]],
        row[["gamlss2 coverage"]],
        row[["gamlss.longitudinal bias"]],
        row[["gamlss.longitudinal RMSE"]],
        row[["gamlss.longitudinal coverage"]],
        sep = " & "
      ) |>
        paste0("\\\\")
    )
    previous_parameter <- current_parameter
  }

  lines <- c(lines, body, "\\bottomrule", "\\end{tabular}", "\\end{table}")
  writeLines(lines, path, useBytes = TRUE)
}

write_fit_latex_table <- function(data, path) {
  lines <- c(
    "\\begin{table}",
    "\\centering",
    "\\caption{\\label{tab:nbi-fit-characteristics}NBI high-signal fit and predictive characteristics for the gamlss2 reference and gamlss.longitudinal model. Values are mean (SD) across replicates where retained. Marginal fixed-effect RMSE excludes the dependence parameter $\\theta$; marginal smooth IRMSE uses only $\\mu$ and $\\sigma$ smooths.}",
    "\\begin{tabular}[t]{lrr}",
    "\\toprule",
    "Metric & gamlss2 & gamlss.longitudinal\\\\",
    "\\midrule"
  )

  body <- apply(data, 1, function(row) {
    paste(latex_escape(row[["Metric"]]), row[["gamlss2"]], row[["gamlss.longitudinal"]], sep = " & ") |>
      paste0("\\\\")
  })

  lines <- c(lines, body, "\\bottomrule", "\\end{tabular}", "\\end{table}")
  writeLines(lines, path, useBytes = TRUE)
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
  filter(model %in% names(model_labels)) %>%
  adjust_gamlss_mu_intercept() %>%
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
  filter(model %in% names(model_labels)) %>%
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

fixed_has_se <- "std_error" %in% names(fixed_by_rep)
fixed_table_long <- fixed_by_rep %>%
  group_by(method, parameter, term) %>%
  summarise(
    Bias = mean(error, na.rm = TRUE),
    RMSE = sqrt(mean(error^2, na.rm = TRUE)),
    Coverage = if (fixed_has_se) {
      mean(
        true_value >= estimate - stats::qnorm(0.975) * std_error &
          true_value <= estimate + stats::qnorm(0.975) * std_error,
        na.rm = TRUE
      )
    } else {
      NA_real_
    },
    .groups = "drop"
  ) %>%
  mutate(
    Method = as.character(method),
    Parameter = as.character(parameter),
    Term = as.character(term)
  ) %>%
  select(Parameter, Term, Method, Bias, RMSE, Coverage)

fixed_wide <- fixed_table_long %>%
  pivot_wider(
    names_from = Method,
    values_from = c(Bias, RMSE, Coverage),
    names_glue = "{Method}_{.value}"
  )

smooth_table_long <- read_result("smooth_integrated_metrics.csv") %>%
  filter(model %in% names(model_labels)) %>%
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
  transmute(
    Parameter,
    Term = "smooth",
    gamlss2_Bias = NA_real_,
    gamlss2_RMSE = NA_real_,
    gamlss2_Coverage = NA_real_,
    `gamlss.longitudinal_Bias` = NA_real_,
    `gamlss.longitudinal_RMSE` = NA_real_,
    `gamlss.longitudinal_Coverage` = NA_real_,
    gamlss2_Smooth_IRMSE,
    `gamlss.longitudinal_Smooth_IRMSE`
  )

recovery_table_csv <- bind_rows(
  fixed_wide %>% mutate(
    gamlss2_Smooth_IRMSE = NA_real_,
    `gamlss.longitudinal_Smooth_IRMSE` = NA_real_
  ),
  smooth_wide
) %>%
  arrange(
    factor(Parameter, levels = parameter_order),
    factor(Term, levels = c("intercept", "time_scaled", "x1", "x2", "smooth"))
  ) %>%
  select(
    Parameter,
    Term,
    gamlss2_Bias,
    gamlss2_RMSE,
    gamlss2_Coverage,
    gamlss2_Smooth_IRMSE,
    `gamlss.longitudinal_Bias`,
    `gamlss.longitudinal_RMSE`,
    `gamlss.longitudinal_Coverage`,
    `gamlss.longitudinal_Smooth_IRMSE`
  )

write.csv(
  recovery_table_csv,
  file.path(comparison_dir, paste0(prefix, "_fixed_parameter_bias_rmse.csv")),
  row.names = FALSE
)

recovery_latex <- recovery_table_csv %>%
  mutate(
    gamlss2_RMSE_or_IRMSE = ifelse(Term == "smooth", gamlss2_Smooth_IRMSE, gamlss2_RMSE),
    `gamlss.longitudinal_RMSE_or_IRMSE` = ifelse(
      Term == "smooth",
      `gamlss.longitudinal_Smooth_IRMSE`,
      `gamlss.longitudinal_RMSE`
    )
  ) %>%
  transmute(
    Parameter,
    Term,
    `gamlss2 bias` = format_num(gamlss2_Bias),
    `gamlss2 RMSE` = format_num(gamlss2_RMSE_or_IRMSE),
    `gamlss2 coverage` = format_num(gamlss2_Coverage),
    `gamlss.longitudinal bias` = format_num(`gamlss.longitudinal_Bias`),
    `gamlss.longitudinal RMSE` = format_num(`gamlss.longitudinal_RMSE_or_IRMSE`),
    `gamlss.longitudinal coverage` = format_num(`gamlss.longitudinal_Coverage`)
  )

write_recovery_latex_table(
  recovery_latex,
  file.path(comparison_dir, paste0(prefix, "_fixed_parameter_bias_rmse.tex"))
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

summarise_optional_metric <- function(data, metric, label, digits = 3) {
  if (!metric %in% names(data) || !any(is.finite(data[[metric]]))) {
    return(data.frame(
      Metric = label,
      Method = method_order,
      Value = "--",
      stringsAsFactors = FALSE
    ))
  }
  summarise_metric(data, metric, label, digits)
}

format_variogram_p <- function(p_value) {
  format(as.numeric(p_value), trim = TRUE, scientific = FALSE)
}

summarise_variogram_metrics <- function(data, digits = 3) {
  if (!"variogram_score" %in% names(data) || !any(is.finite(data$variogram_score))) {
    return(data.frame(
      Metric = "Variogram score (p=0.5)",
      Method = method_order,
      Value = "--",
      stringsAsFactors = FALSE
    ))
  }
  if (!"variogram_p" %in% names(data)) {
    return(summarise_metric(data, "variogram_score", "Variogram score", digits))
  }
  p_numeric <- as.numeric(data$variogram_p)
  p_values <- sort(unique(p_numeric[is.finite(p_numeric)]))
  bind_rows(lapply(p_values, function(p_value) {
    summarise_metric(
      data %>% filter(abs(as.numeric(variogram_p) - p_value) < 1e-8),
      "variogram_score",
      paste0("Variogram score (p=", format_variogram_p(p_value), ")"),
      digits
    )
  }))
}

run_logs <- read_result("nbi_sigma_compare_logs.csv") %>%
  filter(engine %in% names(model_labels)) %>%
  mutate(Method = as_method(engine))

joint_metrics_by_rep <- read_optional_result("joint_distribution_metrics_by_rep.csv") %>%
  filter(model %in% names(model_labels)) %>%
  mutate(Method = as_method(model))

predictive_scores_by_rep <- read_optional_result("predictive_scores_by_rep.csv") %>%
  filter(model %in% names(model_labels)) %>%
  mutate(Method = as_method(model))

predictive_scores_unique <- predictive_scores_by_rep %>%
  distinct(Method, rep, .keep_all = TRUE)

fixed_rmse_by_rep <- fixed_by_rep %>%
  filter(parameter %in% marginal_parameters) %>%
  group_by(Method = method, rep) %>%
  summarise(fixed_rmse = sqrt(mean(error^2, na.rm = TRUE)), .groups = "drop") %>%
  mutate(Method = as.character(Method))

smooth_irmse_by_rep <- read_result("smooth_integrated_metrics.csv") %>%
  filter(model %in% names(model_labels), parameter %in% marginal_parameters) %>%
  mutate(Method = as_method(model)) %>%
  group_by(Method, rep) %>%
  summarise(smooth_irmse = mean(irmse, na.rm = TRUE), .groups = "drop")

fit_table_long <- bind_rows(
  summarise_optional_metric(run_logs, "logLik", "Log-likelihood", 1),
  summarise_optional_metric(predictive_scores_unique, "test_log_score_per_obs", "Test log score / obs.", 3),
  summarise_variogram_metrics(predictive_scores_by_rep, 3),
  summarise_optional_metric(joint_metrics_by_rep, "abs_rosenblatt_lag1_cor", "Abs. Rosenblatt lag-1 cor.", 3),
  summarise_metric(fixed_rmse_by_rep, "fixed_rmse", "Marginal fixed-effect RMSE", 3),
  summarise_metric(smooth_irmse_by_rep, "smooth_irmse", "Marginal smooth IRMSE", 3),
  summarise_metric(run_logs, "elapsed_sec", "Runtime (s)", 1)
)

variogram_metric_order <- unique(fit_table_long$Metric[grepl("^Variogram score", fit_table_long$Metric)])

fit_table_csv <- fit_table_long %>%
  mutate(Method = factor(Method, levels = method_order)) %>%
  pivot_wider(names_from = Method, values_from = Value) %>%
  select(Metric, all_of(method_order)) %>%
  arrange(
    match(
      Metric,
      c(
        "Log-likelihood",
        "Test log score / obs.",
        variogram_metric_order,
        "Abs. Rosenblatt lag-1 cor.",
        "Marginal fixed-effect RMSE",
        "Marginal smooth IRMSE",
        "Runtime (s)"
      )
    )
  )

write.csv(
  fit_table_csv,
  file.path(comparison_dir, paste0(prefix, "_fit_characteristics.csv")),
  row.names = FALSE
)

write_fit_latex_table(
  fit_table_csv,
  file.path(comparison_dir, paste0(prefix, "_fit_characteristics.tex"))
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
