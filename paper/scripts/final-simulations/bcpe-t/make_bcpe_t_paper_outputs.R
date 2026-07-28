#!/usr/bin/env Rscript

suppressPackageStartupMessages({
  library(dplyr)
  library(ggplot2)
  library(knitr)
  library(tidyr)
})

root <- normalizePath(".", winslash = "/", mustWork = TRUE)
result_root <- file.path(root, "results", "jss-exploratory", "01-continuous-bcpe-t")
if (!dir.exists(result_root)) {
  result_root <- file.path(root, "results")
}
rs_joint_dir <- Sys.getenv(
  "BCPE_T_PAPER_RS_JOINT_DIR",
  unset = file.path(result_root, "bcpe_t_current_defaults_rep100_rs_joint")
)
comparison_dir <- Sys.getenv(
  "BCPE_T_PAPER_COMPARISON_DIR",
  unset = file.path(result_root, "bcpe_t_current_defaults_rep100_comparison")
)
prefix <- "paper_simulation_bcpe_t"

model_labels <- c(
  "gamlss.longitudinal" = "gamlss.longitudinal",
  "gamlss2" = "gamlss2"
)

method_colours <- c(
  "gamlss.longitudinal" = "#0072B2",
  "gamlss2" = "#D55E00"
)

method_fills <- c(
  "gamlss.longitudinal" = "#56B4E9",
  "gamlss2" = "#E69F00"
)

param_labels <- c(
  "mu" = "mu",
  "sigma" = "sigma",
  "nu" = "nu",
  "tau" = "tau",
  "theta" = "theta",
  "zeta" = "zeta"
)

parameter_order <- names(param_labels)
marginal_parameters <- c("mu", "sigma", "nu", "tau")
marginal_smooth_parameters <- c("mu", "sigma")

read_result <- function(name) {
  path <- file.path(rs_joint_dir, name)
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
  ifelse(is.na(x), "--", formatC(x, digits = digits, format = "f"))
}

format_signed <- function(x, digits = 3) {
  ifelse(is.na(x), "--", sprintf(paste0("%+.", digits, "f"), x))
}

format_mean_sd <- function(mean, sd, digits = 3) {
  ifelse(
    is.na(mean),
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
    "nu" = "$\\nu$",
    "tau" = "$\\tau$",
    "theta" = "$\\theta$",
    "zeta" = "$\\zeta$",
    .default = latex_escape(x)
  )
}

term_latex <- function(x) {
  recode(
    x,
    "intercept" = "$\\beta_0$",
    "t" = "$t$",
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
    "\\caption{\\label{tab:parameter-recovery}Parameter recovery for the gamlss2 reference and gamlss.longitudinal model. Smooth rows report IRMSE in the RMSE column.}",
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

write_fit_latex_table <- function(data, path, note = NULL) {
  lines <- c(
    "\\begin{table}",
    "\\centering",
    "\\caption{\\label{tab:fit-characteristics}Fit and predictive characteristics for the gamlss2 reference and gamlss.longitudinal model. Values are mean (SD) across replicates. Marginal fixed-effect RMSE excludes copula parameters $\\theta$ and $\\zeta$; marginal smooth IRMSE uses only $\\mu$ and $\\sigma$ smooths.}",
    "\\begin{tabular}[t]{lrr}",
    "\\toprule",
    "Metric & gamlss2 & gamlss.longitudinal\\\\",
    "\\midrule"
  )

  body <- apply(data, 1, function(row) {
    paste(latex_escape(row[["Metric"]]), row[["gamlss2"]], row[["gamlss.longitudinal"]], sep = " & ") |>
      paste0("\\\\")
  })

  lines <- c(lines, body, "\\bottomrule", "\\end{tabular}")
  if (!is.null(note) && nzchar(note)) {
    lines <- c(lines, "\\vspace{2pt}", paste0("\\begin{flushleft}\\footnotesize ", note, "\\end{flushleft}"))
  }
  lines <- c(lines, "\\end{table}")
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

write_paper_latex_files <- function(recovery_table_latex, fit_table_csv, comparison_dir, prefix, fit_note = NULL) {
  recovery_path <- file.path(comparison_dir, paste0(prefix, "_fixed_parameter_bias_rmse.tex"))
  fit_path <- file.path(comparison_dir, paste0(prefix, "_fit_characteristics.tex"))
  write_recovery_latex_table(recovery_table_latex, recovery_path)
  write_fit_latex_table(fit_table_csv, fit_path, note = fit_note)
  invisible(c(recovery = recovery_path, fit = fit_path))
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
  filter(model %in% names(model_labels), parameter %in% parameter_order) %>%
  mutate(
    method = factor(as_method(model), levels = c("gamlss2", "gamlss.longitudinal")),
    parameter = factor(parameter, levels = parameter_order),
    term = factor(term, levels = c("intercept", "t", "x1", "x2")),
    term_math = recode(
      as.character(term),
      "intercept" = "beta[0]",
      "t" = "t",
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
  mutate(
    term_math = factor(term_math, levels = c("beta[0]", "t", "x[1]", "x[2]"))
  )

fixed_limits <- fixed_summary %>%
  group_by(parameter) %>%
  summarise(
    limit = max(abs(c(lower_error, upper_error)), na.rm = TRUE),
    .groups = "drop"
  ) %>%
  mutate(
    term_math = factor("beta[0]", levels = c("beta[0]", "t", "x[1]", "x[2]"))
  ) %>%
  tidyr::pivot_longer(
    cols = limit,
    names_to = "bound",
    values_to = "bias"
  ) %>%
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

fixed_plot <- ggplot(
  fixed_summary,
  aes(x = term_math, y = mean_error, colour = method)
) +
  geom_blank(data = fixed_limits, aes(x = term_math, y = bias), inherit.aes = FALSE) +
  geom_hline(yintercept = 0, linewidth = 0.35, colour = "grey30") +
  geom_pointrange(
    aes(ymin = lower_error, ymax = upper_error),
    position = position_dodge(width = 0.5),
    linewidth = 0.5,
    fatten = 2.5
  ) +
  facet_wrap(~ parameter, scales = "free_y", nrow = 1, labeller = label_parsed) +
  scale_colour_manual(values = method_colours) +
  scale_x_discrete(labels = function(x) parse(text = x)) +
  labs(
    x = NULL,
    y = "Bias"
  ) +
  theme_paper()

ggsave(
  file.path(comparison_dir, paste0(prefix, "_fixed_effect_recovery.png")),
  fixed_plot,
  width = 12,
  height = 3.8,
  dpi = 320,
  bg = "white"
)

smooth_parameters <- c("mu", "sigma", "theta")

smooth_by_rep <- read_result("smooth_estimates_by_rep.csv") %>%
  filter(model %in% names(model_labels), parameter %in% smooth_parameters) %>%
  mutate(
    method = factor(as_method(model), levels = c("gamlss2", "gamlss.longitudinal")),
    parameter = factor(parameter, levels = smooth_parameters)
  )

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
  geom_ribbon(
    aes(ymin = q25_hat, ymax = q75_hat, fill = method),
    colour = NA,
    alpha = 0.32
  ) +
  geom_line(
    data = truth_curves,
    aes(x = s1, y = true),
    inherit.aes = FALSE,
    colour = "black",
    linewidth = 0.75,
    linetype = "dashed"
  ) +
  geom_line(aes(y = mean_hat, colour = method), linewidth = 0.8) +
  facet_grid(
    method ~ parameter,
    scales = "free_y",
    drop = FALSE,
    labeller = labeller(parameter = label_parsed, method = label_value)
  ) +
  scale_colour_manual(values = method_colours, drop = FALSE) +
  scale_fill_manual(values = method_fills, drop = FALSE) +
  labs(
    x = NULL,
    y = "Smooth effect"
  ) +
  theme_paper() +
  theme(
    strip.text.y = element_blank(),
    strip.text.y.right = element_blank()
  )

ggsave(
  file.path(comparison_dir, paste0(prefix, "_smooth_effect_recovery.png")),
  smooth_plot,
  width = 11,
  height = 6.5,
  dpi = 320,
  bg = "white"
)
unlink(file.path(comparison_dir, paste0(prefix, "_gamlss_smooth_effect_recovery.png")))
unlink(file.path(comparison_dir, paste0(prefix, "_gamlss_longitudinal_smooth_effect_recovery.png")))

fixed_table_long <- read_result("fixed_effects_bias_rmse_table.csv") %>%
  filter(model %in% names(model_labels), parameter %in% parameter_order) %>%
  mutate(
    Method = as_method(model),
    Parameter = as.character(parameter),
    Term = term,
    Coverage = ifelse(is.na(rmse), NA_real_, coverage_95)
  ) %>%
  select(Parameter, Term, Method, Bias = bias, RMSE = rmse, Coverage)

fixed_wide <- fixed_table_long %>%
  pivot_wider(
    names_from = Method,
    values_from = c(Bias, RMSE, Coverage),
    names_glue = "{Method}_{.value}"
  ) %>%
  arrange(
    factor(Parameter, levels = parameter_order),
    desc(Term == "intercept"),
    Term
  )

smooth_table_long <- read_result("smooth_integrated_metrics.csv") %>%
  filter(model %in% names(model_labels), parameter %in% parameter_order) %>%
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
    names_glue = "{Method}_Smooth IRMSE"
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
    gamlss2_Smooth_IRMSE = `gamlss2_Smooth IRMSE`,
    `gamlss.longitudinal_Smooth_IRMSE` = `gamlss.longitudinal_Smooth IRMSE`
  )

fixed_wide <- fixed_wide %>%
  mutate(
    gamlss2_Smooth_IRMSE = NA_real_,
    `gamlss.longitudinal_Smooth_IRMSE` = NA_real_
  )

recovery_table_csv <- bind_rows(fixed_wide, smooth_wide) %>%
  arrange(
    factor(Parameter, levels = parameter_order),
    factor(Term, levels = c("intercept", "t", "x1", "x2", "smooth"))
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

recovery_table_latex <- recovery_table_csv %>%
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

summarise_metric <- function(data, metric, label, digits = 3) {
  data %>%
    group_by(Method) %>%
    summarise(
      mean = mean(.data[[metric]], na.rm = TRUE),
      sd = sd(.data[[metric]], na.rm = TRUE),
      .groups = "drop"
    ) %>%
    mutate(
      Metric = label,
      Value = format_mean_sd(mean, sd, digits)
    ) %>%
    select(Metric, Method, Value)
}

format_variogram_p <- function(p_value) {
  format(as.numeric(p_value), trim = TRUE, scientific = FALSE)
}

summarise_variogram_metrics <- function(data, digits = 3) {
  if (!"variogram_score" %in% names(data) || !any(is.finite(data$variogram_score))) {
    return(data.frame(
      Metric = "Variogram score (p=0.5)",
      Method = c("gamlss2", "gamlss.longitudinal"),
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

joint_by_rep <- read_result("joint_distribution_metrics_by_rep.csv") %>%
  filter(model %in% names(model_labels)) %>%
  mutate(Method = as_method(model))

predictive_by_rep <- read_result("predictive_scores_by_rep.csv") %>%
  filter(model %in% names(model_labels)) %>%
  mutate(
    Method = as_method(model),
    variogram_p_num = as.numeric(variogram_p)
  )

predictive_unique <- predictive_by_rep %>%
  distinct(Method, rep, .keep_all = TRUE)

predictive_for_variogram <- predictive_by_rep
fit_table_note <- NULL
trim_candidate <- predictive_by_rep %>%
  filter(Method == "gamlss2", abs(variogram_p_num - 2) < 1e-8, is.finite(variogram_score))
if (nrow(trim_candidate) > 0) {
  trim_q1 <- stats::quantile(trim_candidate$variogram_score, 0.25, na.rm = TRUE)
  trim_q3 <- stats::quantile(trim_candidate$variogram_score, 0.75, na.rm = TRUE)
  trim_iqr <- stats::IQR(trim_candidate$variogram_score, na.rm = TRUE)
  trim_cutoff <- as.numeric(trim_q3 + 3 * trim_iqr)
  trim_exclusions <- trim_candidate %>%
    filter(variogram_score > trim_cutoff) %>%
    arrange(desc(variogram_score)) %>%
    mutate(
      trim_rule = "Tukey outer fence: Q3 + 3 * IQR",
      trim_cutoff = trim_cutoff,
      raw_mean = mean(trim_candidate$variogram_score, na.rm = TRUE),
      raw_sd = stats::sd(trim_candidate$variogram_score, na.rm = TRUE)
    )
  write.csv(
    trim_exclusions,
    file.path(comparison_dir, paste0(prefix, "_gamlss2_variogram_p2_trim_exclusions.csv")),
    row.names = FALSE
  )
  if (nrow(trim_exclusions) > 0) {
    excluded_reps <- trim_exclusions$rep
    predictive_for_variogram <- predictive_by_rep %>%
      filter(!(Method == "gamlss2" & abs(variogram_p_num - 2) < 1e-8 & rep %in% excluded_reps))
    fit_table_note <- paste0(
      "* gamlss2 variogram score \\((p=2)\\) excludes ",
      nrow(trim_exclusions),
      " outer-fence replicate",
      ifelse(nrow(trim_exclusions) == 1, "", "s"),
      " (Tukey \\(Q_3 + 3\\,IQR\\); reps ",
      paste(excluded_reps, collapse = ", "),
      "). Raw value: ",
      format_mean_sd(mean(trim_candidate$variogram_score, na.rm = TRUE), stats::sd(trim_candidate$variogram_score, na.rm = TRUE), 3),
      "."
    )
  }
}

comparison_run_log <- file.path(comparison_dir, "current_defaults_rep100_run_logs.csv")
run_logs <- if (file.exists(comparison_run_log)) {
  read.csv(comparison_run_log, stringsAsFactors = FALSE) %>%
    filter(method %in% c("GAMLSS default", "RS joint")) %>%
    mutate(Method = recode(method, "GAMLSS default" = "gamlss2", "RS joint" = "gamlss.longitudinal"))
} else {
  read_result("fit_run_log.csv") %>%
    filter(model %in% names(model_labels)) %>%
    mutate(Method = as_method(model))
}

fixed_rmse_by_rep <- fixed_by_rep %>%
  filter(parameter %in% marginal_parameters) %>%
  group_by(Method = method, rep) %>%
  summarise(fixed_rmse = sqrt(mean(error^2, na.rm = TRUE)), .groups = "drop") %>%
  mutate(Method = as.character(Method))

smooth_irmse_by_rep <- read_result("smooth_integrated_metrics.csv") %>%
  filter(model %in% names(model_labels), parameter %in% marginal_smooth_parameters) %>%
  mutate(Method = as_method(model)) %>%
  group_by(Method, rep) %>%
  summarise(smooth_irmse = mean(irmse, na.rm = TRUE), .groups = "drop")

fit_table_long <- bind_rows(
  summarise_metric(joint_by_rep, "logLik", "Log-likelihood", 1),
  summarise_metric(predictive_unique, "test_log_score_per_obs", "Test log score / obs.", 3),
  summarise_variogram_metrics(predictive_for_variogram, 3),
  summarise_metric(joint_by_rep, "abs_rosenblatt_lag1_cor", "Abs. Rosenblatt lag-1 cor.", 3),
  summarise_metric(fixed_rmse_by_rep, "fixed_rmse", "Marginal fixed-effect RMSE", 3),
  summarise_metric(smooth_irmse_by_rep, "smooth_irmse", "Marginal smooth IRMSE", 3),
  summarise_metric(run_logs, "elapsed_sec", "Runtime (s)", 1)
)

if (!is.null(fit_table_note)) {
  fit_table_long <- fit_table_long %>%
    mutate(Value = ifelse(
      Method == "gamlss2" & Metric == "Variogram score (p=2)",
      paste0(Value, "*"),
      Value
    ))
}

variogram_metric_order <- unique(fit_table_long$Metric[grepl("^Variogram score", fit_table_long$Metric)])

fit_table_csv <- fit_table_long %>%
  mutate(Method = factor(Method, levels = c("gamlss2", "gamlss.longitudinal"))) %>%
  pivot_wider(names_from = Method, values_from = Value) %>%
  select(Metric, gamlss2, `gamlss.longitudinal`) %>%
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

write_paper_latex_files(recovery_table_latex, fit_table_csv, comparison_dir, prefix, fit_note = fit_table_note)

message("Wrote reduced paper outputs to: ", comparison_dir)
