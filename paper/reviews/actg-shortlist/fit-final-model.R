#!/usr/bin/env Rscript

# Refit the selected ACTG model from the package-shipped full ACTG dataset.
#
# The established review object is not overwritten. This script writes a fresh
# fit to final-model-refit.rds in this directory unless --output= is supplied.

find_script_path <- function() {
  script_argument <- grep(
    "^--file=", commandArgs(trailingOnly = FALSE), value = TRUE
  )
  if (length(script_argument)) {
    return(normalizePath(
      sub("^--file=", "", script_argument[[1L]]),
      mustWork = TRUE
    ))
  }

  # `source()` records the sourced filename in one of the active frames.
  source_paths <- vapply(sys.frames(), function(frame) {
    if (is.null(frame$ofile)) NA_character_ else as.character(frame$ofile)
  }, character(1L))
  source_paths <- source_paths[!is.na(source_paths) & nzchar(source_paths)]
  source_paths <- source_paths[file.exists(source_paths)]
  if (length(source_paths)) {
    return(normalizePath(tail(source_paths, 1L), mustWork = TRUE))
  }

  # Support the RStudio Source button when its API is available.
  if (requireNamespace("rstudioapi", quietly = TRUE) && rstudioapi::isAvailable()) {
    active_path <- rstudioapi::getActiveDocumentContext()$path
    if (nzchar(active_path) && file.exists(active_path)) {
      return(normalizePath(active_path, mustWork = TRUE))
    }
  }

  # Last fallback for a console launched from the package root.
  candidate <- file.path(
    getwd(), "paper", "reviews", "actg-shortlist", "fit-final-model.R"
  )
  if (file.exists(candidate)) return(normalizePath(candidate, mustWork = TRUE))

  stop(
    "Could not determine the location of fit-final-model.R. ",
    "Run it with source('C:/full/path/fit-final-model.R') or Rscript."
  )
}

script_path <- find_script_path()
review_dir <- dirname(script_path)
package_root <- normalizePath(file.path(review_dir, "..", "..", ".."))

argument_value <- function(name, default) {
  prefix <- paste0("--", name, "=")
  hit <- grep(paste0("^", prefix), commandArgs(trailingOnly = TRUE), value = TRUE)
  if (!length(hit)) return(default)
  sub(prefix, "", hit[[1L]], fixed = TRUE)
}
output_path <- normalizePath(
  argument_value("output", file.path(review_dir, "final-model-refit.rds")),
  mustWork = FALSE
)

if (!requireNamespace("devtools", quietly = TRUE)) {
  stop("Install devtools with install.packages('devtools').")
}
devtools::load_all(package_root, quiet = TRUE)

load_package_data <- function(name) {
  environment <- new.env(parent = emptyenv())
  utils::data(
    list = name,
    package = "gamlss.longitudinal",
    envir = environment
  )
  get(name, envir = environment, inherits = FALSE)
}

actg_data <- load_package_data("actg_cd4")
actg_data$week_factor <- factor(actg_data$week)
actg_data$visit1 <- as.integer(actg_data$week == min(actg_data$week))

optimizer_control <- gamlss_longitudinal_control(
  outer_tol = 1e-4,
  max_outer_iter = 50L,
  max_elapsed_sec = 7200,
  rs = list(
    inner_tol = 1e-5,
    max_inner_iter = 100L,
    warm_start_joint = TRUE,
    warm_start_joint_iter = 5L,
    discrete_score_method = "analytical"
  )
)

plot_dist(actg_data,
  margin_dist = gamlss.dist::ZISICHEL(),
  time_var = "week_factor",
  subject_var = "subject",
  response_var= "cd4"
)
# save plot
ggsave(filename = file.path(review_dir, "final-model-plot_dist.png"))

missing = check_missingness(
  data = actg_data,
  time_var = "week_factor",
  subject_var = "subject",
  response_var = "cd4"
)

margin_screen = select_margin(
  data = actg_data,
  time_var = "week_factor",
  subject_var = "subject",
  response_var = "cd4"
)

margin_screen_time = select_margin(
  data = actg_data,
  time_var = "week_factor",
  subject_var = "subject",
  response_var = "cd4",
  time_intercepts = TRUE
)

#Margin selection>
#Criterion: AIC
#Selected: Zero Adjusted BNB (ZABNB)
#
# rank   family      AIC delta_AIC    screen_model converged pooled_AIC
#    1    ZABNB 29455.07     0.000 time_intercepts      TRUE   29478.27
#    2    ZIBNB 29455.07     0.000 time_intercepts      TRUE   29478.27
#    3    ZIPIG 29470.86    15.786 time_intercepts      TRUE   29496.88
#    4    ZAPIG 29470.86    15.786 time_intercepts      TRUE   29496.88
#    5 ZASICHEL 29474.56    19.486 time_intercepts      TRUE   29497.62
#    6 ZISICHEL 29474.56    19.486 time_intercepts      TRUE   29497.62
#    7      BNB 29517.54    62.466 time_intercepts      TRUE   29544.64
#    8   WARING 29577.59   122.520 time_intercepts      TRUE   29606.30
#    9      GPO 29589.97   134.892 time_intercepts      TRUE   29619.14
#   10   SICHEL 29665.88   210.810 time_intercepts      TRUE   29695.94

# ZIBNB, ZIPIG, ZISICHEL are top candidates based on AIC

copula_screen = select_copula(
  data = actg_data,
  time_var = "week_factor",
  subject_var = "subject",
  response_var = "cd4",
  margin_dist = ZIBNB() #best_fit_family(margin_screen_time)
)
#<Copula selection>
#Criterion: AIC | Complete pairs: 1849
#Selected: Student t (t)

# rank family        theta df tau   logLik  AIC       BIC       delta  
# 1    Student t (t) 0.804  4 0.594 945.548 -1887.096 -1876.051   0.000
# 2    Gumbel (G)    2.397 NA 0.583 936.299 -1870.599 -1865.076  16.498
# 3    Gaussian (N)  0.792 NA 0.582 882.022 -1762.043 -1756.521 125.053
# 4    Frank (F)     7.821 NA 0.596 868.196 -1734.393 -1728.870 152.704
# 5    Joe (J)       2.940 NA 0.511 820.720 -1639.440 -1633.917 247.657
# 6    Clayton (C)   1.862 NA 0.482 646.485 -1290.970 -1285.448 596.126
#Tail diagnostics are available in the diagnostic attributes.

# Student t and Gumbel are the top candidates based on AIC.

screen_joint = select_joint_distribution(
  data = actg_data,
  time_var = "week_factor",
  subject_var = "subject",
  response_var = "cd4",
  margin_families=c("ZIBNB", "ZIPIG", "ZISICHEL"),
  copula_families=c("t", "G")
)

plot_margin_fit(
  x = NULL,
  margin_dist = NULL,
  data = NULL,
  fit = NULL,
  response_var = "response",
  bins = 30,
  grid_n = 200,
  response_scale = c("response", "log"),
  time_intercepts = FALSE,
  by_time = FALSE,
  time_var = "time",
  fit_control = gamlss::gamlss.control(n.cyc = 50),
  plot = TRUE,
  ...
)

# Plot fitted marginal screen distributions against

zisichel_plot=plot_margin_fit(
  data = actg_data,
  time_var = "week_factor",
  response_var = "cd4",
  margin_dist = ZISICHEL(),
  time_intercepts=TRUE,
  by_time=TRUE,
  bins=100,
  separate_zero=TRUE
)
zisichel_plot$plot + ggplot2::coord_cartesian(xlim = c(0, 150)) +
  ggplot2::facet_wrap(~split_group, nrow = 1)
#save the plot we just showed with the ggplot adjustments
ggsave(file.path(review_dir, "zisichel_plot_adjusted.png"), plot = zisichel_plot$plot + ggplot2::coord_cartesian(xlim = c(0, 150)) +
  ggplot2::facet_wrap(~split_group, nrow = 1))
nbi_plot=plot_margin_fit(
  data = actg_data,
  time_var = "week_factor",
  response_var = "cd4",
  margin_dist = NBI(),
  time_intercepts=TRUE,
  by_time=TRUE,
  bins=100,
  separate_zero=TRUE
)
nbi_plot$plot + ggplot2::coord_cartesian(xlim = c(0, 150)) + 
ggplot2::facet_wrap(~split_group, nrow = 1)
#save the plot we just showed with the ggplot adjustments
ggsave(file.path(review_dir, "nbi_plot_adjusted.png"), plot = nbi_plot$plot + ggplot2::coord_cartesian(xlim = c(0, 150)) +
  ggplot2::facet_wrap(~split_group, nrow = 1))

zibnb_plot=plot_margin_fit(
  data = actg_data,
  time_var = "week_factor",
  response_var = "cd4",
  margin_dist = ZIBNB(),
  time_intercepts=TRUE,
  by_time=TRUE,
  bins=100,
  separate_zero=TRUE,
  fit_control = gamlss::gamlss.control(n.cyc = 1000),
)
zibnb_plot$plot + ggplot2::coord_cartesian(xlim = c(0, 150)) + 
ggplot2::facet_wrap(~split_group, nrow = 1)
#save the plot we just showed with the ggplot adjustments
ggsave(file.path(review_dir, "zibnb_plot_adjusted.png"), plot = zibnb_plot$plot + ggplot2::coord_cartesian(xlim = c(0, 150)) +
  ggplot2::facet_wrap(~split_group, nrow = 1))


# Overlay fitted gumbel copula contours against the observed data

gumbel_copula_fit=plot_copula_fit(
  data = actg_data,
  time_var = "week_factor",
  response_var = "cd4",
  margin_dist = ZISICHEL(),
  copula_dist = "G",
  by_time = TRUE
)
#save plot
ggsave(file.path(review_dir, "gumbel_copula_fit.png"), plot = gumbel_copula_fit$plot)

t_copula_fit=plot_copula_fit(
  data = actg_data,
  time_var = "week_factor",
  response_var = "cd4",
  margin_dist = ZISICHEL(),
  copula_dist = "t",
  by_time = TRUE
)
#save plot
ggsave(file.path(review_dir, "t_copula_fit.png"), plot = t_copula_fit$plot)

contours <- gamlss.longitudinal:::plot.copula_contour_compare(
  final_model,
  transform = "uniform",
  time_stratified = TRUE,
  contour_bins = 10,
  plot = FALSE
)

empirical_plot <- contours$plots$empirical_contours +
  ggplot2::facet_wrap(~time_pair, nrow = 1)

empirical_plot


# Overlay fitted marginal distributions on the observed data as well as poisson and NBI to compare
final_model <- gamlss_longitudinal(
  data = actg_data,
  time_var = "week_factor",
  subject_var = "subject",
  margin_dist = gamlss.dist::ZISICHEL(),
  copula_dist = "G", # Gumbel
  mu.formula = cd4 ~ treatment + age +
    s(baseline_cd4, bs = "ps", k = 8) +
    week_factor + baseline_cd4:visit1,
  sigma.formula = ~ treatment + baseline_cd4,
  nu.formula = ~ visit1,
  tau.formula = ~ treatment + sex + baseline_cd4,
  theta.formula = ~ treatment +
    s(baseline_cd4, bs = "ps", k = 8) +
    week_factor,
  zeta.formula = ~ 1,
  missingness = "segment",
  method = "RS",
  optimization = "joint",
  optimizer_control = optimizer_control,
  verbose = 3
)

saveRDS(final_model, output_path, compress = "gzip")
print(summary(final_model, include_vcov = FALSE))
cat("\nSaved refitted model to:\n", output_path, "\n", sep = "")

# Add option to load saved RDS above
if (file.exists(output_path)) {
  final_model <- readRDS(output_path)
}

final_model_terms <- plot_terms(final_model)
final_model_plot <- plot(final_model)
final_model_copula_diagnostics <- plot_copula_diagnostics(final_model, type="all",ncol=3)
#save plot
# plot()/plot_copula_diagnostics() dashboards are ggarrange grobs, but plot_terms()$dashboard
# is only layout metadata (plotlist/ncol/nrow), so it must be rearranged with ggarrange to save
ggsave(
  file.path(review_dir, "final_model_terms.png"),
  plot = ggpubr::ggarrange(
    plotlist = final_model_terms$dashboard$plotlist,
    ncol = final_model_terms$dashboard$ncol,
    nrow = final_model_terms$dashboard$nrow
  )
)
ggsave(file.path(review_dir, "final_model_plot.png"), plot = final_model_plot$dashboard)
ggsave(file.path(review_dir, "final_model_copula_diagnostics.png"), plot = final_model_copula_diagnostics$dashboard)

plot_margin_trt=plot_margin_fit(
  fit = final_model,
  data = actg_data,
  by = "treatment",
  bins = 100
)
# Use common axis limits and display treatment facets in one row. The package
# preserves the treatment factor levels for both marginal and copula plots.
plot_margin_trt_2=plot_margin_trt$plot +
  ggplot2::coord_cartesian(xlim = c(0, 150), ylim = c(0, 0.3)) +
  ggplot2::facet_wrap(~split_group, nrow = 1)
# save plot
ggsave(file.path(review_dir, "plot_margin_trt.png"), plot = plot_margin_trt_2)

plot_copula_trt=plot_copula_fit(
  fit = final_model,
  data = actg_data,
  by = "treatment"
)
plot_copula_trt_fin = plot_copula_trt$plot + 
ggplot2::facet_wrap(~split_group, nrow = 1)
ggsave(file.path(review_dir, "plot_copula_trt.png"), plot = plot_copula_trt_fin)


check_model(final_model)

publication_table(final_model)

# The subject-clustered sandwich inference used in the review can be reproduced
# after fitting with:
#
# final_inference <- summary(
#   final_model,
#   include_vcov = TRUE,
#   vcov_method = "sandwich"
# )
# saveRDS(final_inference, file.path(review_dir, "final-model-sandwich-summary.rds"))

check_model(final_model)
