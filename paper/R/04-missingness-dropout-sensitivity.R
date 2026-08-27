jss_run_04_missingness_dropout <- function(settings) {
  if (identical(settings$profile, "paper")) return(jss_run_04_missingness_from_public_inputs(settings))
  out <- file.path(settings$out_dir, "missingness")
  reps <- if (identical(settings$profile, "full")) "20" else if (identical(settings$profile, "paper")) "20" else "1"
  levels <- if (identical(settings$profile, "smoke")) "0,0.2" else "0,0.1,0.2,0.3,0.4,0.5"
  cores <- if (settings$workers > 1L) settings$workers else max(1L, parallel::detectCores() - 2L)
  jss_run_script(
    file.path(settings$root, "paper", "scripts", "final-simulations", "missingness", "run_missingness_study.R"),
    c(
      OUT_DIR = out, N_FITS = reps, N_CORES = as.character(cores),
      MISSING_LEVELS = levels, COMPUTE_SE = if (identical(settings$profile, "smoke")) "0" else "1",
      COMPUTE_PREDICTIVE_SCORES = "0", MAX_OUTER_ITER = if (identical(settings$profile, "smoke")) "3" else "1000",
      MAX_INNER_ITER = if (identical(settings$profile, "smoke")) "3" else "100",
      MAX_ELAPSED_SEC = if (identical(settings$profile, "smoke")) "20" else "180"
    ),
    settings$root
  )
  manuscript_figures <- file.path(
    settings$figures_dir,
    c("fixed_margin_rmse_by_missingness.png", "smooth_selected_recovery_curves.png")
  )
  source_figures <- file.path(out, basename(manuscript_figures))
  copied <- file.copy(source_figures, manuscript_figures, overwrite = TRUE)
  if (!all(copied)) {
    stop(
      "Missingness study completed but manuscript figures could not be staged: ",
      paste(source_figures[!copied], collapse = ", "),
      call. = FALSE
    )
  }
  list(
    module_id = "04-missingness-dropout-sensitivity", status = "regenerated",
    data = list.files(out, pattern = "[.]csv$", full.names = TRUE), tables = character(),
    figures = manuscript_figures
  )
}

jss_run_04_missingness_from_public_inputs <- function(settings) {
  input <- file.path(settings$public_data_dir, "missingness")
  fixed <- utils::read.csv(file.path(input, "fixed_term_summary_by_missingness.csv"), stringsAsFactors = FALSE)
  smooth <- utils::read.csv(file.path(input, "smooth_selected_plot_data.csv"), stringsAsFactors = FALSE)
  fixed$model <- factor(fixed$model, levels = c("gamlss2", "gamlss.longitudinal"))
  p1 <- ggplot2::ggplot(fixed[fixed$parameter %in% c("mu", "sigma", "nu", "tau"), ],
    ggplot2::aes(x = target_missing_pct, y = mean_rmse, colour = model, shape = model, group = model)) +
    ggplot2::geom_line(linewidth = 0.7) + ggplot2::geom_point(size = 1.9) +
    ggplot2::facet_grid(parameter ~ missing_mechanism, scales = "free_y") +
    ggplot2::scale_x_continuous(breaks = seq(0, 50, 10)) +
    ggplot2::labs(x = "Target missingness (%)", y = "Fixed-effect RMSE", colour = NULL, shape = NULL) +
    ggplot2::theme_bw(base_size = 10) + ggplot2::theme(legend.position = "top")
  f1 <- file.path(settings$figures_dir, "fixed_margin_rmse_by_missingness.png")
  ggplot2::ggsave(f1, p1, width = 8.5, height = 5.5, dpi = 180, bg = "white")

  smooth$missing_label <- factor(smooth$missing_label, levels = c("0%", "30%", "50%"))
  smooth$model <- factor(smooth$model, levels = c("gamlss2", "gamlss.longitudinal"))
  p2 <- ggplot2::ggplot(smooth, ggplot2::aes(x = s1)) +
    ggplot2::geom_line(ggplot2::aes(y = smooth_true), colour = "black", linewidth = 0.75) +
    ggplot2::geom_line(ggplot2::aes(y = smooth_median, colour = model), linewidth = 0.75, linetype = "dashed") +
    ggplot2::facet_grid(parameter + missing_mechanism ~ missing_label, scales = "free_y") +
    ggplot2::labs(x = "s1", y = "Centered smooth", colour = NULL) +
    ggplot2::theme_bw(base_size = 10) + ggplot2::theme(legend.position = "top")
  f2 <- file.path(settings$figures_dir, "smooth_selected_recovery_curves.png")
  ggplot2::ggsave(f2, p2, width = 8.5, height = 7, dpi = 180, bg = "white")
  list(module_id = "04-missingness-dropout-sensitivity", status = "regenerated",
    data = c(file.path(input, "fixed_term_summary_by_missingness.csv"), file.path(input, "smooth_selected_plot_data.csv")),
    tables = character(), figures = c(f1, f2))
}
