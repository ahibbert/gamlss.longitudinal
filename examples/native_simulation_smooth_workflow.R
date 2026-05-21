# Native simulation, RS fitting, and diagnostics workflow
#
# This script is intended as the package's main end-to-end example. It uses the
# native copula simulation path, includes factor effects and interactions in the
# simulated truth, fits a joint RS model, and prepares a few diagnostic
# summaries that can be reused in the README or a vignette.

suppressPackageStartupMessages({
  library(gamlss.longitudinal)
  library(gamlss.dist)
  library(ggplot2)
})

set.seed(20260521)

# Keep the example quick enough for interactive use. Increase n_subject for a
# smoother recovery demonstration when you do not need term uncertainty plots.
n_subject <- 300
time_points <- seq(0, 1, length.out = 9)

treatment_effect <- c(
  control = 0,
  low = 0.20,
  high = 0.45
)

truth_mu_age_smooth <- function(x) {
  sim_smooth_bump(x, amplitude = 0.45, location = 0.62, width = 0.18) +
    sim_smooth_sin(x, amplitude = 0.12)
}

truth_sigma_age_smooth <- function(x) {
  sim_smooth_u(x, amplitude = 1.00)
}

truth_theta_age_smooth <- function(x) {
  sim_smooth_bump(x, amplitude = 0.55, location = 0.45, width = 0.22)
}

dat <- simulate_longitudinal_dataset(
  n = n_subject,
  times = time_points,
  margin_dist = GA(mu.link = "log", sigma.link = "log"),
  copula_dist = "C",
  covariates = function(base) {
    simulate_longitudinal_covariates(
      base,
      subject = list(
        treatment = function(subject_data) {
          factor(
            sample(c("control", "low", "high"), nrow(subject_data), replace = TRUE),
            levels = c("control", "low", "high")
          )
        },
        age = function(subject_data) stats::runif(nrow(subject_data), 25, 70)
      ),
      observation = list(
        time_scaled = function(long_data) sim_rescale01(long_data$time)
      )
    )
  },
  margin_params = list(
    mu = function(d) {
      age_scaled <- sim_rescale01(d$age)
      eta_mu <- 1.00 +
        sim_factor_effect(d$treatment, treatment_effect) +
        0.20 * d$time_scaled +
        truth_mu_age_smooth(age_scaled)
      exp(eta_mu)
    },
    sigma = function(d) {
      age_scaled <- sim_rescale01(d$age)
      eta_sigma <- -0.75 +
        0.12 * sim_factor_effect(d$treatment, c(control = 0, low = 1, high = 1)) +
        0.10 * d$time_scaled +
        truth_sigma_age_smooth(age_scaled)
      exp(eta_sigma)
    }
  ),
  copula_params = list(
    theta = function(e) {
      age_scaled <- sim_rescale01(e$age)
      time_left_scaled <- sim_rescale01(e$time_left)
      eta_theta <- 0.05 +
        0.20 * time_left_scaled +
        truth_theta_age_smooth(age_scaled)
      exp(eta_theta)
    }
  ),
  seed = 20260521
)

dat$age_scaled <- sim_rescale01(dat$age)

run_exploratory_plots <- isTRUE(getOption("gamlss.longitudinal.example.plots", interactive()))
run_model_diagnostics <- isTRUE(getOption("gamlss.longitudinal.example.model_diagnostics", FALSE))
run_marginal_family_screen <- isTRUE(getOption("gamlss.longitudinal.example.fitDist", FALSE))

example_dashboard <- function(..., plotlist = NULL, ncol = 1, nrow = NULL) {
  plots <- if (is.null(plotlist)) list(...) else plotlist
  if (is.null(nrow)) {
    nrow <- ceiling(length(plots) / ncol)
  }
  list(plotlist = plots, ncol = ncol, nrow = nrow)
}

print_example_dashboard <- function(x) {
  grid::grid.newpage()
  grid::pushViewport(grid::viewport(layout = grid::grid.layout(x$nrow, x$ncol)))
  on.exit(grid::popViewport(), add = TRUE)
  for (i_plot in seq_along(x$plotlist)) {
    row_id <- ((i_plot - 1) %/% x$ncol) + 1
    col_id <- ((i_plot - 1) %% x$ncol) + 1
    print(
      x$plotlist[[i_plot]],
      vp = grid::viewport(layout.pos.row = row_id, layout.pos.col = col_id)
    )
  }
  invisible(x)
}

sample_subjects <- sample(
  levels(dat$subject),
  size = min(40, length(levels(dat$subject)))
)
raw_data_plot <- ggplot(dat[dat$subject %in% sample_subjects, ], aes(
  x = time,
  y = response,
  group = subject,
  color = treatment
)) +
  geom_line(alpha = 0.35, linewidth = 0.4) +
  geom_point(alpha = 0.45, size = 0.8) +
  labs(
    title = "Observed Longitudinal Responses",
    x = "Time",
    y = "Response"
  ) +
  theme_minimal()

distribution_diagnostics <- NULL
if (run_exploratory_plots) {
  distribution_diagnostics <- plotDist(
    dat,
    GA(mu.link = "log", sigma.link = "log"),
    offdiag_scale = "pseudo"
  )
}

marginal_family_screen <- NULL
if (run_marginal_family_screen) {
  invisible(utils::capture.output({
    marginal_family_screen <- suppressWarnings(suppressMessages(gamlss::fitDist(
      dat$response,
      type = "realplus",
      try.gamlss = FALSE,
      trace = FALSE
    )))
  }))
  print(marginal_family_screen)
} else {
  message(
    "Skipping gamlss::fitDist() in the runnable example. ",
    "Set options(gamlss.longitudinal.example.fitDist = TRUE) to run it."
  )
}

# In a real analysis, pseudo-observations for this step would usually come
# from a preliminary marginal fit. In this simulated example we can use the
# stored truth uniforms to demonstrate the native copula screening step.
copula_family_screen <- select_copula(
  data = dat,
  u_var = "u",
  families = c("N", "C", "F", "G", "J", "t"),
  t_df_grid = c(3, 4, 6, 10, 20)
)
print(copula_family_screen)

truth_by_time <- stats::aggregate(
  cbind(true_mu, true_sigma, true_theta) ~ time + time_scaled + treatment,
  data = dat,
  FUN = function(x) mean(x, na.rm = TRUE)
)

fit_rs_joint <- gamlss.longitudinal(
  dataset = dat,
  margin_dist = GA(mu.link = "log", sigma.link = "log"),
  copula_dist = "C",
  time_var = "time",
  subject_var = "subject",
  mu.formula = response ~ treatment + time_scaled + s(age_scaled, bs = "ps", k = 8),
  sigma.formula = ~ treatment + time_scaled + s(age_scaled, bs = "ps", k = 5),
  theta.formula = ~ time_scaled + s(age_scaled, bs = "ps", k = 8)
)

print(fit_rs_joint$calc_lik_out_end$log_lik)
print(coef(fit_rs_joint))

copula_summary <- copula_time_summary(fit_rs_joint)

copula_truth <- stats::aggregate(
  true_theta ~ time + time_scaled,
  data = dat[!is.na(dat$true_theta), ],
  FUN = mean
)

eta_fit <- calc_eta(
  fit_rs_joint$par,
  fit_rs_joint$model_matrix,
  fit_rs_joint$margin_dist,
  get_copula_dist(fit_rs_joint$copula_dist)$copula_link,
  fit_rs_joint$par_s
)$eta_inv

fitted_margin <- data.frame(
  subject = fit_rs_joint$response_subject,
  time = fit_rs_joint$response_margin,
  fitted_mu = eta_fit$mu,
  fitted_sigma = eta_fit$sigma
)
fitted_margin <- merge(
  fitted_margin,
  dat[, c("subject", "time", "time_scaled", "treatment")],
  by = c("subject", "time"),
  all.x = TRUE,
  sort = FALSE
)
fitted_by_time <- stats::aggregate(
  cbind(fitted_mu, fitted_sigma) ~ time + time_scaled + treatment,
  data = fitted_margin,
  FUN = function(x) mean(x, na.rm = TRUE)
)

fit_time_levels <- sort(unique(fit_rs_joint$response_margin))
theta_rows <- fit_rs_joint$response_margin %in% fit_time_levels[-length(fit_time_levels)]
fitted_theta <- data.frame(
  subject = fit_rs_joint$response_subject[theta_rows],
  time = fit_rs_joint$response_margin[theta_rows],
  fitted_theta = eta_fit$theta
)
fitted_theta <- merge(
  fitted_theta,
  dat[, c("subject", "time", "time_scaled")],
  by = c("subject", "time"),
  all.x = TRUE,
  sort = FALSE
)
fitted_theta_by_time <- stats::aggregate(
  fitted_theta ~ time + time_scaled,
  data = fitted_theta,
  FUN = function(x) mean(x, na.rm = TRUE)
)

truth_mu_plot <- ggplot(truth_by_time, aes(x = time_scaled, y = true_mu, color = treatment)) +
  geom_line(linewidth = 0.8) +
  geom_point(size = 1.8) +
  labs(
    title = "Truth: mu",
    x = "Scaled time",
    y = "mu"
  ) +
  theme_minimal()

truth_sigma_plot <- ggplot(truth_by_time, aes(x = time_scaled, y = true_sigma, color = treatment)) +
  geom_line(linewidth = 0.8) +
  geom_point(size = 1.8) +
  labs(
    title = "Truth: sigma",
    x = "Scaled time",
    y = "sigma"
  ) +
  theme_minimal()

truth_theta_plot <- ggplot(copula_truth, aes(x = time_scaled, y = true_theta)) +
  geom_line(linewidth = 0.8, color = "#B3262E") +
  geom_point(size = 1.8, color = "#B3262E") +
  labs(
    title = "Truth: Clayton theta",
    x = "Scaled time",
    y = "theta"
  ) +
  theme_minimal()

truth_dashboard <- example_dashboard(
  truth_mu_plot,
  truth_sigma_plot,
  truth_theta_plot,
  ncol = 1,
  nrow = 3
)

fitted_mu_plot <- ggplot(fitted_by_time, aes(x = time_scaled, y = fitted_mu, color = treatment)) +
  geom_line(linewidth = 0.8) +
  geom_point(size = 1.8) +
  labs(
    title = "Fitted: mu",
    x = "Scaled time",
    y = "mu"
  ) +
  theme_minimal()

fitted_sigma_plot <- ggplot(fitted_by_time, aes(x = time_scaled, y = fitted_sigma, color = treatment)) +
  geom_line(linewidth = 0.8) +
  geom_point(size = 1.8) +
  labs(
    title = "Fitted: sigma",
    x = "Scaled time",
    y = "sigma"
  ) +
  theme_minimal()

fitted_theta_plot <- ggplot(fitted_theta_by_time, aes(x = time_scaled, y = fitted_theta)) +
  geom_line(linewidth = 0.8, color = "#1F4E79") +
  geom_point(size = 1.8, color = "#1F4E79") +
  labs(
    title = "Fitted: Clayton theta",
    x = "Scaled time",
    y = "theta"
  ) +
  theme_minimal()

fitted_parameter_dashboard <- example_dashboard(
  fitted_mu_plot,
  fitted_sigma_plot,
  fitted_theta_plot,
  ncol = 1,
  nrow = 3
)

truth_term_theme <- function(p) {
  p +
    theme_minimal() +
    theme(
      plot.title = element_text(size = 14),
      axis.title = element_text(size = 11),
      plot.caption = element_text(size = 9)
    )
}

truth_smooth_plot <- function(x, contribution, title, xlab = "time_scaled") {
  plot_df <- data.frame(
    x = x,
    contribution = contribution
  )
  plot_df <- stats::aggregate(contribution ~ x, data = plot_df, FUN = mean)
  plot_df <- plot_df[order(plot_df$x), , drop = FALSE]

  truth_term_theme(
    ggplot(plot_df, aes(x = x, y = contribution)) +
      geom_line(color = "black", linewidth = 2) +
      labs(
        title = title,
        x = xlab,
        y = paste0("smooth(", xlab, ")"),
        caption = "truth"
      )
  )
}

truth_continuous_plot <- function(x, contribution, title, xlab, ylab) {
  plot_df <- data.frame(
    x = x,
    contribution = contribution
  )
  plot_df <- stats::aggregate(contribution ~ x, data = plot_df, FUN = mean)
  plot_df <- plot_df[order(plot_df$x), , drop = FALSE]

  truth_term_theme(
    ggplot(plot_df, aes(x = x, y = contribution)) +
      geom_line(color = "black", linewidth = 2) +
      labs(
        title = title,
        x = xlab,
        y = ylab,
        caption = "truth"
      )
  )
}

truth_factor_plot <- function(effects, title, xlab, ylab) {
  plot_df <- data.frame(
    x = seq_along(effects),
    level = factor(names(effects), levels = names(effects)),
    contribution = as.numeric(effects)
  )

  truth_term_theme(
    ggplot(plot_df, aes(x = x, y = contribution)) +
      geom_hline(yintercept = 0, color = "grey70", linetype = 3) +
      geom_point(color = "black", size = 1.2) +
      geom_errorbar(aes(ymin = contribution, ymax = contribution), color = "red", width = 0.15) +
      scale_x_continuous(breaks = plot_df$x, labels = plot_df$level) +
      labs(
        title = title,
        x = xlab,
        y = ylab,
        caption = "truth"
      )
  )
}

time_grid <- seq(0, 1, length.out = 200)
age_grid <- seq(0, 1, length.out = 200)

truth_term_plots <- list(
  truth_smooth_plot(
    x = age_grid,
    contribution = truth_mu_age_smooth(age_grid),
    title = 'mu: s(age_scaled, bs = "ps", k = 8)',
    xlab = "age_scaled"
  ),
  truth_smooth_plot(
    x = age_grid,
    contribution = truth_sigma_age_smooth(age_grid),
    title = 'sigma: s(age_scaled, bs = "ps", k = 5)',
    xlab = "age_scaled"
  ),
  truth_smooth_plot(
    x = age_grid,
    contribution = truth_theta_age_smooth(age_grid),
    title = 'theta: s(age_scaled, bs = "ps", k = 8)',
    xlab = "age_scaled"
  ),
  truth_factor_plot(
    effects = treatment_effect,
    title = "mu: treatment",
    xlab = "treatment",
    ylab = "fixed contribution: mu.treatment"
  ),
  truth_continuous_plot(
    x = time_grid,
    contribution = 0.20 * time_grid,
    title = "mu: time_scaled",
    xlab = "time_scaled",
    ylab = "fixed contribution: mu.time_scaled"
  ),
  truth_factor_plot(
    effects = c(control = 0, low = 0.12, high = 0.12),
    title = "sigma: treatment",
    xlab = "treatment",
    ylab = "fixed contribution: sigma.treatment"
  ),
  truth_continuous_plot(
    x = time_grid,
    contribution = 0.10 * time_grid,
    title = "sigma: time_scaled",
    xlab = "time_scaled",
    ylab = "fixed contribution: sigma.time_scaled"
  ),
  truth_continuous_plot(
    x = time_grid,
    contribution = 0.20 * time_grid,
    title = "theta: time_scaled",
    xlab = "time_scaled",
    ylab = "fixed contribution: theta.time_scaled"
  )
)

truth_terms_dashboard <- example_dashboard(
  plotlist = truth_term_plots,
  ncol = 4,
  nrow = ceiling(length(truth_term_plots) / 4)
)

copula_plot_data <- merge(
  copula_truth,
  copula_summary$time_summary[, c("time", "theta_fit", "tau_fit")],
  by = "time",
  all.x = TRUE,
  sort = FALSE
)

copula_plot <- ggplot(copula_plot_data, aes(x = time_scaled)) +
  geom_line(aes(y = true_theta), linewidth = 0.8, color = "#B3262E") +
  geom_point(aes(y = theta_fit), size = 1.8, color = "#1F4E79") +
  labs(
    title = "Clayton Dependence: Truth and Fitted Time Summary",
    x = "Scaled time",
    y = "Theta"
  ) +
  theme_minimal()

term_diagnostics <- NULL
model_diagnostics <- NULL
copula_diagnostics <- NULL

if (interactive()) {
  print(raw_data_plot)
  if (!is.null(distribution_diagnostics)) print(distribution_diagnostics)
  print_example_dashboard(truth_dashboard)
  print_example_dashboard(fitted_parameter_dashboard)
  print(copula_plot)
  print_example_dashboard(truth_terms_dashboard)
  term_diagnostics <- plot.terms(fit_rs_joint, data = dat)
  if (run_model_diagnostics) {
    model_diagnostics <- plot(fit_rs_joint, data = dat)
  }
  copula_diagnostics <- plot.copula(fit_rs_joint, data = dat)
} else {
  print_example_dashboard(truth_terms_dashboard)
  term_diagnostics <- plot.terms(fit_rs_joint, data = dat, paginate = TRUE)
  if (run_model_diagnostics) {
    model_diagnostics <- plot(fit_rs_joint, data = dat)
  }
  copula_diagnostics <- plot.copula(fit_rs_joint, data = dat, plot = FALSE)
}

example_objects <- list(
  data = dat,
  fit = fit_rs_joint,
  marginal_family_screen = marginal_family_screen,
  copula_family_screen = copula_family_screen,
  truth_by_time = truth_by_time,
  fitted_by_time = fitted_by_time,
  fitted_theta_by_time = fitted_theta_by_time,
  copula_summary = copula_summary,
  diagnostics = list(
    terms = term_diagnostics,
    model = model_diagnostics,
    copula = copula_diagnostics
  ),
  plots = list(
    truth_dashboard = truth_dashboard,
    truth_terms_dashboard = truth_terms_dashboard,
    truth_terms = truth_term_plots,
    raw_data = raw_data_plot,
    distribution = distribution_diagnostics,
    fitted_parameter_dashboard = fitted_parameter_dashboard,
    copula = copula_plot
  )
)
