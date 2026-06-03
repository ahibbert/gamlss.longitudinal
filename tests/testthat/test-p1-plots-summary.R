test_that("T151 summary output contract is stable", {
  dat <- make_fixture_factor_time_interaction(n_subject = 18L)
  fit <- fit_fixture_model(dat, include_dlcopdpar = TRUE)

  s <- summary(fit, include_vcov = FALSE)

  expect_s3_class(s, "summary.gamlss.longitudinal")
  expect_true(all(c("model", "fit", "coefficients", "smooth_terms") %in% names(s)))
  expect_true(is.data.frame(s$coefficients))
  expect_true(all(c("term", "estimate", "std_error", "p_value", "signif") %in% names(s$coefficients)))
  expect_true(is.finite(s$fit$logLik))
})

test_that("T151b summary print shows tiny p-values as less than threshold", {
  dat <- make_fixture_factor_time_interaction(n_subject = 18L)
  fit <- fit_fixture_model(dat, include_dlcopdpar = TRUE)

  s <- summary(fit, include_vcov = FALSE)
  s$coefficients$p_value[1] <- 1e-8
  s$coefficients$signif[1] <- "***"

  txt <- capture.output(print(s))
  expect_true(any(grepl("<0.00001", txt, fixed = TRUE)))
})

test_that("T152 plot_terms interaction rendering metadata includes factor levels", {
  dat <- make_fixture_factor_time_interaction(n_subject = 16L)
  dat$time_raw <- factor(as.character(dat$time_raw), levels = levels(dat$time_raw), ordered = FALSE)
  fit <- fit_fixture_model(dat, include_dlcopdpar = TRUE)

  suppressPackageStartupMessages(library(grid))

  pt <- suppressWarnings(plot_terms(fit, data = dat, plot_interactions = TRUE))
  pf <- pt$fixed_terms

  mu_entries <- pf$mu
  entry_names <- names(mu_entries)
  interaction_entries <- entry_names[grepl(":", entry_names, fixed = TRUE)]

  expect_equal(length(interaction_entries), 1)
  interaction_entry <- mu_entries[[interaction_entries[1]]]
  expect_setequal(interaction_entry$levels, levels(dat$time_raw))
  expect_setequal(interaction_entry$series, levels(dat$gender))
  expect_equal(nrow(interaction_entry$plot_data), length(levels(dat$time_raw)) * length(levels(dat$gender)))
})

test_that("T153 plot methods smoke test return dashboard structures", {
  dat <- make_fixture_factor_time_interaction(n_subject = 14L)
  fit <- fit_fixture_model(dat, include_dlcopdpar = TRUE)

  suppressPackageStartupMessages(library(grid))

  pt <- suppressWarnings(plot_terms(fit, data = dat))
  expect_true(is.list(pt))
  expect_true(all(c("smooth_terms", "fixed_terms", "dashboard") %in% names(pt)))

  pdiag <- suppressWarnings(plot(fit, data = dat))
  expect_true(is.list(pdiag))
  expect_true(all(c("diagnostics", "forecasts", "dashboard") %in% names(pdiag)))
})

test_that("T154 plot_terms handles no-data fixed-term plots with many time levels", {
  set.seed(42)

  subject_tbl <- data.frame(
    id = seq_len(18L),
    gender = factor(sample(c("F", "M"), 18L, replace = TRUE)),
    age = round(runif(18L, min = 20, max = 70), 1),
    stringsAsFactors = FALSE
  )
  time_levels <- paste0("t", seq_len(8L))
  grid <- expand.grid(
    id = subject_tbl$id,
    time_raw = factor(time_levels, levels = time_levels, ordered = TRUE),
    KEEP.OUT.ATTRS = FALSE,
    stringsAsFactors = FALSE
  )
  dat <- merge(grid, subject_tbl, by = "id", sort = FALSE)
  dat <- dat[order(dat$id, dat$time_raw), ]
  rownames(dat) <- NULL

  t_num <- as.integer(dat$time_raw)
  g_num <- ifelse(dat$gender == "M", 1, 0)
  dat$y <- 1.5 + 0.15 * t_num + 0.35 * g_num + 0.2 * t_num * g_num + 0.01 * dat$age + rnorm(nrow(dat), sd = 0.12)

  fit <- fit_fixture_model(
    dat,
    include_dlcopdpar = TRUE,
    mu_formula = "y ~ time_raw * gender + age",
    sigma_formula = "~ time_raw + gender",
    theta_formula = "~ time_raw",
    max_outer_iter = 2,
    max_inner_iter = 2,
    use_backtracking = TRUE
  )

  pt <- suppressWarnings(plot_terms(fit))
  expect_true(is.list(pt))
  expect_true(length(pt$fixed_terms) > 0)
  mu_entries <- pt$fixed_terms$mu
  expect_false(any(grepl("\\|", names(mu_entries), perl = TRUE)))
})

test_that("T155 ordered factor time is handled like nominal factor in grouped plots", {
  dat <- make_fixture_factor_time_interaction(n_subject = 16L)
  dat$time_raw <- factor(as.character(dat$time_raw), levels = levels(dat$time_raw), ordered = TRUE)
  fit <- fit_fixture_model(dat, include_dlcopdpar = TRUE)

  pt <- suppressWarnings(plot_terms(fit, data = dat, plot_interactions = TRUE))
  mu_entries <- pt$fixed_terms$mu
  interaction_entries <- names(mu_entries)[grepl(":", names(mu_entries), fixed = TRUE)]

  expect_equal(length(interaction_entries), 1)
  interaction_entry <- mu_entries[[interaction_entries[1]]]
  expect_setequal(interaction_entry$levels, levels(dat$time_raw))
  expect_setequal(interaction_entry$series, levels(dat$gender))
})

test_that("T156 transformed smooth covariates keep their x-axis scale", {
  dat <- make_fixture_factor_time_interaction(n_subject = 16L)
  fit <- fit_fixture_model(
    dat,
    include_dlcopdpar = TRUE,
    mu_formula = "y ~ time_raw * gender + s(log(age), bs = 'ps')"
  )

  pt <- suppressWarnings(plot_terms(fit))
  smooth_entries <- pt$smooth_terms$mu
  expect_equal(length(smooth_entries), 1)

  smooth_entry <- smooth_entries[[1]]
  expect_equal(smooth_entry$x, log(dat$age), tolerance = 1e-12)
  expect_equal(smooth_entry$x[1], log(dat$age)[1], tolerance = 1e-12)
  expect_s3_class(smooth_entry$plot, "ggplot")
})

test_that("T157 copula plot wrappers remain available after install", {
  dat <- make_fixture_factor_time_interaction(n_subject = 16L)
  fit <- fit_fixture_model(dat, include_dlcopdpar = TRUE)

  p1 <- suppressWarnings(plot_copula_diagnostics(fit, plot = FALSE))
  p2 <- suppressWarnings(plot.copula_contour_compare(fit, plot = FALSE))

  expect_true(is.list(p1))
  expect_true(all(c("plots", "dashboard", "fit_data", "pair_data", "quartile_summary") %in% names(p1)))
  expect_equal(length(p1$plots), 9L)
  expect_true("rosenblatt_qq" %in% names(p1$plots))
  expect_true(is.list(p2))
  expect_true(all(c("plots", "dashboard", "grid", "metrics") %in% names(p2)))
})

test_that("T158 old plotting entry points delegate with lifecycle warnings", {
  dat <- make_fixture_factor_time_interaction(n_subject = 14L)
  fit <- fit_fixture_model(dat, include_dlcopdpar = TRUE)

  expect_warning(
    old_terms <- plot.terms(fit, data = dat),
    "plot_terms"
  )
  expect_true(is.list(old_terms))

  old_copula <- suppressWarnings(plot.copula(fit, plot = FALSE))
  expect_true(is.list(old_copula))
  expect_true(all(c("plots", "dashboard", "fit_data", "pair_data", "quartile_summary") %in% names(old_copula)))
})

test_that("T159 standard fit inspection plotting helpers are available", {
  dat0 <- make_fixture_factor_time_interaction(n_subject = 12L)
  dat <- data.frame(
    subject = dat0$id,
    time = dat0$time_raw,
    response = dat0$y
  )

  margin_plot <- suppressWarnings(plot_margin_fit(
    dat,
    margin_dist = gamlss.dist::NO(),
    response_var = "response",
    plot = FALSE
  ))
  expect_s3_class(margin_plot$plot, "ggplot")
  expect_true(all(c("plot", "data", "density") %in% names(margin_plot)))

  copula_screen <- select_copula(
    data = dat,
    response_var = "response",
    margin_dist = gamlss.dist::NO(),
    subject_var = "subject",
    time_var = "time",
    families = c("N", "C"),
    min_pairs = 5
  )
  copula_plot <- suppressWarnings(plot_copula_fit(
    data = dat,
    copula_dist = copula_screen,
    response_var = "response",
    margin_dist = gamlss.dist::NO(),
    subject_var = "subject",
    time_var = "time",
    plot = FALSE
  ))
  expect_s3_class(copula_plot$plot, "ggplot")
  expect_true(all(c("plot", "pair_data", "density") %in% names(copula_plot)))

  grid_plot <- suppressWarnings(plot_dist(
    dat,
    margin_dist = gamlss.dist::NO(),
    subject_var = "subject",
    time_var = "time",
    response_var = "response",
    overlay = "margin"
  ))
  expect_false(is.null(grid_plot))

  fit <- fit_fixture_model(dat0, include_dlcopdpar = TRUE)
  raw_grid_plot <- suppressWarnings(plot_dist(
    dat,
    subject_var = "subject",
    time_var = "time",
    response_var = "response"
  ))
  model_margin_plot <- suppressWarnings(plot_margin_fit(data = dat0, fit = fit, plot = FALSE))
  model_copula_plot <- suppressWarnings(plot_copula_fit(data = dat0, fit = fit, plot = FALSE))
  model_copula_overlay_plot <- suppressWarnings(plot_copula_overlay(data = dat0, fit = fit, plot = FALSE))
  model_grid_plot <- suppressWarnings(plot_dist(
    dat,
    fit = fit
  ))

  expect_false(is.null(raw_grid_plot))
  expect_s3_class(model_margin_plot$plot, "ggplot")
  expect_s3_class(model_copula_plot$plot, "ggplot")
  expect_s3_class(model_copula_overlay_plot$plot, "ggplot")
  expect_false(is.null(model_grid_plot))

  expect_error(
    plot_margin_fit(dat, family = gamlss.dist::NO(), response_var = "response", plot = FALSE),
    "margin_dist"
  )
  expect_error(
    plot_copula_fit(
      data = dat,
      copula = copula_screen,
      response_var = "response",
      margin_dist = gamlss.dist::NO(),
      subject_var = "subject",
      time_var = "time",
      plot = FALSE
    ),
    "copula_dist"
  )
  expect_error(
    plot_dist(
      dat,
      dist = gamlss.dist::NO(),
      subject_var = "subject",
      time_var = "time",
      response_var = "response"
    ),
    "margin_dist"
  )
})
