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

test_that("T151c summary includes EDF for smooth terms", {
  dat <- make_fixture_factor_time_interaction(n_subject = 16L)
  fit <- fit_fixture_model(
    dat,
    include_dlcopdpar = TRUE,
    mu_formula = "y ~ time_raw * gender + s(log(age), bs = 'ps')"
  )

  s <- summary(fit, include_vcov = FALSE)

  expect_true(all(c("parameter", "smooth_term", "edf") %in% names(s$smooth_terms)))
  expect_equal(nrow(s$smooth_terms), 1L)
  expect_equal(s$smooth_terms$parameter, "mu")
  expect_equal(s$smooth_terms$smooth_term, "s(log(age), bs = \"ps\")")
  expect_true(is.finite(s$smooth_terms$edf))

  txt <- capture.output(print(s))
  expect_true(any(grepl("\\bedf\\b", txt)))
  expect_true(any(grepl("s\\(log\\(age\\), bs = \"ps\"\\)", txt)))
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
  expect_null(pdiag$forecasts$fitted_quantiles)
  expect_null(pdiag$forecasts$newdata_quantiles)
  expect_null(pdiag$fitted_data)
  expect_null(pdiag$newdata_data)

  pdiag_forecasts <- suppressWarnings(plot(
    fit,
    data = dat,
    include_fitted_quantiles = TRUE,
    include_newdata_quantiles = TRUE
  ))
  expect_s3_class(pdiag_forecasts$forecasts$fitted_quantiles, "ggplot")
  expect_s3_class(pdiag_forecasts$forecasts$newdata_quantiles, "ggplot")
  expect_true(is.data.frame(pdiag_forecasts$fitted_data))
  expect_true(is.data.frame(pdiag_forecasts$newdata_data))
})

test_that("T153b PIT histogram bins stay inside unit interval", {
  dat <- make_fixture_factor_time_interaction(n_subject = 14L)
  fit <- fit_fixture_model(dat, include_dlcopdpar = TRUE)

  pit_plot <- suppressWarnings(pithist(fit, bins = 20, plot = TRUE))
  hist_layer <- ggplot2::ggplot_build(pit_plot)$data[[1]]

  expect_true(all(hist_layer$xmin >= 0))
  expect_true(all(hist_layer$xmax <= 1))
  expect_equal(pit_plot$layers[[1]]$stat_params$breaks, seq(0, 1, length.out = 21L))
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

test_that("T156b smooth plotting tolerates unavailable interval values", {
  dat <- make_fixture_factor_time_interaction(n_subject = 16L)
  fit <- fit_fixture_model(
    dat,
    include_dlcopdpar = TRUE,
    mu_formula = "y ~ time_raw * gender + s(log(age), bs = 'ps')"
  )

  smooth_terms <- suppressWarnings(plot_smooth_terms(
    fit,
    vcov_obj = list(vcov = list(smooth_vcov = NULL, smooth_se = NULL)),
    data = dat,
    setup_mfrow = FALSE,
    show_legend = FALSE,
    grid_n = 20
  ))

  smooth_entry <- smooth_terms$mu[[1]]
  expect_s3_class(smooth_entry$plot, "ggplot")
  expect_true(all(is.na(smooth_entry$ci_lower)))
  expect_true(all(is.na(smooth_entry$ci_upper)))
})

test_that("T156c transformed factor time is grouped into point interval plots", {
  dat <- make_fixture_factor_time_interaction(n_subject = 16L)
  dat$time_raw <- factor(as.character(dat$time_raw), levels = levels(dat$time_raw), ordered = FALSE)
  fit <- fit_fixture_model(
    dat,
    include_dlcopdpar = TRUE,
    mu_formula = "y ~ as.factor(time_raw) + gender + age",
    sigma_formula = "~ as.factor(time_raw) + gender",
    theta_formula = "~ as.factor(time_raw)"
  )

  pt <- suppressWarnings(plot_terms(fit, data = dat))

  for (par_name in c("mu", "sigma", "theta")) {
    entries <- pt$fixed_terms[[par_name]]
    factor_entries <- names(entries)[grepl("as.factor\\(time_covariate\\)|as.factor\\(time_raw\\)", names(entries))]

    expect_equal(length(factor_entries), 1L)
    factor_entry <- entries[[factor_entries[1]]]
    expect_setequal(factor_entry$levels, levels(dat$time_raw))
    expect_equal(factor_entry$fitted[1], 0)
    expect_equal(factor_entry$se[1], 0)

    layer_classes <- vapply(
      factor_entry$plot$layers,
      function(layer) class(layer$geom)[1],
      character(1)
    )
    expect_true("GeomPoint" %in% layer_classes)
    expect_true("GeomErrorbar" %in% layer_classes)
    expect_false("GeomLine" %in% layer_classes)
    expect_false("GeomRibbon" %in% layer_classes)
  }
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

  positive_response <- gamlss.dist::rGG(
    80,
    mu = 3,
    sigma = 0.7,
    nu = 1.2
  )
  positive_margin_plot <- suppressWarnings(plot_margin_fit(
    positive_response,
    margin_dist = gamlss.dist::GG(mu.link = "log", sigma.link = "log", nu.link = "identity"),
    plot = FALSE
  ))
  expect_gt(min(positive_margin_plot$density$response, na.rm = TRUE), 0)
  expect_true(any(is.finite(positive_margin_plot$density$density)))

  set.seed(159)
  gamma_low_end_response <- c(
    0.02,
    gamlss.dist::rGA(79, mu = 4, sigma = 0.35)
  )
  gamma_low_end_plot <- suppressWarnings(plot_margin_fit(
    gamma_low_end_response,
    margin_dist = gamlss.dist::GA(mu.link = "log", sigma.link = "log"),
    plot = FALSE
  ))
  gamma_density <- gamma_low_end_plot$density
  finite_density <- gamma_density[is.finite(gamma_density$density), ]
  gamma_response <- gamma_low_end_plot$data$response
  gamma_response <- gamma_response[is.finite(gamma_response) & gamma_response > 0]
  old_lower_start <- max(
    as.numeric(stats::quantile(gamma_response, probs = 0.05, type = 8, names = FALSE)),
    max(gamma_response) * 0.02
  )
  expect_lt(min(finite_density$response), old_lower_start)
  expect_gt(min(finite_density$response), 0)

  set.seed(159)
  unstable_bcpe_response <- gamlss.dist::rGG(
    540,
    mu = 7.2,
    sigma = 0.96,
    nu = 1.63
  )
  unstable_bcpe_plot <- NULL
  expect_warning(
    unstable_bcpe_plot <- plot_margin_fit(
      unstable_bcpe_response,
      margin_dist = gamlss.dist::BCPE(),
      plot = FALSE
    ),
    "BCPE marginal overlay fit did not converge"
  )
  expect_s3_class(unstable_bcpe_plot$plot, "ggplot")
  expect_equal(nrow(unstable_bcpe_plot$density), 0L)

  unstable_bcpe_time_data <- data.frame(
    response = unstable_bcpe_response,
    time = rep(seq_len(3), length.out = length(unstable_bcpe_response))
  )
  unstable_bcpe_time_plot <- NULL
  expect_warning(
    unstable_bcpe_time_plot <- plot_margin_fit(
      unstable_bcpe_time_data,
      margin_dist = gamlss.dist::BCPE(),
      response_var = "response",
      time_var = "time",
      by_time = TRUE,
      plot = FALSE
    ),
    "BCPE marginal overlay fit did not converge"
  )
  expect_s3_class(unstable_bcpe_time_plot$plot, "ggplot")
  expect_equal(nrow(unstable_bcpe_time_plot$density), 0L)

  time_intercept_margin_plot <- suppressWarnings(plot_margin_fit(
    dat,
    margin_dist = gamlss.dist::NO(),
    response_var = "response",
    time_intercepts = TRUE,
    plot = FALSE
  ))
  time_intercept_facet_plot <- suppressWarnings(plot_margin_fit(
    dat,
    margin_dist = gamlss.dist::NO(),
    response_var = "response",
    time_intercepts = TRUE,
    by_time = TRUE,
    plot = FALSE
  ))
  expect_s3_class(time_intercept_margin_plot$plot, "ggplot")
  expect_equal(unique(time_intercept_margin_plot$density$split_group), "All")
  expect_setequal(unique(time_intercept_facet_plot$density$split_group), as.character(unique(dat$time)))
  expect_true(any(is.finite(time_intercept_facet_plot$density$density)))

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

  copula_time_plot <- suppressWarnings(plot_copula_fit(
    data = dat,
    copula_dist = copula_screen,
    response_var = "response",
    margin_dist = gamlss.dist::NO(),
    subject_var = "subject",
    time_var = "time",
    by_time = TRUE,
    plot = FALSE
  ))
  expect_s3_class(copula_time_plot$plot, "ggplot")
  expect_true("split_group" %in% names(copula_time_plot$pair_data))
  expect_true("split_group" %in% names(copula_time_plot$density))
  expect_gt(length(unique(copula_time_plot$pair_data$split_group)), 1L)

  clayton_dat <- simulate_longitudinal_dataset(
    n = 120,
    times = 1:4,
    margin_dist = gamlss.dist::NO(),
    margin_params = list(mu = 0, sigma = 1),
    copula_dist = "C",
    copula_params = list(theta = 3),
    seed = 159
  )
  clayton_plot <- suppressWarnings(plot_copula_fit(
    data = clayton_dat,
    copula_dist = "C",
    u_var = "u",
    plot = FALSE
  ))
  best_plot <- suppressWarnings(plot_copula_fit(
    data = clayton_dat,
    copula_dist = "best",
    u_var = "u",
    plot = FALSE
  ))
  expect_equal(clayton_plot$copula$family, "C")
  expect_gt(clayton_plot$copula$par, 0.5)
  expect_s3_class(clayton_plot$selection, "copula_selection")
  expect_equal(best_plot$copula$family, best_fit_family(best_plot$selection))
  expect_gt(abs(best_plot$copula$tau), 0.1)

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
