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

test_that("T152 plot.terms interaction rendering metadata includes factor levels", {
  dat <- make_fixture_factor_time_interaction(n_subject = 16L)
  dat$time_raw <- factor(as.character(dat$time_raw), levels = levels(dat$time_raw), ordered = FALSE)
  fit <- fit_fixture_model(dat, include_dlcopdpar = TRUE)

  suppressPackageStartupMessages(library(grid))

  pt <- suppressWarnings(plot.terms(fit, data = dat))
  pf <- pt$fixed_terms

  mu_entries <- pf$mu
  entry_names <- names(mu_entries)
  interaction_entries <- entry_names[grepl(":", entry_names, fixed = TRUE) & grepl("|", entry_names, fixed = TRUE)]

  expect_equal(length(interaction_entries), 2)
  panel_levels <- vapply(mu_entries[interaction_entries], function(e) e$panel_level, character(1))
  expect_setequal(panel_levels, levels(dat$gender))
  expect_true(all(vapply(mu_entries[interaction_entries], function(e) length(e$levels), integer(1)) == length(levels(dat$time_raw))))
})

test_that("T153 plot methods smoke test return dashboard structures", {
  dat <- make_fixture_factor_time_interaction(n_subject = 14L)
  fit <- fit_fixture_model(dat, include_dlcopdpar = TRUE)

  suppressPackageStartupMessages(library(grid))

  pt <- suppressWarnings(plot.terms(fit, data = dat))
  expect_true(is.list(pt))
  expect_true(all(c("smooth_terms", "fixed_terms", "dashboard") %in% names(pt)))

  pdiag <- suppressWarnings(plot(fit, data = dat))
  expect_true(is.list(pdiag))
  expect_true(all(c("diagnostics", "forecasts", "dashboard") %in% names(pdiag)))
})
