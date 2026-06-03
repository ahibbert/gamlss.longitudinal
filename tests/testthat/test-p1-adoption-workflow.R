test_that("T201 predict method returns standard fitted outputs", {
  dat <- make_fixture_factor_time_interaction(n_subject = 14L)
  fit <- fit_fixture_model(dat, include_dlcopdpar = TRUE)

  p_response <- predict(fit, type = "response")
  expect_type(p_response, "double")
  expect_length(p_response, length(fit$response))
  expect_true(all(is.finite(p_response)))

  p_response_se <- predict(
    fit,
    type = "response",
    se.fit = TRUE,
    interval = "confidence",
    vcov_method = "numderiv",
    progress = FALSE
  )
  expect_s3_class(p_response_se, "data.frame")
  expect_true(all(c("fit", "se.fit", "conf.low", "conf.high") %in% names(p_response_se)))
  expect_equal(nrow(p_response_se), length(fit$response))

  p_params <- predict(fit, type = "parameters")
  expect_true(all(c("subject", "time", "response", "mu", "sigma") %in% names(p_params)))
  expect_equal(nrow(p_params), length(fit$response))

  p_quantile <- predict(fit, type = "quantile", probs = c(0.1, 0.5, 0.9))
  expect_true(all(c("q01", "q05", "q09") %in% names(p_quantile)))

  p_cdf <- predict(fit, type = "cdf")
  expect_true(all(c("subject", "time", "response", "q", "cdf") %in% names(p_cdf)))
  expect_true(all(p_cdf$cdf >= 0 & p_cdf$cdf <= 1))

  p_cdf_q <- predict(fit, type = "cdf", q = 2)
  expect_true(all(p_cdf_q$q == 2))
  expect_true(all(p_cdf_q$cdf >= 0 & p_cdf_q$cdf <= 1))

  p_prob <- predict(fit, type = "probability", q = 2, direction = "above")
  expect_true(all(c("q", "direction", "probability") %in% names(p_prob)))
  expect_true(all(p_prob$direction == "above"))
  expect_true(all(p_prob$probability >= 0 & p_prob$probability <= 1))

  p_density_y <- predict(fit, type = "density", y = 2)
  expect_true(all(c("y", "density") %in% names(p_density_y)))
  expect_true(all(p_density_y$y == 2))
  expect_true(all(is.finite(p_density_y$density)))
})

test_that("T202 simulate method returns fitted-length simulation columns", {
  dat <- make_fixture_factor_time_interaction(n_subject = 12L)
  fit <- fit_fixture_model(dat, include_dlcopdpar = TRUE)

  sim_marginal <- simulate(fit, nsim = 2, seed = 123, type = "marginal")
  expect_s3_class(sim_marginal, "data.frame")
  expect_equal(dim(sim_marginal), c(length(fit$response), 2L))
  expect_true(all(is.finite(as.matrix(sim_marginal))))

  sim_copula <- simulate(fit, nsim = 1, seed = 123, type = "copula")
  expect_equal(nrow(sim_copula), length(fit$response))
  expect_true(all(is.finite(sim_copula$sim_1)))
})

test_that("T203 check_model returns decision-oriented diagnostic object", {
  dat <- make_fixture_factor_time_interaction(n_subject = 14L)
  fit <- fit_fixture_model(dat, include_dlcopdpar = TRUE)

  chk <- check_model(fit, include_plots = FALSE)

  expect_s3_class(chk, "gamlss_longitudinal_check")
  expect_true(all(c("model", "fit", "convergence", "scores", "pit", "tail", "checks", "warnings", "recommendation", "decisions") %in% names(chk)))
  expect_true(is.data.frame(chk$scores))
  expect_true(is.data.frame(chk$checks))
  expect_true(is.data.frame(chk$warnings))
  expect_true(is.data.frame(chk$recommendation))
  expect_true(chk$recommendation$role %in% c("primary_candidate", "primary_candidate_with_caveats", "revise_before_primary"))
  expect_true(all(chk$warnings$severity %in% c("concern", "review", "note")))
  expect_true(is.character(chk$decisions))
  expect_true(length(chk$decisions) >= 1L)
  expect_message(plot(chk), "No plot objects stored")
})

test_that("T204 marginal_effects summarizes counterfactual parameter contrasts", {
  dat <- make_fixture_factor_time_interaction(n_subject = 14L)
  fit <- fit_fixture_model(dat, include_dlcopdpar = TRUE)

  eff <- marginal_effects(
    fit,
    newdata = dat,
    variable = "gender",
    values = levels(dat$gender),
    parameter = "mu",
    se.fit = TRUE,
    vcov_method = "numderiv",
    progress = FALSE
  )

  expect_s3_class(eff, "data.frame")
  expect_true(all(c("variable", "value", "parameter", "estimate", "std_error", "reference", "contrast", "conf.low", "conf.high") %in% names(eff)))
  expect_equal(nrow(eff), length(levels(dat$gender)))
  expect_equal(eff$contrast[1], 0)

  eff_alias <- effects(
    fit,
    newdata = dat,
    variable = "gender",
    values = levels(dat$gender),
    parameter = "mu",
    se.fit = FALSE
  )
  expect_equal(eff_alias$estimate, eff[, "estimate"])
})

test_that("T204b copula_time_summary prints dependence summaries", {
  dat <- make_fixture_factor_time_interaction(n_subject = 10L)
  fit <- fit_fixture_model(dat, include_dlcopdpar = TRUE)

  dep <- copula_time_summary(fit)

  expect_s3_class(dep, "copula_time_summary")
  expect_true(all(c("time_summary", "pair_summary", "fit_data", "pair_data") %in% names(dep)))
  expect_output(print(dep), "Copula Dependence Summary")
  expect_true(all(c("theta_fit", "tau_fit") %in% names(dep$time_summary)))
})

test_that("T205 select_margin and screen_margin return ordered candidate tables with best-fit accessors", {
  skip_if_not_installed("gamlss")

  y <- rnorm(30)
  screen <- suppressWarnings(suppressMessages(select_margin(y, type = "realAll", try.gamlss = FALSE, trace = FALSE)))

  expect_s3_class(screen, "margin_selection")
  expect_s3_class(screen, "margin_screen")
  expect_true(all(c("family", "AIC", "type") %in% names(screen)))
  expect_true(nrow(screen) >= 1L)
  expect_equal(attr(screen, "selected"), screen$family[[1L]])
  expect_true(all(c("rank", "delta_AIC", "supported_by_longitudinal") %in% names(screen)))
  expect_equal(screen$best_fit$family_name, screen$family[[1L]])
  expect_equal(best_fit(screen)$family_name, screen$family[[1L]])
  expect_equal(best_fit_family(screen)$family[1], screen$family[[1L]])

  dat <- data.frame(response = y)
  from_data <- suppressWarnings(suppressMessages(screen_margin(dat, response_var = "response", type = "realAll", try.gamlss = FALSE, trace = FALSE)))
  expect_s3_class(from_data, "margin_selection")
  expect_s3_class(from_data, "margin_screen")
  expect_equal(attr(from_data, "selected"), from_data$family[[1L]])
  expect_equal(best_fit_family(from_data)$family[1], from_data$family[[1L]])
})

test_that("T206 gamlss_longitudinal is the single public fitter", {
  dat <- make_fixture_factor_time_interaction(n_subject = 10L)

  fit <- suppressWarnings(gamlss_longitudinal(
    dataset = dat,
    margin_dist = NO(),
    time_var = "time_raw",
    subject_var = "id",
    mu.formula = y ~ time_raw + gender + age,
    sigma.formula = ~time_raw,
    theta.formula = ~1,
    copula_dist = "N",
    max_outer_iter = 2,
    max_inner_iter = 2,
    outer_stop_crit = 1,
    inner_stop_crit = 1,
    verbose = 0
  ))

  expect_s3_class(fit, "gamlss.longitudinal")
  expect_null(fit$workflow)
  expect_s3_class(check_model(fit, include_plots = FALSE), "gamlss_longitudinal_check")
})

test_that("T207 confint returns coefficient intervals users can cite", {
  dat <- make_fixture_factor_time_interaction(n_subject = 10L)
  fit <- fit_fixture_model(dat, include_dlcopdpar = TRUE)

  ci <- confint(fit, parm = names(fit$par)[seq_len(2L)], method = "numderiv", progress = FALSE)

  expect_true(is.matrix(ci))
  expect_equal(ncol(ci), 2L)
  expect_equal(nrow(ci), 2L)
})

test_that("T208 wald_test returns individual and joint hypothesis tests", {
  dat <- make_fixture_factor_time_interaction(n_subject = 10L)
  fit <- fit_fixture_model(dat, include_dlcopdpar = TRUE)
  terms <- names(fit$par)[seq_len(2L)]

  wt <- wald_test(fit, terms = terms, method = "numderiv", progress = FALSE)

  expect_s3_class(wt, "gamlss_longitudinal_wald_test")
  expect_false(isTRUE(attr(wt, "joint")))
  expect_true(all(c("term", "estimate", "std_error", "statistic", "p_value", "method") %in% names(wt)))
  expect_equal(wt$term, terms)
  expect_true(all(is.finite(wt$p_value)))

  joint <- wald_test(fit, terms = terms, joint = TRUE, method = "numderiv", progress = FALSE)
  expect_s3_class(joint, "gamlss_longitudinal_wald_test")
  expect_true(isTRUE(attr(joint, "joint")))
  expect_equal(joint$df, length(terms))
  expect_true(is.finite(joint$p_value))

  L <- matrix(c(1, -1), nrow = 1)
  colnames(L) <- terms
  contrast <- wald_test(fit, L = L, method = "numderiv", progress = FALSE)
  expect_true(isTRUE(attr(contrast, "joint")))
  expect_equal(contrast$df, 1)
})

test_that("T208b wald_test accepts coefficient-name prefixes", {
  dat <- make_fixture_factor_time_interaction(n_subject = 10L)
  fit <- fit_fixture_model(dat, include_dlcopdpar = TRUE)
  prefixed_terms <- names(fit$par)[startsWith(names(fit$par), "mu.gender")]

  expect_gt(length(prefixed_terms), 0L)

  wt <- wald_test(fit, terms = "mu.gender", method = "numderiv", progress = FALSE)
  expect_equal(wt$term, prefixed_terms)

  time_terms <- names(fit$par)[startsWith(names(fit$par), "mu.time_covariate") & !grepl(":", names(fit$par), fixed = TRUE)]
  wt_time <- wald_test(fit, terms = "mu.time_covariate", method = "numderiv", progress = FALSE)
  expect_equal(wt_time$term, time_terms)

  joint <- wald_test(fit, terms = "mu.gender", joint = TRUE, method = "numderiv", progress = FALSE)
  expect_true(isTRUE(attr(joint, "joint")))
  expect_equal(joint$df, length(prefixed_terms))
})

test_that("T209 likelihood_compare returns sequential LR summaries", {
  dat <- make_fixture_factor_time_interaction(n_subject = 10L)
  reduced <- suppressWarnings(gamlss_longitudinal(
    dataset = dat,
    margin_dist = NO(),
    time_var = "time_raw",
    subject_var = "id",
    mu.formula = y ~ time_raw + gender + age,
    sigma.formula = ~1,
    theta.formula = ~1,
    copula_dist = "N",
    compute_vcov = FALSE,
    max_outer_iter = 2,
    max_inner_iter = 2,
    outer_stop_crit = 1,
    inner_stop_crit = 1,
    verbose = 0
  ))
  full <- fit_fixture_model(dat, include_dlcopdpar = TRUE)

  cmp <- likelihood_compare(reduced = reduced, full = full)

  expect_s3_class(cmp, "gamlss_longitudinal_likelihood_compare")
  expect_equal(nrow(cmp), 2L)
  expect_true(all(c("model", "n_obs", "df", "logLik", "AIC", "BIC", "delta_df", "LR_statistic", "p_value") %in% names(cmp)))
  expect_true(all(is.finite(cmp$logLik)))
  expect_true(is.na(cmp$p_value[[1L]]))
})

test_that("T210 bootstrap_inference refits simulated responses", {
  dat <- make_fixture_factor_time_interaction(n_subject = 8L)
  fit <- fit_fixture_model(dat, include_dlcopdpar = TRUE)
  terms <- names(fit$par)[seq_len(2L)]

  boot <- bootstrap_inference(
    fit,
    R = 2,
    terms = terms,
    seed = 123,
    simulation_type = "marginal",
    fit_args = list(
      max_outer_iter = 1,
      max_inner_iter = 1,
      outer_stop_crit = 1,
      inner_stop_crit = 1,
      use_backtracking = FALSE
    )
  )

  expect_s3_class(boot, "gamlss_longitudinal_bootstrap")
  expect_s3_class(boot$summary, "data.frame")
  expect_equal(boot$R, 2L)
  expect_equal(boot$summary$term, terms)
  expect_true(all(c("estimate", "bootstrap_mean", "bootstrap_se", "conf.low", "conf.high", "successful_replicates") %in% names(boot$summary)))
  expect_equal(nrow(boot$replicates), 2L)

  prefix_terms <- names(fit$par)[startsWith(names(fit$par), "mu.gender")]
  boot_prefix <- bootstrap_inference(
    fit,
    R = 1,
    terms = "mu.gender",
    seed = 124,
    simulation_type = "marginal",
    fit_args = list(
      max_outer_iter = 1,
      max_inner_iter = 1,
      outer_stop_crit = 1,
      inner_stop_crit = 1,
      use_backtracking = FALSE
    )
  )
  expect_equal(boot_prefix$summary$term, prefix_terms)
})
