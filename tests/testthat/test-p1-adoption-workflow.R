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

  sim_model <- simulate(fit, nsim = 2, seed = 123)
  expect_s3_class(sim_model, "data.frame")
  expect_equal(dim(sim_model), c(length(fit$response), 2L))
  expect_true(all(is.finite(as.matrix(sim_model))))

  sim_copula_compat <- simulate(fit, nsim = 1, seed = 123, type = "copula")
  expect_equal(nrow(sim_copula_compat), length(fit$response))
  expect_true(all(is.finite(sim_copula_compat$sim_1)))
  expect_error(
    simulate(fit, nsim = 1, seed = 123, type = "marginal"),
    "fitted copula model"
  )
})

test_that("T202a simulate supports newdata panels", {
  dat <- make_fixture_factor_time_interaction(n_subject = 12L)
  fit <- fit_fixture_model(dat, include_dlcopdpar = TRUE)

  nd_new <- dat[seq_len(6), , drop = FALSE]
  nd_new$id <- nd_new$id + 1000L
  nd_new$y <- NA_real_
  nd_new <- nd_new[c(4, 1, 6, 2, 5, 3), , drop = FALSE]

  sim_new <- simulate(fit, nsim = 2, seed = 321, newdata = nd_new)
  expect_s3_class(sim_new, "data.frame")
  expect_equal(dim(sim_new), c(nrow(nd_new), 2L))
  expect_true(all(is.finite(as.matrix(sim_new))))
  expect_equal(
    simulate(fit, nsim = 2, seed = 321, newdata = nd_new),
    sim_new
  )

  nd_bad <- nd_new
  nd_bad$gender <- factor("X", levels = c(levels(dat$gender), "X"))
  expect_error(
    simulate(fit, nsim = 1, seed = 321, newdata = nd_bad),
    "contains level(s) not seen during fitting",
    fixed = TRUE
  )
})

test_that("T202b simulate newdata supports numeric-time extensions", {
  dat <- make_fixture_numeric_time(n_subject = 10L)
  fit <- fit_fixture_model(
    dat,
    include_dlcopdpar = TRUE,
    theta_formula = "~ time_raw",
    max_outer_iter = 2,
    max_inner_iter = 2
  )

  subject_data <- dat[!duplicated(dat$id), c("id", "gender", "age"), drop = FALSE]
  nd_more <- merge(
    expand.grid(
      id = subject_data$id[seq_len(3L)],
      time_raw = 1:4,
      KEEP.OUT.ATTRS = FALSE,
      stringsAsFactors = FALSE
    ),
    subject_data,
    by = "id",
    sort = FALSE
  )
  nd_more <- nd_more[order(nd_more$id, nd_more$time_raw), , drop = FALSE]
  nd_more$y <- NA_real_

  sim_more <- simulate(fit, nsim = 1, seed = 654, newdata = nd_more)
  expect_equal(dim(sim_more), c(nrow(nd_more), 1L))
  expect_true(all(is.finite(sim_more$sim_1)))
})

test_that("T202c simulate newdata breaks dependence across time-grid gaps", {
  diag_data <- list(
    family = "NO",
    params = list(mu = c(0, 0), sigma = c(1, 1))
  )
  fit_data <- data.frame(
    subject = c(1, 1),
    time = c(1, 3),
    theta_fit = c(0.95, NA_real_),
    zeta_fit = c(0, NA_real_)
  )
  object <- list(copula_dist = "N")

  set.seed(202)
  expected <- gamlss.dist::qNO(c(stats::runif(1L), stats::runif(1L)), mu = c(0, 0), sigma = c(1, 1))
  set.seed(202)
  actual <- gamlss.longitudinal:::.gl_simulate_copula_matrix(
    object,
    diag_data,
    nsim = 1,
    fit_data = fit_data,
    time_levels = 1:3
  )

  expect_equal(as.numeric(actual[, 1]), expected)
})

test_that("T203 check_model returns basic-check diagnostic object", {
  dat <- make_fixture_factor_time_interaction(n_subject = 14L)
  fit <- fit_fixture_model(dat, include_dlcopdpar = TRUE)

  check_out <- capture_warnings(check_model(fit, include_plots = FALSE, dependence_cor_cutoff = 0.3))
  chk <- check_out$value

  expect_s3_class(chk, "gamlss_longitudinal_check")
  expect_true(all(c(
    "model", "fit", "convergence", "scores", "pit", "tail",
    "basic_checks", "basic_checks_passed", "basic_checks_result",
    "checks", "warnings"
  ) %in% names(chk)))
  expect_false(any(c("recommendation", "decisions") %in% names(chk)))
  expect_true(is.data.frame(chk$scores))
  expect_true(is.data.frame(chk$basic_checks))
  expect_true(is.data.frame(chk$checks))
  expect_true(is.data.frame(chk$warnings))
  expect_equal(nrow(chk$basic_checks), 5L)
  expect_equal(
    chk$basic_checks$area,
    c("Convergence", "Marginal fit", "Tail fit", "Copula fit", "Variance calculation")
  )
  expect_true(all(chk$basic_checks$status %in% c("PASS", "FAIL", "REVIEW")))
  expect_true(all(chk$warnings$status == "FAIL"))
  expect_true(is.logical(chk$basic_checks_passed))
  expect_true(chk$basic_checks_result %in% c("passed", "review", "failed"))
  expect_equal(chk$residual_dependence$cutoff, 0.3)
  expect_true("n_pairs" %in% names(chk$residual_dependence))

  copula_diag <- suppressWarnings(plot_copula_diagnostics(fit, plot = FALSE, residual_lags = 1))
  expect_equal(
    chk$residual_dependence$normal_score_cor,
    copula_diag$residual_lag_summary$cor_z,
    tolerance = 1e-12
  )
  expect_equal(chk$residual_dependence$n_pairs, copula_diag$residual_lag_summary$n_pairs)

  printed <- capture.output(print(chk))
  expect_lt(grep("Basic Checks", printed), grep("Scores", printed))
  expect_true(any(grepl("Result:", printed, fixed = TRUE)))
  expect_true(any(grepl("broader model diagnostics should also be reviewed", printed, fixed = TRUE)))
  expect_false(any(grepl("Detailed checks", printed, fixed = TRUE)))
  expect_false(any(grepl("Failed checks", printed, fixed = TRUE)))
})

test_that("T203a check_model warns only for failed basic checks", {
  dat <- make_fixture_factor_time_interaction(n_subject = 10L)
  fit <- fit_fixture_model(dat, include_dlcopdpar = TRUE)

  fit_failed <- fit
  fit_failed$convergence$converged <- FALSE
  failed_out <- capture_warnings(check_model(fit_failed, include_plots = FALSE, dependence_cor_cutoff = 0.99))
  chk_failed <- failed_out$value
  expect_true(any(grepl("Basic model checks failed", failed_out$warnings, fixed = TRUE)))
  expect_equal(chk_failed$basic_checks_result, "failed")
  expect_true(any(chk_failed$warnings$area == "Convergence"))

  fit_review <- fit
  fit_review$convergence$converged <- TRUE
  chk_review <- expect_no_warning(
    check_model(fit_review, include_vcov = TRUE, numderiv = TRUE, include_plots = FALSE, dependence_cor_cutoff = 0.99)
  )
  expect_equal(chk_review$basic_checks_result, "review")
  expect_true(any(chk_review$basic_checks$status == "REVIEW"))
  expect_equal(nrow(chk_review$warnings), 0L)

  tail_checks <- gamlss.longitudinal:::.gl_check_table(
    summary_obj = list(
      fit = list(),
      convergence = list(converged = TRUE)
    ),
    scores = data.frame(),
    pit_stats = data.frame(ks_p_value = 0.5),
    tail_stats = data.frame(tail_ratio_max = 2.1),
    lag1_cor = 0,
    dependence_cor_cutoff = 0.25,
    vcov_method = NA_character_
  )
  expect_equal(tail_checks$status[tail_checks$area == "Tail fit"], "FAIL")
  expect_equal(gamlss.longitudinal:::.gl_basic_checks_result(tail_checks), "failed")
})

test_that("T203b check_missingness screens response missingness against observed predictors", {
  dat <- make_fixture_factor_time_interaction(n_subject = 40L)
  set.seed(203)
  miss_prob <- stats::plogis(
    -1.4 +
      0.025 * (dat$age - mean(dat$age)) +
      0.35 * (dat$time_raw == levels(dat$time_raw)[3])
  )
  dat$y[stats::runif(nrow(dat)) < miss_prob] <- NA_real_

  chk <- check_missingness(
    dat,
    response_var = "y",
    time_var = "time_raw",
    subject_var = "id",
    predictors = c("time_raw", "gender", "age")
  )

  expect_s3_class(chk, "gamlss_longitudinal_missingness_check")
  expect_equal(chk$response$n_missing, sum(is.na(dat$y)))
  expect_true(chk$assessment %in% c("covariate_related_missingness", "no_detected_covariate_association"))
  expect_true(all(c("term", "p_value") %in% names(chk$terms)))
  expect_output(print(chk), "Response Missingness Check")

  no_missing <- check_missingness(dat[!is.na(dat$y), ], response_var = "y")
  expect_equal(no_missing$assessment, "no_missing_responses")
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
  expect_false(any(grepl("supported_by_longitudinal", capture.output(print(screen)))))
  expect_equal(screen$best_fit$family_name, screen$family[[1L]])
  expect_equal(best_fit(screen)$family_name, screen$family[[1L]])
  expect_equal(best_fit_family(screen)$family[1], screen$family[[1L]])

  dat <- data.frame(response = y)
  from_data <- suppressWarnings(suppressMessages(screen_margin(dat, response_var = "response", type = "realAll", try.gamlss = FALSE, trace = FALSE)))
  expect_s3_class(from_data, "margin_selection")
  expect_s3_class(from_data, "margin_screen")
  expect_equal(attr(from_data, "selected"), from_data$family[[1L]])
  expect_equal(best_fit_family(from_data)$family[1], from_data$family[[1L]])

  dat_time <- data.frame(
    response = y + rep(c(-0.4, 0.2, 0.7), length.out = length(y)),
    time = rep(1:3, length.out = length(y))
  )
  time_screen <- suppressWarnings(suppressMessages(select_margin(
    dat_time,
    response_var = "response",
    time_var = "time",
    time_intercepts = TRUE,
    type = "realAll",
    families = c("NO", "TF"),
    try.gamlss = FALSE,
    trace = FALSE
  )))
  expect_s3_class(time_screen, "margin_selection")
  expect_true(isTRUE(attr(time_screen, "time_intercepts")))
  expect_equal(attr(time_screen, "time_var"), "time")
  expect_true(all(time_screen$screen_model == "time_intercepts"))
  expect_true(all(c("AIC", "pooled_AIC") %in% names(time_screen)))
  expect_true(all(is.finite(time_screen$AIC)))
  expect_equal(attr(time_screen, "selected"), time_screen$family[[1L]])

  expect_error(
    select_margin(dat_time, response_var = "response", time_intercepts = TRUE, type = "realAll"),
    "time_var"
  )
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
  check_out <- capture_warnings(check_model(fit, include_plots = FALSE))
  expect_s3_class(check_out$value, "gamlss_longitudinal_check")
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

test_that("T208c wald_test resolves factor formula terms to all level coefficients", {
  dat <- make_fixture_factor_time_interaction(n_subject = 12L)
  treatments <- c("active", "control", "placebo")
  dat$treatment <- factor(treatments[(dat$id %% length(treatments)) + 1L], levels = treatments)
  dat$treatmentcontrol_score <- as.numeric(dat$id) / max(dat$id)
  fit <- fit_fixture_model(
    dat,
    include_dlcopdpar = TRUE,
    mu_formula = "y ~ time_raw + treatment + treatmentcontrol_score",
    sigma_formula = "~ time_raw + treatment"
  )

  expected_terms <- paste(
    "mu",
    colnames(fit$model_matrix$x$mu)[attr(fit$model_matrix$x$mu, "assign") == 2L],
    sep = "."
  )
  expect_equal(expected_terms, c("mu.treatmentcontrol", "mu.treatmentplacebo"))

  wt <- wald_test(fit, terms = "mu.treatment", method = "numderiv", progress = FALSE)
  expect_equal(wt$term, expected_terms)
  expect_false("mu.treatmentcontrol_score" %in% wt$term)

  joint <- wald_test(fit, terms = "mu.treatment", joint = TRUE, method = "numderiv", progress = FALSE)
  expect_equal(joint$df, length(expected_terms))
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
  expect_true(all(c("estimate", "bootstrap_mean", "bootstrap_se", "conf.low", "conf.high", "reps") %in% names(boot$summary)))
  expect_equal(nrow(boot$replicates), 2L)

  prefix_terms <- names(fit$par)[startsWith(names(fit$par), "mu.gender")]
  boot_prefix <- bootstrap_inference(
    fit,
    R = 1,
    terms = "mu.gender",
    seed = 124,
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
