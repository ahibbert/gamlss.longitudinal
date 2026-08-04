predict.fake_reporting_fit <- function(object, newdata, type, probs = NULL, q = NULL, direction = NULL, ...) {
  n <- nrow(newdata)
  switch(
    type,
    mean = seq_len(n),
    mu = seq_len(n) + 10,
    median = seq_len(n) + 20,
    quantile = data.frame(
      subject = newdata$subject,
      time = newdata$time,
      response = NA_real_,
      q10 = seq_len(n) + 30,
      q90 = seq_len(n) + 40
    ),
    probability = data.frame(probability = seq_len(n) / 10),
    stop("unexpected prediction type", call. = FALSE)
  )
}
registerS3method("predict", "fake_reporting_fit", predict.fake_reporting_fit, envir = asNamespace("stats"))

test_that("model_spec builds printable fitted-model audit metadata", {
  object <- list(
    response = c(1, NA, Inf),
    response_var = "y",
    subject_var = "id",
    time_var = "visit",
    margin_dist = list(family = c("NO", "Normal"), mu.link = "identity", sigma.link = "log"),
    copula_dist = "N",
    formulas = list(mu = y ~ x, sigma = ~ 1),
    optim_method = "RS",
    convergence = list(converged = TRUE, stop_reason = "tolerance", outer_iterations = 3L, max_outer_iter = 20L),
    vcov_meta = list(precomputed = FALSE),
    vcov = list(hessian_diagnostics = list(rank = 2L))
  )
  class(object) <- "gamlss.longitudinal"

  spec <- model_spec(object)
  printed <- utils::capture.output(out <- print(spec))

  expect_s3_class(spec, "gamlss_longitudinal_model_spec")
  expect_identical(out, spec)
  expect_identical(spec$variables$response, "y")
  expect_identical(spec$distributions$margin, "NO")
  expect_identical(spec$optimisation$method, "RS")
  expect_identical(spec$missingness$n_missing_response, 1L)
  expect_identical(spec$missingness$n_nonfinite_response, 1L)
  expect_equal(spec$margin_links$parameter, c("mu", "sigma"))
  expect_true(any(grepl("GAMLSS Longitudinal Model Specification", printed)))
})

test_that("model_spec rejects non-fitted inputs", {
  expect_error(model_spec(list()), "fitted 'gamlss.longitudinal' object")
})

test_that("reporting_table summarizes predictions overall and by group", {
  object <- structure(list(), class = c("fake_reporting_fit", "gamlss.longitudinal"))
  newdata <- data.frame(
    subject = 1:4,
    time = c(1, 1, 2, 2),
    group = c("a", "a", "b", "b")
  )

  overall <- reporting_table(object, newdata = newdata, probs = c(0.1, 0.9), threshold = 5, direction = "above")
  grouped <- reporting_table(object, newdata = newdata, by = "group", probs = c(0.1, 0.9))

  expect_identical(overall$n, 4L)
  expect_equal(overall$mean, 2.5)
  expect_equal(overall$mu, 12.5)
  expect_equal(overall$median, 22.5)
  expect_equal(overall$q10, 32.5)
  expect_equal(overall$q90, 42.5)
  expect_equal(overall$prob_above_5, 0.25)
  expect_identical(grouped$group, c("a", "b"))
  expect_identical(grouped$n, c(2L, 2L))
  expect_equal(grouped$mean, c(1.5, 3.5))
})

test_that("reporting_table validates inputs before prediction", {
  object <- structure(list(), class = c("fake_reporting_fit", "gamlss.longitudinal"))

  expect_error(reporting_table(list(), data.frame(x = 1)), "fitted 'gamlss.longitudinal' object")
  expect_error(reporting_table(object), "'newdata' is required")
  expect_error(reporting_table(object, data.frame(x = 1), by = "missing"), "'by' column")
})

make_publication_helper_object <- function() {
  object <- structure(
    list(
      response = c(1, 2, NA),
      response_subject = c("a", "a", "b"),
      response_margin = c(1, 2, 1),
      dataset = data.frame(response = c(1, 2, NA)),
      margin_dist = gamlss.dist::NO(),
      copula_dist = "N",
      par = c(
        theta.intercept = 0.2,
        mu.intercept = 1,
        sigma.intercept = -0.1
      ),
      par_s = list(
        mu = list(`s(age)` = c(0.1, 0.2)),
        theta = list()
      ),
      df_s = list(
        mu = list(`s(age)` = 1.5),
        theta = list()
      ),
      calc_lik_out_end = list(
        log_lik = c(marginal = -10, copula = -2, joint = -12)
      )
    ),
    class = "gamlss.longitudinal"
  )
  V <- diag(c(0.04, 0.09, 0.16))
  dimnames(V) <- list(names(object$par), names(object$par))
  object$vcov <- list(
    vcov = list(overall = V),
    se = list(overall = stats::setNames(sqrt(diag(V)), names(object$par))),
    method = "analytical",
    method_requested = "analytical"
  )
  object$vcov_meta <- list(numderiv = FALSE, method = "analytical")
  object$convergence <- list(converged = TRUE, outer_iterations = 2L)
  object
}

test_that("tidy method returns broom-compatible coefficient columns", {
  object <- make_publication_helper_object()

  out <- generics::tidy(object, conf.level = 0.9)

  expect_s3_class(out, "data.frame")
  expect_true(all(c(
    "component", "parameter", "term", "coefficient", "estimate",
    "std.error", "statistic", "p.value", "conf.low", "conf.high"
  ) %in% names(out)))
  expect_equal(out$parameter, c("mu", "sigma", "theta"))
  expect_equal(out$term[1], "(Intercept)")
  expect_equal(out$std.error, c(0.3, 0.4, 0.2))
  expect_equal(out$component, c("margin", "margin", "copula"))
  expect_true(all(is.finite(out$conf.low)))
})

test_that("glance method returns one-row model fit metadata", {
  object <- make_publication_helper_object()

  out <- generics::glance(object, include_vcov = FALSE)

  expect_s3_class(out, "data.frame")
  expect_equal(nrow(out), 1L)
  expect_equal(out$margin_dist, "NO")
  expect_equal(out$copula_dist, "N")
  expect_equal(out$nobs, 2L)
  expect_equal(out$nobs_expanded, 3L)
  expect_equal(out$AIC, 33)
  expect_true(out$converged)
})

test_that("augment method adds fitted values and residuals to supplied data", {
  object <- structure(
    list(response_var = "y"),
    class = c("fake_reporting_fit", "gamlss.longitudinal")
  )
  newdata <- data.frame(y = c(2, 4, 6), x = 1:3)

  out <- generics::augment(object, newdata = newdata, type = "mean")

  expect_s3_class(out, "data.frame")
  expect_equal(out$.fitted, 1:3)
  expect_equal(out$.resid, c(1, 2, 3))
})

test_that("publication_table formats coefficient and model summaries", {
  object <- make_publication_helper_object()

  coefs <- publication_table(object, table = "coefficients", conf.level = 0.9, digits = 2, p_digits = 3)
  model <- publication_table(object, table = "model", include_vcov = FALSE, digits = 2)

  expect_s3_class(coefs, "gamlss_longitudinal_publication_table")
  expect_equal(names(coefs), c("Component", "Parameter", "Term", "Estimate", "Std. Error", "90% CI", "p-value", "Significance"))
  expect_equal(coefs$Estimate[1], "1.00")
  expect_match(coefs$`90% CI`[1], "[0.51, 1.49]", fixed = TRUE)
  expect_equal(coefs$Component, c("Margin", "Margin", "Copula"))

  expect_s3_class(model, "gamlss_longitudinal_publication_table")
  expect_equal(model$Value[model$Statistic == "Observed rows"], "2")
  expect_equal(model$Value[model$Statistic == "AIC"], "33.00")
})

test_that("publication_table formats grouped prediction summaries", {
  object <- structure(list(), class = c("fake_reporting_fit", "gamlss.longitudinal"))
  newdata <- data.frame(
    subject = 1:4,
    time = c(1, 1, 2, 2),
    group = c("a", "a", "b", "b")
  )

  out <- publication_table(
    object,
    table = "predictions",
    newdata = newdata,
    by = "group",
    probs = c(0.1, 0.9),
    threshold = 5,
    digits = 1
  )

  expect_s3_class(out, "gamlss_longitudinal_publication_table")
  expect_true(all(c("Group", "N", "Mean", "Mu", "Median", "Q10", "Q90", "Pr(above 5)") %in% names(out)))
  expect_equal(out$N, c(2L, 2L))
  expect_equal(out$Mean, c("1.5", "3.5"))
})

test_that("publication_table can emit latex output", {
  skip_if_not_installed("knitr")
  object <- make_publication_helper_object()

  out <- publication_table(
    object,
    table = "coefficients",
    conf.level = 0.9,
    digits = 2,
    output = "latex",
    caption = "Coefficient summary"
  )
  latex <- paste(as.character(out), collapse = "\n")

  expect_s3_class(out, "knitr_kable")
  expect_match(latex, "\\\\begin\\{tabular\\}", perl = TRUE)
  expect_match(latex, "\\\\toprule", perl = TRUE)
  expect_match(latex, "Coefficient summary", fixed = TRUE)
  expect_match(latex, "\\[0.51, 1.49\\]", perl = TRUE)
})
