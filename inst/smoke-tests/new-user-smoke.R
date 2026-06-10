if (file.exists("DESCRIPTION") && requireNamespace("pkgload", quietly = TRUE)) {
  pkgload::load_all(quiet = TRUE)
} else {
  suppressPackageStartupMessages(library(gamlss.longitudinal))
}
suppressPackageStartupMessages(library(gamlss.dist))

set.seed(1)
dat <- simulate_longitudinal_dataset(
  n = 20,
  times = 0:2,
  margin_dist = NO(),
  copula_dist = "N",
  margin_params = list(mu = function(d) 1 + 0.2 * d$time, sigma = 0.6),
  copula_params = list(theta = 0.25),
  seed = 1
)

fit_no <- gamlss_longitudinal(
  dataset = dat,
  margin_dist = NO(),
  copula_dist = "N",
  time_var = "time",
  subject_var = "subject",
  mu.formula = response ~ time,
  sigma.formula = ~ 1,
  theta.formula = ~ 1,
  max_outer_iter = 5,
  max_inner_iter = 5,
  compute_vcov = FALSE,
  verbose = 0
)

stopifnot(inherits(fit_no, "gamlss.longitudinal"))
stopifnot(all(is.finite(predict(fit_no, type = "mean"))))
stopifnot(all(is.finite(predict(fit_no, type = "quantile", probs = c(0.1, 0.5, 0.9))$q05)))
stopifnot(inherits(check_model(fit_no, include_vcov = FALSE), "gamlss_longitudinal_check"))
stopifnot(ncol(simulate(fit_no, nsim = 2, seed = 1)) == 2)

dat_pos <- simulate_longitudinal_dataset(
  n = 20,
  times = 0:2,
  margin_dist = GA(),
  copula_dist = "N",
  margin_params = list(mu = function(d) exp(1 + 0.1 * d$time), sigma = 0.4),
  copula_params = list(theta = 0.25),
  seed = 2
)

fit_ga <- gamlss_longitudinal(
  dataset = dat_pos,
  margin_dist = GA(),
  copula_dist = "N",
  time_var = "time",
  subject_var = "subject",
  mu.formula = response ~ time,
  sigma.formula = ~ 1,
  theta.formula = ~ 1,
  max_outer_iter = 5,
  max_inner_iter = 5,
  compute_vcov = FALSE,
  verbose = 0
)

bench <- benchmark_standard_models(
  data = dat,
  formula = response ~ time,
  subject_var = "subject",
  family = "gaussian",
  comparators = "glm",
  fit = fit_no
)

stopifnot(inherits(fit_ga, "gamlss.longitudinal"))
stopifnot(inherits(model_spec(fit_no), "gamlss_longitudinal_model_spec"))
stopifnot(is.data.frame(reporting_table(fit_no, dat, by = "time", threshold = 1)))
stopifnot(inherits(bench, "gamlss_longitudinal_benchmark"))

message("gamlss.longitudinal smoke test completed successfully.")
