#!/usr/bin/env Rscript

script_argument <- grep("^--file=", commandArgs(trailingOnly = FALSE), value = TRUE)
script_path <- normalizePath(sub("^--file=", "", script_argument[[1L]]), winslash = "/", mustWork = TRUE)
root <- normalizePath(file.path(dirname(script_path), "..", "..", ".."), winslash = "/", mustWork = TRUE)
setwd(root)
suppressPackageStartupMessages(devtools::load_all(root, quiet = TRUE))

args <- commandArgs(TRUE)
input_path <- normalizePath(args[[1L]], winslash = "/", mustWork = TRUE)
output_path <- args[[2L]]
repetitions <- if (length(args) >= 3L) as.integer(args[[3L]]) else 3L
fast_zibnb <- if (length(args) >= 4L) as.logical(args[[4L]]) else TRUE
compare_zibnb <- if (length(args) >= 5L) as.logical(args[[5L]]) else FALSE
fast_t <- if (length(args) >= 6L) as.logical(args[[6L]]) else TRUE
fast_zibnb_family <- if (length(args) >= 7L) as.logical(args[[7L]]) else TRUE
options(gamlss.longitudinal.fast_zibnb_cdf_scores = fast_zibnb)
options(gamlss.longitudinal.fast_t_rectangle_scores = fast_t)
options(gamlss.longitudinal.fast_zibnb_family = fast_zibnb_family)

container <- readRDS(input_path)
fit <- if (inherits(container, "gamlss.longitudinal")) container else container$value
stopifnot(inherits(fit, "gamlss.longitudinal"))

dataset <- data.frame(
  response = fit$response,
  time = fit$response_margin,
  subject = fit$response_subject
)
pair_cache <- gamlss.longitudinal:::build_copula_pair_cache(
  dataset$response, dataset$time, dataset$subject
)
copula_link <- gamlss.longitudinal:::get_copula_dist(fit$copula_dist)$copula_link
eta_out <- gamlss.longitudinal:::calc_eta(
  fit$par,
  fit$model_matrix,
  fit$margin_dist,
  copula_link,
  par_s = fit$par_s
)
eta_inv <- eta_out$eta_inv
margin_eval_cache <- gamlss.longitudinal:::.build_margin_eval_cache(
  fit$margin_dist, calc_d2 = FALSE
)
calc_lik <- gamlss.longitudinal:::calc_likelihood_minimal(
  eta_inv,
  mm = fit$model_matrix$x,
  margin_dist = fit$margin_dist,
  copula_dist = fit$copula_dist,
  calc_d2 = FALSE,
  response = dataset$response,
  response_margin = dataset$time,
  response_subject = dataset$subject,
  pair_cache = pair_cache,
  margin_eval_cache = margin_eval_cache
)

time_call <- function(label, fn) {
  timings <- numeric(repetitions)
  value <- NULL
  for (index in seq_len(repetitions)) {
    timing <- system.time(value <- fn())
    timings[[index]] <- unname(timing[["elapsed"]])
  }
  data.frame(
    operation = label,
    repetition = seq_len(repetitions),
    elapsed_sec = timings,
    value_length = length(unlist(value, use.names = FALSE)),
    stringsAsFactors = FALSE
  )
}

score_call <- function(score_names) {
  gamlss.longitudinal:::.calc_discrete_rectangle_scores(
    eta_inv,
    mm = fit$model_matrix$x,
    margin_dist = fit$margin_dist,
    copula_dist = fit$copula_dist,
    response = dataset$response,
    response_margin = dataset$time,
    response_subject = dataset$subject,
    pair_cache = pair_cache,
    calc_lik = calc_lik,
    method = "analytical",
    score_names = score_names
  )
}

parameters <- intersect(
  names(fit$model_matrix$x),
  c("mu", "sigma", "nu", "tau", "theta", "zeta")
)
results <- list(time_call("all", function() score_call(NULL)))
for (parameter in parameters) {
  results[[length(results) + 1L]] <- time_call(parameter, function() score_call(parameter))
}
results[[length(results) + 1L]] <- time_call(
  "copula_theta_raw",
  function() gamlss.longitudinal:::.calc_copula_rectangle_par_derivative(
    calc_lik, fit$copula_dist, "theta"
  )
)
if ("zeta" %in% parameters) {
  results[[length(results) + 1L]] <- time_call(
    "copula_zeta_raw",
    function() gamlss.longitudinal:::.calc_copula_rectangle_par_derivative(
      calc_lik, fit$copula_dist, "zeta"
    )
  )
}

results <- do.call(rbind, results)
utils::write.csv(results, output_path, row.names = FALSE)
print(results)
if (compare_zibnb) {
  likelihood_with_option <- function(fast) {
    options(gamlss.longitudinal.fast_zibnb_family = fast)
    cache <- gamlss.longitudinal:::.build_margin_eval_cache(
      fit$margin_dist, calc_d2 = FALSE
    )
    system.time(value <- gamlss.longitudinal:::calc_likelihood_minimal(
      eta_inv,
      mm = fit$model_matrix$x,
      margin_dist = fit$margin_dist,
      copula_dist = fit$copula_dist,
      calc_d2 = FALSE,
      response = dataset$response,
      response_margin = dataset$time,
      response_subject = dataset$subject,
      pair_cache = pair_cache,
      margin_eval_cache = cache
    )) -> timing
    list(value = value, elapsed = unname(timing[["elapsed"]]))
  }
  reference_likelihood <- likelihood_with_option(FALSE)
  candidate_likelihood <- likelihood_with_option(TRUE)
  print(c(
    reference_logLik = reference_likelihood$value$log_lik[["joint"]],
    candidate_logLik = candidate_likelihood$value$log_lik[["joint"]],
    difference = candidate_likelihood$value$log_lik[["joint"]] -
      reference_likelihood$value$log_lik[["joint"]],
    reference_elapsed = reference_likelihood$elapsed,
    candidate_elapsed = candidate_likelihood$elapsed
  ))
  options(gamlss.longitudinal.fast_zibnb_cdf_scores = FALSE)
  options(gamlss.longitudinal.fast_t_rectangle_scores = FALSE)
  reference <- score_call(NULL)
  options(gamlss.longitudinal.fast_zibnb_cdf_scores = TRUE)
  options(gamlss.longitudinal.fast_t_rectangle_scores = TRUE)
  candidate <- score_call(NULL)
  comparison <- do.call(rbind, lapply(names(reference), function(parameter) {
    difference <- as.numeric(candidate[[parameter]] - reference[[parameter]])
    scale <- pmax(1, abs(reference[[parameter]]))
    data.frame(
      parameter = parameter,
      max_abs_difference = max(abs(difference)),
      max_scaled_difference = max(abs(difference) / scale),
      rms_difference = sqrt(mean(difference^2)),
      stringsAsFactors = FALSE
    )
  }))
  comparison_path <- sub("\\.csv$", "-comparison.csv", output_path, ignore.case = TRUE)
  utils::write.csv(comparison, comparison_path, row.names = FALSE)
  print(comparison)
}
