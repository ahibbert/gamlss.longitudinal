#!/usr/bin/env Rscript

script_argument <- grep("^--file=", commandArgs(trailingOnly = FALSE), value = TRUE)
script_path <- normalizePath(sub("^--file=", "", script_argument[[1L]]), winslash = "/", mustWork = TRUE)
root <- normalizePath(file.path(dirname(script_path), "..", "..", ".."), winslash = "/", mustWork = TRUE)
setwd(root)
suppressPackageStartupMessages(devtools::load_all(root, quiet = TRUE))

args <- commandArgs(TRUE)
fit <- readRDS(args[[1L]])
if (!inherits(fit, "gamlss.longitudinal")) fit <- fit$fit %||% fit$value
coefficient <- args[[2L]]
offsets <- as.numeric(strsplit(args[[3L]], ",", fixed = TRUE)[[1L]])
output_path <- args[[4L]]
stopifnot(coefficient %in% names(fit$par))

copula_link <- get_copula_dist(fit$copula_dist)$copula_link
pair_cache <- build_copula_pair_cache(
  fit$response, fit$response_margin, fit$response_subject
)
evaluate <- function(beta) {
  eta <- calc_eta(beta, fit$model_matrix, fit$margin_dist, copula_link, fit$par_s)
  as.numeric(calc_likelihood_minimal(
    eta$eta_inv,
    mm = fit$model_matrix$x,
    margin_dist = fit$margin_dist,
    copula_dist = fit$copula_dist,
    calc_d2 = FALSE,
    response = fit$response,
    response_margin = fit$response_margin,
    response_subject = fit$response_subject,
    pair_cache = pair_cache,
    calc_margin_deriv = FALSE
  )$log_lik[["joint"]])
}

base <- fit$par
values <- vapply(offsets, function(offset) {
  candidate <- base
  candidate[[coefficient]] <- candidate[[coefficient]] + offset
  evaluate(candidate)
}, numeric(1L))
result <- data.frame(
  coefficient = coefficient,
  base_value = base[[coefficient]],
  offset = offsets,
  value = base[[coefficient]] + offsets,
  logLik = values,
  delta_logLik = values - evaluate(base),
  stringsAsFactors = FALSE
)
utils::write.csv(result, output_path, row.names = FALSE)
print(result)
