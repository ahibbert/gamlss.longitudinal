#!/usr/bin/env Rscript

script_argument <- grep("^--file=", commandArgs(trailingOnly = FALSE), value = TRUE)
script_path <- normalizePath(sub("^--file=", "", script_argument[[1L]]), winslash = "/", mustWork = TRUE)
root <- normalizePath(file.path(dirname(script_path), "..", "..", ".."), winslash = "/", mustWork = TRUE)
setwd(root)
suppressPackageStartupMessages(devtools::load_all(root, quiet = TRUE))

args <- commandArgs(TRUE)
fit <- readRDS(args[[1L]])
if (!inherits(fit, "gamlss.longitudinal")) fit <- fit$fit %||% fit$value
output_path <- args[[2L]]
threshold <- if (length(args) >= 3L) as.numeric(args[[3L]]) else 2
interval_floor <- if (length(args) >= 4L) as.numeric(args[[4L]]) else 4

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

beta <- fit$par
selected <- names(beta)[abs(beta) >= threshold]
trace <- vector("list", length(selected))
for (index in seq_along(selected)) {
  coefficient <- selected[[index]]
  start <- beta[[coefficient]]
  start_loglik <- evaluate(beta)
  radius <- max(interval_floor, abs(start))
  result <- stats::optimize(
    function(value) {
      candidate <- beta
      candidate[[coefficient]] <- value
      evaluate(candidate)
    },
    interval = start + c(-radius, radius),
    maximum = TRUE,
    tol = 1e-3
  )
  if (is.finite(result$objective) && result$objective >= start_loglik) {
    beta[[coefficient]] <- result$maximum
  }
  trace[[index]] <- data.frame(
    coefficient = coefficient,
    start = start,
    end = beta[[coefficient]],
    start_logLik = start_loglik,
    end_logLik = evaluate(beta),
    stringsAsFactors = FALSE
  )
}
fit$par <- beta
attr(fit, "coordinate_rescue_trace") <- do.call(rbind, trace)
saveRDS(fit, output_path, version = 3L)
print(attr(fit, "coordinate_rescue_trace"))
