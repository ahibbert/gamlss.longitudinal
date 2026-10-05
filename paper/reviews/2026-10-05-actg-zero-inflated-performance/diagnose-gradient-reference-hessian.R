#!/usr/bin/env Rscript

script_argument <- grep("^--file=", commandArgs(trailingOnly = FALSE), value = TRUE)
script_path <- normalizePath(sub("^--file=", "", script_argument[[1L]]), winslash = "/", mustWork = TRUE)
root <- normalizePath(file.path(dirname(script_path), "..", "..", ".."), winslash = "/", mustWork = TRUE)
setwd(root)
suppressPackageStartupMessages(devtools::load_all(root, quiet = TRUE))

args <- commandArgs(TRUE)
fit_container <- readRDS(args[[1L]])
fit <- if (inherits(fit_container, "gamlss.longitudinal")) fit_container else fit_container$fit %||% fit_container$value
if (length(args) >= 3L && identical(args[[3L]], "bfgs")) {
  source_container <- readRDS(fit_container$input_path)
  fit <- if (inherits(source_container, "gamlss.longitudinal")) source_container else source_container$fit %||% source_container$value
  unpacked <- .gl_unpack_cg_beta(
    fit_container$bfgs$par,
    par_cov_template = fit$par,
    par_s_template = fit$par_s
  )
  fit$par <- unpacked$par_cov
  fit$par_s <- unpacked$par_s
}
reference <- readRDS(args[[2L]])

copula_link <- get_copula_dist(fit$copula_dist)$copula_link
pair_cache <- build_copula_pair_cache(fit$response, fit$response_margin, fit$response_subject)
evaluate <- function(beta, margin_derivatives = FALSE) {
  eta <- calc_eta(beta, fit$model_matrix, fit$margin_dist, copula_link, fit$par_s)
  likelihood <- calc_likelihood_minimal(
    eta$eta_inv,
    mm = fit$model_matrix$x,
    margin_dist = fit$margin_dist,
    copula_dist = fit$copula_dist,
    calc_d2 = FALSE,
    response = fit$response,
    response_margin = fit$response_margin,
    response_subject = fit$response_subject,
    pair_cache = pair_cache,
    calc_margin_deriv = margin_derivatives
  )
  list(eta = eta, likelihood = likelihood)
}

base <- evaluate(fit$par, TRUE)
log_lik <- as.numeric(base$likelihood$log_lik[["joint"]])
analytical <- .cg_analytical_gradient(
  fit$par, fit$model_matrix, base$eta, base$likelihood,
  fit$margin_dist, fit$copula_dist, isTRUE(fit$include_dlcopdpar),
  fit$response, fit$response_margin, fit$response_subject
)[names(fit$par)]
direct <- vapply(seq_along(fit$par), function(index) {
  step <- 1e-4 * max(1, abs(fit$par[[index]]))
  plus <- minus <- fit$par
  plus[[index]] <- plus[[index]] + step
  minus[[index]] <- minus[[index]] - step
  lp <- tryCatch(as.numeric(evaluate(plus)$likelihood$log_lik[["joint"]]), error = function(e) NA_real_)
  lm <- tryCatch(as.numeric(evaluate(minus)$likelihood$log_lik[["joint"]]), error = function(e) NA_real_)
  if (is.finite(lp) && is.finite(lm)) (lp - lm) / (2 * step) else NA_real_
}, numeric(1L))
names(direct) <- names(fit$par)
scale <- sqrt(reference$information_diag[names(fit$par)])

cat(sprintf(
  "logLik=%.15f direct=%.8g analytical=%.8g agreement=%.8g\n",
  log_lik,
  max(abs(direct) / scale),
  max(abs(analytical) / scale),
  max(abs(direct - analytical) / scale)
))
ranked <- order(abs(analytical) / scale, decreasing = TRUE)
print(utils::head(data.frame(
  coefficient = names(analytical)[ranked],
  analytical = analytical[ranked],
  direct = direct[ranked],
  scale = scale[ranked],
  scaled = abs(analytical[ranked]) / scale[ranked]
), 10L), row.names = FALSE)
