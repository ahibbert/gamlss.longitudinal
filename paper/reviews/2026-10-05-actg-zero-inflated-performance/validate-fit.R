#!/usr/bin/env Rscript

script_argument <- grep("^--file=", commandArgs(trailingOnly = FALSE), value = TRUE)
script_path <- normalizePath(sub("^--file=", "", script_argument[[1L]]), winslash = "/", mustWork = TRUE)
root <- normalizePath(file.path(dirname(script_path), "..", "..", ".."), winslash = "/", mustWork = TRUE)
setwd(root)
suppressPackageStartupMessages(devtools::load_all(root, quiet = TRUE))

args <- commandArgs(TRUE)
fit_path <- normalizePath(args[[1L]], winslash = "/", mustWork = TRUE)
output_path <- args[[2L]]
container <- readRDS(fit_path)
fit <- if (inherits(container, "gamlss.longitudinal")) {
  container
} else if (inherits(container$fit, "gamlss.longitudinal")) {
  container$fit
} else {
  container$value
}
stopifnot(inherits(fit, "gamlss.longitudinal"))

copula_link <- gamlss.longitudinal:::get_copula_dist(fit$copula_dist)$copula_link
pair_cache <- gamlss.longitudinal:::build_copula_pair_cache(
  fit$response, fit$response_margin, fit$response_subject
)
evaluate <- function(beta, margin_derivatives = FALSE) {
  eta <- gamlss.longitudinal:::calc_eta(
    beta, fit$model_matrix, fit$margin_dist, copula_link, fit$par_s
  )
  likelihood <- gamlss.longitudinal:::calc_likelihood_minimal(
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

total_timing <- system.time({
  base_timing <- system.time({
  base_evaluation <- evaluate(fit$par, margin_derivatives = TRUE)
  corrected_loglik <- as.numeric(base_evaluation$likelihood$log_lik[["joint"]])
  })
  analytical_timing <- system.time({
  analytical_gradient <- gamlss.longitudinal:::.cg_analytical_gradient(
    fit$par,
    fit$model_matrix,
    base_evaluation$eta,
    base_evaluation$likelihood,
    fit$margin_dist,
    fit$copula_dist,
    isTRUE(fit$include_dlcopdpar),
    fit$response,
    fit$response_margin,
    fit$response_subject
  )
  })
  direct_timing <- system.time({
  direct_gradient <- vapply(seq_along(fit$par), function(index) {
    step <- 1e-4 * max(1, abs(fit$par[[index]]))
    plus <- minus <- fit$par
    plus[[index]] <- plus[[index]] + step
    minus[[index]] <- minus[[index]] - step
    plus_value <- tryCatch(
      as.numeric(evaluate(plus)$likelihood$log_lik[["joint"]]),
      error = function(error) NA_real_
    )
    minus_value <- tryCatch(
      as.numeric(evaluate(minus)$likelihood$log_lik[["joint"]]),
      error = function(error) NA_real_
    )
    if (is.finite(plus_value) && is.finite(minus_value)) {
      (plus_value - minus_value) / (2 * step)
    } else if (is.finite(plus_value) && is.finite(corrected_loglik)) {
      (plus_value - corrected_loglik) / step
    } else if (is.finite(minus_value) && is.finite(corrected_loglik)) {
      (corrected_loglik - minus_value) / step
    } else {
      NA_real_
    }
  }, numeric(1L))
  names(direct_gradient) <- names(fit$par)
  analytical_gradient <- analytical_gradient[names(fit$par)]
  })
  hessian_timing <- system.time({
  hessian <- gamlss.longitudinal:::calc_analytical_hessian(
    fit, progress = FALSE, h = 1e-4
  )
  })
  information_diag <- diag(-(hessian + t(hessian)) / 2)
  scale <- sqrt(information_diag)
  finite <- is.finite(corrected_loglik) &&
    all(is.finite(analytical_gradient)) && all(is.finite(direct_gradient)) &&
    all(is.finite(scale)) && all(scale > 0)
  direct_scaled <- if (finite) max(abs(direct_gradient) / scale) else Inf
  analytical_scaled <- if (finite) max(abs(analytical_gradient) / scale) else Inf
  agreement_scaled <- if (finite) {
    max(abs(analytical_gradient - direct_gradient) / scale)
  } else Inf
})

result <- list(
  fit_path = fit_path,
  fit_converged = isTRUE(fit$convergence$converged),
  fit_stop_reason = fit$convergence$stop_reason,
  corrected_logLik = corrected_loglik,
  analytical_gradient = analytical_gradient,
  direct_gradient = direct_gradient,
  hessian = hessian,
  information_diag = information_diag,
  direct_scaled_gradient = direct_scaled,
  analytical_scaled_gradient = analytical_scaled,
  gradient_agreement_scaled = agreement_scaled,
  direct_pass = is.finite(direct_scaled) && direct_scaled <= 0.25,
  analytical_pass = is.finite(analytical_scaled) && analytical_scaled <= 0.25,
  strict_pass = is.finite(direct_scaled) && direct_scaled <= 0.05,
  analytical_strict_pass = is.finite(analytical_scaled) && analytical_scaled <= 0.05,
  agreement_pass = is.finite(agreement_scaled) && agreement_scaled <= 0.05,
  phase_elapsed_sec = c(
    base = unname(base_timing[["elapsed"]]),
    analytical_gradient = unname(analytical_timing[["elapsed"]]),
    direct_gradient = unname(direct_timing[["elapsed"]]),
    hessian = unname(hessian_timing[["elapsed"]])
  ),
  elapsed_sec = unname(total_timing[["elapsed"]])
)
saveRDS(result, output_path, version = 3L)
cat(sprintf(
  paste0("converged=%s stop=%s logLik=%.15f direct=%.8g analytical=%.8g ",
         "agreement=%.8g strict=%s analytical_strict=%s agreement_pass=%s elapsed=%.3f\n"),
  result$fit_converged, result$fit_stop_reason, result$corrected_logLik,
  result$direct_scaled_gradient, result$analytical_scaled_gradient,
  result$gradient_agreement_scaled, result$strict_pass,
  result$analytical_strict_pass, result$agreement_pass, result$elapsed_sec
))
