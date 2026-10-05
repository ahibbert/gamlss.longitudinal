#!/usr/bin/env Rscript

# Complete the covariance availability audit without recomputing the exact
# Hessian already produced by validate-fit.R. For ACTG's segmented likelihood,
# the production inference contract uses subject-clustered sandwich covariance;
# this script uses the audited Hessian as its bread and recomputes the same
# finite-difference cluster score matrix used by vcov(..., method = "sandwich").

script_argument <- grep("^--file=", commandArgs(trailingOnly = FALSE), value = TRUE)
script_path <- normalizePath(sub("^--file=", "", script_argument[[1L]]),
                             winslash = "/", mustWork = TRUE)
root <- normalizePath(file.path(dirname(script_path), "..", "..", ".."),
                      winslash = "/", mustWork = TRUE)
setwd(root)
suppressPackageStartupMessages(devtools::load_all(root, quiet = TRUE))

args <- commandArgs(TRUE)
fit_path <- normalizePath(args[[1L]], winslash = "/", mustWork = TRUE)
audit_path <- normalizePath(args[[2L]], winslash = "/", mustWork = TRUE)
output_path <- args[[3L]]

fit_container <- readRDS(fit_path)
fit <- if (inherits(fit_container, "gamlss.longitudinal")) {
  fit_container
} else if (inherits(fit_container$fit, "gamlss.longitudinal")) {
  fit_container$fit
} else {
  fit_container$value
}
audit <- readRDS(audit_path)
stopifnot(inherits(fit, "gamlss.longitudinal"), is.matrix(audit$hessian))

H <- audit$hessian[names(fit$par), names(fit$par), drop = FALSE]
information_diag <- diag(-(H + t(H)) / 2)
gradient <- list(
  gradient = audit$direct_gradient[names(fit$par)],
  steps = gamlss.longitudinal:::.gl_vcov_score_steps(fit$par, 1e-4),
  scaled_max = max(abs(audit$direct_gradient[names(fit$par)]) /
                     sqrt(information_diag))
)

warnings <- character()
timing <- system.time({
  solved <- withCallingHandlers(
    gamlss.longitudinal:::.gl_solve_hessian_vcov(
      H,
      parameter_names = names(fit$par),
      inference = inference_control("standard"),
      source = "audited_analytical",
      gradient = gradient
    ),
    warning = function(warning) {
      warnings <<- c(warnings, conditionMessage(warning))
      invokeRestart("muffleWarning")
    }
  )
  scores <- gamlss.longitudinal:::.gl_cluster_score_matrix(
    fit,
    par_cov = fit$par,
    par_s = fit$par_s,
    cluster = fit$response_subject,
    h = 1e-5
  )
  bread <- solved$vcov[colnames(scores), colnames(scores), drop = FALSE]
  correction <- gamlss.longitudinal:::.gl_sandwich_correction(
    n_obs = length(fit$response),
    n_clusters = nrow(scores),
    n_parameters = ncol(scores),
    adjust = TRUE
  )
  meat <- crossprod(scores) * correction
  robust <- bread %*% meat %*% bread
  rownames(robust) <- colnames(robust) <- colnames(scores)
})

result <- list(
  fit_path = fit_path,
  audit_path = audit_path,
  model_based_vcov = solved$vcov,
  sandwich_vcov = robust,
  model_based_se = sqrt(diag(solved$vcov)),
  sandwich_se = sqrt(diag(robust)),
  hessian_diagnostics = solved$hessian_diagnostics,
  cluster_score_sum_max = max(abs(colSums(scores))),
  joint_logLik_from_clusters = sum(attr(scores, "base_contributions")),
  joint_logLik = attr(scores, "joint_loglik"),
  correction = correction,
  warnings = unique(warnings),
  model_based_available = all(is.finite(solved$vcov)) &&
    all(diag(solved$vcov) > 0),
  sandwich_available = all(is.finite(robust)) && all(diag(robust) > 0),
  elapsed_sec = unname(timing[["elapsed"]])
)
saveRDS(result, output_path, version = 3L)
cat(sprintf(
  paste0("model_based=%s sandwich=%s clusters=%d score_sum=%.8g ",
         "logLik_delta=%.8g warnings=%d elapsed=%.3f\n"),
  result$model_based_available,
  result$sandwich_available,
  nrow(scores),
  result$cluster_score_sum_max,
  result$joint_logLik_from_clusters - result$joint_logLik,
  length(result$warnings),
  result$elapsed_sec
))
