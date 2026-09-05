#!/usr/bin/env Rscript

# Small, checkpointed joint-model comparison. Each dataset gets three models:
# constant non-mean/correlation parameters, the preferred marginal-parameter
# covariates, and those marginal covariates plus prespecified theta covariates.

Sys.setenv(MEDICAL_CANDIDATE_DEFINE_ONLY = "true")
source(file.path(
  "paper", "reviews", "2026-09-05-medical-candidates", "run-screening.R"
))

joint_specs <- list(
  tlc = list(
    missingness = "error", max_elapsed_sec = 300,
    sigma = ~ treatment, nu = ~ treatment,
    theta = ~ treatment
  ),
  actg = list(
    missingness = "segment", max_elapsed_sec = 900,
    sigma = ~ treatment + baseline_log_z, nu = ~ 1,
    theta = ~ treatment + baseline_log_z
  ),
  skin = list(
    missingness = "segment", max_elapsed_sec = 900,
    sigma = ~ age_z, nu = ~ 1,
    theta = ~ treatment + age_z + exposure_z
  ),
  epilepsy = list(
    missingness = "error", max_elapsed_sec = 300,
    sigma = ~ treatment, nu = ~ 1,
    theta = ~ treatment + baseline_log_z
  )
)

joint_model_grid <- function(specification) {
  list(
    constant = list(sigma = ~ 1, nu = ~ 1, theta = ~ 1),
    margin_covariates = list(
      sigma = specification$sigma,
      nu = specification$nu,
      theta = ~ 1
    ),
    margin_and_theta_covariates = list(
      sigma = specification$sigma,
      nu = specification$nu,
      theta = specification$theta
    )
  )
}

fit_one <- function(name, model, force = FALSE) {
  data_spec <- candidate_specs[[name]]
  joint_spec <- joint_specs[[name]]
  screening_path <- file.path(output_dir, paste0(name, "-screening.rds"))
  if (!file.exists(screening_path)) {
    stop("Run the screening analysis first for ", name, ".")
  }
  screening <- readRDS(screening_path)
  data <- data_spec$read()
  formulas <- joint_model_grid(joint_spec)[[model]]
  fit_path <- file.path(output_dir, paste0(name, "-joint-", model, ".rds"))
  if (file.exists(fit_path) && !force) {
    message("Reusing ", fit_path)
    return(readRDS(fit_path))
  }

  message("Fitting ", name, ": ", model)
  control <- gamlss.longitudinal::gamlss_longitudinal_control(
    outer_tol = if (name %in% c("actg", "skin")) 1e-3 else 1e-4,
    max_outer_iter = 80L,
    max_elapsed_sec = joint_spec$max_elapsed_sec,
    rs = list(
      inner_tol = 1e-5,
      max_inner_iter = 150L,
      warm_start_joint = !identical(joint_spec$missingness, "segment"),
      warm_start_joint_iter = if (identical(joint_spec$missingness, "segment")) 0L else 5L
    )
  )
  attempt <- capture_fit(gamlss.longitudinal::gamlss_longitudinal(
    dataset = data,
    margin_dist = family_object(screening$selected_margin),
    copula_dist = screening$selected_copula,
    subject_var = data_spec$subject,
    time_var = data_spec$time,
    missingness = joint_spec$missingness,
    mu.formula = data_spec$mu,
    sigma.formula = formulas$sigma,
    nu.formula = formulas$nu,
    tau.formula = ~ 1,
    theta.formula = formulas$theta,
    zeta.formula = ~ 1,
    method = "RS",
    optimizer_control = control,
    compute_vcov = FALSE,
    verbose = 0
  ))
  result <- list(
    fit = if (inherits(attempt$value, "error")) NULL else attempt$value,
    error = if (inherits(attempt$value, "error")) conditionMessage(attempt$value) else "",
    warnings = attempt$warnings,
    elapsed_sec = attempt$elapsed,
    model = model,
    formulas = formulas,
    selected_margin = screening$selected_margin,
    selected_copula = screening$selected_copula
  )
  saveRDS(result, fit_path)
  result
}

summarize_joint <- function(result, name) {
  failed <- is.null(result$fit)
  fit <- result$fit
  converged <- !failed && isTRUE(fit$convergence$converged)
  coefficients <- if (failed) numeric() else stats::coef(fit)
  log_likelihood <- if (converged) tryCatch(stats::logLik(fit), error = function(e) NULL) else NULL
  stop_reason <- if (failed) {
    "error"
  } else if (is.null(fit$convergence$stop_reason)) {
    "unknown"
  } else {
    fit$convergence$stop_reason
  }
  data.frame(
    dataset = name,
    model = result$model,
    margin = result$selected_margin,
    copula = result$selected_copula,
    converged = converged,
    stop_reason = stop_reason,
    df = if (is.null(log_likelihood)) NA_real_ else attr(log_likelihood, "df"),
    logLik = if (is.null(log_likelihood)) NA_real_ else as.numeric(log_likelihood),
    AIC = if (is.null(log_likelihood)) NA_real_ else stats::AIC(fit),
    BIC = if (is.null(log_likelihood)) NA_real_ else stats::BIC(fit),
    elapsed_sec = result$elapsed_sec,
    theta_coefficients = paste(
      sprintf("%s=%.6g", names(coefficients)[grepl("^theta\\.", names(coefficients))],
        coefficients[grepl("^theta\\.", names(coefficients))]),
      collapse = "; "
    ),
    sigma_coefficients = paste(
      sprintf("%s=%.6g", names(coefficients)[grepl("^sigma\\.", names(coefficients))],
        coefficients[grepl("^sigma\\.", names(coefficients))]),
      collapse = "; "
    ),
    nu_coefficients = paste(
      sprintf("%s=%.6g", names(coefficients)[grepl("^nu\\.", names(coefficients))],
        coefficients[grepl("^nu\\.", names(coefficients))]),
      collapse = "; "
    ),
    warnings = paste(result$warnings, collapse = " | "),
    error = result$error,
    stringsAsFactors = FALSE
  )
}

args <- commandArgs(trailingOnly = TRUE)
force <- "--force" %in% args
requested <- setdiff(args, "--force")
if (!length(requested)) requested <- names(joint_specs)
unknown <- setdiff(requested, names(joint_specs))
if (length(unknown)) stop("Unknown dataset(s): ", paste(unknown, collapse = ", "))

rows <- list()
k <- 0L
for (name in requested) {
  for (model in names(joint_model_grid(joint_specs[[name]]))) {
    result <- fit_one(name, model, force = force)
    k <- k + 1L
    rows[[k]] <- summarize_joint(result, name)
    utils::write.csv(
      do.call(rbind, rows),
      file.path(output_dir, "joint-model-summary.csv"),
      row.names = FALSE
    )
    print(rows[[k]])
  }
}
