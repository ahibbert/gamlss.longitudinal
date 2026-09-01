#' Runs a separately optimised RS fit to provide safe starting values for a joint RS fit 
#' and extracts fit information needed to pass to start for the joint fit.
#'
#' @noRd
.gl_run_joint_warm_start <- function(
    start_from,
    method,
    include_dlcopdpar,
    warm_start_joint,
    warm_start_joint_iter,
    user_supplied_start,
    dataset_original,
    margin_dist,
    copula_dist,
    time_var,
    subject_var,
    mu.formula,
    sigma.formula,
    nu.formula,
    tau.formula,
    theta.formula,
    zeta.formula,
    optimizer_control,
    true_val,
    vcov_method,
    vcov_numderiv,
    use_Rcpp,
    lambda_start,
    lambda_penalty_K,
    verbose,
    fit_fn = NULL,
    warm_start_fit_call_fn = .gl_call_joint_warm_start_fit) {
  warm_start_info <- .gl_default_warm_start_info()
  warm_start_par_s <- NULL

  if (!.gl_should_run_joint_warm_start(
    method = method,
    include_dlcopdpar = include_dlcopdpar,
    warm_start_joint = warm_start_joint,
    warm_start_joint_iter = warm_start_joint_iter,
    user_supplied_start = user_supplied_start
  )) {
    return(list(
      start_from = start_from,
      warm_start_par_s = warm_start_par_s,
      warm_start_info = warm_start_info
    ))
  }

  if (is.null(fit_fn)) {
    fit_fn <- gamlss_longitudinal
  }

  if (verbose > 0) {
    cat(
      "\nRunning separate RS warm-start phase for ",
      warm_start_joint_iter,
      " outer iteration(s) before joint RS fit...\n",
      sep = ""
    )
  }

  warm_call <- warm_start_fit_call_fn(
    fit_fn = fit_fn,
    dataset_original = dataset_original,
    margin_dist = margin_dist,
    copula_dist = copula_dist,
    time_var = time_var,
    subject_var = subject_var,
    mu.formula = mu.formula,
    sigma.formula = sigma.formula,
    nu.formula = nu.formula,
    tau.formula = tau.formula,
    theta.formula = theta.formula,
    zeta.formula = zeta.formula,
    optimizer_control = optimizer_control,
    true_val = true_val,
    method = method,
    warm_start_joint_iter = warm_start_joint_iter,
    vcov_method = vcov_method,
    vcov_numderiv = vcov_numderiv,
    use_Rcpp = use_Rcpp,
    lambda_start = lambda_start,
    lambda_penalty_K = lambda_penalty_K
  )
  warm_fit <- warm_call$fit
  warm_output <- warm_call$output
  warm_warnings <- warm_call$warnings
  warm_err <- warm_call$error

  .gl_validate_joint_warm_start_fit(warm_fit, warm_err)
  warm_start_result <- .gl_joint_warm_start_result(
    warm_fit = warm_fit,
    warm_output = warm_output,
    warm_warnings = warm_warnings,
    warm_start_joint_iter = warm_start_joint_iter
  )

  if (verbose > 1 && length(warm_output) > 0) {
    cat(paste(warm_output, collapse = "\n"), "\n")
  }

  warm_start_result
}
