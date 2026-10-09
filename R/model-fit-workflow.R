#' Model fit workflow is called from model-fit-entrypoint.R and contains the first 
#' of three primary steps of the fitting process, covering pre-fit work:
#' - Preprocessing and validation of input data and formulas (model-preprocess.R)
#' - Construction of model matrices and related components (model-matrix.R)
#' - Initialization of starting parameters and optimizer context (model-fit-setup.R)
#' 
#' @noRd
.gl_initial_optimizer_likelihood_valid <- function(
    optimizer_context,
    mm,
    margin_dist,
    copula_dist,
    dataset,
    likelihood_fn = calc_likelihood_minimal) {
  eta_out <- tryCatch(
    optimizer_context$rs_calc_eta(
      par_cov_current = optimizer_context$par_cov,
      par_s_current = optimizer_context$par_s
    ),
    error = identity
  )
  if (inherits(eta_out, "error")) return(FALSE)
  likelihood <- tryCatch(
    likelihood_fn(
      eta_inv = eta_out$eta_inv,
      mm = mm$x,
      margin_dist = margin_dist,
      copula_dist = copula_dist,
      response = dataset$response,
      response_margin = dataset$time,
      response_subject = dataset$subject,
      pair_cache = optimizer_context$pair_cache,
      margin_eval_cache = optimizer_context$margin_eval_cache,
      calc_margin_deriv = FALSE
    ),
    error = identity
  )
  !inherits(likelihood, "error") && isTRUE(likelihood$valid) &&
    is.finite(likelihood$log_lik[["joint"]])
}

#' @noRd
.gl_prepare_fit_workflow <- function(
    dataset,
    margin_dist,
    copula_dist,
    time_var,
    subject_var,
    missingness,
    mu.formula,
    sigma.formula,
    nu.formula,
    tau.formula,
    theta.formula,
    zeta.formula,
    include_dlcopdpar,
    start_from,
    verbose,
    true_val,
    method,
    optimizer_control,
    vcov_method,
    vcov_numderiv,
    use_Rcpp,
    lambda_start,
    lambda_penalty_K,
    data_fn = .gl_prepare_fit_data,
    matrix_fn = .gl_build_model_matrix_bundle,
    warm_start_fn = .gl_run_joint_warm_start,
    step_control_fn = .gl_normalize_step_controls,
    optimizer_context_fn = .gl_initialize_fit_optimizer_context,
    initial_likelihood_valid_fn = .gl_initial_optimizer_likelihood_valid) {
  rs_control <- optimizer_control$rs
  cg_control <- optimizer_control$cg
  shared_control <- optimizer_control$shared
  vcov_controls <- .gl_normalize_vcov_controls(vcov_method, vcov_numderiv)
  fit_controls <- list(
    method = method,
    user_supplied_start = if (inherits(start_from, "gamlss.longitudinal") || is.list(start_from)) {
      !is.null(start_from$par)
    } else {
      !all(is.na(start_from))
    },
    warm_start_joint = rs_control$warm_start_joint,
    warm_start_joint_iter = rs_control$warm_start_joint_iter,
    backtracking_max_halves = if (identical(method, "RS")) rs_control$backtracking_max_halves else cg_control$backtracking_max_halves,
    cg_max_stall = cg_control$max_stall,
    cg_max_delta = cg_control$max_delta,
    cg_lambda_update_every = cg_control$lambda_update_every,
    cg_max_lambda_updates = cg_control$max_lambda_updates,
    cg_raw_loglik_drop_tol = cg_control$raw_loglik_drop_tol,
    cg_line_search = cg_control$line_search,
    cg_max_line_search_evals = cg_control$max_line_search_evals,
    cg_gradient_method = cg_control$gradient_method,
    discrete_score_method = rs_control$discrete_score_method,
    cg_zeta_hessian = cg_control$zeta_hessian,
    cg_hessian_method = cg_control$hessian_method,
    vcov_method = vcov_controls$vcov_method,
    vcov_numderiv = vcov_controls$vcov_numderiv,
    rs_smooth_trust_radius = rs_control$smooth_trust_radius
  )

  fit_data <- data_fn(
    dataset = dataset,
    time_var = time_var,
    subject_var = subject_var,
    mu.formula = mu.formula,
    sigma.formula = sigma.formula,
    nu.formula = nu.formula,
    tau.formula = tau.formula,
    theta.formula = theta.formula,
    zeta.formula = zeta.formula,
    missingness = missingness,
    verbose = verbose
  )

  margin_dist <- .gl_bind_margin_args(margin_dist, fit_data$dataset)

  matrix_bundle <- matrix_fn(
    formulas_int = fit_data$formulas_int,
    margin_dist = margin_dist,
    copula_dist = copula_dist,
    dataset = fit_data$dataset
  )

  warm_start_args <- list(
    start_from = start_from,
    method = fit_controls$method,
    include_dlcopdpar = include_dlcopdpar,
    warm_start_joint = fit_controls$warm_start_joint,
    warm_start_joint_iter = fit_controls$warm_start_joint_iter,
    user_supplied_start = fit_controls$user_supplied_start,
    dataset_original = fit_data$dataset_original,
    margin_dist = margin_dist,
    copula_dist = copula_dist,
    time_var = time_var,
    subject_var = subject_var,
    missingness = missingness,
    mu.formula = mu.formula,
    sigma.formula = sigma.formula,
    nu.formula = nu.formula,
    tau.formula = tau.formula,
    theta.formula = theta.formula,
    zeta.formula = zeta.formula,
    optimizer_control = optimizer_control,
    true_val = true_val,
    vcov_method = fit_controls$vcov_method,
    vcov_numderiv = fit_controls$vcov_numderiv,
    use_Rcpp = use_Rcpp,
    lambda_start = lambda_start,
    lambda_penalty_K = lambda_penalty_K,
    verbose = verbose
  )
  warm_start_formals <- names(formals(warm_start_fn))
  if ("..." %in% warm_start_formals || "target_mm" %in% warm_start_formals) {
    warm_start_args$target_mm <- matrix_bundle$mm
  }
  warm_start <- do.call(warm_start_fn, warm_start_args)

  step_controls <- step_control_fn(
    method = fit_controls$method,
    include_dlcopdpar = include_dlcopdpar,
    start_step_size = rs_control$start_step_size,
    max_steps = rs_control$max_steps,
    step_adjustment = rs_control$step_adjustment,
    verbose = verbose
  )

  optimizer_context <- optimizer_context_fn(
    start_from = warm_start$start_from,
    warm_start_par_s = warm_start$warm_start_par_s,
    warm_start_lambda_s = warm_start$warm_start_lambda_s,
    mm = matrix_bundle$mm,
    margin_dist = margin_dist,
    copula_dist = copula_dist,
    dataset = fit_data$dataset,
    lambda_start = lambda_start,
    start_step_size = step_controls$start_step_size,
    copula_link = matrix_bundle$copula_link,
    inner_stop_crit = rs_control$inner_tol,
    outer_stop_crit = shared_control$outer_tol,
    cg_grad_tol = cg_control$grad_tol,
    cg_step_tol = cg_control$step_tol,
    method = fit_controls$method,
    verbose = verbose
  )

  projected_fallback <- warm_start$projected_start_fallback
  mapped_fallback <- warm_start$mapped_start_fallback
  alternate_start <- if (!is.null(mapped_fallback)) {
    mapped_fallback
  } else {
    projected_fallback
  }
  if (!is.null(alternate_start) &&
      !initial_likelihood_valid_fn(
        optimizer_context = optimizer_context,
        mm = matrix_bundle$mm,
        margin_dist = margin_dist,
        copula_dist = copula_dist,
        dataset = fit_data$dataset
    )) {
    projected_context <- optimizer_context_fn(
      start_from = alternate_start$par,
      warm_start_par_s = alternate_start$par_s,
      warm_start_lambda_s = alternate_start$lambda_s,
      mm = matrix_bundle$mm,
      margin_dist = margin_dist,
      copula_dist = copula_dist,
      dataset = fit_data$dataset,
      lambda_start = lambda_start,
      start_step_size = step_controls$start_step_size,
      copula_link = matrix_bundle$copula_link,
      inner_stop_crit = rs_control$inner_tol,
      outer_stop_crit = shared_control$outer_tol,
      cg_grad_tol = cg_control$grad_tol,
      cg_step_tol = cg_control$step_tol,
      method = method,
      verbose = verbose
    )
    if (initial_likelihood_valid_fn(
      optimizer_context = projected_context,
      mm = matrix_bundle$mm,
      margin_dist = margin_dist,
      copula_dist = copula_dist,
      dataset = fit_data$dataset
    )) {
      optimizer_context <- projected_context
      warm_start$start_from <- alternate_start$par
      warm_start$warm_start_par_s <- alternate_start$par_s
      warm_start$warm_start_lambda_s <- alternate_start$lambda_s
      warm_start$warm_start_info$projection_used <- is.null(mapped_fallback)
    }
  }

  list(
    margin_dist = margin_dist,
    controls = fit_controls,
    optimizer_control_requested = optimizer_control,
    optimizer_control_effective = .gl_effective_optimizer_control(
      optimizer_control,
      optimizer_context
    ),
    fit_data = fit_data,
    matrix_bundle = matrix_bundle,
    warm_start = warm_start,
    step_controls = step_controls,
    optimizer_context = optimizer_context
  )
}
