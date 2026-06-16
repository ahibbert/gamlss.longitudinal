#' Runs the separate fit used for joint RS warm start
#'
#' @noRd
.gl_call_joint_warm_start_fit <- function(
    fit_fn,
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
    inner_stop_crit,
    outer_stop_crit,
    start_step_size,
    step_adjustment,
    max_steps,
    true_val,
    method,
    warm_start_joint_iter,
    max_inner_iter,
    max_negative_outer_streak,
    max_elapsed_sec,
    use_backtracking,
    backtracking_max_halves,
    cg_max_stall,
    cg_max_delta,
    cg_armijo_c1,
    cg_grad_tol,
    cg_step_tol,
    cg_update_lambda,
    cg_lambda_update_every,
    cg_line_search,
    cg_max_line_search_evals,
    cg_gradient_method,
    discrete_score_method,
    cg_zeta_hessian,
    cg_hessian_method,
    vcov_method,
    vcov_numderiv,
    use_Rcpp,
    lambda_start,
    lambda_penalty_K) {
  warm_fit <- NULL
  warm_output <- NULL
  warm_warnings <- character(0)
  warm_err <- NULL

  tryCatch(
    {
      warm_output <- capture.output(
        withCallingHandlers(
          {
            warm_fit <- fit_fn(
              dataset = dataset_original,
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
              include_dlcopdpar = FALSE,
              check_dlcopdpar_gradient = FALSE,
              inner_stop_crit = inner_stop_crit,
              outer_stop_crit = outer_stop_crit,
              start_step_size = start_step_size,
              step_adjustment = step_adjustment,
              max_steps = max_steps,
              start_from = NA,
              warm_start_joint = FALSE,
              warm_start_joint_iter = 0L,
              verbose = 0,
              plot_results = FALSE,
              true_val = true_val,
              method = method,
              max_outer_iter = warm_start_joint_iter,
              max_inner_iter = max_inner_iter,
              max_negative_outer_streak = max_negative_outer_streak,
              max_elapsed_sec = max_elapsed_sec,
              use_backtracking = use_backtracking,
              backtracking_max_halves = backtracking_max_halves,
              cg_max_stall = cg_max_stall,
              cg_max_delta = cg_max_delta,
              cg_armijo_c1 = cg_armijo_c1,
              cg_grad_tol = cg_grad_tol,
              cg_step_tol = cg_step_tol,
              cg_update_lambda = cg_update_lambda,
              cg_lambda_update_every = cg_lambda_update_every,
              cg_line_search = cg_line_search,
              cg_max_line_search_evals = cg_max_line_search_evals,
              cg_gradient_method = cg_gradient_method,
              discrete_score_method = discrete_score_method,
              cg_zeta_hessian = cg_zeta_hessian,
              cg_hessian_method = cg_hessian_method,
              compute_vcov = FALSE,
              vcov_method = vcov_method,
              vcov_numderiv = vcov_numderiv,
              use_Rcpp = use_Rcpp,
              lambda_start = lambda_start,
              lambda_penalty_K = lambda_penalty_K
            )
          },
          warning = function(w) {
            warm_warnings <<- c(warm_warnings, conditionMessage(w))
            invokeRestart("muffleWarning")
          }
        ),
        type = "output"
      )
    },
    error = function(e) {
      warm_err <<- e
    }
  )

  list(
    fit = warm_fit,
    output = warm_output,
    warnings = warm_warnings,
    error = warm_err
  )
}
