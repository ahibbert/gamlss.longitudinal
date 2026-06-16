#' Collect core fitted-object fields into list
#'
#' @noRd
.gl_build_fit_object_core <- function(
    par_cov,
    log_lik_history,
    par_history,
    calc_lik_out_end,
    mm,
    margin_dist,
    copula_dist,
    include_dlcopdpar,
    dataset,
    dataset_original,
    response_var,
    time_var,
    subject_var,
    formulas,
    formulas_int,
    var_map,
    par_s,
    lambda_s,
    df_s,
    weights_final,
    method,
    warm_start_info,
    convergence_info) {
  return_list <- list(
    par_cov,
    log_lik_history,
    par_history,
    calc_lik_out_end,
    mm,
    margin_dist,
    copula_dist,
    include_dlcopdpar,
    dataset$response,
    dataset$time,
    dataset$subject,
    par_s,
    lambda_s,
    df_s,
    weights_final
  )
  names(return_list) <- c(
    "par",
    "log_lik_history",
    "par_history",
    "calc_lik_out_end",
    "model_matrix",
    "margin_dist",
    "copula_dist",
    "include_dlcopdpar",
    "response",
    "response_margin",
    "response_subject",
    "par_s",
    "lambda_s",
    "df_s",
    "weights"
  )

  return_list$dataset <- dataset
  return_list$dataset_original <- dataset_original
  return_list$response_var <- response_var
  return_list$time_var <- time_var
  return_list$subject_var <- subject_var
  return_list$formulas <- formulas
  return_list$formulas_int <- formulas_int
  return_list$var_map <- var_map
  return_list$optim_method <- method
  return_list$warm_start_joint <- warm_start_info
  return_list$convergence <- convergence_info

  return_list
}

#' Attach optimiser traces to the final fitted object
#' These can be plot later to review the optimisation path and diagnose convergence issues.
#'
#' @noRd
.gl_attach_fit_optimizer_traces <- function(
    return_list,
    method,
    rs_block_trace,
    cg_lambda_trace,
    cg_step_trace) {
  if (!identical(method, "CG")) {
    return_list$rs_block_trace <- if (length(rs_block_trace)) {
      do.call(rbind, rs_block_trace)
    } else {
      data.frame()
    }
  }

  if (identical(method, "CG")) {
    return_list$cg_lambda_trace <- cg_lambda_trace
    return_list$cg_step_trace <- if (length(cg_step_trace)) {
      do.call(rbind, cg_step_trace)
    } else {
      data.frame()
    }
  }

  return_list
}
