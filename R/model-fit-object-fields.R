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
    missingness_contract,
    formulas,
    formulas_int,
    var_map,
    par_s,
    lambda_s,
    df_s,
    weights_final,
    method,
    warm_start_info,
    convergence_info,
    optimizer_control_requested = NULL,
    optimizer_control_effective = NULL) {
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
  return_list$missingness <- missingness_contract
  return_list$likelihood_contract <- list(
    objective = missingness_contract$objective %||% "ordinary",
    between_segment_assumption = missingness_contract$between_segment_assumption %||% "not_applicable",
    criteria_status = missingness_contract$criteria_status %||% "available",
    statement = missingness_contract$statement %||% NA_character_,
    future_support = missingness_contract$future_support %||% NA_character_
  )
  return_list$formulas <- formulas
  return_list$formulas_int <- formulas_int
  return_list$var_map <- var_map
  return_list$optim_method <- method
  family <- .gl_capability_margin_code(margin_dist)
  capability <- .gl_capability_margin_spec(family)
  return_list$capability_registry_version <- .gl_capability_registry_version()
  return_list$capability_route <- list(
    registry_version = .gl_capability_registry_version(),
    margin_family = family,
    copula = copula_dist,
    family_type = capability$family_type %||% NA_character_,
    likelihood_route = capability$likelihood_route %||% NA_character_,
    diagnostics = capability$diagnostics %||% NA_character_
  )
  return_list$warm_start_joint <- warm_start_info
  return_list$convergence <- convergence_info
  return_list$optimizer_control_requested <- optimizer_control_requested
  return_list$optimizer_control_effective <- optimizer_control_effective

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
