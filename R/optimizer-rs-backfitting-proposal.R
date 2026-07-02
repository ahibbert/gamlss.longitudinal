#' Prepare the RS backfitting proposal for one parameter block
#'
#' @noRd
.gl_rs_prepare_backfitting_proposal <- function(
    lambda_s,
    par_name,
    par_s,
    par_cov,
    beta_start,
    K,
    margin_dist,
    copula_dist,
    dataset,
    mm,
    copula_link,
    df_s,
    step_size,
    rs_update_lambda,
    inner_run_counter,
    outer_only_run_counter,
    verbose,
    backfitting_fn,
    lambda_update_fn = .gl_rs_update_smoothing_parameters) {
  num_smooths <- length(lambda_s[[par_name]])
  initial_results <- NULL

  if (num_smooths > 0 && outer_only_run_counter != 1) {
    lambda_s <- lambda_update_fn(
      lambda_s = lambda_s,
      par_name = par_name,
      par_s = par_s,
      par_cov = par_cov,
      beta_start = beta_start,
      K = K,
      margin_dist = margin_dist,
      copula_dist = copula_dist,
      dataset = dataset,
      mm = mm,
      copula_link = copula_link,
      df_s = df_s,
      step_size = step_size,
      rs_update_lambda = rs_update_lambda,
      inner_run_counter = inner_run_counter,
      verbose = verbose,
      backfitting_fn = backfitting_fn
    )
  }

  proposed_results <- backfitting_fn(
    par_s = par_s,
    par_cov = par_cov,
    beta_start = beta_start,
    lambda_s = lambda_s,
    first_inner_run = FALSE,
    K = K,
    margin_dist = margin_dist,
    copula_dist = copula_dist,
    dataset = dataset,
    mm = mm,
    copula_link = copula_link,
    df_s = df_s,
    step_size = step_size,
    par_name = par_name
  )

  list(
    lambda_s = lambda_s,
    proposed_results = proposed_results,
    initial_results = initial_results,
    num_smooths = num_smooths
  )
}
