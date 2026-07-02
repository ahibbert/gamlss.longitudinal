#' Update RS smoothing parameters for one parameter block
#'
#' Runs the existing GAIC-based one-dimensional smoothing-parameter search for
#' each smooth attached to the active parameter block.
#'
#' @noRd
.gl_rs_update_smoothing_parameters <- function(
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
    verbose,
    backfitting_fn,
    optim_fn = stats::optim) {
  for (smooth_name in names(lambda_s[[par_name]])) {
    if (isTRUE(rs_update_lambda) && inner_run_counter == 1) {
      cat(paste("\nOptimising smoothing parameter for", par_name, "-", smooth_name))

      optim_lambda <- function(lambda_val) {
        lambda_s_temp <- lambda_s
        lambda_s_temp[[par_name]][[smooth_name]] <- lambda_val
        backfitting_iteration_results <- backfitting_fn(
          par_s = par_s,
          par_cov = par_cov,
          beta_start = beta_start,
          lambda_s = lambda_s_temp,
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

        loglik <- backfitting_iteration_results$calc_lik_out_end$log_lik["joint"]
        df_total <- sum(unlist(backfitting_iteration_results$df_s[[par_name]]))
        gaic_val <- backfitting_iteration_results$GAIC_lambda_k

        backfitting_iteration_results$GAIC_lambda_k
      }

      optim_lambda_out <- optim_fn(
        par = lambda_s[[par_name]][[smooth_name]],
        fn = optim_lambda,
        method = "L-BFGS-B",
        lower = 0.01,
        upper = 1e6,
        control = list(factr = 1, pgtol = .1)
      )

      lambda_s[[par_name]][[smooth_name]] <- optim_lambda_out$par
      if (verbose > 2) {
        print(paste("Chosen lambda:", round(lambda_s[[par_name]][[smooth_name]], 2), "| Penalty K =", K))
      }
    }
  }

  lambda_s
}
