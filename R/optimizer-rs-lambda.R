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
      if (verbose > 1) {
        cat("\nOptimizing smoothing parameter for ", par_name, " - ", smooth_name, "\n", sep = "")
      }

      invalid_gaic <- 1e100
      optim_lambda <- function(lambda_val) {
        lambda_s_temp <- lambda_s
        lambda_s_temp[[par_name]][[smooth_name]] <- lambda_val
        backfitting_iteration_results <- tryCatch(
          backfitting_fn(
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
          ),
          error = function(error) {
            if (inherits(error, "gamlss.longitudinal_time_limit_error")) {
              stop(error)
            }
            NULL
          }
        )
        if (is.null(backfitting_iteration_results)) return(invalid_gaic)
        gaic <- as.numeric(backfitting_iteration_results$GAIC_lambda_k)
        if (length(gaic) != 1L || !is.finite(gaic)) invalid_gaic else gaic
      }

      optim_lambda_out <- tryCatch(
        optim_fn(
          par = lambda_s[[par_name]][[smooth_name]],
          fn = optim_lambda,
          method = "L-BFGS-B",
          lower = 0.01,
          upper = 1e6,
          control = list(factr = 1, pgtol = .1)
        ),
        error = function(error) {
          if (inherits(error, "gamlss.longitudinal_time_limit_error")) {
            stop(error)
          }
          NULL
        }
      )

      candidate_lambda <- optim_lambda_out$par %||% NA_real_
      candidate_value <- optim_lambda_out$value %||% NA_real_
      valid_candidate <- !is.null(optim_lambda_out) &&
        length(candidate_lambda) == 1L && is.finite(candidate_lambda) &&
        candidate_lambda >= 0.01 && candidate_lambda <= 1e6 &&
        (is.na(candidate_value) ||
          (length(candidate_value) == 1L && is.finite(candidate_value) &&
            candidate_value < invalid_gaic))
      if (valid_candidate) {
        lambda_s[[par_name]][[smooth_name]] <- candidate_lambda
      }
      if (verbose > 2) {
        print(paste("Chosen lambda:", round(lambda_s[[par_name]][[smooth_name]], 2), "| Penalty K =", K))
      }
    }
  }

  lambda_s
}
