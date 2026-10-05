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
    missingness = "error",
    mu.formula,
    sigma.formula,
    nu.formula,
    tau.formula,
    theta.formula,
    zeta.formula,
    optimizer_control,
    true_val,
    method,
    warm_start_joint_iter,
    vcov_method,
    vcov_numderiv,
    use_Rcpp,
    lambda_start,
    lambda_penalty_K) {
  warm_fit <- NULL
  warm_output <- NULL
  warm_warnings <- character(0)
  warm_err <- NULL

  warm_control <- optimizer_control
  warm_control$shared$max_outer_iter <- as.integer(warm_start_joint_iter)
  warm_control$rs$warm_start_joint <- FALSE
  warm_control$rs$warm_start_joint_iter <- 0L
  # The separate fit is only an initializer for the full joint optimizer.
  # Preserve its bounded inner iterations, backtracking, adaptive step, and
  # lambda policy, but do not pay for the endpoint-only joint BFGS/refinement
  # sequence twice. The joint fit retains the user's full refinement controls.
  warm_control$rs$bfgs_refine <- FALSE
  specified <- attr(warm_control, "specified")
  specified$shared[] <- FALSE
  specified$shared[["max_outer_iter"]] <- TRUE
  specified$rs <- c("warm_start_joint", "warm_start_joint_iter", "bfgs_refine")
  specified$cg <- character(0)
  attr(warm_control, "specified") <- specified

  tryCatch(
    {
      warm_output <- capture.output(
        withCallingHandlers(
          {
            warm_fit <- fit_fn(
              data = dataset_original,
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
              optimization = "separate",
              start = NA,
              verbose = 0,
              method = method,
              optimizer_control = warm_control,
              compute_vcov = FALSE,
              vcov_method = vcov_method
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
