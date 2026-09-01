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
  attr(warm_control, "specified") <- list(
    shared = c(outer_tol = FALSE, max_outer_iter = TRUE, max_elapsed_sec = FALSE, stop_on_convergence = FALSE),
    rs = c("warm_start_joint", "warm_start_joint_iter"),
    cg = character(0)
  )

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
              start_from = NA,
              verbose = 0,
              plot_results = FALSE,
              true_val = true_val,
              method = method,
              optimizer_control = warm_control,
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
