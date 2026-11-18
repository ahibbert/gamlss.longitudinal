#' Simulate longitudinal GAMLSS-copula data
#'
#' Simulates a balanced longitudinal dataset from a GAMLSS marginal distribution
#' and a first-order copula dependence model between consecutive time points.
#' Parameter specifications may be constants, time-indexed vectors,
#' subject-by-time matrices, or functions of the generated design data.
#'
#' @param n Number of subjects.
#' @param times Vector of observed time values.
#' @param margin_dist A `gamlss.dist` family object, for example `NO()` or
#'   `GA(mu.link = "log", sigma.link = "log")`.
#' @param copula_dist Copula family code. Supported codes are `"N"`, `"C"`,
#'   `"F"`, `"G"`, `"J"`, and `"t"`.
#' @param margin_params Named list of marginal distribution parameter
#'   specifications. Names should match the quantile function arguments, such
#'   as `mu`, `sigma`, `nu`, and `tau`, including `bd` for binomial trials.
#'   For a discrete margin, callbacks are evaluated once on the design before
#'   any responses are drawn.
#' @param copula_params Named list of copula parameter specifications. Use
#'   `theta` or `par` for the primary copula parameter, `tau` to specify
#'   Kendall's tau instead, and `zeta` or `par2` for the t-copula degrees of
#'   freedom.
#' @param covariates Optional data frame or function returning covariates. A
#'   data frame may have either `n` rows for subject-level covariates or
#'   `n * length(times)` rows for long-format covariates. A function is called
#'   with the base long-format data. It may return only new covariate columns
#'   or the input data with new columns added; columns already present in the
#'   base design are ignored.
#' @param seed Optional random seed.
#' @param subject_var,time_var,response_var Column names for the subject, time,
#'   and response variables.
#' @param include_truth If `TRUE`, include simulated uniforms and true
#'   parameter columns.
#' @param u_bounds Optional length-two numeric vector giving lower and upper
#'   bounds used to clamp simulated uniforms before applying the marginal
#'   quantile function. The default `NULL` leaves uniforms unchanged.
#'
#' @details With a discrete margin, each next count is drawn from the
#'   adjacent copula rectangle probabilities divided by the current count's
#'   marginal mass. Each transition uses a fresh uniform draw within the
#'   current count's CDF interval. With `u_bounds`, the clamp is applied to
#'   the reported probability and count after generating the underlying path.
#'   When `include_truth` is true, `u` is the probability supplied to the
#'   marginal quantile for the reported response.
#'
#' @return A long-format data frame with one row per subject-time observation.
#' @export
simulate_longitudinal_dataset <- function(
    n = 100,
    times = seq_len(3),
    margin_dist = gamlss.dist::NO(),
    copula_dist = "N",
    margin_params = list(mu = 0, sigma = 1),
    copula_params = list(theta = 0),
    covariates = NULL,
    seed = NULL,
    subject_var = "subject",
    time_var = "time",
    response_var = "response",
    include_truth = TRUE,
    u_bounds = NULL) {
  restore_seed <- .sim_set_seed_restore(seed)
  on.exit(restore_seed(), add = TRUE)

  n <- .sim_check_count(n, "n")
  if (length(times) < 2L) {
    stop("times must contain at least two values.", call. = FALSE)
  }
  if (anyDuplicated(times)) {
    stop("times must not contain duplicate values.", call. = FALSE)
  }
  .gl_validate_capability_route(
    margin_dist = margin_dist,
    copula_dist = copula_dist,
    context = "simulation"
  )
  copula_dist <- .copula_family_code(copula_dist)

  long <- .sim_base_long_data(n, times, subject_var, time_var)
  long <- .sim_add_covariates(long, covariates, n, length(times), subject_var)

  discrete <- .is_discrete_margin(margin_dist)
  family <- margin_dist$family[1L]
  evaluate_margin_args <- function(data) {
    params <- lapply(names(margin_params), function(parameter) {
      .sim_eval_long_param(
        margin_params[[parameter]], data, n, length(times),
        paste0("margin_params$", parameter)
      )
    })
    stats::setNames(params, names(margin_params))
  }
  qargs <- if (discrete) evaluate_margin_args(long) else list()

  u_mat <- .sim_copula_uniform_matrix(
    n = n,
    n_time = length(times),
    family = copula_dist,
    copula_params = copula_params,
    long_data = long,
    subject_var = subject_var,
    time_var = time_var,
    discrete_margin = if (discrete) family else NULL,
    margin_args = qargs
  )

  long$.sim_u <- .sim_apply_u_bounds(as.vector(t(u_mat)), u_bounds)
  if (!discrete) qargs <- evaluate_margin_args(long)
  long[[response_var]] <- .gl_simulation_quantiles(
    family = margin_dist$family[1L],
    p = long$.sim_u,
    args = qargs
  )

  if (include_truth) {
    long$u <- long$.sim_u
    for (param_name in names(margin_params)) {
      long[[paste0("true_", param_name)]] <- if (discrete) qargs[[param_name]] else .sim_eval_long_param(
        margin_params[[param_name]],
        long,
        n = n,
        n_time = length(times),
        label = paste0("margin_params$", param_name)
      )
    }
    edge_truth <- attr(u_mat, "edge_truth")
    if (!is.null(edge_truth)) {
      long$true_theta <- NA_real_
      long$true_zeta <- NA_real_
      edge_rows <- long$.sim_time_index > 1L
      long$true_theta[edge_rows] <- as.vector(t(edge_truth$theta))
      long$true_zeta[edge_rows] <- as.vector(t(edge_truth$zeta))
    }
  }

  long$.sim_u <- NULL
  long$.sim_time_index <- NULL
  long$.sim_subject_index <- NULL
  front_cols <- c(subject_var, time_var, response_var)
  long <- long[c(front_cols, setdiff(names(long), front_cols))]
  rownames(long) <- NULL
  long
}
