#' Fit a longitudinal joint regression model
#'
#' This function fits a longitudinal model to a dataset with gamlss margins
#' and a first-order copula dependence structure. Any linear, factor, or smooth
#' covariates can be included in the model formulas for either the margin parameters
#' (up to four parameters depending on the family, mu, sigma, nu, tau) 
#' or the copula parameters theta (up to two depending on the copula, theta, zeta).
#' 
#' Formulas are specified as for standard gamlss models, with the response variable on the left-hand side
#' of the formula for the mu parameter and the right-hand side specifying the predictors for each parameter.
#' The sigma, nu, tau, theta, and zeta formulas may be specified with a left-hand side of `~` to indicate no response variable.
#' e.g. `mu.formula = response ~ x1 + s(x2)`, `sigma.formula = ~ x3`, `theta.formula = ~ time`.
#' 
#' The marginal distribution is specified by a registered `gamlss.dist` family
#' object and the copula by a registered character code. Use
#' [longitudinal_capabilities()] to inspect the tested family/copula routes and
#' their response-domain requirements. The copula dependence structure is first-order, with
#' adjacent copula pairs linking the margins at each time point. The copula parameters
#' are shared across all adjacent pairs but may vary with covariates including time.
#' 
#' The user must specify both the time variable and the subject identifier as these are required for the model structure. 
#' The time variable is used to order the margins and adjacent copula pairs, 
#' and the subject identifier is used to link the repeated measurements for each subject. 
#' The response variable must be included in the `mu.formula` as the left-hand side.
#' 
#' The model is fit using Rigby and Stasinopoulos (RS) optimization against
#' the full joint likelihood by default, which iteratively updates the margin 
#' and copula parameters in turn until convergence. We provide alternative methods, 
#' including a separately optimized RS method which optimizes margin and copula
#' likelihood separately which is faster but slightly less optimal for overall fit,
#' and the Cole and Green method (CG) which updates based on the full first 
#' and second derivative matrix updating all parameters at once which is generally the 
#' slowest method but may provide better fits for some models. 
#'
#' @param data Long-format data frame containing the response, subject, time,
#'   and covariate columns.
#' @param margin_dist One registered `gamlss.dist` family object. A fit uses one
#'   homogeneous margin family; see [longitudinal_capabilities()].
#' @param copula_dist Copula distribution code for a registered route; see
#'   [longitudinal_capabilities()].
#' @param time_var Name of the time variable in `dataset`.
#' @param subject_var Name of the subject identifier in `dataset`.
#' @param missingness Handling of response patterns outside complete panels and
#'   observed prefixes ending in terminal dropout. The default, `"auto"`, uses
#'   the ordinary likelihood where possible and automatically switches to the
#'   segmented likelihood when leading unobserved visits or intermittent gaps
#'   are detected. Use `"error"` to stop instead, or `"segment"` to explicitly
#'   permit segmentation. With a segmented objective, the likelihood is
#'   evaluated within contiguous observed
#'   segments and different segments are treated as independent. The first
#'   observed visit after leading missing visits explicitly starts a segment.
#'   Terminal dropout does not require segmentation. AIC and
#'   BIC remain available for segmented fits under the stated independence
#'   assumption. Standard inference methods automatically use subject-clustered
#'   sandwich covariance instead of ordinary model-Hessian inference for this
#'   composite objective; model-based standard errors remain available for
#'   comparison, and bootstrap sensitivity analysis remains advisable.
#'   Numerical integration across gaps and a distinct delayed-entry
#'   likelihood are not yet implemented.
#' @param mu.formula Formula for the mu parameter of the marginal distribution
#' @param sigma.formula Formula for the sigma parameter of the marginal distribution
#' @param nu.formula Formula for the nu parameter of the marginal distribution
#' @param tau.formula Formula for the tau parameter of the marginal distribution
#' @param theta.formula Formula for the theta parameter of the copula distribution
#' @param zeta.formula Formula for the zeta parameter of the copula distribution
#' @param method Optimization method: `"RS"` for Rigby-Stasinopoulos
#'   backfitting (the default) or `"CG"` for Cole-Green updates.
#' @param optimization Parameter-update strategy. `"joint"` includes the
#'   copula contribution when updating marginal parameters and is the default.
#'   `"separate"` uses separate margin and copula updates and is available only
#'   for `method = "RS"`.
#' @param optimizer_control Optimizer settings created by
#'   [gamlss_longitudinal_control()]. `method` remains a top-level choice.
#' @param start Optional named starting coefficients, or a previously fitted
#'   `gamlss.longitudinal` model. A fitted-model start maps every retained fixed
#'   coefficient by name and also carries compatible smooth coefficients and
#'   smoothing penalties, which is useful for nested-model selection. `NULL`
#'   uses automatic starting values.
#' @param compute_vcov Logical; if `TRUE` (default), compute and store the
#' model variance-covariance output at the end of fitting.
#' @param vcov_method Character; fit-time vcov method when `compute_vcov = TRUE`.
#' One of `"analytical"` or `"numderiv"`. Analytical vcov covers continuous
#' margins and supported exact-discrete rectangle likelihoods, and falls back
#' to the numerical reference path if the analytical Hessian cannot be inverted.
#' @param verbose Logical or integer verbosity level. The default, `1`, reports
#'   fitting stages and outer iterations. Use `0` or `FALSE` for no progress
#'   output; levels `2` and `3` add increasingly detailed optimizer tracing.
#' @param ... Deprecated argument aliases retained for compatibility. New code
#'   should use `data`, `optimization`, `start`, and `optimizer_control`.
#' @param margin_args Named list of fixed marginal-distribution arguments.
#'   Binomial families require `list(bd = ...)`, where `bd` is either one
#'   positive integer used for every row or the name of a data column containing
#'   row-specific trial counts. Other fixed family arguments are not currently
#'   supported.
#'
#' @details
#' The ordered unique values of `time_var` define the scheduled visit order.
#' Adjacent copula pairs link globally adjacent scheduled visits, rather than
#' the next observed record for a subject. Numeric spacing is not used as an
#' elapsed-time distance, so irregularly spaced values affect ordering only.
#' Duplicate subject-time rows are rejected. Pair-level dependence covariates
#' and parameters are evaluated at the earlier (left) scheduled row. With a
#' segmented objective, pairs crossing an intermittent gap are omitted and the
#' resulting observed segments contribute independently.
#'
#' Legacy aliases, including `dataset`, `include_dlcopdpar`, `start_from`, and
#' flat optimizer controls, are accepted temporarily through `...`. Supplying
#' both an alias and its replacement is an error. Inactive optimizer sections
#' in a new control object are also rejected.
#'
#' @export
gamlss_longitudinal <- function(data = NULL,
                                margin_dist,
                                copula_dist,
                                time_var = NA,
                                subject_var = NA,
                                missingness = c("auto", "segment", "error"),
                                mu.formula = response ~ 1,
                                sigma.formula = ~ 1,
                                nu.formula = ~ 1,
                                tau.formula = ~ 1,
                                theta.formula = ~ 1,
                                zeta.formula = ~ 1,
                                method = "RS",
                                optimization = c("joint", "separate"),
                                optimizer_control = NULL,
                                start = NULL,
                                compute_vcov = TRUE,
                                vcov_method = c("analytical", "numderiv"),
                                margin_args = list(),
                                verbose = 1,
                                ...) {
  data_supplied <- !missing(data)
  optimization_supplied <- !missing(optimization)
  start_supplied <- !missing(start)
  legacy <- .gl_resolve_legacy_fit_dots(list(...))
  dataset <- .gl_resolve_fit_data_alias(data, data_supplied, legacy)
  start_from <- .gl_resolve_fit_start_alias(start, start_supplied, legacy)
  verbose <- .gl_normalize_verbose(verbose)
  missingness <- match.arg(missingness)
  if (!missing(margin_args) || is.null(attr(margin_dist, "margin_args_spec", exact = TRUE))) {
    margin_dist <- .gl_attach_margin_args_spec(margin_dist, margin_args)
  }
  optimizer_legacy_names <- unique(c(
    names(.gl_legacy_optimizer_map()),
    "use_backtracking", "backtracking_max_halves"
  ))
  legacy_values <- legacy$values[optimizer_legacy_names]
  legacy_supplied <- legacy$supplied[optimizer_legacy_names]
  optimizer_control <- .gl_resolve_optimizer_control(
    method = method,
    optimizer_control = optimizer_control,
    legacy_values = legacy_values,
    legacy_supplied = legacy_supplied
  )
  method <- attr(optimizer_control, "method")
  optimization <- .gl_resolve_optimization(
    method = method,
    optimization = optimization,
    optimization_supplied = optimization_supplied,
    legacy = legacy
  )
  optimizer_control <- .gl_apply_route_optimizer_defaults(
    control = optimizer_control,
    method = method,
    margin_dist = margin_dist,
    copula_dist = copula_dist,
    optimization = optimization,
    purpose = "fit"
  )
  include_dlcopdpar <- identical(optimization, "joint")

  fit <- .gl_run_gamlss_longitudinal_entrypoint(
    dataset = dataset,
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
    include_dlcopdpar = include_dlcopdpar,
    check_dlcopdpar_gradient = legacy$values$check_dlcopdpar_gradient,
    start_from = start_from,
    verbose = verbose,
    plot_results = legacy$values$plot_results,
    true_val = legacy$values$true_val,
    method = method,
    optimizer_control = optimizer_control,
    compute_vcov = compute_vcov,
    vcov_method = vcov_method,
    vcov_numderiv = legacy$values$vcov_numderiv,
    use_Rcpp = legacy$values$use_Rcpp,
    lambda_start = optimizer_control$shared$lambda_start,
    lambda_penalty_K = optimizer_control$shared$lambda_penalty
  )
  fit$optimization <- optimization
  fit$convergence$optimization <- optimization
  fit
}
gamlss.longitudinal <- gamlss_longitudinal
