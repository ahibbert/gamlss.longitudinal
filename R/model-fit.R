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
#' The marginal distribution is specified by a gamlss family object, e.g. GA(), NO(), PO(), NBI(), etc.,
#' while the copula distribution is specified by a character code for one of the implemented
#' copula families: Gaussian ("N"), Clayton ("C"), Frank ("F"), Gumbel ("G"), Joe ("J"),
#' or Student's t ("t"). The copula dependence structure is first-order, with
#' adjacent copula pairs linking the margins at each time point. The copula parameters
#' are shared across all adjacent pairs but may vary with covariates including time.
#' 
#' The user must specify both the time variable and the subject identifier as these are required for the model structure. 
#' The time variable is used to order the margins and adjacent copula pairs, 
#' and the subject identifier is used to link the repeated measurements for each subject. 
#' The response variable must be included in the `mu.formula` as the left-hand side.
#' 
#' The model is fit using Rigby and Stasinopoulis (RS) optimisation against
#' the full joint likelihood by default, which iteratively updates the margin 
#' and copula parameters in turn until convergence. We provide alternative methods, 
#' including a separately optimised RS method which optimises margin and copula
#' likelihood separately which is faster but slightly less optimal for overall fit,
#' and the Cole and Green method (CG) which updates based on the full first 
#' and second derivative matrix updating all parameters at once which is generally the 
#' slowest method but may provide better fits for some models. 
#'
#' @param dataset Long-format data frame containing the response, subject, time,
#'   and covariate columns.
#' @param margin_dist Marginal distribution specified as a gamlss family object,
#' e.g. GA(), NO(), PO(), NBI(), etc.
#' @param copula_dist Copula distribution code, one of "N", "C", "F", "G", "J", or "t".
#' @param time_var Name of the time variable in `dataset`.
#' @param subject_var Name of the subject identifier in `dataset`.
#' @param mu.formula Formula for the mu parameter of the marginal distribution
#' @param sigma.formula Formula for the sigma parameter of the marginal distribution
#' @param nu.formula Formula for the nu parameter of the marginal distribution
#' @param tau.formula Formula for the tau parameter of the marginal distribution
#' @param theta.formula Formula for the theta parameter of the copula distribution
#' @param zeta.formula Formula for the zeta parameter of the copula distribution
#' @param include_dlcopdpar Include the derivative of the copula likelihood with respect
#' to the margin parameters in the joint likelihood. This is only relevant for method=`RS` fits
#' and is `TRUE` by default. Setting this to `FALSE` can speed up fitting for most models
#' due to separately optimised margin and copula fits at the expense of a slightly less optimal overall fit. 
#' For some models, including the dlcopdpar contribution can be important for good convergence and fit.
#' @param check_dlcopdpar_gradient If `TRUE`, run an optional finite-difference
#' diagnostic for the margin score contribution when `include_dlcopdpar = TRUE`.
#' @param inner_stop_crit Stopping criterion for the inner loop. If `NA` or
#' `NULL`, an automatic data-adaptive value is used. This is based on change in log
#' likelihood in each iteration so `0.1` results in algorithm stopping in the inner loop
#' if likelihood changes by less than 0.1 in that iteration. 
#' Setting this to `0` forces the inner loop to run for the full `max_inner_iter` iterations.
#' Note inner iterations are only relevant for the RS algorithm as CG does not have an inner loop.
#' @param outer_stop_crit Stopping criterion for the outer loop. If `NA` or
#' `NULL`, an automatic data-adaptive value is used. This is based on change 
#' in log likelihood in each iteration so `0.1` results in algorithm stopping 
#' if likelihood changes by less than 0.1 in that outer iteration. 
#' Setting this to `0` forces the outer loop to run for the full `max_outer_iter` iterations.
#' @param start_step_size Initial step size for the backfitting algorithm
#' @param step_adjustment Step size adjustment factor
#' @param max_steps Maximum number of times for reducing the step size
#' @param start_from Starting values for the parameters if needed
#' @param warm_start_joint Logical; if `TRUE` (default), RS joint fits started
#' without explicit `start_from` first run a short separate RS stabilisation
#' phase and use those coefficients as the joint starting values.
#' @param warm_start_joint_iter Integer; number of separate RS outer iterations
#' used for the default joint warm start.
#' @param verbose Level of output to the console 3 = ALL, 0 = Minimal
#' @param plot_results Plot the results of the optimisation (depreciated) 
#' @param true_val True values for the parameters if known for plotting (depreciated)
#' @param method Optimisation method to use, 'RS' for Rigby and Stasinopoulis backfitting (default) or 'CG' for Cole and Green.
#' @param max_outer_iter Maximum number of outer iterations for the optimisation algorithm (both RS and CG).
#' @param max_inner_iter Maximum number of inner iterations for the RS backfitting algorithm. 
#' Not relevant for CG which does not have an inner loop.
#' @param max_negative_outer_streak Maximum number of consecutive negative outer
#' log-likelihood changes allowed before stopping. 
#' (This is a safeguard against very badly formed likelihood steps which can cause the algorithm to diverge and produce NaNs. Setting this to `Inf` disables this stopping criterion.)
#' @param max_elapsed_sec Optional maximum elapsed fitting time in seconds.
#' If finite, the optimiser stops with an error once this budget is exceeded.
#' @param use_backtracking Logical; if `TRUE` (default), apply step-halving
#' backtracking to reject downhill inner updates.
#' @param backtracking_max_halves Integer; maximum number of consecutive
#' step halvings attempted after a rejected update before taking no step.
#' @param cg_max_stall Integer; for `method = "CG"` only. Maximum number of
#' consecutive outer iterations where no improving step is found before CG stops.
#' @param cg_max_delta Numeric; for `method = "CG"` only. Maximum absolute
#' coefficient step size used to limit Newton/trust-region updates.
#' @param cg_armijo_c1 Numeric; for `method = "CG"` only. Minimum improvement
#' threshold used by the line-search acceptance rule.
#' @param cg_grad_tol Numeric; for `method = "CG"` only. Penalized-gradient
#' infinity-norm convergence tolerance. If `NA`, selected from `outer_stop_crit`.
#' @param cg_step_tol Numeric; for `method = "CG"` only. Accepted-step L2
#' convergence tolerance. If `NA`, selected from `outer_stop_crit`.
#' @param cg_update_lambda Logical; for `method = "CG"` only. If `TRUE`, update
#' smoother penalties during CG iterations.
#' @param cg_lambda_update_every Integer; for `method = "CG"` only. When
#' `cg_update_lambda = TRUE`, update each smoother's lambda every this many
#' outer iterations. Use `1` to update every CG iteration.
#' @param cg_max_lambda_updates Integer; for `method = "CG"` only. Maximum
#' number of smoother penalty update rounds. Use `NA` for no cap.
#' @param cg_raw_loglik_drop_tol Numeric; for `method = "CG"` only. Stop CG as
#' not converged if the raw joint log-likelihood drops this far below the best
#' raw joint log-likelihood seen after at least one lambda update. Use `NA` to
#' disable.
#' @param cg_line_search Character; for `method = "CG"` only. `"best"` evaluates
#' candidate steps up to `cg_max_line_search_evals` before taking the largest
#' improvement, while `"first"` accepts the first improving candidate step.
#' @param cg_max_line_search_evals Integer; for `method = "CG"` only. Optional
#' cap on the number of candidate likelihood evaluations per outer iteration.
#' @param cg_gradient_method Character; for `method = "CG"` only.
#' `"analytical"` uses the same score components as RS, `"forward"` uses
#' one-sided finite differences, and `"central"` uses two-sided finite
#' differences.
#' @param discrete_score_method Character. For discrete margins using exact
#' rectangle likelihoods, choose `"analytical"` for vectorised rectangle-score
#' assembly or `"finite"` for slow row-wise finite-difference scores.
#' @param cg_zeta_hessian Character; for `method = "CG"` only. `"analytical"`
#' uses the analytical Hessian for the zeta block, while `"finite"` replaces
#' the zeta-zeta block with central finite differences of the raw joint
#' log-likelihood.
#' @param cg_hessian_method Character; for `method = "CG"` only. `"analytical"`
#' uses the semi-analytical Hessian for Newton steps, `"finite"` uses a full
#' finite-difference Hessian, and `"auto"` tries analytical then falls back to
#' finite differences when needed.
#' @param compute_vcov Logical; if `TRUE` (default), compute and store the
#' model variance-covariance output at the end of fitting.
#' @param vcov_method Character; fit-time vcov method when `compute_vcov = TRUE`.
#' One of `"analytical"` or `"numderiv"`. Analytical vcov covers continuous
#' margins and supported exact-discrete rectangle likelihoods, and falls back
#' to the numerical reference path if the analytical Hessian cannot be inverted.
#' @param vcov_numderiv Logical; passed to `vcov.gamlss.longitudinal()` when
#' `compute_vcov = TRUE`.
#' @param use_Rcpp Use Rcpp for matrix operations (depreciated)
#' @param lambda_start Optional starting value for smooth-term penalties.
#' @param lambda_penalty_K Penalty strength used when updating smooth-term
#'   smoothing parameters.
#' @param rs_update_lambda Logical; for `method = "RS"` only. If `TRUE`,
#' update smoothing parameters by the RS GAIC step; if `FALSE`, keep
#' `lambda_start` fixed.
#' @param rs_smooth_trust_radius Numeric; for `method = "RS"` only. Optional
#' L2 trust radius applied separately to each smooth coefficient block after
#' the RS weighted least-squares proposal. Use `Inf` to disable.
#'
#' @export
gamlss_longitudinal <- function(dataset,
                                margin_dist,
                                copula_dist,
                                time_var = NA,
                                subject_var = NA,
                                mu.formula = ("response ~ 1"),
                                sigma.formula = ("~ 1"),
                                nu.formula = ("~ 1"),
                                tau.formula = ("~ 1"),
                                theta.formula = ("~ 1"),
                                zeta.formula = ("~ 1"),
                                include_dlcopdpar = TRUE,
                                check_dlcopdpar_gradient = FALSE,
                                inner_stop_crit = NA,
                                outer_stop_crit = NA,
                                start_step_size = .5,
                                step_adjustment = NA,
                                max_steps = 5,
                                start_from = NA,
                                warm_start_joint = TRUE,
                                warm_start_joint_iter = 5,
                                verbose = 1,
                                plot_results = FALSE,
                                true_val = NA,
                                method = "RS",
                                max_outer_iter = 100,
                                max_inner_iter = 100,
                                max_negative_outer_streak = 10,
                                max_elapsed_sec = Inf,
                                use_backtracking = TRUE,
                                backtracking_max_halves = 50,
                                cg_max_stall = 5,
                                cg_max_delta = 0.5,
                                cg_armijo_c1 = 1e-4,
                                cg_grad_tol = NA,
                                cg_step_tol = NA,
                                cg_update_lambda = TRUE,
                                cg_lambda_update_every = 10,
                                cg_max_lambda_updates = NA,
                                cg_raw_loglik_drop_tol = 10,
                                cg_line_search = "best",
                                cg_max_line_search_evals = 60,
                                cg_gradient_method = "forward",
                                discrete_score_method = c("analytical", "finite"),
                                cg_zeta_hessian = "analytical",
                                cg_hessian_method = c("analytical", "finite", "auto"),
                                compute_vcov = TRUE,
                                vcov_method = c("analytical", "numderiv"),
                                vcov_numderiv = FALSE,
                                use_Rcpp = FALSE,
                                lambda_start = NA,
                                lambda_penalty_K = 2,
                                rs_update_lambda = TRUE,
                                rs_smooth_trust_radius = Inf) {
  .gl_run_gamlss_longitudinal_entrypoint(
    dataset = dataset,
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
    include_dlcopdpar = include_dlcopdpar,
    check_dlcopdpar_gradient = check_dlcopdpar_gradient,
    inner_stop_crit = inner_stop_crit,
    outer_stop_crit = outer_stop_crit,
    start_step_size = start_step_size,
    step_adjustment = step_adjustment,
    max_steps = max_steps,
    start_from = start_from,
    warm_start_joint = warm_start_joint,
    warm_start_joint_iter = warm_start_joint_iter,
    verbose = verbose,
    plot_results = plot_results,
    true_val = true_val,
    method = method,
    max_outer_iter = max_outer_iter,
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
    cg_max_lambda_updates = cg_max_lambda_updates,
    cg_raw_loglik_drop_tol = cg_raw_loglik_drop_tol,
    cg_line_search = cg_line_search,
    cg_max_line_search_evals = cg_max_line_search_evals,
    cg_gradient_method = cg_gradient_method,
    discrete_score_method = discrete_score_method,
    cg_zeta_hessian = cg_zeta_hessian,
    cg_hessian_method = cg_hessian_method,
    compute_vcov = compute_vcov,
    vcov_method = vcov_method,
    vcov_numderiv = vcov_numderiv,
    use_Rcpp = use_Rcpp,
    lambda_start = lambda_start,
    lambda_penalty_K = lambda_penalty_K,
    rs_update_lambda = rs_update_lambda,
    rs_smooth_trust_radius = rs_smooth_trust_radius
  )
}
gamlss.longitudinal <- gamlss_longitudinal
