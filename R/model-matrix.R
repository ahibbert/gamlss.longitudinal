#' Create model matrices for model fitting
#'
#' This function takes the forumlas for each parameter mu,sigma,nu,tau,theta,zeta
#' and creates a list of model matrices mm with items mm$x and mm$s for
#' fixed and smooth terms respectively, with each of those lists being lists of each parameter
#' and their respective model matrices
#' 
#' For smooth term construction, the function looks for terms in the formula of the form s(x) 
#' and uses mgcv::smoothCon() to construct the model matrix for that term. 
#' 
#' Smooth model matrix is stored as mm$s and the fixed model matrix is stored as mm$x.
#' 
#' @param mu.formula Formula for the mu parameter of the marginal distribution
#' @param sigma.formula Formula for the sigma parameter of the marginal distribution
#' @param nu.formula Formula for the nu parameter of the marginal distribution
#' @param tau.formula Formula for the tau parameter of the marginal distribution
#' @param theta.formula Formula for the theta parameter of the copula distribution
#' @param zeta.formula Formula for the zeta parameter of the copula distribution
#' @param margin.family Marginal distribution specified as a gamlss family object,
#' e.g. GA(), NO(), PO(), NBI(), etc.
#' @param copula.family Copula distribution code, one of "N", "C", "F", "G", "J", or "t".
#' @param copula.link List of link functions for the copula parameters
#' @return Returns a list mm with items mm$x and mm$s for fixed and smooth terms respectively,
#' with each of those lists being lists of each parameter and their respective model matrices
#'
#' @keywords internal
#' @noRd
create_model_matrices <- function(
    mu.formula = ("response ~ 1"),
    sigma.formula = ("1"),
    nu.formula = ("1"),
    tau.formula = ("1"),
    theta.formula = ("1"),
    zeta.formula = ("1"),
    margin.family = NO(),
    copula.family = "N",
    copula.link = NA,
    dataset = NA,
    quiet_gamlss2 = TRUE,
    preserve_factor_levels = FALSE) {
  dataset_mm <- .gl_build_model_matrix_proxy_dataset(dataset)
  included_parameters <- .gl_model_matrix_included_parameters(margin.family, copula.family)
  formulas <- list()

  for (parameter in included_parameters) {
    formulas[[parameter]] <- get(paste(parameter, "formula", sep = "."))
  }

  if (!requireNamespace("mgcv", quietly = TRUE)) {
    stop("Package 'mgcv' is required to construct smooth-term model matrices.")
  }

  formulas[["mu"]] <- as.formula(mu.formula)

  for (parameter in included_parameters[2:length(included_parameters)]) {
    formulas[[parameter]] <- .gl_to_response_formula(formulas[[parameter]], response_name = "response")
  }

  mm_x <- list()
  mm_s <- list()
  smooth_eval_env <- new.env(parent = baseenv())
  smooth_eval_env$s <- mgcv::s

  for (parameter in included_parameters) {
    data_for_par <- .gl_model_matrix_parameter_dataset(dataset_mm, parameter)

    data_for_par <- .gl_sanitize_for_gamlss2(
      data_for_par,
      formulas[[parameter]],
      preserve_factor_levels = preserve_factor_levels
    )

    mm_x[[parameter]] <- .gl_build_fixed_model_matrix(formulas[[parameter]], data_for_par)
    mm_s[[parameter]] <- .gl_build_smooth_model_matrices(formulas[[parameter]], data_for_par, smooth_eval_env)
  }

  mm <- list(mm_x, mm_s)
  names(mm) <- c("x", "s")
  return(mm)
}
