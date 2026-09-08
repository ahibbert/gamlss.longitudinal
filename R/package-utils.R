#' @importFrom rlang .data

#' @importFrom graphics abline contour hist par plot.new

#' @importFrom methods formalArgs

#' @importFrom stats aggregate as.formula ave coef complete.cases contr.treatment contrasts "contrasts<-" cor effects fitted formula model.frame nobs optim qnorm rbeta residuals rgamma rnorm runif setNames terms var

#' @importFrom utils capture.output head

#' @importFrom gamlss gamlss coefAll

#' @importFrom gamlss.dist NO qZISICHEL

#' @importFrom generics augment glance tidy
NULL

########### NEW SIMPLIFIED FUNCTIONS

#' Null-coalescing helper

#'

#' Returns the left-hand value unless it is `NULL`, in which case the

#' right-hand value is used. This keeps optional fitted-object fields and

#' fallback controls concise without changing falsey values such as `FALSE` or

#' `0`.

#'

#' @noRd

`%||%` <- function(a, b) if (!is.null(a)) a else b

utils::globalVariables(c(
  "u1", "u2", "quartile", "tau_emp", "tau_fit", "density", "x_id", "time_pair", "split_group",
  "time", "z", "z_prev", "z_curr", "empirical", "fitted", "threshold", "tail", "probability",
  "emp_copula", "fit_copula", "lag", "cor_z", "n_pairs", "source", "cut_group",
  "x", "y", "X1", "X2", "idx", "response", "time_plot", "person_plot", "group",
  "density_emp", "density_fit", "density_diff", "theta_fit", "zeta_fit",
  "ci_lower", "ci_upper",
  "copula_dist", "dataset", "dlcopdpar", "eta", "margin_deriv_names",
  "nd_cross_m", "nd_impact", "rand_mvt", "row_id1"
))

#' Solve a linear system with numerical fallbacks

#'

#' Attempts Cholesky, ordinary solve, QR solve, and finally a pseudo-inverse.

#' Used by variance-covariance helpers where nearly singular Hessians should

#' degrade to a guarded estimate rather than an immediate low-level error.

#'

#' @noRd

.solve_linear_system <- function(A, b = NULL) {
  A <- as.matrix(A)

  if (is.null(b)) {
    b_mat <- diag(nrow(A))
  } else {
    b_mat <- as.matrix(b)
  }

  # Cholesky is fastest and most stable for positive-definite systems.

  chol_A <- tryCatch(chol(A), error = function(e) NULL)

  if (!is.null(chol_A)) {
    sol <- backsolve(chol_A, forwardsolve(t(chol_A), b_mat))

    if (is.null(b)) {
      return(sol)
    }

    return(sol)
  }

  # General full-rank solve.

  sol <- tryCatch(solve(A, b_mat), error = function(e) NULL)

  if (!is.null(sol)) {
    return(sol)
  }

  # Rank-revealing fallback.

  sol <- tryCatch(qr.solve(A, b_mat), error = function(e) NULL)

  if (!is.null(sol)) {
    return(sol)
  }

  # Last-resort pseudo-inverse for near-singular systems.

  MASS::ginv(A) %*% b_mat
}

#' Detect discrete or count-like GAMLSS margins

#'

#' Identifies margins that need rectangle probabilities or count-tail handling

#' rather than continuous-density approximations.

#'

#' @noRd

.is_discrete_margin <- function(margin_dist) {
  identical(.gl_capability_likelihood_route(margin_dist), "exact_discrete_rectangle")
}

#' Fixed distribution arguments supplied by the longitudinal likelihood
#'
#' @noRd
.gl_margin_fixed_family_args <- function(margin_dist, n) {
  family <- if (is.character(margin_dist)) {
    margin_dist[[1L]]
  } else {
    margin_dist$family[[1L]]
  }
  family <- as.character(family)
  out <- list()

  if (identical(family, "BI")) {
    out$bd <- rep(1, n)
  }

  out
}

#' Check whether a GAMLSS margin parameter has complete link metadata

#'

#' A parameter is considered fit-ready only when the family object exposes the

#' link, inverse link, and derivative functions used by model-matrix and

#' Hessian calculations.

#'

#' @noRd

.margin_parameter_has_link <- function(margin_dist, par_name) {
  all(vapply(
    paste0(par_name, c(".linkfun", ".linkinv", ".dr")),
    function(nm) is.function(margin_dist[[nm]]),
    logical(1)
  ))
}

#' Drop unsupported unlinked GAMLSS family parameters

#'

#' Some `gamlss.dist` family objects expose parameters without the full link

#' interface needed for longitudinal fitting. This helper fixes those at their

#' family defaults and records the fixed values as attributes.

#'

#' @noRd

.normalise_margin_dist_links <- function(margin_dist) {
  parameter_names <- names(margin_dist$parameters)

  if (length(parameter_names) == 0L) {
    return(margin_dist)
  }

  linked <- vapply(parameter_names, .margin_parameter_has_link, logical(1), margin_dist = margin_dist)

  if (all(linked)) {
    return(margin_dist)
  }

  dropped <- parameter_names[!linked]

  qfun <- tryCatch(

    get(paste0("q", margin_dist$family[1]), envir = asNamespace("gamlss.dist"), inherits = FALSE),
    error = function(e) NULL
  )

  fixed_values <- stats::setNames(as.list(rep(NA_real_, length(dropped))), dropped)

  if (!is.null(qfun)) {
    q_formals <- formals(qfun)

    for (par_name in dropped) {
      value <- tryCatch(eval(q_formals[[par_name]], envir = baseenv()), error = function(e) NA_real_)

      fixed_values[[par_name]] <- as.numeric(value)[1L]
    }
  }

  margin_dist$parameters <- margin_dist$parameters[linked]

  margin_dist$nopar <- length(margin_dist$parameters)

  attr(margin_dist, "fixed_unlinked_parameters") <- dropped

  attr(margin_dist, "fixed_unlinked_values") <- fixed_values

  margin_dist
}
