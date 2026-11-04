#' @importFrom rlang .data

#' @importFrom graphics abline contour hist par plot.new

#' @importFrom methods formalArgs

#' @importFrom stats aggregate as.formula ave BIC coef complete.cases contr.treatment contrasts "contrasts<-" cor effects fitted formula model.frame nobs optim qnorm rbeta residuals rgamma rnorm runif setNames terms var

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
  "nd_cross_m", "nd_impact", "row_id1"
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
  out <- attr(margin_dist, "fixed_family_args", exact = TRUE)
  if (is.null(out)) return(list())
  lapply(out, function(value) {
    if (length(value) == 1L) return(rep(value, n))
    if (length(value) != n) {
      stop("Stored marginal family arguments do not match the evaluation rows.", call. = FALSE)
    }
    value
  })
}

.gl_binomial_families <- function() {
  c("BI", "BB", "DBI", "ZABB", "ZABI", "ZIBB", "ZIBI")
}

.gl_attach_margin_args_spec <- function(margin_dist, margin_args = list()) {
  if (is.null(margin_args)) margin_args <- list()
  if (!is.list(margin_args) || (length(margin_args) && (
      is.null(names(margin_args)) || any(!nzchar(names(margin_args))) ||
        anyDuplicated(names(margin_args))))) {
    stop("'margin_args' must be a named list.", call. = FALSE)
  }
  unknown <- setdiff(names(margin_args), "bd")
  if (length(unknown)) {
    stop("Unsupported 'margin_args': ", paste(unknown, collapse = ", "), ".", call. = FALSE)
  }
  attr(margin_dist, "margin_args_spec") <- margin_args
  attr(margin_dist, "fixed_family_args") <- NULL
  margin_dist
}

.gl_resolve_margin_args <- function(margin_dist, dataset, response = NULL) {
  family <- .gl_capability_margin_code(margin_dist)
  spec <- attr(margin_dist, "margin_args_spec", exact = TRUE) %||% list()
  is_binomial <- family %in% .gl_binomial_families()
  if (!is_binomial && length(spec)) {
    stop("'margin_args' is currently available only for binomial GAMLSS families.", call. = FALSE)
  }
  if (!is_binomial) return(list())
  if (!"bd" %in% names(spec)) {
    .gl_capability_stop(
      "gamlss_longitudinal_missing_margin_argument_error",
      paste0(
        "Margin family '", family, "' requires explicit binomial trials. ",
        "Supply margin_args = list(bd = 10) for a constant denominator or ",
        "margin_args = list(bd = 'trials_column') for row-specific denominators."
      ),
      margin_family = family, margin_argument = "bd"
    )
  }
  bd_spec <- spec$bd
  if (is.character(bd_spec) && length(bd_spec) == 1L && !is.na(bd_spec)) {
    if (!is.data.frame(dataset) || !bd_spec %in% names(dataset)) {
      stop("Binomial trials column '", bd_spec, "' was not found in the data.", call. = FALSE)
    }
    bd <- dataset[[bd_spec]]
  } else if (is.numeric(bd_spec) && length(bd_spec) == 1L && is.finite(bd_spec)) {
    bd <- rep(as.numeric(bd_spec), nrow(dataset))
  } else {
    stop("'margin_args$bd' must be one positive integer or one data-column name.", call. = FALSE)
  }
  if (!is.numeric(bd) && !is.integer(bd)) {
    stop("Binomial trials must be numeric integers.", call. = FALSE)
  }
  observed <- if (is.null(response)) rep(TRUE, length(bd)) else !is.na(response)
  invalid_bd <- observed & (!is.finite(bd) | bd <= 0 | abs(bd - round(bd)) > sqrt(.Machine$double.eps))
  if (any(invalid_bd)) {
    stop("Binomial trials must be finite positive integers for every observed response.", call. = FALSE)
  }
  if (!is.null(response)) {
    invalid_y <- observed & (response < 0 | abs(response - round(response)) > sqrt(.Machine$double.eps) |
      response > bd)
    if (any(invalid_y)) {
      .gl_capability_stop(
        "gamlss_longitudinal_response_domain_error",
        "Binomial responses must be non-negative integer successes no greater than their row-specific trials.",
        margin_family = family, response_domain = "bounded_count"
      )
    }
  }
  list(bd = as.numeric(bd))
}

.gl_validate_margin_args_contract <- function(margin_dist, dataset, response = NULL) {
  invisible(.gl_resolve_margin_args(margin_dist, dataset, response))
}

.gl_bind_margin_args <- function(margin_dist, dataset) {
  response <- if (is.data.frame(dataset) && "response" %in% names(dataset)) dataset$response else NULL
  attr(margin_dist, "fixed_family_args") <- .gl_resolve_margin_args(
    margin_dist, dataset, response
  )
  margin_dist
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
