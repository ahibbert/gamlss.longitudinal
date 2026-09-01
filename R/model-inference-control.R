#' Control validation of Hessian-based inference
#'
#' `inference_control()` defines scale-aware checks applied before standard
#' errors, Wald tests, or confidence intervals are made available. It is
#' deliberately separate from optimizer controls: changing these thresholds
#' cannot change the fitted coefficients.
#'
#' @param profile Validation profile. `"strict"` uses tighter numerical
#'   thresholds and requests an analytical-versus-numerical Hessian check.
#' @param symmetry_tol Maximum relative Hessian asymmetry.
#' @param rank_tol Relative numerical-rank tolerance.
#' @param condition_max Maximum condition number of the diagonally scaled
#'   information matrix.
#' @param gradient_tol Maximum curvature-standardized fitted score.
#' @param agreement_tol Maximum relative scaled discrepancy between analytical
#'   and numerical Hessians when agreement is checked.
#' @param gradient_step Relative step used for the fitted-score check.
#'
#' @return An object of class `gamlss_longitudinal_inference_control`.
#'
#' @details The standard and strict thresholds are versioned provisional
#'   defaults (`0.1.0-provisional`) pending Phase-gate simulation calibration;
#'   they are not claimed to be empirically calibrated. `gradient_tol` checks a
#'   post-fit score standardized by observed curvature. It is unrelated to the
#'   CG optimizer stopping tolerance commonly called `cg_grad_tol`.
#' @export
inference_control <- function(
    profile = c("standard", "strict"),
    symmetry_tol = NULL,
    rank_tol = NULL,
    condition_max = NULL,
    gradient_tol = NULL,
    agreement_tol = NULL,
    gradient_step = NULL) {
  profile <- match.arg(profile)
  defaults <- if (identical(profile, "strict")) {
    list(
      symmetry_tol = 1e-9,
      rank_tol = 1e-9,
      condition_max = 1e8,
      gradient_tol = 5e-2,
      agreement_tol = 5e-2,
      gradient_step = 1e-4,
      check_agreement = TRUE
    )
  } else {
    list(
      symmetry_tol = 1e-7,
      rank_tol = 1e-10,
      condition_max = 1e12,
      gradient_tol = 2.5e-1,
      agreement_tol = 2e-1,
      gradient_step = 1e-4,
      check_agreement = FALSE
    )
  }

  supplied <- list(
    symmetry_tol = symmetry_tol,
    rank_tol = rank_tol,
    condition_max = condition_max,
    gradient_tol = gradient_tol,
    agreement_tol = agreement_tol,
    gradient_step = gradient_step
  )
  override <- supplied[!vapply(supplied, is.null, logical(1))]
  values <- defaults
  for (nm in names(override)) values[[nm]] <- override[[nm]]

  positive <- c("symmetry_tol", "rank_tol", "condition_max", "gradient_tol",
                "agreement_tol", "gradient_step")
  for (nm in positive) {
    value <- values[[nm]]
    if (!is.numeric(value) || length(value) != 1L || !is.finite(value) || value <= 0) {
      stop(sprintf("'%s' must be one positive finite number.", nm), call. = FALSE)
    }
  }

  structure(
    c(list(
      profile = profile,
      defaults_version = "0.1.0-provisional"
    ), values, list(expert_override = override)),
    class = "gamlss_longitudinal_inference_control"
  )
}

#' @noRd
.gl_normalize_inference_control <- function(control = NULL) {
  if (is.null(control)) return(inference_control("standard"))
  if (!inherits(control, "gamlss_longitudinal_inference_control")) {
    stop("'inference' must be created by inference_control().", call. = FALSE)
  }
  control
}
