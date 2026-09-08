#' Solve a Hessian into covariance output and diagnostics
#'
#' @noRd
.gl_solve_hessian_vcov <- function(H, parameter_names = colnames(H),
                                   inference = inference_control("standard"),
                                   source = "unknown", fallback = NULL,
                                   gradient = NULL, reference_hessian = NULL) {
  .gl_validate_hessian_inference(
    H = H,
    parameter_names = parameter_names,
    control = inference,
    source = source,
    fallback = fallback,
    gradient = gradient,
    reference_hessian = reference_hessian
  )
}

#' Source analytical Hessian helpers for standalone development use
#'
#' @noRd
.gl_source_analytical_hessian_helpers <- function() {
  if (exists("calc_analytical_hessian", mode = "function")) {
    return(invisible(TRUE))
  }

  hessian_files <- c(
    "hessian-linkinv-derivatives.R",
    "hessian-fd-step.R",
    "hessian-warnings.R",
    "hessian-margin-cdf.R",
    "hessian-copula-pair-result.R",
    "hessian-copula-pair-inputs.R",
    "hessian-copula-pair-eval.R",
    "hessian-copula-parameter-curvature.R",
    "hessian-copula-margin-accumulators.R",
    "hessian-copula-margin-pair-terms.R",
    "hessian-copula-margin-terms.R",
    "hessian-copula.R",
    "hessian-discrete-rectangle.R",
    "hessian-margin-derivatives.R",
    "hessian-assembly-helpers.R",
    "hessian-assembly-copula-blocks.R",
    "hessian-assembly-margin-copula-blocks.R",
    "hessian-assembly-margin-blocks.R",
    "hessian-assembly.R",
    "hessian-analytical.R"
  )

  candidate_dirs <- c(
    dirname(attr(body(vcov.gamlss.longitudinal), "srcfile")$filename %||% ""),
    "R",
    file.path(getwd(), "R")
  )

  loaded <- FALSE

  for (dir in candidate_dirs) {
    candidate_paths <- file.path(dir, hessian_files)

    if (all(file.exists(candidate_paths))) {
      for (cp in candidate_paths) source(cp, local = FALSE)
      loaded <- TRUE
      break
    }
  }

  if (!loaded) stop("Cannot locate hessian helper files. Source them manually or ensure the working directory is the package root.")

  invisible(TRUE)
}
