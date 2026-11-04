.gl_cg_hessian_source_paths <- function(package_root = getwd()) {
  file.path(package_root, "R", c(
    "copula-backend-mixed-arguments.R",
    "copula-backend-parameter-stencil.R",
    "copula-backend-parameter-curvature.R",
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
  ))
}

#' Ensure analytical Hessian helpers are available for CG fitting
#'
#' @noRd
.gl_ensure_cg_hessian_available <- function(
    package_root = getwd(),
    exists_fn = exists,
    source_fn = source,
    file_exists_fn = file.exists) {
  has_hessian <- function() {
    all(vapply(c("calc_analytical_hessian", ".calc_copula_hessian_contributions",
      ".copula_hessian_pair_derivatives", ".copula_density_argument_curvature",
      ".copula_parameter_argument_mixed", ".copula_hessian_parameter_log_curvature",
      ".copula_t_df_derivative_domain"),
      function(name) exists_fn(name, mode = "function"), logical(1L)))
  }

  if (!has_hessian()) {
    hess_paths <- .gl_cg_hessian_source_paths(package_root)

    if (all(file_exists_fn(hess_paths))) {
      for (hess_path in hess_paths) {
        source_fn(hess_path, local = FALSE)
      }
    }
  }

  if (!has_hessian()) {
    stop("CG requires the complete analytical Hessian helpers; source R/hessian-*.R and R/copula-backend-*.R first.")
  }

  invisible(TRUE)
}
