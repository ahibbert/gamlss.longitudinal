.gl_cg_hessian_source_paths <- function(package_root = getwd()) {
  file.path(package_root, "R", c(
    "hessian-linkinv-derivatives.R",
    "hessian-fd-step.R",
    "hessian-warnings.R",
    "hessian-margin-cdf.R",
    "hessian-copula.R",
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
    exists_fn("calc_analytical_hessian", mode = "function")
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
    stop("CG requires calc_analytical_hessian(); source the R/hessian-*.R files first.")
  }

  invisible(TRUE)
}
