#' Collect all model matric components for all parameters, including fixed and smooth terms, and copula links
#' 
#' This is essentially a wrapper that orchestrates the construction of all model matrices and copula links, 
#' It is called internally within the model fitting process and is not intended for direct use by users of the package.
#' 
#' @keywords internal
#' @noRd
.gl_build_model_matrix_bundle <- function(
    formulas_int,
    margin_dist,
    copula_dist,
    dataset) {
  copula_link <- get_copula_dist(copula_dist)$copula_link

  mm <- suppressWarnings(create_model_matrices(
    formulas_int$mu,
    formulas_int$sigma,
    formulas_int$nu,
    formulas_int$tau,
    formulas_int$theta,
    formulas_int$zeta,
    margin.family = margin_dist,
    copula.family = copula_dist,
    copula.link = copula_link,
    dataset = dataset
  ))

  .gl_warn_rank_deficient_model_matrices(mm$x)

  list(
    mm = mm,
    copula_link = copula_link
  )
}


#' Warn about rank-deficient fixed-effect model matrices
#'
#' @noRd
.gl_warn_rank_deficient_model_matrices <- function(mm_x) {
  for (parameter in names(mm_x)) {
    X <- mm_x[[parameter]]
    if (!is.matrix(X) && !is.data.frame(X)) next
    X <- as.matrix(X)
    if (nrow(X) == 0L || ncol(X) < 2L) next
    finite_rows <- stats::complete.cases(X)
    if (sum(finite_rows) < 2L) next
    qr_x <- qr(X[finite_rows, , drop = FALSE])
    if (qr_x$rank < ncol(X)) {
      warning(
        "Fixed-effect model matrix for parameter '", parameter,
        "' is rank deficient; estimates may be non-identifiable.",
        call. = FALSE
      )
    }
  }
  invisible(NULL)
}
