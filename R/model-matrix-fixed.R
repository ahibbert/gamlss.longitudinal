#' Build a fixed-effect model matrix for one parameter (e.g. mu)
#'
#' @keywords internal
#' @noRd
.gl_build_fixed_model_matrix <- function(formula, data_for_par) {
  formula_terms <- stats::terms(formula)
  has_intercept <- as.integer(attr(formula_terms, "intercept")) == 1L
  term_labels <- attr(formula_terms, "term.labels")

  fixed_terms <- term_labels[!grepl("^\\s*s\\(", term_labels)]

  if (length(fixed_terms) > 0 || has_intercept) {
    fixed_formula <- if (length(fixed_terms) == 0L && has_intercept) {
      stats::as.formula("~ 1")
    } else {
      stats::reformulate(termlabels = fixed_terms, intercept = has_intercept)
    }

    X_fixed <- stats::model.matrix(fixed_formula, data = data_for_par)
    fixed_assign <- attr(X_fixed, "assign")
    fixed_term_labels <- attr(stats::terms(fixed_formula), "term.labels")
    colnames(X_fixed) <- sub("^\\(Intercept\\)$", "intercept", colnames(X_fixed))
    colnames(X_fixed) <- .gl_normalize_time_covariate_colnames(colnames(X_fixed))

    out <- as.data.frame(X_fixed, check.names = FALSE)
    attr(out, "assign") <- fixed_assign
    attr(out, "term.labels") <- fixed_term_labels
    return(out)
  }

  data.frame(row.names = seq_len(nrow(data_for_par)))
}
