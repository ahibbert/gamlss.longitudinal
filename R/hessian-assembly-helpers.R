.hessian_assembly_param_of <- function(coef_name) {
  match <- regexpr("^[^.]+", coef_name)
  regmatches(coef_name, match)
}

.hessian_assembly_x_col <- function(coef_name, param_block, mm_x) {
  block <- param_block[coef_name]
  if (!block %in% names(mm_x)) {
    return(NULL)
  }

  X <- mm_x[[block]]
  col_name <- sub(paste0("^", block, "\\."), "", coef_name)
  if (identical(col_name, coef_name)) col_name <- "(Intercept)"

  if (!col_name %in% colnames(X)) {
    # Smooth columns can be appended to mm$x as numbered basis columns while
    # the coefficient names retain the full smooth label.
    basis_col <- sub("^.*\\.([0-9]+)$", "\\1", col_name)
    if (!identical(basis_col, col_name) && basis_col %in% colnames(X)) {
      col_name <- basis_col
    }
  }

  if (!col_name %in% colnames(X)) {
    return(NULL)
  }
  X[, col_name, drop = TRUE]
}

.hessian_assembly_margin_d2_obs <- function(pn1, pn2, margin_d2l, copula_hess, eta_dr) {
  margin_pars <- names(margin_d2l)
  n_obs <- length(eta_dr[[if (pn1 %in% names(eta_dr)) pn1 else pn2]])
  value <- numeric(n_obs)

  if (pn1 %in% margin_pars && pn2 %in% margin_pars) {
    if (!is.null(margin_d2l[[pn1]][[pn2]])) {
      value <- value + margin_d2l[[pn1]][[pn2]]
    }
    if (!is.null(copula_hess$cop_d2l_margin[[pn1]][[pn2]])) {
      value <- value + copula_hess$cop_d2l_margin[[pn1]][[pn2]]
    }
  }

  value
}

.hessian_assembly_margin_d1_obs <- function(pn, margin_d1l, copula_hess, eta_dr) {
  n_obs <- length(eta_dr[[pn]])
  value <- numeric(n_obs)

  if (pn %in% names(margin_d1l)) {
    value <- value + margin_d1l[[pn]]
  }
  if (!is.null(copula_hess$cop_d1l_margin) && pn %in% names(copula_hess$cop_d1l_margin)) {
    value <- value + copula_hess$cop_d1l_margin[[pn]]
  }

  value
}
