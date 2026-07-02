.joint_selection_fit_formulas <- function(
    response_var,
    time_var,
    time_intercepts,
    copula_time_intercepts) {
  mu_formula <- if (isTRUE(time_intercepts)) {
    stats::as.formula(paste(.select_copula_formula_name(response_var), "~ factor(", .select_copula_formula_name(time_var), ")"))
  } else {
    stats::reformulate("1", response = response_var)
  }
  par_formula <- if (isTRUE(time_intercepts)) {
    stats::as.formula(paste("~ factor(", .select_copula_formula_name(time_var), ")"))
  } else {
    ~1
  }
  theta_formula <- if (isTRUE(copula_time_intercepts)) {
    stats::as.formula(paste("~ factor(", .select_copula_formula_name(time_var), ")"))
  } else {
    ~1
  }

  list(
    mu_formula = mu_formula,
    par_formula = par_formula,
    theta_formula = theta_formula
  )
}
