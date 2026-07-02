.select_copula_margin_selection_time_var <- function(margin_dist) {
  if (!inherits(margin_dist, "margin_selection") || !isTRUE(attr(margin_dist, "time_intercepts"))) {
    return(NULL)
  }
  time_var <- attr(margin_dist, "time_var")
  if (is.null(time_var) || !is.character(time_var) || length(time_var) != 1L || !nzchar(time_var)) {
    return(NULL)
  }
  time_var
}

.select_copula_time_response_formula <- function(formula, response_var, time_var) {
  if (!is.null(formula)) {
    return(formula)
  }
  stats::as.formula(paste(.select_copula_formula_name(response_var), "~ factor(", .select_copula_formula_name(time_var), ")"))
}

.select_copula_time_rhs_formula <- function(formula, time_var) {
  if (!is.null(formula)) {
    return(formula)
  }
  stats::as.formula(paste("~ factor(", .select_copula_formula_name(time_var), ")"))
}

.select_copula_formula_name <- function(name) {
  if (make.names(name) == name) {
    return(name)
  }
  paste0("`", gsub("`", "\\\\`", name), "`")
}

.select_copula_response_formula <- function(formula, response_var) {
  if (is.null(formula)) {
    return(stats::as.formula(paste(response_var, "~ 1")))
  }
  formula <- stats::as.formula(formula)
  if (length(formula) == 2L) {
    return(stats::as.formula(paste(response_var, deparse(formula), sep = " ")))
  }
  formula
}

.select_copula_rhs_formula <- function(formula) {
  if (is.null(formula)) {
    return(~1)
  }
  formula <- stats::as.formula(formula)
  if (length(formula) == 3L) {
    return(stats::as.formula(paste("~", deparse(formula[[3L]]))))
  }
  formula
}
