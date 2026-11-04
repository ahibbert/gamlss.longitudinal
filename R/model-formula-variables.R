.gl_formula_variables <- function(...) {
  formulas <- list(...)

  unique(unlist(lapply(formulas, function(fml) all.vars(stats::as.formula(fml))), use.names = FALSE))
}
