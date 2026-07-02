#' @export
best_fit.copula_selection <- function(x, ...) {
  if (nrow(x) == 0L) {
    return(list(
      family = NA_character_,
      criterion = attr(x, "criterion")
    ))
  }
  row <- as.list(as.data.frame(x)[1L, , drop = FALSE])
  row$criterion <- attr(x, "criterion")
  row
}

#' @export
best_fit_family.copula_selection <- function(x, ...) {
  best_fit(x)$family
}

#' @export
`$.copula_selection` <- function(x, name) {
  if (identical(name, "best_fit")) {
    return(best_fit(x))
  }
  .subset2(as.data.frame(x), name, exact = FALSE)
}
