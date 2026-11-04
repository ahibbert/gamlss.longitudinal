#' @export
best_fit.copula_selection <- function(x, ...) {
  if (nrow(x) == 0L) {
    return(list(
      family = NA_character_,
      criterion = attr(x, "criterion")
    ))
  }
  row <- .gl_row_as_plain_list(x, 1L)
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
  if (name %in% c("tail_cooccurrence", "conditional_tail_exceedance", "tail_thresholds")) {
    return(attr(x, name))
  }
  .subset2(as.data.frame(x), name, exact = FALSE)
}

#' @export
print.copula_selection <- function(x, ..., n = 10L, digits = 3L) {
  cat("\n<Copula selection>\n")
  if (!nrow(x)) {
    cat("No candidate families were retained.\n")
    return(invisible(x))
  }
  selected <- attr(x, "selected")
  cat(
    "Criterion: ", attr(x, "criterion"), " | Complete pairs: ", x$n_pairs[[1L]], "\n",
    "Selected: ", .gl_copula_display_label(selected), "\n\n",
    sep = ""
  )
  display <- .gl_selection_display(
    x,
    c("rank", "family", "par", "par2", "tau", "logLik", "AIC", "BIC", "delta"),
    n = n,
    digits = digits
  )
  display$par2[display$family != "t"] <- NA_real_
  display$family <- .gl_selection_family_labels(display$family, "copula")
  names(display)[names(display) == "par"] <- "theta"
  names(display)[names(display) == "par2"] <- "df"
  print(display, row.names = FALSE, right = FALSE)
  .gl_selection_footer(x, n)
  cat("Tail diagnostics are available in the diagnostic attributes.\n")
  invisible(x)
}
