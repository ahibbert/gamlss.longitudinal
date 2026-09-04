#' @export
best_fit.joint_distribution_selection <- function(x, ...) {
  eligible <- which(
      x$rank == 1L &
      x$selection_eligible &
      x$converged %in% TRUE &
      is.na(x$error)
  )
  if (!length(eligible)) {
    return(list(
      margin_family_name = NA_character_,
      margin_family = NULL,
      copula_family = NA_character_,
      criterion = attr(x, "criterion")
    ))
  }
  row <- .gl_row_as_plain_list(x, eligible[[1L]])
  margin_family_name <- row$margin_family
  copula_family <- row$copula_family
  row$margin_family <- NULL
  row$copula_family <- NULL
  c(
    list(
      margin_family_name = margin_family_name,
      margin_family = .margin_family_object(margin_family_name),
      copula_family = copula_family,
      criterion = attr(x, "criterion")
    ),
    row
  )
}

#' @export
best_fit_family.joint_distribution_selection <- function(x, ...) {
  best <- best_fit(x)
  list(
    margin_dist = best$margin_family,
    copula_dist = best$copula_family
  )
}

#' @export
`$.joint_distribution_selection` <- function(x, name) {
  if (identical(name, "best_fit")) {
    return(best_fit(x))
  }
  .subset2(as.data.frame(x), name, exact = FALSE)
}

#' @export
print.joint_distribution_selection <- function(x, ..., n = 10L) {
  cat("\n<Joint distribution selection>\n")
  if (nrow(x) == 0L) {
    cat("No candidate combinations were retained.\n")
    return(invisible(x))
  }
  selected <- attr(x, "selected")
  cat("Criterion: ", attr(x, "criterion"), "\n", sep = "")
  cat("Selected: ", if (length(selected) != 1L || is.na(selected)) "none" else selected, "\n\n", sep = "")
  cols <- intersect(
    c(
      "rank", "margin_family", "copula_family", "logLik", "AIC", "BIC",
      "EDF", "converged", "selection_eligible"
    ),
    names(x)
  )
  display <- utils::head(as.data.frame(x)[cols], n = n)
  print(display, row.names = FALSE)
  .gl_selection_footer(x, n)
  if ("converged" %in% names(display) && any(!display$converged, na.rm = TRUE)) {
    cat("\nOne or more displayed candidates did not converge; inspect `as.data.frame(x)` for details.\n")
  }
  invisible(x)
}
