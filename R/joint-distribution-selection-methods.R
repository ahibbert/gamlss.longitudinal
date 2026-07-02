#' @export
best_fit.joint_distribution_selection <- function(x, ...) {
  if (nrow(x) == 0L || !is.na(x$error[[1L]]) || !is.finite(x$AIC[[1L]])) {
    return(list(
      margin_family_name = NA_character_,
      margin_family = NULL,
      copula_family = NA_character_,
      criterion = attr(x, "criterion")
    ))
  }
  row <- as.list(as.data.frame(x)[1L, , drop = FALSE])
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
  cat("\nJoint Distribution Screen\n")
  cat("-------------------------\n")
  if (nrow(x) == 0L) {
    cat("No candidate combinations were retained.\n")
    return(invisible(x))
  }
  selected <- attr(x, "selected")
  cat("Selected:", if (length(selected) != 1L || is.na(selected)) "none" else selected, "\n")
  cat("Criterion:", attr(x, "criterion"), "\n\n")
  cols <- intersect(
    c(
      "rank", "margin_family", "copula_family", "logLik", "AIC", "BIC",
      "EDF", "converged", "elapsed_sec", "error"
    ),
    names(x)
  )
  display <- utils::head(as.data.frame(x)[cols], n = n)
  print(display, row.names = FALSE)
  if ("converged" %in% names(display) && any(!display$converged, na.rm = TRUE)) {
    cat("\nWarning: one or more printed joint candidate fits did not report convergence.\n")
  }
  invisible(x)
}
