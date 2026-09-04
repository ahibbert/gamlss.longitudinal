#' @rdname select_margin

#' @export

screen_margin <- function(...) {
  .Deprecated("select_margin", package = "gamlss.longitudinal")
  select_margin(...)
}

#' @export

print.margin_screen <- function(x, ..., n = 10L) {
  cat("\n<Margin selection>\n")

  if (nrow(x) == 0L) {
    cat("No candidate families were retained.\n")

    return(invisible(x))
  }

  selected <- attr(x, "selected")
  selected_family <- .margin_family_object(selected)
  selected_label <- if (is.null(selected_family)) selected else .gl_margin_display_label(selected_family)
  cat("Criterion: AIC\nSelected: ", selected_label, "\n\n", sep = "")

  unsupported <- "supported_by_longitudinal" %in% names(x) &&
    any(!x$supported_by_longitudinal[seq_len(min(n, nrow(x)))], na.rm = TRUE)

  display <- .gl_selection_display(
    x,
    c("rank", "family", "AIC", "delta_AIC", "screen_model", "converged", "pooled_AIC"),
    n = n,
    digits = 3L
  )

  display$supported_by_longitudinal <- NULL

  print(display, row.names = FALSE)

  .gl_selection_footer(x, n)

  if (isTRUE(unsupported)) {
    cat("\nNote: one or more printed families are not currently supported by gamlss_longitudinal().\n")
  }

  invisible(x)
}
