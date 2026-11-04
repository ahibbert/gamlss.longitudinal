.gl_copula_labels <- c(
  N = "Gaussian", C = "Clayton", F = "Frank", G = "Gumbel",
  J = "Joe", t = "Student t"
)

.gl_margin_display_label <- function(margin_dist) {
  family <- margin_dist$family %||% margin_dist
  code <- if (length(family)) as.character(family[[1L]]) else NA_character_
  label <- if (length(family) >= 2L) as.character(family[[2L]]) else code
  if (!nzchar(code) || is.na(code)) return("unknown")
  if (!nzchar(label) || is.na(label) || identical(label, code)) return(code)
  paste0(label, " (", code, ")")
}

.gl_copula_display_label <- function(code) {
  code <- if (length(code)) as.character(code[[1L]]) else NA_character_
  if (!nzchar(code) || is.na(code)) return("unknown")
  label <- unname(.gl_copula_labels[[code]])
  if (is.null(label) || is.na(label)) return(code)
  paste0(label, " (", code, ")")
}

.gl_format_console_number <- function(x, digits = max(3L, getOption("digits") - 3L)) {
  if (!length(x) || is.na(x) || !is.finite(x)) return("n/a")
  formatC(x, digits = digits, format = "fg", flag = "#")
}

.gl_selection_family_labels <- function(code, component = c("copula", "margin")) {
  component <- match.arg(component)
  if (identical(component, "copula")) {
    return(vapply(code, .gl_copula_display_label, character(1)))
  }
  as.character(code)
}

.gl_selection_display <- function(x, columns, n, digits) {
  columns <- intersect(columns, names(x))
  out <- utils::head(as.data.frame(x)[columns], n = n)
  numeric_cols <- names(out)[vapply(out, is.numeric, logical(1))]
  for (nm in numeric_cols) out[[nm]] <- round(out[[nm]], digits)
  out
}

.gl_selection_footer <- function(x, n) {
  remaining <- nrow(x) - min(nrow(x), n)
  if (remaining > 0L) {
    cat("\n", remaining, " additional candidate", if (remaining == 1L) "" else "s",
        " not shown. Use `as.data.frame(x)` to inspect all rows.\n", sep = "")
  }
}

.gl_row_as_plain_list <- function(x, row) {
  setNames(lapply(names(x), function(nm) x[[nm]][[row]]), names(x))
}

.gl_progress <- function(verbose, ..., prefix = NULL) {
  if (verbose <= 0) return(invisible(FALSE))
  if (!is.null(prefix)) cat("[", prefix, "] ", sep = "")
  cat(..., "\n", sep = "")
  invisible(TRUE)
}
