.gl_validate_tabular_shape <- function(dataset, context = "dataset") {
  if (!is.data.frame(dataset)) {
    return(invisible(NULL))
  }

  if (nrow(dataset) == 0L) {
    stop("ERROR: ", context, " must contain at least one row.", call. = FALSE)
  }

  if (ncol(dataset) == 0L) {
    stop("ERROR: ", context, " must contain at least one column.", call. = FALSE)
  }


  list_cols <- names(dataset)[vapply(dataset, is.list, logical(1))]

  if (length(list_cols) > 0L) {
    stop(

      "ERROR: ", context, " contains unsupported list-column(s): ",
      paste(list_cols, collapse = ", "),
      ". Flatten or remove list-columns before fitting.",
      call. = FALSE
    )
  }

  matrix_cols <- names(dataset)[vapply(dataset, function(x) is.matrix(x) || is.data.frame(x), logical(1))]
  if (length(matrix_cols) > 0L) {
    stop(

      "ERROR: ", context, " contains unsupported matrix/data-frame column(s): ",
      paste(matrix_cols, collapse = ", "),
      ". Expand these columns into ordinary vector columns before fitting.",
      call. = FALSE
    )
  }

  invisible(NULL)
}
