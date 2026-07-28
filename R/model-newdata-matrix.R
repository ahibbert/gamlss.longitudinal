.gl_align_model_matrix_columns <- function(mm_use, mm_reference) {
  if (is.null(mm_use) || is.null(mm_reference)) {
    return(mm_use)
  }

  for (par_name in intersect(names(mm_reference$x), names(mm_use$x))) {
    ref_cols <- colnames(mm_reference$x[[par_name]])

    use_cols <- colnames(mm_use$x[[par_name]])

    missing_cols <- setdiff(ref_cols, use_cols)

    if (length(missing_cols) > 0L) {
      for (col_name in missing_cols) {
        mm_use$x[[par_name]][[col_name]] <- 0
      }
    }

    extra_cols <- setdiff(colnames(mm_use$x[[par_name]]), ref_cols)

    if (length(extra_cols) > 0L) {
      mm_use$x[[par_name]] <- mm_use$x[[par_name]][

        ,

        setdiff(colnames(mm_use$x[[par_name]]), extra_cols),
        drop = FALSE
      ]
    }

    mm_use$x[[par_name]] <- mm_use$x[[par_name]][, ref_cols, drop = FALSE]
  }

  mm_use
}
