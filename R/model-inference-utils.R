.fixed_term_coefficient_names <- function(object) {
  if (!inherits(object, "gamlss.longitudinal") || is.null(object$model_matrix$x)) {
    return(NULL)
  }

  out <- list()

  for (parameter in names(object$model_matrix$x)) {
    X <- object$model_matrix$x[[parameter]]

    if (is.null(X)) next

    col_names <- colnames(X)

    if (length(col_names) == 0L) next

    assign <- attr(X, "assign")

    term_labels <- attr(X, "term.labels")

    if (length(assign) != length(col_names) || length(term_labels) == 0L) next

    for (term_idx in seq_along(term_labels)) {
      cols <- col_names[assign == term_idx]

      if (length(cols) == 0L) next

      out[[paste(parameter, term_labels[[term_idx]], sep = ".")]] <- paste(parameter, cols, sep = ".")
    }
  }

  out
}

.resolve_coefficient_terms <- function(terms, coefficient_names, arg = "terms", term_map = NULL) {
  if (is.null(terms)) {
    return(seq_along(coefficient_names))
  }

  if (!is.character(terms)) {
    idx <- terms

    if (length(idx) == 0L || any(is.na(idx)) || any(idx < 1L) || any(idx > length(coefficient_names))) {
      stop(sprintf("'%s' contains unknown coefficient names or indices.", arg), call. = FALSE)
    }

    return(as.integer(idx))
  }

  if (length(terms) == 0L || any(is.na(terms)) || any(!nzchar(terms))) {
    stop(sprintf("'%s' contains unknown coefficient names or indices.", arg), call. = FALSE)
  }

  idx_list <- lapply(terms, function(term) {
    exact <- match(term, coefficient_names)

    if (!is.na(exact)) {
      return(exact)
    }

    if (!is.null(term_map) && term %in% names(term_map)) {
      mapped_idx <- match(term_map[[term]], coefficient_names)

      mapped_idx <- mapped_idx[!is.na(mapped_idx)]

      if (length(mapped_idx) > 0L) {
        return(mapped_idx)
      }
    }

    prefix_idx <- which(startsWith(coefficient_names, term))

    if (length(prefix_idx) == 0L) {
      return(prefix_idx)
    }

    suffix <- substring(coefficient_names[prefix_idx], nchar(term) + 1L)

    main_idx <- prefix_idx[!grepl(":", suffix, fixed = TRUE)]

    if (length(main_idx) > 0L) {
      return(main_idx)
    }

    prefix_idx
  })

  missing <- vapply(idx_list, length, integer(1)) == 0L

  if (any(missing)) {
    stop(sprintf("'%s' contains unknown coefficient names or indices.", arg), call. = FALSE)
  }

  unique(as.integer(unlist(idx_list, use.names = FALSE)))
}
