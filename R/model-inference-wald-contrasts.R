#' Build Wald-test contrast matrix
#'
#' @param L Optional user-supplied contrast matrix.
#' @param terms Optional coefficient selectors used when `L` is `NULL`.
#' @param estimates Named fixed-effect coefficient vector.
#' @param object Fitted model object, used for formula-term coefficient maps.
#' @param joint Logical joint-test request.
#' @return List with full contrast matrix `L` and resolved `joint` flag.
#' @noRd
.gl_wald_contrast_matrix <- function(L, terms, estimates, object, joint) {
  if (!is.null(L)) {
    L <- as.matrix(L)

    if (nrow(L) == 0L || ncol(L) == 0L) {
      stop("'L' must contain at least one contrast row and one coefficient column.", call. = FALSE)
    }

    if (!is.null(colnames(L))) {
      idx <- match(colnames(L), names(estimates))

      if (any(is.na(idx))) {
        stop("Column names in 'L' must match coefficient names.", call. = FALSE)
      }

      L_full <- matrix(0, nrow = nrow(L), ncol = length(estimates))

      colnames(L_full) <- names(estimates)

      rownames(L_full) <- rownames(L)

      L_full[, idx] <- L

      L <- L_full
    } else if (ncol(L) != length(estimates)) {
      stop("Unnamed 'L' must have one column per fixed coefficient.", call. = FALSE)
    } else {
      colnames(L) <- names(estimates)
    }

    if (is.null(rownames(L))) {
      rownames(L) <- paste0("H", seq_len(nrow(L)))
    }

    joint <- TRUE
  } else {
    idx <- .resolve_coefficient_terms(
      terms,
      names(estimates),
      arg = "terms",
      term_map = .fixed_term_coefficient_names(object)
    )

    L <- diag(length(estimates))[idx, , drop = FALSE]

    colnames(L) <- names(estimates)

    rownames(L) <- names(estimates)[idx]
  }

  list(L = L, joint = joint)
}
