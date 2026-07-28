.plot_smooth_terms_extract_var <- function(s_name) {
  s_txt <- trimws(s_name)
  s_call <- tryCatch(parse(text = s_txt)[[1]], error = function(e) NULL)

  if (!is.null(s_call) && length(s_call) >= 2) {
    out <- paste(deparse(s_call[[2]]), collapse = " ")
  } else {
    out <- sub("^s\\((.*)\\)$", "\\1", s_txt)
  }

  trimws(gsub("`", "", out, fixed = TRUE))
}

.plot_smooth_terms_eval_x <- function(x_expr, data_frame) {
  if (is.null(data_frame) || !is.data.frame(data_frame)) {
    return(NULL)
  }

  tryCatch(eval(parse(text = x_expr)[[1]], envir = data_frame), error = function(e) NULL)
}

.plot_smooth_terms_x_info <- function(par_name,
                                      s_name,
                                      B,
                                      object,
                                      data = NULL,
                                      fallback_to_index = TRUE) {
  x_var <- .plot_smooth_terms_extract_var(s_name)
  x <- NULL

  x_basis <- attr(B, "smooth_x")
  x_basis_var <- attr(B, "smooth_var")
  if (!is.null(x_basis) && length(x_basis) == nrow(B)) {
    x <- x_basis
    if (!is.null(x_basis_var) && nzchar(x_basis_var)) {
      x_var <- x_basis_var
    }
  }

  if (is.null(x) && !is.null(data) && is.data.frame(data)) {
    data_names <- names(data)
    idx_exact <- which(data_names == x_var)
    idx_ci <- which(tolower(data_names) == tolower(x_var))
    idx_mn <- which(make.names(data_names) == make.names(x_var))
    idx <- c(idx_exact, idx_ci, idx_mn)
    idx <- idx[!duplicated(idx)]

    if (length(idx) > 0) {
      matched_name <- data_names[idx[1]]
      x_candidate <- data[[matched_name]]
      if (length(x_candidate) == nrow(B)) {
        x <- x_candidate
      } else if (!is.null(rownames(B)) && !is.null(rownames(data))) {
        row_idx <- match(rownames(B), rownames(data))
        if (all(!is.na(row_idx))) {
          x <- x_candidate[row_idx]
        }
      }
      x_var <- matched_name
    }

    if (is.null(x)) {
      x_candidate <- .plot_smooth_terms_eval_x(x_var, data)
      if (!is.null(x_candidate) && length(x_candidate) == nrow(B)) {
        x <- x_candidate
      }
    }
  }

  if (is.null(x) &&
    !is.null(object$model_matrix$x[[par_name]]) &&
    x_var %in% colnames(object$model_matrix$x[[par_name]])) {
    x <- object$model_matrix$x[[par_name]][, x_var]
  }

  if (is.null(x) && fallback_to_index) {
    x <- seq_len(nrow(B))
    x_var <- "index"
    warning(
      "Falling back to index for smooth term '", s_name,
      "' because covariate was not found in supplied data/model matrix."
    )
  }

  if (is.null(x)) {
    stop(
      "Could not infer x-axis for smooth term '", s_name,
      "'. Provide 'data' with the smooth covariate columns."
    )
  }

  if (length(x) != nrow(B)) {
    stop(
      "Length mismatch for smooth term '", s_name, "': length(x)=",
      length(x), " but nrow(B)=", nrow(B), "."
    )
  }

  list(x = as.numeric(x), x_var = x_var)
}
