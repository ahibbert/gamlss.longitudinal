.gl_translate_newdata_names <- function(nd, var_map) {
  # Translate user variable names to internal names used by model formulas.

  for (old_name in names(var_map)) {
    new_name <- var_map[[old_name]]

    if (old_name %in% names(nd) && !new_name %in% names(nd)) {
      names(nd)[names(nd) == old_name] <- new_name
    }
  }

  nd
}

.gl_add_newdata_default_columns <- function(nd, object) {
  # Fitting keeps the original user time scale as `time_covariate` for

  # formulas, while `time` is used internally for ordering/pairing. Recreate

  # that column for prediction data supplied with either original or internal

  # names.

  if (!"time_covariate" %in% names(nd) && "time" %in% names(nd)) {
    nd$time_covariate <- nd$time
  }

  if (!"time" %in% names(nd) && "time" %in% names(object$model_matrix$x$mu)) {
    nd$time <- NA
  }

  if (!"subject" %in% names(nd)) {
    nd$subject <- seq_len(nrow(nd))
  }

  if (!"response" %in% names(nd)) {
    nd$response <- NA_real_
  }

  nd
}

.gl_align_newdata_factor_levels <- function(nd, object) {
  if (!is.null(object$dataset)) {
    factor_cols <- names(object$dataset)[vapply(object$dataset, is.factor, logical(1))]

    for (nm in intersect(factor_cols, names(nd))) {
      train_col <- object$dataset[[nm]]

      train_levels <- levels(train_col)

      nd_values <- as.character(nd[[nm]])

      unknown <- setdiff(unique(nd_values[!is.na(nd_values)]), train_levels)

      if (length(unknown) > 0L) {
        stop(

          "newdata column '", nm, "' contains level(s) not seen during fitting: ",
          paste(unknown, collapse = ", "),
          call. = FALSE
        )
      }

      nd[[nm]] <- factor(nd_values, levels = train_levels, ordered = is.ordered(train_col))

      if (is.ordered(train_col) && length(train_levels) > 1L) {
        contr <- contr.treatment(length(train_levels))

        colnames(contr) <- train_levels[-1]

        contrasts(nd[[nm]]) <- contr
      }
    }
  }

  nd
}

.gl_validate_newdata_response <- function(nd, require_response) {
  if (require_response && all(is.na(nd$response))) {
    stop("newdata must include a response column (or mapped response variable) for this operation.")
  }

  invisible(TRUE)
}
