#' Turns ordered factors into factors with a baseline and effects
#'
#' @noRd
.gl_model_matrix_normalize_ordered_factor <- function(col) {
  if (!is.factor(col)) {
    return(col)
  }

  if (!is.ordered(col)) {
    return(col)
  }

  levs <- levels(col)
  col_nom <- factor(as.character(col), levels = levs, ordered = FALSE)

  if (length(levs) > 1) {
    contr <- contr.treatment(length(levs))
    colnames(contr) <- levs[-1]
    contrasts(col_nom) <- contr
  }

  col_nom
}

#' Return NA if all values are NA, otherwise return the most common non-NA value
#'
#' @noRd
.gl_model_matrix_mode_value <- function(x) {
  x_non_na <- x[!is.na(x)]

  if (length(x_non_na) == 0) {
    return(NA)
  }

  tab <- table(x_non_na)

  names(tab)[which.max(tab)]
}

#' Build the NA-free proxy data used only for model-matrix construction
#' This does not alter likelihood calculations (which still use original dataset).
#' This is just used as a helper during data transformations to account for NA values that have been inserted into the dataset for time-varying copula parameters, since those NA values would cause model matrix construction to fail.
#' 
#' For numeric columns, NA values are replaced with the mean of the non-NA values (or 0 if all values are NA). For factor columns, NA values are replaced with the most common level (or "missing" if all values are NA). For other types of columns, NA values are replaced with the most common value (or "missing" if all values are NA).
#' The "time" and "subject" columns are not altered, since they are not used in model matrix construction and we want to preserve their NA values for later steps.
#' If the response variable has NA values, those are also replaced with the mean of the non-NA values (or 0 if all values are NA), since the response variable is needed for model matrix construction and cannot have NA values.
#' This function is used in the create_model_matrices() function to create a proxy dataset for model matrix construction that does not have NA values, which allows the model matrix construction to proceed without errors due to NA values. The original dataset with NA values is still used for likelihood calculations and other steps, so this function does not alter the actual data used for modeling, it just creates a temporary dataset for model matrix construction.
#' This is essentially a workaround to handle the fact that we need to insert NA values into the dataset for time-varying copula parameters, but those NA values would cause model matrix construction to fail. By creating a proxy dataset that fills in those NA values with reasonable defaults, we can allow model matrix construction to proceed without errors, while still using the original dataset with NA values for likelihood calculations and other steps.
#' This function is not intended to be used outside of the model matrix construction process, and the resulting dataset should not be used for any other purpose, since it may have altered values that are not representative of the original data.
#' 
#' @keywords internal
#' @noRd
.gl_build_model_matrix_proxy_dataset <- function(dataset) {
  dataset_mm <- dataset

  # Build model matrices from an NA-free proxy dataset.
  for (nm in names(dataset_mm)) {
    if (!any(is.na(dataset_mm[[nm]]))) next
    if (nm %in% c("time", "subject")) next
    col <- dataset_mm[[nm]]

    if (is.numeric(col) || is.integer(col)) {
      obs <- col[!is.na(col)]
      fill_val <- if (length(obs) > 0) mean(obs) else 0
      col[is.na(col)] <- fill_val
      dataset_mm[[nm]] <- col
    } else if (is.factor(col)) {
      col <- .gl_model_matrix_normalize_ordered_factor(col)
      fill_val <- .gl_model_matrix_mode_value(col)
      if (is.na(fill_val)) {
        fill_val <- if (length(levels(col)) > 0) levels(col)[1] else "missing"
      }
      col_chr <- as.character(col)
      col_chr[is.na(col_chr)] <- fill_val
      dataset_mm[[nm]] <- factor(col_chr, levels = levels(col), ordered = FALSE)
    } else {
      fill_val <- .gl_model_matrix_mode_value(col)
      if (is.na(fill_val)) fill_val <- "missing"
      col[is.na(col)] <- fill_val
      dataset_mm[[nm]] <- col
    }
  }

  if ("response" %in% names(dataset_mm) && any(is.na(dataset_mm$response))) {
    obs_resp <- dataset_mm$response[!is.na(dataset_mm$response)]
    fill_val <- if (length(obs_resp) > 0) mean(obs_resp) else 0
    dataset_mm$response[is.na(dataset_mm$response)] <- fill_val
  }

  dataset_mm
}

#' Prepare the per-parameter data passed to gamlss2/mgcv matrix routines
#'
#' @noRd
.gl_sanitize_for_gamlss2 <- function(data_in, fml, preserve_factor_levels = FALSE) {
  vars_needed <- unique(all.vars(stats::as.formula(fml)))
  vars_needed <- vars_needed[vars_needed %in% names(data_in)]
  data_out <- data_in[, vars_needed, drop = FALSE]

  for (nm in names(data_out)) {
    col <- data_out[[nm]]
    if (is.factor(col)) {
      col <- .gl_model_matrix_normalize_ordered_factor(col)
      if (!isTRUE(preserve_factor_levels)) {
        col <- droplevels(col)
      }

      data_out[[nm]] <- col
    } else if (is.numeric(col) || is.integer(col)) {
      col[!is.finite(col)] <- NA

      if (any(is.na(col))) {
        obs <- col[!is.na(col)]
        fill_val <- if (length(obs) > 0) mean(obs) else 0
        col[is.na(col)] <- fill_val
      }

      data_out[[nm]] <- col
    } else {
      if (any(is.na(col))) {
        x_non_na <- col[!is.na(col)]
        fill_val <- if (length(x_non_na) > 0) {
          tab <- table(x_non_na)
          names(tab)[which.max(tab)]
        } else {
          "missing"
        }
        col[is.na(col)] <- fill_val
      }
      data_out[[nm]] <- col
    }
  }
  data_out
}
