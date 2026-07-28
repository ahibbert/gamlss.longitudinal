.copula_v2_attach_group <- function(pair_data, object, by, data = NULL) {
  if (is.null(by) || (is.character(by) && length(by) == 1 && !nzchar(by))) {
    pair_data$split_group <- factor(pair_data$time_pair)

    return(pair_data)
  }

  if (!is.character(by) || length(by) != 1) {
    stop("'by' must be NULL or a single column name as a character string.")
  }

  if (by %in% c("time", "time_pair")) {
    pair_data$split_group <- factor(pair_data$time_pair)

    return(pair_data)
  }

  if (by %in% c("subject", "lag") && by %in% names(pair_data)) {
    pair_data$split_group <- factor(pair_data[[by]])

    return(pair_data)
  }

  if (is.null(data)) {
    stop("To split plot.copula by '", by, "', provide data= containing that column.")
  }

  df <- as.data.frame(data, stringsAsFactors = FALSE)

  if (!is.null(object$var_map)) {
    for (old_name in names(object$var_map)) {
      new_name <- object$var_map[[old_name]]

      if (old_name %in% names(df) && !new_name %in% names(df)) {
        names(df)[names(df) == old_name] <- new_name
      }
    }
  }

  by_col <- by

  if (!by_col %in% names(df) && !is.null(object$var_map) && by %in% names(object$var_map)) {
    mapped_col <- object$var_map[[by]]

    if (mapped_col %in% names(df)) {
      by_col <- mapped_col
    }
  }

  if (!by_col %in% names(df)) {
    stop("Column '", by, "' not found in provided data after internal name mapping.")
  }

  if (!all(c("subject", "time") %in% names(df))) {
    stop("Provided data must contain subject and time columns (or names mappable via object$var_map) to split by '", by, "'.")
  }

  key_df <- paste(df$subject, as.character(df$time), sep = "::")

  key_pair <- paste(pair_data$subject, as.character(pair_data$time_left), sep = "::")

  matched <- df[[by_col]][match(key_pair, key_df)]

  if (all(is.na(matched))) {
    by_subj <- tapply(df[[by_col]], as.character(df$subject), function(v) {
      vv <- unique(v[!is.na(v)])

      if (length(vv) == 1) vv else NA
    })

    matched <- by_subj[as.character(pair_data$subject)]
  }

  pair_data$split_group <- factor(matched)

  pair_data <- pair_data[!is.na(pair_data$split_group), , drop = FALSE]

  if (nrow(pair_data) == 0) {
    stop("No valid paired rows remained after grouping by '", by, "'.")
  }

  pair_data
}
