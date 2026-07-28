.gl_has_supported_base_column_class <- function(x) {
  if (is.factor(x)) {
    return(TRUE)
  }

  type <- typeof(x)

  if (!type %in% c("double", "integer", "logical", "character")) {
    return(FALSE)
  }

  base_class <- switch(type,
    double = "numeric",
    integer = "integer",
    logical = "logical",
    character = "character"
  )

  identical(class(x), base_class)
}

.gl_validate_fitting_data_policy <- function(dataset, formulas, response_name = "response") {
  formula_vars <- unique(unlist(lapply(formulas, function(fml) all.vars(stats::as.formula(fml))), use.names = FALSE))

  missing_vars <- setdiff(formula_vars, names(dataset))

  if (length(missing_vars) > 0L) {
    stop(

      "ERROR: formula variable(s) not found in dataset after internal name mapping: ",
      paste(missing_vars, collapse = ", "),
      call. = FALSE
    )
  }

  response <- dataset[[response_name]]

  if (!is.numeric(response) && !is.integer(response)) {
    stop("ERROR: response variable must be numeric for gamlss_longitudinal().", call. = FALSE)
  }

  if (!.gl_has_supported_base_column_class(response)) {
    stop(

      "ERROR: response variable has unsupported class: ",
      paste(class(response), collapse = "/"),
      ". Use an ordinary numeric or integer response column.",
      call. = FALSE
    )
  }

  if (any(is.nan(response) | is.infinite(response))) {
    stop("ERROR: response variable contains NaN or Inf values; only finite values or NA are allowed.", call. = FALSE)
  }

  predictor_vars <- setdiff(formula_vars, response_name)
  unsupported <- character(0)
  nonfinite <- character(0)
  missing <- character(0)
  character_predictors <- character(0)

  for (nm in predictor_vars) {
    col <- dataset[[nm]]
    supported <- is.numeric(col) || is.integer(col) || is.logical(col) ||
      is.factor(col) || is.character(col)
    if (!supported || !.gl_has_supported_base_column_class(col)) {
      unsupported <- c(unsupported, nm)
      next
    }

    if (is.numeric(col) || is.integer(col)) {
      if (any(is.nan(col) | is.infinite(col))) {
        nonfinite <- c(nonfinite, nm)
      }
      if (any(is.na(col))) {
        missing <- c(missing, nm)
      }
    } else if (is.factor(col) || is.character(col) || is.logical(col)) {
      if (any(is.na(col))) {
        missing <- c(missing, nm)
      }
      if (is.character(col)) {
        character_predictors <- c(character_predictors, nm)
      }
    }
  }

  if (length(unsupported) > 0L) {
    stop(
      "ERROR: predictor column(s) have unsupported classes: ",
      paste(unsupported, collapse = ", "),
      ". Use numeric, integer, logical, factor, ordered factor, or character columns.",
      call. = FALSE
    )
  }

  if (length(nonfinite) > 0L) {
    stop(
      "ERROR: predictor column(s) contain NaN or Inf values: ",
      paste(unique(nonfinite), collapse = ", "),
      ". Clean non-finite predictor values before fitting.",
      call. = FALSE
    )
  }

  if (length(missing) > 0L) {
    stop(
      "ERROR: predictor column(s) contain missing values in submitted rows: ",
      paste(unique(missing), collapse = ", "),
      ". Missing responses and structurally missing visits are allowed, but predictor values must be observed.",
      call. = FALSE
    )
  }

  if (length(character_predictors) > 0L) {
    warning(
      "Character predictor column(s) will be treated as unordered categorical variables: ",
      paste(unique(character_predictors), collapse = ", "),
      ". Convert to factor to control level ordering explicitly.",
      call. = FALSE
    )
  }

  invisible(NULL)
}
