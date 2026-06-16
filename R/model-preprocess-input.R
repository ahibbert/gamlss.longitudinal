#' Standardises names and validates input for time, subject, and response variables, 
#' and translates user-specified formulas to internal format.
#'
#' @keywords internal
#' @noRd
.gl_normalize_fit_input_columns <- function(
    dataset,
    time_var,
    subject_var,
    mu.formula,
    verbose = 1) {
  dataset_original <- dataset

  if (!is.data.frame(dataset)) {
    dataset <- as.data.frame(dataset, stringsAsFactors = FALSE)
  } else {
    dataset <- as.data.frame(dataset, stringsAsFactors = FALSE)
  }

  .gl_validate_tabular_shape(dataset, context = "dataset")

  if (all(is.na(time_var)) || all(is.na(subject_var))) {
    stop(
      "ERROR: Required input variables not specified.\n",
      "Please specify:\n",
      "  - time_var: column name for time/margin variable (e.g., 'time')\n",
      "  - subject_var: column name for subject ID variable (e.g., 'subject')\n",
      "Example: gamlss_longitudinal(..., time_var='time', subject_var='subject')"
    )
  }

  if (!time_var %in% colnames(dataset)) {
    stop(
      "ERROR: time_var='", time_var, "' not found in dataset.\n",
      "Available columns: ", paste(colnames(dataset), collapse = ", ")
    )
  }
  if (!subject_var %in% colnames(dataset)) {
    stop(
      "ERROR: subject_var='", subject_var, "' not found in dataset.\n",
      "Available columns: ", paste(colnames(dataset), collapse = ", ")
    )
  }

  mu_formula_obj <- as.formula(mu.formula)
  response_var <- all.vars(mu_formula_obj)[1]
  if (verbose > 1) print(paste("Identified response variable:", response_var))

  if (!response_var %in% names(dataset)) {
    stop(
      "ERROR: response variable '", response_var, "' not found in dataset.\n",
      "Available columns: ", paste(names(dataset), collapse = ", ")
    )
  }

  names(dataset)[names(dataset) == time_var] <- "time"
  names(dataset)[names(dataset) == subject_var] <- "subject"
  names(dataset)[names(dataset) == response_var] <- "response"

  time_subject <- .gl_prepare_time_subject_columns(dataset, time_var)

  var_map <- c()
  var_map[[time_var]] <- "time"
  var_map[[subject_var]] <- "subject"
  var_map[[response_var]] <- "response"

  formula_var_map <- var_map
  formula_var_map[[time_var]] <- "time_covariate"

  list(
    dataset_original = dataset_original,
    dataset = time_subject$dataset,
    response_var = response_var,
    var_map = var_map,
    formula_var_map = formula_var_map,
    time_covariate_is_factor = time_subject$time_covariate_is_factor,
    time_covariate_levels = time_subject$time_covariate_levels,
    time_covariate_ordered = time_subject$time_covariate_ordered
  )
}
