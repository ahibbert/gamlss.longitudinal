#' Prepare fitting data for gamlss_longitudinal()
#'
#' Renames user-supplied response, time, and subject columns to the internal names used by the fitting engine, 
#' validates submitted data against various policies, expands structurally missing subject/time combinations 
#' (for missingness checks and model matrix construction), and returns translated formulas for matrix construction.
#'
#' @noRd
.gl_prepare_fit_data <- function(
    dataset,
    time_var,
    subject_var,
    mu.formula,
    sigma.formula,
    nu.formula,
    tau.formula,
    theta.formula,
    zeta.formula,
    verbose = 1) {

  input_columns <- .gl_normalize_fit_input_columns(
    dataset = dataset,
    time_var = time_var,
    subject_var = subject_var,
    mu.formula = mu.formula,
    verbose = verbose
  )
  dataset_original <- input_columns$dataset_original
  dataset <- input_columns$dataset
  response_var <- input_columns$response_var
  var_map <- input_columns$var_map
  formula_var_map <- input_columns$formula_var_map
  time_covariate_is_factor <- input_columns$time_covariate_is_factor
  time_covariate_levels <- input_columns$time_covariate_levels
  time_covariate_ordered <- input_columns$time_covariate_ordered

  formulas_int <- .gl_translate_fit_formulas(
    mu.formula = mu.formula,
    sigma.formula = sigma.formula,
    nu.formula = nu.formula,
    tau.formula = tau.formula,
    theta.formula = theta.formula,
    zeta.formula = zeta.formula,
    formula_var_map = formula_var_map,
    response_name = "response"
  )

  mu.formula.int <- formulas_int$mu
  sigma.formula.int <- formulas_int$sigma
  nu.formula.int <- formulas_int$nu
  tau.formula.int <- formulas_int$tau
  theta.formula.int <- formulas_int$theta
  zeta.formula.int <- formulas_int$zeta

  .gl_validate_fitting_data_policy(
    dataset,
    formulas = list(mu.formula.int, sigma.formula.int, nu.formula.int, tau.formula.int, theta.formula.int, zeta.formula.int),
    response_name = "response"
  )

  if (verbose > 1) {
    cat("Input validation successful.\n")
    cat("Data dimensions:", nrow(dataset), "x", ncol(dataset), "\n")
    cat("Response variable:", response_var, "-> renamed to 'response'\n")
    cat("Time variable:", time_var, "-> internal index 'time' and covariate 'time_covariate'\n")
    cat("Subject variable:", subject_var, "-> renamed to 'subject'\n")
    cat("Time points:", length(unique(dataset$time)), "\n")
    cat("Subjects:", length(unique(dataset$subject)), "\n")
  }

  # Validate that all subject/time combinations are unique
  subject_time_combo <- paste(dataset$subject, dataset$time, sep = "_")
  if (length(subject_time_combo) != length(unique(subject_time_combo))) {
    duplicate_combos <- subject_time_combo[duplicated(subject_time_combo)]
    stop(
      "ERROR: Duplicate subject/time combinations found.\n",
      "Each subject must have exactly one observation per time point.\n",
      "Duplicate combinations (first 10): ",
      paste(unique(duplicate_combos)[1:min(10, length(unique(duplicate_combos)))], collapse = ", ")
    )
  }

  if (verbose > 1) {
    cat("Subject/time uniqueness check passed.\n")
    cat("Unique subject/time combinations:", length(unique(subject_time_combo)), "\n\n")
  }

  expanded_panel <- .gl_expand_fit_panel(
    dataset = dataset,
    time_covariate_is_factor = time_covariate_is_factor,
    time_covariate_levels = time_covariate_levels,
    time_covariate_ordered = time_covariate_ordered,
    verbose = verbose
  )

  dataset <- expanded_panel$dataset
  missingness <- .gl_summarize_fit_missingness(dataset, verbose = verbose)
  miss_by_time <- missingness$miss_by_time
  pair_summary <- missingness$pair_summary
  .gl_validate_fit_missingness_support(miss_by_time, pair_summary)

  list(
    dataset_original = dataset_original,
    dataset = dataset,
    response_var = response_var,
    var_map = var_map,
    formulas_int = list(
      mu = mu.formula.int,
      sigma = sigma.formula.int,
      nu = nu.formula.int,
      tau = tau.formula.int,
      theta = theta.formula.int,
      zeta = zeta.formula.int
    ),
    miss_by_time = miss_by_time,
    pair_summary = pair_summary
  )
}
