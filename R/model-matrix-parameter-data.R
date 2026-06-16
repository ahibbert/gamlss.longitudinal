#' Select rows of model matrix used for theta and zeta:
#' here we are selecting just the first margin of the data for each time point, 
#' since theta and zeta are time-varying but not subject-varying we have to select 
#' one row per time point if they differ in value. Essentially the value for the earlier edge 
#' is relied on for the model.
#'
#' @keywords internal
#' @noRd
.gl_model_matrix_parameter_dataset <- function(dataset_mm, parameter) {
  if (parameter %in% c("theta", "zeta")) {
    times <- unique(dataset_mm$time)
    return(dataset_mm[dataset_mm$time %in% times[1:(length(times) - 1)], , drop = FALSE])
  }

  dataset_mm
}
