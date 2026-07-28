#' Validate optional simulated-uniform clipping bounds
#'
#' @noRd
.sim_validate_u_bounds <- function(u_bounds) {
  if (is.null(u_bounds)) {
    return(NULL)
  }

  if (!is.numeric(u_bounds) || length(u_bounds) != 2L ||
    any(!is.finite(u_bounds)) || u_bounds[1L] < 0 ||
    u_bounds[2L] > 1 || u_bounds[1L] >= u_bounds[2L]) {
    stop("u_bounds must be NULL or a finite increasing length-two vector inside [0, 1].", call. = FALSE)
  }

  u_bounds
}

#' Apply optional clipping bounds to simulated uniforms
#'
#' @noRd
.sim_apply_u_bounds <- function(u, u_bounds) {
  u_bounds <- .sim_validate_u_bounds(u_bounds)
  if (is.null(u_bounds)) {
    return(u)
  }

  pmin(pmax(u, u_bounds[1L]), u_bounds[2L])
}
