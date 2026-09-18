#' Compute the effective degrees of freedom for one RS smooth
#'
#' @noRd
.gl_rs_smooth_term_edf <- function(B, weights, penalty, lambda) {
  if (!is.matrix(B) || !nrow(B) || !ncol(B)) {
    stop("RS smooth EDF requires a non-empty smooth design matrix.", call. = FALSE)
  }
  weights <- as.numeric(weights)
  if (length(weights) != nrow(B) || any(!is.finite(weights))) {
    stop("RS smooth EDF requires one finite working weight per design row.", call. = FALSE)
  }
  if (!is.matrix(penalty) || !identical(dim(penalty), c(ncol(B), ncol(B))) ||
      any(!is.finite(penalty))) {
    stop("RS smooth EDF requires a finite square penalty matrix matching the smooth basis.", call. = FALSE)
  }
  lambda <- as.numeric(lambda)[1L]
  if (!is.finite(lambda) || lambda < 0) {
    stop("RS smooth EDF requires a finite non-negative smoothing parameter.", call. = FALSE)
  }

  information <- crossprod(B, B * weights)
  edf <- sum(.solve_linear_system(information + lambda * penalty) * information)
  if (!is.finite(edf) || edf < -sqrt(.Machine$double.eps)) {
    stop("RS smooth EDF calculation returned an invalid value.", call. = FALSE)
  }
  max(0, as.numeric(edf))
}

#' Refresh all RS smooth EDFs from the final optimizer state
#'
#' An RS block can legitimately require no accepted coefficient update in the
#' last fit. Its working weights and fitted smooth remain valid, but the
#' iteration-local EDF value may then still be the zero used at initialization.
#' Information criteria must instead use EDFs recomputed from the final stored
#' weights, penalties, and smoothing parameters for every smooth block.
#'
#' @noRd
.gl_refresh_rs_smooth_edf <- function(df_s, mm_s, lambda_s, weights_final) {
  if (is.null(mm_s) || !length(mm_s)) return(df_s)
  if (is.null(df_s)) df_s <- list()

  for (par_name in names(mm_s)) {
    smooths <- mm_s[[par_name]]
    if (is.null(smooths) || !length(smooths)) next
    if (is.null(weights_final[[par_name]])) {
      stop("Final RS working weights are missing for smooth parameter '", par_name, "'.", call. = FALSE)
    }
    if (is.null(df_s[[par_name]])) df_s[[par_name]] <- list()

    for (smooth_name in names(smooths)) {
      B <- smooths[[smooth_name]]
      penalty <- attr(B, "penalty")
      lambda <- lambda_s[[par_name]][[smooth_name]]
      if (is.null(lambda)) {
        stop(
          "Final RS smoothing parameter is missing for '", par_name, ": ",
          smooth_name, "'.", call. = FALSE
        )
      }
      value <- .gl_rs_smooth_term_edf(
        B = B,
        weights = weights_final[[par_name]],
        penalty = penalty,
        lambda = lambda
      )
      df_s[[par_name]][[smooth_name]] <- stats::setNames(value, smooth_name)
    }
  }

  df_s
}
