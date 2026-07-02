#' Build RS optimizer design matrix
#'
#' Combines fixed-effect and smooth design matrices for each parameter and
#' records smooth penalty metadata used by the RS backfitting loop.
#'
#' @noRd
.gl_build_rs_design_cache <- function(mm, par_s) {
  rs_design_cache <- setNames(vector("list", length(names(mm$x))), names(mm$x))

  for (pn in names(mm$x)) {
    X_fixed <- as.matrix(mm$x[[pn]])
    fixed_names <- paste(pn, colnames(mm$x[[pn]]), sep = ".")
    X_parts <- list(X_fixed)
    smooth_penalty_meta <- list()

    if (length(mm$s[[pn]]) > 0) {
      start_idx <- ncol(X_fixed) + 1L
      for (s_name in names(mm$s[[pn]])) {
        B <- as.matrix(mm$s[[pn]][[s_name]])
        smooth_names <- names(par_s[[pn]][[s_name]])
        colnames(B) <- smooth_names
        X_parts[[length(X_parts) + 1L]] <- B

        n_B <- ncol(B)
        idx <- start_idx:(start_idx + n_B - 1L)
        pen_attr <- attr(mm$s[[pn]][[s_name]], "penalty")
        if (!is.null(pen_attr) && is.matrix(pen_attr) &&
          nrow(pen_attr) == n_B && ncol(pen_attr) == n_B) {
          S_base <- pen_attr
        } else {
          D <- diff(diag(n_B), differences = 2)
          S_base <- t(D) %*% D
        }

        smooth_penalty_meta[[s_name]] <- list(idx = idx, B = B, S_base = S_base)
        start_idx <- start_idx + n_B
      }
    }

    X_combined <- do.call(cbind, X_parts)
    colnames(X_combined)[seq_along(fixed_names)] <- fixed_names
    rs_design_cache[[pn]] <- list(
      X = X_combined,
      fixed_names = fixed_names,
      smooth_penalty_meta = smooth_penalty_meta
    )
  }

  rs_design_cache
}
