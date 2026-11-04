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
    fixed_names <- if (ncol(X_fixed)) paste(pn, colnames(mm$x[[pn]]), sep = ".") else character()
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

# The same cached constrained-basis penalty is used by WLS, block acceptance,
# and complete fixed-lambda sweeps. K belongs to GAIC, not this quadratic.
.gl_rs_penalty_value <- function(par_s, lambda_s, design_cache) {
  value <- 0
  for (pn in names(design_cache)) {
    for (sn in names(design_cache[[pn]]$smooth_penalty_meta)) {
      meta <- design_cache[[pn]]$smooth_penalty_meta[[sn]]
      ids <- colnames(design_cache[[pn]]$X)[meta$idx]
      b <- par_s[[pn]][[sn]]
      lambda <- lambda_s[[pn]][[sn]]
      if (!is.numeric(b) || !identical(names(b), ids) ||
          !is.numeric(lambda) || length(lambda) != 1L ||
          !is.matrix(meta$S_base) || !identical(dim(meta$S_base), c(length(b), length(b)))) {
        stop("Invalid RS coefficient/penalty coordinates.", call. = FALSE)
      }
      if (any(!is.finite(b)) || !is.finite(lambda) || lambda < 0 ||
          any(!is.finite(meta$S_base))) return(NA_real_)
      value <- value + 0.5 * lambda * sum(b * as.numeric(meta$S_base %*% b))
    }
  }
  as.numeric(value)
}

.gl_rs_objective <- function(likelihood, par_s, lambda_s, design_cache) {
  as.numeric(likelihood$log_lik["joint"]) -
    .gl_rs_penalty_value(par_s, lambda_s, design_cache)
}

.gl_rs_same_lambda <- function(x, y) {
  numeric_values <- function(z) lapply(z, function(block) lapply(block, as.numeric))
  identical(numeric_values(x), numeric_values(y))
}
