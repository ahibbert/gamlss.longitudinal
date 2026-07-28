#' Apply CG trust-region and per-coefficient step limits
#'
#' @noRd
.gl_limit_cg_step <- function(delta, trust_radius, max_delta) {
  dnorm <- sqrt(sum(delta^2))
  if (is.finite(dnorm) && dnorm > trust_radius) delta <- delta * trust_radius / dnorm
  dc <- max(abs(delta))
  if (is.finite(dc) && dc > max_delta) delta <- delta * max_delta / dc
  delta
}

#' Shrink CG trust radius after a conservative event
#'
#' @noRd
.gl_shrink_cg_trust_radius <- function(trust_radius, step_tol) {
  max(trust_radius / 2, step_tol)
}

#' Expand CG trust radius after a boundary-sized accepted step
#'
#' @noRd
.gl_expand_cg_trust_radius <- function(trust_radius, step_l2, step_tol, max_delta) {
  if (is.finite(step_l2) && is.finite(trust_radius) && step_l2 >= 0.8 * trust_radius) {
    return(min(as.numeric(max_delta), max(step_tol, 1.5 * trust_radius)))
  }
  trust_radius
}

#' Build CG candidate steps for line search
#'
#' @noRd
.gl_build_cg_candidate_steps <- function(g_pen, H_pen, trust_radius) {
  candidate_steps <- list()
  grad_norm <- sqrt(sum(g_pen^2))
  if (is.finite(grad_norm) && grad_norm > 0) {
    candidate_steps[[length(candidate_steps) + 1L]] <- as.numeric(trust_radius * g_pen / grad_norm)
  }

  for (ridge in c(0, 1e-8, 1e-6, 1e-4, 1e-2, 1, 10, 100)) {
    d <- tryCatch(-as.numeric(.solve_linear_system(H_pen - diag(ridge, nrow(H_pen)), g_pen)), error = function(e) NULL)
    if (!is.null(d) && all(is.finite(d))) {
      candidate_steps[[length(candidate_steps) + 1L]] <- d
      candidate_steps[[length(candidate_steps) + 1L]] <- -d
    }
  }

  candidate_steps
}
