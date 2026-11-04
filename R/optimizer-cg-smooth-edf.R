.gl_cg_smooth_edf_list <- function(H_obs_current, penalty_current, beta_names, par_s) {
  edf_out <- setNames(lapply(names(par_s), function(x) list()), names(par_s))
  for (pn in names(par_s)) {
    if (length(par_s[[pn]]) == 0) next
    for (sn in names(par_s[[pn]])) {
      idx <- match(names(par_s[[pn]][[sn]]), rownames(H_obs_current))
      idx <- idx[!is.na(idx)]
      if (length(idx) == 0) next
      H_block <- H_obs_current[idx, idx, drop = FALSE]
      P_block <- penalty_current[idx, idx, drop = FALSE]
      info_block <- -0.5 * (H_block + t(H_block))
      if (sum(diag(info_block), na.rm = TRUE) < 0) {
        info_block <- -info_block
      }
      info_block <- tryCatch(
        {
          eg <- eigen(0.5 * (info_block + t(info_block)), symmetric = TRUE)
          eg$values[eg$values < 0] <- 0
          eg$vectors %*% diag(eg$values, nrow = length(eg$values)) %*% t(eg$vectors)
        },
        error = function(e) info_block
      )
      P_block <- 0.5 * (P_block + t(P_block))
      edf_val <- tryCatch(
        {
          k <- nrow(info_block)
          ridge <- max(1e-8, 1e-8 * max(1, max(abs(diag(info_block)), na.rm = TRUE)))
          sum(diag(.solve_linear_system(info_block + P_block + diag(ridge, k), info_block)))
        },
        error = function(e) NA_real_
      )
      if (!is.finite(edf_val)) edf_val <- length(idx)
      edf_out[[pn]][[sn]] <- max(0, min(length(idx), as.numeric(edf_val)))
    }
  }
  edf_out
}

#' Default source files needed by the CG analytical Hessian path
#'
