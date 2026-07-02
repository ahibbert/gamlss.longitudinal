#' Aggregate pair-level copula scores to one RS copula parameter
#'
#' @noRd
.gl_rs_copula_parameter_score <- function(par_name, eta, response, calc_lik_out, copula_derivatives) {
  if (par_name == "theta") {
    derivative <- copula_derivatives$dldth
    col_name <- "dldtheta"
  } else if (par_name == "zeta") {
    derivative <- copula_derivatives$dldz
    col_name <- "dldzeta"
  } else {
    stop("Unexpected copula parameter in optimisation: ", par_name)
  }

  n_par <- length(eta[[par_name]])
  d1_full <- matrix(0, nrow = n_par, ncol = 1)
  row_id1 <- calc_lik_out$copula_row_id1

  if (length(row_id1) > 0) {
    if (n_par == length(response)) {
      par_idx <- row_id1
    } else {
      par_idx <- calc_lik_out$copula_theta_index_map[row_id1]
    }

    valid_idx <- is.finite(par_idx) & par_idx >= 1 & par_idx <= nrow(d1_full)

    if (any(valid_idx)) {
      d1_sum <- rowsum(derivative[valid_idx], par_idx[valid_idx], reorder = FALSE)
      d1_full[as.integer(rownames(d1_sum)), 1] <- d1_sum[, 1]
    }
  }

  colnames(d1_full) <- col_name
  as.matrix(d1_full)
}
