#' First derivative of the log-copula term with respect to margin parameters
#'
#' @keywords internal
#' @noRd
.calc_deriv_copula_wrt_margin_d1 <- function(input, margin_par, par_name) {
  dlcopdpar_row1 <- matrix(0, nrow = nrow(input), ncol = length(margin_par))
  dlcopdpar_row2 <- matrix(0, nrow = nrow(input), ncol = length(margin_par))
  i <- 1
  for (inner_par_name in margin_par) {
    if (inner_par_name == par_name) {
      # Take parameters from input for clarity
      dc_tplus_du_t <- .copula_margin_derivative_numeric_column(input, "dcdu1")
      dc_tplus_du_tplus <- .copula_margin_derivative_numeric_column(input, "dcdu2")
      l_t <- .copula_margin_derivative_numeric_column(input, paste(paste("dld", inner_par_name, sep = ""), ".x", sep = ""))
      l_t_plus <- .copula_margin_derivative_numeric_column(input, paste(paste("dld", inner_par_name, sep = ""), ".y", sep = ""))
      x_t <- .copula_margin_derivative_numeric_column(input, "response.x")
      x_t_plus <- .copula_margin_derivative_numeric_column(input, "response.y")
      f_t <- .copula_margin_derivative_numeric_column(input, "margin_d.x")
      f_t_plus <- .copula_margin_derivative_numeric_column(input, "margin_d.y")
      c_tplus <- .copula_margin_derivative_numeric_column(input, "copula_d")
      mu_t <- .copula_margin_derivative_numeric_column(input, "mu.x")
      mu_t_plus <- .copula_margin_derivative_numeric_column(input, "mu.y")

      F_nd_t <- .copula_margin_derivative_numeric_column(input, "F_nd.x")
      F_nd_t_plus <- .copula_margin_derivative_numeric_column(input, "F_nd.y")

      du_t_dmu <- F_nd_t
      du_t_plus_dmu <- F_nd_t_plus

      # Exact endpoint attribution for pair log-copula derivative:
      # row_id1 gets (dc/du1)*(du1/dpar)/c, row_id2 gets (dc/du2)*(du2/dpar)/c.
      dlogc_row1 <- (dc_tplus_du_t * du_t_dmu) / c_tplus
      dlogc_row2 <- (dc_tplus_du_tplus * du_t_plus_dmu) / c_tplus
      dlogc_row1[!is.finite(dlogc_row1)] <- 0
      dlogc_row2[!is.finite(dlogc_row2)] <- 0

      dlcopdpar_row1[, i] <- dlogc_row1
      dlcopdpar_row2[, i] <- dlogc_row2
    }
    i <- i + 1
  }
  colnames(dlcopdpar_row1) <- paste("dlcopd", margin_par, sep = "")
  colnames(dlcopdpar_row2) <- paste("dlcopd", margin_par, sep = "")

  # Prefer explicit row-index accumulation to avoid merge-order instability.
  if (all(c("row_id1", "row_id2") %in% colnames(input))) {
    n_obs <- suppressWarnings(max(c(
      .copula_margin_derivative_numeric_column(input, "row_id1"),
      .copula_margin_derivative_numeric_column(input, "row_id2")
    ), na.rm = TRUE))
    if (!is.finite(n_obs) || n_obs < 1) {
      stop("Invalid row ids in copula-to-margin derivative assembly.")
    }
    d1_cop <- matrix(0, nrow = as.integer(n_obs), ncol = length(margin_par))
    colnames(d1_cop) <- margin_par

    row_id1 <- as.integer(.copula_margin_derivative_numeric_column(input, "row_id1"))
    row_id2 <- as.integer(.copula_margin_derivative_numeric_column(input, "row_id2"))
    for (j in seq_along(margin_par)) {
      contrib1 <- as.numeric(dlcopdpar_row1[, j])
      contrib2 <- as.numeric(dlcopdpar_row2[, j])
      valid1 <- is.finite(contrib1) & is.finite(row_id1) & row_id1 >= 1 & row_id1 <= n_obs
      valid2 <- is.finite(contrib2) & is.finite(row_id2) & row_id2 >= 1 & row_id2 <= n_obs
      d1_cop[row_id1[valid1], j] <- d1_cop[row_id1[valid1], j] + contrib1[valid1]
      d1_cop[row_id2[valid2], j] <- d1_cop[row_id2[valid2], j] + contrib2[valid2]
    }
    return(d1_cop)
  }

  dlcopdpar <- dlcopdpar_row1 + dlcopdpar_row2

  par_dlcopdpar <- dlcopdpar[, paste("dlcopd", margin_par, sep = "")]
  merged_dlcopdpar <- merge(cbind(input[, c("time1", "time2", "subject1", "subject2")], par_dlcopdpar), cbind(input[, c("time1", "time2", "subject1", "subject2")], par_dlcopdpar), by.x = c("time2", "subject2"), by.y = c("time1", "subject1"), all = TRUE)
  merged_dlcopdpar[is.na(merged_dlcopdpar)] <- 0

  x_comp <- grepl("dlcopd", colnames(merged_dlcopdpar)) & grepl(".x", colnames(merged_dlcopdpar))
  y_comp <- grepl("dlcopd", colnames(merged_dlcopdpar)) & grepl(".y", colnames(merged_dlcopdpar))

  d1_cop <- 0.5 * (merged_dlcopdpar[, x_comp] + merged_dlcopdpar[, y_comp])

  return(d1_cop)
}
