#' Legacy fixed-effect vcov derivative-matrix path
#'
#' This helper preserves the historical derivative-matrix implementation that
#' is not reached by the current public `method` choices. It is kept separate so
#' reviewers can follow the active analytical/numerical vcov workflow in
#' `model-vcov.R` without losing the old reference code.
#'
#' @noRd
.gl_vcov_legacy_derivative_matrix <- function(
    object,
    eta_inv,
    mm,
    margin_dist,
    response,
    response_margin,
    response_subject,
    margin_names,
    num_margins,
    margin_p,
    margin_d,
    copula_d,
    copula_derivatives,
    margin_derivatives,
    nd_impact_F,
    nd_impact_F2,
    eta,
    include_dlcopdpar,
    sep_d2) {
  ####### to delete############

  nd_impact_C2 <- calc_Fx2_derivatives(eta_inv, mm, margin_dist, response, testing = TRUE, response_margin, response_subject)[[2]]

  ### MARGIN LIKELIHOOD DERIVATIVES

  margin_deriv_subnames <- c("m", "d", "v", "t")

  names(margin_deriv_subnames) <- c("mu", "sigma", "nu", "tau")

  margin_par <- names(mm)[names(mm) %in% c("mu", "sigma", "nu", "tau")]

  order_margin <- cbind(object$response_margin, object$response_subject)

  colnames(order_margin) <- c("time", "subject")

  dcdu1 <- copula_derivatives$dcdu1
  dcdu2 <- copula_derivatives$dcdu2
  d2cdu12 <- copula_derivatives$d2cdu12
  d2cdu22 <- copula_derivatives$d2cdu22

  ##################### For each parameter...

  d1_cop <- d2_cop <- matrix(0, nrow = length(response), ncol = length(margin_par))

  colnames(d1_cop) <- colnames(d2_cop) <- margin_par

  pair_cache_diag <- build_copula_pair_cache(response, response_margin, response_subject)

  for (par_name in margin_par) {
    if (object$include_dlcopdpar == TRUE | include_dlcopdpar == TRUE) {
      order_copula <- data.frame()

      for (i in 1:(num_margins - 1)) {
        order_copula <- rbind(order_copula, cbind(order_margin[response_margin == margin_names[i], c("time", "subject")], order_margin[response_margin == margin_names[i + 1], c("time", "subject")]))
      }

      colnames(order_copula) <- c("time1", "subject1", "time2", "subject2")

      margin_deriv_1 <- matrix(0, ncol = length(margin_par), nrow = length(response))

      colnames(margin_deriv_1) <- paste("dld", margin_par, sep = "")

      margin_deriv_1[, paste("dld", par_name, sep = "")] <- margin_derivatives[grepl("dld", names(margin_derivatives))][[which(margin_par == par_name)]]

      # COPULA DERIVS WITH RESPECT TO

      mu <- eta_inv[["mu"]]

      F_nd <- nd_impact_F[[par_name]]

      F_nd2 <- nd_impact_F2[[par_name]]

      c_nd2 <- nd_impact_C2[[par_name]]

      margin_components <- cbind(order_margin, response, margin_p, margin_d, margin_deriv_1, mu, F_nd, F_nd2)

      margin_components_Ft_plus <- margin_components

      margin_components_Ft_plus[, "time"] <- normalize_lag_time(margin_components_Ft_plus[, "time"])

      margin_plus <- merge(margin_components, margin_components_Ft_plus, by = c("time", "subject"), all.x = TRUE)

      copula_components <- cbind(

        order_copula,
        row_id1 = pair_cache_diag$row_id1,
        row_id2 = pair_cache_diag$row_id2,
        dcdu1,
        dcdu2,
        copula_d,
        d2cdu12,
        d2cdu22,
        c_nd2
      )

      copula_merged <- merge(copula_components, margin_plus, by.x = c("time1", "subject1"), by.y = c("time", "subject"), all.x = TRUE)

      # Calculate copula derivative with respect to marginal parameters

      input <- copula_merged

      d1_cop[, par_name] <- calc_deriv_copula_wrt_margin(input, margin_par, par_name, calc_d2 = FALSE)[, which(margin_par == par_name)]

      # OK so let's calcute the numerical d2lcopdpar and pass it through input

      d2_cop[, par_name] <- calc_deriv_copula_wrt_margin(input, margin_par, par_name, calc_d2 = TRUE)[, which(margin_par == par_name)]
    }
  }

  ########### Need d1 and d2 for score function

  m_d1_names <- names(margin_derivatives)[grepl("dld", names(margin_derivatives))]

  c_d1_names <- names(copula_derivatives)[grepl("dld", names(copula_derivatives))]

  m_d2_names <- names(margin_derivatives)[grepl("d2ld", names(margin_derivatives))]

  c_d2_names <- names(copula_derivatives)[grepl("d2ld", names(copula_derivatives))]

  d1_all <- list()
  d2_all <- list()

  i <- 1

  for (par_name in c(m_d1_names)) {
    d1_all[[par_name]] <- c(margin_derivatives[[par_name]])

    if ((object$include_dlcopdpar == TRUE | include_dlcopdpar == TRUE)) {
      d1_all[[par_name]] <- c(margin_derivatives[[par_name]] + d1_cop[, i])

      i <- i + 1
    }
  }

  for (par_name in c(c_d1_names)) {
    d1_all[[par_name]] <- c(copula_derivatives[[par_name]])
  }

  names(d1_all) <- names(mm)

  i <- 1

  for (par_name in c(m_d2_names)) {
    d2_all[[par_name]] <- c(margin_derivatives[[par_name]])

    if ((object$include_dlcopdpar == TRUE | include_dlcopdpar == TRUE) & endsWith(par_name, "2")) {
      d2_all[[par_name]] <- c(margin_derivatives[[par_name]] + d2_cop[, i] * (if (sep_d2 == TRUE) {
        0
      } else {
        1
      }))

      i <- i + 1
    }
  }

  for (par_name in c(c_d2_names)) {
    d2_all[[par_name]] <- c(copula_derivatives[[par_name]])
  }

  d2_all_mean <- rep(0, length = length(d2_all))

  names(d2_all_mean) <- names(d2_all)

  for (deriv_name in names(d2_all)) {
    d2_all_mean[deriv_name] <- mean(d2_all[[deriv_name]])
  }

  d2_mat_diag <- d2_all_mean[endsWith(names(d2_all_mean), "2")]

  d2_mat_cross <- d2_all_mean[!endsWith(names(d2_all_mean), "2")]

  d2_mat <- matrix(nrow = length(eta), ncol = length(eta))

  # print(d2_mat);print(names(eta))

  colnames(d2_mat) <- rownames(d2_mat) <- names(eta)

  copula_deriv_subnames <- c("th", "z")

  names(copula_deriv_subnames) <- c("theta", "zeta")

  all_names <- c(margin_deriv_subnames, copula_deriv_subnames)

  sub_names_in <- all_names[names(eta)]

  print(sub_names_in)

  for (row_name in rownames(d2_mat)) {
    for (col_name in colnames(d2_mat)) {
      if (!row_name == col_name) {
        deriv_name_temp <- paste("d2ld", sub_names_in[row_name], "d", sub_names_in[col_name], sep = "")

        if (is.na(d2_all_mean[deriv_name_temp])) {
          deriv_name_temp <- paste("d2ld", sub_names_in[col_name], "d", sub_names_in[row_name], sep = "")
        }

        deriv_val_temp <- d2_all_mean[deriv_name_temp]

        d2_mat[row_name, col_name] <- deriv_val_temp
      }
    }
  }

  # cop_row=(grepl("theta",rownames(d2_mat))|grepl("zeta",rownames(d2_mat)))

  # d2_mat[!cop_row,!cop_row][upper.tri(d2_mat[!cop_row,!cop_row])]=d2_mat_cross

  diag(d2_mat) <- d2_mat_diag

  d2_mat[is.na(d2_mat)] <- 0

  d2_mat
}
