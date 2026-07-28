calc_Fx_derivatives <- function(eta_inv, mm, margin_dist, response, par_names = NULL) {
  # Allow callers to pass full model matrix object; we only need fixed-effect blocks.
  if (is.list(mm) && all(c("x", "s") %in% names(mm))) {
    mm <- mm$x
  }

  nd_impact <- nd_impact_m <- nd_impact_c <- rep(0, length(names(eta_inv)))
  nd_impact_F <- list() ###### WE DON"T NEED TO BE CALCULATING THIS FOR ALL PARAMETERS EVERY TIME

  margin_par_names <- names(eta_inv)[names(eta_inv) %in% c("mu", "sigma", "nu", "tau")]
  if (!is.null(par_names)) {
    margin_par_names <- intersect(margin_par_names, par_names)
  }

  for (eta_par_names_nd in margin_par_names) {
    adj_fac <- .0001
    change <- change_m <- change_c <- c(0, 0)
    change_F <- matrix(0, nrow = length(eta_inv[[1]]), ncol = 2)
    i <- 1
    for (adj in c(-1 * adj_fac, adj_fac)) {
      eta_inv_adj <- eta_inv
      eta_inv_adj[[eta_par_names_nd]] <- eta_inv_adj[[eta_par_names_nd]] + adj
      # change[i]=calc_likelihood_minimal(eta_inv_adj,mm,margin_dist,copula_dist)$log_lik["joint"]
      # change_m[i]=calc_likelihood_minimal(eta_inv_adj,mm,margin_dist,copula_dist)$log_lik["marginal"]
      # change_c[i]=calc_likelihood_minimal(eta_inv_adj,mm,margin_dist,copula_dist)$log_lik["copula"]
      change_F[, i] <- calc_F_x(eta_inv_adj, mm, margin_dist, response) # calc_likelihood_minimal(eta_inv_adj,mm,margin_dist,copula_dist)$margin_p
      # print(change_F[,i]==calc_F_x(eta_inv_adj,mm,margin_dist))
      i <- i + 1
    }
    # nd_impact[eta_par_names_nd]=(change[2]-change[1])/(2*adj_fac)
    # nd_impact_m[eta_par_names_nd]=(change_m[2]-change_m[1])/(2*adj_fac)
    # nd_impact_c[eta_par_names_nd]=(change_c[2]-change_c[1])/(2*adj_fac)
    nd_impact_F[[eta_par_names_nd]] <- (change_F[, 2] - change_F[, 1]) / (2 * adj_fac)
  }
  return(nd_impact_F)
}

calc_Fx2_derivatives <- function(eta_inv, mm, margin_dist, response, testing = FALSE, response_margin = NA, response_subject = NA) {
  nd_impact <- nd_impact_m <- nd_impact_c <- rep(0, length(names(eta_inv)))
  nd_impact_F <- nd_impact_c <- list() ###### WE DON"T NEED TO BE CALCULATING THIS FOR ALL PARAMETERS EVERY TIME

  margin_par_names <- names(eta_inv)[names(eta_inv) %in% c("mu", "sigma", "nu", "tau")]

  for (eta_par_names_nd in margin_par_names) {
    adj_fac <- .0001
    change_F <- matrix(0, nrow = length(eta_inv[[1]]), ncol = 3)
    change_c <- matrix(0, nrow = length(eta_inv[["theta"]]), ncol = 3)
    i <- 1
    for (adj in c(-1 * adj_fac, adj_fac)) {
      eta_inv_adj <- eta_inv
      eta_inv_adj[[eta_par_names_nd]] <- eta_inv_adj[[eta_par_names_nd]] + adj
      # change[i]=calc_likelihood_minimal(eta_inv_adj,mm,margin_dist,copula_dist)$log_lik["joint"]
      # change_m[i]=calc_likelihood_minimal(eta_inv_adj,mm,margin_dist,copula_dist)$log_lik["marginal"]
      if (testing == TRUE) {
        change_c[, i] <- log(calc_likelihood_minimal(eta_inv_adj, mm, margin_dist, copula_dist, calc_d2 = FALSE, response, response_margin, response_subject)$copula_d)
      }
      change_F[, i] <- calc_F_x(eta_inv_adj, mm, margin_dist, response) # calc_likelihood_minimal(eta_inv_adj,mm,margin_dist,copula_dist)$margin_p
      # print(change_F[,i]==calc_F_x(eta_inv_adj,mm,margin_dist))
      i <- i + 1
    }
    change_F[, 3] <- calc_F_x(eta_inv, mm, margin_dist, response)
    if (testing == TRUE) {
      change_c[, 3] <- log(calc_likelihood_minimal(eta_inv, mm, margin_dist, copula_dist, calc_d2 = FALSE, response, response_margin, response_subject)$copula_d)
    }
    # nd_impact[eta_par_names_nd]=(change[2]-change[1])/(2*adj_fac)
    # nd_impact_m[eta_par_names_nd]=(change_m[2]-change_m[1])/(2*adj_fac)
    nd_impact_c[[eta_par_names_nd]] <- (change_c[, 2] + change_c[, 1] - 2 * change_c[, 3]) / (adj_fac^2)
    nd_impact_F[[eta_par_names_nd]] <- (change_F[, 2] + change_F[, 1] - 2 * change_F[, 3]) / (adj_fac^2)
  }
  if (testing == FALSE) {
    return(nd_impact_F)
  } else {
    return(list(nd_impact_F, nd_impact_c))
  }
}
