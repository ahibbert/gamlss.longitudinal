calc_true_SE_numderiv_only <- function(eta_inv, mm, margin_dist, response, testing = FALSE, response_margin = NA, response_subject = NA) {
  adj_fac <- .001
  nd_impact <- rep(0, length(names(eta_inv)))
  names(nd_impact) <- margin_par_names <- names(eta_inv) # [names(eta_inv) %in% c("mu","sigma","nu","tau")]

  for (eta_par_names_nd in margin_par_names) {
    change <- rep(0, length(names(eta_inv)))
    i <- 1
    for (adj in c(-1 * adj_fac, adj_fac)) {
      eta_inv_adj <- eta_inv
      eta_inv_adj[[eta_par_names_nd]] <- eta_inv_adj[[eta_par_names_nd]] + adj
      change[i] <- calc_likelihood_minimal(eta_inv_adj, mm, margin_dist, copula_dist, calc_d2 = FALSE, response, response_margin, response_subject)$log_lik["joint"]
      i <- i + 1
    }
    change[3] <- calc_likelihood_minimal(eta_inv, mm, margin_dist, copula_dist, calc_d2 = FALSE, response, response_margin, response_subject)$log_lik["joint"]
    nd_impact[eta_par_names_nd] <- (change[2] + change[1] - 2 * change[3]) / (adj_fac^2)

    # print(c(change,nd_impact[eta_par_names_nd]))
  }

  nd_cross <- matrix(0, nrow = length(names(eta_inv)), ncol = length(names(eta_inv)))
  colnames(nd_cross) <- rownames(nd_cross) <- names(eta_inv)
  for (name1 in margin_par_names) {
    for (name2 in margin_par_names) {
      for (adj1 in c(-1 * adj_fac, adj_fac)) {
        for (adj2 in c(-1 * adj_fac, adj_fac)) {
          if (name1 != name2) {
            eta_inv_adj <- eta_inv
            eta_inv_adj[[name1]] <- eta_inv_adj[[name1]] + adj1
            eta_inv_adj[[name2]] <- eta_inv_adj[[name2]] + adj2
            change <- calc_likelihood_minimal(eta_inv_adj, mm, margin_dist, copula_dist, calc_d2 = FALSE, response, response_margin, response_subject)$log_lik["joint"]
            nd_cross[name1, name2] <- nd_cross[name1, name2] + change * if (adj1 == adj2) {
              1
            } else {
              -1
            }
          }
        }
      }
    }
  }
  nd_cross <- nd_cross / (4 * (adj_fac^2))

  nd2 <- (diag(nd_impact) + nd_cross)

  return(nd2)
}

calc_true_SE_numderiv_only_rowwise <- function(eta_inv, mm, margin_dist, response, testing = FALSE, response_margin = NA, response_subject = NA) {
  adj_fac <- .00001
  nd_impact_m <- nd_impact_c <- list()

  for (eta_par_names_nd in margin_par_names) {
    change_m <- change_c <- list()
    i <- 1
    for (adj in c(-1 * adj_fac, adj_fac)) {
      eta_inv_adj <- eta_inv
      eta_inv_adj[[eta_par_names_nd]] <- eta_inv_adj[[eta_par_names_nd]] + adj

      calc_lik_temp <- calc_likelihood_minimal(eta_inv_adj, mm, margin_dist, copula_dist, calc_d2 = FALSE, response, response_margin, response_subject)

      change_m[[i]] <- calc_lik_temp$margin_d
      change_c[[i]] <- calc_lik_temp$copula_d
      i <- i + 1
    }
    calc_lik_temp <- calc_likelihood_minimal(eta_inv, mm, margin_dist, copula_dist, calc_d2 = FALSE, response, response_margin, response_subject)
    change_m[[3]] <- calc_lik_temp$margin_d
    change_c[[3]] <- calc_lik_temp$copula_d
    nd_impact_m[[eta_par_names_nd]] <- (change_m[[3]] + change_m[[1]] - 2 * change_m[[3]]) / (adj_fac^2)
    nd_impact_c[[eta_par_names_nd]] <- (change_c[[3]] + change_c[[1]] - 2 * change_c[[3]]) / (adj_fac^2)
  }
  names(nd_impact_m) <- names(nd_impact_c) <- margin_par_names <- names(eta_inv)

  nd_cross_m
  colnames(nd_cross) <- rownames(nd_cross) <- names(eta_inv)
  for (name1 in margin_par_names) {
    for (name2 in margin_par_names) {
      for (adj1 in c(-1 * adj_fac, adj_fac)) {
        for (adj2 in c(-1 * adj_fac, adj_fac)) {
          if (name1 != name2) {
            eta_inv_adj <- eta_inv
            eta_inv_adj[[name1]] <- eta_inv_adj[[name1]] + adj1
            eta_inv_adj[[name2]] <- eta_inv_adj[[name2]] + adj2
            change <- calc_likelihood_minimal(eta_inv_adj, mm, margin_dist, copula_dist, calc_d2 = FALSE, response, response_margin, response_subject)$log_lik["joint"]
            nd_cross[name1, name2] <- nd_cross[name1, name2] + change * if (adj1 == adj2) {
              1
            } else {
              -1
            }
          }
        }
      }
    }
  }
  nd_cross <- nd_cross / (4 * (adj_fac^2))

  nd2 <- (diag(nd_impact) + nd_cross)

  return(nd2)
}
