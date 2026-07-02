#' @keywords internal

#' @noRd

get_copula_dist <- function(copula_dist) {
  copula_dist <- .copula_family_code(copula_dist)


  if (copula_dist == "C") {
    copula_link <- list(log, exp, dloginv = exp)
    two_par_cop <- FALSE

    parameters <- c("theta")
  } else if (copula_dist == "F") {
    copula_link <- list(identity, identity, function(x) rep(1, length(x)))
    two_par_cop <- FALSE

    parameters <- c("theta")
  } else if (copula_dist == "J") {
    copula_link <- list(log_1plus, log_1plus_inv, dlog_1plus_inv)
    two_par_cop <- FALSE

    parameters <- c("theta")
  } else if (copula_dist == "G") {
    copula_link <- list(gumbel_linkfun, gumbel_linkinv, dgumbel_linkinv)
    two_par_cop <- FALSE

    parameters <- c("theta")
  } else if (copula_dist == "N") {
    copula_link <- list(fisher_z, fisher_z_inv, dfisher_z_inv)
    two_par_cop <- FALSE

    parameters <- c("theta")
  } else if (copula_dist == "t") {
    copula_link <- list(fisher_z, fisher_z_inv, dfisher_z_inv, log_2plus, log_2plus_inv, dlog_2plus_inv)
    two_par_cop <- TRUE

    parameters <- c("theta", "zeta")
  } else {
    stop("ERROR: COPULA DIST LINK FUNCTIONS NOT YET IMPLEMENTED.")
  }


  if (two_par_cop) {
    names(copula_link) <- c("theta.linkfun", "theta.linkinv", "theta.dr", "zeta.linkfun", "zeta.linkinv", "zeta.dr")
  } else {
    names(copula_link) <- c("theta.linkfun", "theta.linkinv", "theta.dr")
  }


  return_list <- list()

  return_list[["copula_link"]] <- copula_link

  return_list[["copula_dist"]] <- copula_dist

  return_list[["parameters"]] <- parameters


  return(return_list)
}
