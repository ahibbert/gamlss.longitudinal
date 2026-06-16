#' Transforms eta scale parameter to parameter through link inverse for starting values.
#' 
#' @keywords internal
#' @noRd
eta_to_par <- function(eta, margin_dist, copula_dist) {
  par <- eta * 0

  for (par_name in names(eta)) {
    if (par_name %in% names(margin_dist$parameters)) {
      FUN <- eval(parse(text = paste(paste(paste("margin_dist$", par_name, sep = ""), "linkinv", sep = "."))))
      par[par_name] <- FUN(eta[par_name])
    }
    if (par_name %in% names(copula_dist$parameters)) {
      FUN <- eval(parse(text = paste(paste(paste("copula_dist$", par_name, sep = ""), "linkinv", sep = "."))))
      par[par_name] <- FUN(eta[par_name])
    }
  }

  return(par)
}

#' Transforms parameter starting values to eta transform through link for starting values. 
#' 
#' @keywords internal
#' @noRd
par_to_eta <- function(par, copula_dist, margin_dist) {
  margin_dist <- .normalise_margin_dist_links(margin_dist)
  margin_par <- par[names(margin_dist$parameters)]
  names(margin_par) <- names(margin_dist$parameters)

  cop_par <- par[get_copula_dist(copula_dist)$parameters]
  names(cop_par) <- get_copula_dist(copula_dist)$parameters

  margin_par_eta <- margin_par
  cop_par_eta <- cop_par

  for (par_name in names(margin_par)) {
    FUN <- eval(parse(text = paste(paste(paste("margin_dist$", par_name, sep = ""), "linkfun", sep = "."))))
    margin_par_eta[par_name] <- FUN(margin_par[par_name])
  }

  for (par_name in names(cop_par)) {
    cop_par_eta[par_name] <- get_copula_dist(copula_dist)$copula_link[[paste(par_name, ".linkfun", sep = "")]](cop_par[par_name])
    names(cop_par_eta) <- names(cop_par)
  }

  return_list <- c(margin_par_eta, cop_par_eta)

  return(return_list)
}
