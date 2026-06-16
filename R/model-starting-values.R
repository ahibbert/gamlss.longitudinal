#' Get starting values for copula and margin parameters
#' We use Kendall's tau to copula parameter conversion for copula starting values, 
#' and either method-of-moments or a brief gamlss fit for margin starting values. 
#' 
#' For t-copulas, we select the zeta starting value by a small grid search over candidate zeta values which can be a bit slow.
#' 
#' @keywords internal
#' @noRd
get_starting_values <- function(copula_dist, margin_dist, dataset, eta_transform = FALSE) {
  margin_dist <- .normalise_margin_dist_links(margin_dist)
  margin_names <- unique(dataset$time)
  num_margins <- length(margin_names)
  finite_response <- dataset$response[is.finite(dataset$response)]

  tau_start <- cor(dataset[dataset$time %in% (margin_names[1:(num_margins - 1)]), "response"],
    dataset[dataset$time %in% (margin_names[2:(num_margins)]), "response"],
    method = "kendall", use = "complete.obs"
  )

  if (!is.finite(tau_start)) {
    warning("Non-finite Kendall tau in get_starting_values(); using tau = 0 for copula initialisation.")
    tau_start <- 0
  }

  tau_start <- max(min(tau_start, 0.9999), -0.9999)
  copula_spec <- get_copula_dist(copula_dist)
  theta_start <- .copula_tau_to_par(
    family = copula_dist,
    tau = tau_start
  )
  margin_start <- .starting_margin_parameter_values(
    margin_dist = margin_dist,
    finite_response = finite_response,
    dataset = dataset
  )

  margin_par <- margin_start$margin_par
  margin_par_already_eta <- margin_start$margin_par_already_eta

  if ("zeta" %in% copula_spec$parameters) {
    # .copula_tau_to_par() returns only theta for t-copula; select zeta by a small grid search.
    zeta_start <- .select_t_copula_zeta_start(
      dataset = dataset,
      margin_dist = margin_dist,
      copula_dist = copula_dist,
      margin_par = margin_par,
      theta_start = theta_start
    )
    cop_par <- c(theta = as.numeric(theta_start)[1], zeta = as.numeric(zeta_start))
  } else {
    cop_par <- c(theta = as.numeric(theta_start)[1])
  }

  if (eta_transform == TRUE) {
    margin_par_eta <- margin_par
    cop_par_eta <- cop_par

    if (!isTRUE(margin_par_already_eta)) {
      for (par_name in names(margin_par)) {
        FUN <- eval(parse(text = paste(paste(paste("margin_dist$", par_name, sep = ""), "linkfun", sep = "."))))
        margin_par_eta[par_name] <- FUN(margin_par[par_name])
      }
    }

    for (par_name in names(cop_par)) {
      cop_par_eta[par_name] <- get_copula_dist(copula_dist)$copula_link[[paste(par_name, ".linkfun", sep = "")]](cop_par[par_name])
    }

    return_list <- c(margin_par_eta, cop_par_eta)
  } else {
    return_list <- c(margin_par, cop_par)
  }

  return(return_list)
}
