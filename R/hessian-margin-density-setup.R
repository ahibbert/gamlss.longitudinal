#' Prepare density function and arguments for margin log-density finite differences
#'
#' @noRd
.hessian_margin_density_setup <- function(eta_inv, mm, margin_dist, response) {
  margin_pars <- names(mm$x)[names(mm$x) %in% c("mu", "sigma", "nu", "tau")]
  dfun_name <- paste0("d", margin_dist$family[1])
  dfun <- eval(parse(text = dfun_name))

  base_args <- list(x = response, log = TRUE)
  for (pn in margin_pars) {
    base_args[[pn]] <- as.numeric(eta_inv[[pn]])
  }

  list(
    margin_pars = margin_pars,
    dfun = dfun,
    base_args = base_args
  )
}
