#' Prepare CDF function and arguments for margin finite differences
#'
#' @noRd
.hessian_margin_cdf_setup <- function(eta_inv, margin_dist, response) {
  margin_pars <- names(eta_inv)[names(eta_inv) %in% c("mu", "sigma", "nu", "tau")]
  pfun_name <- paste0("p", margin_dist$family[1])
  pfun <- tryCatch(eval(parse(text = pfun_name)), error = function(e) NULL)
  if (is.null(pfun)) stop("Cannot find CDF function: ", pfun_name)

  args_base <- list(q = response)
  for (pn in c("mu", "sigma", "nu", "tau")) {
    if (pn %in% names(eta_inv)) args_base[[pn]] <- eta_inv[[pn]]
  }
  args_base <- c(args_base, .gl_margin_fixed_family_args(margin_dist, length(response)))

  list(
    margin_pars = margin_pars,
    pfun = pfun,
    args_base = args_base[names(args_base) %in% names(formals(pfun))]
  )
}
