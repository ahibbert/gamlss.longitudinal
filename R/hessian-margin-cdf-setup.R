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

  base_cdf <- calc_F_x(eta_inv, NULL, margin_dist, response)
  use_survival <- is.finite(response) & is.finite(base_cdf) & base_cdf > 0.5 &
    !.is_discrete_margin(margin_dist) & "lower.tail" %in% names(formals(pfun))

  list(
    margin_pars = margin_pars,
    use_survival = use_survival,
    pfun = pfun,
    args_base = args_base[names(args_base) %in% names(formals(pfun))]
  )
}

#' Evaluate a fixed tail representation for CDF finite differences
#' @noRd
.hessian_margin_cdf_value <- function(eta_inv, mm, margin_dist, response, setup) {
  value <- calc_F_x(eta_inv, mm, margin_dist, response)
  if (any(setup$use_survival)) {
    input <- .gl_margin_cdf_input(eta_inv, margin_dist, response)
    input$lower.tail <- FALSE
    survival <- .gl_evaluate_margin_cdf(input, margin_dist, response)
    # F = 1 - S, so every positive-order derivative of F equals that of -S.
    # Keep the base-point tail choice fixed across all perturbations.
    value[setup$use_survival] <- -survival[setup$use_survival]
  }
  value
}
