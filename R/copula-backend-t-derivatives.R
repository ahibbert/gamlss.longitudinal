.copula_t_deriv_many <- function(u1, u2, par, par2, derivs, log = FALSE) {
  derivs <- unique(as.character(derivs))
  cmp <- .copula_t_components(u1, u2, par, par2)
  log_for_deriv <- function(deriv) {
    if (length(log) == 1L && is.null(names(log))) {
      return(isTRUE(log))
    }
    isTRUE(log[[deriv]])
  }

  out <- list()
  for (deriv in derivs) {
    score <- switch(deriv,
      u1 = {
        dlog_dx <- -(cmp$df + 2) * (cmp$x - cmp$rho * cmp$y) / cmp$denom +
          (cmp$df + 1) * cmp$x / (cmp$df + cmp$x^2)
        dlog_dx / stats::dt(cmp$x, df = cmp$df)
      },
      u2 = {
        dlog_dy <- -(cmp$df + 2) * (cmp$y - cmp$rho * cmp$x) / cmp$denom +
          (cmp$df + 1) * cmp$y / (cmp$df + cmp$y^2)
        dlog_dy / stats::dt(cmp$y, df = cmp$df)
      },
      par = {
        dq_drho <- (-2 * cmp$x * cmp$y * cmp$rho_denom + 2 * cmp$rho * cmp$numerator) / cmp$rho_denom^2
        cmp$rho / cmp$rho_denom - (cmp$df + 2) * dq_drho / (2 * (cmp$df + cmp$q))
      },
      par2 = {
        h <- .copula_t_df_step(cmp$df)
        dcd_df <- (
          .copula_t_pdf(cmp$u1, cmp$u2, cmp$rho, cmp$df + h) -
            .copula_t_pdf(cmp$u1, cmp$u2, cmp$rho, cmp$df - h)
        ) / (2 * h)
        dcd_df / cmp$density
      },
      stop("Unsupported t derivative: ", deriv, call. = FALSE)
    )

    value <- if (log_for_deriv(deriv)) score else cmp$density * score
    value[!is.finite(value)] <- 0
    out[[deriv]] <- value
  }

  attr(out, "density") <- cmp$density
  out
}

.copula_t_deriv <- function(u1, u2, par, par2, deriv, log = FALSE) {
  .copula_t_deriv_many(u1, u2, par, par2, derivs = deriv, log = log)[[deriv]]
}
