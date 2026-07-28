#' Finite-difference derivatives of marginal CDF bounds
#'
#' Computes observation-level derivatives of CDF bounds with respect to natural
#' margin parameters. These feed the discrete rectangle likelihood derivatives.
#'
#' @noRd
.calc_F_bounds_derivatives <- function(eta_inv, mm, margin_dist, response, par_names = NULL, h = 1e-4) {
  if (is.list(mm) && all(c("x", "s") %in% names(mm))) {
    mm <- mm$x
  }

  margin_par_names <- names(eta_inv)[names(eta_inv) %in% c("mu", "sigma", "nu", "tau")]
  if (!is.null(par_names)) {
    margin_par_names <- intersect(margin_par_names, par_names)
  }

  out <- setNames(vector("list", length(margin_par_names)), margin_par_names)
  for (pn in margin_par_names) {
    hp <- h * pmax(1, abs(as.numeric(eta_inv[[pn]])))
    eta_plus <- eta_minus <- eta_inv
    eta_plus[[pn]] <- eta_plus[[pn]] + hp
    eta_minus[[pn]] <- eta_minus[[pn]] - hp
    if (pn %in% c("mu", "sigma", "tau", "zeta")) {
      eta_plus[[pn]] <- pmax(eta_plus[[pn]], 1e-8)
      eta_minus[[pn]] <- pmax(eta_minus[[pn]], 1e-8)
    }
    if (pn == "nu" && identical(as.character(margin_dist$family[1]), "DEL")) {
      eta_plus[[pn]] <- pmin(pmax(eta_plus[[pn]], 1e-8), 1 - 1e-8)
      eta_minus[[pn]] <- pmin(pmax(eta_minus[[pn]], 1e-8), 1 - 1e-8)
    }
    denom <- eta_plus[[pn]] - eta_minus[[pn]]
    denom[!is.finite(denom) | abs(denom) < .Machine$double.eps] <- NA_real_

    Fu_plus <- calc_F_x(eta_plus, mm, margin_dist, response)
    Fu_minus <- calc_F_x(eta_minus, mm, margin_dist, response)
    Fl_plus <- calc_F_x(eta_plus, mm, margin_dist, response - 1)
    Fl_minus <- calc_F_x(eta_minus, mm, margin_dist, response - 1)

    out[[pn]] <- list(
      upper = as.numeric((Fu_plus - Fu_minus) / denom),
      lower = as.numeric((Fl_plus - Fl_minus) / denom)
    )
    out[[pn]]$upper[!is.finite(out[[pn]]$upper)] <- 0
    out[[pn]]$lower[!is.finite(out[[pn]]$lower)] <- 0
  }
  out
}
