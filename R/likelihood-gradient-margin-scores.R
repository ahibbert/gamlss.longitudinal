#' Build natural-parameter margin scores for the CG analytical gradient
#'
#' @keywords internal
#' @noRd
.cg_margin_natural_scores <- function(
    margin_par,
    eta,
    eta_inv,
    mm_cg,
    calc_lik,
    margin_dist,
    copula_derivatives,
    include_dlcopdpar,
    response,
    response_margin,
    response_subject) {
  margin_score_natural <- list()
  if (length(margin_par) > 0) {
    margin_deriv_subnames <- c("m", "d", "v", "t")
    names(margin_deriv_subnames) <- c("mu", "sigma", "nu", "tau")

    for (par_name in margin_par) {
      d_name <- paste0("dld", margin_deriv_subnames[par_name])
      hit <- grep(paste0("^", d_name, "$"), names(calc_lik$margin_deriv))
      if (length(hit) == 0) {
        hit <- grep(d_name, names(calc_lik$margin_deriv))
      }
      if (length(hit) == 0) {
        margin_score_natural[[par_name]] <- rep(0, length(eta[[par_name]]))
      } else {
        margin_score_natural[[par_name]] <- as.numeric(calc_lik$margin_deriv[[hit[1]]])
      }
    }

    if (isTRUE(include_dlcopdpar)) {
      nd_impact_F <- calc_Fx_derivatives(eta_inv, mm_cg$x, margin_dist, response = response)

      order_margin <- data.frame(time = response_margin, subject = response_subject)
      margin_deriv_1 <- matrix(0, ncol = length(margin_par), nrow = length(response))
      colnames(margin_deriv_1) <- paste0("dld", margin_par)
      for (par_name in margin_par) {
        margin_deriv_1[, paste0("dld", par_name)] <- margin_score_natural[[par_name]]
      }

      margin_components_base <- cbind(
        order_margin,
        response = response,
        margin_p = calc_lik$margin_p,
        margin_d = calc_lik$margin_d,
        margin_deriv_1,
        mu = eta_inv[["mu"]]
      )
      names(margin_components_base)[seq_len(ncol(order_margin))] <- c("time", "subject")

      copula_components <- cbind(
        calc_lik$order_copula,
        row_id1 = calc_lik$copula_row_id1,
        row_id2 = calc_lik$copula_row_id2,
        dcdu1 = copula_derivatives$dcdu1,
        dcdu2 = copula_derivatives$dcdu2,
        copula_d = calc_lik$copula_d
      )

      for (par_name in margin_par) {
        margin_components <- cbind(
          margin_components_base,
          F_nd = nd_impact_F[[par_name]]
        )
        margin_components_Ft_plus <- margin_components
        margin_components_Ft_plus$time <- normalize_lag_time(margin_components_Ft_plus$time)
        margin_plus <- merge(
          margin_components,
          margin_components_Ft_plus,
          by = c("time", "subject"),
          all.x = TRUE
        )
        copula_merged <- merge(
          copula_components,
          margin_plus,
          by.x = c("time1", "subject1"),
          by.y = c("time", "subject"),
          all.x = TRUE
        )
        d1_cop <- calc_deriv_copula_wrt_margin(
          copula_merged,
          margin_par,
          par_name,
          calc_d2 = FALSE
        )[, which(margin_par == par_name)]
        n_score <- length(margin_score_natural[[par_name]])
        if (length(d1_cop) >= n_score) {
          margin_score_natural[[par_name]] <- margin_score_natural[[par_name]] + d1_cop[seq_len(n_score)]
        }
      }
    }
  }

  margin_score_natural
}
