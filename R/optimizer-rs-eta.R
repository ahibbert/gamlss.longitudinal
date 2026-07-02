#' Build the RS eta calculator for cached or full matrix paths
#'
#' @noRd
.gl_build_rs_eta_calculator <- function(
    rs_design_cache,
    mm,
    margin_dist,
    copula_link,
    option_fn = getOption,
    cached_eta_fn = .calc_eta_rs_cached,
    full_eta_fn = calc_eta) {
  force(rs_design_cache)
  force(mm)
  force(margin_dist)
  force(copula_link)
  force(option_fn)
  force(cached_eta_fn)
  force(full_eta_fn)

  function(par_cov_current, par_s_current, update_only = NULL, eta_out_current = NULL) {
    if (isTRUE(option_fn("gamlss.longitudinal.fast_rs_eta", TRUE))) {
      cached_eta_fn(
        rs_design_cache = rs_design_cache,
        par_cov = par_cov_current,
        par_s = par_s_current,
        margin_dist = margin_dist,
        copula_link = copula_link,
        update_only = update_only,
        eta_out = eta_out_current
      )
    } else {
      full_eta_fn(
        par_cov_current,
        mm,
        margin_dist,
        copula_link,
        par_s = par_s_current
      )
    }
  }
}

#' Validate RS eta lengths against the response
#'
#' @noRd
.gl_validate_rs_eta_lengths <- function(eta_inv, mm, response, margin_params = c("mu", "sigma", "nu", "tau")) {
  n_resp <- length(response)
  active_margin_params <- intersect(names(mm$x), margin_params)
  bad_lengths <- active_margin_params[
    sapply(active_margin_params, function(pn) length(eta_inv[[pn]]) != n_resp)
  ]

  if (length(bad_lengths) > 0) {
    detail <- paste(sapply(bad_lengths, function(pn) {
      paste0(pn, "=", length(eta_inv[[pn]]), " vs response=", n_resp)
    }), collapse = ", ")

    stop(
      "ERROR: Parameter vector lengths do not match response length. ",
      detail,
      ".\nThis usually indicates model-matrix rows were dropped (often due to NA handling)."
    )
  }

  invisible(TRUE)
}
