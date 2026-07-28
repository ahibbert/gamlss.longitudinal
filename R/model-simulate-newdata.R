.gl_simulation_time_levels <- function(time, reference_time = NULL) {
  if (is.factor(time) || is.factor(reference_time)) {
    lev <- unique(c(levels(reference_time), levels(time)))

    return(lev[!is.na(lev)])
  }

  u <- unique(c(reference_time, time))

  if (is.numeric(u) || is.integer(u)) sort(u) else sort(as.character(u))
}

.gl_simulation_newdata <- function(object, newdata) {
  copula_spec <- get_copula_dist(object$copula_dist)

  nd <- .gl_prepare_newdata_internal(object, newdata, require_response = FALSE)

  nd$.gl_sim_row_id <- seq_len(nrow(nd))

  time_grid <- .gl_simulation_time_grid(nd, object)

  time_for_grid <- time_grid$time_for_grid

  time_levels <- time_grid$time_levels

  nd$.gl_sim_time_idx <- time_grid$time_idx

  nd_eval <- nd[order(nd$subject, nd$.gl_sim_time_idx, nd$.gl_sim_row_id), , drop = FALSE]

  mm_use <- do.call(
    create_model_matrices,
    list(
      mu.formula = object$formulas_int$mu,
      sigma.formula = object$formulas_int$sigma,
      nu.formula = object$formulas_int$nu,
      tau.formula = object$formulas_int$tau,
      theta.formula = object$formulas_int$theta,
      zeta.formula = object$formulas_int$zeta,
      margin.family = object$margin_dist,
      copula.family = object$copula_dist,
      copula.link = copula_spec$copula_link,
      dataset = nd_eval,
      quiet_gamlss2 = TRUE,
      preserve_factor_levels = TRUE
    )
  )

  mm_use <- .gl_align_model_matrix_columns(mm_use, object$model_matrix)

  eta_out <- calc_eta(
    par_cov = object$par,
    mm = mm_use,
    margin_dist = object$margin_dist,
    copula_link = copula_spec$copula_link,
    par_s = object$par_s
  )

  margin_params <- eta_out$eta_inv[names(object$margin_dist$parameters)]

  common_n <- min(
    nrow(nd_eval),
    if (length(margin_params) > 0) min(vapply(margin_params, length, integer(1))) else nrow(nd_eval)
  )

  if (!is.finite(common_n) || common_n < 1L) {
    stop("No newdata rows are available for simulation.", call. = FALSE)
  }

  nd_eval <- nd_eval[seq_len(common_n), , drop = FALSE]

  margin_params <- lapply(margin_params, function(x) x[seq_len(common_n)])

  keep <- rep(TRUE, common_n)

  for (par_name in names(margin_params)) {
    keep <- keep & is.finite(margin_params[[par_name]])
  }

  if (!all(keep)) {
    stop("newdata produced non-finite fitted marginal parameters.", call. = FALSE)
  }

  theta_fit <- .gl_simulation_align_dependence_parameter(
    eta_out$eta_inv$theta %||% numeric(0),
    common_n,
    nd_eval,
    time_levels
  )

  zeta_fit <- .gl_simulation_align_dependence_parameter(
    eta_out$eta_inv$zeta %||% numeric(0),
    common_n,
    nd_eval,
    time_levels
  )

  ord <- order(nd_eval$.gl_sim_row_id)

  margin_params <- lapply(margin_params, function(x) as.numeric(x[ord]))

  list(
    diag_data = list(
      response = nd_eval$response[ord],
      params = margin_params,
      family = object$margin_dist$family[1],
      subject = nd_eval$subject[ord],
      time = nd_eval$time[ord],
      keep_index = nd_eval$.gl_sim_row_id[ord]
    ),
    time_levels = time_levels,
    fit_data = data.frame(
      subject = nd_eval$subject[ord],
      time = time_for_grid[nd_eval$.gl_sim_row_id[ord]],
      theta_fit = theta_fit[ord],
      zeta_fit = zeta_fit[ord],
      stringsAsFactors = FALSE
    )
  )
}
