.copula_v2_fit_data <- function(object, data = NULL) {
  if (!inherits(object, "gamlss.longitudinal")) {
    stop("'object' must be a fitted 'gamlss.longitudinal' object.")
  }

  copula_spec <- get_copula_dist(object$copula_dist)

  copula_family_name <- .copula_family_code(copula_spec$copula_dist)

  if (is.null(data)) {
    mm_use <- object$model_matrix

    response <- object$response

    subject <- object$response_subject

    time <- object$response_margin
  } else {
    nd <- .gl_prepare_newdata_internal(object, data, require_response = TRUE)

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
        dataset = nd,
        quiet_gamlss2 = TRUE,
        preserve_factor_levels = TRUE
      )
    )

    mm_use <- .gl_align_model_matrix_columns(mm_use, object$model_matrix)

    response <- nd$response

    subject <- nd$subject

    time <- nd$time
  }

  eta_out <- calc_eta(
    par_cov = object$par,
    mm = mm_use,
    margin_dist = object$margin_dist,
    copula_link = copula_spec$copula_link,
    par_s = object$par_s
  )

  # Extract only margin parameters that are actually in eta_out$eta_inv

  margin_param_names <- names(object$margin_dist$parameters)

  margin_params <- list()

  for (param_name in margin_param_names) {
    if (param_name %in% names(eta_out$eta_inv)) {
      margin_params[[param_name]] <- eta_out$eta_inv[[param_name]]
    }
  }

  theta_fit <- if ("theta" %in% names(eta_out$eta_inv)) eta_out$eta_inv$theta else numeric(0)

  zeta_fit <- if ("zeta" %in% names(eta_out$eta_inv)) eta_out$eta_inv$zeta else numeric(0)

  # Align response-side vectors to a common leading length.

  margin_min_n <- if (length(margin_params) > 0) {
    min(vapply(margin_params, length, integer(1)))
  } else {
    length(response)
  }

  common_n <- min(length(response), length(subject), length(time), margin_min_n)

  if (!is.finite(common_n) || common_n < 1) {
    stop("No finite fitted observations are available for copula diagnostics.")
  }

  response <- response[seq_len(common_n)]

  subject <- subject[seq_len(common_n)]

  time <- time[seq_len(common_n)]

  margin_params <- lapply(margin_params, function(x) x[seq_len(common_n)])

  align_copula_param <- function(param_vec) {
    n_resp <- common_n

    if (length(param_vec) == 0) {
      return(rep(NA_real_, n_resp))
    }

    # Full-row parameterization.

    if (length(param_vec) == n_resp) {
      return(param_vec)
    }

    # Pair-row parameterization: parameters correspond to times 1:(T-1) only.

    margin_names <- sort(unique(time))

    left_time_rows <- which(time %in% margin_names[seq_len(max(1, length(margin_names) - 1))])

    if (length(param_vec) == length(left_time_rows)) {
      out <- rep(NA_real_, n_resp)

      out[left_time_rows] <- param_vec

      return(out)
    }

    # Fallback for unexpected lengths.

    rep(param_vec, length.out = n_resp)
  }

  theta_fit <- align_copula_param(theta_fit)

  zeta_fit <- align_copula_param(zeta_fit)

  # Filter by finite values

  keep <- is.finite(response)

  for (param_name in names(margin_params)) {
    keep <- keep & is.finite(margin_params[[param_name]])
  }

  response <- response[keep]

  subject <- subject[keep]

  time <- time[keep]

  margin_params <- lapply(margin_params, function(x) x[keep])

  theta_fit <- theta_fit[keep]

  zeta_fit <- zeta_fit[keep]

  if (length(response) == 0) {
    stop("No finite fitted observations are available for copula diagnostics.")
  }

  # Convert margin_dist$family to family name if needed

  family_name <- object$margin_dist$family[1]

  if (!is.character(family_name)) {
    family_name <- object$margin_dist$family[1]$family
  }

  u <- .gl_call_family_fun("p", family_name, response, margin_params)

  u <- .copula_v2_clamp01(u)

  family_num <- tryCatch(
    {
      .copula_family_code(copula_family_name)
    },
    error = function(e) NA_character_
  )

  # Compute tau_fit, suppressing coercion warnings

  tau_fit <- suppressWarnings(
    vapply(seq_along(theta_fit), function(i) {
      .copula_v2_tau_from_par(family_num, theta_fit[i], zeta_fit[i])
    }, numeric(1), USE.NAMES = FALSE)
  )

  data.frame(
    subject = subject,
    time = time,
    response = response,
    u = u,
    theta_fit = theta_fit,
    zeta_fit = zeta_fit,
    tau_fit = tau_fit,
    stringsAsFactors = FALSE
  )
}
