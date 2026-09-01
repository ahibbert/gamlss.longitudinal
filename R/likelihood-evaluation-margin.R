#' Evaluate marginal likelihood pieces for the joint likelihood
#'
#' @keywords internal
#' @noRd
.gl_likelihood_margin_derivative_input <- function(
    eta_inv,
    mm,
    margin_dist,
    response) {
  n_obs <- length(response)
  margin_deriv_input <- list()
  margin_deriv_input[["y"]] <- response
  margin_deriv_input[["q"]] <- response
  margin_deriv_input[["x"]] <- response
  for (par_name in names(mm)) {
    if (par_name %in% c("mu", "sigma", "nu", "tau")) {
      margin_deriv_input[[par_name]] <- eta_inv[[par_name]]
    }
  }
  fixed_unlinked_values <- attr(margin_dist, "fixed_unlinked_values")
  if (length(fixed_unlinked_values) > 0L) {
    for (par_name in names(fixed_unlinked_values)) {
      value <- fixed_unlinked_values[[par_name]]
      if (is.numeric(value) && length(value) == 1L && is.finite(value)) {
        margin_deriv_input[[par_name]] <- rep(value, n_obs)
      }
    }
  }
  # Production callers pass a GAMLSS family object. Some internal optimizer
  # tests deliberately inject only a family name together with a derivative
  # cache; there is no fixed-family metadata to recover in that case.
  fixed_family_args <- if (is.list(margin_dist)) {
    .gl_margin_fixed_family_args(margin_dist, n_obs)
  } else {
    list()
  }
  margin_deriv_input <- c(margin_deriv_input, fixed_family_args)
  margin_deriv_input
}

#' @keywords internal
#' @noRd
.gl_likelihood_evaluate_margin_derivatives <- function(
    margin_deriv_input,
    margin_eval_cache,
    response,
    calc_margin_deriv = TRUE,
    margin_deriv_names = NULL) {
  n_obs <- length(response)
  obs_response <- !is.na(response)
  margin_deriv <- list()
  if (isTRUE(calc_margin_deriv)) {
    margin_deriv_cache <- margin_eval_cache$margin_deriv_cache
    if (!is.null(margin_deriv_names)) {
      margin_deriv_cache <- Filter(
        function(deriv_info) deriv_info$name %in% margin_deriv_names,
        margin_deriv_cache
      )
    }
    for (deriv_info in margin_deriv_cache) {
      FUN_args <- names(margin_deriv_input)[names(margin_deriv_input) %in% deriv_info$args]
      deriv_val <- tryCatch(
        do.call(deriv_info$FUN, args = margin_deriv_input[FUN_args]),
        error = function(e) rep(0, n_obs)
      )
      if (length(deriv_val) == n_obs) {
        deriv_val[!obs_response] <- 0
        deriv_val[!is.finite(deriv_val)] <- 0
      }
      margin_deriv[[deriv_info$name]] <- deriv_val
    }
  }
  margin_deriv
}

#' Evaluate marginal log densities without avoidable density-scale underflow
#'
#' @noRd
.gl_likelihood_evaluate_margin_log_density <- function(
    margin_deriv_input,
    margin_eval_cache,
    margin_dist,
    margin_d,
    discrete_margin) {
  if (!"log" %in% margin_eval_cache$margin_d_args) {
    return(log(margin_d))
  }

  log_input <- margin_deriv_input
  log_input$log <- TRUE
  margin_log_d <- .call_fast_count_family("logd", margin_dist$family[1], log_input)
  if (is.null(margin_log_d)) {
    margin_log_d <- .call_margin_family_cached(
      margin_eval_cache$margin_dFUN,
      log_input,
      names(log_input)[names(log_input) %in% margin_eval_cache$margin_d_args],
      # `log` is a scalar control argument for GAMLSS density functions. The
      # row-cache expands scalar arguments, which is not safe for this flag.
      cacheable = FALSE,
      cache_env = margin_eval_cache$family_call_cache,
      cache_prefix = paste0(margin_dist$family[1], ":logd")
    )
  }
  margin_log_d
}

#' Evaluate marginal likelihood pieces for the joint likelihood
#'
#' @keywords internal
#' @noRd
.gl_likelihood_evaluate_margins <- function(
    eta_inv,
    mm,
    margin_dist,
    margin_eval_cache,
    response,
    discrete_margin,
    calc_margin_deriv = TRUE,
    margin_deriv_names = NULL) {
  n_obs <- length(response)
  obs_response <- !is.na(response)

  margin_deriv_input <- .gl_likelihood_margin_derivative_input(
    eta_inv = eta_inv,
    mm = mm,
    margin_dist = margin_dist,
    response = response
  )
  margin_deriv <- .gl_likelihood_evaluate_margin_derivatives(
    margin_deriv_input = margin_deriv_input,
    margin_eval_cache = margin_eval_cache,
    response = response,
    calc_margin_deriv = calc_margin_deriv,
    margin_deriv_names = margin_deriv_names
  )

  FUN_args <- names(margin_deriv_input)[names(margin_deriv_input) %in% margin_eval_cache$margin_p_args]
  margin_p <- .call_fast_count_family("p", margin_dist$family[1], margin_deriv_input)
  if (is.null(margin_p)) {
    margin_p <- .call_margin_family_cached(
      margin_eval_cache$margin_pFUN,
      margin_deriv_input,
      FUN_args,
      cacheable = discrete_margin,
      cache_env = margin_eval_cache$family_call_cache,
      cache_prefix = paste0(margin_dist$family[1], ":p:upper")
    )
  }
  margin_p[!obs_response] <- NA
  margin_p[!is.finite(margin_p)] <- NA

  margin_p_lower <- NULL
  margin_log_survival <- NULL
  margin_log_survival_lower <- NULL
  likelihood_type <- "continuous_density"
  if (discrete_margin) {
    margin_deriv_input_lower <- margin_deriv_input
    lower_response <- response - 1
    margin_deriv_input_lower[["q"]] <- lower_response
    margin_deriv_input_lower[["x"]] <- lower_response
    margin_p_lower <- rep(NA_real_, length(response))
    negative_lower <- obs_response & is.finite(lower_response) & lower_response < 0
    margin_p_lower[negative_lower] <- 0
    valid_lower <- obs_response & is.finite(lower_response) & lower_response >= 0
    if (any(valid_lower)) {
      lower_call_args <- lapply(margin_deriv_input_lower, function(value) {
        if (length(value) == n_obs) value[valid_lower] else value
      })
      lower_eval <- .call_fast_count_family("p", margin_dist$family[1], lower_call_args)
      if (is.null(lower_eval)) {
        lower_eval <- .call_margin_family_cached(
          margin_eval_cache$margin_pFUN,
          lower_call_args,
          names(lower_call_args)[names(lower_call_args) %in% margin_eval_cache$margin_p_args],
          cacheable = TRUE,
          cache_env = margin_eval_cache$family_call_cache,
          cache_prefix = paste0(margin_dist$family[1], ":p:lower")
        )
      }
      margin_p_lower[valid_lower] <- lower_eval
    }
    margin_p_lower[!obs_response] <- NA
    margin_p_lower[!is.finite(margin_p_lower)] <- NA

    if (all(c("lower.tail", "log.p") %in% margin_eval_cache$margin_p_args)) {
      survival_input <- margin_deriv_input
      survival_input$lower.tail <- FALSE
      survival_input$log.p <- TRUE
      margin_log_survival <- .call_margin_family_cached(
        margin_eval_cache$margin_pFUN,
        survival_input,
        names(survival_input)[names(survival_input) %in% margin_eval_cache$margin_p_args],
        cacheable = FALSE,
        cache_env = margin_eval_cache$family_call_cache,
        cache_prefix = paste0(margin_dist$family[1], ":log-survival:upper")
      )
      lower_survival_input <- margin_deriv_input_lower
      lower_survival_input$lower.tail <- FALSE
      lower_survival_input$log.p <- TRUE
      margin_log_survival_lower <- .call_margin_family_cached(
        margin_eval_cache$margin_pFUN,
        lower_survival_input,
        names(lower_survival_input)[names(lower_survival_input) %in% margin_eval_cache$margin_p_args],
        cacheable = FALSE,
        cache_env = margin_eval_cache$family_call_cache,
        cache_prefix = paste0(margin_dist$family[1], ":log-survival:lower")
      )
      margin_log_survival[!obs_response] <- NA_real_
      margin_log_survival_lower[!obs_response] <- NA_real_
    }
    likelihood_type <- "discrete_rectangle"
  }

  FUN_args <- names(margin_deriv_input)[names(margin_deriv_input) %in% margin_eval_cache$margin_d_args]
  margin_d <- .call_fast_count_family("d", margin_dist$family[1], margin_deriv_input)
  if (is.null(margin_d)) {
    margin_d <- .call_margin_family_cached(
      margin_eval_cache$margin_dFUN,
      margin_deriv_input,
      FUN_args,
      cacheable = discrete_margin,
      cache_env = margin_eval_cache$family_call_cache,
      cache_prefix = paste0(margin_dist$family[1], ":d")
    )
  }
  margin_log_d <- .gl_likelihood_evaluate_margin_log_density(
    margin_deriv_input = margin_deriv_input,
    margin_eval_cache = margin_eval_cache,
    margin_dist = margin_dist,
    margin_d = margin_d,
    discrete_margin = discrete_margin
  )
  margin_d[!obs_response] <- NA
  margin_log_d[!obs_response] <- NA

  list(
    margin_d = margin_d,
    margin_log_d = margin_log_d,
    margin_p = margin_p,
    margin_p_lower = margin_p_lower,
    margin_log_survival = margin_log_survival,
    margin_log_survival_lower = margin_log_survival_lower,
    margin_deriv = margin_deriv,
    likelihood_type = likelihood_type
  )
}
