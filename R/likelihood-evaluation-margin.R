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
.gl_likelihood_subset_row_input <- function(input, keep, n_obs) {
  lapply(input, function(value) {
    if (length(value) == n_obs) value[keep] else value
  })
}

#' Derive discrete CDF companions from one upper CDF and the point mass
#'
#' Uses F(y - 1) = F(y) - P(Y = y) and
#' S(y - 1) = S(y) + P(Y = y). A directly evaluated upper-tail probability can
#' be supplied for observations where subtracting the CDF from one is unstable.
#'
#' @noRd
.gl_discrete_cdf_companions <- function(
    upper,
    mass,
    log_mass = log(mass),
    direct_log_survival = NULL,
    tolerance = 1e-12) {
  upper <- as.numeric(upper)
  mass <- as.numeric(mass)
  log_mass <- as.numeric(log_mass)
  valid <- is.finite(upper) & is.finite(mass) &
    upper >= -tolerance & upper <= 1 + tolerance & mass >= -tolerance
  lower <- rep(NA_real_, length(upper))
  lower[valid] <- upper[valid] - mass[valid]
  lower[valid & lower > -tolerance & lower < 0] <- 0
  lower[valid & lower > 1 & lower < 1 + tolerance] <- 1
  valid_lower <- valid & lower >= 0 & lower <= 1

  log_survival <- rep(NA_real_, length(upper))
  log_survival_lower <- rep(NA_real_, length(upper))
  valid_upper <- valid & upper >= 0 & upper <= 1
  log_survival[valid_upper] <- log1p(-upper[valid_upper])
  if (!is.null(direct_log_survival)) {
    direct_log_survival <- as.numeric(direct_log_survival)
    use_direct <- valid & !is.na(direct_log_survival) & direct_log_survival <= 0
    log_survival[use_direct] <- direct_log_survival[use_direct]
  }

  valid_survival_lower <- valid_lower & !is.na(log_survival) &
    !is.na(log_mass) & log_mass <= 0
  if (any(valid_survival_lower)) {
    maximum <- pmax(
      log_survival[valid_survival_lower],
      log_mass[valid_survival_lower]
    )
    both_zero <- is.infinite(maximum) & maximum < 0
    combined <- rep(-Inf, length(maximum))
    if (any(!both_zero)) {
      combined[!both_zero] <- maximum[!both_zero] + log(
        exp(log_survival[valid_survival_lower][!both_zero] - maximum[!both_zero]) +
          exp(log_mass[valid_survival_lower][!both_zero] - maximum[!both_zero])
      )
    }
    log_survival_lower[valid_survival_lower] <- combined
  }
  list(
    lower = lower,
    log_survival = log_survival,
    log_survival_lower = log_survival_lower
  )
}

.gl_numerical_margin_log_density_derivative <- function(
    par_name,
    margin_deriv_input,
    margin_eval_cache,
    margin_dist,
    response,
    h = 1e-5) {
  n_obs <- length(response)
  observed <- is.finite(response)
  out <- numeric(n_obs)
  if (!any(observed) || !par_name %in% names(margin_deriv_input)) return(out)

  input <- .gl_likelihood_subset_row_input(
    margin_deriv_input, observed, n_obs
  )
  value <- input[[par_name]]
  linkfun <- margin_dist[[paste0(par_name, ".linkfun")]]
  linkinv <- margin_dist[[paste0(par_name, ".linkinv")]]
  if (!is.function(linkfun) || !is.function(linkinv)) {
    out[observed] <- NA_real_
    return(out)
  }
  eta <- tryCatch(as.numeric(linkfun(value)), error = function(e) rep(NA_real_, length(value)))
  step <- h * pmax(1, abs(eta))
  plus_value <- tryCatch(as.numeric(linkinv(eta + step)), error = function(e) rep(NA_real_, length(value)))
  minus_value <- tryCatch(as.numeric(linkinv(eta - step)), error = function(e) rep(NA_real_, length(value)))

  evaluate_log_density <- function(call_input) {
    if ("log" %in% margin_eval_cache$margin_d_args) {
      call_input$log <- TRUE
      fast <- .call_fast_count_family("logd", margin_dist$family[1], call_input)
      if (!is.null(fast)) return(as.numeric(fast))
    }
    args <- names(call_input)[names(call_input) %in% margin_eval_cache$margin_d_args]
    density <- .call_margin_family_cached(
      margin_eval_cache$margin_dFUN,
      call_input,
      args,
      cacheable = FALSE,
      cache_env = margin_eval_cache$family_call_cache,
      cache_prefix = paste0(margin_dist$family[1], ":numerical-logd"),
      force_scalar = margin_eval_cache$force_scalar_family
    )
    if ("log" %in% margin_eval_cache$margin_d_args) density else log(density)
  }

  plus_input <- input
  minus_input <- input
  plus_input[[par_name]] <- plus_value
  minus_input[[par_name]] <- minus_value
  log_plus <- evaluate_log_density(plus_input)
  log_minus <- evaluate_log_density(minus_input)
  derivative <- (log_plus - log_minus) / (plus_value - minus_value)
  derivative[!is.finite(derivative)] <- NA_real_
  out[observed] <- derivative
  out
}

#' @keywords internal
#' @noRd
.gl_likelihood_evaluate_margin_derivatives <- function(
    margin_deriv_input,
    margin_eval_cache,
    margin_dist,
    mm,
    response,
    calc_margin_deriv = TRUE,
    margin_deriv_names = NULL) {
  n_obs <- length(response)
  obs_response <- is.finite(response)
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
      observed_input <- .gl_likelihood_subset_row_input(
        margin_deriv_input, obs_response, n_obs
      )
      FUN_args <- names(observed_input)[names(observed_input) %in% deriv_info$args]
      observed_val <- .call_margin_family_cached(
        deriv_info$FUN,
        observed_input,
        FUN_args,
        cacheable = FALSE,
        cache_env = margin_eval_cache$family_call_cache,
        cache_prefix = paste0(margin_dist$family[1], ":", deriv_info$name),
        force_scalar = margin_eval_cache$force_scalar_family
      )
      if (length(observed_val) == 1L) {
        observed_val <- rep(observed_val, sum(obs_response))
      }
      deriv_val <- rep(NA_real_, n_obs)
      if (length(observed_val) == sum(obs_response)) {
        deriv_val[obs_response] <- observed_val
      }
      deriv_val[!obs_response] <- 0
      margin_deriv[[deriv_info$name]] <- deriv_val
    }

    suffix <- c(mu = "m", sigma = "d", nu = "v", tau = "t")
    parameters <- intersect(names(mm), names(suffix))
    expected_names <- paste0("dld", suffix[parameters])
    if (!is.null(margin_deriv_names)) {
      keep <- expected_names %in% margin_deriv_names
      parameters <- parameters[keep]
      expected_names <- expected_names[keep]
    }
    for (ii in seq_along(parameters)) {
      name <- expected_names[[ii]]
      current <- margin_deriv[[name]]
      if (is.null(current) || length(current) != n_obs) {
        current <- rep(NA_real_, n_obs)
        current[!obs_response] <- 0
      }
      failed <- obs_response & !is.finite(current)
      if (any(failed)) {
        numerical <- .gl_numerical_margin_log_density_derivative(
          par_name = parameters[[ii]],
          margin_deriv_input = margin_deriv_input,
          margin_eval_cache = margin_eval_cache,
          margin_dist = margin_dist,
          response = response
        )
        current[failed] <- numerical[failed]
      }
      current[!is.finite(current)] <- 0
      margin_deriv[[name]] <- current
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
      cache_prefix = paste0(margin_dist$family[1], ":logd"),
      force_scalar = margin_eval_cache$force_scalar_family
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
  obs_response <- is.finite(response)

  margin_deriv_input <- .gl_likelihood_margin_derivative_input(
    eta_inv = eta_inv,
    mm = mm,
    margin_dist = margin_dist,
    response = response
  )
  margin_deriv <- .gl_likelihood_evaluate_margin_derivatives(
    margin_deriv_input = margin_deriv_input,
    margin_eval_cache = margin_eval_cache,
    margin_dist = margin_dist,
    mm = mm,
    response = response,
    calc_margin_deriv = calc_margin_deriv,
    margin_deriv_names = margin_deriv_names
  )

  observed_input <- .gl_likelihood_subset_row_input(margin_deriv_input, obs_response, n_obs)
  margin_p <- .gl_evaluate_margin_cdf(
    input = margin_deriv_input,
    margin_dist = margin_dist,
    response = response,
    FUN = margin_eval_cache$margin_pFUN,
    FUN_args = margin_eval_cache$margin_p_args,
    cacheable = discrete_margin,
    cache_env = margin_eval_cache$family_call_cache,
    cache_prefix = paste0(margin_dist$family[1], ":p:upper"),
    use_fast = discrete_margin,
    strict = FALSE
  )

  margin_d <- rep(NA_real_, n_obs)
  margin_log_d <- rep(NA_real_, n_obs)
  if (any(obs_response)) {
    FUN_args <- names(observed_input)[names(observed_input) %in% margin_eval_cache$margin_d_args]
    observed_d <- .call_fast_count_family("d", margin_dist$family[1], observed_input)
    if (is.null(observed_d)) {
      observed_d <- .call_margin_family_cached(
        margin_eval_cache$margin_dFUN,
        observed_input,
        FUN_args,
        cacheable = discrete_margin,
        cache_env = margin_eval_cache$family_call_cache,
        cache_prefix = paste0(margin_dist$family[1], ":d"),
        force_scalar = margin_eval_cache$force_scalar_family
      )
    }
    observed_log_d <- .gl_likelihood_evaluate_margin_log_density(
      margin_deriv_input = observed_input,
      margin_eval_cache = margin_eval_cache,
      margin_dist = margin_dist,
      margin_d = observed_d,
      discrete_margin = discrete_margin
    )
    margin_d[obs_response] <- observed_d
    margin_log_d[obs_response] <- observed_log_d
  }

  margin_p_lower <- NULL
  margin_log_survival <- NULL
  margin_log_survival_lower <- NULL
  likelihood_type <- "continuous_density"
  if (discrete_margin) {
    direct_log_survival <- rep(NA_real_, n_obs)
    stable_tail <- obs_response & is.finite(margin_p) & margin_p > 0.9
    if (any(stable_tail) &&
        all(c("lower.tail", "log.p") %in% margin_eval_cache$margin_p_args)) {
      survival_input <- .gl_likelihood_subset_row_input(
        margin_deriv_input, stable_tail, n_obs
      )
      survival_input$lower.tail <- FALSE
      survival_input$log.p <- TRUE
      survival_args <- names(survival_input)[
        names(survival_input) %in% margin_eval_cache$margin_p_args
      ]
      evaluated_survival <- .call_fast_count_family(
        "p", margin_dist$family[1], survival_input
      )
      if (is.null(evaluated_survival)) {
        evaluated_survival <- .call_margin_family_cached(
          margin_eval_cache$margin_pFUN,
          survival_input,
          survival_args,
          cacheable = FALSE,
          cache_env = margin_eval_cache$family_call_cache,
          cache_prefix = paste0(margin_dist$family[1], ":log-survival:upper"),
          force_scalar = margin_eval_cache$force_scalar_family
        )
      }
      direct_log_survival[stable_tail] <- evaluated_survival
    }
    companions <- .gl_discrete_cdf_companions(
      margin_p,
      margin_d,
      log_mass = margin_log_d,
      direct_log_survival = direct_log_survival
    )
    margin_p_lower <- companions$lower
    margin_log_survival <- companions$log_survival
    margin_log_survival_lower <- companions$log_survival_lower
    likelihood_type <- "discrete_rectangle"
  }

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
