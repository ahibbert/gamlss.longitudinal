.select_copula_u_from_fit <- function(object) {
  copula_link <- get_copula_dist(object$copula_dist)$copula_link
  eta_out <- calc_eta(
    object$par,
    object$model_matrix,
    object$margin_dist,
    copula_link,
    object$par_s
  )
  u <- calc_F_x(
    eta_out$eta_inv,
    object$model_matrix$x,
    object$margin_dist,
    object$response
  )
  .copula_clamp01(u)
}

.select_copula_pit_from_fit <- function(object, pit_seed = 1L) {
  diag_data <- .gl_fitted_distribution(object, newdata = NULL, require_response = TRUE)
  upper <- .gl_call_family_fun(
    "p",
    diag_data$family,
    diag_data$response,
    diag_data$params
  )
  upper <- pmin(pmax(as.numeric(upper), 0), 1)
  route <- .gl_capability_likelihood_route(object$margin_dist)
  if (!identical(route, "exact_discrete_rectangle")) {
    return(list(
      pit = upper,
      upper = NULL,
      lower = NULL,
      log_mass = NULL,
      log_survival = NULL,
      subject = diag_data$subject,
      time = diag_data$time,
      likelihood = "copula_density",
      randomized = FALSE,
      pit_seed = NA_integer_
    ))
  }

  lower <- .gl_call_family_fun(
    "p",
    diag_data$family,
    diag_data$response - 1,
    diag_data$params
  )
  lower <- pmin(pmax(as.numeric(lower), 0), 1)
  log_mass <- .gl_call_family_fun(
    "d",
    diag_data$family,
    diag_data$response,
    diag_data$params,
    extra_args = list(log = TRUE)
  )
  log_survival <- tryCatch(
    .gl_call_family_fun(
      "p",
      diag_data$family,
      diag_data$response,
      diag_data$params,
      extra_args = list(lower.tail = FALSE, log.p = TRUE)
    ),
    error = function(e) rep(NA_real_, length(upper))
  )
  pit <- .gl_with_preserved_rng(
    pit_seed,
    lower + stats::runif(length(upper)) * pmax(upper - lower, 0)
  )
  list(
    pit = .copula_clamp01(pit),
    upper = upper,
    lower = lower,
    log_mass = as.numeric(log_mass),
    log_survival = as.numeric(log_survival),
    subject = diag_data$subject,
    time = diag_data$time,
    likelihood = "exact_discrete_rectangle",
    randomized = TRUE,
    pit_seed = as.integer(pit_seed)
  )
}
