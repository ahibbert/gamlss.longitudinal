.gl_pit <- function(object, randomize = FALSE) {
  diag_data <- .gl_fitted_distribution(object, newdata = NULL, require_response = TRUE)

  y <- diag_data$response

  params <- diag_data$params

  pit_upper <- .gl_call_family_fun("p", diag_data$family, y, params)

  pit <- pit_upper

  route <- .gl_capability_likelihood_route(object$margin_dist)

  if (randomize && identical(route, "exact_discrete_rectangle")) {
    family <- .gl_capability_margin_code(object$margin_dist)
    spec <- .gl_capability_margin_spec(family)
    if (!isTRUE(spec$randomized_pit)) {
      .gl_capability_stop(
        "gamlss_longitudinal_diagnostic_capability_error",
        paste0("Randomized PIT diagnostics are not registered for margin family '", family, "'."),
        margin_family = family
      )
    }
    pit_lower <- .gl_call_family_fun("p", diag_data$family, y - 1, params)

    pit_lower <- pmin(pmax(as.numeric(pit_lower), 0), 1)

    pit_upper <- pmin(pmax(as.numeric(pit_upper), 0), 1)

    interval_width <- pmax(pit_upper - pit_lower, 0)

    pit <- pit_lower + stats::runif(length(pit_upper)) * interval_width
  }

  pit <- pmin(pmax(as.numeric(pit), 0), 1)

  list(diag = diag_data, pit = pit)
}
