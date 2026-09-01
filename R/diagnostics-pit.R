.gl_with_preserved_rng <- function(seed = NULL, code) {
  had_seed <- exists(".Random.seed", envir = .GlobalEnv, inherits = FALSE)
  old_seed <- if (had_seed) get(".Random.seed", envir = .GlobalEnv, inherits = FALSE) else NULL
  on.exit({
    if (had_seed) {
      assign(".Random.seed", old_seed, envir = .GlobalEnv)
    } else if (exists(".Random.seed", envir = .GlobalEnv, inherits = FALSE)) {
      rm(".Random.seed", envir = .GlobalEnv)
    }
  }, add = TRUE)
  if (!is.null(seed)) {
    if (!is.numeric(seed) || length(seed) != 1L || !is.finite(seed)) {
      stop("'seed' must be NULL or one finite number.", call. = FALSE)
    }
    set.seed(as.integer(seed))
  }
  force(code)
}

.gl_pit <- function(object, randomize = NULL, seed = 1L) {
  diag_data <- .gl_fitted_distribution(object, newdata = NULL, require_response = TRUE)

  y <- diag_data$response

  params <- diag_data$params

  pit_upper <- .gl_call_family_fun("p", diag_data$family, y, params)

  pit <- pit_upper

  route <- .gl_capability_likelihood_route(object$margin_dist)
  if (is.null(randomize)) {
    randomize <- identical(route, "exact_discrete_rectangle")
  } else {
    randomize <- isTRUE(randomize)
  }

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

    pit <- .gl_with_preserved_rng(
      seed,
      pit_lower + stats::runif(length(pit_upper)) * interval_width
    )
  }

  pit <- pmin(pmax(as.numeric(pit), 0), 1)

  list(
    diag = diag_data,
    pit = pit,
    randomized = isTRUE(randomize) && identical(route, "exact_discrete_rectangle"),
    seed = if (isTRUE(randomize) && identical(route, "exact_discrete_rectangle") && !is.null(seed)) as.integer(seed) else NA_integer_
  )
}
