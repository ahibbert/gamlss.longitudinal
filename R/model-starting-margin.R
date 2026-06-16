.starting_margin_parameter_values <- function(margin_dist, finite_response, dataset) {
  margin_par_already_eta <- FALSE

  # Method of moments method for starting values where it's simple, otherwise a brief GAMLSS starting fit. 
  if (margin_dist$family[1] == "GA" | margin_dist$family[1] == "EXP") {
    margin_par <- c(
      mean(finite_response),
      stats::sd(finite_response) / mean(finite_response),
      .starting_moment_skewness(finite_response),
      .starting_moment_kurtosis(finite_response)
    )
  } else if (margin_dist$family[1] == "NO") {
    margin_par <- c(
      mean(finite_response),
      stats::sd(finite_response)
    )
  } else if (margin_dist$family[1] == "PO") {
    margin_par <- c(
      mean(finite_response)
    )
  } else if (margin_dist$family[1] == "NBI") {
    margin_par <- c(
      mean(finite_response),
      stats::sd(finite_response) / mean(finite_response)
    )
  } else {
    cat("Fitting initial GAMLSS model for margin to obtain starting values...\n")

    # Deliberately low-iteration startup fit; silence expected convergence warnings.
    start_fit <- tryCatch(
      {
        fit <- NULL
        invisible(utils::capture.output({
          fit <- suppressWarnings(suppressMessages(
            gamlss(

              dataset$response ~ 1,
              family = margin_dist,
              control = gamlss::gamlss.control(n.cyc = 5, trace = FALSE)
            )
          ))
        }))

        fit
      },
      error = function(e) NULL
    )

    margin_par <- vapply(names(margin_dist$parameters), function(parameter) {
      cf <- tryCatch(stats::coef(start_fit, what = parameter), error = function(e) numeric(0))

      if (length(cf) > 0L && is.finite(cf[[1L]])) {
        return(as.numeric(cf[[1L]]))
      }

      qfun <- get(paste0("q", margin_dist$family[1]), envir = asNamespace("gamlss.dist"), inherits = FALSE)

      value <- tryCatch(eval(formals(qfun)[[parameter]], envir = baseenv()), error = function(e) NA_real_)

      if (!is.numeric(value) || length(value) != 1L || !is.finite(value)) {
        value <- switch(parameter,
          mu = if (.is_discrete_margin(margin_dist)) 3 else mean(finite_response),
          sigma = stats::sd(finite_response) / max(abs(mean(finite_response)), 1e-8),
          nu = 0.5,
          tau = 2,
          0
        )
      }

      linkfun <- margin_dist[[paste0(parameter, ".linkfun")]]

      if (is.null(linkfun)) {
        return(as.numeric(value))
      }

      as.numeric(linkfun(value))[1L]
    }, numeric(1))

    names(margin_par) <- names(margin_dist$parameters)

    margin_par_already_eta <- TRUE
  }

  names(margin_par) <- names(margin_dist$parameters)

  margin_par <- margin_par[!is.na(names(margin_par))]

  list(
    margin_par = margin_par,
    margin_par_already_eta = margin_par_already_eta
  )
}
