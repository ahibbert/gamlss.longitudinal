.select_margin_time_intercept_fits <- function(response, time, families, trace = FALSE) {
  response <- as.numeric(response)

  time <- as.character(time)

  keep <- is.finite(response) & !is.na(time)

  if (sum(keep) < 3L) {
    stop("Need at least three finite response values with non-missing time values to screen time-intercept margins.", call. = FALSE)
  }

  fit_data <- data.frame(
    response = response[keep],
    time_intercept = factor(time[keep], levels = unique(time[keep])),
    stringsAsFactors = FALSE
  )

  has_time_contrast <- length(unique(fit_data$time_intercept)) > 1L

  mu_formula <- if (has_time_contrast) {
    stats::as.formula("response ~ time_intercept")
  } else {
    stats::as.formula("response ~ 1")
  }

  par_formula <- if (has_time_contrast) {
    stats::as.formula("~ time_intercept")
  } else {
    stats::as.formula("~ 1")
  }

  rows <- lapply(families, function(family_name) {
    family <- .margin_family_object(family_name)

    if (is.null(family) || is.null(family$parameters)) {
      return(data.frame(
        family = family_name,
        time_intercept_AIC = NA_real_,
        converged = FALSE,
        stringsAsFactors = FALSE
      ))
    }

    fit_args <- list(
      formula = mu_formula,
      family = family,
      data = fit_data,
      trace = trace,
      control = gamlss::gamlss.control(n.cyc = 100, trace = trace)
    )

    for (par_name in setdiff(names(family$parameters), "mu")) {
      fit_args[[paste0(par_name, ".formula")]] <- par_formula
    }

    fit <- tryCatch(
      {
        if (isTRUE(trace)) {
          do.call(gamlss::gamlss, fit_args)
        } else {
          fit <- NULL

          invisible(utils::capture.output({
            fit <- suppressWarnings(suppressMessages(do.call(gamlss::gamlss, fit_args)))
          }))

          fit
        }
      },
      error = function(e) NULL
    )

    if (is.null(fit)) {
      return(data.frame(
        family = family_name,
        time_intercept_AIC = NA_real_,
        converged = FALSE,
        stringsAsFactors = FALSE
      ))
    }

    aic <- tryCatch(stats::AIC(fit), error = function(e) NA_real_)

    data.frame(
      family = family_name,
      time_intercept_AIC = as.numeric(aic)[1L],
      converged = !identical(fit$converged, FALSE),
      stringsAsFactors = FALSE
    )
  })

  do.call(rbind, rows)
}
