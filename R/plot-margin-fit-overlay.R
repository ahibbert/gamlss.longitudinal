#' @keywords internal
#' @noRd
.plot_margin_constant_params <- function(y, family, warn = TRUE, fit_control = gamlss::gamlss.control(n.cyc = 50)) {
  y <- as.numeric(y)
  y <- y[is.finite(y)]
  if (length(y) < 3L) {
    stop("Need at least three finite response values to fit a marginal overlay.", call. = FALSE)
  }

  fit <- NULL
  invisible(utils::capture.output({
    fit <- suppressWarnings(suppressMessages(
      gamlss::gamlss(y ~ 1, family = family, trace = FALSE, control = fit_control)
    ))
  }))
  if (!.plot_margin_check_fit(fit, family, warn = warn)) {
    return(NULL)
  }
  params <- lapply(names(family$parameters), function(par_name) {
    as.numeric(stats::fitted(fit, what = par_name))
  })
  names(params) <- names(family$parameters)
  params
}

#' @keywords internal
#' @noRd
.plot_margin_time_intercept_params <- function(y, time, family, warn = TRUE, fit_control = gamlss::gamlss.control(n.cyc = 50)) {
  y <- as.numeric(y)
  time <- as.character(time)
  keep <- is.finite(y) & !is.na(time)
  if (sum(keep) < 3L) {
    stop("Need at least three finite response values with non-missing time values to fit a time-intercept marginal overlay.", call. = FALSE)
  }

  fit_data <- data.frame(
    y = y[keep],
    time_intercept = factor(time[keep], levels = unique(time[keep])),
    stringsAsFactors = FALSE
  )
  formula_time <- stats::as.formula("y ~ time_intercept")
  fit_args <- list(formula = formula_time, family = family, data = fit_data, trace = FALSE, control = fit_control)
  for (par_name in setdiff(names(family$parameters), "mu")) {
    fit_args[[paste0(par_name, ".formula")]] <- stats::as.formula("~ time_intercept")
  }
  fit <- NULL
  invisible(utils::capture.output({
    fit <- suppressWarnings(suppressMessages(do.call(gamlss::gamlss, fit_args)))
  }))
  if (!.plot_margin_check_fit(fit, family, warn = warn)) {
    return(NULL)
  }

  params <- lapply(names(family$parameters), function(par_name) {
    out <- rep(NA_real_, length(y))
    out[keep] <- as.numeric(stats::fitted(fit, what = par_name))
    out
  })
  names(params) <- names(family$parameters)
  params
}

#' @keywords internal
#' @noRd
.plot_margin_check_fit <- function(fit, family, warn = TRUE) {
  if (is.null(fit)) {
    stop("The marginal overlay fit failed.", call. = FALSE)
  }
  if (identical(fit$converged, FALSE)) {
    if (isTRUE(warn)) {
      family_name <- .plot_margin_family_name(family)
      warning(
        "The ", family_name, " marginal overlay fit did not converge; no fitted density was drawn.",
        call. = FALSE
      )
    }
    return(FALSE)
  }
  TRUE
}
