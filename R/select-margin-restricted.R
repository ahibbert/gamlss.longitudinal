.select_margin_restricted_fits <- function(
    response,
    families,
    try.gamlss = FALSE,
    trace = FALSE,
    ...) {
  if (!is.character(families) || !length(families) || anyNA(families) ||
      any(!nzchar(families))) {
    stop("'families' must contain one or more family names.", call. = FALSE)
  }
  families <- unique(families)
  dots <- list(...)
  k <- if (!is.null(dots$k)) dots$k else 2
  dots$k <- NULL
  fit_data <- data.frame(.response = response)

  aic <- vapply(families, function(family) {
    fit_ml <- function() {
      do.call(
        gamlss::gamlssML,
        c(
          list(formula = .response ~ 1, family = family, data = fit_data),
          dots
        )
      )
    }
    fit <- if (isTRUE(trace)) {
      try(fit_ml(), silent = TRUE)
    } else {
      suppressWarnings(suppressMessages(try(fit_ml(), silent = TRUE)))
    }

    if (inherits(fit, "try-error") && isTRUE(try.gamlss)) {
      fit_gamlss <- function() {
        do.call(
          gamlss::gamlss,
          c(
            list(
              formula = .response ~ 1,
              family = family,
              data = fit_data,
              trace = trace
            ),
            dots
          )
        )
      }
      fit <- if (isTRUE(trace)) {
        try(fit_gamlss(), silent = TRUE)
      } else {
        suppressWarnings(suppressMessages(try(fit_gamlss(), silent = TRUE)))
      }
    }
    if (inherits(fit, "try-error")) return(NA_real_)
    as.numeric(stats::AIC(fit, k = k))[1L]
  }, numeric(1))

  aic[is.finite(aic)]
}
