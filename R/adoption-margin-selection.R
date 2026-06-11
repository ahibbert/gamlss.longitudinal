#' Extract the best-fitting candidate from a selection result

#'

#' @param x A selection result, such as from [select_margin()] or

#'   [select_copula()].

#' @param ... Reserved for methods.

#'

#' @return `best_fit()` returns a list of selected-fit metadata.

#'   `best_fit_family()` returns the family value that can be supplied directly

#'   to fitting helpers.

#' @export

best_fit <- function(x, ...) {

  UseMethod("best_fit")

}


#' @rdname best_fit

#' @export

best_fit_family <- function(x, ...) {

  UseMethod("best_fit_family")

}


.margin_family_object <- function(family_name) {

  if (length(family_name) != 1L || is.na(family_name) || !nzchar(family_name)) {

    return(NULL)

  }

  family_fun <- tryCatch(

    get(family_name, envir = asNamespace("gamlss.dist"), mode = "function", inherits = FALSE),

    error = function(e) NULL

  )

  if (is.null(family_fun)) {

    return(NULL)

  }

  tryCatch(do.call(family_fun, list()), error = function(e) NULL)

}


#' @export

best_fit.margin_selection <- function(x, ...) {

  if (nrow(x) == 0L) {

    return(list(

      family_name = NA_character_,

      family = NULL,

      rank = NA_integer_,

      AIC = NA_real_,

      type = attr(x, "response_type")

    ))

  }

  row <- as.data.frame(x)[1L, , drop = FALSE]

  family_name <- as.character(row$family[[1L]])

  row_meta <- as.list(row)

  row_meta$family <- NULL

  c(

    list(

      family_name = family_name,

      family = .margin_family_object(family_name)

    ),

    row_meta

  )

}


#' @export

best_fit.margin_screen <- best_fit.margin_selection


#' @export

best_fit_family.margin_selection <- function(x, ...) {

  best_fit(x)$family

}


#' @export

best_fit_family.margin_screen <- best_fit_family.margin_selection


#' @export

`$.margin_selection` <- function(x, name) {

  if (identical(name, "best_fit")) {

    return(best_fit(x))

  }

  .subset2(as.data.frame(x), name, exact = FALSE)

}


#' @export

`$.margin_screen` <- `$.margin_selection`


#' Select candidate marginal distributions

#'

#' `select_margin()` is a lightweight wrapper around [gamlss::fitDist()] for

#' the recommended longitudinal workflow: choose a plausible marginal family,

#' then fit dependence with [gamlss_longitudinal()].

#'

#' @param response Numeric response vector. For the common

#'   `select_margin(dat, response_var = "y")` call, a data frame supplied here

#'   is treated as `data`.

#' @param data Optional data frame containing the response.

#' @param response_var Optional response column name in `data`.

#' @param type Optional `gamlss::fitDist()` type. If `NULL`, a simple heuristic

#'   uses `"counts"` for non-negative integer responses, `"realplus"` for

#'   positive continuous responses, and `"realAll"` otherwise.

#' @param families Optional character vector used to filter the returned table.

#' @param time_intercepts Logical; if `TRUE`, rank retained families by a

#'   GAMLSS fit with factor time intercepts for each distribution parameter.

#'   This is useful for longitudinal screening where the final model is allowed

#'   to vary marginal location, scale, or shape over time. Finite-AIC fits are

#'   retained even if the temporary GAMLSS fit did not report convergence; check

#'   the returned `converged` column before selecting a final marginal family.

#' @param time_var Optional time column name in `data`, required when

#'   `time_intercepts = TRUE`.

#' @param try.gamlss,trace,... Passed to [gamlss::fitDist()].

#'

#' @return A data frame ordered by AIC and class `margin_selection`. With

#'   `time_intercepts = TRUE`, the table includes a `converged` column for the

#'   temporary time-intercept marginal fit. The selected family is also stored in

#'   the `"selected"` attribute for backward compatibility. Use [best_fit()] or

#'   [best_fit_family()] to extract the selected family in a fitting-friendly

#'   form.

#' @export

select_margin <- function(

  response = NULL,

  data = NULL,

  response_var = NULL,

  type = NULL,

  families = NULL,

  time_intercepts = FALSE,

  time_var = NULL,

  try.gamlss = FALSE,

  trace = FALSE,

  ...

) {

  if (is.null(data) && is.data.frame(response) && !is.null(response_var)) {

    data <- response

    response <- NULL

  }


  if (!is.null(data)) {

    data <- as.data.frame(data, stringsAsFactors = FALSE)

    if (is.null(response_var) || !is.character(response_var) || length(response_var) != 1L) {

      stop("'response_var' must be a single column name when 'data' is supplied.", call. = FALSE)

    }

    if (!response_var %in% names(data)) {

      stop("response_var='", response_var, "' not found in 'data'.", call. = FALSE)

    }

    response <- data[[response_var]]

  }


  if (is.null(response)) {

    stop("Supply either 'response' or both 'data' and 'response_var'.", call. = FALSE)

  }

  response_all <- as.numeric(response)

  if (isTRUE(time_intercepts)) {

    if (is.null(data)) {

      stop("'time_intercepts = TRUE' requires 'data' and 'time_var'.", call. = FALSE)

    }

    if (is.null(time_var) || !is.character(time_var) || length(time_var) != 1L) {

      stop("'time_var' must be a single column name when 'time_intercepts = TRUE'.", call. = FALSE)

    }

    if (!time_var %in% names(data)) {

      stop("time_var='", time_var, "' not found in 'data'.", call. = FALSE)

    }

  }

  response <- response_all[is.finite(response_all)]

  if (length(response) < 3L) {

    stop("Need at least three finite response values to screen margins.", call. = FALSE)

  }


  if (is.null(type)) {

    is_count <- all(response >= 0) && all(abs(response - round(response)) < .Machine$double.eps^0.5)

    type <- if (is_count) {

      "counts"

    } else if (all(response > 0)) {

      "realplus"

    } else {

      "realAll"

    }

  }


  fit <- NULL

  if (isTRUE(trace)) {

    fit <- gamlss::fitDist(

      response,

      type = type,

      try.gamlss = try.gamlss,

      trace = trace,

      ...

    )

  } else {

    invisible(utils::capture.output({

      fit <- suppressWarnings(suppressMessages(gamlss::fitDist(

        response,

        type = type,

        try.gamlss = try.gamlss,

        trace = trace,

        ...

      )))

    }))

  }


  if (is.null(fit$fits)) {

    stop("gamlss::fitDist() did not return a fits table.", call. = FALSE)

  }


  out <- data.frame(

    family = names(fit$fits),

    AIC = as.numeric(fit$fits),

    type = type,

    stringsAsFactors = FALSE,

    row.names = NULL

  )

  out <- out[is.finite(out$AIC), , drop = FALSE]

  if (!is.null(families)) {

    out <- out[out$family %in% families, , drop = FALSE]

  }


  out$screen_model <- "pooled_intercept"


  if (isTRUE(time_intercepts) && nrow(out) > 0L) {

    time_fits <- .select_margin_time_intercept_fits(

      response = response_all,

      time = data[[time_var]],

      families = out$family,

      trace = trace

    )

    out$pooled_AIC <- out$AIC

    out <- merge(

      out,

      time_fits,

      by = "family",

      all.x = TRUE,

      sort = FALSE

    )

    out$AIC <- out$time_intercept_AIC

    out$time_intercept_AIC <- NULL

    out$screen_model <- "time_intercepts"

    out <- out[is.finite(out$AIC), , drop = FALSE]

    if (nrow(out) > 0L && any(!out$converged, na.rm = TRUE)) {

      not_converged <- out$family[!out$converged]

      warning(

        "Time-intercept marginal screen retained finite-AIC fit(s) without confirmed convergence: ",

        paste(not_converged, collapse = ", "),

        ". Review the 'converged' column before selecting a final marginal family.",

        call. = FALSE

      )

    }

  }


  out <- out[order(out$AIC), , drop = FALSE]

  rownames(out) <- NULL

  out$rank <- seq_len(nrow(out))

  out$delta_AIC <- if (nrow(out) > 0L) out$AIC - min(out$AIC, na.rm = TRUE) else numeric(0)

  out$supported_by_longitudinal <- vapply(out$family, function(family) {

    all(vapply(

      paste0(c("d", "p", "q"), family),

      exists,

      logical(1),

      envir = asNamespace("gamlss.dist"),

      inherits = FALSE

    ))

  }, logical(1), USE.NAMES = FALSE)


  attr(out, "selected") <- if (nrow(out) > 0L) out$family[[1L]] else NA_character_

  attr(out, "response_type") <- type

  attr(out, "time_intercepts") <- isTRUE(time_intercepts)

  attr(out, "time_var") <- if (isTRUE(time_intercepts)) time_var else NULL

  class(out) <- c("margin_selection", "margin_screen", class(out))

  out

}


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


#' @rdname select_margin

#' @export

screen_margin <- function(...) {

  select_margin(...)

}


#' @export

print.margin_screen <- function(x, ..., n = 10L) {

  cat("\nMarginal Distribution Screen\n")

  cat("----------------------------\n")

  if (nrow(x) == 0L) {

    cat("No candidate families were retained.\n")

    return(invisible(x))

  }

  cat("Selected:", attr(x, "selected"), "\n\n")

  display <- utils::head(as.data.frame(x), n = n)

  unsupported <- "supported_by_longitudinal" %in% names(display) &

    any(!display$supported_by_longitudinal, na.rm = TRUE)

  display$supported_by_longitudinal <- NULL

  print(display, row.names = FALSE)

  if ("converged" %in% names(display) && any(!display$converged, na.rm = TRUE)) {

    cat("\nWarning: one or more printed time-intercept marginal fits did not report convergence.\n")

  }

  if (isTRUE(unsupported)) {

    cat("\nNote: one or more printed families are not currently supported by gamlss_longitudinal().\n")

  }

  invisible(x)

}


