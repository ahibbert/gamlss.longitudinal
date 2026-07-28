#' @keywords internal
#' @noRd
.plot_margin_empty_density <- function(group = "All") {
  data.frame(
    response = numeric(),
    density = numeric(),
    split_group = as.character(group)[0L],
    stringsAsFactors = FALSE
  )
}

#' @keywords internal
#' @noRd
.plot_margin_support_bounds <- function(family, params = NULL) {
  family_name <- .plot_margin_family_name(family)
  qfun <- tryCatch(
    get(paste0("q", family_name), envir = asNamespace("gamlss.dist"), inherits = FALSE),
    error = function(e) NULL
  )
  if (!is.function(qfun)) {
    return(c(lower = -Inf, upper = Inf))
  }

  q_args <- list()
  for (par_name in names(family$parameters)) {
    par_value <- NULL
    if (!is.null(params) && par_name %in% names(params)) {
      par_vec <- as.numeric(params[[par_name]])
      par_vec <- par_vec[is.finite(par_vec)]
      if (length(par_vec) > 0L) {
        par_value <- stats::median(par_vec)
      }
    }
    if (is.null(par_value)) {
      par_value <- tryCatch(
        eval(formals(qfun)[[par_name]], envir = baseenv()),
        error = function(e) NA_real_
      )
      par_value <- as.numeric(par_value)[1L]
    }
    if (is.finite(par_value)) {
      q_args[[par_name]] <- par_value
    }
  }

  bounds <- tryCatch(
    do.call(qfun, c(list(p = c(0, 1)), q_args)),
    error = function(e) {
      tryCatch(
        do.call(qfun, c(list(p = c(.Machine$double.eps, 1 - .Machine$double.eps)), q_args)),
        error = function(e2) c(-Inf, Inf)
      )
    }
  )
  bounds <- as.numeric(bounds)
  if (length(bounds) < 2L) {
    return(c(lower = -Inf, upper = Inf))
  }
  c(
    lower = if (is.finite(bounds[1L])) bounds[1L] else -Inf,
    upper = if (is.finite(bounds[2L])) bounds[2L] else Inf
  )
}
