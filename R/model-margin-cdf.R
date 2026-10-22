#' Classed marginal-CDF evaluation failure
#'
#' @keywords internal
#' @noRd
.gl_margin_cdf_error <- function(message, margin_dist, rows = integer()) {
  family <- as.character(margin_dist$family[1L])
  structure(
    list(
      message = paste0("Marginal CDF evaluation failed for family ", family, ": ", message),
      call = NULL,
      family = family,
      rows = as.integer(rows)
    ),
    class = c("gamlss_longitudinal_margin_cdf_error", "error", "condition")
  )
}

#' Build complete marginal-CDF arguments
#'
#' @keywords internal
#' @noRd
.gl_margin_cdf_input <- function(eta_inv, margin_dist, response) {
  n <- length(response)
  input <- list(y = response, q = response, x = response)

  for (par_name in intersect(names(eta_inv), c("mu", "sigma", "nu", "tau"))) {
    value <- eta_inv[[par_name]]
    if (!length(value) %in% c(1L, n)) {
      stop(.gl_margin_cdf_error(
        paste0("parameter '", par_name, "' has length ", length(value),
               " but expected 1 or ", n, "."),
        margin_dist
      ))
    }
    input[[par_name]] <- value
  }

  fixed_unlinked_values <- attr(margin_dist, "fixed_unlinked_values")
  if (length(fixed_unlinked_values) > 0L) {
    for (par_name in names(fixed_unlinked_values)) {
      value <- fixed_unlinked_values[[par_name]]
      if (is.numeric(value) && length(value) == 1L && is.finite(value)) {
        input[[par_name]] <- rep(value, n)
      }
    }
  }

  c(input, .gl_margin_fixed_family_args(margin_dist, n))
}

#' Retry a marginal CDF one observed row at a time
#'
#' @keywords internal
#' @noRd
.gl_margin_cdf_scalar_eval <- function(FUN, input, FUN_args, n) {
  call_args <- input[FUN_args]
  lengths <- vapply(call_args, length, integer(1))
  if (length(lengths) && !all(lengths %in% c(1L, n))) {
    return(rep(NA_real_, n))
  }
  expanded <- lapply(call_args, rep, length.out = n)
  vapply(seq_len(n), function(ii) {
    scalar_args <- lapply(expanded, `[`, ii)
    value <- tryCatch(do.call(FUN, scalar_args), error = function(e) NA_real_)
    as.numeric(value)[1L]
  }, numeric(1))
}

#' Evaluate a marginal CDF without passing missing responses to family code
#'
#' Missing response rows remain `NA`. For a discrete margin, finite negative
#' responses represent lower rectangle endpoints and are returned as zero.
#'
#' @keywords internal
#' @noRd
.gl_evaluate_margin_cdf <- function(
    input,
    margin_dist,
    response,
    FUN = NULL,
    FUN_args = NULL,
    cacheable = FALSE,
    cache_env = NULL,
    cache_prefix = "",
    use_fast = FALSE) {
  n <- length(response)
  if (is.null(FUN)) {
    FUN <- get(
      paste0("p", margin_dist$family[1L]),
      envir = asNamespace("gamlss.dist"),
      mode = "function",
      inherits = FALSE
    )
  }
  if (is.null(FUN_args)) FUN_args <- formalArgs(FUN)

  out <- rep(NA_real_, n)
  observed <- is.finite(response)
  negative_discrete <- observed & .is_discrete_margin(margin_dist) & response < 0
  out[negative_discrete] <- 0
  evaluate <- observed & !negative_discrete
  if (!any(evaluate)) return(out)

  observed_input <- lapply(input, function(value) {
    if (length(value) == n) value[evaluate] else value
  })
  observed_args <- names(observed_input)[names(observed_input) %in% FUN_args]
  n_evaluate <- sum(evaluate)

  value <- NULL
  if (isTRUE(use_fast)) {
    value <- .call_fast_count_family("p", margin_dist$family[1L], observed_input)
  }
  if (is.null(value)) {
    value <- .call_margin_family_cached(
      FUN,
      observed_input,
      observed_args,
      cacheable = cacheable,
      cache_env = cache_env,
      cache_prefix = cache_prefix,
      force_scalar = .gl_margin_requires_scalar_calls(margin_dist$family[1L])
    )
  }
  value <- as.numeric(value)

  invalid_length <- length(value) != n_evaluate
  invalid_value <- invalid_length || any(!is.finite(value)) ||
    any(value < -100 * .Machine$double.eps | value > 1 + 100 * .Machine$double.eps)
  if (invalid_value) {
    value <- .gl_margin_cdf_scalar_eval(FUN, observed_input, observed_args, n_evaluate)
  }

  bad <- !is.finite(value) |
    value < -100 * .Machine$double.eps |
    value > 1 + 100 * .Machine$double.eps
  if (length(value) != n_evaluate || any(bad)) {
    evaluated_rows <- which(evaluate)
    bad_rows <- if (length(value) == n_evaluate) evaluated_rows[bad] else evaluated_rows
    stop(.gl_margin_cdf_error(
      paste0(
        "nonfinite or out-of-range values remained on observed row(s) ",
        paste(utils::head(bad_rows, 10L), collapse = ", "),
        if (length(bad_rows) > 10L) ", ..." else "."
      ),
      margin_dist,
      bad_rows
    ))
  }

  out[evaluate] <- pmin(pmax(value, 0), 1)
  out
}

#' Validate a CDF derivative while zeroing only unobserved rows
#'
#' @keywords internal
#' @noRd
.gl_validate_margin_cdf_derivative <- function(
    value,
    response,
    margin_dist,
    context,
    parameter = NULL) {
  value <- as.numeric(value)
  observed <- is.finite(response)
  if (length(value) != length(response)) {
    stop(.gl_margin_cdf_error(
      paste0(context, " returned ", length(value), " values for ",
             length(response), " response rows."),
      margin_dist
    ))
  }
  bad <- observed & !is.finite(value)
  if (any(bad)) {
    label <- if (is.null(parameter)) "" else paste0(" for parameter '", parameter, "'")
    stop(.gl_margin_cdf_error(
      paste0(context, label, " was nonfinite on observed row(s) ",
             paste(utils::head(which(bad), 10L), collapse = ", "),
             if (sum(bad) > 10L) ", ..." else "."),
      margin_dist,
      which(bad)
    ))
  }
  value[!observed] <- 0
  value
}

#' Marginal CDF evaluated on the model's natural parameter scale
#'
#' @keywords internal
#' @noRd
calc_F_x <- function(eta_inv, mm, margin_dist, response) {
  input <- .gl_margin_cdf_input(eta_inv, margin_dist, response)
  .gl_evaluate_margin_cdf(
    input = input,
    margin_dist = margin_dist,
    response = response,
    cacheable = .is_discrete_margin(margin_dist),
    cache_prefix = paste0(margin_dist$family[1L], ":cdf")
  )
}
