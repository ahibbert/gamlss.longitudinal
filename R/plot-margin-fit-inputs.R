#' @keywords internal
#' @noRd
.plot_margin_resolve_family <- function(margin_dist) {
  if (inherits(margin_dist, "margin_selection") || inherits(margin_dist, "margin_screen")) {
    margin_dist <- best_fit_family(margin_dist)
  }
  if (is.character(margin_dist) && length(margin_dist) == 1L) {
    margin_dist <- .margin_family_object(margin_dist)
  }
  if (is.null(margin_dist) || is.null(margin_dist$family) || is.null(margin_dist$parameters)) {
    stop("'margin_dist' must be a gamlss.dist family object, family name, or margin_selection result.", call. = FALSE)
  }
  margin_dist
}

#' @keywords internal
#' @noRd
.plot_reject_old_args <- function(dots, old_args) {
  dot_names <- names(dots)
  old_used <- intersect(dot_names[!is.na(dot_names) & nzchar(dot_names)], old_args)
  if (length(old_used) > 0L) {
    replacements <- c(
      family = "margin_dist",
      dist = "margin_dist",
      copula = "copula_dist",
      selected_fit = "copula_dist"
    )
    msg <- vapply(old_used, function(arg) {
      paste0("'", arg, "' has been removed; use '", replacements[[arg]], "' instead")
    }, character(1), USE.NAMES = FALSE)
    stop(paste(msg, collapse = "; "), call. = FALSE)
  }
  if (length(dots) > 0L) {
    stop("Unused argument(s): ", paste(dot_names, collapse = ", "), call. = FALSE)
  }
  invisible(NULL)
}

#' @keywords internal
#' @noRd
.plot_reject_old_call_args <- function(call, old_args) {
  call_names <- names(as.list(call)[-1L])
  old_used <- intersect(call_names[!is.na(call_names) & nzchar(call_names)], old_args)
  if (length(old_used) == 0L) {
    return(invisible(NULL))
  }
  replacements <- c(
    family = "margin_dist",
    dist = "margin_dist",
    copula = "copula_dist",
    selected_fit = "copula_dist"
  )
  msg <- vapply(old_used, function(arg) {
    paste0("'", arg, "' has been removed; use '", replacements[[arg]], "' instead")
  }, character(1), USE.NAMES = FALSE)
  stop(paste(msg, collapse = "; "), call. = FALSE)
}

#' @keywords internal
#' @noRd
.plot_margin_family_name <- function(family) {
  if (is.character(family$family)) {
    family$family[1]
  } else {
    as.character(family$family[[1L]])
  }
}

#' @keywords internal
#' @noRd
.plot_margin_response <- function(x = NULL, data = NULL, response_var = "response") {
  if (!is.null(data)) {
    data <- as.data.frame(data, stringsAsFactors = FALSE)
    if (!is.character(response_var) || length(response_var) != 1L || !response_var %in% names(data)) {
      stop("'response_var' must name a response column in 'data'.", call. = FALSE)
    }
    return(as.numeric(data[[response_var]]))
  }
  if (is.data.frame(x)) {
    if (!is.character(response_var) || length(response_var) != 1L || !response_var %in% names(x)) {
      stop("'response_var' must name a response column in 'x'.", call. = FALSE)
    }
    return(as.numeric(x[[response_var]]))
  }
  as.numeric(x)
}
