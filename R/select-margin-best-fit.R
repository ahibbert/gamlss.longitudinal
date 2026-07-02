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
