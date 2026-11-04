.select_copula_margin_family <- function(
    margin_dist = NULL,
    object = NULL,
    margin_selection = NULL) {
  if (!is.null(object) && !is.null(object$margin_dist)) {
    return(.gl_capability_margin_code(object$margin_dist))
  }
  if (inherits(margin_dist, "margin_selection")) {
    selected <- best_fit(margin_dist)$family_name
    return(if (is.null(selected) || is.na(selected)) NA_character_ else selected)
  }
  if (!is.null(margin_dist)) {
    return(.gl_capability_margin_code(margin_dist))
  }
  if (inherits(margin_selection, "margin_selection")) {
    selected <- best_fit(margin_selection)$family_name
    return(if (is.null(selected) || is.na(selected)) NA_character_ else selected)
  }
  NA_character_
}
