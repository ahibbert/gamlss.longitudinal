.gl_prepare_newdata_internal <- function(object, newdata, require_response = FALSE) {
  if (is.null(newdata)) {
    return(NULL)
  }

  if (is.null(object$formulas_int) || is.null(object$var_map)) {
    stop("newdata prediction requires a model fit created with stored formulas/variable map. Refit with the current package version.")
  }

  nd <- as.data.frame(newdata, stringsAsFactors = FALSE)

  nd <- .gl_translate_newdata_names(nd, object$var_map)

  nd <- .gl_add_newdata_default_columns(nd, object)

  nd <- .gl_align_newdata_factor_levels(nd, object)

  .gl_validate_newdata_response(nd, require_response)

  nd
}
