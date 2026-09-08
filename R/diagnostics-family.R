.gl_get_family_fun <- function(family_name, prefix) {
  fun_name <- paste0(prefix, family_name)

  if (exists(fun_name, envir = asNamespace("gamlss.dist"), mode = "function", inherits = FALSE)) {
    return(get(fun_name, envir = asNamespace("gamlss.dist"), mode = "function", inherits = FALSE))
  }

  if (exists(fun_name, mode = "function")) {
    return(get(fun_name, mode = "function"))
  }

  stop(
    "Distribution function '", fun_name,
    "' is not available in gamlss.dist or the current session."
  )
}

.gl_call_family_fun <- function(prefix, family_name, x, params, extra_args = list()) {
  fun <- .gl_get_family_fun(family_name, prefix)

  arg_name <- switch(prefix,
    p = "q",
    d = "x",
    q = "p",
    "x"
  )

  family_args <- attr(params, "family_args", exact = TRUE) %||% list()
  args <- c(stats::setNames(list(x), arg_name), extra_args, family_args, params)

  args <- args[names(args) %in% formalArgs(fun)]

  do.call(fun, args)
}

.gl_subset_family_params <- function(params, index) {
  family_args <- attr(params, "family_args", exact = TRUE) %||% list()
  out <- lapply(params, function(value) value[index])
  family_args <- lapply(family_args, function(value) value[index])
  attr(out, "family_args") <- family_args
  out
}
