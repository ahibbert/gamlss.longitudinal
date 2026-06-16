#' Standardise formula input to a formula object, ensuring it has the expected structure
#'
#' @noRd
.gl_normalize_formula_input <- function(fml, response_name = "response", require_lhs = FALSE) {
  if (inherits(fml, "formula")) {
    return(fml)
  }

  if (!is.character(fml) || length(fml) != 1 || is.na(fml) || nchar(trimws(fml)) == 0) {
    stop("ERROR: Invalid formula input: ", deparse(fml))
  }

  txt <- trimws(fml)

  if (!grepl("~", txt, fixed = TRUE)) {
    txt <- if (require_lhs) paste0(response_name, " ~ ", txt) else paste0("~ ", txt)
  } else if (require_lhs) {
    parts <- strsplit(txt, "~", fixed = TRUE)[[1]]
    lhs <- trimws(parts[1])
    rhs <- trimws(parts[2])
    if (nchar(lhs) == 0) txt <- paste0(response_name, " ~ ", rhs)
  }

  as.formula(txt, env = parent.frame())
}

#' Translate formula variable names from user columns to internal columns
#'
#' @noRd
.gl_translate_formula_vars <- function(fml, var_map, response_name = "response", require_lhs = FALSE) {
  f_obj <- .gl_normalize_formula_input(fml, response_name = response_name, require_lhs = require_lhs)
  f_txt <- paste(deparse(f_obj), collapse = " ")

  for (old_name in names(var_map)) {
    new_name <- var_map[[old_name]]
    if (!identical(old_name, new_name)) {
      old_esc <- gsub("([][{}()+*^$|\\\\.?])", "\\\\\\1", old_name)
      f_txt <- gsub(paste0("`", old_esc, "`"), new_name, f_txt, perl = TRUE)
      f_txt <- gsub(paste0("\\b", old_esc, "\\b"), new_name, f_txt, perl = TRUE)
    }
  }

  as.formula(f_txt, env = environment(f_obj))
}

#' Translate all fitting formulas for the internal fitting process
#'
#' @noRd
.gl_translate_fit_formulas <- function(
    mu.formula,
    sigma.formula,
    nu.formula,
    tau.formula,
    theta.formula,
    zeta.formula,
    formula_var_map,
    response_name = "response") {
  list(
    mu = .gl_translate_formula_vars(mu.formula, formula_var_map, response_name = response_name, require_lhs = TRUE),
    sigma = .gl_translate_formula_vars(sigma.formula, formula_var_map, response_name = response_name, require_lhs = FALSE),
    nu = .gl_translate_formula_vars(nu.formula, formula_var_map, response_name = response_name, require_lhs = FALSE),
    tau = .gl_translate_formula_vars(tau.formula, formula_var_map, response_name = response_name, require_lhs = FALSE),
    theta = .gl_translate_formula_vars(theta.formula, formula_var_map, response_name = response_name, require_lhs = FALSE),
    zeta = .gl_translate_formula_vars(zeta.formula, formula_var_map, response_name = response_name, require_lhs = FALSE)
  )
}
