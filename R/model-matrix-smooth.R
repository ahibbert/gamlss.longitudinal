#' Build smooth-term model matrices for one parameter(e.g. mu)
#'
#' @keywords internal
#' @noRd
.gl_build_smooth_model_matrices <- function(formula, data_for_par, smooth_eval_env) {
  formula_terms <- stats::terms(formula)
  term_labels <- attr(formula_terms, "term.labels")
  smooth_terms <- term_labels[grepl("^\\s*s\\(", term_labels)]

  if (length(smooth_terms) == 0) {
    return(NULL)
  }

  smooth_matrices <- list()
  for (s_label in smooth_terms) {
    s_txt <- trimws(s_label)
    s_call <- tryCatch(parse(text = s_txt)[[1]], error = function(e) NULL)
    s_obj <- eval(parse(text = s_txt), envir = smooth_eval_env)
    s_con <- mgcv::smoothCon(s_obj, data = data_for_par, knots = NULL, absorb.cons = TRUE)

    if (length(s_con) > 0 && !is.null(s_con[[1]]$X)) {
      B_s <- s_con[[1]]$X

      if (!is.null(s_con[[1]]$S) && length(s_con[[1]]$S) > 0) {
        attr(B_s, "penalty") <- s_con[[1]]$S[[1]]
      }

      if (!is.null(s_call) && length(s_call) >= 2) {
        x_expr <- s_call[[2]]
        x_var <- trimws(gsub("`", "", paste(deparse(x_expr), collapse = " "), fixed = TRUE))
        x_value <- tryCatch(eval(x_expr, envir = data_for_par, enclos = parent.frame()), error = function(e) NULL)
        if (!is.null(x_value) && length(x_value) == nrow(B_s)) {
          attr(B_s, "smooth_x") <- as.numeric(x_value)
          attr(B_s, "smooth_var") <- x_var
        } else if (nzchar(x_var)) {
          attr(B_s, "smooth_var") <- x_var
        }
      }

      smooth_matrices[[s_txt]] <- B_s
    }
  }

  if (length(smooth_matrices) == 0) {
    return(NULL)
  }

  smooth_matrices
}
