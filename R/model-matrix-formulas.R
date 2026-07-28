#' Fix time-covariate column names produced by contrast expansion
#' 
#' This function normalizes column names for time covariates that are produced by contrast expansion in model matrices.
#' The function looks for column names that match the pattern "time_covariate.L", "time_covariate.Q", "time_covariate.C" 
#' and replaces them with "time_covariate.1", "time_covariate.2", "time_covariate.3" respectively. 
#' 
#' This is essentially a bug fix for the naming continuing to show the ordering in a human readable way.
#'
#' @noRd
.gl_normalize_time_covariate_colnames <- function(nms) {
  if (length(nms) == 0) {
    return(nms)
  }

  out <- nms
  suffix_map <- c(L = "1", Q = "2", C = "3")

  for (sx in names(suffix_map)) {
    out <- gsub(

      paste0("time_covariate\\.", sx, "\\b"),
      paste0("time_covariate.", suffix_map[[sx]]),
      out,
      perl = TRUE
    )
  }

  # contr.poly names can appear as ^4, ^5, ...; normalize to .4, .5, ...
  out <- gsub("time_covariate\\^([0-9]+)", "time_covariate.\\1", out, perl = TRUE)

  out
}

#' Convert a one-sided parameter formula to a response formula
#' 
#' This is so the user doesn't have to specify a response variable 
#' in the formula for parameters other than mu, since the response variable 
#' is not actually used in model matrix construction. It's just needed 
#' for how we initialise the formula object for model matrix construction.
#'
#' @noRd
.gl_to_response_formula <- function(fml, response_name = "response") {
  if (inherits(fml, "formula")) {
    rhs_txt <- if (length(fml) == 3L) {
      paste(deparse(fml[[3]]), collapse = " ")
    } else {
      paste(deparse(fml[[2]]), collapse = " ")
    }
  } else if (is.character(fml) && length(fml) == 1L) {
    txt <- trimws(fml)

    if (grepl("~", txt, fixed = TRUE)) {
      parts <- strsplit(txt, "~", fixed = TRUE)[[1]]

      rhs_txt <- trimws(parts[length(parts)])
    } else {
      rhs_txt <- txt
    }
  } else {
    stop("Invalid formula input: ", deparse(fml))
  }

  as.formula(paste(response_name, "~", rhs_txt), env = parent.frame())
}

#' Determines if we have a one or two parameter copula based on the copula family, 
#' and returns the names of the parameters to be included in model matrix construction.
#' 
#' We can already extract the margin parameters from the margin.family object, 
#' so we just need to determine if we have one or two copula parameters based on the copula family.
#'
#' @noRd
.gl_model_matrix_included_parameters <- function(margin.family, copula.family) {
  if (copula.family %in% c("t", "T", "Student")) {
    two_par_cop <- TRUE
  } else {
    two_par_cop <- FALSE
  }

  c(names(margin.family$parameters), if (two_par_cop) c("theta", "zeta") else c("theta"))
}
