#' Solve a Hessian into covariance output and diagnostics
#'
#' @noRd
.gl_solve_hessian_vcov <- function(H) {
  vc <- -solve(H)

  se <- sqrt(abs(diag(solve(H))))

  if (!is.matrix(vc) || any(!is.finite(vc)) || any(!is.finite(se))) {
    stop("Hessian inversion produced non-finite variance-covariance values.", call. = FALSE)
  }

  eig <- tryCatch(eigen((H + t(H)) / 2, symmetric = TRUE, only.values = TRUE)$values, error = function(e) NA_real_)

  list(
    vcov = vc,
    se = se,
    hessian_diagnostics = list(
      condition_number = tryCatch(kappa(H), error = function(e) NA_real_),
      min_abs_eigen = if (any(is.finite(eig))) min(abs(eig[is.finite(eig)])) else NA_real_,
      max_abs_eigen = if (any(is.finite(eig))) max(abs(eig[is.finite(eig)])) else NA_real_
    )
  )
}

#' Source analytical Hessian helpers for standalone development use
#'
#' @noRd
.gl_source_analytical_hessian_helpers <- function() {
  if (exists("calc_analytical_hessian", mode = "function")) {
    return(invisible(TRUE))
  }

  hessian_files <- c(
    "hessian-linkinv-derivatives.R",
    "hessian-fd-step.R",
    "hessian-warnings.R",
    "hessian-margin-cdf.R",
    "hessian-copula.R",
    "hessian-margin-derivatives.R",
    "hessian-assembly.R",
    "hessian-analytical.R"
  )

  candidate_dirs <- c(
    dirname(attr(body(vcov.gamlss.longitudinal), "srcfile")$filename %||% ""),
    "R",
    file.path(getwd(), "R")
  )

  loaded <- FALSE

  for (dir in candidate_dirs) {
    candidate_paths <- file.path(dir, hessian_files)

    if (all(file.exists(candidate_paths))) {
      for (cp in candidate_paths) source(cp, local = FALSE)
      loaded <- TRUE
      break
    }
  }

  if (!loaded) stop("Cannot locate hessian helper files. Source them manually or ensure the working directory is the package root.")

  invisible(TRUE)
}
