#' Runs a separately optimised RS fit to provide safe starting values for a joint RS fit 
#' and extracts fit information needed to pass to start for the joint fit.
#'
#' @noRd
.gl_fitted_parameter_eta <- function(fit, par_name) {
  old_x <- fit$model_matrix$x[[par_name]]
  if (is.null(old_x)) return(NULL)
  old_x <- as.matrix(old_x)
  fixed_names <- paste(par_name, colnames(old_x), sep = ".")
  if (!all(fixed_names %in% names(fit$par))) return(NULL)
  eta <- as.numeric(old_x %*% as.numeric(fit$par[fixed_names]))
  old_s <- fit$model_matrix$s[[par_name]] %||% list()
  old_beta_s <- fit$par_s[[par_name]] %||% list()
  for (smooth_name in intersect(names(old_s), names(old_beta_s))) {
    basis <- as.matrix(old_s[[smooth_name]])
    beta <- as.numeric(old_beta_s[[smooth_name]])
    if (nrow(basis) != length(eta) || ncol(basis) != length(beta)) return(NULL)
    eta <- eta + as.numeric(basis %*% beta)
  }
  eta
}

#' Project a fitted linear predictor into a candidate model design
#'
#' Backward selection can remove a factor whose coefficients offset a large
#' intercept. Copying only retained coefficient names then creates an extreme,
#' sometimes invalid starting state. Least-squares projection preserves as much
#' of the fitted parent predictor as the reduced candidate design can represent.
#'
#' @noRd
.gl_project_fitted_start <- function(fit, target_mm, lm_fit_fn = stats::lm.fit) {
  if (!inherits(fit, "gamlss.longitudinal") ||
      is.null(fit$model_matrix$x) || is.null(target_mm$x)) {
    return(list(
      par = fit$par,
      par_s = fit$par_s,
      lambda_s = fit$lambda_s,
      projected = character(),
      rmse = numeric()
    ))
  }

  projected_par <- fit$par
  projected_par_s <- fit$par_s %||% list()
  projected <- character()
  projection_rmse <- numeric()

  for (par_name in intersect(names(target_mm$x), names(fit$model_matrix$x))) {
    target_x <- as.matrix(target_mm$x[[par_name]])
    target_s <- target_mm$s[[par_name]] %||% list()
    old_x <- as.matrix(fit$model_matrix$x[[par_name]])
    old_s <- fit$model_matrix$s[[par_name]] %||% list()
    unchanged <- identical(colnames(target_x), colnames(old_x)) &&
      identical(names(target_s), names(old_s)) &&
      all(vapply(names(target_s), function(name) {
        identical(dim(target_s[[name]]), dim(old_s[[name]]))
      }, logical(1L)))
    if (isTRUE(unchanged)) next

    old_eta <- .gl_fitted_parameter_eta(fit, par_name)
    design_parts <- c(list(target_x), lapply(target_s, as.matrix))
    target_design <- do.call(cbind, design_parts)
    if (is.null(old_eta) || nrow(target_design) != length(old_eta) ||
        !all(is.finite(old_eta)) || !all(is.finite(target_design))) next

    fitted_projection <- tryCatch(
      lm_fit_fn(x = target_design, y = old_eta),
      error = identity
    )
    if (inherits(fitted_projection, "error") ||
        length(fitted_projection$coefficients) != ncol(target_design)) next
    coefficients <- as.numeric(fitted_projection$coefficients)
    coefficients[!is.finite(coefficients)] <- 0

    fixed_count <- ncol(target_x)
    fixed_names <- paste(par_name, colnames(target_x), sep = ".")
    projected_par[fixed_names] <- coefficients[seq_len(fixed_count)]
    offset <- fixed_count
    projected_par_s[[par_name]] <- list()
    for (smooth_name in names(target_s)) {
      width <- ncol(target_s[[smooth_name]])
      indices <- offset + seq_len(width)
      beta <- coefficients[indices]
      names(beta) <- paste(par_name, smooth_name, seq_len(width), sep = ".")
      projected_par_s[[par_name]][[smooth_name]] <- beta
      offset <- offset + width
    }
    fitted_eta <- as.numeric(target_design %*% coefficients)
    projection_rmse[[par_name]] <- sqrt(mean((old_eta - fitted_eta)^2))
    projected <- c(projected, par_name)
  }

  list(
    par = projected_par,
    par_s = projected_par_s,
    lambda_s = fit$lambda_s,
    projected = projected,
    rmse = projection_rmse
  )
}

#' @noRd
.gl_run_joint_warm_start <- function(
    start_from,
    method,
    include_dlcopdpar,
    warm_start_joint,
    warm_start_joint_iter,
    user_supplied_start,
    dataset_original,
    margin_dist,
    copula_dist,
    time_var,
    subject_var,
    missingness,
    mu.formula,
    sigma.formula,
    nu.formula,
    tau.formula,
    theta.formula,
    zeta.formula,
    optimizer_control,
    true_val,
    vcov_method,
    vcov_numderiv,
    use_Rcpp,
    lambda_start,
    lambda_penalty_K,
    verbose,
    target_mm = NULL,
    fit_fn = NULL,
    warm_start_fit_call_fn = .gl_call_joint_warm_start_fit) {
  warm_start_info <- .gl_default_warm_start_info()
  warm_start_par_s <- NULL
  warm_start_lambda_s <- NULL
  projected_start_fallback <- NULL

  if (inherits(start_from, "gamlss.longitudinal") ||
      (is.list(start_from) && !is.null(start_from$par))) {
    fitted_start <- start_from
    projected <- if (inherits(start_from, "gamlss.longitudinal") &&
        !is.null(target_mm)) {
      .gl_project_fitted_start(start_from, target_mm)
    } else {
      NULL
    }
    if (!is.null(projected) && length(projected$projected)) {
      projected_start_fallback <- projected
      warm_start_info$projected_parameters <- projected$projected
      warm_start_info$projection_rmse <- projected$rmse
      warm_start_info$projection_used <- FALSE
    }
    warm_start_par_s <- fitted_start$par_s
    warm_start_lambda_s <- fitted_start$lambda_s
    start_from <- fitted_start$par
  }

  if (!.gl_should_run_joint_warm_start(
    method = method,
    include_dlcopdpar = include_dlcopdpar,
    warm_start_joint = warm_start_joint,
    warm_start_joint_iter = warm_start_joint_iter,
    user_supplied_start = user_supplied_start
  )) {
    return(list(
      start_from = start_from,
      warm_start_par_s = warm_start_par_s,
      warm_start_lambda_s = warm_start_lambda_s,
      projected_start_fallback = projected_start_fallback,
      warm_start_info = warm_start_info
    ))
  }

  if (is.null(fit_fn)) {
    fit_fn <- gamlss_longitudinal
  }

  if (verbose > 0) {
    cat(
      "\nRunning separate RS warm-start phase for ",
      warm_start_joint_iter,
      " outer iteration(s) before joint RS fit...\n",
      sep = ""
    )
  }

  warm_call <- warm_start_fit_call_fn(
    fit_fn = fit_fn,
    dataset_original = dataset_original,
    margin_dist = margin_dist,
    copula_dist = copula_dist,
    time_var = time_var,
    subject_var = subject_var,
    missingness = missingness,
    mu.formula = mu.formula,
    sigma.formula = sigma.formula,
    nu.formula = nu.formula,
    tau.formula = tau.formula,
    theta.formula = theta.formula,
    zeta.formula = zeta.formula,
    optimizer_control = optimizer_control,
    true_val = true_val,
    method = method,
    warm_start_joint_iter = warm_start_joint_iter,
    vcov_method = vcov_method,
    vcov_numderiv = vcov_numderiv,
    use_Rcpp = use_Rcpp,
    lambda_start = lambda_start,
    lambda_penalty_K = lambda_penalty_K
  )
  warm_fit <- warm_call$fit
  warm_output <- warm_call$output
  warm_warnings <- warm_call$warnings
  warm_err <- warm_call$error

  .gl_validate_joint_warm_start_fit(warm_fit, warm_err)
  warm_start_result <- .gl_joint_warm_start_result(
    warm_fit = warm_fit,
    warm_output = warm_output,
    warm_warnings = warm_warnings,
    warm_start_joint_iter = warm_start_joint_iter
  )

  if (verbose > 1 && length(warm_output) > 0) {
    cat(paste(warm_output, collapse = "\n"), "\n")
  }

  warm_start_result
}
