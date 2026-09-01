#' Print a summary output for the completed fit, 
#' including covariate estimates and model selection criteria.
#'
#' @noRd
.gl_print_fit_summary <- function(par_cov, aics, dataset, margin_dist, copula_dist, total_fit_time) {
  cat("\n\n############ MODEL FIT ############\n")
  cat(paste("\nMargin distribution:", margin_dist$family[2]))
  cat(paste("\nCopula distribution:", copula_dist))
  cat("\n")
  cat(paste("\nParameter count:", length(par_cov)))
  cat(paste("\nObservations:", nrow(dataset)))
  cat(paste("\nMargins:", length(unique(dataset$time))))
  cat("\n")
  cat(paste("\nTotal time (seconds):", round(total_fit_time, 2)))
  cat("\n\n")

  par_mat_out_temp <- t(t((par_cov)))
  colnames(par_mat_out_temp) <- c("estimate")
  print(par_mat_out_temp)

  cat("\n")
  cat("Model Selection Criteria:")
  cat("\n")
  print(aics)
  cat("\n####################################\n")
}

#' Kickoff variance / covariance calculation and attach to fitted object if requested (it's on by default)
#'
#' @noRd
.gl_attach_fit_vcov <- function(return_list, compute_vcov, vcov_numderiv, vcov_method, verbose) {
  return_list$vcov <- NULL
  return_list$vcov_meta <- list(
    precomputed = FALSE,
    numderiv = isTRUE(vcov_numderiv),
    method = vcov_method,
    cache_version = .gl_vcov_cache_version()
  )

  objective <- return_list$likelihood_contract$objective %||% "ordinary"
  if (identical(objective, "segmented") &&
      !identical(vcov_method, "sandwich")) {
    return_list$vcov_meta$inference_status <- "unavailable_for_segmented_objective"
    return_list$vcov_meta$hessian_diagnostics <- list(
      status = "unavailable",
      reason = "segmented_objective",
      recommendation = "Use cluster-sandwich or bootstrap sensitivity analysis."
    )
    return(return_list)
  }

  if (isTRUE(compute_vcov)) {
    if (verbose > 0) {
      cat("Calculating variance-covariance matrix at fit completion...\n")
    }

    vcov_cached <- NULL
    inference_failure <- NULL
    vcov_cached <- tryCatch(
      {
        vcov.gamlss.longitudinal(
          return_list,
          numderiv = isTRUE(vcov_numderiv),
          method = vcov_method,
          details = TRUE,
          progress = isTRUE(verbose > 0)
        )
      },
      error = function(e) {
        if (inherits(e, "gamlss_longitudinal_inference_unavailable")) {
          inference_failure <<- e$diagnostics
          warning(structure(
            list(message = conditionMessage(e), call = NULL,
                 diagnostics = e$diagnostics),
            class = c("gamlss_longitudinal_inference_unavailable",
                      "warning", "condition")
          ))
        } else {
          warning(
            "Could not precompute variance-covariance matrix at fit completion: ",
            conditionMessage(e),
            call. = FALSE
          )
        }
        NULL
      }
    )

    if (!is.null(vcov_cached)) {
      return_list$vcov <- vcov_cached
      return_list$vcov_meta$precomputed <- TRUE

      if (!is.null(vcov_cached$method)) {
        return_list$vcov_meta$method_used <- vcov_cached$method
      }

      if (!is.null(vcov_cached$method_requested)) {
        return_list$vcov_meta$method <- vcov_cached$method_requested
      }
      return_list$vcov_meta$inference_contract <- vcov_cached$inference_contract
      return_list$vcov_meta$smooth_inference_contract <- vcov_cached$smooth_inference_contract
    } else if (!is.null(inference_failure)) {
      return_list$vcov_meta$inference_status <- "unavailable"
      return_list$vcov_meta$hessian_diagnostics <- inference_failure
    }
  }

  return_list
}
