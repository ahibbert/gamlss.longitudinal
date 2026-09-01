#' Calculate smooth-term variance components for vcov output
#'
#' @noRd
.gl_compute_smooth_vcov <- function(object, eta_inv, response) {
  smooth_vcov_list <- list()

  smooth_se_list <- list()

  # Extract residual variance estimate (using reciprocal of mean weights as proxy for sigma^2)

  if (!is.null(object$weights) && length(object$weights) > 0 && is.numeric(object$weights)) {
    sigma2_est <- 1 / mean(object$weights, na.rm = TRUE)
  } else {
    # Fallback: estimate from residuals if weights not available

    fitted_response <- eta_inv[["mu"]]

    residuals <- response - fitted_response

    sigma2_est <- var(residuals, na.rm = TRUE)
  }

  # Ensure scalar numeric to avoid deprecated array recycling warnings.

  sigma2_est <- as.numeric(sigma2_est)[1]

  if (!is.finite(sigma2_est)) {
    sigma2_est <- 1
  }

  # Process each parameter that has smooth terms

  for (par_name in names(object$par_s)) {
    if (length(object$par_s[[par_name]]) > 0) {
      smooth_vcov_list[[par_name]] <- list()

      smooth_se_list[[par_name]] <- list()

      # Process each smooth term for this parameter

      for (s_name in names(object$par_s[[par_name]])) {
        # Get the B-spline basis matrix

        B <- object$model_matrix$s[[par_name]][[s_name]]

        # Get the smoothing parameter

        lambda <- object$lambda_s[[par_name]][[s_name]]

        # Use the mgcv-generated penalty stored on the basis matrix; fall back to

        # a generic second-difference penalty only when unavailable.

        k <- ncol(B)

        pen_attr <- attr(B, "penalty")

        if (!is.null(pen_attr) && is.matrix(pen_attr) &&

          nrow(pen_attr) == k && ncol(pen_attr) == k) {
          P <- pen_attr
        } else if (k > 2) {
          D2 <- diff(diag(k), differences = 2)

          P <- t(D2) %*% D2
        } else {
          P <- diag(k)
        }

        # Get per-parameter IRLS working weights. object$weights is a named list

        # keyed by parameter name; fall back to unit weights if not available.

        w_par <- object$weights[[par_name]]

        if (!is.null(w_par) && is.numeric(w_par) && length(w_par) == nrow(B)) {
          w_diag <- as.vector(w_par)
        } else {
          w_diag <- rep(1, nrow(B))
        }

        W <- diag(w_diag)

        # Per-parameter sigma2: scale consistent with IRLS, 1/mean(w)

        sigma2_par <- if (all(w_diag > 0)) 1 / mean(w_diag) else sigma2_est

        # Calculate the penalized precision matrix: X'WX + lambda*P

        XWX <- t(B) %*% W %*% B

        penalized_precision <- XWX + lambda * P

        # Variance-covariance matrix for this smooth: (X'WX + lambda*P)^(-1) * sigma^2

        tryCatch(
          {
            smooth_vcov <- solve(penalized_precision) * sigma2_par

            smooth_se <- .gl_sqrt_derived_variance(
              diag(smooth_vcov), "smooth coefficient covariance",
              allow_zero = FALSE
            )

            # Store results

            smooth_vcov_list[[par_name]][[s_name]] <- smooth_vcov

            smooth_se_list[[par_name]][[s_name]] <- smooth_se

            # Also calculate the smoother matrix for fitted values variance

            # A = X(X'WX + lambda*P)^(-1)X'W

            smoother_matrix <- B %*% solve(penalized_precision) %*% t(B) %*% W

            fitted_se <- .gl_sqrt_derived_variance(
              as.vector(diag(smoother_matrix)) * sigma2_par,
              "smooth fitted covariance",
              allow_zero = FALSE
            )

            cat(sprintf("\nSmooth term variance estimates for %s:%s\n", par_name, s_name))

            cat(sprintf(
              "  Basis coefficients SE: min=%.4f, max=%.4f, mean=%.4f\n",
              min(smooth_se), max(smooth_se), mean(smooth_se)
            ))

            cat(sprintf(
              "  Fitted values SE: min=%.4f, max=%.4f, mean=%.4f\n",
              min(fitted_se), max(fitted_se), mean(fitted_se)
            ))

            cat(sprintf(
              "  Effective DF: %.2f (trace of smoother matrix)\n",
              sum(diag(smoother_matrix))
            ))

            cat(sprintf("  Smoothing parameter lambda: %.4f\n", lambda))
          },
          error = function(e) {
            warning(sprintf(
              "Could not calculate variance for smooth %s:%s - %s",
              par_name, s_name, e$message
            ))

            smooth_vcov_list[[par_name]][[s_name]] <- NULL

            smooth_se_list[[par_name]][[s_name]] <- NULL
          }
        )
      }
    }
  }

  list(
    smooth_vcov = smooth_vcov_list,
    smooth_se = smooth_se_list
  )
}
