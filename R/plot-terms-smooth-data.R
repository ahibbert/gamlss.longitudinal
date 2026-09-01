.plot_smooth_terms_index <- function(object) {
  smooth_index <- list()

  for (par_name in names(object$par_s)) {
    if (length(object$par_s[[par_name]]) == 0) next

    for (s_name in names(object$par_s[[par_name]])) {
      B <- object$model_matrix$s[[par_name]][[s_name]]
      beta_s <- object$par_s[[par_name]][[s_name]]
      if (is.null(B) || is.null(beta_s)) next

      smooth_index[[length(smooth_index) + 1]] <- list(par_name = par_name, s_name = s_name)
    }
  }

  smooth_index
}

.plot_smooth_terms_fit_se <- function(B, smooth_vcov = NULL, smooth_se = NULL) {
  if (!is.null(smooth_vcov) && all(dim(smooth_vcov) == c(ncol(B), ncol(B)))) {
    return(.gl_sqrt_derived_variance(
      diag(B %*% smooth_vcov %*% t(B)),
      "smooth-term covariance", allow_zero = TRUE
    ))
  }

  if (!is.null(smooth_se) && length(smooth_se) == ncol(B)) {
    beta_var_diag <- as.numeric(smooth_se)^2
    return(.gl_sqrt_derived_variance(
      rowSums((B^2) * rep(beta_var_diag, each = nrow(B))),
      "smooth-term diagonal covariance", allow_zero = TRUE
    ))
  }

  rep(NA_real_, nrow(B))
}

.plot_smooth_terms_plot_df <- function(x,
                                       fitted_smooth,
                                       ci_lower,
                                       ci_upper,
                                       sort_x = TRUE,
                                       even_grid = TRUE,
                                       grid_n = 200) {
  if (isTRUE(even_grid)) {
    x_ok <- is.finite(x)
    df_obs <- data.frame(
      x = x[x_ok],
      fitted = fitted_smooth[x_ok],
      ci_lower = ci_lower[x_ok],
      ci_upper = ci_upper[x_ok]
    )

    if (nrow(df_obs) >= 2 && length(unique(df_obs$x)) >= 2) {
      agg_df <- stats::aggregate(
        df_obs[, c("fitted", "ci_lower", "ci_upper")],
        by = list(x = df_obs$x),
        FUN = mean
      )
      agg_df <- agg_df[order(agg_df$x), , drop = FALSE]
      n_grid_use <- max(20, as.integer(grid_n))
      x_grid <- seq(min(agg_df$x), max(agg_df$x), length.out = n_grid_use)

      safe_approx <- function(y) {
        ok <- is.finite(agg_df$x) & is.finite(y)
        if (sum(ok) >= 2 && length(unique(agg_df$x[ok])) >= 2) {
          stats::approx(agg_df$x[ok], y[ok], xout = x_grid, method = "linear", rule = 2)$y
        } else if (sum(ok) == 1) {
          rep(y[ok][1], length(x_grid))
        } else {
          rep(NA_real_, length(x_grid))
        }
      }

      return(data.frame(
        x = x_grid,
        fitted = safe_approx(agg_df$fitted),
        ci_lower = safe_approx(agg_df$ci_lower),
        ci_upper = safe_approx(agg_df$ci_upper)
      ))
    }
  }

  ord <- if (sort_x) order(x) else seq_along(x)
  data.frame(
    x = x[ord],
    fitted = fitted_smooth[ord],
    ci_lower = ci_lower[ord],
    ci_upper = ci_upper[ord]
  )
}
