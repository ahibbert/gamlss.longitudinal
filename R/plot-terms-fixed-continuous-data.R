.plot_fixed_terms_interaction_factor_data <- function(spec, par_name, object, V, z) {
  col_name <- spec$col_name
  coef_name <- spec$coef_name
  X <- object$model_matrix$x[[par_name]]
  x_raw <- as.numeric(X[, col_name])
  coef_info <- .plot_fixed_terms_coef_info(coef_name, object, V)
  x_levels <- sort(unique(x_raw[is.finite(x_raw)]))

  if (length(x_levels) == 0L) {
    return(NULL)
  }

  fitted_term <- x_levels * coef_info$estimate
  term_se <- abs(x_levels) * .gl_sqrt_derived_variance(
    coef_info$variance, "fixed-term covariance", allow_zero = FALSE
  )
  ci_lower <- fitted_term - z * term_se
  ci_upper <- fitted_term + z * term_se
  x_labels <- as.character(signif(x_levels, 6))
  x_plot <- seq_along(x_levels)

  list(
    col_name = col_name,
    coef_name = coef_name,
    x_levels = x_levels,
    x_labels = x_labels,
    x_plot = x_plot,
    fitted = fitted_term,
    se = term_se,
    ci_lower = ci_lower,
    ci_upper = ci_upper,
    y_lim = .plot_fixed_terms_y_limits(fitted_term, ci_lower, ci_upper),
    plot_df = data.frame(
      x = x_plot,
      fitted = fitted_term,
      ci_lower = ci_lower,
      ci_upper = ci_upper,
      stringsAsFactors = FALSE
    )
  )
}

.plot_fixed_terms_continuous_data <- function(spec,
                                              par_name,
                                              object,
                                              V,
                                              z,
                                              sort_x = TRUE,
                                              fallback_to_index = TRUE) {
  col_name <- spec$col_name
  coef_name <- spec$coef_name
  X <- object$model_matrix$x[[par_name]]
  x_raw <- as.numeric(X[, col_name])
  coef_info <- .plot_fixed_terms_coef_info(coef_name, object, V)

  fitted_term <- x_raw * coef_info$estimate
  term_se <- .gl_sqrt_derived_variance(
    (x_raw^2) * coef_info$variance,
    "fixed-term prediction covariance",
    allow_zero = TRUE
  )
  ci_lower <- fitted_term - z * term_se
  ci_upper <- fitted_term + z * term_se

  if (length(unique(x_raw)) <= 1 && fallback_to_index) {
    x_plot <- seq_along(x_raw)
    xlab_text <- paste(col_name, "(index)")
    ord <- seq_along(x_plot)
  } else {
    x_plot <- x_raw
    xlab_text <- col_name
    ord <- if (sort_x) order(x_plot) else seq_along(x_plot)
  }

  list(
    col_name = col_name,
    coef_name = coef_name,
    x_plot = x_plot,
    xlab_text = xlab_text,
    ord = ord,
    fitted = fitted_term,
    se = term_se,
    ci_lower = ci_lower,
    ci_upper = ci_upper,
    y_lim = .plot_fixed_terms_y_limits(fitted_term, ci_lower, ci_upper),
    plot_df = data.frame(
      x = x_plot[ord],
      fitted = fitted_term[ord],
      ci_lower = ci_lower[ord],
      ci_upper = ci_upper[ord]
    )
  )
}
