.plot_fixed_terms_y_limits <- function(...) {
  y_vals <- c(...)
  y_vals <- y_vals[is.finite(y_vals)]
  if (length(y_vals) == 0L) {
    return(NULL)
  }

  y_rng <- range(y_vals)
  y_pad <- 0.05 * max(1e-8, diff(y_rng))
  c(y_rng[1] - y_pad, y_rng[2] + y_pad)
}

.plot_fixed_terms_dashboard <- function(plot_objects, ncol = NULL) {
  if (length(plot_objects) == 0L) {
    return(NULL)
  }
  if (is.null(ncol)) {
    ncol <- min(2, length(plot_objects))
  }

  list(
    plotlist = plot_objects,
    ncol = ncol,
    nrow = ceiling(length(plot_objects) / ncol)
  )
}

.plot_fixed_terms_has_coef <- function(coef_name, object, V) {
  coef_name %in% names(object$par) &&
    coef_name %in% rownames(V) &&
    coef_name %in% colnames(V)
}

.plot_fixed_terms_coef_info <- function(coef_name, object, V) {
  if (!.plot_fixed_terms_has_coef(coef_name, object, V)) {
    return(list(estimate = NA_real_, variance = NA_real_, se = NA_real_))
  }

  variance <- as.numeric(V[coef_name, coef_name])
  list(
    estimate = as.numeric(object$par[coef_name]),
    variance = variance,
    se = sqrt(pmax(0, variance))
  )
}

.plot_fixed_terms_ci_caption <- function(prefix, ci_level) {
  paste(prefix, "/", round(ci_level * 100), "% CI")
}
