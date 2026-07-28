.plot_fixed_terms_factor_data <- function(fg, par_name, object, V, z) {
  levs <- fg$levels
  x_plot <- seq_along(levs)
  fitted_term <- rep(NA_real_, length(levs))
  term_se <- rep(NA_real_, length(levs))

  for (j in seq_along(levs)) {
    lev <- levs[j]
    if (identical(lev, fg$ref_level)) {
      fitted_term[j] <- 0
      term_se[j] <- 0
    } else if (lev %in% names(fg$level_col_map)) {
      coef_name_lev <- paste(par_name, fg$level_col_map[[lev]], sep = ".")
      coef_info <- .plot_fixed_terms_coef_info(coef_name_lev, object, V)
      fitted_term[j] <- coef_info$estimate
      term_se[j] <- coef_info$se
    }
  }

  keep <- is.finite(fitted_term) & is.finite(term_se)
  ci_lower <- fitted_term - z * term_se
  ci_upper <- fitted_term + z * term_se
  plot_df <- data.frame(
    x = x_plot,
    fitted = fitted_term,
    ci_lower = ci_lower,
    ci_upper = ci_upper,
    keep = keep
  )

  list(
    x_plot = x_plot,
    fitted = fitted_term,
    se = term_se,
    ci_lower = ci_lower,
    ci_upper = ci_upper,
    keep = keep,
    y_lim = .plot_fixed_terms_y_limits(fitted_term[keep], ci_lower[keep], ci_upper[keep]),
    plot_df = plot_df
  )
}

.plot_fixed_terms_factor_factor_data <- function(ig, par_name, object, V, z) {
  pg <- ig$panel_group
  og <- ig$other_group
  panel_levels <- pg$levels
  other_levels <- og$levels
  x_plot <- seq_along(panel_levels)
  plot_rows <- list()

  for (other_lev in other_levels) {
    fitted_term <- rep(NA_real_, length(panel_levels))
    term_se <- rep(NA_real_, length(panel_levels))

    for (j in seq_along(panel_levels)) {
      panel_lev <- panel_levels[j]
      if (identical(panel_lev, pg$ref_level) || identical(other_lev, og$ref_level)) {
        fitted_term[j] <- 0
        term_se[j] <- 0
      } else if (panel_lev %in% names(ig$interaction_col_map) && other_lev %in% names(ig$interaction_col_map[[panel_lev]])) {
        coef_name_lev <- paste(par_name, ig$interaction_col_map[[panel_lev]][[other_lev]], sep = ".")
        coef_info <- .plot_fixed_terms_coef_info(coef_name_lev, object, V)
        fitted_term[j] <- coef_info$estimate
        term_se[j] <- coef_info$se
      }
    }

    keep <- is.finite(fitted_term) & is.finite(term_se)
    ci_lower <- fitted_term - z * term_se
    ci_upper <- fitted_term + z * term_se
    plot_rows[[length(plot_rows) + 1]] <- data.frame(
      x = x_plot,
      group = factor(rep(other_lev, length(panel_levels)), levels = other_levels),
      fitted = fitted_term,
      se = term_se,
      ci_lower = ci_lower,
      ci_upper = ci_upper,
      keep = keep,
      stringsAsFactors = FALSE
    )
  }

  plot_df <- do.call(rbind, plot_rows)
  list(
    panel_levels = panel_levels,
    other_levels = other_levels,
    x_plot = x_plot,
    x_labels = panel_levels,
    plot_df = plot_df,
    y_lim = .plot_fixed_terms_y_limits(
      plot_df$fitted[plot_df$keep],
      plot_df$ci_lower[plot_df$keep],
      plot_df$ci_upper[plot_df$keep]
    )
  )
}
