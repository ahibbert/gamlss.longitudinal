.plot_fixed_terms_plot_specs <- function(object,
                                         V,
                                         data_for_terms = NULL,
                                         include_intercept = FALSE,
                                         plot_interactions = FALSE) {
  plot_specs <- list()

  for (par_name in names(object$model_matrix$x)) {
    X <- object$model_matrix$x[[par_name]]
    if (is.null(X) || ncol(X) == 0) next

    factor_groups <- .plot_fixed_terms_factor_groups(X, data_for_terms)
    interaction_groups <- if (plot_interactions) {
      .plot_fixed_terms_factor_interaction_groups(X, factor_groups)
    } else {
      list()
    }

    grouped_cols <- unique(unlist(lapply(factor_groups, function(g) g$matched_cols), use.names = FALSE))
    if (length(grouped_cols) == 0) grouped_cols <- character(0)

    grouped_interaction_cols <- unique(unlist(lapply(interaction_groups, function(g) g$matched_cols), use.names = FALSE))
    if (length(grouped_interaction_cols) == 0) grouped_interaction_cols <- character(0)

    for (var_name in names(factor_groups)) {
      fg <- factor_groups[[var_name]]
      has_valid_coef <- FALSE

      for (lev in names(fg$level_col_map)) {
        coef_name <- paste(par_name, fg$level_col_map[[lev]], sep = ".")
        if (.plot_fixed_terms_has_coef(coef_name, object, V)) {
          has_valid_coef <- TRUE
          break
        }
      }

      if (has_valid_coef) {
        plot_specs[[length(plot_specs) + 1]] <- list(
          type = "factor",
          par_name = par_name,
          var_name = var_name,
          group = fg
        )
      }
    }

    for (inter_name in names(interaction_groups)) {
      ig <- interaction_groups[[inter_name]]
      has_valid_coef <- FALSE

      for (panel_lev in names(ig$interaction_col_map)) {
        for (other_lev in names(ig$interaction_col_map[[panel_lev]])) {
          coef_name <- paste(par_name, ig$interaction_col_map[[panel_lev]][[other_lev]], sep = ".")
          if (.plot_fixed_terms_has_coef(coef_name, object, V)) {
            has_valid_coef <- TRUE
            break
          }
        }
        if (has_valid_coef) break
      }

      if (has_valid_coef) {
        plot_specs[[length(plot_specs) + 1]] <- list(
          type = "interaction_factor_factor",
          par_name = par_name,
          group = ig
        )
      }
    }

    for (col_name in colnames(X)) {
      if (col_name %in% grouped_cols) next
      if (col_name %in% grouped_interaction_cols) next
      if (!include_intercept && col_name == "intercept") next
      if (!plot_interactions && grepl(":", col_name, fixed = TRUE)) next

      coef_name <- paste(par_name, col_name, sep = ".")
      if (!.plot_fixed_terms_has_coef(coef_name, object, V)) next

      term_type <- if (grepl(":", col_name, fixed = TRUE)) "interaction_factor" else "continuous"
      plot_specs[[length(plot_specs) + 1]] <- list(
        type = term_type,
        par_name = par_name,
        col_name = col_name,
        coef_name = coef_name
      )
    }
  }

  plot_specs
}
