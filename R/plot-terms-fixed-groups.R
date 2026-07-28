#' Add a ggplot component while keeping compatibility with current ggplot2
#'
#' @noRd
.plot_fixed_terms_gg_add <- function(plot, object, object_name = "") {
  ggplot2::ggplot_add(object, plot, object_name)
}

#' Build factor main-effect groups for fixed-term plots
#'
#' @noRd
.plot_fixed_terms_factor_groups <- function(X, data) {
  groups <- list()

  if (is.null(X) || ncol(X) == 0) {
    return(groups)
  }

  x_cols <- colnames(X)

  assign <- attr(X, "assign")

  term_labels <- attr(X, "term.labels")

  if (!is.null(assign) && !is.null(term_labels) && length(assign) == length(x_cols)) {
    for (term_idx in seq_along(term_labels)) {
      term_name <- term_labels[term_idx]

      if (grepl(":", term_name, fixed = TRUE)) next

      term_cols <- x_cols[assign == term_idx]

      if (length(term_cols) == 0) next

      term_var <- .plot_fixed_terms_clean_expr_name(term_name)

      levs <- .plot_fixed_terms_levels_from_data(term_name, term_var, data)

      if (is.null(levs)) {
        levs <- .plot_fixed_terms_levels_from_factor_expr(term_name, term_cols)
      }

      if (is.null(levs) || length(levs) < 2) next

      fg <- .plot_fixed_terms_make_factor_group_from_cols(term_name, term_cols, levs)

      if (!is.null(fg)) {
        groups[[term_name]] <- fg
      }
    }
  }

  if (is.null(data) || !is.data.frame(data)) {
    return(groups)
  }

  for (var_name in names(data)) {
    if (var_name %in% names(groups)) next

    fg <- .plot_fixed_terms_fallback_factor_group(var_name, data, x_cols, groups)

    if (!is.null(fg)) {
      groups[[var_name]] <- fg
    }
  }

  groups
}

#' Build factor-by-factor interaction groups for fixed-term plots
#'
#' @noRd
.plot_fixed_terms_factor_interaction_groups <- function(X, factor_groups) {
  groups <- list()

  if (is.null(X) || ncol(X) == 0 || length(factor_groups) == 0) {
    return(groups)
  }

  x_cols <- colnames(X)

  fg_names <- names(factor_groups)

  if (length(fg_names) < 2) {
    return(groups)
  }

  for (i in seq_len(length(fg_names) - 1)) {
    for (j in (i + 1):length(fg_names)) {
      g1 <- factor_groups[[fg_names[i]]]

      g2 <- factor_groups[[fg_names[j]]]

      if (length(g1$level_col_map) == 0 || length(g2$level_col_map) == 0) next

      n1 <- length(g1$levels)

      n2 <- length(g2$levels)

      is_gender_1 <- grepl("gender|sex", g1$var_name, ignore.case = TRUE)

      is_gender_2 <- grepl("gender|sex", g2$var_name, ignore.case = TRUE)

      is_time_1 <- grepl("time", g1$var_name, ignore.case = TRUE)

      is_time_2 <- grepl("time", g2$var_name, ignore.case = TRUE)

      if (n1 > n2 || (n1 == n2 && is_time_1 && !is_time_2) || (n1 == n2 && !is_gender_1 && is_gender_2)) {
        panel_group <- g1

        other_group <- g2
      } else {
        panel_group <- g2

        other_group <- g1
      }

      panel_level_col_map <- panel_group$level_col_map

      other_level_col_map <- other_group$level_col_map

      interaction_col_map <- list()

      matched_cols <- character(0)

      for (panel_lev in names(panel_level_col_map)) {
        panel_col <- panel_level_col_map[[panel_lev]]

        interaction_col_map[[panel_lev]] <- list()

        for (other_lev in names(other_level_col_map)) {
          other_col <- other_level_col_map[[other_lev]]

          candidates <- c(
            paste0(other_col, ":", panel_col),
            paste0(panel_col, ":", other_col)
          )

          hit <- candidates[candidates %in% x_cols]

          if (length(hit) > 0) {
            interaction_col_map[[panel_lev]][[other_lev]] <- hit[1]

            matched_cols <- c(matched_cols, hit[1])
          }
        }
      }

      if (length(matched_cols) > 0) {
        interaction_name <- paste(other_group$var_name, panel_group$var_name, sep = ":")

        groups[[interaction_name]] <- list(
          interaction_name = interaction_name,
          panel_group = panel_group,
          other_group = other_group,
          interaction_col_map = interaction_col_map,
          matched_cols = unique(matched_cols)
        )
      }
    }
  }

  groups
}
