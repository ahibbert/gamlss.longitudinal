.plot_terms_clean_expr_name <- function(x) {
  x <- trimws(x)
  x <- gsub("`", "", x, fixed = TRUE)
  factor_match <- regexec("^(?:as\\.)?factor\\(([^)]+)\\)$", x)
  matched <- regmatches(x, factor_match)[[1]]
  if (length(matched) >= 2) {
    x <- trimws(matched[2])
  }
  x
}

.plot_terms_is_factor_term <- function(term_name, data_for_terms = NULL) {
  term_var <- .plot_terms_clean_expr_name(term_name)

  if (grepl("^(?:as\\.)?factor\\(", term_name)) {
    return(TRUE)
  }

  if (is.null(data_for_terms) || !is.data.frame(data_for_terms)) {
    return(FALSE)
  }

  candidates <- unique(c(
    term_var,
    make.names(term_var),
    if (identical(term_var, "time_covariate")) {
      names(data_for_terms)[grepl("time", names(data_for_terms), ignore.case = TRUE)]
    } else {
      character(0)
    }
  ))

  any(candidates %in% names(data_for_terms) & vapply(candidates, function(candidate) {
    candidate %in% names(data_for_terms) && is.factor(data_for_terms[[candidate]])
  }, logical(1)))
}

.plot_terms_count <- function(obj,
                              data = NULL,
                              include_intercept = FALSE,
                              plot_interactions = FALSE) {
  n_smooth <- 0
  n_fixed <- 0
  data_for_terms <- data

  if ((is.null(data_for_terms) || !is.data.frame(data_for_terms)) && !is.null(obj$dataset)) {
    data_for_terms <- obj$dataset
  }

  for (par_name in names(obj$par_s)) {
    if (length(obj$par_s[[par_name]]) > 0) {
      n_smooth <- n_smooth + length(obj$par_s[[par_name]])
    }
  }

  for (par_name in names(obj$model_matrix$x)) {
    X <- obj$model_matrix$x[[par_name]]

    if (!is.null(X) && ncol(X) > 0) {
      x_cols <- colnames(X)
      assign <- attr(X, "assign")
      term_labels <- attr(X, "term.labels")
      grouped_cols <- character(0)

      if (!is.null(assign) && !is.null(term_labels) && length(assign) == length(x_cols)) {
        for (term_idx in seq_along(term_labels)) {
          term_name <- term_labels[term_idx]
          if (grepl(":", term_name, fixed = TRUE)) next

          term_cols <- x_cols[assign == term_idx]
          if (length(term_cols) == 0 || !.plot_terms_is_factor_term(term_name, data_for_terms)) next

          coef_names <- paste(par_name, term_cols, sep = ".")
          if (any(coef_names %in% names(obj$par))) {
            n_fixed <- n_fixed + 1
            grouped_cols <- c(grouped_cols, term_cols)
          }
        }
      }

      keep_cols <- !(x_cols == "intercept" & !include_intercept)
      keep_cols <- keep_cols & !(x_cols %in% grouped_cols)
      coef_names <- paste(par_name, x_cols[keep_cols], sep = ".")
      if (!plot_interactions) {
        coef_names <- coef_names[!grepl(":", coef_names, fixed = TRUE)]
      }
      n_fixed <- n_fixed + sum(coef_names %in% names(obj$par))
    }
  }

  list(smooth = n_smooth, fixed = n_fixed, total = n_smooth + n_fixed)
}

.plot_terms_collect_plot_objects <- function(smooth_results = list(), fixed_results = list()) {
  plot_objects <- list()

  if (!is.null(smooth_results$plots)) {
    plot_objects <- c(plot_objects, smooth_results$plots)
  }
  if (!is.null(fixed_results$plots)) {
    plot_objects <- c(plot_objects, fixed_results$plots)
  }

  plot_objects
}

.plot_terms_dashboard_layout <- function(plot_objects, ncol = 4, paginate = FALSE) {
  if (length(plot_objects) == 0) {
    return(NULL)
  }

  if (isTRUE(paginate)) {
    return(list(plotlist = plot_objects, paginate = TRUE))
  }

  if (is.null(ncol)) {
    ncol <- min(2, length(plot_objects))
  }

  list(
    plotlist = plot_objects,
    ncol = ncol,
    nrow = ceiling(length(plot_objects) / ncol),
    paginate = FALSE
  )
}

.plot_terms_render_dashboard <- function(plot_objects, ncol = 4, paginate = FALSE) {
  if (length(plot_objects) == 0) {
    return(NULL)
  }

  if (length(plot_objects) > 16 && !isTRUE(paginate)) {
    warning(
      "More than 16 charts (", length(plot_objects), ") were generated. ",
      "Rendering all charts at once may fail in some environments. ",
      "Use paginate=TRUE to view one chart at a time.",
      call. = FALSE
    )
  }

  dashboard <- .plot_terms_dashboard_layout(plot_objects, ncol = ncol, paginate = paginate)

  if (isTRUE(paginate)) {
    for (i_plot in seq_along(plot_objects)) {
      grid::grid.newpage()
      print(plot_objects[[i_plot]])
      if (i_plot < length(plot_objects) && interactive()) {
        invisible(readline(prompt = sprintf("Press [Enter] for next chart (%d/%d)... ", i_plot, length(plot_objects))))
      }
    }
  } else {
    grid::grid.newpage()
    grid::pushViewport(grid::viewport(layout = grid::grid.layout(dashboard$nrow, dashboard$ncol)))
    for (i_plot in seq_along(plot_objects)) {
      r <- ((i_plot - 1) %/% dashboard$ncol) + 1
      c <- ((i_plot - 1) %% dashboard$ncol) + 1
      print(plot_objects[[i_plot]], vp = grid::viewport(layout.pos.row = r, layout.pos.col = c))
    }
    grid::popViewport()
  }

  dashboard
}
