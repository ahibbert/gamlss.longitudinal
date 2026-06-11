#' Plot all smooth terms with confidence bands
#'
#' This utility plots every smooth term in a fitted `gamlss.longitudinal` object
#' and computes pointwise confidence bands using the smooth coefficient
#' covariance matrices returned by `vcov.gamlss.longitudinal()`.
#'
#' @param object A fitted `gamlss.longitudinal` object.
#' @param vcov_obj Optional output from `vcov(object, ...)`. If `NULL`, this is
#' computed internally with the analytical vcov path.
#' @param data Optional data frame containing original covariates used for the
#' x-axis variable of each smooth.
#' @param ci_level Confidence level for pointwise intervals.
#' @param ncol Number of plot columns (defaults to 2 or fewer if needed).
#' @param ci_col Color for confidence bands.
#' @param fit_col Color for fitted smooth line.
#' @param ci_lty Line type for confidence bands.
#' @param fit_lwd Line width for fitted smooth.
#' @param sort_x Logical; sort points by x before plotting lines.
#' @param even_grid Logical; if TRUE, plot smooths on an evenly spaced x-grid
#' built over observed x-range.
#' @param grid_n Number of grid points when `even_grid = TRUE`.
#' @param fallback_to_index Logical; if x variable cannot be inferred, plot
#' against row index.
#' @param setup_mfrow Logical; if TRUE (default), configure par(mfrow) inside
#' this function. Set FALSE when caller configures layout.
#' @param show_legend Logical; if TRUE, draw a small legend in each panel.
#'
#' @return Invisibly returns a nested list with x, fitted values, standard
#' errors, and confidence limits for each smooth term.
#' @export
plot_smooth_terms = function(

  object,

  vcov_obj = NULL,

  data = NULL,

  ci_level = 0.95,

  ncol = NULL,

  ci_col = "red",

  fit_col = "black",

  ci_lty = 2,

  fit_lwd = 2,

  sort_x = TRUE,

  even_grid = TRUE,

  grid_n = 200,

  fallback_to_index = TRUE,

  setup_mfrow = TRUE,

  show_legend = TRUE

) {

  if(!inherits(object, "gamlss.longitudinal")) {

    stop("'object' must be of class 'gamlss.longitudinal'.")

  }


  if(is.null(vcov_obj)) {

    vcov_obj = .resolve_vcov(object, numderiv = FALSE, extra_args = list(method = "analytical"))

  }


  if(!is.list(vcov_obj) || is.null(vcov_obj$vcov)) {

    stop("'vcov_obj' must be the output of vcov.gamlss.longitudinal().")

  }


  smooth_vcov_list = vcov_obj$vcov$smooth_vcov

  smooth_se_list = vcov_obj$vcov$smooth_se


  extract_smooth_var = function(s_name) {

    s_txt = trimws(s_name)

    s_call = tryCatch(parse(text = s_txt)[[1]], error = function(e) NULL)

    if(!is.null(s_call) && length(s_call) >= 2) {

      out = paste(deparse(s_call[[2]]), collapse = " ")

    } else {

      out = sub("^s\\((.*)\\)$", "\\1", s_txt)

    }

    out = trimws(gsub("`", "", out, fixed = TRUE))

    out

  }


  eval_smooth_x = function(x_expr, data_frame) {

    if(is.null(data_frame) || !is.data.frame(data_frame)) {

      return(NULL)

    }

    tryCatch(eval(parse(text = x_expr)[[1]], envir = data_frame), error = function(e) NULL)

  }


  get_x_for_smooth = function(par_name, s_name, B) {

    x_var = extract_smooth_var(s_name)

    x = NULL


    # Prefer x saved with the smooth basis because it is guaranteed row-aligned.

    x_basis = attr(B, "smooth_x")

    x_basis_var = attr(B, "smooth_var")

    if(!is.null(x_basis) && length(x_basis) == nrow(B)) {

      x = x_basis

      if(!is.null(x_basis_var) && nzchar(x_basis_var)) {

        x_var = x_basis_var

      }

    }


    if(is.null(x) && !is.null(data) && is.data.frame(data)) {

      data_names = names(data)

      # Prefer exact match, then case-insensitive match, then make.names match.

      idx_exact = which(data_names == x_var)

      idx_ci = which(tolower(data_names) == tolower(x_var))

      idx_mn = which(make.names(data_names) == make.names(x_var))

      idx = c(idx_exact, idx_ci, idx_mn)

      idx = idx[!duplicated(idx)]

      if(length(idx) > 0) {

        matched_name = data_names[idx[1]]

        x_candidate = data[[matched_name]]

        if(length(x_candidate) == nrow(B)) {

          x = x_candidate

        } else if(!is.null(rownames(B)) && !is.null(rownames(data))) {

          row_idx = match(rownames(B), rownames(data))

          if(all(!is.na(row_idx))) {

            x = x_candidate[row_idx]

          }

        }

        x_var = matched_name

      }


      if(is.null(x)) {

        x_candidate = eval_smooth_x(x_var, data)

        if(!is.null(x_candidate) && length(x_candidate) == nrow(B)) {

          x = x_candidate

        }

      }

    }


    if(is.null(x) && !is.null(object$model_matrix$x[[par_name]]) && x_var %in% colnames(object$model_matrix$x[[par_name]])) {

      x = object$model_matrix$x[[par_name]][, x_var]

    }


    if(is.null(x) && fallback_to_index) {

      x = seq_len(nrow(B))

      x_var = "index"

      warning("Falling back to index for smooth term '", s_name,

              "' because covariate was not found in supplied data/model matrix.")

    }


    if(is.null(x)) {

      stop("Could not infer x-axis for smooth term '", s_name, "'. Provide 'data' with the smooth covariate columns.")

    }


    if(length(x) != nrow(B)) {

      stop("Length mismatch for smooth term '", s_name, "': length(x)=", length(x), " but nrow(B)=", nrow(B), ".")

    }


    list(x = as.numeric(x), x_var = x_var)

  }


  z = qnorm((1 + ci_level) / 2)

  smooth_index = list()


  for(par_name in names(object$par_s)) {

    if(length(object$par_s[[par_name]]) == 0) next

    for(s_name in names(object$par_s[[par_name]])) {

      B = object$model_matrix$s[[par_name]][[s_name]]

      beta_s = object$par_s[[par_name]][[s_name]]

      if(is.null(B) || is.null(beta_s)) next


      smooth_index[[length(smooth_index) + 1]] = list(par_name = par_name, s_name = s_name)

    }

  }


  n_plots = length(smooth_index)

  if(n_plots == 0) {

    warning("No smooth terms found to plot.")

    return(invisible(list()))

  }


  out = list()

  plot_objects = list()

  for(i in seq_len(n_plots)) {

    par_name = smooth_index[[i]]$par_name

    s_name = smooth_index[[i]]$s_name


    B = object$model_matrix$s[[par_name]][[s_name]]

    beta_s = object$par_s[[par_name]][[s_name]]

    x_info = get_x_for_smooth(par_name, s_name, B)

    x = x_info$x


    fitted_smooth = as.numeric(B %*% beta_s)


    smooth_vcov = NULL

    smooth_se = NULL

    if(!is.null(smooth_vcov_list) && !is.null(smooth_vcov_list[[par_name]])) {

      smooth_vcov = smooth_vcov_list[[par_name]][[s_name]]

    }

    if(!is.null(smooth_se_list) && !is.null(smooth_se_list[[par_name]])) {

      smooth_se = smooth_se_list[[par_name]][[s_name]]

    }


    if(!is.null(smooth_vcov) && all(dim(smooth_vcov) == c(ncol(B), ncol(B)))) {

      smooth_fit_se = sqrt(pmax(0, diag(B %*% smooth_vcov %*% t(B))))

    } else if(!is.null(smooth_se) && length(smooth_se) == ncol(B)) {

      beta_var_diag = as.numeric(smooth_se)^2

      smooth_fit_se = sqrt(pmax(0, rowSums((B^2) * rep(beta_var_diag, each = nrow(B)))))

    } else {

      smooth_fit_se = rep(NA_real_, nrow(B))

    }


    ci_lower = fitted_smooth - z * smooth_fit_se

    ci_upper = fitted_smooth + z * smooth_fit_se

    main_title = paste(par_name, s_name, sep = ": ")

    ylab_text = paste("smooth(", x_info$x_var, ")", sep = "")


    if(isTRUE(even_grid)) {

      x_ok = is.finite(x)

      df_obs = data.frame(

        x = x[x_ok],

        fitted = fitted_smooth[x_ok],

        ci_lower = ci_lower[x_ok],

        ci_upper = ci_upper[x_ok]

      )


      if(nrow(df_obs) >= 2 && length(unique(df_obs$x)) >= 2) {

        agg_df = stats::aggregate(df_obs[, c("fitted", "ci_lower", "ci_upper")], by = list(x = df_obs$x), FUN = mean)

        agg_df = agg_df[order(agg_df$x), , drop = FALSE]

        n_grid_use = max(20, as.integer(grid_n))

        x_grid = seq(min(agg_df$x), max(agg_df$x), length.out = n_grid_use)

        safe_approx = function(y) {

          ok = is.finite(agg_df$x) & is.finite(y)

          if(sum(ok) >= 2 && length(unique(agg_df$x[ok])) >= 2) {

            stats::approx(agg_df$x[ok], y[ok], xout = x_grid, method = "linear", rule = 2)$y

          } else if(sum(ok) == 1) {

            rep(y[ok][1], length(x_grid))

          } else {

            rep(NA_real_, length(x_grid))

          }

        }


        plot_df = data.frame(

          x = x_grid,

          fitted = safe_approx(agg_df$fitted),

          ci_lower = safe_approx(agg_df$ci_lower),

          ci_upper = safe_approx(agg_df$ci_upper)

        )

      } else {

        ord = if(sort_x) order(x) else seq_along(x)

        plot_df = data.frame(

          x = x[ord],

          fitted = fitted_smooth[ord],

          ci_lower = ci_lower[ord],

          ci_upper = ci_upper[ord]

        )

      }

    } else {

      ord = if(sort_x) order(x) else seq_along(x)

      plot_df = data.frame(

        x = x[ord],

        fitted = fitted_smooth[ord],

        ci_lower = ci_lower[ord],

        ci_upper = ci_upper[ord]

      )

    }


    y_vals = c(plot_df$fitted, plot_df$ci_lower, plot_df$ci_upper)

    y_vals = y_vals[is.finite(y_vals)]

    y_lim = NULL

    if(length(y_vals) > 0) {

      y_rng = range(y_vals)

      y_pad = 0.05 * max(1e-8, diff(y_rng))

      y_lim = c(y_rng[1] - y_pad, y_rng[2] + y_pad)

    }


    p = ggplot2::ggplot(plot_df, ggplot2::aes(x = x, y = fitted)) +

      ggplot2::geom_ribbon(ggplot2::aes(ymin = ci_lower, ymax = ci_upper), fill = ci_col, alpha = 0.16) +

      ggplot2::geom_line(color = fit_col, linewidth = fit_lwd) +

      ggplot2::labs(title = main_title, x = x_info$x_var, y = ylab_text)


    if(!is.null(y_lim)) {

      p = p + ggplot2::coord_cartesian(ylim = y_lim)

    }


    if(show_legend) {

      p = p + ggplot2::labs(caption = paste("fit /", round(ci_level * 100), "% CI"))

    }


    p = p + ggplot2::theme_minimal()


    plot_objects[[length(plot_objects) + 1]] = p


    if(is.null(out[[par_name]])) out[[par_name]] = list()

    out[[par_name]][[s_name]] = list(

      x = x,

      fitted = fitted_smooth,

      se = smooth_fit_se,

      ci_lower = ci_lower,

      ci_upper = ci_upper,

      plot = p

    )

  }


  if(length(plot_objects) > 0) {

    if(is.null(ncol)) {

      ncol = min(2, n_plots)

    }

    nrow = ceiling(length(plot_objects) / ncol)

    dashboard = list(plotlist = plot_objects, ncol = ncol, nrow = nrow)

    if(setup_mfrow) {

      grid::grid.newpage()

      grid::pushViewport(grid::viewport(layout = grid::grid.layout(nrow, ncol)))

      for(i_plot in seq_along(plot_objects)) {

        r = ((i_plot - 1) %/% ncol) + 1

        c = ((i_plot - 1) %% ncol) + 1

        print(plot_objects[[i_plot]], vp = grid::viewport(layout.pos.row = r, layout.pos.col = c))

      }

      grid::popViewport()

    }

    out$plots = plot_objects

    out$dashboard = dashboard

  }


  invisible(out)

}


#' Plot all fixed terms with confidence bands

#'

#' This utility plots fixed-effect term contributions for a fitted

#' `gamlss.longitudinal` object using coefficient uncertainty from

#' `vcov.gamlss.longitudinal()`.

#'

#' For each fixed-effect design-matrix column \eqn{x_j}, it plots

#' \eqn{x_j \hat{\beta}_j} with pointwise confidence bands

#' \eqn{x_j \hat{\beta}_j \pm z_{\alpha/2}\sqrt{x_j^2 \mathrm{Var}(\hat{\beta}_j)}}.

#'

#' @param object A fitted `gamlss.longitudinal` object.

#' @param vcov_obj Optional output from `vcov(object, ...)`. If `NULL`, this is

#' computed internally with the analytical vcov path.

#' @param ci_level Confidence level for pointwise intervals.

#' @param ncol Number of plot columns (defaults to 2 or fewer if needed).

#' @param include_intercept Logical; include intercept columns in plots.

#' @param plot_interactions Logical; include interaction columns in plots.

#' @param ci_col Color for confidence bands.

#' @param fit_col Color for fitted fixed-term line.

#' @param ci_lty Line type for confidence bands.

#' @param fit_lwd Line width for fitted fixed-term line.

#' @param sort_x Logical; sort x-values before drawing lines.

#' @param fallback_to_index Logical; if x has one unique value, use index on x-axis.

#' @param setup_mfrow Logical; if TRUE (default), configure par(mfrow) inside

#' this function. Set FALSE when caller configures layout.

#' @param data Optional data frame used to detect factor columns and show

#' factor levels on x-axis for categorical fixed terms. Factor terms are grouped

#' into one point-and-interval plot per model term and parameter.

#' @param factor_pch Point symbol for factor-level estimates.

#' @param factor_cex Point size for factor-level estimates.

#' @param show_legend Logical; if TRUE, draw a small legend in each panel.

#'

#' @return Invisibly returns a nested list with x, fitted values, standard

#' errors, and confidence limits for each fixed term.

#' @export

plot_fixed_terms = function(

  object,

  vcov_obj = NULL,

  ci_level = 0.95,

  ncol = NULL,

  include_intercept = FALSE,

  plot_interactions = FALSE,

  ci_col = "red",

  fit_col = "black",

  ci_lty = 2,

  fit_lwd = 2,

  sort_x = TRUE,

  fallback_to_index = TRUE,

  setup_mfrow = TRUE,

  data = NULL,

  factor_pch = 16,

  factor_cex = 1.2,

  show_legend = TRUE

) {

  if(!inherits(object, "gamlss.longitudinal")) {

    stop("'object' must be of class 'gamlss.longitudinal'.")

  }


  if(is.null(vcov_obj)) {

    vcov_obj = .resolve_vcov(object, numderiv = FALSE, extra_args = list(method = "analytical"))

  }


  if(!is.list(vcov_obj) || is.null(vcov_obj$vcov) || is.null(vcov_obj$vcov$overall)) {

    stop("'vcov_obj' must be the output of vcov.gamlss.longitudinal() with vcov$overall present.")

  }


  V = vcov_obj$vcov$overall

  if(is.null(rownames(V)) || is.null(colnames(V))) {

    stop("vcov$overall must have row and column names matching fixed coefficients.")

  }


  z = qnorm((1 + ci_level) / 2)

  gg_add = function(plot, object, object_name = "") {

    ggplot2::ggplot_add(object, plot, object_name)

  }

  data_for_terms = data

  if((is.null(data_for_terms) || !is.data.frame(data_for_terms)) && !is.null(object$dataset)) {

    data_for_terms = object$dataset

  }


  build_factor_groups = function(X, data) {

    groups = list()

    if(is.null(X) || ncol(X) == 0) return(groups)


    x_cols = colnames(X)

    assign = attr(X, "assign")

    term_labels = attr(X, "term.labels")


    clean_expr_name = function(x) {

      x = trimws(x)

      x = gsub("`", "", x, fixed = TRUE)

      factor_match = regexec("^(?:as\\.)?factor\\(([^)]+)\\)$", x)

      matched = regmatches(x, factor_match)[[1]]

      if(length(matched) >= 2) {

        x = trimws(matched[2])

      }

      x

    }


    make_factor_group_from_cols = function(term_name, term_cols, levs) {

      if(length(levs) < 2 || length(term_cols) == 0) return(NULL)

      level_col_map = list()

      matched_cols = character(0)


      for(lev in levs[-1]) {

        lev_plain = as.character(lev)

        lev_mn = make.names(lev_plain)

        hits = term_cols[

          endsWith(term_cols, lev_plain) |

            endsWith(term_cols, lev_mn) |

            grepl(paste0("(^|[^[:alnum:]_.])", gsub("([][{}()+*^$|\\\\?.])", "\\\\\\1", lev_plain), "$"), term_cols)

        ]

        if(length(hits) == 0 && length(term_cols) == length(levs) - 1L) {

          hits = term_cols[seq_along(levs[-1]) == which(levs[-1] == lev)]

        }

        if(length(hits) > 0) {

          level_col_map[[lev]] = hits[1]

          matched_cols = c(matched_cols, hits[1])

        }

      }


      if(length(level_col_map) == 0) return(NULL)

      list(

        var_name = term_name,

        levels = levs,

        ref_level = levs[1],

        level_col_map = level_col_map,

        matched_cols = unique(matched_cols)

      )

    }


    if(!is.null(assign) && !is.null(term_labels) && length(assign) == length(x_cols)) {

      for(term_idx in seq_along(term_labels)) {

        term_name = term_labels[term_idx]

        if(grepl(":", term_name, fixed = TRUE)) next

        term_cols = x_cols[assign == term_idx]

        if(length(term_cols) == 0) next


        term_var = clean_expr_name(term_name)

        levs = NULL

        if(!is.null(data) && is.data.frame(data)) {

          data_candidates = unique(c(

            term_var,

            make.names(term_var),

            if(identical(term_var, "time_covariate")) names(data)[grepl("time", names(data), ignore.case = TRUE)] else character(0)

          ))

          for(candidate in data_candidates) {

            if(candidate %in% names(data) && is.factor(data[[candidate]])) {

              levs = levels(data[[candidate]])

              break

            } else if(candidate %in% names(data) && grepl("^(?:as\\.)?factor\\(", term_name)) {

              levs = levels(as.factor(data[[candidate]]))

              break

            }

          }

        }

        if(is.null(levs) && grepl("^(?:as\\.)?factor\\(", term_name)) {

          levs = sub(paste0("^", gsub("([][{}()+*^$|\\\\?.])", "\\\\\\1", term_name)), "", term_cols)

          levs = c(sub("\\).*", ")", term_name), levs)

        }

        if(is.null(levs) || length(levs) < 2) next


        fg = make_factor_group_from_cols(term_name, term_cols, levs)

        if(!is.null(fg)) {

          groups[[term_name]] = fg

        }

      }

    }


    if(is.null(data) || !is.data.frame(data)) return(groups)


    for(var_name in names(data)) {

      if(var_name %in% names(groups)) next

      v = data[[var_name]]

      if(!is.factor(v)) next

      levs = levels(v)

      if(length(levs) < 2) next


      var_tokens = strsplit(var_name, "_", fixed = TRUE)[[1]]

      var_prefixes = unique(c(

        var_name,

        make.names(var_name),

        if(length(var_tokens) > 0) var_tokens[1] else character(0),

        if(length(var_tokens) > 0) make.names(var_tokens[1]) else character(0)

      ))

      # Internal fitting often renames the user time variable to time_covariate.

      if(grepl("time", var_name, ignore.case = TRUE) && any(grepl("^time_covariate", x_cols))) {

        var_prefixes = unique(c(var_prefixes, "time_covariate", make.names("time_covariate")))

      }

      var_prefixes = var_prefixes[nzchar(var_prefixes)]


      level_col_map = list()

      matched_cols = character(0)

      for(lev in levs[-1]) {

        lev_plain = as.character(lev)

        lev_mn = make.names(lev_plain)

        candidates = unique(unlist(lapply(var_prefixes, function(pref) {

          c(

            paste0(pref, lev_plain),

            paste0(pref, lev_mn),

            paste0(pref, "_", lev_plain),

            paste0(pref, "_", lev_mn)

          )

        }), use.names = FALSE))

        hit = candidates[candidates %in% x_cols]

        if(length(hit) > 0) {

          level_col_map[[lev]] = hit[1]

          matched_cols = c(matched_cols, hit[1])

        }

      }


      if(length(matched_cols) == 0 && length(levs) == 2 && var_name %in% x_cols) {

        level_col_map[[levs[2]]] = var_name

        matched_cols = var_name

      }


      if(length(level_col_map) > 0 && !any(matched_cols %in% unlist(lapply(groups, function(g) g$matched_cols), use.names = FALSE))) {

        groups[[var_name]] = list(

          var_name = var_name,

          levels = levs,

          ref_level = levs[1],

          level_col_map = level_col_map,

          matched_cols = unique(matched_cols)

        )

      }

    }

    groups

  }


  build_factor_interaction_groups = function(X, factor_groups) {

    groups = list()

    if(is.null(X) || ncol(X) == 0 || length(factor_groups) == 0) return(groups)


    x_cols = colnames(X)

    fg_names = names(factor_groups)

    if(length(fg_names) < 2) return(groups)


    for(i in seq_len(length(fg_names) - 1)) {

      for(j in (i + 1):length(fg_names)) {

        g1 = factor_groups[[fg_names[i]]]

        g2 = factor_groups[[fg_names[j]]]


        if(length(g1$level_col_map) == 0 || length(g2$level_col_map) == 0) next


        n1 = length(g1$levels)

        n2 = length(g2$levels)

        is_gender_1 = grepl("gender|sex", g1$var_name, ignore.case = TRUE)

        is_gender_2 = grepl("gender|sex", g2$var_name, ignore.case = TRUE)

        is_time_1 = grepl("time", g1$var_name, ignore.case = TRUE)

        is_time_2 = grepl("time", g2$var_name, ignore.case = TRUE)


        if(n1 > n2 || (n1 == n2 && is_time_1 && !is_time_2) || (n1 == n2 && !is_gender_1 && is_gender_2)) {

          panel_group = g1

          other_group = g2

        } else {

          panel_group = g2

          other_group = g1

        }


        panel_level_col_map = panel_group$level_col_map

        other_level_col_map = other_group$level_col_map


        interaction_col_map = list()

        matched_cols = character(0)

        for(panel_lev in names(panel_level_col_map)) {

          panel_col = panel_level_col_map[[panel_lev]]

          interaction_col_map[[panel_lev]] = list()


          for(other_lev in names(other_level_col_map)) {

            other_col = other_level_col_map[[other_lev]]

            candidates = c(

              paste0(other_col, ":", panel_col),

              paste0(panel_col, ":", other_col)

            )

            hit = candidates[candidates %in% x_cols]

            if(length(hit) > 0) {

              interaction_col_map[[panel_lev]][[other_lev]] = hit[1]

              matched_cols = c(matched_cols, hit[1])

            }

          }

        }


        if(length(matched_cols) > 0) {

          interaction_name = paste(other_group$var_name, panel_group$var_name, sep = ":")

          groups[[interaction_name]] = list(

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


  plot_specs = list()

  for(par_name in names(object$model_matrix$x)) {

    X = object$model_matrix$x[[par_name]]

    if(is.null(X) || ncol(X) == 0) next


    factor_groups = build_factor_groups(X, data_for_terms)

    interaction_groups = if(plot_interactions) build_factor_interaction_groups(X, factor_groups) else list()

    grouped_cols = unique(unlist(lapply(factor_groups, function(g) g$matched_cols), use.names = FALSE))

    if(length(grouped_cols) == 0) grouped_cols = character(0)

    grouped_interaction_cols = unique(unlist(lapply(interaction_groups, function(g) g$matched_cols), use.names = FALSE))

    if(length(grouped_interaction_cols) == 0) grouped_interaction_cols = character(0)


    for(var_name in names(factor_groups)) {

      fg = factor_groups[[var_name]]

      has_valid_coef = FALSE

      for(lev in names(fg$level_col_map)) {

        coef_name = paste(par_name, fg$level_col_map[[lev]], sep = ".")

        if(coef_name %in% names(object$par) && coef_name %in% rownames(V) && coef_name %in% colnames(V)) {

          has_valid_coef = TRUE

          break

        }

      }

      if(has_valid_coef) {

        plot_specs[[length(plot_specs) + 1]] = list(

          type = "factor",

          par_name = par_name,

          var_name = var_name,

          group = fg

        )

      }

    }


    for(inter_name in names(interaction_groups)) {

      ig = interaction_groups[[inter_name]]

      has_valid_coef = FALSE

      for(panel_lev in names(ig$interaction_col_map)) {

        for(other_lev in names(ig$interaction_col_map[[panel_lev]])) {

          coef_name = paste(par_name, ig$interaction_col_map[[panel_lev]][[other_lev]], sep = ".")

          if(coef_name %in% names(object$par) && coef_name %in% rownames(V) && coef_name %in% colnames(V)) {

            has_valid_coef = TRUE

            break

          }

        }

        if(has_valid_coef) break

      }

      if(has_valid_coef) {

        plot_specs[[length(plot_specs) + 1]] = list(

          type = "interaction_factor_factor",

          par_name = par_name,

          group = ig

        )

      }

    }


    for(col_name in colnames(X)) {

      if(col_name %in% grouped_cols) next

      if(col_name %in% grouped_interaction_cols) next

      if(!include_intercept && col_name == "intercept") next

      if(!plot_interactions && grepl(":", col_name, fixed = TRUE)) next

      coef_name = paste(par_name, col_name, sep = ".")

      if(!coef_name %in% names(object$par)) next

      if(!coef_name %in% rownames(V) || !coef_name %in% colnames(V)) next


      term_type = if(grepl(":", col_name, fixed = TRUE)) "interaction_factor" else "continuous"


      plot_specs[[length(plot_specs) + 1]] = list(

        type = term_type,

        par_name = par_name,

        col_name = col_name,

        coef_name = coef_name

      )

    }

  }


  n_plots = length(plot_specs)

  if(n_plots == 0) {

    warning("No fixed terms found to plot with matching vcov entries.")

    return(invisible(list()))

  }


  out = list()

  plot_objects = list()

  for(i in seq_len(n_plots)) {

    spec = plot_specs[[i]]

    par_name = spec$par_name


    if(identical(spec$type, "factor")) {

      fg = spec$group

      levs = fg$levels

      x_plot = seq_along(levs)

      fitted_term = rep(NA_real_, length(levs))

      term_se = rep(NA_real_, length(levs))


      for(j in seq_along(levs)) {

        lev = levs[j]

        if(identical(lev, fg$ref_level)) {

          fitted_term[j] = 0

          term_se[j] = 0

        } else if(lev %in% names(fg$level_col_map)) {

          col_name_lev = fg$level_col_map[[lev]]

          coef_name_lev = paste(par_name, col_name_lev, sep = ".")

          if(coef_name_lev %in% names(object$par) && coef_name_lev %in% rownames(V) && coef_name_lev %in% colnames(V)) {

            fitted_term[j] = as.numeric(object$par[coef_name_lev])

            term_se[j] = sqrt(pmax(0, as.numeric(V[coef_name_lev, coef_name_lev])))

          }

        }

      }


      keep = is.finite(fitted_term) & is.finite(term_se)

      ci_lower = fitted_term - z * term_se

      ci_upper = fitted_term + z * term_se


      y_vals = c(fitted_term[keep], ci_lower[keep], ci_upper[keep])

      if(length(y_vals) > 0) {

        y_rng = range(y_vals)

        y_pad = 0.05 * max(1e-8, diff(y_rng))

        y_lim = c(y_rng[1] - y_pad, y_rng[2] + y_pad)

      } else {

        y_lim = NULL

      }


      plot_df = data.frame(

        x = x_plot,

        fitted = fitted_term,

        ci_lower = ci_lower,

        ci_upper = ci_upper,

        keep = keep

      )


      p = ggplot2::ggplot(plot_df[plot_df$keep, , drop = FALSE], ggplot2::aes(x = x, y = fitted))

      p = gg_add(p, ggplot2::geom_hline(yintercept = 0, color = "grey70", linetype = 3), "geom_hline")

      p = gg_add(p, ggplot2::geom_point(color = fit_col, size = factor_cex), "geom_point")

      p = gg_add(p, ggplot2::geom_errorbar(ggplot2::aes(ymin = ci_lower, ymax = ci_upper), color = ci_col, width = 0.15), "geom_errorbar")

      p = gg_add(p, ggplot2::scale_x_continuous(breaks = x_plot, labels = levs), "scale_x_continuous")

      p = gg_add(

        p,

        ggplot2::labs(

          title = paste(par_name, fg$var_name, sep = ": "),

          x = fg$var_name,

          y = paste("fixed contribution:", paste(par_name, fg$var_name, sep = "."))

        ),

        "labs"

      )

      p = gg_add(p, ggplot2::theme_minimal(), "theme_minimal")


      if(!is.null(y_lim)) {

        p = gg_add(p, ggplot2::coord_cartesian(ylim = y_lim), "coord_cartesian")

      }


      if(show_legend) {

        p = gg_add(p, ggplot2::labs(caption = paste("estimate /", round(ci_level * 100), "% CI")), "labs")

      }


      if(is.null(out[[par_name]])) out[[par_name]] = list()

      out[[par_name]][[fg$var_name]] = list(

        coefficient = paste(par_name, fg$var_name, sep = "."),

        x = x_plot,

        levels = levs,

        fitted = fitted_term,

        se = term_se,

        ci_lower = ci_lower,

        ci_upper = ci_upper,

        plot = p

      )

      plot_objects[[length(plot_objects) + 1]] = p

    } else if(identical(spec$type, "interaction_factor_factor")) {

      ig = spec$group

      pg = ig$panel_group

      og = ig$other_group


      panel_levels = pg$levels

      other_levels = og$levels

      x_plot = seq_along(panel_levels)

      x_labels = panel_levels


      plot_rows = list()

      for(other_lev in other_levels) {

        fitted_term = rep(NA_real_, length(panel_levels))

        term_se = rep(NA_real_, length(panel_levels))


        for(j in seq_along(panel_levels)) {

          panel_lev = panel_levels[j]

          if(identical(panel_lev, pg$ref_level) || identical(other_lev, og$ref_level)) {

            fitted_term[j] = 0

            term_se[j] = 0

          } else if(panel_lev %in% names(ig$interaction_col_map) && other_lev %in% names(ig$interaction_col_map[[panel_lev]])) {

            col_name_lev = ig$interaction_col_map[[panel_lev]][[other_lev]]

            coef_name_lev = paste(par_name, col_name_lev, sep = ".")

            if(coef_name_lev %in% names(object$par) && coef_name_lev %in% rownames(V) && coef_name_lev %in% colnames(V)) {

              fitted_term[j] = as.numeric(object$par[coef_name_lev])

              term_se[j] = sqrt(pmax(0, as.numeric(V[coef_name_lev, coef_name_lev])))

            }

          }

        }


        keep = is.finite(fitted_term) & is.finite(term_se)

        ci_lower = fitted_term - z * term_se

        ci_upper = fitted_term + z * term_se


        plot_rows[[length(plot_rows) + 1]] = data.frame(

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


      plot_df = do.call(rbind, plot_rows)

      y_vals = c(plot_df$fitted[plot_df$keep], plot_df$ci_lower[plot_df$keep], plot_df$ci_upper[plot_df$keep])

      if(length(y_vals) > 0) {

        y_rng = range(y_vals)

        y_pad = 0.05 * max(1e-8, diff(y_rng))

        y_lim = c(y_rng[1] - y_pad, y_rng[2] + y_pad)

      } else {

        y_lim = NULL

      }


      p = ggplot2::ggplot(plot_df[plot_df$keep, , drop = FALSE], ggplot2::aes(x = x, y = fitted, color = group, group = group))

      p = gg_add(p, ggplot2::geom_hline(yintercept = 0, color = "grey70", linetype = 3), "geom_hline")

      p = gg_add(p, ggplot2::geom_line(linewidth = 0.8), "geom_line")

      p = gg_add(p, ggplot2::geom_point(size = factor_cex), "geom_point")

      p = gg_add(p, ggplot2::geom_errorbar(ggplot2::aes(ymin = ci_lower, ymax = ci_upper), width = 0.15), "geom_errorbar")

      p = gg_add(p, ggplot2::scale_x_continuous(breaks = x_plot, labels = x_labels), "scale_x_continuous")

      p = gg_add(p, ggplot2::scale_color_discrete(name = og$var_name), "scale_color_discrete")

      p = gg_add(

        p,

        ggplot2::labs(

          title = paste(par_name, ig$interaction_name, sep = ": "),

          x = pg$var_name,

          y = paste("fixed contribution:", paste(par_name, ig$interaction_name, sep = "."))

        ),

        "labs"

      )

      p = gg_add(p, ggplot2::theme_minimal(), "theme_minimal")


      if(!is.null(y_lim)) {

        p = gg_add(p, ggplot2::coord_cartesian(ylim = y_lim), "coord_cartesian")

      }


      if(show_legend) {

        p = gg_add(p, ggplot2::labs(caption = paste("estimate /", round(ci_level * 100), "% CI")), "labs")

      }


      if(is.null(out[[par_name]])) out[[par_name]] = list()

      out[[par_name]][[ig$interaction_name]] = list(

        coefficient = paste(par_name, ig$interaction_name, sep = "."),

        x = x_plot,

        levels = panel_levels,

        series = other_levels,

        fitted = plot_df$fitted,

        se = plot_df$se,

        ci_lower = plot_df$ci_lower,

        ci_upper = plot_df$ci_upper,

        plot_data = plot_df,

        plot = p

      )

      plot_objects[[length(plot_objects) + 1]] = p

    } else if(identical(spec$type, "interaction_factor")) {

      col_name = spec$col_name

      coef_name = spec$coef_name

      X = object$model_matrix$x[[par_name]]

      x_raw = as.numeric(X[, col_name])

      beta_hat = as.numeric(object$par[coef_name])

      var_beta = as.numeric(V[coef_name, coef_name])


      x_levels = sort(unique(x_raw[is.finite(x_raw)]))

      if(length(x_levels) == 0) next


      fitted_term = x_levels * beta_hat

      term_se = abs(x_levels) * sqrt(pmax(0, var_beta))

      ci_lower = fitted_term - z * term_se

      ci_upper = fitted_term + z * term_se


      x_labels = as.character(signif(x_levels, 6))

      x_plot = seq_along(x_levels)


      y_vals = c(fitted_term, ci_lower, ci_upper)

      y_vals = y_vals[is.finite(y_vals)]

      if(length(y_vals) > 0) {

        y_rng = range(y_vals)

        y_pad = 0.05 * max(1e-8, diff(y_rng))

        y_lim = c(y_rng[1] - y_pad, y_rng[2] + y_pad)

      } else {

        y_lim = NULL

      }


      plot_df = data.frame(

        x = x_plot,

        fitted = fitted_term,

        ci_lower = ci_lower,

        ci_upper = ci_upper,

        stringsAsFactors = FALSE

      )


      p = ggplot2::ggplot(plot_df, ggplot2::aes(x = x, y = fitted))

      p = gg_add(p, ggplot2::geom_hline(yintercept = 0, color = "grey70", linetype = 3), "geom_hline")

      p = gg_add(p, ggplot2::geom_point(color = fit_col, size = factor_cex), "geom_point")

      p = gg_add(p, ggplot2::geom_errorbar(ggplot2::aes(ymin = ci_lower, ymax = ci_upper), color = ci_col, width = 0.15), "geom_errorbar")

      p = gg_add(p, ggplot2::scale_x_continuous(breaks = x_plot, labels = x_labels), "scale_x_continuous")

      p = gg_add(

        p,

        ggplot2::labs(

          title = paste(par_name, col_name, sep = ": "),

          x = paste(col_name, "(interaction level)"),

          y = paste("fixed contribution:", coef_name)

        ),

        "labs"

      )

      p = gg_add(p, ggplot2::theme_minimal(), "theme_minimal")


      if(!is.null(y_lim)) {

        p = gg_add(p, ggplot2::coord_cartesian(ylim = y_lim), "coord_cartesian")

      }


      if(show_legend) {

        p = gg_add(p, ggplot2::labs(caption = paste("estimate /", round(ci_level * 100), "% CI")), "labs")

      }


      if(is.null(out[[par_name]])) out[[par_name]] = list()

      out[[par_name]][[col_name]] = list(

        coefficient = coef_name,

        x = x_levels,

        levels = x_labels,

        fitted = fitted_term,

        se = term_se,

        ci_lower = ci_lower,

        ci_upper = ci_upper,

        plot = p

      )

      plot_objects[[length(plot_objects) + 1]] = p

    } else {

      col_name = spec$col_name

      coef_name = spec$coef_name

      X = object$model_matrix$x[[par_name]]

      x_raw = as.numeric(X[, col_name])

      beta_hat = as.numeric(object$par[coef_name])

      var_beta = as.numeric(V[coef_name, coef_name])


      fitted_term = x_raw * beta_hat

      term_se = sqrt(pmax(0, (x_raw^2) * var_beta))

      ci_lower = fitted_term - z * term_se

      ci_upper = fitted_term + z * term_se


      if(length(unique(x_raw)) <= 1 && fallback_to_index) {

        x_plot = seq_along(x_raw)

        xlab_text = paste(col_name, "(index)")

        ord = seq_along(x_plot)

      } else {

        x_plot = x_raw

        xlab_text = col_name

        ord = if(sort_x) order(x_plot) else seq_along(x_plot)

      }


      main_title = paste(par_name, col_name, sep = ": ")

      ylab_text = paste("fixed contribution:", coef_name)


      y_vals = c(fitted_term, ci_lower, ci_upper)

      y_vals = y_vals[is.finite(y_vals)]

      if(length(y_vals) > 0) {

        y_rng = range(y_vals)

        y_pad = 0.05 * max(1e-8, diff(y_rng))

        y_lim = c(y_rng[1] - y_pad, y_rng[2] + y_pad)

      } else {

        y_lim = NULL

      }


      plot_df = data.frame(

        x = x_plot[ord],

        fitted = fitted_term[ord],

        ci_lower = ci_lower[ord],

        ci_upper = ci_upper[ord]

      )


      p = ggplot2::ggplot(plot_df, ggplot2::aes(x = x, y = fitted))

      p = gg_add(p, ggplot2::geom_ribbon(ggplot2::aes(ymin = ci_lower, ymax = ci_upper), fill = ci_col, alpha = 0.16), "geom_ribbon")

      p = gg_add(p, ggplot2::geom_line(color = fit_col, linewidth = fit_lwd), "geom_line")

      p = gg_add(p, ggplot2::labs(title = main_title, x = xlab_text, y = ylab_text), "labs")

      p = gg_add(p, ggplot2::theme_minimal(), "theme_minimal")


      if(!is.null(y_lim)) {

        p = gg_add(p, ggplot2::coord_cartesian(ylim = y_lim), "coord_cartesian")

      }


      if(show_legend) {

        p = gg_add(p, ggplot2::labs(caption = paste("fit /", round(ci_level * 100), "% CI")), "labs")

      }


      if(is.null(out[[par_name]])) out[[par_name]] = list()

      out[[par_name]][[col_name]] = list(

        coefficient = coef_name,

        x = x_plot,

        fitted = fitted_term,

        se = term_se,

        ci_lower = ci_lower,

        ci_upper = ci_upper,

        plot = p

      )

      plot_objects[[length(plot_objects) + 1]] = p

    }

  }


  if(length(plot_objects) > 0) {

    if(is.null(ncol)) {

      ncol = min(2, n_plots)

    }

    nrow = ceiling(length(plot_objects) / ncol)

    dashboard = list(plotlist = plot_objects, ncol = ncol, nrow = nrow)

    if(setup_mfrow) {

      grid::grid.newpage()

      grid::pushViewport(grid::viewport(layout = grid::grid.layout(nrow, ncol)))

      for(i_plot in seq_along(plot_objects)) {

        r = ((i_plot - 1) %/% ncol) + 1

        c = ((i_plot - 1) %% ncol) + 1

        print(plot_objects[[i_plot]], vp = grid::viewport(layout.pos.row = r, layout.pos.col = c))

      }

      grid::popViewport()

    }

    out$plots = plot_objects

    out$dashboard = dashboard

  }


  invisible(out)

}


#' Plot term effects for a fitted longitudinal model

#'

#' @param x A fitted `gamlss.longitudinal` object.

#' @param y Unused; included for compatibility with older calls.

#' @param data Optional data frame used to recover factor levels, transformed

#'   covariate scales, and interaction plotting metadata. Factor fixed terms

#'   are grouped into one point-and-interval plot per model term and parameter.

#' @param ci_level Confidence level for pointwise intervals.

#' @param ncol Number of columns in the combined dashboard.

#' @param include_intercept Logical; include intercept terms in fixed-effect

#'   plots.

#' @param plot_interactions Logical; include fixed-effect interaction terms.

#' @param ci_col,fit_col Colours for interval and fitted-term layers.

#' @param show_legend Logical; include plot captions/legends where available.

#' @param smooth_even_grid Logical; draw smooth terms on an evenly spaced grid.

#' @param smooth_grid_n Number of grid points for smooth-term plots.

#' @param paginate Logical; print one chart at a time for large dashboards.

#' @param ... Additional arguments reserved for future use.

#'

#' @return Invisibly returns a list with smooth-term, fixed-term, and dashboard

#'   plot objects.

#' @export

plot_terms <- function(

  x,

  y = NULL,

  data = NULL,

  ci_level = 0.95,

  ncol = 4,

  include_intercept = FALSE,

  plot_interactions = FALSE,

  ci_col = "red",

  fit_col = "black",

  show_legend = TRUE,

  smooth_even_grid = TRUE,

  smooth_grid_n = 200,

  paginate = FALSE,

  ...

) {

  .plot_terms_gamlss_longitudinal(

    x = x,

    y = y,

    data = data,

    ci_level = ci_level,

    ncol = ncol,

    include_intercept = include_intercept,

    plot_interactions = plot_interactions,

    ci_col = ci_col,

    fit_col = fit_col,

    show_legend = show_legend,

    smooth_even_grid = smooth_even_grid,

    smooth_grid_n = smooth_grid_n,

    paginate = paginate,

    ...

  )

}


#' @rdname plot_terms

#' @usage NULL

#' @rawNamespace export(plot.terms)

plot.terms <- function(x, ...) {

  .Deprecated("plot_terms", package = "gamlss.longitudinal")

  plot_terms(x, ...)

}


.plot_terms_gamlss_longitudinal <- function(

  x,

  y = NULL,

  data = NULL,

  ci_level = 0.95,

  ncol = 4,

  include_intercept = FALSE,

  plot_interactions = FALSE,

  ci_col = "red",

  fit_col = "black",

  show_legend = TRUE,

  smooth_even_grid = TRUE,

  smooth_grid_n = 200,

  paginate = FALSE,

  ...

) {

  if(!inherits(x, "gamlss.longitudinal")) {

    stop("'x' must be a fitted 'gamlss.longitudinal' object.")

  }


  cat("\n=== Plotting term effects for gamlss.longitudinal object ===\n")


  count_plot_terms = function(obj) {

    n_smooth = 0

    n_fixed = 0

    data_for_terms = data

    if((is.null(data_for_terms) || !is.data.frame(data_for_terms)) && !is.null(obj$dataset)) {

      data_for_terms = obj$dataset

    }


    clean_expr_name = function(x) {

      x = trimws(x)

      x = gsub("`", "", x, fixed = TRUE)

      factor_match = regexec("^(?:as\\.)?factor\\(([^)]+)\\)$", x)

      matched = regmatches(x, factor_match)[[1]]

      if(length(matched) >= 2) {

        x = trimws(matched[2])

      }

      x

    }


    is_factor_term = function(term_name) {

      term_var = clean_expr_name(term_name)

      if(grepl("^(?:as\\.)?factor\\(", term_name)) return(TRUE)

      if(is.null(data_for_terms) || !is.data.frame(data_for_terms)) return(FALSE)

      candidates = unique(c(

        term_var,

        make.names(term_var),

        if(identical(term_var, "time_covariate")) names(data_for_terms)[grepl("time", names(data_for_terms), ignore.case = TRUE)] else character(0)

      ))

      any(candidates %in% names(data_for_terms) & vapply(candidates, function(candidate) {

        candidate %in% names(data_for_terms) && is.factor(data_for_terms[[candidate]])

      }, logical(1)))

    }


    for(par_name in names(obj$par_s)) {

      if(length(obj$par_s[[par_name]]) > 0) {

        n_smooth = n_smooth + length(obj$par_s[[par_name]])

      }

    }


    for(par_name in names(obj$model_matrix$x)) {

      X = obj$model_matrix$x[[par_name]]

      if(!is.null(X) && ncol(X) > 0) {

        x_cols = colnames(X)

        assign = attr(X, "assign")

        term_labels = attr(X, "term.labels")

        grouped_cols = character(0)

        if(!is.null(assign) && !is.null(term_labels) && length(assign) == length(x_cols)) {

          for(term_idx in seq_along(term_labels)) {

            term_name = term_labels[term_idx]

            if(grepl(":", term_name, fixed = TRUE)) next

            term_cols = x_cols[assign == term_idx]

            if(length(term_cols) == 0 || !is_factor_term(term_name)) next

            coef_names = paste(par_name, term_cols, sep = ".")

            if(any(coef_names %in% names(obj$par))) {

              n_fixed = n_fixed + 1

              grouped_cols = c(grouped_cols, term_cols)

            }

          }

        }


        keep_cols = !(x_cols == "intercept" & !include_intercept)

        keep_cols = keep_cols & !(x_cols %in% grouped_cols)

        coef_names = paste(par_name, x_cols[keep_cols], sep = ".")

        if(!plot_interactions) {

          coef_names = coef_names[!grepl(":", coef_names, fixed = TRUE)]

        }

        n_fixed = n_fixed + sum(coef_names %in% names(obj$par))

      }

    }


    list(smooth = n_smooth, fixed = n_fixed, total = n_smooth + n_fixed)

  }


  counts = count_plot_terms(x)

  cat(sprintf("Found %d smooth terms and %d fixed terms (total: %d plots).\n\n",

              counts$smooth, counts$fixed, counts$total))


  if(counts$total == 0) {

    warning("No term plots to display.")

    return(invisible(list(smooth_terms = list(), fixed_terms = list())))

  }


  vcov_obj = .resolve_vcov(x, numderiv = FALSE, extra_args = list(method = "analytical"))


  smooth_results = list()

  fixed_results = list()

  plot_objects = list()


  if(counts$smooth > 0) {

    smooth_results = plot_smooth_terms(

      object = x,

      vcov_obj = vcov_obj,

      data = data,

      ci_level = ci_level,

      ncol = ncol,

      ci_col = ci_col,

      fit_col = fit_col,

      even_grid = smooth_even_grid,

      grid_n = smooth_grid_n,

      setup_mfrow = FALSE,

      show_legend = show_legend

    )

    if(!is.null(smooth_results$plots)) {

      plot_objects = c(plot_objects, smooth_results$plots)

    }

  }


  if(counts$fixed > 0) {

    fixed_results = plot_fixed_terms(

      object = x,

      vcov_obj = vcov_obj,

      ci_level = ci_level,

      ncol = ncol,

      include_intercept = include_intercept,

      plot_interactions = plot_interactions,

      ci_col = ci_col,

      fit_col = fit_col,

      setup_mfrow = FALSE,

      data = data,

      show_legend = show_legend

    )

    if(!is.null(fixed_results$plots)) {

      plot_objects = c(plot_objects, fixed_results$plots)

    }

  }


  dashboard = NULL

  if(length(plot_objects) > 0) {

    if(length(plot_objects) > 16 && !isTRUE(paginate)) {

      warning(

        "More than 16 charts (", length(plot_objects), ") were generated. ",

        "Rendering all charts at once may fail in some environments. ",

        "Use paginate=TRUE to view one chart at a time.",

        call. = FALSE

      )

    }


    if(isTRUE(paginate)) {

      dashboard = list(plotlist = plot_objects, paginate = TRUE)

      for(i_plot in seq_along(plot_objects)) {

        grid::grid.newpage()

        print(plot_objects[[i_plot]])

        if(i_plot < length(plot_objects) && interactive()) {

          invisible(readline(prompt = sprintf("Press [Enter] for next chart (%d/%d)... ", i_plot, length(plot_objects))))

        }

      }

    } else {

      if(is.null(ncol)) {

        ncol = min(2, length(plot_objects))

      }

      nrow = ceiling(length(plot_objects) / ncol)

      dashboard = list(plotlist = plot_objects, ncol = ncol, nrow = nrow, paginate = FALSE)


      grid::grid.newpage()

      grid::pushViewport(grid::viewport(layout = grid::grid.layout(nrow, ncol)))

      for(i_plot in seq_along(plot_objects)) {

        r = ((i_plot - 1) %/% ncol) + 1

        c = ((i_plot - 1) %% ncol) + 1

        print(plot_objects[[i_plot]], vp = grid::viewport(layout.pos.row = r, layout.pos.col = c))

      }

      grid::popViewport()

    }

  }


  invisible(list(

    smooth_terms = smooth_results,

    fixed_terms = fixed_results,

    dashboard = dashboard

  ))

}


