source(file.path("paper", "R", "08-simulation-sensitivity-correlation-misspecification", "standard-model-benchmarking", "00-benchmark-setup.R"))

run_dir <- bmk_env(
  "GAMLSS_LONGITUDINAL_BENCHMARK_REPORT_DIR",
  bmk_env(
    "GAMLSS_LONGITUDINAL_BENCHMARK_RUN_DIR",
    {
      latest_combined <- file.path(bmk_output_root, "latest_combined_review_dir.txt")
      latest <- file.path(bmk_output_root, "latest_run_dir.txt")
      if (file.exists(latest_combined)) {
        trimws(readLines(latest_combined, warn = FALSE)[1L])
      } else if (file.exists(latest)) {
        trimws(readLines(latest, warn = FALSE)[1L])
      } else {
        ""
      }
    }
  )
)
if (!nzchar(run_dir)) {
  stop("No run directory supplied. Set GAMLSS_LONGITUDINAL_BENCHMARK_REPORT_DIR or GAMLSS_LONGITUDINAL_BENCHMARK_RUN_DIR.", call. = FALSE)
}
run_dir <- normalizePath(run_dir, winslash = "/", mustWork = TRUE)
out_dir <- file.path(run_dir, "story_tables")
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)
unlink(c(
  Sys.glob(file.path(out_dir, "table_*.csv")),
  Sys.glob(file.path(out_dir, "table_*.tex")),
  Sys.glob(file.path(out_dir, "appendix_table_*.csv")),
  Sys.glob(file.path(out_dir, "appendix_table_*.tex")),
  file.path(out_dir, "appendix_t50_tables.tex"),
  file.path(out_dir, "story_tables.tex")
), force = TRUE)

read_csv <- function(name, required = TRUE) {
  path <- file.path(run_dir, name)
  if (!file.exists(path)) {
    if (required) stop("Missing required file: ", path, call. = FALSE)
    return(data.frame())
  }
  utils::read.csv(path, stringsAsFactors = FALSE, check.names = FALSE)
}

read_sandwich_csv <- function(base_name, sandwich_name, required = TRUE) {
  sandwich_dir <- file.path(
    run_dir,
    bmk_env("GAMLSS_LONGITUDINAL_STORY_TABLE_SANDWICH_DIR", "sandwich_t20_grid")
  )
  sandwich_path <- file.path(sandwich_dir, sandwich_name)
  if (file.exists(sandwich_path)) {
    return(utils::read.csv(sandwich_path, stringsAsFactors = FALSE, check.names = FALSE))
  }
  read_csv(base_name, required = required)
}

add_n_time <- function(x) {
  if (nrow(x) > 0L && !"n_time" %in% names(x)) x$n_time <- 4L
  x
}

method_order <- c("rs_joint", "rs_joint_sandwich", "gee_exchangeable", "gee_ar1", "gee_unstructured", "glm")
method_labels <- c(
  rs_joint = "gamlss.longitudinal",
  rs_joint_sandwich = "gamlss.longitudinal robust",
  gee_exchangeable = "geepack exchangeable",
  gee_ar1 = "geepack AR(1)",
  gee_unstructured = "geepack unstructured",
  glm = "glm"
)
method_column_labels <- c(
  rs_joint = "gamlss.long",
  rs_joint_sandwich = "robust",
  gee_exchangeable = "exch.",
  gee_ar1 = "AR(1)",
  gee_unstructured = "unstr.",
  glm = "glm"
)
scenario_labels <- c(
  external_exchangeable_moderate = "External exchangeable, moderate",
  external_exchangeable_high = "External exchangeable, high",
  external_ar1_moderate = "External AR(1), moderate",
  external_ar1_high = "External AR(1), high",
  internal_time_varying_high = "Internal time-varying adjacent",
  internal_covariate_dependent_high = "Internal covariate-dependent adjacent"
)
scenario_ids <- c(
  external_exchangeable_moderate = "EXCH-M",
  external_exchangeable_high = "EXCH-H",
  external_ar1_moderate = "AR1-M",
  external_ar1_high = "AR1-H",
  internal_time_varying_high = "TIME-COR",
  internal_covariate_dependent_high = "COV-COR"
)
dependence_labels <- c(
  exchangeable = "Exchangeable",
  ar1 = "AR(1)",
  time_varying_adjacent = "Time-varying adjacent",
  covariate_dependent_adjacent = "Covariate-dependent adjacent"
)
family_order <- bmk_env_vector("GAMLSS_LONGITUDINAL_STORY_TABLE_FAMILIES", c("gaussian", "gamma", "binary"))
family_order <- intersect(family_order, c("gaussian", "gamma", "binary", "poisson"))
family_labels <- c(gaussian = "Normal", gamma = "Gamma", binary = "Binary", poisson = "Poisson")
n_subjects <- bmk_env_int(
  "GAMLSS_LONGITUDINAL_STORY_TABLE_N",
  bmk_env_int("GAMLSS_LONGITUDINAL_BENCHMARK_N", 120L)
)

base_scenario <- function(x) sub("_t[0-9]+$", "", as.character(x))
scenario_code <- function(scenario) {
  base <- base_scenario(scenario)
  code <- unname(scenario_ids[base])
  code[is.na(code)] <- base[is.na(code)]
  code
}
scenario_display <- function(scenario, n_time) {
  base <- base_scenario(scenario)
  label <- unname(scenario_labels[base])
  label[is.na(label)] <- base[is.na(label)]
  paste0("T=", n_time, ": ", label)
}
fmt_num <- function(x, digits = 3L) ifelse(is.na(x) | !is.finite(x), "--", formatC(x, digits = digits, format = "f"))
fmt_int <- function(x) ifelse(is.na(x) | !is.finite(x), "--", as.character(as.integer(round(x))))
ok <- function(x) x %in% c(TRUE, "TRUE", "True", "true", "1")
fmt_range <- function(x, digits = 1L) {
  x <- x[is.finite(x)]
  if (length(x) == 0L) return("--")
  r <- range(x)
  if (isTRUE(all.equal(r[[1L]], r[[2L]]))) return(fmt_num(r[[1L]], digits))
  paste0(fmt_num(r[[1L]], digits), "-", fmt_num(r[[2L]], digits))
}
df_range_note_by_method <- function(df, methods, labels, digits = 1L) {
  parts <- vapply(methods, function(method) {
    method_df <- df[df$method == method, , drop = FALSE]
    paste0(
      labels[[method]], " ",
      fmt_range(method_df$total_df_median, digits), "/",
      fmt_range(method_df$dependence_df_median, digits)
    )
  }, character(1L))
  paste0(" Median df ranges (total/dependence): ", paste(parts, collapse = "; "), ".")
}
esc <- function(x) {
  x <- as.character(x)
  x <- gsub("\\\\", "\\\\textbackslash{}", x)
  x <- gsub("_", "\\\\_", x, fixed = TRUE)
  x <- gsub("%", "\\\\%", x, fixed = TRUE)
  x <- gsub("&", "\\\\&", x, fixed = TRUE)
  x
}
latex_tabular <- function(df, align = NULL, resize = TRUE, escape = TRUE) {
  if (!is.data.frame(df)) {
    stop("Expected a data frame for LaTeX table output.", call. = FALSE)
  }
  n_rows <- dim(df)[1L]
  if (is.na(n_rows)) {
    stop("Expected a two-dimensional data frame for LaTeX table output.", call. = FALSE)
  }
  if (is.null(align)) align <- paste0("l", paste(rep("r", ncol(df) - 1L), collapse = ""))
  line_end <- paste0(" ", "\\", "\\")
  cell <- if (isTRUE(escape)) esc else as.character
  tabular <- c(
    paste0("\\begin{tabular}{", align, "}"),
    "\\toprule",
    paste(cell(names(df)), collapse = " & ") |> paste0(line_end),
    "\\midrule"
  )
  for (i in seq_len(n_rows)) {
    tabular <- c(tabular, paste(cell(unlist(df[i, ], use.names = FALSE)), collapse = " & ") |> paste0(line_end))
  }
  tabular <- c(tabular, "\\bottomrule", "\\end{tabular}")
  if (!resize) return(tabular)
  c("\\resizebox{\\textwidth}{!}{%", tabular, "}")
}

write_table <- function(df, stem, caption, label, align = NULL, resize = TRUE, escape = TRUE) {
  csv_path <- file.path(out_dir, paste0(stem, ".csv"))
  tex_path <- file.path(out_dir, paste0(stem, ".tex"))
  utils::write.csv(df, csv_path, row.names = FALSE, na = "")
  tex <- c(
    paste0("\\begin{table}[htbp]"),
    "\\centering",
    paste0("\\caption{", caption, "}"),
    paste0("\\label{", label, "}"),
    latex_tabular(df, align = align, resize = resize, escape = escape),
    "\\end{table}"
  )
  writeLines(tex, tex_path, useBytes = TRUE)
  invisible(list(csv = csv_path, tex = tex_path))
}

pivot_metric_table <- function(df, row_cols, metrics, methods, method_labels_for_columns, digits = 3L) {
  keys <- unique(df[row_cols])
  out <- vector("list", length = nrow(keys) * length(metrics))
  k <- 0L
  for (i in seq_len(nrow(keys))) {
    row_match <- rep(TRUE, nrow(df))
    for (col in row_cols) {
      row_match <- row_match & df[[col]] == keys[[col]][[i]]
    }
    row_df <- df[row_match, , drop = FALSE]
    for (metric_name in names(metrics)) {
      k <- k + 1L
      item <- keys[i, , drop = FALSE]
      item$Metric <- unname(metrics[[metric_name]])
      for (method in methods) {
        values <- row_df[[metric_name]][row_df$method == method]
        value <- if (length(values) == 0L) NA_real_ else values[[1L]]
        item[[unname(method_labels_for_columns[[method]])]] <- fmt_num(value, digits)
      }
      out[[k]] <- item
    }
  }
  do.call(rbind, out)
}

pivot_metric_time_table <- function(df, row_cols, time_col, metrics, methods, method_labels_for_columns, times = NULL, digits = 3L) {
  if (is.null(times)) times <- sort(unique(df[[time_col]]))
  keys <- unique(df[row_cols])
  out <- vector("list", length = nrow(keys) * length(metrics))
  k <- 0L
  for (i in seq_len(nrow(keys))) {
    key_match <- rep(TRUE, nrow(df))
    for (col in row_cols) {
      key_match <- key_match & df[[col]] == keys[[col]][[i]]
    }
    key_df <- df[key_match, , drop = FALSE]
    for (metric_name in names(metrics)) {
      k <- k + 1L
      item <- keys[i, , drop = FALSE]
      item$Metric <- unname(metrics[[metric_name]])
      for (method in methods) {
        method_label <- unname(method_labels_for_columns[[method]])
        for (time in times) {
          values <- key_df[[metric_name]][key_df$method == method & key_df[[time_col]] == time]
          value <- if (length(values) == 0L) NA_real_ else values[[1L]]
          item[[paste0(method_label, " T=", time)]] <- fmt_num(value, digits)
        }
      }
      out[[k]] <- item
    }
  }
  do.call(rbind, out)
}

latex_time_pivot_tabular <- function(df, row_cols, methods, method_labels_for_columns, times, resize = TRUE) {
  fixed_cols <- row_cols
  n_fixed <- length(fixed_cols)
  n_times <- length(times)
  n_cols <- n_fixed + length(methods) * n_times
  align <- paste0(paste(rep("l", n_fixed), collapse = ""), paste(rep("r", length(methods) * n_times), collapse = ""))
  line_end <- paste0(" ", "\\", "\\")
  is_gee <- grepl("^gee_", methods)
  if (n_times == 1L && any(is_gee)) {
    is_gamlss <- methods %in% c("rs_joint", "rs_joint_sandwich")
    has_gamlss_pair <- all(c("rs_joint", "rs_joint_sandwich") %in% methods)
    first_gee <- which(is_gee)[1L]
    last_gee <- utils::tail(which(is_gee), 1L)
    first_gamlss <- if (any(is_gamlss)) which(is_gamlss)[1L] else NA_integer_
    last_gamlss <- if (any(is_gamlss)) utils::tail(which(is_gamlss), 1L) else NA_integer_
    method_group_header <- character()
    i <- 1L
    while (i <= length(methods)) {
      if (isTRUE(has_gamlss_pair) && is.finite(first_gamlss) && i == first_gamlss) {
        method_group_header <- c(method_group_header, paste0("\\multicolumn{", sum(is_gamlss), "}{c}{gamlss.long}"))
        i <- last_gamlss + 1L
      } else if (i == first_gee) {
        method_group_header <- c(method_group_header, paste0("\\multicolumn{", sum(is_gee), "}{c}{geepack (robust)}"))
        i <- last_gee + 1L
      } else {
        method_group_header <- c(method_group_header, "")
        i <- i + 1L
      }
    }
    top_header <- c(
      rep("", n_fixed),
      method_group_header
    )
    sub_header <- c(
      rep("", n_fixed),
      vapply(seq_along(methods), function(i) {
        if (isTRUE(has_gamlss_pair) && identical(methods[[i]], "rs_joint")) {
          "naive"
        } else if (isTRUE(has_gamlss_pair) && identical(methods[[i]], "rs_joint_sandwich")) {
          "robust"
        } else {
          esc(method_labels_for_columns[[methods[[i]]]])
        }
      }, character(1L))
    )
    group_rules <- paste0("\\cmidrule(lr){", n_fixed + first_gee, "-", n_fixed + last_gee, "}")
    if (isTRUE(has_gamlss_pair) && is.finite(first_gamlss)) {
      group_rules <- c(
        paste0("\\cmidrule(lr){", n_fixed + first_gamlss, "-", n_fixed + last_gamlss, "}"),
        group_rules
      )
    }
  } else {
    top_header <- c(esc(fixed_cols), vapply(methods, function(method) {
      if (n_times == 1L) {
        esc(method_labels_for_columns[[method]])
      } else {
        paste0("\\multicolumn{", n_times, "}{c}{", esc(method_labels_for_columns[[method]]), "}")
      }
    }, character(1L)))
    sub_header <- NULL
    group_rules <- NULL
  }
  time_header <- c(rep("", n_fixed), rep(paste0("T=", times), times = length(methods)))
  cmidrules <- vapply(seq_along(methods), function(i) {
    start <- n_fixed + ((i - 1L) * n_times) + 1L
    end <- start + n_times - 1L
    paste0("\\cmidrule(lr){", start, "-", end, "}")
  }, character(1L))
  tabular <- c(
    paste0("\\begin{tabular}{", align, "}"),
    "\\toprule",
    paste(top_header, collapse = " & ") |> paste0(line_end),
    if (is.null(group_rules)) character() else paste(group_rules, collapse = " "),
    if (is.null(sub_header)) character() else paste(sub_header, collapse = " & ") |> paste0(line_end),
    if (n_times == 1L) character() else paste(cmidrules, collapse = " "),
    if (n_times == 1L) character() else paste(esc(time_header), collapse = " & ") |> paste0(line_end),
    "\\midrule"
  )
  value_cols <- setdiff(names(df), c(row_cols, "Metric"))
  for (metric in unique(df$Metric)) {
    tabular <- c(
      tabular,
      paste0("\\addlinespace[0.25em]"),
      paste0("\\multicolumn{", n_cols, "}{l}{\\textit{", esc(metric), "}}", line_end)
    )
    metric_df <- df[df$Metric == metric, c(row_cols, value_cols), drop = FALSE]
    for (i in seq_len(nrow(metric_df))) {
      tabular <- c(tabular, paste(esc(unlist(metric_df[i, ], use.names = FALSE)), collapse = " & ") |> paste0(line_end))
    }
  }
  tabular <- c(tabular, "\\bottomrule", "\\end{tabular}")
  if (!resize) return(tabular)
  c("\\resizebox{\\textwidth}{!}{%", tabular, "}")
}

write_time_pivot_table <- function(df, stem, caption, label, row_cols, methods, method_labels_for_columns, times, resize = TRUE) {
  csv_path <- file.path(out_dir, paste0(stem, ".csv"))
  tex_path <- file.path(out_dir, paste0(stem, ".tex"))
  utils::write.csv(df, csv_path, row.names = FALSE, na = "")
  tex <- c(
    paste0("\\begin{table}[htbp]"),
    "\\centering",
    paste0("\\caption{", caption, "}"),
    paste0("\\label{", label, "}"),
    latex_time_pivot_tabular(
      df,
      row_cols = row_cols,
      methods = methods,
      method_labels_for_columns = method_labels_for_columns,
      times = times,
      resize = resize
    ),
    "\\end{table}"
  )
  writeLines(tex, tex_path, useBytes = TRUE)
  invisible(list(csv = csv_path, tex = tex_path))
}

pivot_metric_family_time_table <- function(df, row_cols, family_col, time_col, metrics, methods, method_labels_for_columns, families, times = NULL, digits = 3L) {
  if (is.null(times)) times <- sort(unique(df[[time_col]]))
  keys <- unique(df[row_cols])
  out <- vector("list", length = nrow(keys) * length(metrics))
  k <- 0L
  for (i in seq_len(nrow(keys))) {
    key_match <- rep(TRUE, nrow(df))
    for (col in row_cols) {
      key_match <- key_match & df[[col]] == keys[[col]][[i]]
    }
    key_df <- df[key_match, , drop = FALSE]
    for (metric_name in names(metrics)) {
      k <- k + 1L
      item <- keys[i, , drop = FALSE]
      item$Metric <- unname(metrics[[metric_name]])
      for (family in families) {
        for (method in methods) {
          method_label <- unname(method_labels_for_columns[[method]])
          for (time in times) {
            values <- key_df[[metric_name]][
              key_df[[family_col]] == family &
                key_df$method == method &
                key_df[[time_col]] == time
            ]
            value <- if (length(values) == 0L) NA_real_ else values[[1L]]
            item[[paste0(family, " ", method_label, " T=", time)]] <- fmt_num(value, digits)
          }
        }
      }
      out[[k]] <- item
    }
  }
  do.call(rbind, out)
}

latex_family_time_pivot_tabular <- function(df, row_cols, families, methods, method_labels_for_columns, times, resize = TRUE) {
  fixed_cols <- row_cols
  n_fixed <- length(fixed_cols)
  n_times <- length(times)
  n_method_cols <- length(methods) * n_times
  n_cols <- n_fixed + length(families) * n_method_cols
  align <- paste0(paste(rep("l", n_fixed), collapse = ""), paste(rep("r", length(families) * n_method_cols), collapse = ""))
  line_end <- paste0(" ", "\\", "\\")
  family_header <- c(esc(fixed_cols), vapply(families, function(family) {
    paste0("\\multicolumn{", n_method_cols, "}{c}{", esc(family), "}")
  }, character(1L)))
  method_header <- c(rep("", n_fixed), unlist(lapply(families, function(family) {
    if (n_times == 1L) {
      esc(unname(method_labels_for_columns[methods]))
    } else {
      vapply(methods, function(method) {
        paste0("\\multicolumn{", n_times, "}{c}{", esc(method_labels_for_columns[[method]]), "}")
      }, character(1L))
    }
  }), use.names = FALSE))
  time_header <- c(rep("", n_fixed), rep(rep(paste0("T=", times), times = length(methods)), times = length(families)))
  family_rules <- vapply(seq_along(families), function(i) {
    start <- n_fixed + ((i - 1L) * n_method_cols) + 1L
    end <- start + n_method_cols - 1L
    paste0("\\cmidrule(lr){", start, "-", end, "}")
  }, character(1L))
  method_rules <- unlist(lapply(seq_along(families), function(family_i) {
    vapply(seq_along(methods), function(method_i) {
      start <- n_fixed + ((family_i - 1L) * n_method_cols) + ((method_i - 1L) * n_times) + 1L
      end <- start + n_times - 1L
      paste0("\\cmidrule(lr){", start, "-", end, "}")
    }, character(1L))
  }), use.names = FALSE)
  tabular <- c(
    paste0("\\begin{tabular}{", align, "}"),
    "\\toprule",
    paste(family_header, collapse = " & ") |> paste0(line_end),
    paste(family_rules, collapse = " "),
    paste(method_header, collapse = " & ") |> paste0(line_end),
    if (n_times == 1L) character() else paste(method_rules, collapse = " "),
    if (n_times == 1L) character() else paste(esc(time_header), collapse = " & ") |> paste0(line_end),
    "\\midrule"
  )
  value_cols <- setdiff(names(df), c(row_cols, "Metric"))
  for (metric in unique(df$Metric)) {
    tabular <- c(
      tabular,
      paste0("\\addlinespace[0.25em]"),
      paste0("\\multicolumn{", n_cols, "}{l}{\\textit{", esc(metric), "}}", line_end)
    )
    metric_df <- df[df$Metric == metric, c(row_cols, value_cols), drop = FALSE]
    for (i in seq_len(nrow(metric_df))) {
      tabular <- c(tabular, paste(esc(unlist(metric_df[i, ], use.names = FALSE)), collapse = " & ") |> paste0(line_end))
    }
  }
  tabular <- c(tabular, "\\bottomrule", "\\end{tabular}")
  if (!resize) return(tabular)
  c("\\resizebox{\\textwidth}{!}{%", tabular, "}")
}

write_family_time_pivot_table <- function(df, stem, caption, label, row_cols, families, methods, method_labels_for_columns, times, resize = TRUE) {
  csv_path <- file.path(out_dir, paste0(stem, ".csv"))
  tex_path <- file.path(out_dir, paste0(stem, ".tex"))
  utils::write.csv(df, csv_path, row.names = FALSE, na = "")
  tex <- c(
    paste0("\\begin{table}[htbp]"),
    "\\centering",
    paste0("\\caption{", caption, "}"),
    paste0("\\label{", label, "}"),
    latex_family_time_pivot_tabular(
      df,
      row_cols = row_cols,
      families = families,
      methods = methods,
      method_labels_for_columns = method_labels_for_columns,
      times = times,
      resize = resize
    ),
    "\\end{table}"
  )
  writeLines(tex, tex_path, useBytes = TRUE)
  invisible(list(csv = csv_path, tex = tex_path))
}

pivot_external_correlation_table <- function(df, correlation_groups, levels, time_col, metrics, methods, method_labels_for_columns, time, digits = 3L) {
  out <- vector("list", length(levels) * length(metrics))
  k <- 0L
  for (metric_name in names(metrics)) {
    for (level_name in names(levels)) {
      k <- k + 1L
      item <- data.frame(Level = unname(levels[[level_name]]), Metric = unname(metrics[[metric_name]]), stringsAsFactors = FALSE)
      for (group_name in names(correlation_groups)) {
        scenario_base <- unname(correlation_groups[[group_name]][[level_name]])
        row_df <- df[base_scenario(df$scenario) == scenario_base & df[[time_col]] == time, , drop = FALSE]
        for (method in methods) {
          method_label <- unname(method_labels_for_columns[[method]])
          values <- row_df[[metric_name]][row_df$method == method]
          value <- if (length(values) == 0L) NA_real_ else values[[1L]]
          item[[paste0(group_name, " ", method_label)]] <- fmt_num(value, digits)
        }
      }
      out[[k]] <- item
    }
  }
  do.call(rbind, out)
}

latex_external_correlation_tabular <- function(df, correlation_groups, methods, method_labels_for_columns, resize = TRUE) {
  n_groups <- length(correlation_groups)
  n_methods <- length(methods)
  n_cols <- 1L + n_groups * n_methods
  group_align <- paste(vapply(seq_len(n_groups), function(i) {
    paste(rep("r", n_methods), collapse = "")
  }, character(1L)), collapse = "@{\\hspace{1em}}")
  align <- paste0("l", group_align)
  line_end <- paste0(" ", "\\", "\\")
  top_header <- c(
    "",
    vapply(names(correlation_groups), function(group_name) {
      paste0("\\multicolumn{", n_methods, "}{c}{", esc(group_name), "}")
    }, character(1L))
  )
  top_rules <- vapply(seq_along(correlation_groups), function(i) {
    start <- 2L + ((i - 1L) * n_methods)
    end <- start + n_methods - 1L
    paste0("\\cmidrule(lr){", start, "-", end, "}")
  }, character(1L))
  method_header <- c(
    "",
    unlist(lapply(seq_along(correlation_groups), function(group_i) {
      is_gee <- grepl("^gee_", methods)
      is_gamlss <- methods %in% c("rs_joint", "rs_joint_sandwich")
      has_gamlss_pair <- all(c("rs_joint", "rs_joint_sandwich") %in% methods)
      first_gee <- which(is_gee)[1L]
      last_gee <- utils::tail(which(is_gee), 1L)
      first_gamlss <- if (any(is_gamlss)) which(is_gamlss)[1L] else NA_integer_
      last_gamlss <- if (any(is_gamlss)) utils::tail(which(is_gamlss), 1L) else NA_integer_
      out <- character()
      i <- 1L
      while (i <= length(methods)) {
        if (isTRUE(has_gamlss_pair) && is.finite(first_gamlss) && i == first_gamlss) {
          out <- c(out, paste0("\\multicolumn{", sum(is_gamlss), "}{c}{gamlss.long}"))
          i <- last_gamlss + 1L
        } else if (i == first_gee) {
          out <- c(out, paste0("\\multicolumn{", sum(is_gee), "}{c}{geepack (robust)}"))
          i <- last_gee + 1L
        } else {
          out <- c(out, "")
          i <- i + 1L
        }
      }
      out
    }), use.names = FALSE)
  )
  method_rules <- unlist(lapply(seq_along(correlation_groups), function(group_i) {
    is_gee <- grepl("^gee_", methods)
    is_gamlss <- methods %in% c("rs_joint", "rs_joint_sandwich")
    has_gamlss_pair <- all(c("rs_joint", "rs_joint_sandwich") %in% methods)
    first_gee <- which(is_gee)[1L]
    last_gee <- utils::tail(which(is_gee), 1L)
    rules <- character()
    if (isTRUE(has_gamlss_pair)) {
      first_gamlss <- which(is_gamlss)[1L]
      last_gamlss <- utils::tail(which(is_gamlss), 1L)
      start <- 1L + ((group_i - 1L) * n_methods) + first_gamlss
      end <- 1L + ((group_i - 1L) * n_methods) + last_gamlss
      rules <- c(rules, paste0("\\cmidrule(lr){", start, "-", end, "}"))
    }
    start <- 1L + ((group_i - 1L) * n_methods) + first_gee
    end <- 1L + ((group_i - 1L) * n_methods) + last_gee
    c(rules, paste0("\\cmidrule(lr){", start, "-", end, "}"))
  }), use.names = FALSE)
  subtype_header <- c(
    "",
    unlist(lapply(seq_along(correlation_groups), function(group_i) {
      vapply(seq_along(methods), function(i) {
        has_gamlss_pair <- all(c("rs_joint", "rs_joint_sandwich") %in% methods)
        if (isTRUE(has_gamlss_pair) && identical(methods[[i]], "rs_joint")) {
          "naive"
        } else if (isTRUE(has_gamlss_pair) && identical(methods[[i]], "rs_joint_sandwich")) {
          "robust"
        } else {
          esc(method_labels_for_columns[[methods[[i]]]])
        }
      }, character(1L))
    }), use.names = FALSE)
  )
  tabular <- c(
    paste0("\\begin{tabular}{", align, "}"),
    "\\toprule",
    paste(top_header, collapse = " & ") |> paste0(line_end),
    paste(top_rules, collapse = " "),
    paste(method_header, collapse = " & ") |> paste0(line_end),
    paste(method_rules, collapse = " "),
    paste(subtype_header, collapse = " & ") |> paste0(line_end),
    "\\midrule"
  )
  value_cols <- setdiff(names(df), c("Level", "Metric"))
  for (metric in unique(df$Metric)) {
    tabular <- c(
      tabular,
      paste0("\\addlinespace[0.25em]"),
      paste0("\\multicolumn{", n_cols, "}{l}{\\textit{", esc(metric), "}}", line_end)
    )
    metric_df <- df[df$Metric == metric, c("Level", value_cols), drop = FALSE]
    for (i in seq_len(nrow(metric_df))) {
      vals <- esc(unlist(metric_df[i, ], use.names = FALSE))
      tabular <- c(tabular, paste(vals, collapse = " & ") |> paste0(line_end))
    }
  }
  tabular <- c(tabular, "\\bottomrule", "\\end{tabular}")
  if (!resize) return(tabular)
  c("\\resizebox{\\textwidth}{!}{%", tabular, "}")
}

write_external_correlation_table <- function(df, stem, caption, label, correlation_groups, methods, method_labels_for_columns, resize = FALSE) {
  csv_path <- file.path(out_dir, paste0(stem, ".csv"))
  tex_path <- file.path(out_dir, paste0(stem, ".tex"))
  utils::write.csv(df, csv_path, row.names = FALSE, na = "")
  tex <- c(
    paste0("\\begin{table}[htbp]"),
    "\\centering",
    paste0("\\caption{", caption, "}"),
    paste0("\\label{", label, "}"),
    latex_external_correlation_tabular(
      df,
      correlation_groups = correlation_groups,
      methods = methods,
      method_labels_for_columns = method_labels_for_columns,
      resize = resize
    ),
    "\\end{table}"
  )
  writeLines(tex, tex_path, useBytes = TRUE)
  invisible(list(csv = csv_path, tex = tex_path))
}

scenario_sort <- function(x) {
  order(
    x$n_time,
    match(base_scenario(x$scenario), names(scenario_labels)),
    x$scenario,
    na.last = TRUE
  )
}

review <- add_n_time(read_sandwich_csv("full_summary_review_table.csv", "augmented_full_summary_review_table.csv"))
dependence <- add_n_time(read_sandwich_csv("dependence_recovery_summary.csv", "augmented_dependence_recovery_summary.csv", required = FALSE))
complexity <- add_n_time(read_sandwich_csv("fit_complexity_summary.csv", "augmented_fit_complexity_summary.csv", required = FALSE))
status <- add_n_time(read_sandwich_csv("primary_status_summary.csv", "augmented_primary_status_summary.csv", required = FALSE))
benchmark_by_rep <- add_n_time(read_csv("benchmark_results_by_rep.csv", required = FALSE))

review <- review[review$method %in% method_order, , drop = FALSE]
dependence <- dependence[dependence$method %in% method_order, , drop = FALSE]
complexity <- complexity[complexity$method %in% method_order, , drop = FALSE]
review <- review[review$family %in% family_order, , drop = FALSE]
dependence <- dependence[dependence$family %in% family_order, , drop = FALSE]
complexity <- complexity[complexity$family %in% family_order, , drop = FALSE]
status <- status[status$family %in% family_order, , drop = FALSE]
benchmark_by_rep <- benchmark_by_rep[benchmark_by_rep$family %in% family_order, , drop = FALSE]

# Table 1: simulation design.
scenario_specs <- bmk_scenario_specs()
included <- unique(review[c("n_time", "scenario")])
included$scenario_base <- base_scenario(included$scenario)
included <- included[included$scenario_base %in% names(scenario_specs), , drop = FALSE]
if (nrow(included) == 0L) {
  stop("No recognised simulation scenarios found in the selected run directory.", call. = FALSE)
}
included_base <- unique(included["scenario_base"])
included_base$sort_order <- match(included_base$scenario_base, names(scenario_labels))
included_base <- included_base[order(included_base$sort_order), , drop = FALSE]
scenario_rho_label <- function(spec) {
  if (identical(spec$generator, "external")) return(fmt_num(spec$rho, 2L))
  if (!is.null(spec$tau_edges)) {
    return(paste0(fmt_num(min(spec$tau_edges), 2L), "--", fmt_num(max(spec$tau_edges), 2L)))
  }
  if (!is.null(spec$tau_base) && !is.null(spec$tau_effect)) {
    return(paste0(
      "$\\operatorname{logit}^{-1}\\{\\operatorname{logit}(",
      fmt_num(spec$tau_base, 2L),
      ") + ",
      fmt_num(spec$tau_effect, 2L),
      "x\\}$"
    ))
  }
  "--"
}
scenario_rows <- lapply(seq_len(nrow(included_base)), function(i) {
  scenario_base <- included_base$scenario_base[[i]]
  spec <- scenario_specs[[scenario_base]]
  t_values <- sort(unique(included$n_time[included$scenario_base == scenario_base]))
  dependence <- unname(dependence_labels[spec$correlation])
  if (is.na(dependence)) dependence <- spec$correlation
  data.frame(
    Scenario = scenario_ids[[scenario_base]],
    n = n_subjects,
    T = paste(t_values, collapse = ", "),
    Dependence = dependence,
    `$\\rho$ / $\\tau$` = scenario_rho_label(spec),
    check.names = FALSE,
    stringsAsFactors = FALSE
  )
})
table_1 <- do.call(rbind, scenario_rows)
write_table(
  table_1,
  "table_1_scenario_design",
  "Simulation scenarios included in the benchmark. All scenarios used Normal, Gamma, and Binary response margins with a Gaussian copula; Poisson margins are excluded from these paper-facing tables. External scenarios are parameterised by Gaussian copula correlation $\\rho$, while internal time-varying and covariate-dependent scenarios are parameterised on Kendall's $\\tau$ scale.",
  "tab:scenario-design",
  align = "lrrll",
  resize = FALSE,
  escape = FALSE
)

# Tables 2-4: combined performance and dependence recovery, split by generated correlation structure.
main_methods <- c("rs_joint", "rs_joint_sandwich", "gee_exchangeable", "gee_ar1", "gee_unstructured", "glm")
appendix_methods <- c("rs_joint", "gee_exchangeable", "gee_ar1", "gee_unstructured", "glm")
dependence_methods <- c("rs_joint", "rs_joint_sandwich", "gee_exchangeable", "gee_ar1", "gee_unstructured")
appendix_dependence_methods <- c("rs_joint", "gee_exchangeable", "gee_ar1", "gee_unstructured")
robust_caption_note <- paste(
  "The robust column uses the same gamlss.longitudinal fit as gamlss.long",
  "but replaces model-based coefficient SEs with subject-cluster sandwich SEs;",
  "fit and dependence metrics are therefore unchanged by design."
)
combined_metrics <- c(
  benchmark_mean_rmse = "Mean RMSE (lower is better)",
  x_se_ratio = "Slope SE ratio (target: 1.0)",
  x_ci_coverage = "Slope coverage (target: 0.95)",
  tau_mae_median = "Tau MAE",
  tau_rmse_median = "Tau RMSE"
)

scenario_method_grid <- merge(
  unique(review[c("n_time", "scenario")]),
  data.frame(method = main_methods, stringsAsFactors = FALSE),
  by = NULL
)
main_summary <- aggregate(
  cbind(benchmark_mean_rmse, x_se_ratio, x_ci_coverage) ~ n_time + scenario + method,
  review[review$method %in% main_methods, , drop = FALSE],
  FUN = function(x) mean(x, na.rm = TRUE)
)
dep_df <- dependence[is.finite(dependence$tau_mae_median) | is.finite(dependence$tau_rmse_median), , drop = FALSE]
dep_summary <- aggregate(
  cbind(tau_mae_median, tau_rmse_median) ~ n_time + scenario + method,
  dep_df,
  FUN = function(x) mean(x, na.rm = TRUE)
)
fit_summary <- aggregate(
  cbind(total_df_median, dependence_df_median) ~ n_time + scenario + method,
  complexity,
  FUN = function(x) mean(x, na.rm = TRUE)
)
combined_summary <- merge(scenario_method_grid, main_summary, by = c("n_time", "scenario", "method"), all.x = TRUE)
combined_summary <- merge(combined_summary, dep_summary, by = c("n_time", "scenario", "method"), all.x = TRUE)
combined_summary <- merge(combined_summary, fit_summary, by = c("n_time", "scenario", "method"), all.x = TRUE)
combined_summary$Scenario <- scenario_code(combined_summary$scenario)
combined_summary$T <- combined_summary$n_time
combined_summary$method <- factor(combined_summary$method, levels = method_order)
combined_summary <- combined_summary[
  order(combined_summary$n_time, match(base_scenario(combined_summary$scenario), names(scenario_labels)), combined_summary$method),
  ,
  drop = FALSE
]

write_combined_correlation_table <- function(
    scenario_bases,
    stem,
    caption,
    label,
    times,
    time_caption = NULL,
    metrics = combined_metrics,
    methods = main_methods,
    dependence_methods_for_note = dependence_methods,
    extra_caption_note = NULL) {
  tab_df <- combined_summary[base_scenario(combined_summary$scenario) %in% scenario_bases, , drop = FALSE]
  tab_df <- tab_df[tab_df$T %in% times, , drop = FALSE]
  tab_df <- tab_df[order(match(base_scenario(tab_df$scenario), scenario_bases), tab_df$T, tab_df$method), , drop = FALSE]
  table <- pivot_metric_time_table(
    tab_df,
    row_cols = c("Scenario"),
    time_col = "T",
    metrics = metrics,
    methods = methods,
    method_labels_for_columns = method_column_labels,
    times = times
  )
  df_note <- df_range_note_by_method(
    tab_df[tab_df$method %in% dependence_methods_for_note, , drop = FALSE],
    methods = dependence_methods_for_note,
    labels = c(
      rs_joint = "gamlss.longitudinal",
      rs_joint_sandwich = "gamlss.longitudinal robust",
      gee_exchangeable = "geepack exch.",
      gee_ar1 = "geepack AR(1)",
      gee_unstructured = "geepack unstr."
    )
  )
  time_note <- if (is.null(time_caption)) "" else paste0(" ", time_caption)
  robust_note <- if (is.null(extra_caption_note)) "" else paste0(" ", extra_caption_note)
  write_time_pivot_table(
    table,
    stem,
    paste0(caption, time_note, " Values are averaged over response families. Scenario labels are defined in Table~\\ref{tab:scenario-design}.", robust_note, df_note),
    label,
    row_cols = c("Scenario"),
    methods = methods,
    method_labels_for_columns = method_column_labels,
    times = times,
    resize = FALSE
  )
}

external_correlation_groups <- list(
  `AR(1)` = c(
    moderate = "external_ar1_moderate",
    high = "external_ar1_high"
  ),
  Exchangeable = c(
    moderate = "external_exchangeable_moderate",
    high = "external_exchangeable_high"
  )
)
external_levels <- c(moderate = "Moderate", high = "High")

external_t20 <- combined_summary[
  base_scenario(combined_summary$scenario) %in% unlist(external_correlation_groups) &
    combined_summary$T == 20L,
  ,
  drop = FALSE
]
external_table_20 <- pivot_external_correlation_table(
  external_t20,
  correlation_groups = external_correlation_groups,
  levels = external_levels,
  time_col = "T",
  metrics = combined_metrics,
  methods = main_methods,
  method_labels_for_columns = method_column_labels,
  time = 20L
)
external_df_note_20 <- df_range_note_by_method(
  external_t20[external_t20$method %in% dependence_methods, , drop = FALSE],
  methods = dependence_methods,
  labels = c(
    rs_joint = "gamlss.longitudinal",
    rs_joint_sandwich = "gamlss.longitudinal robust",
    gee_exchangeable = "geepack exch.",
    gee_ar1 = "geepack AR(1)",
    gee_unstructured = "geepack unstr."
  )
)
write_external_correlation_table(
  external_table_20,
  "table_2_external_correlation",
  paste0("Performance and all-pair dependence recovery under generated AR(1) and exchangeable correlation. Results are shown for T=20. Values are averaged over response families. Scenario definitions are given in Table~\\ref{tab:scenario-design}. ", robust_caption_note, external_df_note_20),
  "tab:external-correlation",
  correlation_groups = external_correlation_groups,
  methods = main_methods,
  method_labels_for_columns = method_column_labels,
  resize = FALSE
)
write_combined_correlation_table(
  c("internal_time_varying_high", "internal_covariate_dependent_high"),
  "table_3_flexible_correlation",
  "Performance and all-pair dependence recovery under generated time-varying and covariate-dependent correlation.",
  "tab:flexible-correlation",
  times = 20L,
  time_caption = "Results are shown for T=20.",
  extra_caption_note = robust_caption_note
)

external_t50 <- combined_summary[
  base_scenario(combined_summary$scenario) %in% unlist(external_correlation_groups) &
    combined_summary$T == 50L,
  ,
  drop = FALSE
]
external_table_50 <- pivot_external_correlation_table(
  external_t50,
  correlation_groups = external_correlation_groups,
  levels = external_levels,
  time_col = "T",
  metrics = if (any(is.finite(external_t50$tau_mae_median) | is.finite(external_t50$tau_rmse_median))) combined_metrics else combined_metrics[!names(combined_metrics) %in% c("tau_mae_median", "tau_rmse_median")],
  methods = appendix_methods,
  method_labels_for_columns = method_column_labels,
  time = 50L
)
external_caption_50 <- if (any(is.finite(external_t50$tau_mae_median) | is.finite(external_t50$tau_rmse_median))) {
  "Appendix performance and all-pair dependence recovery under generated AR(1) and exchangeable correlation."
} else {
  "Appendix performance under generated AR(1) and exchangeable correlation."
}
external_df_note_50 <- df_range_note_by_method(
  external_t50[external_t50$method %in% appendix_dependence_methods, , drop = FALSE],
  methods = appendix_dependence_methods,
  labels = c(
    rs_joint = "gamlss.longitudinal",
    gee_exchangeable = "geepack exch.",
    gee_ar1 = "geepack AR(1)",
    gee_unstructured = "geepack unstr."
  )
)
write_external_correlation_table(
  external_table_50,
  "appendix_table_1_external_correlation_t50",
  paste0(external_caption_50, " Results are shown for T=50. Values are averaged over response families. Scenario definitions are given in Table~\\ref{tab:scenario-design}.", external_df_note_50),
  "tab:appendix-external-correlation-t50",
  correlation_groups = external_correlation_groups,
  methods = appendix_methods,
  method_labels_for_columns = method_column_labels,
  resize = FALSE
)
write_combined_correlation_table(
  c("internal_time_varying_high", "internal_covariate_dependent_high"),
  "appendix_table_2_flexible_correlation_t50",
  if (any(is.finite(combined_summary$tau_mae_median[combined_summary$T == 50L & base_scenario(combined_summary$scenario) %in% c("internal_time_varying_high", "internal_covariate_dependent_high")]) |
          is.finite(combined_summary$tau_rmse_median[combined_summary$T == 50L & base_scenario(combined_summary$scenario) %in% c("internal_time_varying_high", "internal_covariate_dependent_high")]))) {
    "Appendix performance and all-pair dependence recovery under generated time-varying and covariate-dependent correlation."
  } else {
    "Appendix performance under generated time-varying and covariate-dependent correlation."
  },
  "tab:appendix-flexible-correlation-t50",
  times = 50L,
  time_caption = "Results are shown for T=50.",
  metrics = if (any(is.finite(combined_summary$tau_mae_median[combined_summary$T == 50L & base_scenario(combined_summary$scenario) %in% c("internal_time_varying_high", "internal_covariate_dependent_high")]) |
                   is.finite(combined_summary$tau_rmse_median[combined_summary$T == 50L & base_scenario(combined_summary$scenario) %in% c("internal_time_varying_high", "internal_covariate_dependent_high")]))) combined_metrics else combined_metrics[!names(combined_metrics) %in% c("tau_mae_median", "tau_rmse_median")]
  ,
  methods = appendix_methods,
  dependence_methods_for_note = appendix_dependence_methods
)

combined <- c(
  "% Automatically generated by 07-write-story-tables.R",
  paste0("% Source run directory: ", run_dir),
  "",
  readLines(file.path(out_dir, "table_1_scenario_design.tex"), warn = FALSE),
  "",
  readLines(file.path(out_dir, "table_2_external_correlation.tex"), warn = FALSE),
  "",
  readLines(file.path(out_dir, "table_3_flexible_correlation.tex"), warn = FALSE)
)
writeLines(combined, file.path(out_dir, "story_tables.tex"), useBytes = TRUE)

appendix_t50 <- c(
  "% Automatically generated by 07-write-story-tables.R",
  paste0("% Source run directory: ", run_dir),
  "",
  readLines(file.path(out_dir, "appendix_table_1_external_correlation_t50.tex"), warn = FALSE),
  "",
  readLines(file.path(out_dir, "appendix_table_2_flexible_correlation_t50.tex"), warn = FALSE)
)
writeLines(appendix_t50, file.path(out_dir, "appendix_t50_tables.tex"), useBytes = TRUE)

message("Story tables written to: ", out_dir)
