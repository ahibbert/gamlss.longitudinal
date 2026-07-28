#' Clean a model-matrix term label for fixed-term factor matching
#'
#' @noRd
.plot_fixed_terms_clean_expr_name <- function(x) {
  x <- trimws(x)
  x <- gsub("`", "", x, fixed = TRUE)

  factor_match <- regexec("^(?:as\\.)?factor\\(([^)]+)\\)$", x)
  matched <- regmatches(x, factor_match)[[1]]
  if (length(matched) >= 2) {
    x <- trimws(matched[2])
  }

  x
}

#' Escape a string for use inside a regular expression
#'
#' @noRd
.plot_fixed_terms_regex_escape <- function(x) {
  gsub("([][{}()+*^$|\\\\?.])", "\\\\\\1", x)
}

#' Build a factor group from known design-matrix columns
#'
#' @noRd
.plot_fixed_terms_make_factor_group_from_cols <- function(term_name, term_cols, levs) {
  if (length(levs) < 2 || length(term_cols) == 0) {
    return(NULL)
  }

  level_col_map <- list()
  matched_cols <- character(0)

  for (lev in levs[-1]) {
    lev_plain <- as.character(lev)
    lev_mn <- make.names(lev_plain)
    hits <- term_cols[
      endsWith(term_cols, lev_plain) |
        endsWith(term_cols, lev_mn) |
        grepl(
          paste0("(^|[^[:alnum:]_.])", .plot_fixed_terms_regex_escape(lev_plain), "$"),
          term_cols
        )
    ]

    if (length(hits) == 0 && length(term_cols) == length(levs) - 1L) {
      hits <- term_cols[seq_along(levs[-1]) == which(levs[-1] == lev)]
    }

    if (length(hits) > 0) {
      level_col_map[[lev]] <- hits[1]
      matched_cols <- c(matched_cols, hits[1])
    }
  }

  if (length(level_col_map) == 0) {
    return(NULL)
  }

  list(
    var_name = term_name,
    levels = levs,
    ref_level = levs[1],
    level_col_map = level_col_map,
    matched_cols = unique(matched_cols)
  )
}

#' Candidate data columns for a fixed-term factor expression
#'
#' @noRd
.plot_fixed_terms_data_candidates <- function(term_var, data) {
  unique(c(
    term_var,
    make.names(term_var),
    if (identical(term_var, "time_covariate")) {
      names(data)[grepl("time", names(data), ignore.case = TRUE)]
    } else {
      character(0)
    }
  ))
}

#' Recover factor levels from the original data for a model-matrix term
#'
#' @noRd
.plot_fixed_terms_levels_from_data <- function(term_name, term_var, data) {
  if (is.null(data) || !is.data.frame(data)) {
    return(NULL)
  }

  data_candidates <- .plot_fixed_terms_data_candidates(term_var, data)

  for (candidate in data_candidates) {
    if (candidate %in% names(data) && is.factor(data[[candidate]])) {
      return(levels(data[[candidate]]))
    } else if (candidate %in% names(data) && grepl("^(?:as\\.)?factor\\(", term_name)) {
      return(levels(as.factor(data[[candidate]])))
    }
  }

  NULL
}

#' Recover factor levels directly from factor() design columns
#'
#' @noRd
.plot_fixed_terms_levels_from_factor_expr <- function(term_name, term_cols) {
  if (!grepl("^(?:as\\.)?factor\\(", term_name)) {
    return(NULL)
  }

  levs <- sub(paste0("^", .plot_fixed_terms_regex_escape(term_name)), "", term_cols)
  c(sub("\\).*", ")", term_name), levs)
}

#' Candidate design-matrix prefixes for a data factor
#'
#' @noRd
.plot_fixed_terms_var_prefixes <- function(var_name, x_cols) {
  var_tokens <- strsplit(var_name, "_", fixed = TRUE)[[1]]
  var_prefixes <- unique(c(
    var_name,
    make.names(var_name),
    if (length(var_tokens) > 0) var_tokens[1] else character(0),
    if (length(var_tokens) > 0) make.names(var_tokens[1]) else character(0)
  ))

  # Internal fitting often renames the user time variable to time_covariate.
  if (grepl("time", var_name, ignore.case = TRUE) && any(grepl("^time_covariate", x_cols))) {
    var_prefixes <- unique(c(var_prefixes, "time_covariate", make.names("time_covariate")))
  }

  var_prefixes[nzchar(var_prefixes)]
}

#' Match factor levels to design columns using prefix candidates
#'
#' @noRd
.plot_fixed_terms_factor_level_col_map <- function(levs, var_prefixes, x_cols) {
  level_col_map <- list()
  matched_cols <- character(0)

  for (lev in levs[-1]) {
    lev_plain <- as.character(lev)
    lev_mn <- make.names(lev_plain)
    candidates <- unique(unlist(lapply(var_prefixes, function(pref) {
      c(
        paste0(pref, lev_plain),
        paste0(pref, lev_mn),
        paste0(pref, "_", lev_plain),
        paste0(pref, "_", lev_mn)
      )
    }), use.names = FALSE))

    hit <- candidates[candidates %in% x_cols]
    if (length(hit) > 0) {
      level_col_map[[lev]] <- hit[1]
      matched_cols <- c(matched_cols, hit[1])
    }
  }

  list(
    level_col_map = level_col_map,
    matched_cols = matched_cols
  )
}

#' Build a fallback factor group by scanning original data columns
#'
#' @noRd
.plot_fixed_terms_fallback_factor_group <- function(var_name, data, x_cols, groups) {
  v <- data[[var_name]]
  if (!is.factor(v)) {
    return(NULL)
  }

  levs <- levels(v)
  if (length(levs) < 2) {
    return(NULL)
  }

  var_prefixes <- .plot_fixed_terms_var_prefixes(var_name, x_cols)
  matched <- .plot_fixed_terms_factor_level_col_map(levs, var_prefixes, x_cols)
  level_col_map <- matched$level_col_map
  matched_cols <- matched$matched_cols

  if (length(matched_cols) == 0 && length(levs) == 2 && var_name %in% x_cols) {
    level_col_map[[levs[2]]] <- var_name
    matched_cols <- var_name
  }

  existing_cols <- unlist(lapply(groups, function(g) g$matched_cols), use.names = FALSE)
  if (length(level_col_map) == 0 || any(matched_cols %in% existing_cols)) {
    return(NULL)
  }

  list(
    var_name = var_name,
    levels = levs,
    ref_level = levs[1],
    level_col_map = level_col_map,
    matched_cols = unique(matched_cols)
  )
}
